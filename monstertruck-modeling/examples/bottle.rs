//! Modeling a bottle.
//!
//! This is a technical indicator for comparing with Open CASCADE Technology, a great senior.
//! We want to reproduce the bottle made in the [OCCT tutorial].
//!
//! When the `fillet` feature is enabled, the body edges are filleted just like
//! the OCCT tutorial (`BRepFilletAPI_MakeFillet` at radius = thickness / 12).
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.
//!
//! ```bash
//! cargo run -p monstertruck-modeling --features fillet --example bottle
//! ```
//!
//! [OCCT tutorial]: https://dev.opencascade.org/doc/overview/html/occt__tutorial.html

#[cfg(feature = "fillet")]
use monstertruck_modeling::{Curve, Edge, FilletOptions, RadiusSpec, fillet_edges};
use monstertruck_modeling::{EuclideanSpace, Point3, Rad, Shell, Solid, Vector3, builder};
use std::f64::consts::PI;

type ExampleResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn body_shell(bottom: f64, height: f64, width: f64, thickness: f64) -> ExampleResult<Shell> {
    let vertex0 = builder::vertex(Point3::new(-width / 2.0, bottom, thickness / 4.0));
    let vertex1 = builder::vertex(Point3::new(width / 2.0, bottom, thickness / 4.0));
    let transit = Point3::new(0.0, bottom, thickness / 2.0);
    let arc0 = builder::circle_arc(&vertex0, &vertex1, transit)?;
    let arc1 = builder::rotated(&arc0, Point3::origin(), Vector3::unit_y(), Rad(PI));
    let face = builder::homotopy(&arc0, &arc1.inverse());
    let solid = builder::extrude(&face, Vector3::new(0.0, height, 0.0));
    solid
        .into_boundaries()
        .pop()
        .ok_or_else(|| std::io::Error::other("bottle body extrusion produced no shell").into())
}

fn cylinder(bottom: f64, height: f64, radius: f64) -> ExampleResult<Shell> {
    let vertex = builder::vertex(Point3::new(0.0, bottom, radius));
    let circle = builder::revolve(
        &vertex,
        Point3::origin(),
        Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    )?;
    let disk = builder::try_attach_plane([circle])?;
    let solid = builder::extrude(&disk, Vector3::new(0.0, height, 0.0));
    solid
        .into_boundaries()
        .pop()
        .ok_or_else(|| std::io::Error::other("neck extrusion produced no shell").into())
}

fn glue_body_neck(body: &mut Shell, neck: Shell) -> ExampleResult<()> {
    let neck_face = neck
        .first()
        .ok_or_else(|| std::io::Error::other("neck shell has no attachment face"))?;
    let wire = neck_face
        .boundaries()
        .into_iter()
        .next()
        .ok_or_else(|| std::io::Error::other("neck shell has no attachment boundary"))?;
    let body_ceiling = body
        .last_mut()
        .ok_or_else(|| std::io::Error::other("bottle body has no ceiling face"))?;
    body_ceiling.add_boundary(wire)?;
    body.extend(neck.into_iter().skip(1));
    Ok(())
}

#[cfg(feature = "fillet")]
fn vertical_line_edges(shell: &Shell, height: f64) -> Vec<Edge> {
    let mut edge_face_count = std::collections::HashMap::new();
    for face in shell.face_iter() {
        for wire in face.boundaries() {
            for edge in wire.edge_iter() {
                *edge_face_count.entry(edge.id()).or_insert(0u32) += 1;
            }
        }
    }
    let mut seen = std::collections::HashSet::new();
    shell
        .edge_iter()
        .filter(|edge| seen.insert(edge.id()))
        .filter(|edge| {
            edge_face_count
                .get(&edge.id())
                .is_some_and(|&count| count == 2)
        })
        .filter(|edge| {
            matches!(edge.oriented_curve(), Curve::Line(_)) && {
                let p0 = edge.front().point();
                let p1 = edge.back().point();
                (p0.y.abs() < 0.01 && (p1.y - height).abs() < 0.01)
                    || (p1.y.abs() < 0.01 && (p0.y - height).abs() < 0.01)
            }
        })
        .collect()
}

fn bottle(height: f64, width: f64, thickness: f64) -> ExampleResult<Solid> {
    let mut body = body_shell(0.0, height, width, thickness)?;

    #[cfg(feature = "fillet")]
    {
        let edges = vertical_line_edges(&body, height);
        let opts = FilletOptions {
            radius: RadiusSpec::Constant(thickness / 12.0),
            ..Default::default()
        };
        fillet_edges(&mut body, &edges, Some(&opts))?;
    }

    let neck = cylinder(height, height / 10.0, thickness / 4.0)?;
    glue_body_neck(&mut body, neck)?;

    let eps = height / 50.0;
    let inner_height = 2.0f64.mul_add(-eps, height);
    let inner_width = 2.0f64.mul_add(-eps, width);
    let inner_thickness = 2.0f64.mul_add(-eps, thickness);
    let mut inner_body = body_shell(eps, inner_height, inner_width, inner_thickness)?;

    #[cfg(feature = "fillet")]
    {
        let edges = vertical_line_edges(&inner_body, inner_height);
        let opts = FilletOptions {
            radius: RadiusSpec::Constant(inner_thickness / 12.0),
            ..Default::default()
        };
        fillet_edges(&mut inner_body, &edges, Some(&opts))?;
    }

    let inner_neck = cylinder(height - eps, height / 10.0 + eps, thickness / 4.0 - eps)?;
    glue_body_neck(&mut inner_body, inner_neck)?;

    let inner_hat = inner_body
        .pop()
        .ok_or_else(|| std::io::Error::other("inner bottle body has no cap"))?;
    let wire = inner_hat
        .into_boundaries()
        .into_iter()
        .next()
        .ok_or_else(|| std::io::Error::other("inner bottle cap has no boundary"))?
        .inverse();
    body.last_mut()
        .ok_or_else(|| std::io::Error::other("bottle body has no final face"))?
        .add_boundary(wire)?;
    body.extend(inner_body.into_iter().map(|face| face.inverse()));
    Ok(Solid::new_unchecked(vec![body]))
}

fn main() -> ExampleResult<()> {
    let bottle = bottle(1.4, 1.0, 0.6)?;
    let json = serde_json::to_vec_pretty(&bottle)?;
    std::fs::write("bottle.json", json)?;
    Ok(())
}
