//! Modeling a unit cube with a hole through it.
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::builder::SweepAngle;
use monstertruck_modeling::{Edge, EuclideanSpace, Point3, Rad, Solid, Vector3, Wire, builder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let v = builder::vertex(Point3::new(-0.5, -0.5, -0.5));
    let edge = builder::extrude(&v, Vector3::unit_x());
    let mut face = builder::extrude(&edge, Vector3::unit_y());
    let v = builder::vertex(Point3::new(0.2, 0.0, -0.5));
    let edge0 = builder::extrude(&v, Vector3::new(-0.2, 0.2, 0.0));
    let wire1 = builder::revolve(
        edge0.back(),
        Point3::origin(),
        Vector3::unit_z(),
        SweepAngle::Partial(Rad(std::f64::consts::PI / 2.0)),
        2,
    )?;
    let back_vertex = wire1
        .back_vertex()
        .ok_or_else(|| std::io::Error::other("first quarter-turn produced no back vertex"))?;
    let edge2 = builder::extrude(back_vertex, Vector3::new(0.2, -0.2, 0.0));
    let mut wire3 = builder::revolve(
        edge2.back(),
        Point3::origin(),
        Vector3::unit_z(),
        SweepAngle::Partial(Rad(std::f64::consts::PI / 2.0)),
        2,
    )?;
    let back_edge = wire3
        .pop_back()
        .ok_or_else(|| std::io::Error::other("second quarter-turn produced no closing edge"))?;
    let tmp = Edge::new_unchecked(back_edge.front(), edge0.front(), back_edge.curve());
    wire3.push_back(tmp);
    let mut wire = Wire::from_iter(
        std::iter::once(edge0)
            .chain(wire1)
            .chain(std::iter::once(edge2))
            .chain(wire3),
    );
    wire.invert();
    face.add_boundary(wire)?;
    let shape: Solid = builder::extrude(&face, Vector3::unit_z());
    let json = serde_json::to_vec_pretty(&shape)?;
    std::fs::write("punched-cube.json", json)?;
    Ok(())
}
