//! A cube punched by a torus.
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::{
    EuclideanSpace, Face, ParametricSurface3D, Point3, Rad, Result, Shell, Solid, Tolerance,
    Vector3, Wire, builder,
};
use std::f64::consts::PI;

type ExampleResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn cube_shell() -> ExampleResult<Shell> {
    let v = builder::vertex(Point3::origin());
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    let s = builder::extrude(&f, Vector3::unit_z());
    s.into_boundaries()
        .pop()
        .ok_or_else(|| std::io::Error::other("cube extrusion produced no shell").into())
}

fn torus_shell() -> Result<Shell> {
    let v = builder::vertex(Point3::new(0.5, 0.0, 0.25));
    let w = builder::revolve(
        &v,
        Point3::new(0.5, 0.0, 0.5),
        -Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    )?;
    builder::revolve(
        &w,
        Point3::origin(),
        Vector3::unit_z(),
        builder::SweepAngle::Partial(Rad(PI / 2.0)),
        1,
    )
}

fn find_cube_face(cube_shell: &mut Shell, normal: Vector3) -> Option<&mut Face> {
    cube_shell.iter_mut().find(|face| {
        let surface = face.oriented_surface();
        surface.normal(0.5, 0.5).near(&normal)
    })
}

fn find_torus_boundary(bdds: &[Wire], idx: usize) -> Option<&Wire> {
    bdds.iter().find(|wire| {
        wire.front_vertex()
            .is_some_and(|vertex| vertex.point()[idx].near(&0.0))
    })
}

fn main() -> ExampleResult<()> {
    let mut shell = cube_shell()?;
    let torus = torus_shell()?;
    let bdds = torus.extract_boundaries();
    let face = find_cube_face(&mut shell, -Vector3::unit_y())
        .ok_or_else(|| std::io::Error::other("cube has no -Y face"))?;
    let bdd = find_torus_boundary(&bdds, 1)
        .ok_or_else(|| std::io::Error::other("torus has no Y-plane boundary"))?;
    face.add_boundary(bdd.inverse())?;
    let face = find_cube_face(&mut shell, -Vector3::unit_x())
        .ok_or_else(|| std::io::Error::other("cube has no -X face"))?;
    let bdd = find_torus_boundary(&bdds, 0)
        .ok_or_else(|| std::io::Error::other("torus has no X-plane boundary"))?;
    face.add_boundary(bdd.inverse())?;
    shell.extend(torus);
    let solid = Solid::new_unchecked(vec![shell]);
    let json = serde_json::to_vec_pretty(&solid)?;
    std::fs::write("torus-punched-cube.json", json)?;
    Ok(())
}
