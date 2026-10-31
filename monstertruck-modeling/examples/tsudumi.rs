//! Modeling a one-leaf hyperboloid.
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::{EuclideanSpace, Point3, Solid, Vector3, builder};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let v0 = builder::vertex(Point3::new(1.0, 1.0, 0.0));
    let v1 = builder::vertex(Point3::new(0.0, -1.0, 1.0));
    let line = builder::line(&v0, &v1);
    let mut shell = builder::revolve(
        &line,
        Point3::origin(),
        Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    )?;
    let wires = shell.extract_boundaries();
    shell.push(builder::try_attach_plane([wires[0].clone()])?.inverse());
    shell.push(builder::try_attach_plane([wires[1].clone()])?.inverse());
    let solid = Solid::new_unchecked(vec![shell]);
    let json = serde_json::to_vec_pretty(&solid)?;
    std::fs::write("tsudumi.json", json)?;
    Ok(())
}
