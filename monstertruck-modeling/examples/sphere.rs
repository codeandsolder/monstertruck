//! Modeling a sphere
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::{Point3, Solid, Vector3, Wire, builder};
use std::f64::consts::PI;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let v0 = builder::vertex(Point3::new(0.0, 0.5, 0.0));
    let wire: Wire = builder::revolve(
        &v0,
        Point3::origin(),
        Vector3::x(),
        builder::SweepAngle::Partial(PI),
        3,
    )?;
    let shell = builder::revolve_wire(
        &wire,
        Point3::origin(),
        Vector3::y(),
        builder::SweepAngle::Closed,
        4,
    )?;
    let sphere = Solid::new_unchecked(vec![shell]);
    let json = serde_json::to_vec_pretty(&sphere)?;
    std::fs::write("sphere.json", json)?;
    Ok(())
}
