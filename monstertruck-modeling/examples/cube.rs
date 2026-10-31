//! Modeling a unit cube by three sweeps.
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::{Point3, Solid, Vector3, builder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let v = builder::vertex(Point3::new(-0.5, -0.5, -0.5));
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    let cube: Solid = builder::extrude(&f, Vector3::unit_z());
    let json = serde_json::to_vec_pretty(&cube)?;
    std::fs::write("cube.json", json)?;
    Ok(())
}
