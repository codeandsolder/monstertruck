//! Modeling a cylinder by two sweeps.
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::{EuclideanSpace, Point3, Result, Solid, Vector3, builder};

fn cylinder(height: f64, radius: f64) -> Result<Solid> {
    let vertex = builder::vertex(Point3::new(0.0, -height / 2.0, radius));
    let circle = builder::revolve(
        &vertex,
        Point3::origin(),
        Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    )?;
    let disk = builder::try_attach_plane([circle])?;

    Ok(builder::extrude(&disk, Vector3::new(0.0, height, 0.0)))
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let cylinder = cylinder(1.0, 0.5)?;
    let json = serde_json::to_vec_pretty(&cylinder)?;
    std::fs::write("cylinder.json", json)?;
    Ok(())
}
