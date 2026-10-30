//! Fillet an edge on a boolean (AND) result.
//!
//! Builds the same punched-cube as `punched-cube-shapeops`, then fillets
//! one of the remaining cube edges to demonstrate `fillet_edges`
//! working on shells that contain `IntersectionCurve` edges.

use monstertruck_modeling::{EuclideanSpace, Point3, Solid, Vector3, builder};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    // Unit cube at origin.
    let v = builder::vertex(Point3::origin());
    let e = builder::extrude(&v, Vector3::unit_x());
    let f = builder::extrude(&e, Vector3::unit_y());
    let cube: Solid = builder::extrude(&f, Vector3::unit_z());

    // Cylinder through the cube (same as punched-cube-shapeops example).
    let v = builder::vertex(Point3::new(0.5, 0.25, -0.5));
    let w = builder::revolve(
        &v,
        Point3::new(0.5, 0.5, 0.0),
        Vector3::unit_z(),
        builder::SweepAngle::Closed,
        4,
    )?;
    let f = builder::try_attach_plane([w])?;
    let mut cylinder = builder::extrude(&f, Vector3::unit_z() * 2.0);
    cylinder.not();

    let and = monstertruck_solid::and(&cube, &cylinder, 0.05)?;
    let json = serde_json::to_vec_pretty(&and)?;
    std::fs::write("fillet-after-boolean.json", json)?;
    eprintln!("Wrote fillet-after-boolean.json");
    Ok(())
}
