//! Demonstrates `solid_from_planar_profile` with a complex profile containing
//! nested holes.
//!
//! Creates a rectangular plate with two rectangular holes punched through it,
//! showing how the profile pipeline automatically classifies outer vs. hole
//! loops and normalizes winding directions.

use monstertruck_modeling::{Point3, Solid, Vector3, Wire, builder, profile};

fn rect_wire(x0: f64, y0: f64, x1: f64, y1: f64) -> Wire {
    let v0 = builder::vertex(Point3::new(x0, y0, 0.0));
    let v1 = builder::vertex(Point3::new(x1, y0, 0.0));
    let v2 = builder::vertex(Point3::new(x1, y1, 0.0));
    let v3 = builder::vertex(Point3::new(x0, y1, 0.0));
    vec![
        builder::line(&v0, &v1),
        builder::line(&v1, &v2),
        builder::line(&v2, &v3),
        builder::line(&v3, &v0),
    ]
    .into()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let outer = rect_wire(-5.0, -3.0, 5.0, 3.0);
    let hole_left = rect_wire(-4.0, -2.0, -1.0, 2.0);
    let hole_right = rect_wire(1.0, -2.0, 4.0, 2.0);

    let solid: Solid = profile::solid_from_planar_profile(
        vec![outer, hole_left, hole_right],
        Vector3::new(0.0, 0.0, 1.0),
    )?;

    assert!(solid.is_geometric_consistent());
    assert_eq!(solid.boundaries()[0].len(), 14);

    let json = serde_json::to_vec_pretty(&solid)?;
    std::fs::write("profile-with-holes.json", json)?;
    Ok(())
}
