//! An example of the solid with several boundaries

use monstertruck_modeling::{Point3, Solid, Vector3, builder};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let v = builder::vertex(Point3::origin());
    let e = builder::extrude(&v, Vector3::z());
    let f = builder::extrude(&e, Vector3::x());
    let cube0 = builder::extrude(&f, Vector3::y());

    let v = builder::vertex(Point3::new(0.25, 0.25, 0.25));
    let e = builder::extrude(&v, Vector3::x() * 0.5);
    let f = builder::extrude(&e, Vector3::z() * 0.5);
    let cube1 = builder::extrude(&f, Vector3::y() * 0.5);

    let mut boundaries = cube0.into_boundaries();
    boundaries.extend(cube1.into_boundaries());

    let solid = Solid::new_unchecked(boundaries);
    let json = serde_json::to_string_pretty(&solid)?;
    std::fs::write("cube-in-cube.json", json)?;
    Ok(())
}
