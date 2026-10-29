//! Modeling a torus by two sweeps.
//!
//! Generated json file can be visualized by `simple-shape-viewer`, an example of `monstertruck-render`.

use monstertruck_modeling::*;

fn modeling(radius0: f64, radius1: f64) -> Result<Solid> {
    let v = builder::vertex(Point3::new(radius0, 0.0, radius1));
    let w = builder::revolve(
        &v,
        Point3::new(radius0, 0.0, 0.0),
        Vector3::unit_y(),
        builder::SweepAngle::Closed,
        2,
    )?;
    let shell = builder::revolve(
        &w,
        Point3::origin(),
        Vector3::unit_z(),
        builder::SweepAngle::Closed,
        2,
    )?;
    Ok(Solid::new_unchecked(vec![shell]))
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args = std::env::args().collect::<Vec<_>>();
    let (radius0, radius1, filename) = match args.len() {
        1 => (0.75, 0.25, "torus.json".to_string()),
        3 => (args[1].parse()?, args[2].parse()?, "torus.json".to_string()),
        4 => (args[1].parse()?, args[2].parse()?, args[3].clone()),
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "expected 0, 2, or 3 arguments",
            )
            .into());
        }
    };
    let torus = modeling(radius0, radius1)?;
    let json = serde_json::to_vec_pretty(&torus)?;
    std::fs::write(filename, json)?;
    Ok(())
}
