//! Adds smooth normals to and quadrangulate the famous teapot.
//!
//! - Input: teapot.obj
//! - Output: `quaded_pot.obj`

use monstertruck_mesh::{TOLERANCE, obj};
use monstertruck_meshing::filters::{NormalFilters, OptimizingFilter, StructuringFilter};

fn main() -> anyhow::Result<()> {
    const PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/obj/teapot.obj");
    std::fs::copy(PATH, "teapot.obj")?;
    let file = std::fs::File::open(PATH)?;
    let mut mesh = obj::read(file)?;

    mesh.put_together_same_attrs(TOLERANCE * 2.0)
        .add_smooth_normals(std::f64::consts::PI / 3.0, true)
        .quadrangulate(0.1, 1.0);
    let file = std::fs::File::create("quaded_pot.obj")?;
    obj::write(&mesh, file)?;
    Ok(())
}
