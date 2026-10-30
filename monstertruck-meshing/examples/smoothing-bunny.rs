//! Adds smooth normals to the stanford bunny.
//!
//! - Input: bunny.obj
//! - Output: `smooth_bunny.obj`

use monstertruck_mesh::obj;
use monstertruck_meshing::filters::NormalFilters;

fn main() -> anyhow::Result<()> {
    const PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/obj/bunny.obj");
    std::fs::copy(PATH, "bunny.obj")?;
    let file = std::fs::File::open(PATH)?;
    let mut mesh = obj::read(file)?;
    mesh.add_smooth_normals(std::f64::consts::PI / 3.0, true);

    let file = std::fs::File::create("smooth_bunny.obj")?;
    obj::write(&mesh, file)?;
    Ok(())
}
