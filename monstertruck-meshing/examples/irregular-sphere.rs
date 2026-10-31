//! Add the normal to the sphere containing
//! the irregular normal generated from the NURBS containing the critical point.
//!
//! - Input: `irregular_sphere.obj`
//! - Output: `regular_sphere.obj`

use monstertruck_mesh::{KERNEL_EPSILON, obj};
use monstertruck_meshing::filters::{NormalFilters, OptimizingFilter};

fn main() -> anyhow::Result<()> {
    const PATH: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../resources/obj/irregular_sphere.obj"
    );
    std::fs::copy(PATH, "irregular_shpere.obj")?;
    let file = std::fs::File::open(PATH)?;
    let mut mesh = obj::read(file)?;
    mesh.normalize_normals()
        .remove_unused_attrs()
        .put_together_same_attrs(KERNEL_EPSILON * 2.0)
        .remove_degenerate_faces()
        .add_smooth_normals(std::f64::consts::PI / 6.0, true)
        .remove_unused_attrs();
    let file = std::fs::File::create("regular_sphere.obj")?;
    obj::write(&mesh, file)?;
    Ok(())
}
