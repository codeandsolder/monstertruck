//! An experiment to decompose a mesh into elements for future NURBS-shape approximation of the mesh.
//! Contains `doc(hidden)` methods.
//!
//! - Input: `filleted_cube.obj`
//! - Output: planes.obj, lower.obj, upper.obj

use monstertruck_mesh::{KERNEL_EPSILON, obj};
use monstertruck_meshing::{
    analyzers::{ExperimentalSplitters, Splitting},
    filters::{NormalFilters, OptimizingFilter, StructuringFilter},
};

fn main() -> anyhow::Result<()> {
    const PATH: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../resources/obj/filleted_cube.obj"
    );
    std::fs::copy(PATH, "filleted_cube.obj")?;
    let file = std::fs::File::open(PATH)?;
    let mut mesh = obj::read(file)?;
    mesh.put_together_same_attrs(KERNEL_EPSILON * 2.0)
        .quadrangulate(0.1, 1.0)
        .add_smooth_normals(std::f64::consts::PI / 3.0, true);

    let (planes, others) = mesh.extract_planes(0.01);
    let file = std::fs::File::create("planes.obj")?;
    obj::write(&mesh.create_mesh_by_face_indices(&planes), file)?;

    let mesh = mesh.create_mesh_by_face_indices(&others);
    let (upper, lower) = mesh.clustering_faces_by_gcurvature(0.1, false);
    let file = std::fs::File::create("lower.obj")?;
    obj::write(&mesh.create_mesh_by_face_indices(&lower), file)?;
    let file = std::fs::File::create("upper.obj")?;
    obj::write(&mesh.create_mesh_by_face_indices(&upper), file)?;
    Ok(())
}
