//! An experiment to decompose a mesh into elements for future NURBS-shape approximation of the mesh.
//!
//! - Input: sample.obj
//! - Output: `planes_parts_#.obj`, `others_parts_#.obj`

use monstertruck_mesh::obj;
use monstertruck_meshing::{analyzers::Splitting, filters::NormalFilters};

fn main() -> anyhow::Result<()> {
    const PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/obj/sample.obj");
    std::fs::copy(PATH, "sample.obj")?;
    let file = std::fs::File::open(PATH)?;
    let mut mesh = obj::read(file)?;
    mesh.add_smooth_normals(std::f64::consts::PI / 3.0, true);

    let (planes, others) = mesh.extract_planes(0.01);
    let planes = mesh.create_mesh_by_face_indices(&planes);
    let others = mesh.create_mesh_by_face_indices(&others);
    let planes_parts = planes.components(true);
    let others_parts = others.components(true);

    std::fs::create_dir_all("output")?;
    for (i, faces) in planes_parts.into_iter().enumerate() {
        let mesh = planes.create_mesh_by_face_indices(&faces);
        let file = std::fs::File::create(format!("output/planes_parts_{i}.obj"))?;
        obj::write(&mesh, file)?;
    }
    for (i, faces) in others_parts.into_iter().enumerate() {
        let mesh = others.create_mesh_by_face_indices(&faces);
        let file = std::fs::File::create(format!("output/others_parts_{i}.obj"))?;
        obj::write(&mesh, file)?;
    }
    Ok(())
}
