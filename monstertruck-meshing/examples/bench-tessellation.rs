//! Benchmark tessellation of all shape JSON files.
//!
//! ```bash
//! cargo run -p monstertruck-meshing --example bench-tessellation
//! ```

use anyhow::Context as _;
use monstertruck_mesh::KERNEL_EPSILON;
use monstertruck_meshing::{
    analyzers::Topology,
    filters::OptimizingFilter,
    tessellation::{MeshableShape, MeshedShape},
};
use monstertruck_modeling::{
    Point3,
    geometry::{Curve, Surface},
};
use monstertruck_topology::compress::{CompressedShell, CompressedSolid};
use std::{fs::DirEntry, time::Instant};

type CShell = CompressedShell<Point3, Curve, Surface>;
type CSolid = CompressedSolid<Point3, Curve, Surface>;

fn tessellate_shape(name: &str, json: &str) {
    let instant = Instant::now();
    let mut poly = if let Ok(solid) = serde_json::from_str::<CSolid>(json) {
        solid.triangulation(0.005).to_polygon()
    } else if let Ok(shell) = serde_json::from_str::<CShell>(json) {
        shell.triangulation(0.005).to_polygon()
    } else {
        eprintln!("  {name}: failed to parse");
        return;
    };
    let elapsed = instant.elapsed();
    poly.put_together_same_attrs(KERNEL_EPSILON * 10.0)
        .remove_degenerate_faces()
        .remove_unused_attrs();
    let tri_count = poly.tri_faces().len();
    let condition = poly.shell_condition();
    println!("  {name}: {elapsed:.1?}, {tri_count} triangles, condition={condition:?}");
}

fn main() -> anyhow::Result<()> {
    let shapes_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../resources/shape");
    println!("=== Shape tessellation benchmarks ===");
    let mut entries = std::fs::read_dir(shapes_dir)?.collect::<Result<Vec<_>, _>>()?;
    entries.retain(|entry| entry.path().extension().is_some_and(|ext| ext == "json"));
    entries.sort_by_key(DirEntry::file_name);

    let total = Instant::now();
    for entry in &entries {
        let path = entry.path();
        let name = path
            .file_stem()
            .with_context(|| format!("JSON path has no file stem: {}", path.display()))?
            .to_string_lossy();
        let json = std::fs::read_to_string(&path)
            .with_context(|| format!("failed to read {}", path.display()))?;
        tessellate_shape(&name, &json);
    }
    println!("  total: {:.1?}", total.elapsed());
    Ok(())
}
