//! Analyze real OBJ meshes at multiple tessellation levels without depending on
//! known CAD primitive shapes. Emits one JSON record per input path.
//!
//! Run using managed Cargo, for example:
//! cargo run -p monstertruck-modeling --example mesh_recovery_corpus -- mesh.obj

use monstertruck_mesh::{PolygonMesh, obj};
use monstertruck_modeling::mesh_recovery::{
    EvidenceOptions, FeatureEdgeKind, MeshEvidence, analyze_mesh_evidence,
};
use serde_json::{Value, json};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn best_axes(evidence: &MeshEvidence) -> Vec<Value> {
    let mut indices = (0..evidence.parallel_bundles.len()).collect::<Vec<_>>();
    indices.sort_unstable_by_key(|&id| {
        std::cmp::Reverse(evidence.parallel_bundles[id].placements.len())
    });
    indices
        .into_iter()
        .take(6)
        .map(|id| {
            let group = &evidence.parallel_bundles[id];
            let max_residual = group
                .placements
                .iter()
                .map(|placement| placement.angular_residual_radians)
                .fold(0.0_f64, f64::max);
            let direction = group.direction;
            let span = group
                .placements
                .iter()
                .map(|placement| {
                    placement.longitudinal_interval[1] - placement.longitudinal_interval[0]
                })
                .sum::<f64>();
            json!({
                "axis": [direction.x, direction.y, direction.z],
                "members": group.placements.len(),
                "summed_length": span,
                "max_angle_residual": max_residual,
            })
        })
        .collect()
}

fn analyze_file(path: &Path) -> Value {
    let start = Instant::now();
    let base = json!({
        "file": path.to_string_lossy(),
        "model": path.parent()
            .and_then(|dir| dir.file_name())
            .map(|name| name.to_string_lossy().to_string()),
        "lod": path.file_stem().map(|name| name.to_string_lossy().to_string()),
    });
    let loaded = File::open(path)
        .map(BufReader::new)
        .map_err(|error| error.to_string())
        .and_then(|reader| obj::read(reader).map_err(|error| error.to_string()));
    let mesh: PolygonMesh = match loaded {
        Ok(mesh) => mesh,
        Err(error) => {
            return json!({"source": base, "stage": "obj_read",
                "error": error, "seconds": start.elapsed().as_secs_f64()});
        }
    };
    let loaded_seconds = start.elapsed().as_secs_f64();
    let faces = mesh.face_iter().count();
    let vertices = mesh.positions().len();
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default());
    match evidence {
        Ok(evidence) => {
            let boundaries = evidence
                .edges
                .iter()
                .filter(|edge| edge.kind == FeatureEdgeKind::Boundary)
                .count();
            let dihedrals = evidence
                .edges
                .iter()
                .filter(|edge| edge.kind == FeatureEdgeKind::Dihedral)
                .count();
            let total_segments = evidence
                .chains
                .iter()
                .map(|chain| chain.edge_indices.len())
                .sum::<usize>();
            json!({
                "source": base,
                "stage": "ok",
                "source_faces": faces,
                "source_vertices": vertices,
                "canonical_vertices": evidence.positions.len(),
                "feature_edges": evidence.edges.len(),
                "boundary_edges": boundaries,
                "dihedral_edges": dihedrals,
                "suppressed_coplanar_edges": evidence.suppressed_interior_edge_count,
                "non_manifold_edges": evidence.non_manifold_edge_count,
                "winding_conflicts": evidence.inconsistent_winding_edge_count,
                "chains": evidence.chains.len(),
                "chain_segment_total": total_segments,
                "parallel_groups": evidence.parallel_bundles.len(),
                "top_parallel_axes": best_axes(&evidence),
                "read_seconds": loaded_seconds,
                "evidence_seconds": start.elapsed().as_secs_f64() - loaded_seconds,
                "total_seconds": start.elapsed().as_secs_f64(),
            })
        }
        Err(error) => json!({
            "source": base,
            "stage": "evidence",
            "source_faces": faces,
            "source_vertices": vertices,
            "error": error.to_string(),
            "read_seconds": loaded_seconds,
            "evidence_seconds": start.elapsed().as_secs_f64() - loaded_seconds,
        }),
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let paths = std::env::args_os()
        .skip(1)
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    if paths.is_empty() {
        return Err("expected one or more OBJ paths".into());
    }
    let stdout = io::stdout();
    let mut output = io::BufWriter::new(stdout.lock());
    for path in paths {
        let result = analyze_file(&path);
        writeln!(output, "{result}")?;
        output.flush()?;
    }
    Ok(())
}
