//! Analyze real OBJ meshes at multiple tessellation levels without depending on
//! known CAD primitive shapes. Emits one JSON record per input path.
//!
//! Run using managed Cargo, for example:
//! `cargo run -p monstertruck-modeling --example mesh_recovery_corpus -- mesh.obj`

use monstertruck_mesh::{PolygonMesh, obj};
use monstertruck_modeling::mesh_recovery::{
    EvidenceOptions, FeatureEdgeKind, MeshEvidence, analyze_mesh_evidence,
};
use serde_json::{Value, json};
use std::collections::HashSet;
use std::error::Error;
use std::fs::{self, File};
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

/// Optional independent labels used for evaluation, never geometry extraction.
struct NativeFaceLabels {
    /// Source CAD face index for each input polygon.
    polygon_face_ids: Vec<usize>,
    /// Original surface family for each CAD face.
    surface_types: Vec<String>,
}

fn support_interface_summary(evidence: &MeshEvidence) -> Value {
    let mut paths = 0;
    let mut loops = 0;
    let mut branches = 0;
    let mut source_edge_count = 0;
    let mut curved_multi_segment_interfaces = 0;
    for interface in &evidence.support_interfaces {
        match interface.topology {
            monstertruck_modeling::mesh_recovery::SupportBoundaryTopology::OpenPath => paths += 1,
            monstertruck_modeling::mesh_recovery::SupportBoundaryTopology::ClosedLoop => loops += 1,
            monstertruck_modeling::mesh_recovery::SupportBoundaryTopology::Branched => {
                branches += 1
            }
        }
        source_edge_count += interface.source_edges.len();
        if interface.source_edges.len() >= 3 && interface.max_seed_direction_spread_radians > 1.0e-4
        {
            curved_multi_segment_interfaces += 1;
        }
    }
    json!({
        "connected_components": evidence.support_interfaces.len(),
        "source_contact_edges": source_edge_count,
        "open_paths": paths,
        "closed_loops": loops,
        "branches": branches,
        "curved_multi_segment_interfaces": curved_multi_segment_interfaces,
        "max_region_contact_degree": evidence.support_interface_indices
            .iter().map(Vec::len).max().unwrap_or(0),
        "regions_degree_2_to_8": evidence.support_interface_indices
            .iter().filter(|ids| (2..=8).contains(&ids.len())).count(),
        "regions_degree_over_16": evidence.support_interface_indices
            .iter().filter(|ids| ids.len() > 16).count(),
        "regions_degree_over_64": evidence.support_interface_indices
            .iter().filter(|ids| ids.len() > 64).count(),
    })
}

/// Evaluate native-face disagreement for source-interface components without
/// assuming a primitive surface or consuming truth during region extraction.
fn native_interface_error_examples(evidence: &MeshEvidence, labels: &[usize]) -> Vec<Value> {
    let selected_edges = evidence
        .edges
        .iter()
        .filter(|edge| edge.kind == FeatureEdgeKind::Dihedral)
        .map(|edge| edge.vertices)
        .collect::<HashSet<_>>();
    let mut examples = Vec::new();
    for interface in &evidence.support_interfaces {
        let [a, b] = interface.supports;
        let region_a = &evidence.planar_supports[a];
        let region_b = &evidence.planar_supports[b];
        if region_a.source_face_indices.len() < 3 && region_b.source_face_indices.len() < 3 {
            continue;
        }
        let mut false_edges = 0;
        let mut true_edges = 0;
        for (&edge, &[left, right]) in interface
            .source_edges
            .iter()
            .zip(&interface.source_polygon_pairs)
        {
            if !selected_edges.contains(&edge) {
                continue;
            }
            if labels[left] == labels[right] {
                false_edges += 1;
            } else {
                true_edges += 1;
            }
        }
        if false_edges == 0 {
            continue;
        }
        examples.push(json!({
            "false_edges": false_edges,
            "true_edges": true_edges,
            "support_sizes": [region_a.source_face_indices.len(), region_b.source_face_indices.len()],
            "support_areas": [region_a.area, region_b.area],
            "boundary_length": interface.boundary_length,
            "turn_mean_degrees": interface.mean_turn_radians.to_degrees(),
            "turn_min_degrees": interface.min_turn_radians.to_degrees(),
            "turn_max_degrees": interface.max_turn_radians.to_degrees(),
            "max_line_offset": interface.max_seed_line_offset,
            "max_tangent_spread_degrees": interface.max_seed_direction_spread_radians.to_degrees(),
            "source_segment_count": interface.source_edges.len(),
            "topology": format!("{:?}", interface.topology),
        }));
    }
    examples.sort_unstable_by(|a, b| {
        b["false_edges"]
            .as_u64()
            .cmp(&a["false_edges"].as_u64())
            .then_with(|| {
                b["source_segment_count"]
                    .as_u64()
                    .cmp(&a["source_segment_count"].as_u64())
            })
    });
    examples.truncate(12);
    examples
}

/// Compare geometric normal-rotation hypotheses with independent native CAD
/// face IDs. The geometry code never sees the labels or these study thresholds.
fn native_normal_rotation_reference(evidence: &MeshEvidence, labels: &[usize]) -> Value {
    let native_per_support = evidence
        .planar_supports
        .iter()
        .map(|support| {
            let first = support
                .source_face_indices
                .first()
                .map(|&face| labels[face])?;
            support
                .source_face_indices
                .iter()
                .all(|&face| labels[face] == first)
                .then_some(first)
        })
        .collect::<Vec<_>>();
    let mut compared = 0_usize;
    let mut all_reference = [0_usize; 3]; // Same source face / other faces / mixed.
    let mut aligned = 0;
    let mut coherent = 0;
    let mut rotating = [0_usize; 3]; // All / same CAD face / different CAD faces.
    let mut small_turn = [0_usize; 3];
    for (center, interfaces) in evidence.support_interface_indices.iter().enumerate() {
        for (i, &first) in interfaces.iter().enumerate() {
            for &second in &interfaces[i + 1..] {
                let Some(measured) = evidence.support_normal_rotation(center, first, second) else {
                    continue;
                };
                compared += 1;
                let triplet = [
                    native_per_support[center],
                    native_per_support[measured.neighbor_supports[0]],
                    native_per_support[measured.neighbor_supports[1]],
                ];
                let valid = triplet.iter().all(Option::is_some);
                let same_source = valid && triplet[0] == triplet[1] && triplet[1] == triplet[2];
                if same_source {
                    all_reference[0] += 1;
                } else if valid {
                    all_reference[1] += 1;
                } else {
                    all_reference[2] += 1;
                }
                if measured.interface_axis_residual_radians > 1_f64.to_radians()
                    || measured.max_interface_segment_spread_radians > 1_f64.to_radians()
                {
                    continue;
                }
                aligned += 1;
                if !measured.coherent_rotation {
                    continue;
                }
                coherent += 1;
                let max_turn = measured.turn_radians[0].max(measured.turn_radians[1]);
                if max_turn < 1.0e-6
                    || measured.turn_magnitude_mismatch_radians > max_turn * 0.2
                    || measured.max_normal_axis_dot > 0.1
                {
                    continue;
                }
                rotating[0] += 1;
                if same_source {
                    rotating[1] += 1;
                } else if valid {
                    rotating[2] += 1;
                }
                if max_turn < 30_f64.to_radians() {
                    small_turn[0] += 1;
                    if same_source {
                        small_turn[1] += 1;
                    } else if valid {
                        small_turn[2] += 1;
                    }
                }
            }
        }
    }
    json!({
        "all_local_pairs": compared,
        "all_native_reference": {
            "within_one_native_face": all_reference[0],
            "crossing_native_faces": all_reference[1],
            "mixed_source_region_labels": all_reference[2],
        },
        "near_parallel_interface_pairs": aligned,
        "aligned_coherent_rotation_pairs": coherent,
        "near_parallel_equal_normal_turn": {
            "all": rotating[0], "within_one_native_face": rotating[1],
            "crossing_native_faces": rotating[2],
            "mixed_source_region_labels": rotating[0] - rotating[1] - rotating[2],
        },
        "same_but_turn_below_30_deg": {
            "all": small_turn[0], "within_one_native_face": small_turn[1],
            "crossing_native_faces": small_turn[2],
            "mixed_source_region_labels": small_turn[0] - small_turn[1] - small_turn[2],
        },
    })
}

/// Load independently supplied CAD face identities; missing sidecars mean
/// the caller is analyzing an ordinary unlabeled mesh.
fn read_native_labels(
    path: &Path,
    polygon_count: usize,
) -> Result<Option<NativeFaceLabels>, String> {
    let labels_path = path.with_extension("labels");
    if !labels_path.exists() {
        return Ok(None);
    }
    let text = fs::read_to_string(&labels_path).map_err(|error| error.to_string())?;
    let labels = text
        .lines()
        .map(str::parse::<usize>)
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| error.to_string())?;
    if labels.len() != polygon_count {
        return Err(format!(
            "expected {polygon_count} native-face labels, got {}",
            labels.len()
        ));
    }
    let surfaces_path = path
        .parent()
        .ok_or("missing source parent directory")?
        .join("surfaces.json");
    let surfaces: Value = serde_json::from_str(
        &fs::read_to_string(&surfaces_path).map_err(|error| error.to_string())?,
    )
    .map_err(|error| error.to_string())?;
    let kinds = surfaces
        .as_array()
        .ok_or("surfaces.json must be an array")?
        .iter()
        .map(|surface| {
            surface
                .get("kind")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .ok_or_else(|| "missing native CAD surface type".to_owned())
        })
        .collect::<Result<Vec<_>, _>>()?;
    if labels.iter().any(|&label| label >= kinds.len()) {
        return Err("out-of-bounds source surface label".to_owned());
    }
    Ok(Some(NativeFaceLabels {
        polygon_face_ids: labels,
        surface_types: kinds,
    }))
}

/// Compare a candidate geometry observation with optional independent native
/// CAD-face identities. Sidecars are solely benchmark references and never
/// influence extraction or construction.
fn evaluate_native_face_labels(
    path: &Path,
    evidence: &MeshEvidence,
) -> Result<Option<Value>, String> {
    let Some(NativeFaceLabels {
        polygon_face_ids: labels,
        surface_types: kinds,
    }) = read_native_labels(path, evidence.source_planar_support.len())?
    else {
        return Ok(None);
    };
    let mut native_edges = 0;
    let mut changed_surface_types = 0;
    let mut supported_candidates = 0;
    let mut candidates_on_native = 0;
    let mut candidates_on_changed_type = 0;
    let angle_thresholds_degrees = [5.0_f64, 15.0, 30.0];
    let mut baseline = [[0_usize; 3]; 3];
    let mut union = [0_usize; 3];
    for (index, edge) in evidence.edges.iter().enumerate() {
        let Some([first_region, second_region]) = evidence.edge_planar_support_pair(index) else {
            continue;
        };
        let [first_source_face, second_source_face] = edge.source_face_indices.as_slice() else {
            continue;
        };
        let first_native_face = labels[*first_source_face];
        let second_native_face = labels[*second_source_face];
        let native_boundary = first_native_face != second_native_face;
        let changed_type = kinds[first_native_face] != kinds[second_native_face];
        native_edges += usize::from(native_boundary);
        changed_surface_types += usize::from(changed_type);
        let backed = first_region != second_region
            && [first_region, second_region]
                .iter()
                .any(|&region| evidence.planar_supports[region].source_face_indices.len() >= 3);
        if backed {
            supported_candidates += 1;
            candidates_on_native += usize::from(native_boundary);
            candidates_on_changed_type += usize::from(changed_type);
        }
        let angle = edge.dihedral_angle_radians.unwrap_or_default();
        for (threshold, counts) in angle_thresholds_degrees.iter().zip(&mut baseline) {
            if angle >= threshold.to_radians() {
                counts[0] += 1;
                counts[1] += usize::from(native_boundary);
                counts[2] += usize::from(changed_type);
            }
        }
        if backed || angle >= 30_f64.to_radians() {
            union[0] += 1;
            union[1] += usize::from(native_boundary);
            union[2] += usize::from(changed_type);
        }
    }
    let baseline = angle_thresholds_degrees
        .iter()
        .zip(baseline)
        .map(|(angle_degrees, counts)| {
            json!({
                "threshold_degrees": angle_degrees,
                "selected": counts[0],
                "true_native": counts[1],
                "true_surface_type": counts[2],
            })
        })
        .collect::<Vec<_>>();
    Ok(Some(json!({
        "scope": "retained two-sided dihedral edges only",
        "source_native_face_boundaries": native_edges,
        "source_surface_type_boundaries": changed_surface_types,
        "candidate_edges": supported_candidates,
        "candidate_true_native": candidates_on_native,
        "candidate_false_native": supported_candidates - candidates_on_native,
        "candidate_true_surface_type": candidates_on_changed_type,
        "angle_baselines": baseline,
        "false_candidate_interface_examples": native_interface_error_examples(evidence, &labels),
        "normal_rotation_reference": native_normal_rotation_reference(evidence, &labels),
        "planar_or_30_degree_union": {
            "selected": union[0],
            "true_native": union[1],
            "true_surface_type": union[2],
        },
    })))
}

/// Count candidate edge support from the same observed planar-region graph
/// used by the labelled native-CAD confusion matrix.
fn planar_candidate_counts(evidence: &MeshEvidence) -> [usize; 3] {
    let mut supported_seam_candidates = 0;
    let mut bilateral_plane_candidates = 0;
    let mut supported_region_internal_edges = 0;
    for index in 0..evidence.edges.len() {
        let Some([a, b]) = evidence.edge_planar_support_pair(index) else {
            continue;
        };
        let supported_a = evidence.planar_supports[a].source_face_indices.len() >= 3;
        let supported_b = evidence.planar_supports[b].source_face_indices.len() >= 3;
        if a == b {
            supported_region_internal_edges += 1;
        } else if supported_a || supported_b {
            supported_seam_candidates += 1;
            if supported_a && supported_b {
                bilateral_plane_candidates += 1;
            }
        }
    }
    [
        supported_seam_candidates,
        bilateral_plane_candidates,
        supported_region_internal_edges,
    ]
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
            let [
                supported_seam_candidates,
                bilateral_plane_candidates,
                supported_region_internal_edges,
            ] = planar_candidate_counts(&evidence);
            let ground_truth = match evaluate_native_face_labels(path, &evidence) {
                Ok(truth) => truth,
                Err(error) => {
                    return json!({"source": base, "stage": "ground_truth",
                        "error": error});
                }
            };
            json!({
                "source": base,
                "stage": "ok",
                "ground_truth": ground_truth,
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
                "support_interfaces": support_interface_summary(&evidence),
                "planar_regions": evidence.planar_supports.len(),
                "planar_regions_multi_facet": evidence
                    .planar_supports.iter().filter(|region| region.source_face_indices.len() >= 3).count(),
                "planar_supported_seam_candidates": supported_seam_candidates,
                "bilateral_planar_seam_candidates": bilateral_plane_candidates,
                "same_planar_region_feature_edges": supported_region_internal_edges,
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
