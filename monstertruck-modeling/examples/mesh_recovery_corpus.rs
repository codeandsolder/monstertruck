//! Analyze real OBJ meshes at multiple tessellation levels without depending on
//! known CAD primitive shapes. Emits one JSON record per input path.
//!
//! Run using managed Cargo, for example:
//! `cargo run -p monstertruck-modeling --example mesh_recovery_corpus -- mesh.obj`

use monstertruck_mesh::{PolygonMesh, obj};
use monstertruck_modeling::mesh_recovery::{
    CurvatureGrowthIndex, CurvatureGrowthOptions, CurvatureTrack, CurvatureTrackOptions,
    EvidenceOptions, FeatureEdgeKind, MeshEvidence, SpatialTrackEvidence, analyze_mesh_evidence,
    measure_curvature_spatial, measure_curvature_spatial_join, trace_curvature_tracks,
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
    /// Independently recorded exact original STEP axis, when meaningful:
    /// Cylinders, cones and toroids. Source labels are only benchmark truth.
    analytic_axes: Vec<Option<[f64; 3]>>,
    /// Original centers/axis locations and radii, never used for inference.
    analytic_centers: Vec<Option<[f64; 3]>>,
    native_radii: Vec<Option<f64>>,
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
                branches += 1;
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

/// Original source CAD identities touched by a geometry-only track,
/// collected exclusively in the independent native-CAD benchmark path.
fn track_native_face_ids(
    evidence: &MeshEvidence,
    track: &CurvatureTrack,
    labels: &[usize],
) -> std::collections::BTreeSet<usize> {
    track
        .support_indices
        .iter()
        .flat_map(|&support| {
            evidence.planar_supports[support]
                .source_face_indices
                .iter()
                .map(|&polygon| labels[polygon])
        })
        .collect()
}

/// Reference only: compare whole recovered geometric tracks against source
/// native CAD-face labels without ever giving the labels to track discovery.
fn native_curvature_track_reference(
    evidence: &MeshEvidence,
    tracks: &[CurvatureTrack],
    labels: &[usize],
    kinds: &[String],
) -> Value {
    let mut pure = 0;
    let mut mixed = 0;
    let mut native_kind = std::collections::BTreeMap::<String, [usize; 2]>::new();
    let mut long_face_ids = std::collections::BTreeSet::<usize>::new();
    let mut longest_pure = 0;
    let minimum_distinct_supports = [4, 6, 8, 12, 16, 24];
    // All tracks / pure within one native STEP face / crossed native faces.
    let mut length_bins = [[0_usize; 3]; 6];
    for track in tracks {
        let distinct_supports = track.distinct_support_count;
        for (minimum, counts) in minimum_distinct_supports.iter().zip(&mut length_bins) {
            if distinct_supports >= *minimum {
                counts[0] += 1;
            }
        }
        let native = track
            .support_indices
            .iter()
            .try_fold(None, |face, &support| {
                let region = &evidence.planar_supports[support];
                region
                    .source_face_indices
                    .iter()
                    .try_fold(face, |id, &source_face| {
                        let native = labels[source_face];
                        if id.is_none_or(|id| id == native) {
                            Some(Some(native))
                        } else {
                            None
                        }
                    })
            });
        let Some(Some(face)) = native else {
            mixed += 1;
            for (minimum, counts) in minimum_distinct_supports.iter().zip(&mut length_bins) {
                if distinct_supports >= *minimum {
                    counts[2] += 1;
                }
            }
            continue;
        };
        pure += 1;
        for (minimum, counts) in minimum_distinct_supports.iter().zip(&mut length_bins) {
            if distinct_supports >= *minimum {
                counts[1] += 1;
            }
        }
        if distinct_supports >= 8 {
            long_face_ids.insert(face);
        }
        longest_pure = longest_pure.max(distinct_supports);
        let counts = native_kind.entry(kinds[face].clone()).or_default();
        counts[0] += 1;
        counts[1] += distinct_supports;
    }
    let top_track_native = tracks
        .iter()
        .take(8)
        .map(|track| {
            let source_faces = track_native_face_ids(evidence, track, labels);
            let mut types = std::collections::BTreeMap::<String, usize>::new();
            for &face in &source_faces {
                *types.entry(kinds[face].clone()).or_default() += 1;
            }
            json!({"native_face_count":source_faces.len(),"source_type_face_counts":types})
        })
        .collect::<Vec<_>>();
    let length_bins = minimum_distinct_supports
        .iter()
        .zip(length_bins)
        .map(|(&minimum, counts)| {
            json!({
                "min_distinct_supports": minimum,
                "all": counts[0],
                "pure_native_face": counts[1],
                "crossing_or_mixed": counts[2],
            })
        })
        .collect::<Vec<_>>();
    json!({
        "track_count": tracks.len(),
        "purity_by_minimum_length": length_bins,
        "top_track_native": top_track_native,
        "unique_native_faces_in_pure_long_tracks": long_face_ids.len(),
        "pure_native_face": pure,
        "crossing_or_mixed_native_faces": mixed,
        "longest_pure_track_supports": longest_pure,
        "native_kinds": native_kind,
    })
}

fn curvature_track_summary(tracks: &[CurvatureTrack], elapsed: f64) -> Value {
    let long = tracks
        .iter()
        .filter(|track| track.distinct_support_count >= 8)
        .count();
    let closed = tracks.iter().filter(|track| track.closed).count();
    json!({
        "track_count": tracks.len(),
        "at_least_8_supports": long,
        "closed_tracks": closed,
        "discovery_seconds": elapsed,
        "top_tracks": tracks.iter().take(8).map(|track| {
            json!({
                "supports": track.support_indices.len(),
                "distinct_supports": track.distinct_support_count,
                "interfaces": track.interface_indices.len(),
                "closed": track.closed,
                "seed_axis": [track.seed_axis.x, track.seed_axis.y, track.seed_axis.z],
                "max_axis_deviation_degrees": track.max_seed_axis_error_radians.to_degrees(),
                "normal_projection_range": track.normal_axis_projection_spread,
                "normal_axis_projection_mean": track.mean_normal_axis_projection,
                "accumulated_normal_turn_degrees": track.accumulated_normal_turn_radians.to_degrees(),
                "max_step_turn_mismatch_degrees": track.max_adjacent_turn_mismatch_radians.to_degrees(),
            })
        }).collect::<Vec<_>>(),
    })
}

/// Summarize source-only independent spatial evidence on longer tracks.
/// Shorter hypotheses are still retained in the production track graph.
fn spatial_track_summary(
    tracks: &[CurvatureTrack],
    fits: &[Option<SpatialTrackEvidence>],
    elapsed: f64,
) -> Value {
    let candidate_count = tracks
        .iter()
        .filter(|track| track.distinct_support_count >= 8)
        .count();
    let measured = fits.iter().filter(|fit| fit.is_some()).count();
    let three_station = fits
        .iter()
        .flatten()
        .filter(|fit| fit.axial_stations.len() >= 3)
        .count();
    let thresholds = [0.01_f64, 0.02, 0.05, 0.1];
    let bounded = thresholds
        .iter()
        .map(|&threshold| {
            let radial = fits
                .iter()
                .flatten()
                .filter(|fit| fit.relative_max_sample_radial_residual <= threshold)
                .count();
            let radial_and_planes = fits
                .iter()
                .flatten()
                .filter(|fit| {
                    fit.relative_max_sample_radial_residual <= threshold
                        && fit.max_support_plane_offset_residual / fit.mean_sample_radius
                            <= threshold
                })
                .count();
            json!({"maximum_relative_residual":threshold,
            "radial_bounded":radial,
            "radial_and_support_planes_bounded":radial_and_planes})
        })
        .collect::<Vec<_>>();
    let top = fits.iter().take(8).map(|fit| fit.as_ref().map(|fit| {
        json!({
            "transverse_normal_rank_ratio":fit.transverse_normal_rank_ratio,
            "source_vertices":fit.unique_source_vertex_count,
            "axis_origin":[fit.axis_origin.x,fit.axis_origin.y,fit.axis_origin.z],
            "signed_plane_offset":fit.signed_support_plane_offset,
            "max_plane_offset_error":fit.max_support_plane_offset_residual,
            "relative_plane_offset_error":fit.max_support_plane_offset_residual / fit.mean_sample_radius,
            "axial_span":fit.axial_coordinate_span,
            "axial_station_count":fit.axial_stations.len(),
            "axial_station_resolution":fit.axial_station_resolution,
            "maximum_meridian_slope_change":fit.max_meridian_slope_change,
            "mean_vertex_radius":fit.mean_sample_radius,
            "vertex_radius_at_origin":fit.fitted_sample_radius_at_origin,
            "vertex_radius_slope":fit.fitted_sample_radius_slope,
            "normal_implied_radius_slope":fit.normal_implied_radius_slope,
            "max_vertex_radial_error":fit.max_sample_radial_residual,
            "rms_vertex_radial_error":fit.rms_sample_radial_residual,
            "relative_vertex_radial_error":fit.relative_max_sample_radial_residual,
        })
    })).collect::<Vec<_>>();
    json!({
        "candidates_at_least_8_supports":candidate_count,
        "spatially_identifiable":measured,
        "unidentifiable":candidate_count-measured,
        "spatial_fit_with_three_or_more_resolved_axial_stations":three_station,
        "spatial_fit_with_fewer_than_three_axial_stations":measured-three_station,
        "measurement_seconds":elapsed,
        "exploratory_residual_thresholds":bounded,
        "top_tracks":top,
    })
}

/// Grow representative source-connected candidates without original native
/// labels. One seed is one competing observation, not a CAD face hypothesis
/// selected by comparison against the original source STEP surface classes.
fn growth_track_summary(
    evidence: &MeshEvidence,
    tracks: &[CurvatureTrack],
) -> Result<Value, String> {
    let start = Instant::now();
    let index = CurvatureGrowthIndex::new(evidence, tracks, CurvatureGrowthOptions::default())
        .map_err(|error| error.to_string())?;
    let indexed_seconds = start.elapsed().as_secs_f64();
    let mut seed_ids = (0..tracks.len()).collect::<Vec<_>>();
    seed_ids.sort_unstable_by_key(|&id| (std::cmp::Reverse(tracks[id].distinct_support_count), id));
    let mut grown = seed_ids
        .into_iter()
        .filter_map(|id| index.grow(id))
        .take(12)
        .collect::<Vec<_>>();
    let sampled_seed_count = grown.len();
    let expanded_seeds = grown
        .iter()
        .filter(|hypothesis| !hypothesis.steps.is_empty())
        .count();
    let maximum_source_supports = grown
        .iter()
        .map(|hypothesis| hypothesis.support_indices.len())
        .max()
        .unwrap_or(0);
    let maximum_tracks_merged = grown
        .iter()
        .map(|hypothesis| hypothesis.track_indices.len())
        .max()
        .unwrap_or(0);
    let maximum_source_polygons = grown
        .iter()
        .map(|hypothesis| hypothesis.source_polygon_count)
        .max()
        .unwrap_or(0);
    grown.sort_unstable_by_key(|h| std::cmp::Reverse(h.support_indices.len()));
    let top = grown
        .iter()
        .take(6)
        .map(|h| {
            json!({
                "seed_track_index":h.seed_track_index,
                "seed_support_count":tracks[h.seed_track_index].distinct_support_count,
                "source_track_ids":h.track_indices,
                "source_supports":h.support_indices.len(),
                "source_polygons":h.source_polygon_count,
                "source_vertex_count":h.spatial.unique_source_vertex_count,
                "radial_relative_error":h.spatial.relative_max_sample_radial_residual,
                "support_plane_relative_error":
                    h.spatial.max_support_plane_offset_residual / h.spatial.mean_sample_radius,
                "meridian_stations":h.spatial.axial_stations.len(),
                "meridian_max_slope_change":h.spatial.max_meridian_slope_change,
                "maximum_seed_axis_deviation_radians":h.max_seed_axis_deviation_radians,
                "maximum_relative_axis_line_separation":h.max_relative_input_axis_separation,
                "growth_steps":h.steps.iter().map(|s| json!({
                    "track_id":s.track_index,
                    "new_supports":s.new_supports,
                    "cumulative_supports":s.cumulative_supports,
                    "cumulative_radial_relative_error":s.cumulative_relative_radial_residual,
                    "cumulative_support_plane_relative_error":s.cumulative_relative_plane_residual,
                })).collect::<Vec<_>>(),
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "eligible_source_tracks":index.eligible_track_count(),
        "sampled_independent_seed_count":sampled_seed_count,
        "seed_policy":"up to 12 locally identifiable seeds, descending original support count, no original STEP metadata",
        "expanded_seed_count":expanded_seeds,
        "maximum_source_support_count":maximum_source_supports,
        "maximum_tracks_in_one_region":maximum_tracks_merged,
        "maximum_original_polygons_in_one_region":maximum_source_polygons,
        "index_seconds":indexed_seconds,
        "index_and_sample_growth_seconds":start.elapsed().as_secs_f64(),
        "samples":top,
    }))
}

/// Compare a recovered geometric normal-space axis with independently
/// observed native STEP analytic axes, including across distinct CAD faces.
/// Source labels are used only after all tracks have been discovered.
fn native_curvature_axis_reference(
    evidence: &MeshEvidence,
    tracks: &[CurvatureTrack],
    labels: &[usize],
    kinds: &[String],
    native_axes: &[Option<[f64; 3]>],
) -> Value {
    let length_thresholds = [4, 6, 8, 12, 16, 24];
    // Total tracks / native analytic-axis evaluable / source axes consistent
    // / inferred axis agrees with every relevant source analytic axis.
    let mut counts = [[0_usize; 4]; 6];
    let mut by_type = std::collections::BTreeMap::<String, [usize; 2]>::new();
    let mut multi_face_axis_matches = 0;
    let mut source_axis_disagreement = 0;
    let mut max_matching_track = 0;
    let threshold = 2_f64.to_radians();
    for track in tracks {
        let distinct = track.distinct_support_count;
        let native_ids = track_native_face_ids(evidence, track, labels);
        let axes = native_ids
            .iter()
            .map(|&id| native_axes[id])
            .collect::<Option<Vec<_>>>();
        let evaluable = axes.as_ref().is_some_and(|axes| !axes.is_empty());
        let [mut consistent, mut matching] = [false, false];
        if let Some(axes) = &axes
            && let Some(&reference) = axes.first()
        {
            let angle = |a: [f64; 3], b: [f64; 3]| {
                let dot = a.iter().zip(b).map(|(&x, y)| x * y).sum::<f64>().abs();
                dot.clamp(0.0, 1.0).acos()
            };
            consistent = axes.iter().all(|&axis| angle(reference, axis) <= threshold);
            let inferred = [track.seed_axis.x, track.seed_axis.y, track.seed_axis.z];
            matching = axes.iter().all(|&axis| angle(inferred, axis) <= threshold);
            if matching && native_ids.len() > 1 {
                multi_face_axis_matches += 1;
            }
            if !consistent {
                source_axis_disagreement += 1;
            }
            if matching {
                max_matching_track = max_matching_track.max(distinct);
            }
            let native_kinds = native_ids
                .iter()
                .map(|&id| kinds[id].as_str())
                .collect::<std::collections::BTreeSet<_>>();
            if native_kinds.len() == 1
                && let Some(&kind) = native_kinds.first()
            {
                let entry = by_type.entry(kind.to_owned()).or_default();
                entry[0] += 1;
                entry[1] += usize::from(matching);
            }
        }
        for (&minimum, bin) in length_thresholds.iter().zip(&mut counts) {
            if distinct >= minimum {
                bin[0] += 1;
                bin[1] += usize::from(evaluable);
                bin[2] += usize::from(consistent);
                bin[3] += usize::from(matching);
            }
        }
    }
    json!({
        "axis_error_threshold_degrees": 2.0,
        "by_minimum_distinct_supports": length_thresholds
            .iter()
            .zip(counts)
            .map(|(&minimum, bin)| json!({
                "minimum": minimum,
                "total_tracks": bin[0],
                "analytic_axes_evaluable": bin[1],
                "original_axes_mutually_consistent": bin[2],
                "inferred_axis_matches_all_original_axes": bin[3],
            }))
            .collect::<Vec<_>>(),
        "by_single_source_surface_family": by_type,
        "matching_axis_across_multiple_native_faces": multi_face_axis_matches,
        "inconsistent_source_analytic_axes": source_axis_disagreement,
        "longest_track_with_original_axis_match": max_matching_track,
    })
}

/// Decode original STEP analytic-vector reference metadata, used only by
/// the native-CAD benchmark evaluator.
fn native_vector3(surface: &Value, field: &str) -> Result<Option<[f64; 3]>, String> {
    let Some(value) = surface.get(field) else {
        return Ok(None);
    };
    let numbers = value
        .as_array()
        .ok_or_else(|| format!("native {field} must be an array"))?;
    let [x, y, z] = numbers.as_slice() else {
        return Err(format!("native {field} must contain exactly three numbers"));
    };
    let coordinates = [x, y, z].map(Value::as_f64);
    let [Some(x), Some(y), Some(z)] = coordinates else {
        return Err(format!("native {field} must contain numbers"));
    };
    if ![x, y, z].into_iter().all(f64::is_finite) {
        return Err(format!("native {field} must contain finite numbers"));
    }
    Ok(Some([x, y, z]))
}

/// Check whether a geometry-only candidate supported by *original cylindrical
/// STEP faces* independently matches the actual native CAD cylinder axes and
/// radii. Original source face truth is never available in production.
fn native_cylinder_spatial_error(
    fit: &SpatialTrackEvidence,
    native_face_ids: &std::collections::BTreeSet<usize>,
    centers: &[Option<[f64; 3]>],
    axes: &[Option<[f64; 3]>],
    radii: &[Option<f64>],
) -> Option<(f64, f64)> {
    let mut max_axis_line_error = 0.0_f64;
    let mut max_radius_error = 0.0_f64;
    for &native in native_face_ids {
        let [cx, cy, cz] = centers[native]?;
        let [ax, ay, az] = axes[native]?;
        let radius = radii[native]?;
        let dx = fit.axis_origin.x - cx;
        let dy = fit.axis_origin.y - cy;
        let dz = fit.axis_origin.z - cz;
        // A point-to-original-axis distance: invariant to the arbitrary
        // longitudinal location of the fitted axis origin.
        let cross_x = dz.mul_add(-ay, dy * az);
        let cross_y = dx.mul_add(-az, dz * ax);
        let cross_z = dy.mul_add(-ax, dx * ay);
        let distance = cross_x.hypot(cross_y).hypot(cross_z);
        max_axis_line_error = max_axis_line_error.max(distance / radius);
        max_radius_error = max_radius_error.max((fit.mean_sample_radius - radius).abs() / radius);
    }
    Some((max_axis_line_error, max_radius_error))
}

/// Compare local-geometry-only joint fits on rare held-out mixed-source
/// candidate tracks which are spatially plausible over a narrow source strip.
/// Native CAD labels select *diagnostic cases*, never influence geometry fits.
fn diagnose_mixed_surface_track_expansions(
    evidence: &MeshEvidence,
    tracks: &[CurvatureTrack],
    fits: &[Option<SpatialTrackEvidence>],
    truth: &NativeFaceLabels,
) -> Value {
    let mut cases = Vec::new();
    for (id, (track, maybe_fit)) in tracks.iter().zip(fits).enumerate() {
        let Some(fit) = maybe_fit else {
            continue;
        };
        let native = track_native_face_ids(evidence, track, &truth.polygon_face_ids);
        if native_track_surface_family(&native, &truth.surface_types)
            != "MixedNativeSurfaceFamilies"
            || fit.relative_max_sample_radial_residual > 0.01
            || fit.max_support_plane_offset_residual / fit.mean_sample_radius > 0.01
        {
            continue;
        }
        let mut source_connected = 0;
        let mut genuinely_larger = 0;
        let mut joint_within_one_percent = 0;
        let mut smallest_joint_radial_error = f64::INFINITY;
        let mut smallest_joint_plane_error = f64::INFINITY;
        let mut longest_union_supports = track.distinct_support_count;
        for (other_id, neighbor) in tracks.iter().enumerate() {
            if id == other_id || neighbor.distinct_support_count < 8 {
                continue;
            }
            let Some(joint) = measure_curvature_spatial_join(evidence, track, neighbor) else {
                continue;
            };
            source_connected += 1;
            if joint.union_support_count <= track.distinct_support_count {
                continue;
            }
            genuinely_larger += 1;
            longest_union_supports = longest_union_supports.max(joint.union_support_count);
            let relative_plane = joint.combined.max_support_plane_offset_residual
                / joint.combined.mean_sample_radius;
            let relative_radial = joint.combined.relative_max_sample_radial_residual;
            smallest_joint_radial_error = smallest_joint_radial_error.min(relative_radial);
            smallest_joint_plane_error = smallest_joint_plane_error.min(relative_plane);
            joint_within_one_percent +=
                usize::from(relative_radial <= 0.01 && relative_plane <= 0.01);
        }
        cases.push(json!({
            "source_original_faces":native.iter().copied().collect::<Vec<_>>(),
            "source_distinct_supports":track.distinct_support_count,
            "source_vertex_count":fit.unique_source_vertex_count,
            "source_axial_stations":fit.axial_stations.len(),
            "source_relative_radial_error":fit.relative_max_sample_radial_residual,
            "source_relative_plane_error":
                fit.max_support_plane_offset_residual / fit.mean_sample_radius,
            "adjacent_candidate_tracks":source_connected,
            "adjacent_joins_adding_original_supports":genuinely_larger,
            "joint_extensions_with_both_errors_below_one_percent":joint_within_one_percent,
            "smallest_extended_radial_error":smallest_joint_radial_error.is_finite()
                .then_some(smallest_joint_radial_error),
            "smallest_extended_plane_error":smallest_joint_plane_error.is_finite()
                .then_some(smallest_joint_plane_error),
            "largest_extended_original_support_count":longest_union_supports,
        }));
    }
    json!(cases)
}

/// CAD surface-family provenance of an entire original source-support track.
/// Identical original families can still represent different actual faces.
fn native_track_surface_family(
    native_ids: &std::collections::BTreeSet<usize>,
    families: &[String],
) -> String {
    let Some(&first) = native_ids.first() else {
        return "EmptyNativeSource".to_owned();
    };
    if native_ids.iter().all(|&id| families[id] == families[first]) {
        families[first].clone()
    } else {
        "MixedNativeSurfaceFamilies".to_owned()
    }
}

/// Inspect the rare cross-family candidate whose *local* spatial errors
/// are small, so we can test its source-connected neighboring geometry.
fn mixed_spatial_counterexample(
    track: &CurvatureTrack,
    fit: &SpatialTrackEvidence,
    native: &std::collections::BTreeSet<usize>,
    kinds: &[String],
) -> Value {
    json!({
        "distinct_supports": track.distinct_support_count,
        "original_native_face_count": native.len(),
        "original_native_faces": native.iter().copied().collect::<Vec<_>>(),
        "original_surface_families": native.iter().map(|&id| &kinds[id])
            .collect::<std::collections::BTreeSet<_>>(),
        "mean_radius":fit.mean_sample_radius,
        "axial_span":fit.axial_coordinate_span,
        "axial_stations":fit.axial_stations.len(),
        "seed_axis":[fit.axis.x,fit.axis.y,fit.axis.z],
        "axis_origin":[fit.axis_origin.x,fit.axis_origin.y,fit.axis_origin.z],
        "max_plane_offset_error_relative":
            fit.max_support_plane_offset_residual / fit.mean_sample_radius,
        "max_source_radius_error_relative":fit.relative_max_sample_radial_residual,
        "max_meridian_slope_change":fit.max_meridian_slope_change,
        "source_polygon_vertices":fit.unique_source_vertex_count,
    })
}

/// Record resolved meridional slope-change thresholds per original native
/// source family. Three distinct axial stations are necessary for this cue.
fn count_native_meridian_slope_bins(
    histogram: &mut std::collections::BTreeMap<String, [usize; 7]>,
    family: &str,
    change: Option<f64>,
    thresholds: &[f64; 7],
) {
    if let Some(change) = change {
        let counts = histogram.entry(family.to_owned()).or_default();
        for (&threshold, total) in thresholds.iter().zip(counts) {
            *total += usize::from(change <= threshold);
        }
    }
}

/// Format source-native family summaries separately from counting fitted
/// evidence so the benchmark keeps numeric comparisons independently testable.
fn native_family_metric_report(
    groups: std::collections::BTreeMap<String, [usize; 4]>,
    slope_histogram: &std::collections::BTreeMap<String, [usize; 7]>,
    thresholds: &[f64; 7],
) -> std::collections::BTreeMap<String, Value> {
    groups
        .into_iter()
        .map(|(kind, counts)| {
            let bounded = slope_histogram.get(&kind).copied().unwrap_or_default();
            let curvature_bins = thresholds
                .iter()
                .zip(bounded)
                .map(|(&threshold, count)| {
                    json!({
                        "maximum_meridian_slope_change":threshold,
                        "count_at_or_below":count,
                    })
                })
                .collect::<Vec<_>>();
            (
                kind,
                json!({
                    "measured":counts[0],
                    "radial_within_1pct":counts[1],
                    "radial_and_support_planes_within_1pct":counts[2],
                    "at_least_three_axial_stations":counts[3],
                    "meridian_slope_change_histogram":curvature_bins,
                }),
            )
        })
        .collect()
}

/// Native-CAD-only reference checks and deliberately exploratory filtering
/// statistics for independently fitted spatial tracks.
fn native_spatial_track_reference(
    evidence: &MeshEvidence,
    tracks: &[CurvatureTrack],
    fits: &[Option<SpatialTrackEvidence>],
    native_labels: &NativeFaceLabels,
) -> Value {
    let mut groups = std::collections::BTreeMap::<String, [usize; 4]>::new();
    let slope_change_thresholds = [1.0e-6_f64, 0.001, 0.01, 0.05, 0.1, 0.5, 1.0];
    let mut slope_change_bounded = std::collections::BTreeMap::<String, [usize; 7]>::new();
    let mut cylinder_evaluable = 0;
    let mut cylinder_axis_placed = 0;
    let mut cylinder_radius_matched = 0;
    let mut cylinder_jointly_matched = 0;
    let mut original_cylinder_cross_face = 0;
    let mut cylinder_error_examples = Vec::<Value>::new();
    let mut passing_mixed_surface_examples = Vec::<Value>::new();
    for (track, evidence_fit) in tracks.iter().zip(fits) {
        let Some(fit) = evidence_fit else {
            continue;
        };
        let native = track_native_face_ids(evidence, track, &native_labels.polygon_face_ids);
        let kind = native_track_surface_family(&native, &native_labels.surface_types);
        let radial_consistent = fit.relative_max_sample_radial_residual <= 0.01;
        let plane_consistent =
            fit.max_support_plane_offset_residual / fit.mean_sample_radius <= 0.01;
        let row = groups.entry(kind.clone()).or_default();
        row[0] += 1;
        row[1] += usize::from(radial_consistent);
        row[2] += usize::from(radial_consistent && plane_consistent);
        row[3] += usize::from(fit.axial_stations.len() >= 3);
        if kind == "MixedNativeSurfaceFamilies" && radial_consistent && plane_consistent {
            passing_mixed_surface_examples.push(mixed_spatial_counterexample(
                track,
                fit,
                &native,
                &native_labels.surface_types,
            ));
        }
        count_native_meridian_slope_bins(
            &mut slope_change_bounded,
            &kind,
            fit.max_meridian_slope_change,
            &slope_change_thresholds,
        );
        if kind != "Cylinder" {
            continue;
        }
        let Some((axis_error, radius_error)) = native_cylinder_spatial_error(
            fit,
            &native,
            &native_labels.analytic_centers,
            &native_labels.analytic_axes,
            &native_labels.native_radii,
        ) else {
            continue;
        };
        cylinder_evaluable += 1;
        let axis_good = axis_error <= 0.01;
        let radius_good = radius_error <= 0.01;
        cylinder_axis_placed += usize::from(axis_good);
        cylinder_radius_matched += usize::from(radius_good);
        cylinder_jointly_matched += usize::from(axis_good && radius_good);
        original_cylinder_cross_face += usize::from(native.len() > 1 && axis_good && radius_good);
        if cylinder_error_examples.len() < 24 {
            cylinder_error_examples.push(json!({
                "distinct_supports":track.distinct_support_count,
                "original_native_cylinder_faces":native.len(),
                "original_axis_offset_relative_to_radius":axis_error,
                "original_radius_error_relative":radius_error,
                "source_vertex_radial_error_relative":fit.relative_max_sample_radial_residual,
                "plane_offset_error_relative":
                    fit.max_support_plane_offset_residual / fit.mean_sample_radius,
            }));
        }
    }
    let family_report =
        native_family_metric_report(groups, &slope_change_bounded, &slope_change_thresholds);
    json!({
        "minimum_distinct_supports_measured":8,
        "exploratory_relative_residual_threshold":0.01,
        "source_family_track_counts_and_residuals":family_report,
        "cylinder_original_axis_radius_evaluable":cylinder_evaluable,
        "cylinder_original_axis_placed_within_1pct_radius":cylinder_axis_placed,
        "cylinder_original_radius_matched_within_1pct":cylinder_radius_matched,
        "cylinder_joint_axis_and_radius_matches":cylinder_jointly_matched,
        "cylinder_joint_matches_across_multiple_native_faces":original_cylinder_cross_face,
        "sample_source_cylinder_fit_errors":cylinder_error_examples,
        "mixed_surface_tracks_with_small_spatial_error":passing_mixed_surface_examples,
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
    let native_surfaces = surfaces
        .as_array()
        .ok_or("surfaces.json must be an array")?;
    let kinds = native_surfaces
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
    let analytic_axes = native_surfaces
        .iter()
        .map(|surface| -> Result<Option<[f64; 3]>, String> {
            let Some(axis) = native_vector3(surface, "axis")? else {
                return Ok(None);
            };
            let norm_sq = axis.iter().map(|x| x * x).sum::<f64>();
            if !norm_sq.is_finite() || norm_sq <= 0.0 {
                return Err("native axis must have finite nonzero length".to_owned());
            }
            Ok(Some(axis.map(|x| x / norm_sq.sqrt())))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let analytic_centers = native_surfaces
        .iter()
        .map(|surface| native_vector3(surface, "center"))
        .collect::<Result<Vec<_>, _>>()?;
    let native_radii = native_surfaces
        .iter()
        .map(|surface| -> Result<Option<f64>, String> {
            let Some(value) = surface.get("radius") else {
                return Ok(None);
            };
            let radius = value
                .as_f64()
                .ok_or_else(|| "native CAD radius must be a number".to_owned())?;
            if !radius.is_finite() || radius <= 0.0 {
                return Err("native CAD radius must be positive and finite".to_owned());
            }
            Ok(Some(radius))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(NativeFaceLabels {
        polygon_face_ids: labels,
        surface_types: kinds,
        analytic_axes,
        analytic_centers,
        native_radii,
    }))
}

/// Compare a candidate geometry observation with optional independent native
/// CAD-face identities. Sidecars are solely benchmark references and never
/// influence extraction or construction.
fn evaluate_native_face_labels(
    path: &Path,
    evidence: &MeshEvidence,
    tracks: &[CurvatureTrack],
    spatial: &[Option<SpatialTrackEvidence>],
) -> Result<Option<Value>, String> {
    let Some(native) = read_native_labels(path, evidence.source_planar_support.len())? else {
        return Ok(None);
    };
    let spatial_reference = native_spatial_track_reference(evidence, tracks, spatial, &native);
    let suspicious_joint_expansions =
        diagnose_mixed_surface_track_expansions(evidence, tracks, spatial, &native);
    let NativeFaceLabels {
        polygon_face_ids: labels,
        surface_types: kinds,
        analytic_axes,
        analytic_centers: _,
        native_radii: _,
    } = native;
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
        "curvature_track_reference": native_curvature_track_reference(evidence, tracks, &labels, &kinds),
        "curvature_axis_reference": native_curvature_axis_reference(
            evidence, tracks, &labels, &kinds, &analytic_axes,
        ),
        "spatial_track_reference": spatial_reference,
        "suspicious_joint_expansions": suspicious_joint_expansions,
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

/// Primitive-independent evidence plus optional original CAD benchmark truth.
struct CorpusCurvatureAnalysis {
    tracks: Value,
    spatial: Value,
    growth: Value,
    native_ground_truth: Option<Value>,
}

/// Keep source-independent curvature-track discovery distinct from the
/// optional original CAD-face reference evaluation.
fn discover_and_evaluate_curvature(
    path: &Path,
    evidence: &MeshEvidence,
) -> Result<CorpusCurvatureAnalysis, (&'static str, String)> {
    let started = Instant::now();
    let tracks = trace_curvature_tracks(evidence, CurvatureTrackOptions::default())
        .map_err(|error| ("curvature", error.to_string()))?;
    let summary = curvature_track_summary(&tracks, started.elapsed().as_secs_f64());
    let spatial_start = Instant::now();
    // Explicitly study long candidates first rather than paying the source
    // vertex fitting cost for every local 4-support observation.
    let spatial = tracks
        .iter()
        .map(|track| {
            (track.distinct_support_count >= 8)
                .then(|| measure_curvature_spatial(evidence, track))
                .flatten()
        })
        .collect::<Vec<_>>();
    let spatial_summary =
        spatial_track_summary(&tracks, &spatial, spatial_start.elapsed().as_secs_f64());
    let growth = growth_track_summary(evidence, &tracks).map_err(|error| ("growth", error))?;
    let ground_truth = evaluate_native_face_labels(path, evidence, &tracks, &spatial)
        .map_err(|error| ("ground_truth", error))?;
    Ok(CorpusCurvatureAnalysis {
        tracks: summary,
        spatial: spatial_summary,
        growth,
        native_ground_truth: ground_truth,
    })
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
            let analysis = match discover_and_evaluate_curvature(path, &evidence) {
                Ok(result) => result,
                Err((stage, error)) => {
                    return json!({"source": base, "stage": stage, "error": error});
                }
            };
            json!({
                "source": base,
                "stage": "ok",
                "ground_truth": analysis.native_ground_truth,
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
                "curvature_tracks": analysis.tracks,
                "spatial_tracks": analysis.spatial,
                "spatial_growth": analysis.growth,
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
