//! Measure direct projected-station recovery independently of strip unions.

use super::{PartitionCarriers, partition_contact_metrics};

use monstertruck_modeling::mesh_recovery::{
    CircleAxisOptions, CurvatureGrowthHypothesis, CurvatureTrack, MeridianGraphOptions,
    MeridianStationGraph, MeshEvidence, RegionBoundaryIndex, RegionContactGraph,
    SharedTrimGeometry, SharedTrimGraph, SharedTrimOptions, SpatialTrackEvidence,
    StationLinearCandidate, StationMeridianCandidate, SurfaceCarrier, SurfaceRegionBudget,
    SurfaceRegionProposal, discover_circle_axes, measure_curvature_spatial,
};
use monstertruck_modeling::{InnerSpace, Point3, Vector3};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

const RELATIVE_BUDGET: f64 = 1.0e-4;
const AXIS_RESOLUTION: f64 = 1.0e-8;
// Local chord approximation envelope. For a cylindrical chord the
// conservative squared-radius bound is r*sin²(azimuth_step/2); a quarter
// radius corresponds to a 60 degree source step. This is a declared
// sampling-quality constraint, not a claim that every inscribed polygon
// at finer resolution was originally analytic.
const LINEAR_CHORD_FRACTION: f64 = 0.25;
// A quarter-turn circular chord has sagitta r*(1-cos(pi/4)).
// Limit the absolute surface deviation to this fraction of the observed
// meridian scale; exact vertex interpolation alone cannot justify a surface
// with coarser deviation. Continuous bounds are refined at this budget.
const CIRCULAR_CHORD_FRACTION: f64 = 1.0 - std::f64::consts::FRAC_1_SQRT_2;

struct PlacedAxis {
    origin: Point3,
    axis: Vector3,
    source_track: Option<usize>,
    source_circle_group: Option<usize>,
    circular_seed: bool,
}

fn independent_axes(
    evidence: &MeshEvidence,
    spatial: &[Option<SpatialTrackEvidence>],
    local_tracks: &[CurvatureTrack],
) -> Vec<PlacedAxis> {
    let Some(&reference) = evidence.positions.first() else {
        return Vec::new();
    };
    let extent = evidence
        .positions
        .iter()
        .map(|&p| (p - reference).magnitude())
        .fold(0.0_f64, f64::max);
    if extent <= 0.0 {
        return Vec::new();
    }
    let mut tracks = local_tracks
        .iter()
        .enumerate()
        .filter_map(|(id, track)| {
            let fit = spatial
                .get(id)?
                .clone()
                .or_else(|| measure_curvature_spatial(evidence, track))?;
            Some((id, fit))
        })
        .filter(|(_, fit)| {
            fit.max_sample_radial_residual <= AXIS_RESOLUTION * fit.mean_sample_radius
                && fit.max_support_plane_offset_residual <= AXIS_RESOLUTION * fit.mean_sample_radius
        })
        .collect::<Vec<_>>();
    tracks.sort_unstable_by(|(a_id, a), (b_id, b)| {
        a.relative_max_sample_radial_residual
            .total_cmp(&b.relative_max_sample_radial_residual)
            .then_with(|| {
                b.transverse_normal_rank_ratio
                    .total_cmp(&a.transverse_normal_rank_ratio)
            })
            .then(a_id.cmp(b_id))
    });
    let mut result: Vec<PlacedAxis> = Vec::new();
    for (id, fit) in tracks {
        let axis = fit.axis;
        let origin = fit.axis_origin + axis * (reference - fit.axis_origin).dot(axis);
        if result.iter().any(|placed| {
            axis.cross(placed.axis).magnitude() <= AXIS_RESOLUTION
                && (origin - placed.origin).magnitude() <= AXIS_RESOLUTION * extent
        }) {
            continue;
        }
        result.push(PlacedAxis {
            origin,
            axis,
            source_track: Some(id),
            source_circle_group: None,
            circular_seed: true,
        });
    }
    result
}

fn extend_with_circle_axes(evidence: &MeshEvidence, axes: &mut Vec<PlacedAxis>) -> Value {
    let mut circle_summary = Value::Null;
    if let Some(&reference) = evidence.positions.first() {
        let extent = evidence
            .positions
            .iter()
            .map(|&point| (point - reference).magnitude())
            .fold(0.0_f64, f64::max);
        let resolution = extent * AXIS_RESOLUTION;
        if let Some(discovery) = discover_circle_axes(
            evidence,
            CircleAxisOptions {
                vertex_tolerance: resolution,
                relative_tolerance: RELATIVE_BUDGET,
                axis_position_tolerance: 8.0 * resolution,
                axis_angle_tolerance_radians: 1.0e-7,
            },
        ) {
            let mut added = 0;
            for (id, group) in discovery.axes.iter().enumerate() {
                if let Some(placed) = axes.iter_mut().find(|placed| {
                    group.axis.cross(placed.axis).magnitude() <= AXIS_RESOLUTION
                        && (group.origin - placed.origin).magnitude() <= resolution
                }) {
                    placed.circular_seed |= group.stations.len() >= 4;
                    continue;
                }
                axes.push(PlacedAxis {
                    origin: group.origin,
                    axis: group.axis,
                    source_track: None,
                    source_circle_group: Some(id),
                    circular_seed: group.stations.len() >= 4,
                });
                added += 1;
            }
            circle_summary = json!({"source_triples":discovery.source_triples,"circles":discovery.circles.len(),"repeated_axes":discovery.axes.len(),"axes_added":added});
        }
    }
    circle_summary
}

fn selected_partition(
    evidence: &MeshEvidence,
    planar_baseline: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    selected: &[&StationMeridianCandidate],
) -> Result<Value, String> {
    let proposals = selected
        .iter()
        .filter_map(|c| curved_proposal(boundary, c))
        .collect::<Vec<_>>();
    let partition = boundary
        .select_surface_partition(&proposals, planar_baseline, 100_000)
        .map_err(|error| error.to_string())?;
    let carriers =
        PartitionCarriers::new(evidence, partition.carriers, &partition.region_by_polygon);
    Ok(json!({
        "selected_regions":partition.selected_proposals.len(),
        "source_supports_consolidated":partition.consolidated_supports,
        "contact_graph":partition_contact_metrics(&partition.contacts,evidence,Some(&carriers)),
    }))
}

fn curved_proposal<'a>(
    boundary: &RegionBoundaryIndex<'a>,
    candidate: &StationMeridianCandidate,
) -> Option<SurfaceRegionProposal<'a>> {
    let scale = candidate.fit.geometric_scale();
    boundary.measure_surface_region(
        &candidate.supports,
        SurfaceCarrier::from_circular_meridian(&candidate.fit),
        SurfaceRegionBudget {
            vertex: RELATIVE_BUDGET * scale,
            chord: candidate.chord_tolerance,
        },
    )
}

fn candidate_summary(evidence: &MeshEvidence, candidate: &StationMeridianCandidate) -> Value {
    let fit = &candidate.fit;
    let polygons = candidate
        .supports
        .iter()
        .flat_map(|&id| {
            evidence.planar_supports[id]
                .source_face_indices
                .iter()
                .copied()
        })
        .collect::<Vec<_>>();
    json!({
        "supports":candidate.supports.len(),
        "source_polygons":polygons,
        "source_vertices":fit.source_vertices,
        "resolved_meridian_samples":fit.resolved_meridian_samples,
        "repeated_stations":candidate.repeated_stations,
        "axial_center":[fit.axial_center.x,fit.axial_center.y,fit.axial_center.z],
        "axis":[fit.axis.x,fit.axis.y,fit.axis.z],
        "radial_center":fit.radial_center,
        "profile_radius":fit.profile_radius,
        "observed_meridian_span":fit.observed_meridian_span,
        "relative_max_vertex_residual":fit.max_vertex_residual/fit.profile_radius,
        "max_normal_angle_radians":fit.max_normal_angle_radians,
        "boundary_loops":candidate.boundary_loops,
        "source_chord_witness":candidate.source_chords.witnessed_max,
        "source_chord_upper_bound":candidate.source_chords.upper_bound,"chord_tolerance":candidate.chord_tolerance,
        "relative_source_chord_witness":candidate.source_chords.witnessed_max/fit.geometric_scale(),
        "relative_source_chord_upper_bound":candidate.source_chords.upper_bound/fit.geometric_scale(),
    })
}

fn graph_options(evidence: &MeshEvidence, placed: &PlacedAxis) -> MeridianGraphOptions {
    let extent = evidence
        .positions
        .iter()
        .map(|&p| (p - placed.origin).magnitude())
        .fold(0.0_f64, f64::max);
    let resolution = extent * AXIS_RESOLUTION;
    MeridianGraphOptions {
        station_resolution: resolution,
        vertex_tolerance: 4.0 * resolution,
        relative_tolerance: RELATIVE_BUDGET,
        max_linear_relative_chord_error: LINEAR_CHORD_FRACTION,
        max_circular_relative_chord_error: CIRCULAR_CHORD_FRACTION,
        normal_tolerance_radians: 1.0e-5,
    }
}

pub fn summarize(
    evidence: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    spatial: &[Option<SpatialTrackEvidence>],
    grown: &[CurvatureGrowthHypothesis],
    local_tracks: &[CurvatureTrack],
) -> Result<Value, String> {
    let mut axes = independent_axes(evidence, spatial, local_tracks);
    let circle_summary = extend_with_circle_axes(evidence, &mut axes);
    let mut refined = None;
    let mut history = Vec::new();
    loop {
        let active = refined.as_ref().unwrap_or(evidence);
        let split = axes
            .iter()
            .filter_map(|placed| {
                MeridianStationGraph::new(
                    active,
                    placed.origin,
                    placed.axis,
                    graph_options(active, placed),
                )
            })
            .flat_map(|graph| graph.mixed_support_refinements())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        history.push(
            json!({"support_count":active.planar_supports.len(),"split_supports":split.len()}),
        );
        if split.is_empty() {
            break;
        }
        // Monotone refinement, finite because each split makes singletons.
        refined = Some(
            active
                .refine_planar_supports(&split)
                .map_err(|error| error.to_string())?,
        );
    }
    let mut result = if let Some(active) = &refined {
        let active_boundary =
            RegionBoundaryIndex::new(active).map_err(|error| error.to_string())?;
        evaluate_axes(active, evidence, &active_boundary, &axes, grown)?
    } else {
        evaluate_axes(evidence, evidence, boundary, &axes, grown)?
    };
    result["circle_axis_discovery"] = circle_summary;
    result["support_refinement"] = json!({"passes":history,"original_supports":evidence.planar_supports.len(),
        "final_supports":refined.as_ref().map_or(evidence.planar_supports.len(),|active|active.planar_supports.len()),"placed_axes":axes.len()});
    Ok(result)
}

fn evaluate_axes(
    evidence: &MeshEvidence,
    planar_baseline: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    axes: &[PlacedAxis],
    grown: &[CurvatureGrowthHypothesis],
) -> Result<Value, String> {
    let mut results = Vec::new();
    let mut candidates = BTreeMap::<Vec<usize>, StationMeridianCandidate>::new();
    let mut linear = BTreeMap::<Vec<usize>, StationLinearCandidate>::new();
    for placed in axes {
        let options = graph_options(evidence, placed);
        let Some(graph) = MeridianStationGraph::new(evidence, placed.origin, placed.axis, options)
        else {
            continue;
        };
        let linear_candidates = graph.discover_linear(boundary);
        let linear_count = linear_candidates.len();
        for candidate in linear_candidates {
            let entry = linear
                .entry(candidate.supports.clone())
                .or_insert_with(|| candidate.clone());
            if candidate.max_vertex_residual < entry.max_vertex_residual {
                *entry = candidate;
            }
        }
        if !placed.circular_seed {
            if linear_count > 0 {
                results.push(json!({"source_track":placed.source_track,"source_circle_group":placed.source_circle_group,
                    "axis":[placed.axis.x,placed.axis.y,placed.axis.z],"origin":[placed.origin.x,placed.origin.y,placed.origin.z],
                    "stations":graph.stations().len(),"linear_candidates":linear_count,"candidates":0}));
            }
            continue;
        }
        let found = graph.discover_circular(boundary);
        results.push(json!({
            "source_track":placed.source_track,
            "source_circle_group":placed.source_circle_group,
            "axis":[placed.axis.x,placed.axis.y,placed.axis.z],
            "origin":[placed.origin.x,placed.origin.y,placed.origin.z],
            "stations":graph.stations().len(),
            "repeated_transitions":graph.transitions().iter().filter(|transition|transition.repeated).count(),
            "seeds":found.seed_count,
            "profiles":found.profile_count,
            "candidates":found.candidates.len(),
            "linear_candidates":linear_count,
        }));
        for candidate in found.candidates {
            let entry = candidates
                .entry(candidate.supports.clone())
                .or_insert_with(|| candidate.clone());
            if candidate.fit.max_vertex_residual < entry.fit.max_vertex_residual {
                *entry = candidate;
            }
        }
    }
    let mut candidates = candidates.into_values().collect::<Vec<_>>();
    candidates.sort_unstable_by(|a, b| {
        b.source_polygons
            .cmp(&a.source_polygons)
            .then_with(|| {
                a.fit
                    .max_vertex_residual
                    .total_cmp(&b.fit.max_vertex_residual)
            })
            .then_with(|| a.supports.cmp(&b.supports))
    });
    let mut used = BTreeSet::new();
    let mut selected = Vec::new();
    for candidate in &candidates {
        // Provisional admissibility: even a perfectly interpolating surface
        // must remain within the observed feature scale along source chords.
        // Retain failures as alternatives, not as selected analytic regions.
        if candidate.source_chords.upper_bound > candidate.chord_tolerance {
            continue;
        }
        if candidate.supports.iter().any(|id| used.contains(id)) {
            continue;
        }
        used.extend(candidate.supports.iter().copied());
        selected.push(candidate);
    }
    let partition = selected_partition(evidence, planar_baseline, boundary, &selected)?;
    let linear = linear.into_values().collect::<Vec<_>>();
    let joint = joint_partition(
        evidence,
        planar_baseline,
        boundary,
        &candidates,
        &linear,
        grown,
    )?;
    Ok(
        json!({"axes":results,"candidates":candidates.iter().map(|c|candidate_summary(evidence,c)).collect::<Vec<_>>(),"selected":selected.iter().map(|c|candidate_summary(evidence,c)).collect::<Vec<_>>(),"applied_partition":partition,"joint_partition":joint}),
    )
}

fn trim_options(evidence: &MeshEvidence) -> SharedTrimOptions {
    let extent = evidence.positions.first().map_or(0.0, |&reference| {
        evidence
            .positions
            .iter()
            .map(|&point| (point - reference).magnitude())
            .fold(0.0_f64, f64::max)
    });
    SharedTrimOptions {
        absolute_tolerance: extent * AXIS_RESOLUTION * 8.0,
        relative_tolerance: RELATIVE_BUDGET,
        max_interval_edges: 128,
    }
}

fn joint_partition(
    evidence: &MeshEvidence,
    planar_baseline: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    curved: &[StationMeridianCandidate],
    linear: &[StationLinearCandidate],
    grown: &[CurvatureGrowthHypothesis],
) -> Result<Value, String> {
    let mut proposals = Vec::new();
    let mut seen = BTreeSet::new();
    proposals.extend(curved.iter().filter_map(|c| curved_proposal(boundary, c)));
    let curved_count = proposals.len();
    for candidate in linear {
        if let Some(proposal) = boundary.measure_surface_region(
            &candidate.supports,
            candidate.carrier,
            SurfaceRegionBudget {
                vertex: candidate.vertex_tolerance,
                chord: candidate.chord_tolerance,
            },
        ) {
            seen.insert(candidate.supports.clone());
            proposals.push(proposal);
        }
    }
    let direct_linear_count = proposals.len() - curved_count;
    for candidate in grown {
        let fit = &candidate.spatial;
        let scale = fit.mean_sample_radius.min(fit.axial_coordinate_span);
        if scale <= 0.0
            || fit.max_sample_radial_residual > RELATIVE_BUDGET * scale
            || fit.max_support_plane_offset_residual > RELATIVE_BUDGET * fit.mean_sample_radius
            || fit.support_normal_axial_projection_spread > RELATIVE_BUDGET
        {
            continue;
        }
        let supports = candidate
            .support_indices
            .iter()
            .flat_map(|&id| &planar_baseline.planar_supports[id].source_face_indices)
            .map(|&polygon| evidence.source_planar_support[polygon])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        if !seen.insert(supports.clone()) {
            continue;
        }
        let Some(carrier) = SurfaceCarrier::from_spatial_fit(fit) else {
            continue;
        };
        if let Some(proposal) = boundary.measure_surface_region(
            &supports,
            carrier,
            SurfaceRegionBudget {
                vertex: RELATIVE_BUDGET * scale,
                chord: LINEAR_CHORD_FRACTION
                    * fit
                        .axial_stations
                        .iter()
                        .map(|s| s.radius_mean)
                        .fold(f64::INFINITY, f64::min),
            },
        ) {
            proposals.push(proposal);
        }
    }
    let partition = boundary
        .select_surface_partition(&proposals, planar_baseline, 100_000)
        .map_err(|error| error.to_string())?;
    let chosen = partition
        .selected_proposals
        .iter()
        .map(|&id| proposal_summary(evidence, &proposals[id]))
        .collect::<Vec<_>>();
    let selected_curved = partition
        .selected_proposals
        .iter()
        .filter(|&&id| id < curved_count)
        .count();
    let region_count = partition.carriers.iter().flatten().count();
    let trim_graph = boundary
        .assemble_shared_trims(&partition, trim_options(evidence))
        .map_err(|error| error.to_string())?;
    let carriers =
        PartitionCarriers::new(evidence, partition.carriers, &partition.region_by_polygon);
    Ok(json!({
        "curved_proposals":curved_count,
        "linear_proposals":proposals.len()-curved_count,
        "direct_linear_proposals":direct_linear_count,
        "direct_linear_hypotheses":linear.len(),
        "selected_curved":selected_curved,
        "selected_linear":chosen.len()-selected_curved,
        "selected":chosen,
        "consolidated_supports":partition.consolidated_supports,
        "region_count":region_count,
        "original_planar_supports":planar_baseline.planar_supports.len(),
        "exact_overlap_components":partition.exact_components,
        "limited_overlap_components":partition.limited_components,
        "contact_graph":partition_contact_metrics(&partition.contacts,evidence,Some(&carriers)),
        "analytic_trims":trim_summary(&trim_graph, &partition.contacts,&carriers),
    }))
}

fn proposal_summary(evidence: &MeshEvidence, proposal: &SurfaceRegionProposal<'_>) -> Value {
    let polygons = proposal
        .supports()
        .iter()
        .flat_map(|&id| {
            evidence.planar_supports[id]
                .source_face_indices
                .iter()
                .copied()
        })
        .collect::<Vec<_>>();
    let carrier = match proposal.carrier() {
        SurfaceCarrier::Plane { origin, normal } => {
            json!({"kind":"plane","origin":[origin.x,origin.y,origin.z],"normal":[normal.x,normal.y,normal.z]})
        }
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } => {
            json!({"kind":"linear","origin":[origin.x,origin.y,origin.z],"axis":[axis.x,axis.y,axis.z],"radius":radius_at_origin,"slope":slope})
        }
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } => {
            json!({"kind":"circular","origin":[axial_center.x,axial_center.y,axial_center.z],"axis":[axis.x,axis.y,axis.z],"radial_center":radial_center,"profile_radius":profile_radius})
        }
    };
    json!({"carrier":carrier,"source_polygons":polygons,"supports":proposal.supports().len(),"relative_vertex_error":proposal.relative_vertex_error(),"relative_chord_bound":proposal.relative_chord_bound(),"independently_tested_profile":proposal.independently_tested_profile()})
}

fn trim_summary(
    trims: &SharedTrimGraph,
    graph: &RegionContactGraph,
    carriers: &PartitionCarriers,
) -> Value {
    let mut counts = [0_usize; 3];
    let mut edges = [0_usize; 3];
    let mut path_counts = vec![[0_usize; 2]; graph.paths.len()];
    let mut circle_examples = Vec::new();
    let mut unresolved_examples = Vec::new();
    let mut active_revolution_edges = 0;
    let mut unresolved_revolution_edges = 0;
    for trim in &trims.trims {
        let sides = std::iter::once(trim.regions.0)
            .chain(trim.regions.1)
            .filter_map(|id| carriers.carriers[id])
            .collect::<Vec<_>>();
        let nonplanar = sides
            .iter()
            .any(|side| !matches!(side, SurfaceCarrier::Plane { .. }));
        let pairs = || {
            trim.source_edges
                .iter()
                .map(|edge| [Some(edge.interior_polygon), edge.exterior_polygon])
                .collect::<Vec<_>>()
        };
        let kind = match trim.geometry {
            SharedTrimGeometry::Line(_) => 0,
            SharedTrimGeometry::Circle(circle) => {
                if circle_examples.len() < 128 {
                    circle_examples.push(json!({"regions":trim.regions,"source_polygon_pairs":pairs(),
                        "radius":circle.radius,"center":[circle.center.x,circle.center.y,circle.center.z],"normal":[circle.normal.x,circle.normal.y,circle.normal.z],"signed_sweep":circle.signed_sweep_radians}));
                }
                1
            }
            SharedTrimGeometry::Unresolved => {
                if unresolved_examples.len() < 128 {
                    unresolved_examples.push(json!({"regions":trim.regions,"source_polygon_pairs":pairs(),
                        "source_edges":trim.source_edges.iter().map(|edge|edge.vertices).collect::<Vec<_>>(),"revolution":nonplanar}));
                }
                2
            }
        };
        let length = trim.source_edges.len();
        counts[kind] += 1;
        edges[kind] += length;
        path_counts[trim.contact_path][0] += 1;
        path_counts[trim.contact_path][1] += usize::from(kind == 2);
        if nonplanar {
            active_revolution_edges += length;
            if kind == 2 {
                unresolved_revolution_edges += length;
            }
        }
    }
    let loop_count = trims
        .regions
        .iter()
        .filter_map(|region| region.boundary.as_ref().ok())
        .map(|boundary| boundary.loops.len())
        .sum::<usize>();
    let pinched_regions = trims
        .regions
        .iter()
        .filter(|region| {
            region
                .boundary
                .as_ref()
                .is_ok_and(|boundary| !boundary.pinched_vertices.is_empty())
        })
        .count();
    let boundary_errors = trims
        .regions
        .iter()
        .filter_map(|region| {
            region
                .boundary
                .as_ref()
                .err()
                .map(|error| json!({"region":region.region,"error":error.to_string()}))
        })
        .collect::<Vec<_>>();
    json!({"line_pieces":counts[0],"circle_pieces":counts[1],"unresolved_pieces":counts[2],"line_edges":edges[0],"circle_edges":edges[1],"unresolved_edges":edges[2],
        "fully_supported_paths":path_counts.iter().filter(|counts|counts[1]==0).count(),"split_paths":path_counts.iter().filter(|counts|counts[0]>1).count(),"limited_paths":trims.limited_paths,
        "revolution_contact_edges":active_revolution_edges,"unresolved_revolution_edges":unresolved_revolution_edges,
        "circle_examples_first_128":circle_examples,"unresolved_examples_first_128":unresolved_examples,
        "region_loop_count":loop_count,"pinched_regions":pinched_regions,"boundary_errors":boundary_errors,
        "ambiguous_edges":trims.ambiguous_edges.len()})
}
