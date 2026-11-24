//! Measure direct projected-station recovery independently of strip unions.

use super::{PartitionCarriers, partition_contact_metrics};

use monstertruck_modeling::mesh_recovery::{
    CircleAxisOptions, MeridianGraphOptions, MeridianStationGraph, MeshEvidence,
    RegionBoundaryIndex, SpatialTrackEvidence, StationMeridianCandidate, SurfaceCarrier,
    discover_circle_axes,
};
use monstertruck_modeling::{InnerSpace, Point3, Vector3};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

const RELATIVE_BUDGET: f64 = 1.0e-4;
const AXIS_RESOLUTION: f64 = 1.0e-8;

struct PlacedAxis {
    origin: Point3,
    axis: Vector3,
    source_track: Option<usize>,
    source_circle_group: Option<usize>,
}

fn independent_axes(
    evidence: &MeshEvidence,
    spatial: &[Option<SpatialTrackEvidence>],
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
    let mut tracks = spatial
        .iter()
        .enumerate()
        .filter_map(|(id, fit)| Some((id, fit.as_ref()?)))
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
                if group.stations.len() < 4 {
                    continue;
                }
                if axes.iter().any(|placed| {
                    group.axis.cross(placed.axis).magnitude() <= AXIS_RESOLUTION
                        && (group.origin - placed.origin).magnitude() <= resolution
                }) {
                    continue;
                }
                axes.push(PlacedAxis {
                    origin: group.origin,
                    axis: group.axis,
                    source_track: None,
                    source_circle_group: Some(id),
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
    boundary: &RegionBoundaryIndex<'_>,
    selected: &[&StationMeridianCandidate],
) -> Result<Value, String> {
    let mut region_by_support = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
    let mut carriers = evidence
        .planar_supports
        .iter()
        .map(|support| Some(SurfaceCarrier::from_planar_support(support)))
        .collect::<Vec<_>>();
    for candidate in selected {
        let Some(&region) = candidate.supports.first() else {
            continue;
        };
        for &support in &candidate.supports {
            region_by_support[support] = region;
        }
        carriers[region] = Some(SurfaceCarrier::from_circular_meridian(&candidate.fit));
    }
    let region_by_polygon = evidence
        .source_planar_support
        .iter()
        .map(|&support| region_by_support[support])
        .collect::<Vec<_>>();
    let graph = boundary
        .partition_contacts(&region_by_polygon)
        .map_err(|error| error.to_string())?;
    let carriers = PartitionCarriers::new(evidence, carriers, &region_by_polygon);
    Ok(json!({
        "selected_regions": selected.len(),
        "source_supports_consolidated": selected.iter().map(|candidate|candidate.supports.len()-1).sum::<usize>(),
        "contact_graph": partition_contact_metrics(&graph,evidence,Some(&carriers)),
    }))
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
        "source_chord_upper_bound":candidate.source_chords.upper_bound,
        "relative_source_chord_witness":candidate.source_chords.witnessed_max/fit.geometric_scale(),
        "relative_source_chord_upper_bound":candidate.source_chords.upper_bound/fit.geometric_scale(),
    })
}

pub fn summarize(
    evidence: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    spatial: &[Option<SpatialTrackEvidence>],
) -> Result<Value, String> {
    let mut axes = independent_axes(evidence, spatial);
    let circle_summary = extend_with_circle_axes(evidence, &mut axes);
    let mut results = Vec::new();
    let mut candidates = BTreeMap::<Vec<usize>, StationMeridianCandidate>::new();
    for placed in &axes {
        let extent = evidence
            .positions
            .iter()
            .map(|&p| (p - placed.origin).magnitude())
            .fold(0.0_f64, f64::max);
        let resolution = extent * AXIS_RESOLUTION;
        let options = MeridianGraphOptions {
            station_resolution: resolution,
            vertex_tolerance: 4.0 * resolution,
            relative_tolerance: RELATIVE_BUDGET,
            normal_tolerance_radians: 1.0e-5,
        };
        let Some(graph) = MeridianStationGraph::new(evidence, placed.origin, placed.axis, options)
        else {
            continue;
        };
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
        if candidate.source_chords.upper_bound > candidate.fit.geometric_scale() {
            continue;
        }
        if candidate.supports.iter().any(|id| used.contains(id)) {
            continue;
        }
        used.extend(candidate.supports.iter().copied());
        selected.push(candidate);
    }
    let partition = selected_partition(evidence, boundary, &selected)?;
    Ok(
        json!({"circle_axis_discovery":circle_summary,"axes":results,"candidates":candidates.iter().map(|c|candidate_summary(evidence,c)).collect::<Vec<_>>(),"selected":selected.iter().map(|c|candidate_summary(evidence,c)).collect::<Vec<_>>(),"applied_partition":partition}),
    )
}
