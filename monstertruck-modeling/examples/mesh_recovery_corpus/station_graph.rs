//! Measure direct projected-station recovery independently of strip unions.

use super::{PartitionCarriers, partition_contact_metrics};

use monstertruck_modeling::mesh_recovery::{
    AnalyticRecoveryOptions, CurvatureGrowthHypothesis, CurvatureTrack, MeshEvidence,
    NativeTrimGraph, RecoveredAnalyticMesh, RegionBoundaryIndex, RegionContactGraph,
    SelectedSurfaceProposal, SharedTrimGeometry, SharedTrimGraph, SpatialTrackEvidence,
    StationMeridianCandidate, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionProposal,
    construct_analytic_recovery, discover_analytic_regions, recover_analytic_evidence,
};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

fn recovery_options() -> AnalyticRecoveryOptions {
    AnalyticRecoveryOptions {
        max_planar_displacement: 0.0055,
        ..Default::default()
    }
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
        .select_surface_partition(
            &proposals,
            planar_baseline,
            recovery_options().max_partition_search_nodes,
        )
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
            vertex: recovery_options().relative_surface_tolerance * scale,
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

pub fn summarize(
    evidence: &MeshEvidence,
    _boundary: &RegionBoundaryIndex<'_>,
    spatial: &[Option<SpatialTrackEvidence>],
    grown: &[CurvatureGrowthHypothesis],
    local_tracks: &[CurvatureTrack],
) -> Result<Value, String> {
    let options = recovery_options();
    let discovered = discover_analytic_regions(evidence, spatial, grown, local_tracks, options)
        .map_err(|e| e.to_string())?;
    let recovery = construct_analytic_recovery(discovered, options).map_err(|e| e.to_string())?;
    let active = &recovery.evidence;
    let stats = &recovery.discovery;
    let boundary = RegionBoundaryIndex::new(active).map_err(|e| e.to_string())?;
    let mut used = BTreeSet::new();
    let mut selected = Vec::new();
    for candidate in &stats.circular_candidates {
        if candidate.source_chords.upper_bound > candidate.chord_tolerance
            || candidate.supports.iter().any(|id| used.contains(id))
        {
            continue;
        }
        used.extend(candidate.supports.iter().copied());
        selected.push(candidate);
    }
    let baseline = selected_partition(active, evidence, &boundary, &selected)?;
    let axes=stats.axes.iter().map(|a|json!({
        "source_track":a.source_track,"source_circle_group":a.source_circle_group,
        "axis":[a.axis.x,a.axis.y,a.axis.z],"origin":[a.origin.x,a.origin.y,a.origin.z],
        "stations":a.stations,"repeated_transitions":a.repeated_transitions,
        "seeds":a.seed_count,"profiles":a.profile_count,"candidates":a.circular_candidates,"linear_candidates":a.linear_candidates,
    })).collect::<Vec<_>>();
    let circle_summary=stats.circle_axes.map(|c|json!({"source_triples":c.source_triples,"circles":c.circles,"repeated_axes":c.repeated_axes,"axes_added":c.axes_added}));
    let history = stats
        .refinement_passes
        .iter()
        .map(|p| json!({"support_count":p.support_count,"split_supports":p.split_supports}))
        .collect::<Vec<_>>();
    let mut result = json!({
        "axes":axes,"candidates":stats.circular_candidates.iter().map(|c|candidate_summary(active,c)).collect::<Vec<_>>(),
        "selected":selected.iter().map(|c|candidate_summary(active,c)).collect::<Vec<_>>(),
        "applied_partition":baseline,"joint_partition":joint_summary(&recovery),
        "circle_axis_discovery":circle_summary,
        "support_refinement":{"passes":history,"original_supports":stats.original_supports,"final_supports":active.planar_supports.len(),"placed_axes":stats.placed_axes},
    });
    if evidence.non_manifold_edge_count > 0 {
        let normalized = recover_analytic_evidence(
            evidence,
            AnalyticRecoveryOptions {
                normalize_coincident_sheets: true,
                ..options
            },
        )
        .map_err(|e| e.to_string())?;
        result["normalized_recovery"] = json!({
            "canceled_sheets":normalized.evidence.provenance.canceled_sheets.iter().map(|s|json!({"sides":s.sides,"resolved_edges":s.resolved_edges,"max_plane_residual":s.max_plane_residual})).collect::<Vec<_>>(),
            "original_polygon_ids":normalized.evidence.provenance.original_polygon_ids,
            "joint_partition":joint_summary(&normalized),
            "complete_shell":normalized.complete_shell().is_some(),
        });
    }
    Ok(result)
}

fn joint_summary(recovery: &RecoveredAnalyticMesh) -> Value {
    let evidence = &recovery.evidence;
    let stats = &recovery.discovery;
    let repaired = &recovery.reconciled;
    let partition = &repaired.partition;
    let carriers = PartitionCarriers::new(
        evidence,
        partition.carriers.clone(),
        &partition.region_by_polygon,
    );
    let repair_summary = json!({"accepted":repaired.planes.len(),"rejected":repaired.rejected,"max_displacement":repaired.max_displacement,
        "planes":repaired.planes.iter().map(|p|json!({"region":p.region,"supports":p.supports,"max_plane_residual":p.max_plane_residual,"normal_shift":p.normal_shift,
        "source_polygons":p.supports.iter().flat_map(|&id|&evidence.planar_supports[id].source_face_indices).collect::<Vec<_>>() })).collect::<Vec<_>>()});
    json!({
        "curved_proposals":stats.curved_proposals,"expanded_circular_proposals":stats.expanded_circular_proposals,"linear_proposals":stats.linear_proposals,
        "direct_linear_proposals":stats.direct_linear_proposals,"rim_linear_proposals":stats.rim_linear_proposals,"compared_native_neighborhoods":stats.compared_native_neighborhoods,"expanded_native_neighborhoods":stats.expanded_native_neighborhoods,"direct_linear_hypotheses":stats.linear_candidates.len(),
        "selected_curved":stats.selected_curved,"selected_linear":stats.selected.len()-stats.selected_curved,
        "selected":stats.selected.iter().map(|p|proposal_summary(evidence,p)).collect::<Vec<_>>(),
        "consolidated_supports":partition.consolidated_supports,"region_count":partition.carriers.iter().flatten().count(),
        "original_planar_supports":stats.original_supports,"exact_overlap_components":partition.exact_components,"limited_overlap_components":partition.limited_components,
        "contact_graph":partition_contact_metrics(&partition.contacts,evidence,Some(&carriers)),
        "reconciliation":repair_summary,"native":native_summary(&recovery.native,partition),
        "analytic_trims":trim_summary(&recovery.trims,&partition.contacts,&carriers),
    })
}

fn native_summary(
    native: &NativeTrimGraph,
    partition: &monstertruck_modeling::mesh_recovery::SurfaceRegionPartition,
) -> Value {
    let mut errors = BTreeMap::<String, usize>::new();
    let mut faces = [0_usize; 3];
    let mut built = [0_usize; 3];
    let mut chart = [0_usize; 3];
    let mut failed = Vec::new();
    for region in &native.regions {
        let kind = match partition.carriers[region.region] {
            Some(SurfaceCarrier::Plane { .. }) => 0,
            Some(SurfaceCarrier::LinearRevolution { .. }) => 1,
            _ => 2,
        };
        faces[kind] += 1;
        if let Err(error) = &region.faces {
            *errors.entry(error.to_string()).or_default() += 1;
            if failed.len() < 128 {
                failed.push(json!({"region":region.region,"kind":kind,"error":error.to_string()}));
            }
        } else {
            built[kind] += 1;
            chart[kind] += usize::from(region.requires_chart_cut);
        }
    }
    let mut trim_errors = BTreeMap::<String, usize>::new();
    let mut edges = 0;
    for trim in &native.trims {
        match trim {
            Ok(trim) => edges += trim.len(),
            Err(error) => *trim_errors.entry(error.to_string()).or_default() += 1,
        }
    }
    let shell = native
        .regions
        .iter()
        .filter_map(|region| region.faces.as_ref().ok())
        .flatten()
        .cloned()
        .collect::<monstertruck_modeling::Shell>();
    json!({"faces_by_kind":faces,"native_faces_by_kind":built,"requires_chart_cut_by_kind":chart,
        "face_errors":errors,"failed_faces_first_128":failed,"trim_errors":trim_errors,"native_edges":edges,
        "max_vertex_displacement":native.max_vertex_displacement,"chart_seams":native.regions.iter().map(|r|r.chart_seams).sum::<usize>(),"native_face_count":shell.len(),"shell_condition":format!("{:?}",shell.shell_condition())})
}

fn proposal_summary(evidence: &MeshEvidence, proposal: &SelectedSurfaceProposal) -> Value {
    let polygons = proposal
        .supports
        .iter()
        .flat_map(|&id| {
            evidence.planar_supports[id]
                .source_face_indices
                .iter()
                .copied()
        })
        .collect::<Vec<_>>();
    let carrier = match proposal.carrier {
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
    json!({"carrier":carrier,"source_polygons":polygons,"supports":proposal.supports.len(),"relative_vertex_error":proposal.relative_vertex_error,"relative_chord_bound":proposal.relative_chord_bound,"independently_tested_profile":proposal.independently_tested_profile})
}

fn trim_summary(
    trims: &SharedTrimGraph,
    graph: &RegionContactGraph,
    carriers: &PartitionCarriers,
) -> Value {
    let mut counts = [0_usize; 6];
    let mut edges = [0_usize; 6];
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
            SharedTrimGeometry::Conic(_) => 3,
            SharedTrimGeometry::ToroidalSection(_) => 4,
            SharedTrimGeometry::LinearSection(_) => 5,
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
    json!({"toroidal_section_pieces":counts[4],"toroidal_section_edges":edges[4],"linear_section_pieces":counts[5],"linear_section_edges":edges[5],"conic_pieces":counts[3],"conic_edges":edges[3],"line_pieces":counts[0],"circle_pieces":counts[1],"unresolved_pieces":counts[2],"line_edges":edges[0],"circle_edges":edges[1],"unresolved_edges":edges[2],
        "fully_supported_paths":path_counts.iter().filter(|counts|counts[1]==0).count(),"split_paths":path_counts.iter().filter(|counts|counts[0]>1).count(),"limited_paths":trims.limited_paths,
        "revolution_contact_edges":active_revolution_edges,"unresolved_revolution_edges":unresolved_revolution_edges,
        "circle_examples_first_128":circle_examples,"unresolved_examples_first_128":unresolved_examples,
        "region_loop_count":loop_count,"pinched_regions":pinched_regions,"boundary_errors":boundary_errors,
        "ambiguous_edges":trims.ambiguous_edges.len()})
}
