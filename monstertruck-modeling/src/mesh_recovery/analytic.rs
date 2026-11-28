//! Shared analytic region recognition and native reconstruction.
//! Source topology and approximation budgets are distinct. Partial results
//! retain every failed region; only complete oriented shells are returned
//! by the explicit shell extraction method.

use super::evidence::RimLinearOptions;
use super::{
    CircleAxisOptions, CurvatureGrowthHypothesis, CurvatureGrowthIndex, CurvatureGrowthOptions,
    CurvatureTrack, CurvatureTrackOptions, EvidenceOptions, MeridianGraphOptions,
    MeridianStationGraph, MeshEvidence, MeshRecoveryError, NativeTrimError, NativeTrimGraph,
    NativeTrimOptions, ReconciledPartition, ReconciliationOptions, RegionBoundaryError,
    RegionBoundaryIndex, SharedTrimGraph, SharedTrimOptions, SourceCircleObservation,
    SpatialTrackEvidence, StationLinearCandidate, StationMeridianCandidate, SurfaceCarrier,
    SurfaceRegionBudget, SurfaceRegionPartition, SurfaceRegionProposal, analyze_mesh_evidence,
    discover_circle_axes, measure_curvature_spatial, trace_curvature_tracks,
};
use crate::{InnerSpace, KERNEL_EPSILON, Point3, Shell, Vector3};
use monstertruck_mesh::PolygonMesh;
use monstertruck_topology::shell::ShellCondition;
use std::collections::BTreeSet;
mod coplanar_contacts;
pub use coplanar_contacts::CoplanarSupportReassignment;
mod joint_recovery;
pub use super::evidence::{JointFitMethod, JointFitRejection};
pub use joint_recovery::{JointFitAttempt, JointFitReport};
mod carrier_alternatives;
mod native_selection;
use super::carrier_identity::{GeometryBuckets, append_distinct_proposals};
use carrier_alternatives::CarrierAlternatives;
use native_selection::{PartitionAlternative, select_native};

/// Shared policy for evidence, surface recognition and native assembly.
#[derive(Clone, Copy, Debug)]
pub struct AnalyticRecoveryOptions {
    /// Source incidence and planar-support discovery.
    pub evidence: EvidenceOptions,
    /// Local curvature track hypotheses.
    pub tracks: CurvatureTrackOptions,
    /// Spatial track growth hypotheses.
    pub growth: CurvatureGrowthOptions,
    /// Surface fit error relative to the observed feature scale.
    pub relative_surface_tolerance: f64,
    /// Axis placement and station resolution relative to source extent.
    pub axis_relative_resolution: f64,
    /// Angular resolution for independently placed axes.
    pub axis_angle_tolerance_radians: f64,
    /// Normal consistency tolerance for projected station profiles.
    pub profile_normal_tolerance_radians: f64,
    /// Maximum continuous source chord error relative to local radius.
    pub max_linear_relative_chord_error: f64,
    /// Maximum continuous source chord error relative to circular profile scale.
    pub max_circular_relative_chord_error: f64,
    /// Search-node budget per overlapping region proposal component.
    pub max_partition_search_nodes: usize,
    /// Maximum source intervals tested by mixed-trim segmentation.
    pub max_trim_interval_edges: usize,
    /// Maximum additional source-supported displacement in joint reconstruction.
    /// Zero disables this local recovery stage. Not a native incidence tolerance.
    pub max_reconstruction_adjustment: f64,
    /// Known positional uncertainty of the input, when supplied by the caller.
    /// Zero uses local source evidence alone.
    pub source_position_uncertainty: f64,
    /// Maximum failed-contact neighborhood transactions to consider.
    pub max_joint_fit_trials: usize,
    /// Maximum alternative-carrier native assembly trials per starting partition.
    pub max_carrier_trials: usize,
    /// Maximum source displacement for planar reconciliation, in source
    /// coordinate units. Zero disables geometric planar reconciliation.
    pub max_planar_displacement: f64,
    /// Maximum oriented source normal angle of a reconciled plane.
    pub max_planar_normal_angle: f64,
    /// Explicitly allow cancellation of coincident reverse source sheets.
    /// This cannot recover the original profile of a collapsed feature.
    pub normalize_coincident_sheets: bool,
}

impl Default for AnalyticRecoveryOptions {
    fn default() -> Self {
        Self {
            evidence: EvidenceOptions::default(),
            tracks: CurvatureTrackOptions::default(),
            growth: CurvatureGrowthOptions::default(),
            relative_surface_tolerance: 1.0e-4,
            axis_relative_resolution: 1.0e-8,
            axis_angle_tolerance_radians: 1.0e-7,
            profile_normal_tolerance_radians: 1.0e-5,
            // A linear revolution chord at 60 degrees is bounded by r/4.
            max_linear_relative_chord_error: 0.25,
            // A quarter-turn circular chord has this relative sagitta.
            max_circular_relative_chord_error: 1.0 - std::f64::consts::FRAC_1_SQRT_2,
            max_partition_search_nodes: 100_000,
            max_trim_interval_edges: 128,
            max_reconstruction_adjustment: 0.01,
            source_position_uncertainty: 0.0,
            max_joint_fit_trials: 8,
            max_carrier_trials: 8,
            max_planar_displacement: 0.0,
            max_planar_normal_angle: 0.01,
            normalize_coincident_sheets: false,
        }
    }
}

impl AnalyticRecoveryOptions {
    fn validate(self) -> Result<(), AnalyticRecoveryError> {
        let positive = [
            self.relative_surface_tolerance,
            self.axis_relative_resolution,
            self.axis_angle_tolerance_radians,
            self.profile_normal_tolerance_radians,
            self.max_linear_relative_chord_error,
            self.max_circular_relative_chord_error,
        ];
        if positive
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0 || *v >= 1.0)
            || !self.max_reconstruction_adjustment.is_finite()
            || self.max_reconstruction_adjustment < 0.0
            || !self.source_position_uncertainty.is_finite()
            || self.source_position_uncertainty < 0.0
            || !self.max_planar_displacement.is_finite()
            || self.max_planar_displacement < 0.0
            || !self.max_planar_normal_angle.is_finite()
            || self.max_planar_normal_angle <= self.axis_relative_resolution
            || self.max_planar_normal_angle >= std::f64::consts::FRAC_PI_2
            || self.max_partition_search_nodes == 0
            || self.max_trim_interval_edges == 0
        {
            return Err(AnalyticRecoveryError::InvalidOptions);
        }
        Ok(())
    }

    /// Shared trim budget derived from the same source extent and relative
    /// accuracy as surface recognition. Kernel intersection accuracy is
    /// separately limited by the kernel's numerical resolution.
    #[must_use = "the returned options define source trim acceptance"]
    pub fn trim_options(self, evidence: &MeshEvidence) -> SharedTrimOptions {
        let extent = evidence.positions.first().map_or(0.0, |&reference| {
            evidence
                .positions
                .iter()
                .map(|&p| (p - reference).magnitude())
                .fold(0.0_f64, f64::max)
        });
        SharedTrimOptions {
            absolute_tolerance: extent * self.axis_relative_resolution * 8.0,
            relative_tolerance: self.relative_surface_tolerance,
            max_interval_edges: self.max_trim_interval_edges,
        }
    }
}

/// Failure of the recovery workflow before local native region diagnostics.
#[derive(Debug, thiserror::Error)]
pub enum AnalyticRecoveryError {
    /// Invalid numerical options or inconsistent supplied discovery evidence.
    #[error("invalid analytic recovery options or source references")]
    InvalidOptions,
    /// Invalid source geometry.
    #[error(transparent)]
    Source(#[from] MeshRecoveryError),
    /// Inconsistent source partition or boundary construction.
    #[error(transparent)]
    Boundary(#[from] RegionBoundaryError),
    /// Invalid native construction inputs.
    #[error(transparent)]
    Native(#[from] NativeTrimError),
}

/// Counts for independent source circle-axis observations.
#[derive(Clone, Copy, Debug)]
pub struct CircleAxisSummary {
    /// Source triples considered.
    pub source_triples: usize,
    /// Independently sampled circles.
    pub circles: usize,
    /// Groups of circles sharing a placed axis.
    pub repeated_axes: usize,
    /// Distinct axes added to normal-track proposals.
    pub axes_added: usize,
}

/// A placed axis and the candidate evidence discovered on it.
#[derive(Clone, Copy, Debug)]
pub struct AxisDiscoverySummary {
    /// Local curvature track supplying this axis, if any.
    pub source_track: Option<usize>,
    /// Independent source-circle group supplying it, if any.
    pub source_circle_group: Option<usize>,
    /// Unit placed axis direction.
    pub axis: Vector3,
    /// Point on the axis.
    pub origin: Point3,
    /// Projected source stations.
    pub stations: usize,
    /// Repeated station transitions for circular discovery.
    pub repeated_transitions: usize,
    /// Circular profile seeds.
    pub seed_count: usize,
    /// Distinct circular profiles tested.
    pub profile_count: usize,
    /// Circular surface candidates.
    pub circular_candidates: usize,
    /// Linear surface candidates.
    pub linear_candidates: usize,
}

/// One monotone source-support refinement pass.
#[derive(Clone, Copy, Debug)]
pub struct SupportRefinementPass {
    /// Number of supports before splitting.
    pub support_count: usize,
    /// Supports split into original polygons.
    pub split_supports: usize,
}

/// Evidence from a selected independently measured surface proposal.
#[derive(Clone, Debug)]
pub struct SelectedSurfaceProposal {
    /// Sorted original active support IDs.
    pub supports: Vec<usize>,
    /// Selected analytic surface.
    pub carrier: SurfaceCarrier,
    /// Measured vertex residual relative to its bound.
    pub relative_vertex_error: f64,
    /// Continuous chord residual relative to its bound.
    pub relative_chord_bound: f64,
    /// Whether the profile is tested beyond minimal interpolation data.
    pub independently_tested_profile: bool,
}

impl SelectedSurfaceProposal {
    fn from_measured(proposal: &SurfaceRegionProposal<'_>) -> Self {
        Self {
            supports: proposal.supports().to_vec(),
            carrier: proposal.carrier(),
            relative_vertex_error: proposal.relative_vertex_error(),
            relative_chord_bound: proposal.relative_chord_bound(),
            independently_tested_profile: proposal.independently_tested_profile(),
        }
    }
}

/// Candidate discovery and joint selection diagnostics without source labels.
#[derive(Debug)]
pub struct AnalyticDiscoveryStats {
    /// Independent source-circle axis discovery counts.
    pub circle_axes: Option<CircleAxisSummary>,
    /// All distinct placed axes before support refinement.
    pub placed_axes: usize,
    /// Original planar support count.
    pub original_supports: usize,
    /// Monotone support split history.
    pub refinement_passes: Vec<SupportRefinementPass>,
    /// Candidate discovery observations for each contributing axis.
    pub axes: Vec<AxisDiscoverySummary>,
    /// Circular station-region hypotheses.
    pub circular_candidates: Vec<StationMeridianCandidate>,
    /// Direct linear station-region hypotheses.
    pub linear_candidates: Vec<StationLinearCandidate>,
    /// Admissible circular region proposals.
    pub curved_proposals: usize,
    /// Fixed independently measured circular cores expanded over irregular supports.
    pub expanded_circular_proposals: usize,
    /// Fixed station-linear cores expanded across irregular lateral facets.
    pub expanded_linear_proposals: usize,
    /// Admissible linear region proposals (direct and grown).
    pub linear_proposals: usize,
    /// Direct station-graph linear proposals.
    pub direct_linear_proposals: usize,
    /// Additional local circle-rim linear proposals.
    pub rim_linear_proposals: usize,
    /// Measured local circle-rim hypotheses, including alternatives whose
    /// incident trims are not yet reconstructable.
    pub rim_linear_candidates: Vec<SelectedSurfaceProposal>,
    /// Neighborhoods compared using actual native assembly.
    pub compared_native_neighborhoods: usize,
    /// Neighborhoods using the expanded analytic candidate set.
    pub expanded_native_neighborhoods: usize,
    /// Native transactions testing geometrically distinct carriers on identical supports.
    pub carrier_alternative_trials: usize,
    /// Alternative carrier transactions that repaired native geometry.
    pub carrier_alternative_acceptances: usize,
    /// Number of selected circular proposals.
    pub selected_curved: usize,
    /// Selected mutually nonoverlapping measured proposals.
    pub selected: Vec<SelectedSurfaceProposal>,
}

/// Surface partition and its refined source evidence, before native trims.
#[derive(Debug)]
pub struct AnalyticRegionDiscovery {
    /// Active source evidence, including refined planar support ownership.
    pub evidence: MeshEvidence,
    /// Selected region partition.
    pub partition: SurfaceRegionPartition,
    /// Recognition evidence and diagnostics.
    pub stats: AnalyticDiscoveryStats,
    core_alternative: Option<PartitionAlternative>,
    carrier_alternatives: CarrierAlternatives,
}

struct CandidateDiscovery {
    axes: Vec<AxisDiscoverySummary>,
    circular: Vec<StationMeridianCandidate>,
    linear: Vec<StationLinearCandidate>,
}

struct ProposalSelection {
    partition: SurfaceRegionPartition,
    selected: Vec<SelectedSurfaceProposal>,
    curved_proposals: usize,
    expanded_circular_proposals: usize,
    expanded_linear_proposals: usize,
    direct_linear_proposals: usize,
    rim_linear_proposals: usize,
    rim_linear_candidates: Vec<SelectedSurfaceProposal>,
    core_alternative: Option<PartitionAlternative>,
    carrier_alternatives: CarrierAlternatives,
    linear_proposals: usize,
    selected_curved: usize,
}

/// Complete recovery workflow result, including individually failed regions.
#[derive(Debug)]
pub struct RecoveredAnalyticMesh {
    /// Refined active source evidence. Original position IDs retained.
    pub evidence: MeshEvidence,
    /// Evidence from surface recognition and partition selection.
    pub discovery: AnalyticDiscoveryStats,
    /// Final partition and separate reconciled position table.
    pub reconciled: ReconciledPartition,
    /// Shared source-backed analytic trim graph.
    pub trims: SharedTrimGraph,
    /// Kernel-native shared edges and chart faces, with local failures.
    pub native: NativeTrimGraph,
    /// Exact plane ownership corrections applied after native assembly.
    pub coplanar_reassignments: Vec<CoplanarSupportReassignment>,
    /// Accepted uncertainty-bounded joint analytic recovery transactions.
    pub joint_fits: Vec<JointFitReport>,
    /// Rejected local reconstruction attempts with actionable diagnostics.
    pub joint_fit_rejections: Vec<JointFitAttempt>,
    pub(in crate::mesh_recovery) ownership_alternative: Option<NativeAssembly>,
}

impl RecoveredAnalyticMesh {
    /// Extract a complete, oriented, geometrically consistent native shell.
    /// Disconnected components may remain; closure, positive volume and
    /// absence of self-intersection are separate requirements for a solid.
    #[must_use = "the returned shell exists only when the native reconstruction is complete"]
    pub fn complete_shell(&self) -> Option<Shell> {
        native_complete_shell(&self.trims, &self.native)
    }
}

fn native_complete_shell(trims: &SharedTrimGraph, native: &NativeTrimGraph) -> Option<Shell> {
    if !trims.ambiguous_edges.is_empty() || native.trims.iter().any(Result::is_err) {
        return None;
    }
    let mut shell = Shell::new();
    for region in &native.regions {
        if region.requires_chart_cut {
            return None;
        }
        shell.extend(region.faces.as_ref().ok()?.iter().cloned());
    }
    (!shell.is_empty()
        && matches!(
            shell.shell_condition(),
            ShellCondition::Closed | ShellCondition::Oriented
        )
        && shell.is_geometric_consistent())
    .then_some(shell)
}

/// Discover a region partition from already measured track and growth
/// hypotheses. Corpus diagnostics can reuse their source-independent
/// measurements rather than rerunning recognition.
/// # Errors
/// Invalid options, source references or partition construction.
pub fn discover_analytic_regions(
    evidence: &MeshEvidence,
    spatial: &[Option<SpatialTrackEvidence>],
    grown: &[CurvatureGrowthHypothesis],
    tracks: &[CurvatureTrack],
    options: AnalyticRecoveryOptions,
) -> Result<AnalyticRegionDiscovery, AnalyticRecoveryError> {
    options.validate()?;
    if evidence.source_planar_support.is_empty()
        || spatial.len() != tracks.len()
        || grown.iter().any(|g| {
            g.support_indices
                .iter()
                .any(|&id| id >= evidence.planar_supports.len())
        })
    {
        return Err(AnalyticRecoveryError::InvalidOptions);
    }
    let mut axes = independent_axes(evidence, spatial, tracks, options);
    let (circle_axes, circles) = extend_with_circle_axes(evidence, &mut axes, options);
    let (active, refinement_passes, candidates) = refine_and_discover(evidence, &axes, options)?;
    let boundary = RegionBoundaryIndex::new(&active)?;
    let selection = choose_partition(
        &active,
        evidence,
        &boundary,
        &candidates,
        grown,
        &circles,
        options,
    )?;
    let stats = AnalyticDiscoveryStats {
        circle_axes,
        placed_axes: axes.len(),
        original_supports: evidence.planar_supports.len(),
        refinement_passes,
        axes: candidates.axes,
        circular_candidates: candidates.circular,
        linear_candidates: candidates.linear,
        curved_proposals: selection.curved_proposals,
        expanded_circular_proposals: selection.expanded_circular_proposals,
        expanded_linear_proposals: selection.expanded_linear_proposals,
        linear_proposals: selection.linear_proposals,
        direct_linear_proposals: selection.direct_linear_proposals,
        rim_linear_proposals: selection.rim_linear_proposals,
        rim_linear_candidates: selection.rim_linear_candidates,
        compared_native_neighborhoods: 0,
        expanded_native_neighborhoods: 0,
        carrier_alternative_trials: 0,
        carrier_alternative_acceptances: 0,
        selected_curved: selection.selected_curved,
        selected: selection.selected,
    };
    Ok(AnalyticRegionDiscovery {
        evidence: active,
        partition: selection.partition,
        stats,
        core_alternative: selection.core_alternative,
        carrier_alternatives: selection.carrier_alternatives,
    })
}

fn refine_and_discover(
    evidence: &MeshEvidence,
    axes: &[PlacedAxis],
    options: AnalyticRecoveryOptions,
) -> Result<(MeshEvidence, Vec<SupportRefinementPass>, CandidateDiscovery), AnalyticRecoveryError> {
    let mut active = evidence.clone();
    let mut passes = Vec::new();
    // Refinement changes support membership, never source vertex positions.
    // Axes lacking two sampled stations remain impossible on every later pass.
    let mut unsampled = vec![None; axes.len()];
    loop {
        let boundary = RegionBoundaryIndex::new(&active)?;
        let mut split = BTreeSet::new();
        let mut candidates = AxisCandidates::default();
        let coincidence_tolerance = options
            .trim_options(&active)
            .absolute_tolerance
            .min(KERNEL_EPSILON * 0.1);
        for (placed, cached) in axes.iter().zip(&mut unsampled) {
            if let Some(count) = *cached {
                if split.is_empty() {
                    candidates.record_unsampled(placed, count);
                }
                continue;
            }
            let Some((count, graph)) = MeridianStationGraph::for_discovery(
                &active,
                placed.origin,
                placed.axis,
                graph_options(&active, placed, options),
            ) else {
                continue;
            };
            let Some(graph) = graph else {
                *cached = Some(count);
                if split.is_empty() {
                    candidates.record_unsampled(placed, count);
                }
                continue;
            };
            split.extend(graph.mixed_support_refinements());
            // On the final (unchanged) pass, the exact same graph also
            // supplies discovery. Once any refinement is requested,
            // earlier proposals refer to obsolete support IDs and all
            // subsequent discovery is skipped until the next pass.
            if split.is_empty() {
                candidates.discover(placed, &graph, &boundary, &active, coincidence_tolerance);
            }
        }
        passes.push(SupportRefinementPass {
            support_count: active.planar_supports.len(),
            split_supports: split.len(),
        });
        if split.is_empty() {
            return Ok((active, passes, candidates.finish()));
        }
        // Monotone, finite refinement: each split makes original singletons.
        active = active.refine_planar_supports(&split.into_iter().collect::<Vec<_>>())?;
    }
}

/// Construct native shared curves and faces from a discovered partition.
/// # Errors
/// Invalid options, inconsistent partition or native construction inputs.
pub fn construct_analytic_recovery(
    discovered: AnalyticRegionDiscovery,
    options: AnalyticRecoveryOptions,
) -> Result<RecoveredAnalyticMesh, AnalyticRecoveryError> {
    options.validate()?;
    let AnalyticRegionDiscovery {
        evidence,
        partition,
        mut stats,
        core_alternative,
        carrier_alternatives,
    } = discovered;
    let mut ownership = None;
    let mut ownership_rejection = None;
    let mut ownership_alternative = None;
    let mut assembly = if let Some(core) = core_alternative {
        let expanded = PartitionAlternative {
            partition,
            selected: std::mem::take(&mut stats.selected),
        };
        let selection = select_native(&evidence, core, expanded, &carrier_alternatives, options)?;
        stats.carrier_alternative_trials = selection.carrier_trials;
        stats.carrier_alternative_acceptances = selection.carrier_acceptances;
        stats.compared_native_neighborhoods = selection.compared_neighborhoods;
        stats.expanded_native_neighborhoods = selection.expanded_neighborhoods;
        stats.selected = selection.alternative.selected;
        stats.selected_curved = stats
            .selected
            .iter()
            .filter(|proposal| {
                matches!(proposal.carrier, SurfaceCarrier::CircularRevolution { .. })
            })
            .count();
        ownership = selection.ownership;
        ownership_rejection = selection.ownership_rejection;
        ownership_alternative = selection.pending_ownership;
        selection.assembly
    } else {
        let mut alternative = PartitionAlternative {
            partition,
            selected: std::mem::take(&mut stats.selected),
        };
        let mut assembly = assemble_native(&evidence, alternative.partition.clone(), options)?;
        let result =
            carrier_alternatives.repair(&evidence, &mut alternative, &mut assembly, options)?;
        stats.carrier_alternative_trials = result.trials;
        stats.carrier_alternative_acceptances = result.accepted;
        stats.selected = alternative.selected;
        stats.selected_curved = stats
            .selected
            .iter()
            .filter(|p| matches!(p.carrier, SurfaceCarrier::CircularRevolution { .. }))
            .count();
        assembly
    };
    let coplanar_reassignments = coplanar_contacts::repair(&evidence, &mut assembly, options)?;
    let (mut joint_fits, mut joint_fit_rejections) = joint_recovery::repair(
        &evidence,
        &mut assembly,
        options,
        ownership.as_ref().map(|fit| &fit.frozen),
    )?;
    if let Some(fit) = ownership {
        joint_fits.insert(0, fit.report);
    }
    if let Some(rejection) = ownership_rejection {
        joint_fit_rejections.insert(0, rejection);
    }
    let trim = options.trim_options(&evidence);
    RegionBoundaryIndex::new(&evidence)?.repair_planar_crossings(
        &assembly.reconciled.partition,
        &assembly.reconciled.positions,
        &mut assembly.trims,
        &mut assembly.native,
        NativeTrimOptions {
            max_vertex_displacement: trim.absolute_tolerance,
            intersection_tolerance: KERNEL_EPSILON * 0.1,
            angular_resolution: options.axis_relative_resolution,
        },
    )?;
    Ok(RecoveredAnalyticMesh {
        evidence,
        discovery: stats,
        reconciled: assembly.reconciled,
        trims: assembly.trims,
        native: assembly.native,
        coplanar_reassignments,
        joint_fits,
        joint_fit_rejections,
        ownership_alternative,
    })
}

#[derive(Debug)]
pub(in crate::mesh_recovery) struct NativeAssembly {
    pub reconciled: ReconciledPartition,
    pub trims: SharedTrimGraph,
    pub native: NativeTrimGraph,
}

fn assemble_native(
    evidence: &MeshEvidence,
    partition: SurfaceRegionPartition,
    options: AnalyticRecoveryOptions,
) -> Result<NativeAssembly, AnalyticRecoveryError> {
    let boundary = RegionBoundaryIndex::new(evidence)?;
    let trim = options.trim_options(evidence);
    let reconciled = if options.max_planar_displacement > 0.0 {
        boundary.reconcile_planes(
            partition,
            trim,
            ReconciliationOptions {
                max_displacement: options.max_planar_displacement,
                max_normal_angle: options.max_planar_normal_angle,
                parallel_angle_tolerance: options.axis_relative_resolution,
                intersection_tolerance: (trim.absolute_tolerance * 0.1).min(KERNEL_EPSILON * 0.01),
            },
        )?
    } else {
        ReconciledPartition {
            partition,
            positions: evidence.positions.clone(),
            planes: Vec::new(),
            max_displacement: 0.0,
            rejected: 0,
        }
    };
    let trims =
        boundary.assemble_shared_trims(&reconciled.partition, &reconciled.positions, trim)?;
    let native = boundary.construct_native_trims(
        &reconciled.partition,
        &trims,
        &reconciled.positions,
        NativeTrimOptions {
            max_vertex_displacement: trim.absolute_tolerance,
            intersection_tolerance: KERNEL_EPSILON * 0.1,
            angular_resolution: options.axis_relative_resolution,
        },
    )?;
    Ok(NativeAssembly {
        reconciled,
        trims,
        native,
    })
}

/// Recover analytic regions and shared native topology from a polygon mesh.
/// # Errors
/// Invalid source geometry, options or inconsistent partition data.
pub fn recover_analytic_mesh(
    mesh: &PolygonMesh,
    options: AnalyticRecoveryOptions,
) -> Result<RecoveredAnalyticMesh, AnalyticRecoveryError> {
    options.validate()?;
    let evidence = analyze_mesh_evidence(mesh, options.evidence)?;
    recover_analytic_evidence(&evidence, options)
}

/// Reuse existing source evidence for the complete analytic recovery path.
/// # Errors
/// Invalid options or source references.
pub fn recover_analytic_evidence(
    evidence: &MeshEvidence,
    options: AnalyticRecoveryOptions,
) -> Result<RecoveredAnalyticMesh, AnalyticRecoveryError> {
    options.validate()?;
    let normalization = if options.normalize_coincident_sheets {
        evidence.normalize_coincident_sheets(options.evidence)?
    } else {
        None
    };
    let active = normalization.as_ref().unwrap_or(evidence);
    let tracks = trace_curvature_tracks(active, options.tracks)?;
    let spatial = tracks
        .iter()
        .map(|t| {
            (t.distinct_support_count >= 8)
                .then(|| measure_curvature_spatial(active, t))
                .flatten()
        })
        .collect::<Vec<_>>();
    let growth = CurvatureGrowthIndex::new(active, &tracks, options.growth)?;
    let mut grown = (0..tracks.len())
        .filter_map(|id| growth.grow(id))
        .collect::<Vec<_>>();
    grown.sort_unstable_by(|a, b| {
        b.support_indices
            .len()
            .cmp(&a.support_indices.len())
            .then_with(|| {
                a.spatial
                    .relative_max_sample_radial_residual
                    .total_cmp(&b.spatial.relative_max_sample_radial_residual)
            })
            .then(a.seed_track_index.cmp(&b.seed_track_index))
    });
    let discovered = discover_analytic_regions(active, &spatial, &grown, &tracks, options)?;
    construct_analytic_recovery(discovered, options)
}

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
    options: AnalyticRecoveryOptions,
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
            fit.max_sample_radial_residual
                <= options.axis_relative_resolution * fit.mean_sample_radius
                && fit.max_support_plane_offset_residual
                    <= options.axis_relative_resolution * fit.mean_sample_radius
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
            axis.cross(placed.axis).magnitude() <= options.axis_relative_resolution
                && (origin - placed.origin).magnitude() <= options.axis_relative_resolution * extent
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

fn extend_with_circle_axes(
    evidence: &MeshEvidence,
    axes: &mut Vec<PlacedAxis>,
    options: AnalyticRecoveryOptions,
) -> (Option<CircleAxisSummary>, Vec<SourceCircleObservation>) {
    let mut circle_summary = None;
    let mut circles = Vec::new();
    if let Some(&reference) = evidence.positions.first() {
        let extent = evidence
            .positions
            .iter()
            .map(|&point| (point - reference).magnitude())
            .fold(0.0_f64, f64::max);
        let resolution = extent * options.axis_relative_resolution;
        if let Some(discovery) = discover_circle_axes(
            evidence,
            CircleAxisOptions {
                vertex_tolerance: resolution,
                relative_tolerance: options.relative_surface_tolerance,
                axis_position_tolerance: 8.0 * resolution,
                axis_angle_tolerance_radians: options.axis_angle_tolerance_radians,
            },
        ) {
            let mut added = 0;
            for (id, group) in discovery.axes.iter().enumerate() {
                if let Some(placed) = axes.iter_mut().find(|placed| {
                    group.axis.cross(placed.axis).magnitude() <= options.axis_relative_resolution
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
            circle_summary = Some(CircleAxisSummary {
                source_triples: discovery.source_triples,
                circles: discovery.circles.len(),
                repeated_axes: discovery.axes.len(),
                axes_added: added,
            });
            circles = discovery.circles;
        }
    }
    (circle_summary, circles)
}

fn graph_options(
    evidence: &MeshEvidence,
    placed: &PlacedAxis,
    options: AnalyticRecoveryOptions,
) -> MeridianGraphOptions {
    let extent = evidence
        .positions
        .iter()
        .map(|&p| (p - placed.origin).magnitude())
        .fold(0.0_f64, f64::max);
    let resolution = extent * options.axis_relative_resolution;
    MeridianGraphOptions {
        station_resolution: resolution,
        vertex_tolerance: 4.0 * resolution,
        relative_tolerance: options.relative_surface_tolerance,
        max_linear_relative_chord_error: options.max_linear_relative_chord_error,
        max_circular_relative_chord_error: options.max_circular_relative_chord_error,
        normal_tolerance_radians: options.profile_normal_tolerance_radians,
    }
}

fn curved_proposal<'a>(
    boundary: &RegionBoundaryIndex<'a>,
    candidate: &StationMeridianCandidate,
    options: AnalyticRecoveryOptions,
) -> Option<SurfaceRegionProposal<'a>> {
    let scale = candidate.fit.geometric_scale();
    boundary.measure_surface_region(
        &candidate.supports,
        SurfaceCarrier::from_circular_meridian(&candidate.fit),
        SurfaceRegionBudget {
            vertex: options.relative_surface_tolerance * scale,
            chord: candidate.chord_tolerance,
        },
    )
}

#[derive(Default)]
struct AxisCandidates {
    axes: Vec<AxisDiscoverySummary>,
    circular: GeometryBuckets<StationMeridianCandidate>,
    linear: GeometryBuckets<StationLinearCandidate>,
}

impl AxisCandidates {
    fn record_unsampled(&mut self, placed: &PlacedAxis, stations: usize) {
        // Non-circular axes without candidates were never reported.
        if placed.circular_seed {
            self.axes.push(AxisDiscoverySummary {
                source_track: placed.source_track,
                source_circle_group: placed.source_circle_group,
                axis: placed.axis,
                origin: placed.origin,
                stations,
                linear_candidates: 0,
                circular_candidates: 0,
                repeated_transitions: 0,
                seed_count: 0,
                profile_count: 0,
            });
        }
    }

    fn discover(
        &mut self,
        placed: &PlacedAxis,
        graph: &MeridianStationGraph<'_>,
        boundary: &RegionBoundaryIndex<'_>,
        evidence: &MeshEvidence,
        tolerance: f64,
    ) {
        let linear_candidates = graph.discover_linear(boundary);
        let linear_count = linear_candidates.len();
        for candidate in linear_candidates {
            self.linear.insert(candidate, evidence, tolerance);
        }
        if !placed.circular_seed {
            if linear_count > 0 {
                self.axes.push(AxisDiscoverySummary {
                    source_track: placed.source_track,
                    source_circle_group: placed.source_circle_group,
                    axis: placed.axis,
                    origin: placed.origin,
                    stations: graph.stations().len(),
                    linear_candidates: linear_count,
                    circular_candidates: 0,
                    repeated_transitions: 0,
                    seed_count: 0,
                    profile_count: 0,
                });
            }
            return;
        }
        let found = graph.discover_circular(boundary);
        self.axes.push(AxisDiscoverySummary {
            source_track: placed.source_track,
            source_circle_group: placed.source_circle_group,
            axis: placed.axis,
            origin: placed.origin,
            stations: graph.stations().len(),
            repeated_transitions: graph
                .transitions()
                .iter()
                .filter(|transition| transition.repeated)
                .count(),
            seed_count: found.seed_count,
            profile_count: found.profile_count,
            circular_candidates: found.candidates.len(),
            linear_candidates: linear_count,
        });
        for candidate in found.candidates {
            self.circular.insert(candidate, evidence, tolerance);
        }
    }

    fn finish(self) -> CandidateDiscovery {
        let mut circular = self.circular.into_candidates();
        circular.sort_unstable_by(|a, b| {
            b.source_polygons
                .cmp(&a.source_polygons)
                .then_with(|| {
                    a.fit
                        .max_vertex_residual
                        .total_cmp(&b.fit.max_vertex_residual)
                })
                .then_with(|| a.supports.cmp(&b.supports))
        });
        CandidateDiscovery {
            axes: self.axes,
            circular,
            linear: self.linear.into_candidates(),
        }
    }
}

fn choose_partition(
    evidence: &MeshEvidence,
    planar_baseline: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    candidates: &CandidateDiscovery,
    grown: &[CurvatureGrowthHypothesis],
    circles: &[SourceCircleObservation],
    options: AnalyticRecoveryOptions,
) -> Result<ProposalSelection, AnalyticRecoveryError> {
    let absolute_tolerance = options.trim_options(evidence).absolute_tolerance;
    let (mut proposals, expanded) =
        circular_proposals(boundary, &candidates.circular, options, absolute_tolerance);
    let coincidence_tolerance = absolute_tolerance.min(KERNEL_EPSILON * 0.1);
    let curved_count = proposals.len();
    let measured_linear =
        append_station_linear_proposals(boundary, &candidates.linear, &mut proposals);
    let direct_linear_count = proposals.len() - curved_count;
    let growth_proposals =
        measure_growth_proposals(evidence, planar_baseline, boundary, grown, options);
    append_distinct_proposals(
        evidence,
        growth_proposals,
        coincidence_tolerance,
        &mut proposals,
    );
    let core_count = proposals.len();
    append_distinct_proposals(evidence, expanded, coincidence_tolerance, &mut proposals);
    let expanded_circular_count = proposals.len() - core_count;
    let expanded_linear = boundary.expand_station_linear_regions(
        &measured_linear,
        RimLinearOptions {
            absolute_tolerance,
            relative_tolerance: options.relative_surface_tolerance,
            angular_tolerance: options.profile_normal_tolerance_radians,
            chord_fraction: options.max_linear_relative_chord_error,
        },
    );
    append_distinct_proposals(
        evidence,
        expanded_linear,
        coincidence_tolerance,
        &mut proposals,
    );
    let expanded_linear_count = proposals.len() - core_count - expanded_circular_count;
    let rim_linear_candidates = append_rim_proposals(
        evidence,
        boundary,
        circles,
        options,
        absolute_tolerance,
        &mut proposals,
    );
    let rim_linear_count =
        proposals.len() - core_count - expanded_circular_count - expanded_linear_count;
    let partition = boundary.select_surface_partition(
        &proposals,
        planar_baseline,
        options.max_partition_search_nodes,
    )?;
    let selected_for = |partition: &SurfaceRegionPartition| {
        partition
            .selected_proposals
            .iter()
            .map(|&id| SelectedSurfaceProposal::from_measured(&proposals[id]))
            .collect::<Vec<_>>()
    };
    let chosen = selected_for(&partition);
    let core_alternative = if partition
        .selected_proposals
        .iter()
        .any(|&id| id >= core_count)
    {
        let core = boundary.select_surface_partition(
            &proposals[..core_count],
            planar_baseline,
            options.max_partition_search_nodes,
        )?;
        let selected = selected_for(&core);
        Some(PartitionAlternative {
            partition: core,
            selected,
        })
    } else {
        None
    };
    let selected_curved = chosen
        .iter()
        .filter(|proposal| matches!(proposal.carrier, SurfaceCarrier::CircularRevolution { .. }))
        .count();
    Ok(ProposalSelection {
        partition,
        selected: chosen,
        curved_proposals: curved_count + expanded_circular_count,
        expanded_circular_proposals: expanded_circular_count,
        expanded_linear_proposals: expanded_linear_count,
        direct_linear_proposals: direct_linear_count,
        rim_linear_proposals: rim_linear_count,
        rim_linear_candidates,
        core_alternative,
        carrier_alternatives: CarrierAlternatives::new(&proposals),
        linear_proposals: proposals.len() - curved_count - expanded_circular_count,
        selected_curved,
    })
}

fn append_station_linear_proposals<'a, 'b>(
    boundary: &RegionBoundaryIndex<'a>,
    candidates: &'b [StationLinearCandidate],
    proposals: &mut Vec<SurfaceRegionProposal<'a>>,
) -> Vec<&'b StationLinearCandidate> {
    let mut measured_linear = Vec::new();
    for candidate in candidates {
        if let Some(proposal) = boundary.measure_surface_region(
            &candidate.supports,
            candidate.carrier,
            SurfaceRegionBudget {
                vertex: candidate.vertex_tolerance,
                chord: candidate.chord_tolerance,
            },
        ) {
            measured_linear.push(candidate);
            proposals.push(proposal);
        }
    }
    measured_linear
}

fn circular_proposals<'a>(
    boundary: &RegionBoundaryIndex<'a>,
    candidates: &[StationMeridianCandidate],
    options: AnalyticRecoveryOptions,
    absolute_tolerance: f64,
) -> (
    Vec<SurfaceRegionProposal<'a>>,
    Vec<SurfaceRegionProposal<'a>>,
) {
    let mut core = Vec::new();
    let mut expanded = Vec::new();
    for candidate in candidates {
        if let Some(proposal) = curved_proposal(boundary, candidate, options) {
            let budget = SurfaceRegionBudget {
                vertex: options.relative_surface_tolerance * candidate.fit.geometric_scale(),
                chord: candidate.chord_tolerance,
            };
            expanded.extend(boundary.expand_circular_region(
                &proposal,
                budget,
                options.profile_normal_tolerance_radians,
            ));
            core.push(proposal);
        }
    }
    // Spheres already have a placed-carrier-independent core/growth path.
    core.extend(boundary.discover_spherical_regions(
        options.relative_surface_tolerance,
        absolute_tolerance,
        options.max_circular_relative_chord_error,
    ));
    (core, expanded)
}

fn append_rim_proposals<'a>(
    evidence: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'a>,
    circles: &[SourceCircleObservation],
    options: AnalyticRecoveryOptions,
    absolute_tolerance: f64,
    proposals: &mut Vec<SurfaceRegionProposal<'a>>,
) -> Vec<SelectedSurfaceProposal> {
    let rim_candidates = boundary.discover_rim_linear_regions(
        circles,
        RimLinearOptions {
            absolute_tolerance,
            relative_tolerance: options.relative_surface_tolerance,
            angular_tolerance: options.profile_normal_tolerance_radians,
            chord_fraction: options.max_linear_relative_chord_error,
        },
    );
    let first = proposals.len();
    append_distinct_proposals(
        evidence,
        rim_candidates,
        absolute_tolerance.min(KERNEL_EPSILON * 0.1),
        proposals,
    );
    proposals[first..]
        .iter()
        .map(SelectedSurfaceProposal::from_measured)
        .collect()
}

fn measure_growth_proposals<'a>(
    evidence: &MeshEvidence,
    planar_baseline: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'a>,
    grown: &[CurvatureGrowthHypothesis],
    options: AnalyticRecoveryOptions,
) -> Vec<SurfaceRegionProposal<'a>> {
    let mut proposals = Vec::new();
    for candidate in grown {
        let fit = &candidate.spatial;
        let scale = fit.mean_sample_radius.min(fit.axial_coordinate_span);
        if scale <= 0.0
            || fit.max_sample_radial_residual > options.relative_surface_tolerance * scale
            || fit.max_support_plane_offset_residual
                > options.relative_surface_tolerance * fit.mean_sample_radius
            || fit.support_normal_axial_projection_spread > options.relative_surface_tolerance
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
        let Some(carrier) = SurfaceCarrier::from_spatial_fit(fit) else {
            continue;
        };
        if let Some(proposal) = boundary.measure_surface_region(
            &supports,
            carrier,
            SurfaceRegionBudget {
                vertex: options.relative_surface_tolerance * scale,
                chord: options.max_linear_relative_chord_error
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
    proposals
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_mesh::{Faces, StandardAttributes};

    #[test]
    fn recovers_closed_cube_without_supplied_carriers() -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                    Point3::new(1.0, 0.0, 1.0),
                    Point3::new(1.0, 1.0, 1.0),
                    Point3::new(0.0, 1.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([
                [3, 2, 1, 0],
                [4, 5, 6, 7],
                [0, 1, 5, 4],
                [1, 2, 6, 5],
                [2, 3, 7, 6],
                [3, 0, 4, 7],
            ]),
        )?;
        let recovered = recover_analytic_mesh(&mesh, AnalyticRecoveryOptions::default())?;
        let shell = recovered.complete_shell().ok_or("complete cube")?;
        assert_eq!(shell.shell_condition(), ShellCondition::Closed);
        assert_eq!(shell.len(), 6);
        let volume = recovered
            .evidence
            .closed_source_volume()
            .ok_or("source volume")?;
        assert!((volume - 1.0).abs() < 1.0e-12);
        assert_eq!(
            recovered.evidence.provenance.original_polygon_ids,
            (0..6).collect::<Vec<_>>()
        );
        Ok(())
    }

    pub(super) fn sampled_cylinder_mesh() -> Result<PolygonMesh, Box<dyn std::error::Error>> {
        let n = 32_usize;
        let mut positions = Vec::new();
        for z in [0.0, 2.0] {
            for i in 0..n {
                let theta = std::f64::consts::TAU * i as f64 / n as f64;
                positions.push(Point3::new(theta.cos(), theta.sin(), z));
            }
        }
        let mut polygons = (0..n)
            .map(|i| {
                let j = (i + 1) % n;
                vec![i, j, j + n, i + n]
            })
            .collect::<Vec<_>>();
        polygons.push((0..n).rev().collect());
        polygons.push((n..2 * n).collect());
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(polygons),
        )?)
    }

    #[test]
    fn recovers_sampled_cylinder_without_supplied_axis() -> Result<(), Box<dyn std::error::Error>> {
        let mesh = sampled_cylinder_mesh()?;
        let recovered = recover_analytic_mesh(&mesh, AnalyticRecoveryOptions::default())?;
        let shell = recovered.complete_shell().ok_or("complete cylinder")?;
        assert_eq!(shell.shell_condition(), ShellCondition::Closed);
        assert_eq!(shell.len(), 4);
        assert_eq!(recovered.discovery.selected.len(), 1);
        Ok(())
    }
}
