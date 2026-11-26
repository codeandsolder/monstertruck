//! Source-connected, seeded curvature region growth. Local curvature-track
//! compatibility is not transitive: every proposed extension is independently
//! re-fitted over the entire original source polygon-vertex union.

use super::super::MeshRecoveryError;
use super::join::line_distance;
use super::spatial::{SpatialTrackEvidence, fit_spatial_support_set};
use super::{CurvatureTrack, MeshEvidence, measure_curvature_spatial};
use crate::{InnerSpace, Vector3};
use rustc_hash::FxHashSet as HashSet;

/// Hypothesis budgets, never implicit original CAD primitive classifications.
#[derive(Clone, Copy, Debug)]
pub struct CurvatureGrowthOptions {
    /// Minimum actual distinct source supports for a candidate source track.
    pub min_track_supports: usize,
    /// Minimum new source supports required for each extension.
    pub min_added_supports: usize,
    /// Maximum undirected angle against the unchanged first track axis.
    pub max_seed_axis_angle_radians: f64,
    /// Maximum independently placed axis line separation / source radius.
    pub max_relative_axis_line_separation: f64,
    /// Maximum original-vertex radial fit error / mean source radius.
    pub max_relative_vertex_radial_residual: f64,
    /// Maximum supporting-plane fit error / mean source radius.
    pub max_relative_support_plane_residual: f64,
    /// Maximum spread of original source-support normal projections onto the
    /// *joint* placed axis. The normal latitude is constant for one cylinder
    /// or cone, even though a short toroidal blend can fit their vertex
    /// positions extremely closely. This is an independent necessary check
    /// for one globally constant normal-cone surface, not a torus rejection.
    pub max_support_normal_axial_projection_spread: f64,
}

impl Default for CurvatureGrowthOptions {
    fn default() -> Self {
        Self {
            min_track_supports: 8,
            min_added_supports: 1,
            max_seed_axis_angle_radians: 0.05,
            max_relative_axis_line_separation: 0.05,
            max_relative_vertex_radial_residual: 0.01,
            max_relative_support_plane_residual: 0.01,
            max_support_normal_axial_projection_spread: 0.06,
        }
    }
}

impl CurvatureGrowthOptions {
    fn valid(self) -> bool {
        self.min_track_supports >= 3
            && self.min_added_supports >= 1
            && self.max_seed_axis_angle_radians.is_finite()
            && (0.0..std::f64::consts::FRAC_PI_2).contains(&self.max_seed_axis_angle_radians)
            && [
                self.max_relative_axis_line_separation,
                self.max_relative_vertex_radial_residual,
                self.max_relative_support_plane_residual,
                self.max_support_normal_axial_projection_spread,
            ]
            .into_iter()
            .all(|x| x.is_finite() && x >= 0.0)
    }

    fn accepts(self, fit: &SpatialTrackEvidence) -> bool {
        fit.relative_max_sample_radial_residual <= self.max_relative_vertex_radial_residual
            && fit.max_support_plane_offset_residual / fit.mean_sample_radius
                <= self.max_relative_support_plane_residual
            && fit.support_normal_axial_projection_spread
                <= self.max_support_normal_axial_projection_spread
    }
}

/// One accepted extension measured against the entire global original source.
#[derive(Clone, Debug)]
pub struct CurvatureGrowthStep {
    /// Original local candidate track ID.
    pub track_index: usize,
    /// Newly covered source supports in this extension.
    pub new_supports: usize,
    /// Cumulative source-support coverage after refitting.
    pub cumulative_supports: usize,
    /// Largest original-vertex radial error relative to radius.
    pub cumulative_relative_radial_residual: f64,
    /// Largest supporting-plane placement error relative to radius.
    pub cumulative_relative_plane_residual: f64,
    /// Original support-normal latitude spread after cumulative refitting.
    pub cumulative_support_normal_axial_projection_spread: f64,
    /// Undirected angle from the unchanging initial seed axis.
    pub seed_axis_angle_radians: f64,
}

/// Candidate multi-track spatial surface region, NOT an analytic CAD-family
/// certificate. All other seed candidates remain independently available.
#[derive(Clone, Debug)]
pub struct CurvatureGrowthHypothesis {
    /// Fixed first local track, anchoring all direction comparisons.
    pub seed_track_index: usize,
    /// Original local track IDs in deterministic growth order.
    pub track_indices: Vec<usize>,
    /// Exactly retained original source support IDs, sorted.
    pub support_indices: Vec<usize>,
    /// Complete original polygon count in those source supports.
    pub source_polygon_count: usize,
    /// Maximum member divergence from initial seed axis.
    pub max_seed_axis_deviation_radians: f64,
    /// Maximum fitted input axis-line separation / source radius.
    pub max_relative_input_axis_separation: f64,
    /// Joint placed axis/radial/plane model fitted over all original vertices.
    pub spatial: SpatialTrackEvidence,
    /// Individual extensions with their *cumulative* source fit metrics.
    pub steps: Vec<CurvatureGrowthStep>,
}

struct GrowthState {
    hypothesis: CurvatureGrowthHypothesis,
    selected: HashSet<usize>,
    axis_weighted_sum: Vector3,
    max_axis_uncertainty: f64,
}

struct GrowthProposal {
    track_index: usize,
    new_supports: usize,
    supports: Vec<usize>,
    spatial: SpatialTrackEvidence,
    axis_weighted_sum: Vector3,
    axis_uncertainty: f64,
    seed_angle: f64,
    relative_line_distance: f64,
    score: f64,
}

/// Index source-support contacts once; visit only source-adjacent local tracks,
/// not a quadratic all-pairs matrix. Preserves every original seed alternative.
pub struct CurvatureGrowthIndex<'a> {
    evidence: &'a MeshEvidence,
    tracks: &'a [CurvatureTrack],
    options: CurvatureGrowthOptions,
    track_supports: Vec<Vec<usize>>,
    fits: Vec<Option<SpatialTrackEvidence>>,
    tracks_by_support: Vec<Vec<usize>>,
}

impl<'a> CurvatureGrowthIndex<'a> {
    /// Prepare independently supported local curvature evidence for seed growth.
    ///
    /// # Errors
    /// Invalid numerical bounds return the `InvalidEvidenceOptions` error variant.
    pub fn new(
        evidence: &'a MeshEvidence,
        tracks: &'a [CurvatureTrack],
        options: CurvatureGrowthOptions,
    ) -> Result<Self, MeshRecoveryError> {
        if !options.valid() {
            return Err(MeshRecoveryError::InvalidEvidenceOptions);
        }
        let mut track_supports = Vec::with_capacity(tracks.len());
        let mut fits = Vec::with_capacity(tracks.len());
        let mut tracks_by_support = vec![Vec::new(); evidence.planar_supports.len()];
        for (track_id, track) in tracks.iter().enumerate() {
            let mut supports = track.support_indices.clone();
            supports.sort_unstable();
            supports.dedup();
            let spatial = (supports.len() >= options.min_track_supports
                && supports
                    .iter()
                    .all(|&id| id < evidence.planar_supports.len()))
            .then(|| measure_curvature_spatial(evidence, track))
            .flatten()
            .filter(|fit| options.accepts(fit));
            if spatial.is_some() {
                for &support in &supports {
                    tracks_by_support[support].push(track_id);
                }
            }
            track_supports.push(supports);
            fits.push(spatial);
        }
        Ok(Self {
            evidence,
            tracks,
            options,
            track_supports,
            fits,
            tracks_by_support,
        })
    }

    /// Independently eligible local hypotheses, not source CAD surface count.
    #[must_use]
    pub fn eligible_track_count(&self) -> usize {
        self.fits.iter().filter(|fit| fit.is_some()).count()
    }

    /// Grow one anchored candidate using cumulative fit at each source contact.
    /// Different seeds preserve competing hypotheses instead of mutating the
    /// source mesh or accepting a connected component via transitive closure.
    #[must_use]
    pub fn grow(&self, seed_index: usize) -> Option<CurvatureGrowthHypothesis> {
        let local = self.fits.get(seed_index)?.as_ref()?.clone();
        let supports = self.track_supports.get(seed_index)?.clone();
        let axis = self.tracks.get(seed_index)?.seed_axis;
        let mut state = GrowthState {
            selected: HashSet::from_iter([seed_index]),
            axis_weighted_sum: axis * crate::usize_to_f64(supports.len()),
            max_axis_uncertainty: self.tracks[seed_index].max_seed_axis_error_radians,
            hypothesis: CurvatureGrowthHypothesis {
                seed_track_index: seed_index,
                track_indices: vec![seed_index],
                source_polygon_count: self.source_polygon_count(&supports),
                support_indices: supports,
                max_seed_axis_deviation_radians: 0.0,
                max_relative_input_axis_separation: 0.0,
                spatial: local,
                steps: Vec::new(),
            },
        };
        while let Some(proposal) = self.best_extension(&state) {
            state.selected.insert(proposal.track_index);
            state.hypothesis.steps.push(CurvatureGrowthStep {
                track_index: proposal.track_index,
                new_supports: proposal.new_supports,
                cumulative_supports: proposal.supports.len(),
                cumulative_relative_radial_residual: proposal
                    .spatial
                    .relative_max_sample_radial_residual,
                cumulative_relative_plane_residual: proposal
                    .spatial
                    .max_support_plane_offset_residual
                    / proposal.spatial.mean_sample_radius,
                cumulative_support_normal_axial_projection_spread: proposal
                    .spatial
                    .support_normal_axial_projection_spread,
                seed_axis_angle_radians: proposal.seed_angle,
            });
            state.hypothesis.track_indices.push(proposal.track_index);
            state.hypothesis.max_seed_axis_deviation_radians = state
                .hypothesis
                .max_seed_axis_deviation_radians
                .max(proposal.seed_angle);
            state.hypothesis.max_relative_input_axis_separation = state
                .hypothesis
                .max_relative_input_axis_separation
                .max(proposal.relative_line_distance);
            state.hypothesis.source_polygon_count = self.source_polygon_count(&proposal.supports);
            state.hypothesis.support_indices = proposal.supports;
            state.hypothesis.spatial = proposal.spatial;
            state.axis_weighted_sum = proposal.axis_weighted_sum;
            state.max_axis_uncertainty = proposal.axis_uncertainty;
        }
        Some(state.hypothesis)
    }

    fn source_polygon_count(&self, supports: &[usize]) -> usize {
        supports
            .iter()
            .map(|&id| self.evidence.planar_supports[id].source_face_indices.len())
            .sum()
    }

    fn frontier(&self, state: &GrowthState) -> Vec<usize> {
        let mut reachable = HashSet::default();
        for &support in &state.hypothesis.support_indices {
            for &track in &self.tracks_by_support[support] {
                reachable.insert(track);
            }
            for &interface in &self.evidence.support_interface_indices[support] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let neighbor = if a == support { b } else { a };
                for &track in &self.tracks_by_support[neighbor] {
                    reachable.insert(track);
                }
            }
        }
        let mut result = reachable
            .into_iter()
            .filter(|track| !state.selected.contains(track))
            .collect::<Vec<_>>();
        result.sort_unstable();
        result
    }

    fn evaluate_extension(&self, state: &GrowthState, id: usize) -> Option<GrowthProposal> {
        let local = self.fits.get(id)?.as_ref()?;
        let seed = self.tracks[state.hypothesis.seed_track_index].seed_axis;
        let incoming = self.tracks[id].seed_axis;
        let dot = seed.dot(incoming);
        let seed_angle = seed.cross(incoming).magnitude().atan2(dot.abs());
        if !seed_angle.is_finite() || seed_angle > self.options.max_seed_axis_angle_radians {
            return None;
        }
        let scale = state
            .hypothesis
            .spatial
            .mean_sample_radius
            .min(local.mean_sample_radius);
        let relative_line_distance = line_distance(
            state.hypothesis.spatial.axis_origin,
            state.hypothesis.spatial.axis,
            local.axis_origin,
            local.axis,
        ) / scale;
        if !relative_line_distance.is_finite()
            || relative_line_distance > self.options.max_relative_axis_line_separation
        {
            return None;
        }
        let mut supports = state.hypothesis.support_indices.clone();
        let mut newly_added = 0;
        for &support in &self.track_supports[id] {
            if state
                .hypothesis
                .support_indices
                .binary_search(&support)
                .is_err()
            {
                supports.push(support);
                newly_added += 1;
            }
        }
        if newly_added < self.options.min_added_supports {
            return None;
        }
        supports.sort_unstable();
        let oriented = if dot < 0.0 { -incoming } else { incoming };
        let axis_weighted_sum =
            state.axis_weighted_sum + oriented * crate::usize_to_f64(newly_added);
        let axis = axis_weighted_sum.normalize();
        let axis_uncertainty = state
            .max_axis_uncertainty
            .max(self.tracks[id].max_seed_axis_error_radians)
            .max(seed_angle);
        let spatial = fit_spatial_support_set(self.evidence, &supports, axis, axis_uncertainty)?;
        if !self.options.accepts(&spatial) {
            return None;
        }
        let score = (spatial.support_normal_axial_projection_spread
            / self
                .options
                .max_support_normal_axial_projection_spread
                .max(f64::MIN_POSITIVE))
        .max(
            spatial.relative_max_sample_radial_residual
                / self
                    .options
                    .max_relative_vertex_radial_residual
                    .max(f64::MIN_POSITIVE),
        )
        .max(
            (spatial.max_support_plane_offset_residual / spatial.mean_sample_radius)
                / self
                    .options
                    .max_relative_support_plane_residual
                    .max(f64::MIN_POSITIVE),
        );
        Some(GrowthProposal {
            track_index: id,
            new_supports: newly_added,
            supports,
            spatial,
            axis_weighted_sum,
            axis_uncertainty,
            seed_angle,
            relative_line_distance,
            score,
        })
    }

    fn best_extension(&self, state: &GrowthState) -> Option<GrowthProposal> {
        let mut best: Option<GrowthProposal> = None;
        for id in self.frontier(state) {
            let Some(candidate) = self.evaluate_extension(state, id) else {
                continue;
            };
            let better = best.as_ref().is_none_or(|existing| {
                candidate
                    .score
                    .total_cmp(&existing.score)
                    .then_with(|| existing.new_supports.cmp(&candidate.new_supports))
                    .then_with(|| candidate.track_index.cmp(&existing.track_index))
                    .is_lt()
            });
            if better {
                best = Some(candidate);
            }
        }
        best
    }
}
