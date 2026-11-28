//! Lazy source-connected compatibility between independently observed
//! normal-space tracks. A single secant strip may fit a cone perfectly; the
//! union must instead explain *all* original source vertices with one placed
//! axis and radial profile to support a larger analytic-surface hypothesis.

use super::spatial::{SpatialTrackEvidence, measure_spatial_support_set};
use super::{CurvatureTrack, MeshEvidence, measure_curvature_spatial};
use crate::{Point3, Vector3};
use rustc_hash::FxHashSet as HashSet;

/// Evidence for whether two genuinely source-connected curvature hypotheses
/// could share one globally placed radial source surface.
///
/// Normal directions, axis positions and full original source-vertex geometry
/// are rechecked for the union. A small combined residual is still not a
/// proof of a specific CAD primitive or a consistent material boundary.
#[derive(Clone, Debug)]
pub struct SpatialTrackJoinEvidence {
    /// Canonical original planar-support IDs present in both source paths.
    pub shared_support_count: usize,
    /// Original manifold support-interface components linking the paths.
    pub connecting_interface_count: usize,
    /// Distinct canonical original source supports in the combined candidate.
    pub union_support_count: usize,
    /// Undirected angular divergence of the two recovered local axes.
    pub axis_angle_radians: f64,
    /// Minimum separation of the infinite geometric axes, independent of
    /// their arbitrarily chosen longitudinal reference points.
    pub axis_line_separation: f64,
    /// Larger normalized original-vertex radial residual of the two inputs.
    pub max_individual_relative_radial_error: f64,
    /// Larger normalized support-plane offset residual of the two inputs.
    pub max_individual_relative_plane_error: f64,
    /// New placed surface fit over every original source vertex of the union,
    /// not an average of independently fitted radii or axis locations.
    pub combined: SpatialTrackEvidence,
}

fn source_contact_count(
    evidence: &MeshEvidence,
    left: &HashSet<usize>,
    right: &HashSet<usize>,
) -> usize {
    let mut interfaces = HashSet::default();
    for &support in left {
        for &index in evidence
            .support_interface_indices
            .get(support)
            .into_iter()
            .flatten()
        {
            let Some(interface) = evidence.support_interfaces.get(index) else {
                continue;
            };
            let [a, b] = interface.supports;
            if (left.contains(&a) && right.contains(&b))
                || (left.contains(&b) && right.contains(&a))
            {
                interfaces.insert(index);
            }
        }
    }
    interfaces.len()
}

pub(super) fn line_distance(
    left_origin: Point3,
    left_axis: Vector3,
    right_origin: Point3,
    right_axis: Vector3,
) -> f64 {
    let delta = right_origin - left_origin;
    let cross = left_axis.cross(&right_axis);
    let sine = cross.norm();
    if sine < 1.0e-8 {
        delta.cross(&left_axis).norm()
    } else {
        delta.dot(&cross).abs() / sine
    }
}

/// Measure whether two source-connected normal-rotation tracks are spatially
/// compatible after refitting their *union* under one placed radial profile.
///
/// Both shared source-support IDs and true regular source interfaces count as
/// original topological contact; proximity alone cannot connect separate
/// bodies. The direction is a sign-aligned average of the observed axes,
/// with its angular disagreement reported rather than silently discarded.
///
/// A failed spatial fit or absence of original contact returns None. All
/// original tracks and source evidence remain available to other hypotheses.
#[must_use = "the joint fit evidence determines whether the regions can be combined"]
pub fn measure_curvature_spatial_join(
    evidence: &MeshEvidence,
    first: &CurvatureTrack,
    second: &CurvatureTrack,
) -> Option<SpatialTrackJoinEvidence> {
    let left = first
        .support_indices
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    let right = second
        .support_indices
        .iter()
        .copied()
        .collect::<HashSet<_>>();
    let shared_support_count = left.intersection(&right).count();
    let connecting_interface_count = source_contact_count(evidence, &left, &right);
    if shared_support_count == 0 && connecting_interface_count == 0 {
        return None;
    }
    let first_fit = measure_curvature_spatial(evidence, first)?;
    let second_fit = measure_curvature_spatial(evidence, second)?;
    let cross = first.seed_axis.cross(&second.seed_axis).norm();
    let dot = first.seed_axis.dot(&second.seed_axis);
    let axis_angle_radians = cross.atan2(dot.abs());
    let second_axis = if dot < 0.0 {
        -second.seed_axis
    } else {
        second.seed_axis
    };
    let mean_axis = (first.seed_axis + second_axis).normalize();
    let axis_line_separation = line_distance(
        first_fit.axis_origin,
        first.seed_axis,
        second_fit.axis_origin,
        second.seed_axis,
    );
    let max_individual_relative_radial_error = first_fit
        .relative_max_sample_radial_residual
        .max(second_fit.relative_max_sample_radial_residual);
    let max_individual_relative_plane_error = (first_fit.max_support_plane_offset_residual
        / first_fit.mean_sample_radius)
        .max(second_fit.max_support_plane_offset_residual / second_fit.mean_sample_radius);
    let mut union = left.union(&right).copied().collect::<Vec<_>>();
    union.sort_unstable();
    let uncertainty = first
        .max_seed_axis_error_radians
        .max(second.max_seed_axis_error_radians)
        .max(0.5 * axis_angle_radians);
    let combined = measure_spatial_support_set(evidence, &union, mean_axis, uncertainty)?;
    Some(SpatialTrackJoinEvidence {
        shared_support_count,
        connecting_interface_count,
        union_support_count: union.len(),
        axis_angle_radians,
        axis_line_separation,
        max_individual_relative_radial_error,
        max_individual_relative_plane_error,
        combined,
    })
}
