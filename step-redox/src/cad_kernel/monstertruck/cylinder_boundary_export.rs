//! Export a regular crossing-cylinder branch without reintroducing marching
//! chord error. Recognition of the rational carriers is independent of the
//! evaluated SSI points. The branch is selected using the complete directed
//! polyline, and its new cubic leader has an independent analytic error bound.

use super::{validated_step_parameter_curve, validated_step_surface_support};
use monstertruck_core::tolerance::{
    DEFAULT_NUMERICAL_TOLERANCES, DEFAULT_STEP_RECOVERY_TOLERANCES,
};
use monstertruck_geometry::{
    analytic_cylinder_intersection::{BoundedSourceIntersection, CylinderCarrier},
    prelude::TryIntoHomogeneousBsplineSurface,
};
use monstertruck_modeling::{
    BoundedCurve, BsplineCurve, Curve, InnerSpace, KnotVector, NurbsCurve, Surface, SurfaceCurve,
    algo,
};

fn carrier(surface: &Surface) -> Option<CylinderCarrier> {
    if !surface.homogeneous_net_covers_domain() {
        return None;
    }
    CylinderCarrier::from_homogeneous_surface(surface.try_into_homogeneous_bspline_surface()?)
}

/// This operates only on the exchange copy: the original evaluated topology,
/// surface charts, and all directed vertex positions are retained. PCURVEs
/// have independent parameter ranges and must still agree geometrically.
pub(super) fn bounded_crossing_cylinder_step_guide(curve: &Curve, tolerance: f64) -> Option<Curve> {
    let Curve::Intersection(intersection) = curve else {
        return None;
    };
    let carriers = [
        carrier(intersection.surface0())?,
        carrier(intersection.surface1())?,
    ];
    let proof_tolerance = DEFAULT_STEP_RECOVERY_TOLERANCES
        .solid_geometry_mm
        .min(tolerance);
    // Use a finer sampling budget than the requested export accuracy, since
    // the whole-polyline bound amplifies carrier errors near tangency.
    let (parameters, points) =
        algo::curve::parameter_division(intersection, intersection.range_tuple(), tolerance / 16.0);
    if points.len() < 2 || parameters.len() != points.len() {
        return None;
    }
    let knots = KnotVector::from(
        std::iter::once(*parameters.first()?)
            .chain(parameters)
            .chain(std::iter::once(intersection.range_tuple().1))
            .collect::<Vec<_>>(),
    );
    let guide = NurbsCurve::from(BsplineCurve::new(knots, points).ok()?);
    let mut selected = None::<BoundedSourceIntersection>;
    for (primary, secondary) in [(carriers[0], carriers[1]), (carriers[1], carriers[0])] {
        let Some(bound) = BoundedSourceIntersection::from_spline(
            primary,
            secondary,
            &guide,
            proof_tolerance,
            tolerance / 64.0,
        ) else {
            continue;
        };
        if bound.separation_upper <= tolerance
            && selected.is_none_or(|old| bound.separation_upper < old.separation_upper)
        {
            selected = Some(bound);
        }
    }
    let leader = selected?.proposal.arc.cubic_leader(
        proof_tolerance.min(DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_structural_length) * 0.5,
        DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_leader_max_segments,
    )?;
    let mut cubic = leader.curve;
    let start = guide.front();
    let end = guide.back();
    let endpoint_error = (cubic.front() - start)
        .magnitude()
        .max((cubic.back() - end).magnitude());
    if !endpoint_error.is_finite() || leader.maximum_error + endpoint_error > proof_tolerance {
        return None;
    }
    // The clamped polynomial basis is a nonnegative partition of unity.
    // Endpoint control displacement therefore consumes at most this much
    // of the whole-curve approximation budget.
    let last = cubic.control_points().len().checked_sub(1)?;
    *cubic.control_point_mut(0) = start;
    *cubic.control_point_mut(last) = end;
    let leader = Curve::NurbsCurve(NurbsCurve::from(cubic));
    let (_, samples) =
        algo::curve::parameter_division(&leader, leader.range_tuple(), tolerance / 16.0);
    if ![intersection.surface0(), intersection.surface1()]
        .iter()
        .all(|carrier| validated_step_surface_support(carrier, &samples, proof_tolerance))
        || ![intersection.boundary0(), intersection.boundary1()]
            .into_iter()
            .flatten()
            .all(|boundary| validated_step_parameter_curve(boundary, &samples, tolerance))
    {
        return None;
    }
    let mut exported = SurfaceCurve::with_boundaries(
        intersection.surface0().clone(),
        intersection.surface1().clone(),
        Box::new(leader),
        intersection.boundary0().cloned(),
        intersection.boundary1().cloned(),
    );
    let (front, back) = intersection.endpoint_anchors();
    exported.set_endpoint_anchors(front, back);
    Some(Curve::Intersection(exported))
}
