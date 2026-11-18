//! Shared bounded-curve-in-surface verification for polynomial and rational
//! splines. The solver is provided by each surface because the B-spline
//! implementation explicitly presearches, while NURBS uses its optimized
//! search trait (including its cached presearch path).

use super::{KERNEL_EPSILON, KnotVector, ParametricCurve, ParametricSurface, Tolerance};

/// Check curve inclusion using the caller's exact search policy.
///
/// In particular, keep the initial seed generation separate from later
/// seeded searches, retain knot-interval sample order and tolerance tests,
/// and never accept an out-of-domain surface projection.
pub(super) fn includes_curve<S, C, F>(
    surface: &S,
    curve: &C,
    first_point: S::Point,
    knot_vector: &KnotVector,
    surface_knots: (&KnotVector, &KnotVector),
    curve_degree: usize,
    mut search: F,
) -> bool
where
    S: ParametricSurface,
    S::Point: Copy + Tolerance,
    C: ParametricCurve<Point = S::Point>,
    F: FnMut(S::Point, Option<(f64, f64)>) -> Option<(f64, f64)>,
{
    let Some(mut hint) = search(first_point, None) else {
        return false;
    };
    // Preserve the original check against the raw knot vector bounds rather
    // than the (possibly narrower) active parameter range.
    let (u_knots, v_knots) = surface_knots;
    let u_start = u_knots[0];
    let v_start = v_knots[0];
    let u_span = u_knots.range_length();
    let v_span = v_knots.range_length();
    let degree = curve_degree * 6;
    let (knots, _) = knot_vector.to_single_multi();
    for i in 1..knots.len() {
        for j in 1..=degree {
            let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
            let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
            let point = curve.evaluate(t);
            hint = match search(point, Some(hint)) {
                Some(found) => found,
                None => return false,
            };
            if !surface.evaluate(hint.0, hint.1).near(&point)
                || hint.0 < u_start - KERNEL_EPSILON
                || hint.0 - u_start > u_span + KERNEL_EPSILON
                || hint.1 < v_start - KERNEL_EPSILON
                || hint.1 - v_start > v_span + KERNEL_EPSILON
            {
                return false;
            }
        }
    }
    true
}
