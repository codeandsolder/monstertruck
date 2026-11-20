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
    let (knots, _) = knot_vector.to_single_multi();
    // A successful projection of the first point is not proof that an empty
    // or degenerate knot-span sequence is included in the surface.
    if knots.len() < 2 {
        return false;
    }
    // Preserve the original check against the raw knot vector bounds rather
    // than the (possibly narrower) active parameter range.
    let (u_knots, v_knots) = surface_knots;
    let u_start = u_knots[0];
    let v_start = v_knots[0];
    let u_span = u_knots.range_length();
    let v_span = v_knots.range_length();
    let includes_projected_point = |point: S::Point, (u, v): (f64, f64)| {
        surface.evaluate(u, v).near(&point)
            && u >= u_start - KERNEL_EPSILON
            && u - u_start <= u_span + KERNEL_EPSILON
            && v >= v_start - KERNEL_EPSILON
            && v - v_start <= v_span + KERNEL_EPSILON
    };
    let Some(mut hint) = search(first_point, None) else {
        return false;
    };
    // The initial search seed must also satisfy the same projection and
    // parameter-domain checks as the subsequent curve samples.
    if !includes_projected_point(first_point, hint) {
        return false;
    }
    // Degree-zero B-splines are piecewise constant: test the interior of
    // every span, not just its right knot (which belongs to the next span).
    // Positive-degree curves keep exactly their original 6×degree density.
    let samples_per_span = (curve_degree * 6).max(2);
    for i in 1..knots.len() {
        for j in 1..=samples_per_span {
            let p = crate::usize_to_f64(j) / crate::usize_to_f64(samples_per_span);
            let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
            let point = curve.evaluate(t);
            hint = match search(point, Some(hint)) {
                Some(found) => found,
                None => return false,
            };
            if !includes_projected_point(point, hint) {
                return false;
            }
        }
    }
    true
}
