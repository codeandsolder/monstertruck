//! Bounds the maximum distance of an entire spline curve from a plane.
//!
//! Every positive knot interval is converted to a homogeneous Bezier span.
//! Positive weights make its Euclidean locus a convex combination of the
//! dehomogenized controls, giving a global distance upper bound. Recursive
//! de Casteljau subdivision tightens that bound without sampling a grid.

use crate::analytic_bezier_partition::proof_bezier_spans;
use crate::analytic_curve_incidence::valid;
use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

const MAX_SUBDIVISIONS: usize = 32_768;
const MAX_SUBDIVISION_DEPTH: usize = 48;

/// Bounds the maximum absolute signed distance of a complete spline to a plane.
///
/// `lower` is attained at one of the evaluated subdivision endpoints;
/// `upper` bounds every point on every knot span through its positive-weight
/// Bernstein hull. Both include an explicit conservative floating-point
/// conditioning allowance. This is not a certified IEEE outward-rounded
/// interval and does not establish CAD boundary equivalence on its own.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneDeviationBounds {
    pub lower: f64,
    pub upper: f64,
    pub arithmetic_margin: f64,
    pub subdivisions: usize,
}

impl PlaneDeviationBounds {
    #[must_use]
    pub fn within(self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.upper <= tolerance
    }

    #[must_use]
    pub fn exceeds(self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.lower > tolerance
    }
}

type HomogeneousControls = Vec<[f64; 4]>;
type SplitControls = (HomogeneousControls, HomogeneousControls);

struct Span {
    controls: Vec<[f64; 4]>,
    upper: f64,
    depth: usize,
    sequence: usize,
}

impl PartialEq for Span {
    fn eq(&self, other: &Self) -> bool {
        self.sequence == other.sequence && self.upper.total_cmp(&other.upper).is_eq()
    }
}

impl Eq for Span {}

impl PartialOrd for Span {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}

impl Ord for Span {
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper
            .total_cmp(&other.upper)
            .then(self.sequence.cmp(&other.sequence))
    }
}

fn signed_distance(control: &[f64; 4], origin: Point3, normal: Vector3) -> Option<f64> {
    if control[3] <= 0.0 || !control.iter().all(|value| value.is_finite()) {
        return None;
    }
    let x = control[0] / control[3] - origin.x;
    let y = control[1] / control[3] - origin.y;
    let z = control[2] / control[3] - origin.z;
    let distance = normal.x.mul_add(x, normal.y.mul_add(y, normal.z * z));
    distance.is_finite().then_some(distance)
}

fn upper_bound(controls: &[[f64; 4]], origin: Point3, normal: Vector3) -> Option<f64> {
    controls.iter().try_fold(0.0_f64, |max, point| {
        Some(max.max(signed_distance(point, origin, normal)?.abs()))
    })
}

/// Subdivide homogeneous Bernstein controls at exactly the parametric midpoint.
fn subdivide(controls: &[[f64; 4]]) -> Option<SplitControls> {
    let mut row = controls.to_vec();
    let mut left = Vec::with_capacity(row.len());
    let mut right = Vec::with_capacity(row.len());
    while let (Some(first), Some(last)) = (row.first(), row.last()) {
        left.push(*first);
        right.push(*last);
        if row.len() == 1 {
            break;
        }
        row = row
            .windows(2)
            .map(|pair| std::array::from_fn(|i| f64::midpoint(pair[0][i], pair[1][i])))
            .collect();
        if row.iter().flatten().any(|value| !value.is_finite()) {
            return None;
        }
    }
    right.reverse();
    Some((left, right))
}

fn arithmetic_margin(
    controls: &[[f64; 4]],
    origin: Point3,
    degree: usize,
    depth: usize,
) -> Option<f64> {
    let mut min_weight = f64::INFINITY;
    let mut max_weight = 0.0_f64;
    let mut scale = [origin.x.abs(), origin.y.abs(), origin.z.abs()]
        .into_iter()
        .fold(1.0_f64, f64::max);
    for control in controls {
        let weight = control[3];
        min_weight = min_weight.min(weight);
        max_weight = max_weight.max(weight);
        for value in &control[..3] {
            scale = scale.max((value / weight).abs());
        }
    }
    let ratio = max_weight / min_weight;
    let degree = crate::usize_to_f64(degree + 1);
    let depth = crate::usize_to_f64(depth + 1);
    let margin = 256.0 * f64::EPSILON * scale * ratio * degree * degree * depth;
    margin.is_finite().then_some(margin)
}

/// Return a global bracket on maximum whole-curve plane distance.
///
/// `target_gap` requests a maximum width of the geometric hull bracket,
/// before the arithmetic margin. If work/depth limits are reached, the best
/// still-valid (possibly wider) bracket is returned. A missing bracket means
/// malformed geometry, a non-unit normal, or numerically invalid arithmetic.
#[must_use]
pub fn bound_nurbs_plane_deviation(
    curve: &NurbsCurve<Vector4>,
    origin: Point3,
    normal: Vector3,
    target_gap: f64,
) -> Option<PlaneDeviationBounds> {
    if !target_gap.is_finite()
        || target_gap <= 0.0
        || ![origin.x, origin.y, origin.z].iter().all(|v| v.is_finite())
        || !normal.magnitude2().is_finite()
        || (normal.magnitude2() - 1.0).abs() > DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
        || !valid(curve.non_rationalized())
    {
        return None;
    }
    let spans = proof_bezier_spans(curve.non_rationalized())?;
    let degree = spans.first()?.controls.len().checked_sub(1)?;
    let mut nodes = BinaryHeap::with_capacity(spans.len());
    let mut lower = 0.0_f64;
    let mut max_depth = 0_usize;
    let mut max_margin = 0.0_f64;
    for (sequence, span) in spans.into_iter().enumerate() {
        let controls = span.controls;
        let front = signed_distance(controls.first()?, origin, normal)?.abs();
        let back = signed_distance(controls.last()?, origin, normal)?.abs();
        lower = lower.max(front).max(back);
        max_margin = max_margin.max(arithmetic_margin(&controls, origin, degree, 0)?);
        nodes.push(Span {
            upper: upper_bound(&controls, origin, normal)?,
            controls,
            depth: 0,
            sequence,
        });
    }
    let mut next_sequence = nodes.len();

    let mut subdivisions = 0;
    loop {
        let largest = nodes.peek()?;
        if largest.upper - lower <= target_gap
            || largest.depth >= MAX_SUBDIVISION_DEPTH
            || subdivisions >= MAX_SUBDIVISIONS
        {
            break;
        }
        let previous = nodes.pop()?;
        let depth = previous.depth + 1;
        let (left, right) = subdivide(&previous.controls)?;
        // The last control of the left piece is the *actual* curve midpoint.
        lower = lower.max(signed_distance(left.last()?, origin, normal)?.abs());
        for controls in [left, right] {
            let upper = upper_bound(&controls, origin, normal)?;
            nodes.push(Span {
                controls,
                upper,
                depth,
                sequence: next_sequence,
            });
            next_sequence += 1;
        }
        subdivisions += 1;
        max_depth = max_depth.max(depth);
    }
    // The final numerical allowance accounts for dehomogenization, dot products
    // and up to max_depth rounds of homogeneous midpoint interpolation.
    let margin = max_margin * crate::usize_to_f64(max_depth + 1);
    if !margin.is_finite() {
        return None;
    }
    let upper = nodes.peek()?.upper + margin;
    let lower = (lower - margin).max(0.0);
    if !lower.is_finite() || !upper.is_finite() || lower > upper {
        return None;
    }
    Some(PlaneDeviationBounds {
        lower,
        upper,
        arithmetic_margin: margin,
        subdivisions,
    })
}

/// Polynomial B-spline specialization using positive unit weights.
#[must_use]
pub fn bound_bspline_plane_deviation(
    curve: &BsplineCurve<Point3>,
    origin: Point3,
    normal: Vector3,
    target_gap: f64,
) -> Option<PlaneDeviationBounds> {
    bound_nurbs_plane_deviation(&NurbsCurve::from(curve.clone()), origin, normal, target_gap)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Point3 { Point3::new(x, y, z) }

    #[test]
    fn quadratic_peak_is_located_between_endpoints() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![p(0.0, 0.0, 0.0), p(0.5, 0.0, 2.0), p(1.0, 0.0, 0.0)],
        );
        let bound =
            bound_bspline_plane_deviation(&curve, p(0.0, 0.0, 0.0), Vector3::unit_z(), 1.0e-10)
                .expect("valid polynomial");
        assert!(bound.lower <= 1.0 && bound.upper >= 1.0);
        assert!(bound.upper - bound.lower < 1.0e-9, "{bound:?}");
        assert!(bound.within(1.000001));
        assert!(bound.exceeds(0.999999));
    }

    #[test]
    fn rational_and_multispan_curve_stay_globally_bounded() {
        let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 1.0),
                Vector4::new(0.25, 0.0, 0.125, 0.5),
                Vector4::new(1.0, 0.0, 0.0, 1.0),
                Vector4::new(1.5, 0.0, -0.25, 0.5),
            ],
        ));
        let bound =
            bound_nurbs_plane_deviation(&curve, p(0.0, 0.0, 0.0), Vector3::unit_z(), 1.0e-8)
                .expect("valid rational spline");
        assert!(bound.lower > 0.05, "{bound:?}");
        assert!(bound.upper - bound.lower < 1.0e-6, "{bound:?}");
        assert!(bound.exceeds(0.01));
    }

    #[test]
    fn malformed_and_nonunit_inputs_do_not_produce_certificates() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)],
        );
        assert!(
            bound_bspline_plane_deviation(
                &curve,
                p(0.0, 0.0, 0.0),
                Vector3::unit_z() * 2.0,
                1.0e-8,
            )
            .is_none()
        );
        assert!(
            bound_bspline_plane_deviation(&curve, p(0.0, 0.0, 0.0), Vector3::unit_z(), f64::NAN,)
                .is_none()
        );
    }
}
