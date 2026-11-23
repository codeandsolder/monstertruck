//! Bounds the maximum distance of an entire spline curve from a plane.
//!
//! Every positive knot interval is converted to a homogeneous Bezier span.
//! Positive weights make its Euclidean locus a convex combination of the
//! dehomogenized controls, giving a global distance upper bound. Recursive
//! de Casteljau subdivision tightens that bound without sampling a grid.

use crate::analytic_bezier_partition::proof_bezier_spans;
use crate::analytic_bezier_refinement::refine_maximum;
use crate::analytic_curve_incidence::valid;
use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

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
    #[must_use = "The boolean result of `within` must be observed to use the predicate."]
    pub fn within(self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.upper <= tolerance
    }

    #[must_use = "The boolean result of `exceeds` must be observed to use the predicate."]
    pub fn exceeds(self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.lower > tolerance
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
#[must_use = "Ignoring this result discards a global bracket on maximum whole-curve plane distance."]
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
    let max_margin = spans.iter().try_fold(0.0_f64, |margin, span| {
        Some(margin.max(arithmetic_margin(&span.controls, origin, degree, 0)?))
    })?;
    let bracket = refine_maximum(spans, target_gap, |controls| {
        let first = signed_distance(controls.first()?, origin, normal)?.abs();
        let last = signed_distance(controls.last()?, origin, normal)?.abs();
        Some((first.max(last), upper_bound(controls, origin, normal)?))
    })?;
    // The final numerical allowance accounts for dehomogenization, dot products
    // and up to max_depth rounds of homogeneous midpoint interpolation.
    let margin = max_margin * crate::usize_to_f64(bracket.max_depth + 1);
    if !margin.is_finite() {
        return None;
    }
    let upper = bracket.upper + margin;
    let lower = (bracket.lower - margin).max(0.0);
    if !lower.is_finite() || !upper.is_finite() || lower > upper {
        return None;
    }
    Some(PlaneDeviationBounds {
        lower,
        upper,
        arithmetic_margin: margin,
        subdivisions: bracket.subdivisions,
    })
}

/// A lower bound on the directed point-set distance from one entire curve
/// to another, obtained from their distances to the same reference plane.
///
/// Both curves may have unrelated knot vectors, parameterizations, degrees,
/// and numbers of spans. A large source-plane excursion and a small target-
/// plane excursion separate their geometric loci by the reverse triangle
/// inequality. Zero is inconclusive, not an equivalence certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneDirectedSeparationBounds {
    pub source: PlaneDeviationBounds,
    pub target: PlaneDeviationBounds,
    /// A lower bound on sup_{p in source} inf_{q in target} |p - q|.
    pub minimum_source_to_target_mm: f64,
}

/// Bound source-to-target geometric disagreement without matching parameters.
///
/// The target only needs a proven global bound to the reference plane; it
/// need not have an exact planar representation. Both splines are checked over
/// their entire active knot domains by the same Bernstein partitioner.
/// This does not establish B-rep topology, edge trim, or face orientation.
#[must_use = "The optional bound nurbs curves separation by plane outcome must be inspected to distinguish a value from absence."]
pub fn bound_nurbs_curves_separation_by_plane(
    source: &NurbsCurve<Vector4>,
    target: &NurbsCurve<Vector4>,
    plane_origin: Point3,
    plane_normal: Vector3,
    target_gap: f64,
) -> Option<PlaneDirectedSeparationBounds> {
    let source = bound_nurbs_plane_deviation(source, plane_origin, plane_normal, target_gap)?;
    let target = bound_nurbs_plane_deviation(target, plane_origin, plane_normal, target_gap)?;
    let minimum_source_to_target_mm = (source.lower - target.upper).max(0.0);
    Some(PlaneDirectedSeparationBounds {
        source,
        target,
        minimum_source_to_target_mm,
    })
}

/// Polynomial B-spline specialization using positive unit weights.
#[must_use = "The optional bound bspline plane deviation outcome must be inspected to distinguish a value from absence."]
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
    fn planar_separator_proves_real_locus_disagreement_across_parameterizations() {
        let source = NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.3, 0.0, 0.0002),
                p(0.7, 0.0, 0.0002),
                p(1.0, 0.0, 0.0),
            ],
        ));
        let target = NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::from(vec![10.0, 10.0, 11.0, 12.0, 12.0]),
            vec![p(0.0, 0.0, 0.0), p(0.5, 0.0, 0.0), p(1.0, 0.0, 0.0)],
        ));
        let origin = p(0.0, 0.0, 0.0);
        let normal = Vector3::unit_z();
        let bound =
            bound_nurbs_curves_separation_by_plane(&source, &target, origin, normal, 1.0e-10)
                .expect("finite complete curves");
        assert!(bound.minimum_source_to_target_mm > 0.00014999, "{bound:?}");
        assert!(bound.target.upper < 1.0e-8, "{bound:?}");

        let reverse =
            bound_nurbs_curves_separation_by_plane(&target, &source, origin, normal, 1.0e-10)
                .expect("reverse remains defined");
        assert_eq!(reverse.minimum_source_to_target_mm, 0.0);
    }

    #[test]
    fn same_plane_is_inconclusive_even_for_different_loci() {
        let source = NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)],
        ));
        let target = NurbsCurve::from(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![p(0.0, 5.0, 0.0), p(1.0, 5.0, 0.0)],
        ));
        let bound = bound_nurbs_curves_separation_by_plane(
            &source,
            &target,
            p(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            1.0e-10,
        )
        .expect("same plane");
        assert_eq!(bound.minimum_source_to_target_mm, 0.0);
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
