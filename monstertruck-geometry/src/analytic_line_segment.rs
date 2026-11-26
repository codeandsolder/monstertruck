//! Whole-curve geometric locus coverage of a finite straight segment.
//!
//! Uses positive rational Bernstein convex-hull containment, not samples.
//! Connected pieces (or pieces with certified small junction gaps) that stay
//! in a segment's tolerance tube and reach both endpoints cover its locus.

use crate::analytic_bezier_partition::{ProofBezierSpan, proof_bezier_spans};
use crate::analytic_bezier_refinement::refine_maximum;
use crate::analytic_curve_incidence::valid;
use crate::prelude::*;

/// Whether the entire curve and an intended line segment have a mutually
/// bounded spatial point-set distance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegmentLocus {
    /// The Hausdorff distance of the two point sets is within tolerance,
    /// subject to the documented floating-point arithmetic margin.
    WithinTolerance,
    /// The supplied control hull does not establish that distance bound.
    NotProven,
    /// Invalid geometry or numerically unresolved control-hull bounds.
    Unsupported,
}

fn finite(point: Point3) -> bool { [point.x, point.y, point.z].into_iter().all(f64::is_finite) }

fn length(a: Point3, b: Point3) -> f64 { (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z) }

fn cartesian(p: [f64; 4]) -> Option<Point3> {
    if !p.iter().all(|x| x.is_finite()) || p[3] <= 0.0 {
        return None;
    }
    let point = Point3::new(p[0] / p[3], p[1] / p[3], p[2] / p[3]);
    finite(point).then_some(point)
}

fn distance_to_segment(point: Point3, start: Point3, end: Point3) -> Option<f64> {
    let (x, y, z) = (end.x - start.x, end.y - start.y, end.z - start.z);
    let denominator = x.mul_add(x, y.mul_add(y, z * z));
    if !denominator.is_finite() || denominator <= 0.0 {
        return None;
    }
    let position = (x.mul_add(
        point.x - start.x,
        y.mul_add(point.y - start.y, z * (point.z - start.z)),
    ) / denominator)
        .clamp(0.0, 1.0);
    if !position.is_finite() {
        return None;
    }
    let nearest = Point3::new(
        x.mul_add(position, start.x),
        y.mul_add(position, start.y),
        z.mul_add(position, start.z),
    );
    let distance = length(point, nearest);
    distance.is_finite().then_some(distance)
}

// Conservative arithmetic-conditioning margin, not a fully outward-rounded
// interval error certificate. Semantic equality must not follow from it alone.
fn resolution_floor(points: &[[f64; 4]], start: Point3, end: Point3) -> Option<f64> {
    let min_weight = points.iter().map(|p| p[3]).fold(f64::INFINITY, f64::min);
    let max_weight = points.iter().map(|p| p[3]).fold(0.0_f64, f64::max);
    if min_weight <= 0.0 || !min_weight.is_finite() {
        return None;
    }
    let ratio = max_weight / min_weight;
    let scale = points
        .iter()
        .flat_map(|p| p[..3].iter().map(move |x| x.abs() / p[3]))
        .chain([start.x.abs(), start.y.abs(), start.z.abs()])
        .chain([end.x.abs(), end.y.abs(), end.z.abs()])
        .fold(1.0_f64, f64::max);
    let degree = crate::usize_to_f64(points.len().saturating_sub(1));
    let floor = 256.0 * (degree + 1.0).powi(2) * f64::EPSILON * scale * ratio;
    floor.is_finite().then_some(floor)
}

struct SegmentEvidence {
    pieces: Vec<ProofBezierSpan>,
    floor: f64,
    hull: f64,
    gap: f64,
    endpoint_error: f64,
}

fn span_distance(points: &[[f64; 4]], start: Point3, end: Point3) -> Option<(f64, f64)> {
    let distance = |point| distance_to_segment(cartesian(point)?, start, end);
    let witness = distance(*points.first()?)?.max(distance(*points.last()?)?);
    let upper = points
        .iter()
        .try_fold(witness, |bound, &point| Some(bound.max(distance(point)?)))?;
    Some((witness, upper))
}

fn segment_evidence(
    curve: &NurbsCurve<Vector4>,
    start: Point3,
    end: Point3,
) -> Option<SegmentEvidence> {
    if !finite(start) || !finite(end) || !valid(curve.non_rationalized()) {
        return None;
    }
    let chord = length(start, end);
    if !chord.is_finite() || chord <= 0.0 {
        return None;
    }
    let pieces = proof_bezier_spans(curve.non_rationalized())?;
    let mut previous_end = None;
    let mut floor = 0.0_f64;
    let mut hull = 0.0_f64;
    let mut gap = 0.0_f64;
    for piece in &pieces {
        let points = &piece.controls;
        floor = floor.max(resolution_floor(points, start, end)?);
        let front = cartesian(*points.first()?)?;
        let back = cartesian(*points.last()?)?;
        if let Some(previous) = previous_end {
            gap = gap.max(length(previous, front));
        }
        previous_end = Some(back);
        hull = hull.max(span_distance(points, start, end)?.1);
    }
    let front = cartesian(*pieces.first()?.controls.first()?)?;
    let back = previous_end?;
    let endpoint_error = length(front, start)
        .max(length(back, end))
        .min(length(front, end).max(length(back, start)));
    [floor, hull, gap, endpoint_error]
        .into_iter()
        .all(f64::is_finite)
        .then_some(SegmentEvidence {
            pieces,
            floor,
            hull,
            gap,
            endpoint_error,
        })
}

impl SegmentEvidence {
    fn upper_bound(&self, hull: f64) -> f64 {
        // Each continuous piece covers its complete projected interval.
        // A missing interval is bounded by the adjacent junction gap.
        // Outside the covered endpoint interval, distance to the nearest
        // actual curve endpoint is bounded by its endpoint error.
        (hull + self.gap).max(self.endpoint_error) + self.floor
    }
}

/// Measure a continuous two-sided geometric deviation upper bound between
/// a complete rational spline and a finite, nonzero source segment.
///
/// Refines the Bernstein control hull if it exceeds the requested accuracy
/// target. The returned bound can exceed that target; the caller decides
/// admissibility. Positive rational weights, all original positive knot
/// intervals and endpoint/junction coverage are checked. The numerical
/// margin is not a directed-rounding interval arithmetic proof.
#[must_use = "the continuous locus bound is needed to accept a geometric approximation"]
pub fn bound_nurbs_segment_locus(
    curve: &NurbsCurve<Vector4>,
    start: Point3,
    end: Point3,
    target: f64,
) -> Option<f64> {
    if !target.is_finite() || target <= 0.0 {
        return None;
    }
    let evidence = segment_evidence(curve, start, end)?;
    let coarse = evidence.upper_bound(evidence.hull);
    if !coarse.is_finite() {
        return None;
    }
    if coarse <= target || evidence.endpoint_error + evidence.floor > target {
        return Some(coarse);
    }
    let upper_padding = evidence.gap;
    let endpoint_error = evidence.endpoint_error;
    let floor = evidence.floor;
    let refined = refine_maximum(evidence.pieces, (target * 0.01).max(floor), |points| {
        span_distance(points, start, end)
    })?;
    let bound = (refined.upper + upper_padding).max(endpoint_error) + floor;
    bound.is_finite().then_some(bound.min(coarse))
}

/// Certify that a complete rational spline and a line segment have point-set
/// Hausdorff distance no more than the requested geometric tolerance.
///
/// All knot intervals are converted into homogeneous Bezier pieces. Positive
/// weights make each piece a convex combination of its Euclidean controls.
/// If every control is close to the target segment, the whole curve is close
/// to it. Certified junction gaps and endpoint coverage give the reverse
/// inclusion, independently of parameter speed or repeated traversal.
///
/// This is a tolerance and arithmetic-conditioning claim, not a proof that
/// two CAD solids enclose the same occupied set. In particular it does not
/// certify the simplicity of an enclosing face loop.
#[must_use = "The segment-locus certificate records the source spline's bounded global deviation."]
pub fn certify_nurbs_segment_locus(
    curve: &NurbsCurve<Vector4>,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> SegmentLocus {
    if !tolerance.is_finite() || tolerance <= 0.0 || length(start, end) <= 2.0 * tolerance {
        return SegmentLocus::Unsupported;
    }
    let Some(evidence) = segment_evidence(curve, start, end) else {
        return SegmentLocus::Unsupported;
    };
    if evidence.floor >= tolerance {
        return SegmentLocus::Unsupported;
    }
    let bound = evidence.upper_bound(evidence.hull);
    if !bound.is_finite() {
        SegmentLocus::Unsupported
    } else if bound <= tolerance {
        SegmentLocus::WithinTolerance
    } else {
        SegmentLocus::NotProven
    }
}

/// Weight-one embedding of a polynomial spline with the same global proof.
#[must_use = "The segment-locus certificate records the original polynomial spline's global agreement."]
pub fn certify_bspline_segment_locus(
    curve: &BsplineCurve<Point3>,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> SegmentLocus {
    certify_nurbs_segment_locus(&NurbsCurve::from(curve.clone()), start, end, tolerance)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-7;

    fn p(x: f64, y: f64, z: f64) -> Point3 { Point3::new(x, y, z) }

    fn polynomial(points: Vec<Point3>) -> BsplineCurve<Point3> {
        BsplineCurve::new_unchecked(KnotVector::bezier_knot(points.len() - 1), points)
    }

    #[test]
    fn adaptive_bound_resolves_control_hull_overestimate() {
        let start = p(0.0, 0.0, 0.0);
        let end = p(2.0, 0.0, 0.0);
        let curve = NurbsCurve::from(polynomial(vec![start, p(1.0, 0.2, 0.0), end]));
        let bound = bound_nurbs_segment_locus(&curve, start, end, 0.12).unwrap();
        assert!((0.1..=0.12).contains(&bound));
        let tighter = bound_nurbs_segment_locus(&curve, start, end, 0.09).unwrap();
        assert!(tighter >= 0.1);
        // Unlike a carrier-identity test, a source approximation may use
        // a budget greater than half the source segment length.
        assert!(bound_nurbs_segment_locus(&curve, start, end, 2.0).unwrap() < 2.0);
    }

    #[test]
    fn curved_parameter_speed_and_nonmonotone_controls_keep_locus() {
        let curve = polynomial(vec![
            p(0.0, 0.0, 0.0),
            p(0.8, 0.0, 0.0),
            p(0.2, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
        ]);
        assert_eq!(
            certify_bspline_segment_locus(&curve, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), TOL),
            SegmentLocus::WithinTolerance,
        );
    }

    #[test]
    fn overshoot_and_between_endpoint_bump_are_not_proven() {
        let start = p(0.0, 0.0, 0.0);
        let end = p(1.0, 0.0, 0.0);
        for control in [p(4.0, 0.0, 0.0), p(0.5, 0.1, 0.0)] {
            let curve = polynomial(vec![start, control, end]);
            assert_eq!(
                certify_bspline_segment_locus(&curve, start, end, TOL),
                SegmentLocus::NotProven,
            );
        }
    }

    #[test]
    fn positive_weight_rational_line_segment_is_proven() {
        let start = p(0.0, 0.0, 0.0);
        let end = p(2.0, 0.0, 0.0);
        let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 1.0),
                Vector4::new(0.5, 0.0, 0.0, 0.5),
                Vector4::new(2.0, 0.0, 0.0, 1.0),
            ],
        ));
        assert_eq!(
            certify_nurbs_segment_locus(&curve, start, end, TOL),
            SegmentLocus::WithinTolerance,
        );
    }

    #[test]
    fn line_segment_locus_is_independent_of_curve_orientation() {
        let curve = polynomial(vec![p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0)]);
        assert_eq!(
            certify_bspline_segment_locus(&curve, p(1.0, 0.0, 0.0), p(0.0, 0.0, 0.0), TOL),
            SegmentLocus::WithinTolerance,
        );
    }

    #[test]
    fn discontinuous_knot_jump_cannot_claim_complete_segment_coverage() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.25, 0.0, 0.0),
                p(0.75, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
            ],
        );
        assert_eq!(
            certify_bspline_segment_locus(&curve, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), TOL,),
            SegmentLocus::NotProven,
        );
    }

    #[test]
    fn zero_gap_full_multiplicity_knot_is_supported() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.5, 0.0, 0.0),
                p(0.5, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
            ],
        );
        assert_eq!(
            certify_bspline_segment_locus(&curve, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), TOL),
            SegmentLocus::WithinTolerance,
        );
    }

    #[test]
    fn a_small_actual_jump_can_be_proven_within_the_same_tolerance_budget() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.5, 0.0, 0.0),
                p(0.5 + 1.0e-8, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
            ],
        );
        assert_eq!(
            certify_bspline_segment_locus(&curve, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), TOL),
            SegmentLocus::WithinTolerance,
        );
    }

    #[test]
    fn hull_deviation_and_knot_jump_share_one_error_budget() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.5, 0.0, 0.0),
                p(0.56, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
            ],
        );
        // The 0.06 gap is resolvable by a 0.1 tolerance when the
        // control hull itself lies exactly on the line.
        assert_eq!(
            certify_bspline_segment_locus(&curve, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), 0.1,),
            SegmentLocus::WithinTolerance,
        );

        let displaced = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.5, 0.08, 0.0),
                p(0.56, 0.08, 0.0),
                p(1.0, 0.0, 0.0),
            ],
        );
        // A valid bound must account for BOTH displacement and missing
        // junction coverage, not grant each a separate 0.1 allowance.
        assert_eq!(
            certify_bspline_segment_locus(&displaced, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), 0.1,),
            SegmentLocus::NotProven,
        );
    }

    #[test]
    fn malformed_weights_and_unresolved_precision_fail_closed() {
        let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 0.0),
                Vector4::new(1.0, 0.0, 0.0, 1.0),
            ],
        ));
        assert_eq!(
            certify_nurbs_segment_locus(&curve, p(0.0, 0.0, 0.0), p(1.0, 0.0, 0.0), TOL),
            SegmentLocus::Unsupported,
        );
        let huge = 1.0e12;
        let curve = polynomial(vec![p(huge, 0.0, 0.0), p(huge + 1.0, 0.0, 0.0)]);
        assert_eq!(
            certify_bspline_segment_locus(
                &curve,
                p(huge, 0.0, 0.0),
                p(huge + 1.0, 0.0, 0.0),
                1.0e-9,
            ),
            SegmentLocus::Unsupported,
        );
    }
}
