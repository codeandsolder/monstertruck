//! Whole-curve geometric locus coverage of a finite straight segment.
//!
//! Uses positive rational Bernstein convex-hull containment, not samples.
//! Connected pieces (or pieces with certified small junction gaps) that stay
//! in a segment's tolerance tube and reach both endpoints cover its locus.

use crate::analytic_bezier_partition::proof_bezier_spans;
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
#[must_use]
pub fn certify_nurbs_segment_locus(
    curve: &NurbsCurve<Vector4>,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> SegmentLocus {
    if !tolerance.is_finite()
        || tolerance <= 0.0
        || !finite(start)
        || !finite(end)
        || !valid(curve.non_rationalized())
    {
        return SegmentLocus::Unsupported;
    }
    let chord = length(start, end);
    if !chord.is_finite() || chord <= 2.0 * tolerance {
        return SegmentLocus::Unsupported;
    }
    let Some(pieces) = proof_bezier_spans(curve.non_rationalized()) else {
        return SegmentLocus::Unsupported;
    };
    let mut endpoints = [None, None];
    let mut previous_end = None;
    let mut max_floor = 0.0_f64;
    let mut max_hull_deviation = 0.0_f64;
    let mut max_junction_gap = 0.0_f64;
    for (index, piece) in pieces.iter().enumerate() {
        let points = &piece.controls;
        let Some(floor) = resolution_floor(points, start, end) else {
            return SegmentLocus::Unsupported;
        };
        if floor >= tolerance {
            return SegmentLocus::Unsupported;
        }
        max_floor = max_floor.max(floor);
        let Some(front) = points.first().copied().and_then(cartesian) else {
            return SegmentLocus::Unsupported;
        };
        let Some(back) = points.last().copied().and_then(cartesian) else {
            return SegmentLocus::Unsupported;
        };
        if let Some(previous) = previous_end {
            let gap = length(previous, front);
            if !gap.is_finite() {
                return SegmentLocus::Unsupported;
            }
            max_junction_gap = max_junction_gap.max(gap);
        } else {
            endpoints[0] = Some(front);
        }
        previous_end = Some(back);
        if index + 1 == pieces.len() {
            endpoints[1] = Some(back);
        }
        for &point in points {
            let Some(point) = cartesian(point) else {
                return SegmentLocus::Unsupported;
            };
            let Some(distance) = distance_to_segment(point, start, end) else {
                return SegmentLocus::Unsupported;
            };
            max_hull_deviation = max_hull_deviation.max(distance);
        }
    }
    // Every piece stays within the segment's convex tolerance tube. Each
    // connected piece covers all intermediate projected positions; if two
    // consecutive pieces jump, their projected coverage can miss only a gap
    // bounded by the actual spatial endpoint distance. Triangle inequality
    // then bounds any missing target point by hull deviation + junction gap.
    let total_bound = max_hull_deviation + max_junction_gap + max_floor;
    if !total_bound.is_finite() {
        return SegmentLocus::Unsupported;
    }
    if total_bound > tolerance {
        return SegmentLocus::NotProven;
    }
    let [Some(front), Some(back)] = endpoints else {
        return SegmentLocus::Unsupported;
    };
    let margin = tolerance - max_floor;
    let forward = length(front, start) <= margin && length(back, end) <= margin;
    let reversed = length(front, end) <= margin && length(back, start) <= margin;
    if forward || reversed {
        SegmentLocus::WithinTolerance
    } else {
        SegmentLocus::NotProven
    }
}

/// Weight-one embedding of a polynomial spline with the same global proof.
#[must_use]
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
