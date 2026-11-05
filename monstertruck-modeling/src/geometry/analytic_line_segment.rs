//! Native modeling adapter for whole-curve line-segment locus certificates.

use super::{BsplineCurve, Curve, KnotVector, Line, Point3};
use monstertruck_geometry::analytic_line_segment::{
    SegmentLocus, certify_bspline_segment_locus, certify_nurbs_segment_locus,
};

/// Certifies a whole native curve has the same point-set locus as a segment,
/// up to the stated geometric tolerance in both Hausdorff directions.
///
/// The result is not an oriented edge topology, simple face trim or solid
/// occupied-set equivalence certificate.
#[must_use]
pub fn certify_native_curve_segment_locus(
    curve: &Curve,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> SegmentLocus {
    match curve {
        Curve::Line(line) => {
            let Line(first, last) = *line;
            let spline = BsplineCurve::new_unchecked(KnotVector::bezier_knot(1), vec![first, last]);
            certify_bspline_segment_locus(&spline, start, end, tolerance)
        }
        Curve::BsplineCurve(spline) => certify_bspline_segment_locus(spline, start, end, tolerance),
        Curve::NurbsCurve(spline) => certify_nurbs_segment_locus(spline, start, end, tolerance),
        // These variants evaluate a composed/projected spatial curve; the
        // leader alone is not a certificate for its actual geometric image.
        Curve::Intersection(_) | Curve::ParameterCurve(_) => SegmentLocus::Unsupported,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_line_segment_locus_accepts_both_endpoint_orders() {
        let start = Point3::new(0.0, 0.0, 0.0);
        let end = Point3::new(1.0, 0.0, 0.0);
        let line = Curve::Line(Line(start, end));
        for (a, b) in [(start, end), (end, start)] {
            assert_eq!(
                certify_native_curve_segment_locus(&line, a, b, 1.0e-7),
                SegmentLocus::WithinTolerance,
            );
        }
    }

    #[test]
    fn native_parallel_displaced_line_is_not_same_locus() {
        let start = Point3::new(0.0, 0.0, 0.0);
        let end = Point3::new(1.0, 0.0, 0.0);
        let shifted = Curve::Line(Line(
            Point3::new(0.0, 0.02, 0.0),
            Point3::new(1.0, 0.02, 0.0),
        ));
        assert_eq!(
            certify_native_curve_segment_locus(&shifted, start, end, 1.0e-7),
            SegmentLocus::NotProven,
        );
    }
}
