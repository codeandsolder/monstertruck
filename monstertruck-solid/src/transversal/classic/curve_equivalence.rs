use monstertruck_core::tolerance::Tolerance;
use monstertruck_traits::{BoundedCurve, CurveParameter, ParametricCurve, SearchNearestParameter};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum CurveCoincidence {
    Distinct,
    /// Both directions passed the bounded sampling/projection policy.
    ///
    /// This is validation evidence, not analytic curve equality.
    SampledValidation,
}

pub(super) fn validate_sampled_coincidence<P, C>(lhs: &C, rhs: &C) -> CurveCoincidence
where
    P: Tolerance + Clone,
    C: ParametricCurve<Point = P>
        + BoundedCurve<Point = P>
        + SearchNearestParameter<CurveParameter, Point = P>, {
    let samples_on = |curve: &C, other: &C| {
        let (t0, t1) = curve.range_tuple();
        [0.0, 0.25, 0.5, 0.75, 1.0].into_iter().all(|fraction| {
            let point = curve.evaluate(t0 + (t1 - t0) * fraction);
            other
                .search_nearest_parameter(point.clone(), None, 100)
                .is_some_and(|parameter| other.evaluate(parameter).near(&point))
        })
    };
    if samples_on(lhs, rhs) && samples_on(rhs, lhs) {
        CurveCoincidence::SampledValidation
    } else {
        CurveCoincidence::Distinct
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_core::cgmath64::Point3;
    use monstertruck_meshing::prelude::PolylineCurve;

    #[test]
    fn shared_endpoints_and_midpoint_do_not_imply_curve_coincidence() {
        let lhs = PolylineCurve(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
        ]);
        let rhs = PolylineCurve(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, -1.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
        ]);

        assert_eq!(
            validate_sampled_coincidence::<Point3, _>(&lhs, &rhs),
            CurveCoincidence::Distinct
        );
    }
}
