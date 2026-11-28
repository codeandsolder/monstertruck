use super::tolerance::BooleanTolerance;
use monstertruck_core::{math::MetricSpace, tolerance::DEFAULT_NUMERICAL_TOLERANCES};
use monstertruck_traits::{BoundedCurve, CurveParameter, ParametricCurve, SearchNearestParameter};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum CurveCoincidence {
    Distinct,
    /// Both directions passed the bounded sampling/projection policy.
    ///
    /// This is validation evidence, not analytic curve equality.
    SampledValidation,
}

pub(super) fn validate_sampled_coincidence<P, C>(
    lhs: &C,
    rhs: &C,
    tolerances: BooleanTolerance,
) -> CurveCoincidence
where
    P: MetricSpace<Metric = f64> + Clone,
    C: ParametricCurve<Point = P>
        + BoundedCurve<Point = P>
        + SearchNearestParameter<CurveParameter, Point = P>,
{
    let samples_on = |curve: &C, other: &C| {
        let (t0, t1) = curve.range_tuple();
        DEFAULT_NUMERICAL_TOLERANCES
            .boolean_coincidence_fractions
            .into_iter()
            .all(|fraction| {
                let point = curve.evaluate(f64::mul_add(t1 - t0, fraction, t0));
                other
                    .search_nearest_parameter(point.clone(), None, 100)
                    .is_some_and(|parameter| {
                        other.evaluate(parameter).distance(point.clone()) <= tolerances.length
                    })
            })
    };
    if samples_on(lhs, rhs) && samples_on(rhs, lhs) {
        CurveCoincidence::SampledValidation
    } else {
        CurveCoincidence::Distinct
    }
}

/// Validate that the entire trimmed `candidate` curve lies on `container`.
///
/// Unlike full coincidence this intentionally permits a partitioned sub-edge of
/// the same geometric carrier. It is used only where topology has already shown
/// that both curves contain the same imprint point; distinct crossings therefore
/// remain fail-closed.
pub(super) fn validate_sampled_containment<P, C>(
    container: &C,
    candidate: &C,
    tolerances: BooleanTolerance,
) -> CurveCoincidence
where
    P: MetricSpace<Metric = f64> + Clone,
    C: ParametricCurve<Point = P>
        + BoundedCurve<Point = P>
        + SearchNearestParameter<CurveParameter, Point = P>,
{
    let (t0, t1) = candidate.range_tuple();
    let all_samples_contained = DEFAULT_NUMERICAL_TOLERANCES
        .boolean_containment_fractions
        .into_iter()
        .all(|fraction| {
            let point = candidate.evaluate(f64::mul_add(t1 - t0, fraction, t0));
            container
                .search_nearest_parameter(point.clone(), None, 100)
                .is_some_and(|parameter| {
                    container.evaluate(parameter).distance(point.clone()) <= tolerances.length
                })
        });
    if all_samples_contained {
        CurveCoincidence::SampledValidation
    } else {
        CurveCoincidence::Distinct
    }
}

/// Certify an edge pair whose *topological* endpoints have already been
/// matched. Curve inverse solvers need not converge at a hard trim endpoint;
/// those locations are checked directly against the canonical vertices.
/// Symmetric interior projections still have to pass the geometric tolerance.
pub(super) fn validate_shared_endpoint_interior_coincidence<P, C>(
    lhs: &C,
    rhs: &C,
    endpoints: &[P; 2],
    tolerances: BooleanTolerance,
) -> CurveCoincidence
where
    P: MetricSpace<Metric = f64> + Clone,
    C: ParametricCurve<Point = P>
        + BoundedCurve<Point = P>
        + SearchNearestParameter<CurveParameter, Point = P>,
{
    if endpoints[0].clone().distance(endpoints[1].clone()) <= tolerances.length {
        return CurveCoincidence::Distinct;
    }
    let endpoint_evidence = |curve: &C| {
        let (start, end) = curve.range_tuple();
        if !start.is_finite() || !end.is_finite() || start >= end {
            return false;
        }
        let first = curve.evaluate(start);
        let last = curve.evaluate(end);
        (first.clone().distance(endpoints[0].clone()) <= tolerances.length
            && last.clone().distance(endpoints[1].clone()) <= tolerances.length)
            || (first.distance(endpoints[1].clone()) <= tolerances.length
                && last.distance(endpoints[0].clone()) <= tolerances.length)
    };
    if !endpoint_evidence(lhs) || !endpoint_evidence(rhs) {
        return CurveCoincidence::Distinct;
    }
    let check_interior = |source: &C, target: &C| {
        let (start, end) = source.range_tuple();
        DEFAULT_NUMERICAL_TOLERANCES
            .boolean_containment_fractions
            .into_iter()
            .filter(|&fraction| fraction > 0.0 && fraction < 1.0)
            .all(|fraction| {
                let point = source.evaluate((end - start).mul_add(fraction, start));
                target
                    .search_nearest_parameter(point.clone(), None, 100)
                    .is_some_and(|parameter| {
                        target.evaluate(parameter).distance(point.clone()) <= tolerances.length
                    })
            })
    };
    if check_interior(lhs, rhs) && check_interior(rhs, lhs) {
        CurveCoincidence::SampledValidation
    } else {
        CurveCoincidence::Distinct
    }
}

/// Validate two curves by normalized parameter samples in either traversal.
///
/// The caller must already have established that the two edge endpoints are the
/// same canonical topological pair. This fallback exists for curve types
/// (notably SSI curves) whose inverse/nearest-parameter solver can fail even
/// when their forward evaluations agree. It is sampled validation, not analytic
/// equality.
pub(super) fn validate_endpoint_aligned_coincidence<P, C>(
    lhs: &C,
    rhs: &C,
    tolerances: BooleanTolerance,
) -> CurveCoincidence
where
    P: MetricSpace<Metric = f64>,
    C: ParametricCurve<Point = P> + BoundedCurve<Point = P>,
{
    let (l0, l1) = lhs.range_tuple();
    let (r0, r1) = rhs.range_tuple();
    let fractions = DEFAULT_NUMERICAL_TOLERANCES.boolean_containment_fractions;
    let aligned = fractions.into_iter().all(|fraction| {
        lhs.evaluate(f64::mul_add(l1 - l0, fraction, l0))
            .distance(rhs.evaluate(f64::mul_add(r1 - r0, fraction, r0)))
            <= tolerances.length
    });
    let reversed = fractions.into_iter().all(|fraction| {
        lhs.evaluate(f64::mul_add(l1 - l0, fraction, l0))
            .distance(rhs.evaluate(f64::mul_add(r1 - r0, -fraction, r1)))
            <= tolerances.length
    });
    if aligned || reversed {
        CurveCoincidence::SampledValidation
    } else {
        CurveCoincidence::Distinct
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_core::math::Point3;
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
            validate_sampled_coincidence::<Point3, _>(
                &lhs,
                &rhs,
                BooleanTolerance::new(DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length)
            ),
            CurveCoincidence::Distinct
        );
        assert_eq!(
            validate_endpoint_aligned_coincidence::<Point3, _>(
                &lhs,
                &rhs,
                BooleanTolerance::new(DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length)
            ),
            CurveCoincidence::Distinct
        );
    }
    #[test]
    fn sampled_containment_accepts_subcurve_but_rejects_crossing_curve() {
        let full = PolylineCurve(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ]);
        let sub = PolylineCurve(vec![
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 0.0),
        ]);
        let crossing = PolylineCurve(vec![
            Point3::new(1.0, -1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ]);

        assert_eq!(
            validate_sampled_containment::<Point3, _>(
                &full,
                &sub,
                BooleanTolerance::new(DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length)
            ),
            CurveCoincidence::SampledValidation
        );
        assert_eq!(
            validate_sampled_containment::<Point3, _>(
                &sub,
                &full,
                BooleanTolerance::new(DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length)
            ),
            CurveCoincidence::Distinct
        );
        assert_eq!(
            validate_sampled_containment::<Point3, _>(
                &full,
                &crossing,
                BooleanTolerance::new(DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length)
            ),
            CurveCoincidence::Distinct
        );
    }

    #[test]
    fn certified_common_endpoints_do_not_hide_distinct_curve_interiors() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let lhs = PolylineCurve(vec![a, Point3::new(0.1, 0.0, 0.0), b]);
        let same = PolylineCurve(vec![a, Point3::new(0.9, 0.0, 0.0), b]);
        let bowed = PolylineCurve(vec![a, Point3::new(0.5, 0.03, 0.0), b]);
        let tol = BooleanTolerance::new(1.0e-6);
        assert_eq!(
            validate_shared_endpoint_interior_coincidence(&lhs, &same, &[a, b], tol),
            CurveCoincidence::SampledValidation
        );
        assert_eq!(
            validate_shared_endpoint_interior_coincidence(&lhs, &bowed, &[a, b], tol),
            CurveCoincidence::Distinct
        );
        assert_eq!(
            validate_shared_endpoint_interior_coincidence(&lhs, &same, &[a, a], tol),
            CurveCoincidence::Distinct
        );
    }

    #[test]
    fn sampled_coincidence_distinguishes_sub_epsilon_offset_curves() {
        let gap = 5.0e-7;
        let lhs = PolylineCurve(vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]);
        let rhs = PolylineCurve(vec![Point3::new(0.0, 0.0, gap), Point3::new(1.0, 0.0, gap)]);
        assert_eq!(
            validate_sampled_coincidence::<Point3, _>(&lhs, &rhs, BooleanTolerance::new(1.0e-8)),
            CurveCoincidence::Distinct,
        );
        assert_eq!(
            validate_sampled_coincidence::<Point3, _>(&lhs, &rhs, BooleanTolerance::new(1.0e-6)),
            CurveCoincidence::SampledValidation,
        );
    }

    #[test]
    fn endpoint_aligned_validation_accepts_reversed_parameterization() {
        let lhs = PolylineCurve(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ]);
        let rhs = PolylineCurve(lhs.iter().rev().copied().collect());

        assert_eq!(
            validate_endpoint_aligned_coincidence::<Point3, _>(
                &lhs,
                &rhs,
                BooleanTolerance::new(DEFAULT_NUMERICAL_TOLERANCES.boolean_default_numeric_length)
            ),
            CurveCoincidence::SampledValidation
        );
    }
}
