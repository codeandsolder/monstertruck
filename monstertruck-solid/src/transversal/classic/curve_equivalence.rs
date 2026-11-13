use super::tolerance::BooleanTolerance;
use monstertruck_core::{cgmath64::MetricSpace, tolerance::DEFAULT_NUMERICAL_TOLERANCES};
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
