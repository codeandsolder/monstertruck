//! A shared 1D hint resolver must not change either spline family's
//! nearest or exact parameter solver, or the precedence of supplied hints.
use super::*;
use monstertruck_core::math::{Point2, Vector2, Vector3};
use monstertruck_traits::{
    CurveParameter, SearchNearestParameter, SearchParameter, SearchParameterHint1D,
};

fn verify_search_hints<C>(curve: &C)
where C: BoundedCurve<Point = Point2>
        + SearchParameter<CurveParameter, Point = Point2>
        + SearchNearestParameter<CurveParameter, Point = Point2> {
    let parameter = 0.35;
    let point = curve.evaluate(parameter);
    let nearest_point = point + Vector2::new(0.0, 0.25);

    // Explicit parameters are never clamped or pre-searched, even when
    // outside the curve's natural domain.
    assert_eq!(
        resolve_curve_search_hint(curve, point, SearchParameterHint1D::Parameter(2.5)),
        2.5
    );
    for hint in [
        SearchParameterHint1D::Parameter(0.5),
        SearchParameterHint1D::Range(0.2, 0.65),
        SearchParameterHint1D::None,
    ] {
        let seed = resolve_curve_search_hint(curve, point, hint);
        assert!(seed.is_finite(), "{hint:?}");
        match hint {
            SearchParameterHint1D::Parameter(value) => assert_eq!(seed, value),
            SearchParameterHint1D::Range(_, _) | SearchParameterHint1D::None => {
                assert!(
                    (seed - parameter).abs() < 0.03,
                    "incorrect sampled seed {seed} from {hint:?}"
                );
            }
        }

        let exact = SearchParameter::<CurveParameter>::search_parameter(curve, point, hint, 100)
            .expect("exact projection");
        assert!(
            (exact - parameter).abs() < 1.0e-7,
            "wrong exact parameter {exact} from {hint:?}"
        );
        let nearest = SearchNearestParameter::<CurveParameter>::search_nearest_parameter(
            curve,
            nearest_point,
            hint,
            100,
        )
        .expect("nearest projection");
        assert!(
            (nearest - parameter).abs() < 1.0e-7,
            "wrong nearest parameter {nearest} from {hint:?}"
        );
    }

    // A restricted range must be respected when choosing a presearch seed.
    let restricted =
        resolve_curve_search_hint(curve, point, SearchParameterHint1D::Range(0.7, 0.9));
    assert!(
        (0.7..=0.9).contains(&restricted),
        "seed {restricted} escaped explicit search range"
    );
}

#[test]
fn bspline_curve_search_hints_are_shared_across_exact_and_nearest_searches() {
    let curve = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)],
    );
    verify_search_hints(&curve);
}

#[test]
fn nurbs_curve_search_hints_are_shared_across_exact_and_nearest_searches() {
    // Different positive weights make the parameterization rational,
    // exercising the projected point type rather than homogeneous points.
    let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Vector3::new(0.0, 0.0, 1.0), Vector3::new(2.0, 0.0, 2.0)],
    ));
    verify_search_hints(&curve);
}
