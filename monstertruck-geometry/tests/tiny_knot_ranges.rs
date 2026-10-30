use monstertruck_geometry::prelude::*;

fn tiny_curve() -> BsplineCurve<Point3> {
    let a = 0.006_289_815_176_407_898;
    let b = 0.006_290_495_768_656_085;
    BsplineCurve::new(
        KnotVector::from(vec![a, a, a, a, b, b, b, b]),
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.5, 0.0),
            Point3::new(2.0, -0.5, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
    )
    .expect("a finite nonzero knot interval is a valid parameter domain")
}

#[test]
fn sub_tolerance_knot_range_is_not_zero_range() {
    let curve = tiny_curve();
    let (a, b) = curve.range_tuple();
    assert!(b - a < TOLERANCE);
    for fraction in [0.0, 0.17, 0.5, 0.83, 1.0] {
        let point = curve.evaluate(a + (b - a) * fraction);
        assert!(point.x.is_finite() && point.y.is_finite() && point.z.is_finite());
    }
}

#[test]
fn knot_normalization_preserves_tiny_range_curve_geometry() {
    let curve = tiny_curve();
    let (a, b) = curve.range_tuple();
    let mut normalized = curve.clone();
    normalized.knot_normalize();
    for fraction in [0.0, 0.13, 0.37, 0.71, 1.0] {
        let before = curve.evaluate(a + (b - a) * fraction);
        let after = normalized.evaluate(fraction);
        assert!(
            before.distance(after) <= 1.0e-10,
            "{fraction}: {before:?} != {after:?}"
        );
    }
}

#[test]
fn sub_tolerance_surface_knot_ranges_are_valid() {
    let ua = 3.0;
    let ub = 3.0 + 4.0e-7;
    let va = -2.0;
    let vb = -2.0 + 7.0e-7;
    let surface = BsplineSurface::new(
        (
            KnotVector::from(vec![ua, ua, ub, ub]),
            KnotVector::from(vec![va, va, vb, vb]),
        ),
        vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
    )
    .expect("finite nonzero surface parameter intervals are valid");
    let point = surface.evaluate((ua + ub) * 0.5, (va + vb) * 0.5);
    assert!(point.distance(Point3::new(0.5, 0.5, 0.0)) <= 1.0e-8);
}

#[test]
fn exactly_zero_knot_range_is_still_rejected() {
    let knots = KnotVector::from(vec![1.0, 1.0, 1.0, 1.0]);
    let points = vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)];
    assert!(matches!(
        BsplineCurve::new(knots, points),
        Err(Error::ZeroRange)
    ));
}
