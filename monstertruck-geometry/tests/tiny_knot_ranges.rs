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
    assert!(b - a < KERNEL_EPSILON);
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

#[test]
fn exact_multiplicity_and_compression_preserve_close_distinct_knots() {
    let a = 0.25_f64;
    let b = a + 1.0e-12;
    let c = b + 1.0e-12;
    let original = vec![a, a, a, b, c, 0.75, 1.0, 1.0, 1.0];
    let knots = KnotVector::from(original.clone());
    assert_eq!(knots.multiplicity(0), 3);
    assert_eq!(knots.multiplicity(3), 1);
    assert_eq!(knots.multiplicity(4), 1);
    let (distinct, counts) = knots.to_single_multi();
    assert_eq!(counts, vec![3, 1, 1, 1, 3]);
    let restored = KnotVector::from_single_multi(distinct, counts).unwrap();
    assert!(
        restored
            .iter()
            .map(|knot| knot.to_bits())
            .eq(original.iter().map(|knot| knot.to_bits()))
    );
    let zero = KnotVector::from(vec![-0.0, 0.0, 1.0]);
    assert_eq!(zero.multiplicity(0), 2);
    assert_eq!(zero.to_single_multi().1, vec![2, 1]);
}

#[test]
fn short_parameter_span_knot_insertion_and_cut_preserve_geometry() {
    let curve = tiny_curve();
    let (a, b) = curve.range_tuple();
    let first = a + (b - a) * 0.25;
    let second = a + (b - a) * 0.5;
    let mut left = curve.clone();
    left.add_knot(first);
    left.add_knot(second);
    let right = left.cut(second);
    assert_eq!(left.range_tuple().1.to_bits(), second.to_bits());
    assert_eq!(right.range_tuple().0.to_bits(), second.to_bits());
    for i in 0..=32 {
        let t = a + (b - a) * f64::from(i) / 32.0;
        let split = if t <= second { &left } else { &right };
        assert!(
            curve.evaluate(t).distance(split.evaluate(t)) < 1.0e-8,
            "{t:?}"
        );
    }
}
