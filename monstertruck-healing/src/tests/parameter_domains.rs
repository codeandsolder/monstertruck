use super::*;
use monstertruck_geometry::prelude::{BsplineSurface, KnotVector, Matrix4, Plane, Processor};

#[test]
fn healing_keeps_planar_trim_coordinates_outside_synthetic_unit_square() {
    let plane = Plane::xy();
    assert_eq!(
        plane.try_range_tuple(),
        (Some((0.0, 1.0)), Some((0.0, 1.0)))
    );
    assert_eq!(plane.parameter_bounds_are_intrinsic(), (false, false));
    let sample = Point2::new(3.9, -2.1);
    assert_eq!(
        normalize_surface_parameter(&plane, sample, None),
        Some(sample)
    );
    let previous = Point2::new(3.8, -2.2);
    assert_eq!(
        normalize_surface_parameter(&plane, sample, Some(previous)),
        Some(sample),
    );
    let mut processor: Processor<_, Matrix4> = Processor::new(plane);
    assert_eq!(
        normalize_surface_parameter(&processor, sample, None),
        Some(sample)
    );
    processor.invert();
    assert_eq!(
        normalize_surface_parameter(&processor, sample, None),
        Some(sample)
    );
}

#[test]
fn healing_snaps_only_roundoff_at_genuine_bounded_knot_domain() {
    let knots = KnotVector::bezier_knot(1);
    let surface = BsplineSurface::new_unchecked(
        (knots.clone(), knots),
        vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
    );
    assert_eq!(surface.parameter_bounds_are_intrinsic(), (true, true));
    let near = Point2::new(1.0 + 0.5e-6, -0.5e-6);
    assert_eq!(
        normalize_surface_parameter(&surface, near, None),
        Some(Point2::new(1.0, 0.0))
    );
    assert_eq!(
        normalize_surface_parameter(&surface, near, Some(Point2::new(0.8, 0.1))),
        Some(Point2::new(1.0, 0.0))
    );
    assert_eq!(
        normalize_surface_parameter(&surface, Point2::new(1.02, 0.5), None),
        None,
        "A geometric 2% displacement must never be disguised as endpoint roundoff"
    );
    assert_eq!(
        normalize_surface_parameter(&surface, Point2::new(0.5, 0.5), None),
        Some(Point2::new(0.5, 0.5))
    );
}

#[test]
fn healing_native_small_period_is_not_a_zero_length_chart() {
    let period = 1e-9;
    let range = Some((0.0, period));
    let wrapped = normalize_parameter_axis(1.5 * period, None, Some(period), range)
        .expect("positive finite period");
    assert!((wrapped - 0.5 * period).abs() <= 4.0 * f64::EPSILON * period);
    let shifted =
        normalize_parameter_axis(1000.5 * period, Some(0.5 * period), Some(period), range)
            .expect("periodic continuity");
    assert!((shifted - 0.5 * period).abs() <= 2048.0 * f64::EPSILON * period);
    assert_eq!(
        normalize_parameter_axis(0.5, None, Some(0.0), Some((0.0, 1.0))),
        None
    );
    assert_eq!(
        normalize_parameter_axis(0.5, None, Some(f64::NAN), Some((0.0, 1.0))),
        None
    );
    assert_eq!(
        normalize_parameter_axis(0.5, Some(f64::INFINITY), None, None),
        None
    );
    assert_eq!(
        normalize_parameter_axis(0.5, None, None, Some((1.0, 0.0))),
        None
    );
    assert_eq!(
        normalize_parameter_axis(1.0, None, None, Some((1.0, 1.0))),
        Some(1.0)
    );
    assert_eq!(
        normalize_parameter_axis(1.0 + 1e-8, None, None, Some((1.0, 1.0))),
        None
    );
}

#[test]
fn healing_clamp_scales_with_the_actual_parameter_span() {
    for scale in [1.0e-9, 1.0, 1.0e9] {
        let range = (0.0, scale);
        let policy = DEFAULT_NUMERICAL_TOLERANCES;
        let admissible = scale * policy.healing_projection_parameter_relative * 0.5;
        let inadmissible = scale * policy.healing_projection_parameter_relative * 10.0;
        assert_eq!(
            normalize_intrinsic_parameter(scale + admissible, range),
            Some(scale)
        );
        assert_eq!(
            normalize_intrinsic_parameter(scale + inadmissible, range),
            None
        );
    }
}

#[test]
fn healing_trim_cut_accepts_tiny_genuine_native_parameter_intervals() {
    let limit = 1.0e-9;
    let curve = BsplineCurve::new_unchecked(
        KnotVector::from(vec![0.0, 0.0, limit, limit]),
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
    );
    let first = curve.evaluate(0.25 * limit);
    let last = curve.evaluate(0.75 * limit);
    let verts = [first, last];
    let segment = CompressedEdge {
        vertices: (0, 1),
        curve: Line(first, last),
    };
    let edge_use = CompressedEdgeUse {
        index: 0,
        orientation: true,
        trim_curve: Some(curve),
    };
    let split = oriented_trim_segment(&edge_use, true, &segment, &verts, 1.0e-12)
        .expect("small native parameter interval still has an interior");
    let (lo, hi) = split.range_tuple();
    assert!((lo - 0.25 * limit).abs() <= 1.0e-17);
    assert!((hi - 0.75 * limit).abs() <= 1.0e-17);
}

#[test]
fn native_trim_split_tolerance_is_relative_to_native_parameter_span() {
    let policy = DEFAULT_NUMERICAL_TOLERANCES;
    for width in [1.0e-10, 1.0e-7, 1.0, 1.0e8] {
        let margin = trim_parameter_slack((0.0, width)).unwrap();
        assert!(margin >= width * policy.trim_parameter_relative);
        assert!(margin < width * policy.trim_parameter_relative * 1.0000001);
    }
    assert!(trim_parameter_slack((2.0, 1.0)).is_none());
    assert!(trim_parameter_slack((1.0, 1.0)).is_none());
    assert!(trim_parameter_slack((0.0, f64::INFINITY)).is_none());
    let large_origin = trim_parameter_slack((1.0e8, 1.0e8 + 1.0)).unwrap();
    let zero_origin = trim_parameter_slack((0.0, 1.0)).unwrap();
    assert!(large_origin > zero_origin);
}
