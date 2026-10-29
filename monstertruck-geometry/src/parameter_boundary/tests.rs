//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;
use std::f64::consts::FRAC_PI_4;

#[test]
fn plane_line_pcurve_converts_to_bspline_curve() {
    let plane = Plane::xy();
    let line = Line(Point2::new(0.25, 0.5), Point2::new(1.25, -0.5));
    let curve = (BsplineCurve::<Point3>::try_from(ParameterCurve::new(line, plane)))
        .expect("plane line parameter curve should convert to a `BsplineCurve`.");

    assert!(
        curve
            .evaluate(0.0)
            .near(&plane.evaluate(line.0.x, line.0.y))
    );
    assert!(
        curve
            .evaluate(1.0)
            .near(&plane.evaluate(line.1.x, line.1.y))
    );
    assert_eq!(curve.control_points().len(), 2);
}

#[test]
fn revolution_surface_exact_boundary_recovers_internal_iso_angle_line() {
    let profile = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 0.0, 1.0)],
    );
    let surface = RevolutionSurface::by_revolution(profile, Point3::origin(), Vector3::unit_z());
    let angle = FRAC_PI_4;
    let curve = Line(
        Point3::new(angle.cos(), angle.sin(), 0.25),
        Point3::new(angle.cos(), angle.sin(), 0.75),
    );

    let boundary = (curve.exact_parameter_boundary_2d(&surface))
        .expect("iso-angle line must have an exact surface parameter boundary");
    let line = boundary.curve();

    assert!((line.0.x - 0.25).abs() <= TOLERANCE);
    assert!((line.1.x - 0.75).abs() <= TOLERANCE);
    assert!((line.0.y - angle).abs() <= TOLERANCE);
    assert!((line.1.y - angle).abs() <= TOLERANCE);
}

#[test]
fn revolution_line_profile_exact_boundary_extends_beyond_placeholder_range() {
    let profile = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0));
    let surface = RevolutionSurface::by_revolution(profile, Point3::origin(), Vector3::unit_z());
    let angle = FRAC_PI_4;
    let curve = Line(surface.subs(-3.0, angle), surface.subs(4.0, angle));

    let boundary = curve
        .exact_parameter_boundary_2d(&surface)
        .expect("analytic line profile must extend exactly beyond its placeholder range");
    let line = boundary.curve();
    assert!((line.0.x + 3.0).abs() <= TOLERANCE);
    assert!((line.1.x - 4.0).abs() <= TOLERANCE);
    assert!((line.0.y - angle).abs() <= TOLERANCE);
    assert!((line.1.y - angle).abs() <= TOLERANCE);
}

#[test]
fn revolution_surface_exact_boundary_recovers_iso_angle_nurbs_with_fixed_endpoint() {
    let profile = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![Point3::origin(), Point3::new(1.0, 0.0, 1.0)],
    );
    let surface = RevolutionSurface::by_revolution(profile, Point3::origin(), Vector3::unit_z());
    let angle = FRAC_PI_4;
    let curve = NurbsCurve::<Vector4>::from(BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(1),
        vec![surface.evaluate(0.0, angle), surface.evaluate(1.0, angle)],
    ));

    let boundary = (curve.exact_parameter_boundary_2d(&surface))
        .expect("iso-angle NURBS with a fixed endpoint must have an exact boundary");
    let line = boundary.curve();

    assert!(line.0.x.abs() <= TOLERANCE);
    assert!((line.1.x - 1.0).abs() <= TOLERANCE);
    assert!((line.0.y - angle).abs() <= TOLERANCE);
    assert!((line.1.y - angle).abs() <= TOLERANCE);
}

#[test]
fn exact_bspline_boundary_rejects_five_checkpoint_false_positive() {
    let surface = BsplineSurface::new(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
        vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
    );
    // z(t) = t(t-1/4)(t-1/2)(t-3/4)(t-1).  Its degree-5 Bernstein
    // coefficients below make the curve coincide with the y=0 surface edge at
    // exactly the five checkpoints used by the removed sampled "exact" test,
    // while the nonzero interior control values prove it is not that edge.
    let z = [
        0.0,
        3.0 / 160.0,
        -13.0 / 320.0,
        13.0 / 320.0,
        -3.0 / 160.0,
        0.0,
    ];
    let curve = BsplineCurve::new(
        KnotVector::bezier_knot(5),
        z.into_iter()
            .enumerate()
            .map(|(index, z)| Point3::new(index as f64 / 5.0, 0.0, z))
            .collect(),
    );

    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        assert!(curve.subs(t).z.abs() <= TOLERANCE);
    }
    assert!(curve.subs(0.125).z.abs() > TOLERANCE);
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn nurbs_extrusion_exact_boundary_recovers_generator_line() {
    let vector = Vector3::new(0.0, 0.0, 2.0);
    let control_points = vec![
        vec![
            Vector4::new(1.0, 0.0, 0.0, 1.0),
            Vector4::new(1.0, 0.0, 2.0, 1.0),
        ],
        vec![
            Vector4::new(1.0, 1.0, 0.0, 1.0),
            Vector4::new(1.0, 1.0, 2.0, 1.0),
        ],
        vec![
            Vector4::new(0.0, 1.0, 0.0, 1.0),
            Vector4::new(0.0, 1.0, 2.0, 1.0),
        ],
    ];
    let surface = NurbsSurface::new(BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
        control_points,
    ));
    let u = 0.5;
    let v0 = 0.25;
    let v1 = 0.75;
    let curve = Line(surface.evaluate(u, v0), surface.evaluate(u, v1));

    assert!((curve.1 - curve.0).near(&Vector3::new(0.0, 0.0, vector.z * (v1 - v0))));

    let boundary = (curve.exact_parameter_boundary_2d(&surface))
        .expect("extrusion generator line must have an exact surface parameter boundary");
    let line = boundary.curve();

    assert!(line.0.near(&Point2::new(u, v0)));
    assert!(line.1.near(&Point2::new(u, v1)));
}
