use monstertruck_geometry::prelude::*;

#[test]
fn canonical_periodicity_is_directionally_correct() {
    let torus = Torus::new_unchecked(Point3::origin(), 3.0, 1.0);
    assert_eq!(torus.period_u(), Some(2.0 * std::f64::consts::PI));
    assert_eq!(torus.period_v(), Some(2.0 * std::f64::consts::PI));

    let sphere = Sphere::new(Point3::origin(), 2.0);
    assert_eq!(sphere.period_u(), None);
    assert_eq!(sphere.period_v(), Some(2.0 * std::f64::consts::PI));
}

#[test]
fn canonical_plane_axes_are_directionally_correct() {
    let plane = Plane::new(
        Point3::new(1.0, 2.0, 3.0),
        Point3::new(5.0, 2.0, 3.0),
        Point3::new(1.0, 9.0, 3.0),
    );
    assert_eq!(plane.axis_u(), Vector3::new(4.0, 0.0, 0.0));
    assert_eq!(plane.axis_v(), Vector3::new(0.0, 7.0, 0.0));
}

#[test]
fn canonical_sectional_curves_follow_parameter_direction() {
    let surface = BsplineSurface::new_unchecked(
        (KnotVector::bezier_knot(1), KnotVector::bezier_knot(2)),
        vec![
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 1.0),
                Point3::new(2.0, 0.0, 2.0),
            ],
            vec![
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 1.0, 1.0),
                Point3::new(2.0, 1.0, 2.0),
            ],
        ],
    );
    let along_u = surface.curve_u(1);
    assert_eq!(along_u.knot_vector(), &KnotVector::bezier_knot(1));
    assert_eq!(along_u.control_points().len(), 2);
    let along_v = surface.curve_v(1);
    assert_eq!(along_v.knot_vector(), &KnotVector::bezier_knot(2));
    assert_eq!(along_v.control_points().len(), 3);
}
