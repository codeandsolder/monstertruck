use monstertruck_geometry::prelude::{rbf_surface::RadiusFunction, *};
use std::f64::consts::PI;

#[test]
fn contact_circle_as_curve() {
    let line = Line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let (plane0, plane1) = (Plane::xy(), Plane::zx());
    let (r, t) = (0.5, 0.75);
    let rfs = RollingBallFilletSurface::new(line, plane0, plane1, r);
    let cc = (rfs.contact_circle(t)).unwrap();

    assert_near!(cc.evaluate(0.0), cc.contact_point0().point);
    assert_near!(cc.evaluate(1.0), cc.contact_point1().point);

    let nurbs: NurbsCurve<Vector4> = cc.to_same_geometry();
    assert_near!(cc.evaluate(0.0), nurbs.evaluate(0.0));
    assert_near!(cc.evaluate(0.5), nurbs.evaluate(0.5));
    assert_near!(cc.evaluate(1.0), nurbs.evaluate(1.0));
    assert_near!((nurbs.evaluate(0.2) - cc.center()).magnitude(), r);

    const EPS: f64 = 1.0e-4;
    for i in 1..=9 {
        let s = f64::from(i) / 10.0;
        assert_near!(cc.evaluate(s).distance2(cc.center()), r * r);
        let cal_der = (cc.evaluate(s + EPS) - cc.evaluate(s - EPS)) / (2.0 * EPS);
        assert!((cal_der - cc.derivative(s)).magnitude() < EPS);
        let cal_der2 = (cc.derivative(s + EPS) - cc.derivative(s - EPS)) / (2.0 * EPS);
        assert!((cal_der2 - cc.derivative_2(s)).magnitude() < EPS);
    }
}

#[test]
fn fillet_between_two_spheres() {
    let sphere0 = Sphere::new(Point3::new(0.0, 0.0, 1.0), 2.0);
    let sphere1 = Sphere::new(Point3::new(0.0, 0.0, -1.0), 2.0);
    let edge_circle = Processor::with_transform(
        UnitCircle::<Point3>::new(),
        Matrix4::from_scale(f64::sqrt(3.0)),
    );

    #[derive(Clone, Copy, Debug)]
    struct Radius;
    impl RadiusFunction for Radius {
        fn derivative_n(&self, n: usize, t: f64) -> f64 {
            let o = if n == 0 { 1.0 } else { 0.0 };
            let x = match n % 4 {
                0 => f64::cos(t),
                1 => -f64::sin(t),
                2 => -f64::cos(t),
                _ => f64::sin(t),
            };
            0.2f64.mul_add(x, o)
        }
    }

    let fillet = RollingBallFilletSurface::new(edge_circle, sphere0, sphere1, Radius);
    let cp_curve0 = fillet.contact_curve0();
    let cp_curve1 = fillet.contact_curve1();

    let uc = UnitCircle::<Point3>::new();
    const N: usize = 20;
    for i in 0..=N {
        let t = 2.0 * PI * i as f64 / N as f64;
        let cc = (fillet.contact_circle(t)).unwrap();

        let r = Radius.evaluate(t);
        let center_radius = (r + 2.0).mul_add(r + 2.0, -1.0).sqrt();
        assert_near!(cc.center(), center_radius * uc.evaluate(t));

        let contact_radius = 2.0 / (2.0 + r) * center_radius;
        let contact_z = r / (2.0 + r);
        let contact_point0 = contact_radius * uc.evaluate(t) + contact_z * Vector3::unit_z();
        assert_near!(cp_curve0.evaluate(t), contact_point0);
        let contact_point1 = contact_radius * uc.evaluate(t) - contact_z * Vector3::unit_z();
        assert_near!(cp_curve1.evaluate(t), contact_point1);

        let t0 = (cp_curve0.search_parameter(cc.contact_point0().point, None, 10)).unwrap();
        assert_near!(cp_curve0.evaluate(t0), cc.contact_point0().point);
        let t1 = (cp_curve1.search_parameter(cc.contact_point1().point, None, 10)).unwrap();
        assert_near!(cp_curve1.evaluate(t1), cc.contact_point1().point);

        for j in 0..=N {
            let (u, v) = (j as f64 / N as f64, t);

            let eps = 1.0e-4;

            let ders = fillet.derivatives(2, u, v);

            let ders_plus = fillet.derivatives(1, u + eps, v);
            let ders_minus = fillet.derivatives(1, u - eps, v);
            let uder_approx = (ders_plus[0][0] - ders_minus[0][0]) / (2.0 * eps);
            assert!((ders[1][0] - uder_approx).magnitude() < eps);
            let uvder_approx = (ders_plus[0][1] - ders_minus[0][1]) / (2.0 * eps);
            assert!((ders[1][1] - uvder_approx).magnitude() < eps);
            let uuder_approx = (ders_plus[1][0] - ders_minus[1][0]) / (2.0 * eps);
            assert!((ders[2][0] - uuder_approx).magnitude() < eps);

            let ders_plus = fillet.derivatives(1, u, v + eps);
            let ders_minus = fillet.derivatives(1, u, v - eps);
            let vder_approx = (ders_plus[0][0] - ders_minus[0][0]) / (2.0 * eps);
            assert!((ders[0][1] - vder_approx).magnitude() < eps);
            let uvder_approx = (ders_plus[1][0] - ders_minus[1][0]) / (2.0 * eps);
            assert!((ders[1][1] - uvder_approx).magnitude() < eps);
            let vvder_approx = (ders_plus[0][1] - ders_minus[0][1]) / (2.0 * eps);
            assert!((ders[0][2] - vvder_approx).magnitude() < eps);

            let p = cc.evaluate(u);
            let (u0, v0) = fillet
                .search_parameter(p, None, 10)
                .unwrap_or_else(|| panic!("{:?}", (u, v)));
            assert_near!(fillet.evaluate(u0, v0), p);
        }
    }
}

#[test]
fn rolling_ball_fillet_supports_third_order_mixed_derivatives() {
    let sphere0 = Sphere::new(Point3::new(0.0, 0.0, 1.0), 2.0);
    let sphere1 = Sphere::new(Point3::new(0.0, 0.0, -1.0), 2.0);
    let edge_circle = Processor::with_transform(
        UnitCircle::<Point3>::new(),
        Matrix4::from_scale(f64::sqrt(3.0)),
    );

    #[derive(Clone, Copy, Debug)]
    struct Radius;
    impl RadiusFunction for Radius {
        fn derivative_n(&self, n: usize, t: f64) -> f64 {
            let constant = if n == 0 { 1.0 } else { 0.0 };
            let oscillation = match n % 4 {
                0 => f64::cos(t),
                1 => -f64::sin(t),
                2 => -f64::cos(t),
                _ => f64::sin(t),
            };
            0.2f64.mul_add(oscillation, constant)
        }
    }

    let fillet = RollingBallFilletSurface::new(edge_circle, sphere0, sphere1, Radius);
    let (u, v) = (0.37, 1.1);
    let derivatives = fillet.derivatives(3, u, v);
    let epsilon = 1.0e-4;

    let plus_v = fillet.derivatives(2, u, v + epsilon);
    let minus_v = fillet.derivatives(2, u, v - epsilon);
    for (m, n) in [(0, 3), (1, 2), (2, 1)] {
        let approximation = (plus_v[m][n - 1] - minus_v[m][n - 1]) / (2.0 * epsilon);
        assert!(
            (derivatives[m][n] - approximation).magnitude() < 2.0e-3,
            "third derivative ({m}, {n}) differs: exact={:?}, finite_difference={approximation:?}",
            derivatives[m][n],
        );
    }

    let plus_u = fillet.derivatives(2, u + epsilon, v);
    let minus_u = fillet.derivatives(2, u - epsilon, v);
    let u_approximation = (plus_u[2][0] - minus_u[2][0]) / (2.0 * epsilon);
    assert!(
        (derivatives[3][0] - u_approximation).magnitude() < 2.0e-3,
        "third derivative (3, 0) differs: exact={:?}, finite_difference={u_approximation:?}",
        derivatives[3][0],
    );
}
