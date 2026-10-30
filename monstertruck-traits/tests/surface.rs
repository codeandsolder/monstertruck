#![allow(
    clippy::suboptimal_flops,
    reason = "reference polynomial formulas stay expanded so the oracle is readable and does not change rounding via FMA"
)]
#![allow(
    clippy::many_single_char_names,
    reason = "surface tests use canonical u/v parameter and compact algebraic notation"
)]
#![allow(
    clippy::similar_names,
    reason = "paired u/v derivative and division names are canonical surface notation"
)]
#![cfg(feature = "polynomial")]

use monstertruck_core::assert_near2;
use monstertruck_core::cgmath64::{
    EuclideanSpace, InnerSpace, MetricSpace, Point2, Point3, Vector2, Vector3,
};
use monstertruck_core::tolerance::{Origin, Tolerance};
use monstertruck_traits::algo::surface;
use monstertruck_traits::polynomial::{PolynomialCurve, PolynomialSurface};
use monstertruck_traits::{ParametricCurve, ParametricSurface, ParametricSurface3D};
use rand::{RngExt, SeedableRng, rngs::StdRng};

/// Deterministic generator for the randomised search trials below.
///
/// These tests explore a numerical search over random inputs and assert that it
/// succeeds on most of them. With an unseeded generator that is a coin flip on
/// every run: the suite fails on its own schedule, for reasons no diff explains,
/// and a genuine regression is indistinguishable from bad luck. Seeding fixes
/// the sample, so the trial set is identical on every machine and every run --
/// the coverage is the same, the verdict is reproducible.
///
/// Changing the seed re-rolls the sample and is a real change to what is
/// exercised. The minimum success rates below are the explicit robustness
/// contract; they are not copied from the current solver's exact failure count.
fn trial_rng() -> StdRng { StdRng::seed_from_u64(0x5445_5354_5F52_4E47) }

#[test]
fn polysurface() {
    let coef0 = vec![
        Vector3::new(1.0, 0.0, 1.0),
        Vector3::new(3.0, 1.0, 0.0),
        Vector3::new(2.0, 0.0, 0.0),
    ];
    let coef1 = vec![
        Vector3::new(2.0, 1.0, 0.0),
        Vector3::new(-6.0, 0.0, 1.0),
        Vector3::new(4.0, 0.0, 0.0),
    ];
    let curve0 = PolynomialCurve::<Point3>(coef0);
    let curve1 = PolynomialCurve::<Point3>(coef1);
    let poly = PolynomialSurface::by_tensor(curve0, curve1);
    for i in 0..5 {
        let u = f64::from(i);
        for j in 0..5 {
            let v = f64::from(j);
            assert_eq!(
                poly.evaluate(u, v)[0],
                (2.0 * u * u + 3.0 * u + 1.0) * (4.0 * v * v - 6.0 * v + 2.0)
            );
            assert_eq!(
                poly.derivative_u(u, v)[0],
                (4.0 * u + 3.0) * (4.0 * v * v - 6.0 * v + 2.0)
            );
            assert_eq!(
                poly.derivative_v(u, v)[0],
                (2.0 * u * u + 3.0 * u + 1.0) * (8.0 * v - 6.0)
            );
            assert_eq!(
                poly.derivative_uu(u, v)[0],
                4.0 * (4.0 * v * v - 6.0 * v + 2.0)
            );
            assert_eq!(
                poly.derivative_uv(u, v)[0],
                (4.0 * u + 3.0) * (8.0 * v - 6.0)
            );
            assert_eq!(
                poly.derivative_vv(u, v)[0],
                (2.0 * u * u + 3.0 * u + 1.0) * 8.0
            );
            assert!(poly.normal(u, v).dot(poly.derivative_u(u, v)).so_small());
            assert!(poly.normal(u, v).dot(poly.derivative_v(u, v)).so_small());

            let eps = 1.0e-4;
            let normal_uder_approx =
                (poly.normal(u + eps, v) - poly.normal(u - eps, v)) / (2.0 * eps);
            let normal_uder = poly.normal_uder(u, v);
            assert!((normal_uder - normal_uder_approx).magnitude() < eps);
            let normal_vder_approx =
                (poly.normal(u, v + eps) - poly.normal(u, v - eps)) / (2.0 * eps);
            let normal_vder = poly.normal_vder(u, v);
            assert!((normal_vder - normal_vder_approx).magnitude() < eps);
        }
    }
}

#[test]
fn polysurface_presearch() {
    let coef0 = vec![Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 0.0, 0.0)];
    let coef1 = vec![Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 1.0, 0.0)];
    let curve0 = PolynomialCurve::<Point3>(coef0);
    let curve1 = PolynomialCurve::<Point3>(coef1);
    let poly = PolynomialSurface::by_tensor(curve0, curve1);
    let pt = Point3::new(0.2, 0.3, 1.0);
    let (u, v) = surface::presearch(&poly, pt, ((0.0, 1.0), (0.0, 1.0)), 100);
    assert_eq!(u, 0.2);
    assert_eq!(v, 0.3);
}

fn exec_polysurface_snp_on_surface(rng: &mut StdRng) -> bool {
    let coef0 = vec![
        Vector3::new(0.0, 1.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(1.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
    ];
    let coef1 = vec![
        Vector3::new(1.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 1.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
    ];
    let curve0 = PolynomialCurve::<Point3>(coef0);
    let curve1 = PolynomialCurve::<Point3>(coef1);
    let poly = PolynomialSurface::by_tensor(curve0, curve1);
    let u = 10.0 * rng.random::<f64>() - 5.0;
    let v = 10.0 * rng.random::<f64>() - 5.0;
    let pt = poly.evaluate(u, v);
    let u0 = u + 0.2 * rng.random::<f64>() - 0.1;
    let v0 = v + 0.2 * rng.random::<f64>() - 0.1;
    let answer = (u, v);
    let hint = (u0, v0);
    surface::search_nearest_parameter(&poly, pt, hint, 100).map_or_else(
        || {
            eprintln!("not converge\npolynomial: {poly:?}\nanswer: {answer:?}\nhint: {hint:?}");
            false
        },
        |res| {
            if poly.evaluate(res.0, res.1).near(&pt) {
                true
            } else {
                eprintln!(
                    "wrong answer\npolynomial: {poly:?}\nanswer: {answer:?}\nhint: {hint:?}\nresult: {res:?}",
                );
                false
            }
        },
    )
}

#[test]
fn polysurface_snp_on_surface() {
    let mut rng = trial_rng();
    let trials = 200;
    let successes = (0..trials)
        .filter(|_| exec_polysurface_snp_on_surface(&mut rng))
        .count();
    assert!(
        successes >= 180,
        "nearest-parameter search must reconstruct the on-surface point in at least 90% of the deterministic trials: {successes}/{trials} succeeded",
    );
}

fn exec_polysurface_sp_on_surface(rng: &mut StdRng) -> bool {
    let coef0 = vec![
        Vector3::new(0.0, 1.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(1.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
    ];
    let coef1 = vec![
        Vector3::new(1.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 1.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
        Vector3::new(0.0, 0.0, 3.0 * rng.random::<f64>() - 1.5),
    ];
    let curve0 = PolynomialCurve::<Point3>(coef0);
    let curve1 = PolynomialCurve::<Point3>(coef1);
    let poly = PolynomialSurface::by_tensor(curve0, curve1);
    let u = 10.0 * rng.random::<f64>() - 5.0;
    let v = 10.0 * rng.random::<f64>() - 5.0;
    let pt = poly.evaluate(u, v);
    let u0 = u + 2.0 * rng.random::<f64>() - 1.0;
    let v0 = v + 2.0 * rng.random::<f64>() - 1.0;
    let answer = (u, v);
    let hint = (u0, v0);
    surface::search_parameter(&poly, pt, hint, 100).map_or_else(
        || {
            eprintln!("not converge\npolynomial: {poly:?}\nanswer: {answer:?}\nhint: {hint:?}");
            false
        },
        |res| {
            if poly.evaluate(res.0, res.1).near(&pt) {
                true
            } else {
                eprintln!(
                    "wrong answer\npolynomial: {poly:?}\nanswer: {answer:?}\nhint: {hint:?}\nresult: {res:?}",
                );
                false
            }
        },
    )
}

#[test]
fn polysurface_sp_on_surface() {
    let mut rng = trial_rng();
    let trials = 200;
    let successes = (0..trials)
        .filter(|_| exec_polysurface_sp_on_surface(&mut rng))
        .count();
    assert!(
        successes >= 180,
        "parameter search must reconstruct the on-surface point in at least 90% of the deterministic trials: {successes}/{trials} succeeded",
    );
}

fn exec_polysurface_intersection_point(rng: &mut StdRng) -> bool {
    let (a, b) = (rng.random::<f64>(), rng.random::<f64>());
    let coef0 = vec![
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let coef1 = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let curve0 = PolynomialCurve::<Point3>(coef0);
    let curve1 = PolynomialCurve::<Point3>(coef1);
    let surface = PolynomialSurface::by_tensor(curve0, curve1);
    let coef = vec![Vector3::new(a, b, 0.0), Vector3::new(0.0, 0.0, 1.0)];
    let curve = PolynomialCurve::<Point3>(coef);

    match surface::search_intersection_parameter(&surface, (0.5, 0.5), &curve, 0.0, 100) {
        Some(((x, y), z)) => {
            let (p, q) = (surface.evaluate(x, y), curve.evaluate(z));
            p.near(&q) && x.near(&a) && y.near(&b) && p.z.near(&z)
        }
        None => false,
    }
}

#[test]
fn polysurface_intersection_point() {
    let mut rng = trial_rng();
    let count = (0..10)
        .filter(|_| exec_polysurface_intersection_point(&mut rng))
        .count();
    assert_eq!(
        count,
        10,
        "intersection_point regressed: {} of 10 failed",
        10 - count
    );
}

fn exec_polysurface_division(rng: &mut StdRng) -> bool {
    let coef0 = vec![
        Vector3::new(0.0, 1.0, 10.0 * rng.random::<f64>() - 5.0),
        Vector3::new(1.0, 0.0, 10.0 * rng.random::<f64>() - 5.0),
        Vector3::new(0.0, 0.0, 10.0 * rng.random::<f64>() - 5.0),
    ];
    let coef1 = vec![
        Vector3::new(1.0, 0.0, 10.0 * rng.random::<f64>() - 5.0),
        Vector3::new(0.0, 1.0, 10.0 * rng.random::<f64>() - 5.0),
        Vector3::new(0.0, 0.0, 10.0 * rng.random::<f64>() - 5.0),
    ];
    let curve0 = PolynomialCurve::<Point3>(coef0);
    let curve1 = PolynomialCurve::<Point3>(coef1);
    let poly = PolynomialSurface::by_tensor(curve0, curve1);
    let (udiv, vdiv) = surface::parameter_division(&poly, ((-1.0, 1.0), (-1.0, 1.0)), 0.1);
    for (i, u) in udiv
        .windows(2)
        .flat_map(move |u| (1..3).map(move |i| (i, u)))
    {
        for (j, v) in vdiv
            .windows(2)
            .flat_map(move |v| (1..3).map(move |j| (j, v)))
        {
            let p = f64::from(i) / 3.0;
            let q = f64::from(j) / 3.0;
            let pt0 = Point3::from_vec(
                poly.evaluate(u[0], v[0]).to_vec() * (1.0 - p) * (1.0 - q)
                    + poly.evaluate(u[1], v[0]).to_vec() * p * (1.0 - q)
                    + poly.evaluate(u[0], v[1]).to_vec() * (1.0 - p) * q
                    + poly.evaluate(u[1], v[1]).to_vec() * p * q,
            );
            let u = u[0] * (1.0 - p) + u[1] * p;
            let v = v[0] * (1.0 - q) + v[1] * q;
            let pt1 = poly.evaluate(u, v);
            if pt0.distance(pt1) > 0.2 {
                return false;
            }
        }
    }
    true
}

#[test]
fn polysurface_division() {
    let mut rng = trial_rng();
    let count = (0..10)
        .filter(|_| exec_polysurface_division(&mut rng))
        .count();
    assert_eq!(count, 10, "division regressed: {} of 10 failed", 10 - count);
}

#[test]
fn test_composite() {
    let curve_vec = vec![
        Vector2::new(1.0, 0.0),
        Vector2::new(0.0, 2.0),
        Vector2::new(1.0, 0.0),
    ];
    let curve = PolynomialCurve::<Point2>(curve_vec);

    let surface_vec = vec![
        vec![
            Vector2::new(0.0, 0.0),
            Vector2::new(0.0, 1.0),
            Vector2::new(1.0, 0.0),
        ],
        vec![Vector2::new(0.0, 1.0), Vector2::new(1.0, 0.0)],
    ];
    let surface = PolynomialSurface::<Point2>(surface_vec);

    let res = surface.composite(&curve);
    assert_near2!(res.0[0], Vector2::new(0.0, 1.0));
    assert_near2!(res.0[1], Vector2::new(2.0, 2.0));
    assert_near2!(res.0[2], Vector2::new(4.0, 1.0));
    assert_near2!(res.0[3], Vector2::new(2.0, 0.0));
    res.0.iter().skip(4).for_each(|&p| assert!(p.so_small2()));
}
