use monstertruck_geometry::prelude::*;
use proptest::prelude::*;
use std::f64::consts::PI;
use std::ops::RangeBounds;

proptest! {
    #[test]
    fn surface(
        small_radius in 1f64..=5f64,
        radius_diff in 0.1f64..=5f64,
        (u, v) in (0f64..=2.0 * PI, 0f64..=2.0 * PI),
        deform in -0.4f64..=0.4
    ) {
        const EPS: f64 = 1.0e-2;
        let large_radius = 2.0f64.mul_add(small_radius, radius_diff);
        let torus = Torus::new_unchecked(Point3::origin(), large_radius, small_radius);

        let p = torus.evaluate(u, v);
        let r = large_radius * Point3::new(f64::cos(u), f64::sin(u), 0.0);
        prop_assert_near!(p.distance(r), small_radius);
        prop_assert_near!(p.z, small_radius * f64::sin(v));

        let uder0 = torus.derivative_u(u, v);
        let uder1 = (torus.evaluate(u + EPS, v) - torus.evaluate(u - EPS, v)) / (2.0 * EPS);
        prop_assert!((uder0 - uder1).magnitude() < EPS, "{:?} {:?}", uder0, uder1);

        let vder0 = torus.derivative_v(u, v);
        let vder1 = (torus.evaluate(u, v + EPS) - torus.evaluate(u, v - EPS)) / (2.0 * EPS);
        prop_assert!((vder0 - vder1).magnitude() < EPS, "{:?} {:?}", vder0, vder1);

        let uuder0 = torus.derivative_uu(u, v);
        let uuder1 = (torus.derivative_u(u + EPS, v) - torus.derivative_u(u - EPS, v)) / (2.0 * EPS);
        prop_assert!(
            (uuder0 - uuder1).magnitude() < EPS,
            "{:?} {:?}",
            uuder0,
            uuder1
        );

        let uvder0 = torus.derivative_uv(u, v);
        let uvder1 = (torus.derivative_v(u + EPS, v) - torus.derivative_v(u - EPS, v)) / (2.0 * EPS);
        prop_assert!(
            (uvder0 - uvder1).magnitude() < EPS,
            "{:?} {:?}",
            uvder0,
            uvder1
        );

        let vvder0 = torus.derivative_vv(u, v);
        let vvder1 = (torus.derivative_v(u, v + EPS) - torus.derivative_v(u, v - EPS)) / (2.0 * EPS);
        prop_assert!(
            (vvder0 - vvder1).magnitude() < EPS,
            "{:?} {:?}",
            vvder0,
            vvder1
        );

        let n0 = torus.normal(u, v);
        let n1 = torus.derivative_u(u, v).cross(torus.derivative_v(u, v)).normalize();
        prop_assert_near!(n0, n1);

        let (u0, v0) = torus.search_parameter(p, None, 1).unwrap();
        let (urange, vrange) = torus.parameter_range();
        prop_assert!(urange.contains(&u0) && vrange.contains(&v0), "{u0}, {v0}");
        prop_assert_near!(torus.evaluate(u0, v0), p);

        let deform = deform * small_radius;
        let q = p + deform * n0;
        let (u0, v0) = torus.search_nearest_parameter(q, None, 1).unwrap();
        let (urange, vrange) = torus.parameter_range();
        prop_assert!(urange.contains(&u0) && vrange.contains(&v0), "{u0}, {v0}");
        prop_assert_near!(torus.evaluate(u0, v0), p);
    }

    #[test]
    fn test_der_mn(
        (u, v) in (0f64..=PI, 0f64..=2.0 * PI),
        (m, n) in (0usize..=4, 0usize..=4),
        center in prop::array::uniform3(-100f64..=100.0),
        large_radius in 6.0f64..=10.0,
        small_radius in 0.1f64..=5.0,
        u_derivate in prop::bool::ANY,
    ) {
        let torus = Torus::new_unchecked(Point3::from(center), large_radius, small_radius);

        const EPS: f64 = 1.0e-4;
        let (der0, der1) = if u_derivate {
            let der0 = torus.derivative_mn(m + 1, n, u, v);
            let der1 = (torus.derivative_mn(m, n, u + EPS, v) - torus.derivative_mn(m, n, u - EPS, v)) / (2.0 * EPS);
            (der0, der1)
        } else {
            let der0 = torus.derivative_mn(m, n + 1, u, v);
            let der1 = (torus.derivative_mn(m, n, u, v + EPS) - torus.derivative_mn(m, n, u, v - EPS)) / (2.0 * EPS);
            (der0, der1)
        };
        prop_assert!((der0 - der1).magnitude() < 0.01 * der0.magnitude());
    }
}

/// Ring, horn, and spindle tori all need a sound inverse map.
///
/// The spindle used to lose roughly a quarter of its own parameter domain
/// because the inverse derived `v` from unsigned cylindrical radius and could
/// not represent the negative signed-radial branch.  The solver now enumerates
/// the meridian branches explicitly; every point produced by the forward map
/// must therefore round-trip, and nearest-point recovery of an on-surface point
/// must land on that point rather than on the other sheet.
#[test]
fn torus_parameter_recovery_is_sound_in_ring_horn_and_spindle_regimes() {
    let cases: [(&str, f64, f64); 5] = [
        ("Rocky_House witness spindle", 0.633_974_596_215_563, 1.0),
        ("NissanGT-R spindle |major|", 20.852_342_325_613_467, 29.5),
        ("ring", 3.0, 1.0),
        ("exact horn (Ai-14R fillet)", 3.0, 3.0),
        (
            "fp near-horn (Pi fillet)",
            0.099_999_999_992_725,
            0.099_999_999_999_987_88,
        ),
    ];
    const N: usize = 24;
    for (label, large_radius, small_radius) in cases {
        let torus = Torus::new_unchecked(Point3::origin(), large_radius, small_radius);
        let tol = 1.0e-9 * (large_radius + small_radius);
        for i in 0..N {
            // Offsets keep the samples off periodic seams and singular axis points.
            let u = 2.0 * PI * (i as f64 + 0.37) / N as f64;
            for j in 0..N {
                let v = 2.0 * PI * (j as f64 + 0.21) / N as f64;
                let point = torus.evaluate(u, v);
                let (eu, ev) = torus
                    .search_parameter(point, Some((u, v)), 100)
                    .unwrap_or_else(|| panic!("{label}: exact inverse refused ({u}, {v})"));
                assert!(
                    torus.evaluate(eu, ev).distance(point) <= tol,
                    "{label}: exact inverse landed elsewhere at ({u}, {v}) -> ({eu}, {ev})",
                );
                let (nu, nv) = torus
                    .search_nearest_parameter(point, Some((u, v)), 100)
                    .unwrap_or_else(|| panic!("{label}: nearest inverse refused ({u}, {v})"));
                assert!(
                    torus.evaluate(nu, nv).distance(point) <= tol,
                    "{label}: nearest inverse landed elsewhere at ({u}, {v}) -> ({nu}, {nv})",
                );
            }
        }
    }
}

/// The analytic nearest-point solver is global even in the spindle regime.
/// A dense meridian sweep is an intentionally independent lower-resolution
/// oracle: the exact solver may beat it, but must never be measurably worse.
#[test]
fn spindle_torus_nearest_parameter_matches_dense_meridian_reference() {
    let torus = Torus::new_unchecked(Point3::new(0.2, -0.4, 0.7), 0.7, 1.3);
    let queries = [
        Point3::new(0.2, -0.4, 0.7),
        Point3::new(0.8, -0.1, 0.2),
        Point3::new(-1.1, 0.5, 1.9),
        Point3::new(2.7, -0.4, 0.7),
        Point3::new(0.2, -0.4, 3.1),
        Point3::new(-0.3, -0.9, 0.65),
    ];
    const SAMPLES: usize = 16_384;
    for point in queries {
        let (u, v) = torus
            .search_nearest_parameter(point, None, 100)
            .expect("spindle nearest-point solver must return a candidate");
        let got = torus.evaluate(u, v).distance2(point);
        let relative = point - torus.center();
        let rho = (relative.x * relative.x + relative.y * relative.y).sqrt();
        let mut sampled = f64::INFINITY;
        for i in 0..SAMPLES {
            let sv = 2.0 * PI * i as f64 / SAMPLES as f64;
            let q = torus.large_radius() + torus.small_radius() * sv.cos();
            let dz = relative.z - torus.small_radius() * sv.sin();
            let dr = rho - q.abs();
            sampled = sampled.min(dr * dr + dz * dz);
        }
        let scale = (torus.large_radius() + torus.small_radius() + relative.magnitude()).max(1.0);
        assert!(
            got <= sampled + (2.0 * PI / SAMPLES as f64 * scale).powi(2),
            "nearest point worse than dense meridian reference: got {got}, sampled {sampled}, point {point:?}",
        );
    }
}

/// **A seam point must be spelled the way the caller asks for** -- spec 012 W1.
///
/// A torus is exactly `2 pi`-periodic in both parameters, so the point on the
/// `u = 0` seam is `u = 0` and `u = 2 pi` equally. `search_parameter` picks by
/// the sign of `y`, and at the seam that sign is the sign of a ZERO: a
/// negative-zero `y` yields exactly `2 pi`, a positive-zero `y` exactly `0`.
/// Both are correct in isolation and only the caller knows which one its
/// parameter LOOP is written in -- so the caller's hint has to decide, and
/// before this it was discarded outright (`_: H`).
///
/// Measured consequence, ap224 with spec 012's analytic-torus routing on: two
/// fillet tori touching one seam from opposite sides both reported their `u`
/// trim extent as the WHOLE ring (`(0.0645, 6.2832)` and `(0, 6.2187)`) instead
/// of `(0, pi)` and `(pi, 2 pi)`. The SSI then traced the right curves against
/// the box cutter's planes and the trim filter discarded every one of them.
/// Ledger class C4.
///
/// Both directions are asserted, because a rule that only ever answers `0`
/// fixes one of those two faces and breaks the other -- which is exactly what
/// the first attempt at this did.
#[test]
fn a_seam_point_takes_the_periodic_spelling_its_hint_asks_for() {
    let torus = Torus::new_unchecked(Point3::new(1.0, -2.0, 3.0), 7.0, 2.0);
    // Straddle the seam so the answers must NOT all collapse to one end.
    for v in [0.0, 0.7, PI, 5.0] {
        let seam = torus.evaluate(0.0, v);
        let below = torus
            .search_parameter(seam, Some((0.05, v)), 100)
            .expect("a point on the surface must be recovered");
        assert!(
            below.0.abs() < 0.1,
            "a hint just above the seam must keep the low spelling, got {below:?}",
        );
        let above = torus
            .search_parameter(seam, Some((2.0f64.mul_add(PI, -0.05), v)), 100)
            .expect("a point on the surface must be recovered");
        assert!(
            2.0f64.mul_add(-PI, above.0).abs() < 0.1,
            "a hint just below 2 pi must keep the high spelling, got {above:?}",
        );
        // Whichever spelling comes back, it is the SAME 3D point. That is what
        // makes the re-spelling exact rather than a tolerance trade.
        for (u, v) in [below, above] {
            assert!(
                torus.evaluate(u, v).distance(seam) < 1.0e-9,
                "re-spelled ({u}, {v}) must still evaluate to the seam point",
            );
        }
        // The nearest solver is the twin `project_onto_surface_domain`
        // alternates with, so it has to agree.
        let nearest_below = torus
            .search_nearest_parameter(seam, Some((0.05, v)), 100)
            .expect("nearest must answer for an on-surface point");
        let nearest_above = torus
            .search_nearest_parameter(seam, Some((2.0f64.mul_add(PI, -0.05), v)), 100)
            .expect("nearest must answer for an on-surface point");
        assert!(
            nearest_below.0.abs() < 0.1 && 2.0f64.mul_add(-PI, nearest_above.0).abs() < 0.1,
            "the nearest solver must disambiguate the same way: {nearest_below:?} / \
             {nearest_above:?}",
        );
    }
}

/// Without a hint, parameter recovery must return a member of the declared
/// periodic domain that reconstructs the original on-surface point.
///
/// Seam points have multiple valid parameter representations, so no particular
/// branch's numeric answer is part of the contract.
#[test]
fn unhinted_search_recovers_point_in_declared_period() {
    let torus = Torus::new_unchecked(Point3::origin(), 7.0, 2.0);
    for u in [0.0, 0.3, 1.9, PI, 4.4, 2.0f64.mul_add(PI, -0.3)] {
        for v in [0.0, 0.4, 2.2, PI, 5.9] {
            let point = torus.evaluate(u, v);
            let (recovered_u, recovered_v) = torus
                .search_parameter(point, None, 100)
                .expect("an on-surface point must be recovered");
            // Semantic guarantee: stay in the declared periodic domain and
            // reproduce the original point.
            assert!(
                (0.0..=2.0 * PI).contains(&recovered_u) && (0.0..=2.0 * PI).contains(&recovered_v),
                "unhinted ({u}, {v}) left the declared range: \
                 ({recovered_u}, {recovered_v})",
            );
            assert!(
                torus.evaluate(recovered_u, recovered_v).distance(point) < 1.0e-9,
                "unhinted ({u}, {v}) did not reproduce its point",
            );
        }
    }
}
