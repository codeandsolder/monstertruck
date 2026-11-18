#![expect(
    clippy::similar_names,
    reason = "paired u/v coefficient variables reflect geometric parameter axes"
)]

use monstertruck_geometry::prelude::*;
use monstertruck_traits::polynomial::*;
use proptest::prelude::*;

// The explicit Faà di Bruno derivative and polynomial coefficient composition
// evaluate algebraically identical derivatives through different floating-point
// operation sequences. A fixed absolute kernel-space tolerance does not scale
// to fifth derivatives whose magnitudes can be thousands; compare each
// component with a strict relative budget while retaining the absolute floor.
fn prop_assert_derivative_agreement(
    left: Vector3,
    right: Vector3,
) -> std::result::Result<(), TestCaseError> {
    const RELATIVE_TOLERANCE: f64 = 1.0e-9;
    for axis in 0..3 {
        let lhs = left[axis];
        let rhs = right[axis];
        let tolerance = KERNEL_EPSILON.max(RELATIVE_TOLERANCE * lhs.abs().max(rhs.abs()));
        prop_assert!(
            lhs.is_finite() && rhs.is_finite() && (lhs - rhs).abs() <= tolerance,
            "derivative component {axis} differs: left={lhs:?}, right={rhs:?}, tolerance={tolerance:?}",
        );
    }
    Ok(())
}

fn exec_pcurve_derivation(
    curve_coef: &[[f64; 2]],
    surface_ucoef: &[[f64; 3]],
    surface_vcoef: &[[f64; 3]],
    n: usize,
    t: f64,
) -> std::result::Result<(), TestCaseError> {
    let curve_coef = curve_coef.iter().map(|&p| Vector2::from(p)).collect();
    let curve = PolynomialCurve::<Point2>(curve_coef);

    let surface_ucoef = surface_ucoef.iter().map(|&p| Vector3::from(p)).collect();
    let surface_vcoef = surface_vcoef.iter().map(|&p| Vector3::from(p)).collect();
    let ucurve = PolynomialCurve::<Point3>(surface_ucoef);
    let vcurve = PolynomialCurve::<Point3>(surface_vcoef);
    let surface = PolynomialSurface::by_tensor(ucurve, vcurve);

    let pcurve0 = ParameterCurve::new(curve, surface);
    let pcurve1 = pcurve0.surface().composite(pcurve0.curve());

    prop_assert_derivative_agreement(pcurve0.derivative_n(n, t), pcurve1.derivative_n(n, t))?;

    let ders0 = (0..=n)
        .map(|i| pcurve0.derivative_n(i, t))
        .collect::<Vec<_>>();

    let ders1 = pcurve0.derivatives(n, t);

    prop_assert_eq!(ders0.len(), ders1.len());

    ders0
        .into_iter()
        .zip(&*ders1)
        .try_for_each(|(left, right)| prop_assert_derivative_agreement(left, *right))?;
    Ok(())
}

proptest! {
    #[test]
    fn pcurve_derivation(
        curve_coef in prop::array::uniform6(prop::array::uniform2(-1f64..=1.0)),
        curve_degree in 1usize..6,
        surface_ucoef in prop::array::uniform6(prop::array::uniform3(-1f64..=1.0)),
        surface_vcoef in prop::array::uniform6(prop::array::uniform3(-1f64..=1.0)),
        surface_udegree in 1usize..6,
        surface_vdegree in 1usize..6,
        n in 0usize..=5,
        t in -1f64..=1.0,
    ) {
        exec_pcurve_derivation(
            &curve_coef[..=curve_degree],
            &surface_ucoef[..=surface_udegree],
            &surface_vcoef[..=surface_vdegree],
            n,
            t,
        )?;
    }
}

#[test]
fn high_order_pcurve_derivative_regression_from_ci() {
    // Shrunk from the CI proptest failure on 2026-11-18.
    let curve = [
        [-0.017_093_239_736_548_603, 0.066_779_720_106_981_45],
        [0.0, 0.210_174_389_360_958_97],
        [-0.942_841_532_951_020_7, 0.838_479_701_226_748_2],
        [0.0, -0.955_215_861_886_636_2],
        [0.825_498_531_480_303_9, -0.883_044_902_082_598_1],
        [-0.529_182_018_438_617_5, 0.412_348_547_486_633_47],
    ];
    let u = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -0.543_720_736_235_002],
    ];
    let v = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -0.942_648_615_904_186_1],
    ];
    exec_pcurve_derivation(&curve, &u, &v, 5, 0.969_411_026_851_627_3)
        .expect("composite and explicit parameter-curve derivatives must agree");
}
