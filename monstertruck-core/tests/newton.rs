use monstertruck_core::{math::*, newton::*, prop_assert_near, tolerance::*};
use proptest::prelude::*;
use std::f64::consts::PI;

fn cubic(a: [f64; 4], root: f64) -> impl Fn(f64) -> CalcOutput<f64, f64> {
    let poly =
        move |x: f64| (a[3] * x * x).mul_add(x, (a[2] * x).mul_add(x, f64::mul_add(a[1], x, a[0])));
    move |x| CalcOutput {
        value: poly(x) - poly(root),
        derivation: (3.0 * a[3] * x).mul_add(x, f64::mul_add(2.0 * a[2], x, a[1])),
    }
}

proptest! {
    #[test]
    fn test_newton1(
        a in prop::array::uniform4(-10.0f64..=10.0f64),
        x0 in -10.0f64..=10.0f64,
        delta in -0.5f64..=0.5f64,
    ) {
        let function = cubic(a, x0);
        // Unconstrained cubic Newton iterates need not converge, even
        // when a root exists and no iterate has a singular derivative.
        // Success must nevertheless be an actual root.
        if let Ok(res) = solve(&function, x0 + delta, 100, KERNEL_EPSILON) {
            prop_assert_near!(function(res).value, 0.0);
        }
    }

    #[test]
    fn newton_converges_in_certified_local_basin(
        root in -10.0f64..=10.0f64,
        quadratic in -0.25f64..=0.25f64,
        cubic in -0.25f64..=0.25f64,
        delta in -0.5f64..=0.5f64,
    ) {
        // With e=x-root and |e|<=1/2, f=e+b*e²+c*e³ has
        // |f'|>=9/16 and the exact next error
        // |e'|=|e²(b+2ce)/f'|<= (8/9)e² <= (4/9)|e|.
        // This supplies a genuine local convergence guarantee, not a
        // claim of global convergence for arbitrary cubic coefficients.
        let function = |x: f64| {
            let e = x - root;
            CalcOutput {
                value: e * e.mul_add(e.mul_add(cubic, quadratic), 1.0),
                derivation: e.mul_add(e.mul_add(3.0 * cubic, 2.0 * quadratic), 1.0),
            }
        };
        let result = solve(&function, root + delta, 100, KERNEL_EPSILON)
            .map_err(|log| TestCaseError::fail(log.to_string()))?;
        prop_assert_near!(result, root);
        prop_assert_near!(function(result).value, 0.0);
    }

    #[test]
    fn test_newton2(
        n in prop::array::uniform2(-10.0f64..=10.0f64),
        delta in prop::array::uniform2(-0.5f64..=0.5f64),
    ) {
        let n = Vector2::from(n);
        if n.so_small() {
            return Ok(());
        }
        let n = n.normalize();
        let function = |vec: Vector2| CalcOutput {
            value: Vector2::new(vec.norm_squared() - 1.0, vec.dot(&n)),
            derivation: Matrix2::new(2.0 * vec.x, 2.0 * vec.y, n.x, n.y),
        };
        let hint = Rotation2::new(-(PI / 2.0)).into_inner() * n + Vector2::from(delta);
        match solve(function, hint, 10, KERNEL_EPSILON) {
            Ok(res) => prop_assert_near!(function(res).value, Vector2::zero()),
            Err(log) => prop_assert!(log.degenerate(), "{log}"),
        }
    }
}

#[test]
fn cubic_newton_may_exhaust_trials_without_singular_jacobian() {
    let function = cubic(
        [
            0.0,
            -2.408_548_620_740_545_3,
            -8.443_626_344_016_24,
            -7.455_990_877_269_495,
        ],
        0.0,
    );
    let failure = solve(&function, -0.435_058_868_007_126_73, 100, KERNEL_EPSILON)
        .expect_err("this local Newton trajectory does not converge");
    assert!(!failure.degenerate());
    let result = solve(&function, 0.01, 100, KERNEL_EPSILON).expect("seed in the root's basin");
    assert!(function(result).value.so_small());
}
