use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

// Each term describes three numerator knot indices and two denominator
// intervals: [i, j, k, dl0, dl1, dr0, dr1]. Both middle spans share the
// evaluator, but use different choices of affine numerator factors.
const MIDDLE_SPAN_TERMS: [[[usize; 7]; 3]; 2] = [
    [
        [0, 0, 2, 3, 0, 2, 0],
        [0, 1, 3, 3, 0, 3, 1],
        [1, 1, 4, 4, 1, 3, 1],
    ],
    [
        [0, 3, 3, 3, 0, 3, 1],
        [1, 4, 3, 4, 1, 3, 1],
        [2, 4, 4, 4, 1, 4, 2],
    ],
];

/// Cubic five-knot basis and its first two analytic derivatives.
///
/// The derivative order is a compile-time constant so the common knot-span
/// and numerator descriptors are shared without a runtime derivative switch.
/// For repeated knots, terms with a zero-width denominator have zero Cox–de
/// Boor weight and must be omitted instead of producing 0/0.
#[inline]
fn cubic_basis<const ORDER: usize>(u: f64, a: &[f64]) -> f64 {
    debug_assert!(ORDER <= 2);
    if u < a[0] || a[4] <= u {
        return 0.0;
    }

    if u < a[1] || a[3] <= u {
        let leading = u < a[1];
        let d = if leading { u - a[0] } else { a[4] - u };
        let denominator = if leading {
            (a[3] - a[0]) * (a[2] - a[0]) * (a[1] - a[0])
        } else {
            (a[4] - a[1]) * (a[4] - a[2]) * (a[4] - a[3])
        };
        let numerator = match ORDER {
            0 => d * d * d,
            1 => (if leading { 3.0 } else { -3.0 }) * d * d,
            _ => 6.0 * d,
        };
        return numerator / denominator;
    }

    let front = u < a[2];
    let (terms, width, sign) = if front {
        (&MIDDLE_SPAN_TERMS[0], a[2] - a[1], 1.0)
    } else {
        (&MIDDLE_SPAN_TERMS[1], a[3] - a[2], -1.0)
    };
    let sum = terms
        .iter()
        .fold(0.0, |sum, &[i, j, k, dl0, dl1, dr0, dr1]| {
            let denom = (a[dl0] - a[dl1]) * (a[dr0] - a[dr1]);
            if denom == 0.0 {
                // Knot multiplicity collapses the term's Cox–de Boor interval.
                return sum;
            }

            let f0 = u - a[i];
            let f1 = sign * (u - a[j]);
            let f2 = a[k] - u;
            let numerator = match ORDER {
                0 => f0 * f1 * f2,
                1 => f0.mul_add(-f1, (sign * f0).mul_add(f2, f1 * f2)),
                _ => 2.0 * sign.mul_add(f2 - f0, -f1),
            };
            sum + numerator / denom
        });
    sum / width
}

/// Evaluates the cubic B-spline basis for a local, nondecreasing five-knot vector.
///
/// The support is half-open: `[a[0], a[4])`.
#[inline]
pub(super) fn basis_function(u: f64, a: &[f64]) -> f64 { cubic_basis::<0>(u, a) }

/// Evaluates the given derivative order of the cubic B-spline basis.
#[inline]
pub(super) fn basis_function_der(u: f64, a: &[f64], der_order: usize) -> f64 {
    match der_order {
        0 => cubic_basis::<0>(u, a),
        1 => cubic_basis::<1>(u, a),
        2 => cubic_basis::<2>(u, a),
        _ => {
            // Preserve the existing finite-difference fallback above order 2.
            let h = DEFAULT_NUMERICAL_TOLERANCES.finite_difference_step;
            (basis_function_der(u + h, a, der_order - 1)
                - basis_function_der(u - h, a, der_order - 1))
                / (2.0 * h)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{basis_function, basis_function_der};

    // Independent Cox–de Boor recurrence for the mathematical cubic B-spline
    // basis and derivatives. Zero-width intervals contribute zero.
    fn reference(knots: &[f64; 5], u: f64, degree: usize, index: usize, order: usize) -> f64 {
        if degree == 0 {
            return f64::from(order == 0 && knots[index] <= u && u < knots[index + 1]);
        }
        let left = knots[index + degree] - knots[index];
        let right = knots[index + degree + 1] - knots[index + 1];
        let degree_float = [0.0, 1.0, 2.0, 3.0][degree];
        let left_term = if left > 0.0 {
            if order == 0 {
                (u - knots[index]) / left * reference(knots, u, degree - 1, index, 0)
            } else {
                degree_float / left * reference(knots, u, degree - 1, index, order - 1)
            }
        } else {
            0.0
        };
        let right_term = if right > 0.0 {
            if order == 0 {
                (knots[index + degree + 1] - u) / right
                    * reference(knots, u, degree - 1, index + 1, 0)
            } else {
                -degree_float / right * reference(knots, u, degree - 1, index + 1, order - 1)
            }
        } else {
            0.0
        };
        left_term + right_term
    }

    #[test]
    fn cubic_basis_and_derivatives_agree_with_cox_de_boor() {
        let knot_vectors = [
            [0.0, 1.0, 2.0, 3.0, 4.0],
            [0.0, 0.02, 0.3, 0.35, 1.0],
            [0.0, 0.0, 0.2, 0.7, 1.0],
            [0.0, 0.0, 0.0, 0.5, 1.0],
            [0.0, 0.0, 0.5, 0.5, 1.0],
            [0.0, 0.2, 0.2, 0.2, 1.0],
            [0.0, 0.4, 0.4, 1.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0, 1.0, 1.0],
        ];
        for knots in &knot_vectors {
            for interval in knots.windows(2).filter(|pair| pair[0] < pair[1]) {
                for fraction in [0.17f64, 0.37, 0.63, 0.87] {
                    let u = fraction.mul_add(interval[1] - interval[0], interval[0]);
                    for order in 0..=2 {
                        let expected = reference(knots, u, 3, 0, order);
                        let actual = basis_function_der(u, knots, order);
                        let tol = 1.0e-9 * expected.abs().max(1.0);
                        assert!(
                            actual.is_finite() && (actual - expected).abs() <= tol,
                            "{knots:?}: u={u}, order={order}, expected={expected}, got={actual}"
                        );
                    }
                }
            }
            for u in [knots[0] - 1.0, knots[4], knots[4] + 1.0] {
                assert_eq!(basis_function(u, knots), 0.0);
            }
        }
    }
}
