//! Small dense numerical operations for CAD parameter and transform analysis.
use nalgebra::{Matrix3, Matrix4, linalg::Schur};

/// Factor a finite full-rank 3D homogeneous transform into orthogonal and
/// upper-triangular factors with nonnegative diagonal scales.
#[must_use]
pub fn affine_qr3(matrix: Matrix3<f64>) -> Option<(Matrix3<f64>, Matrix3<f64>)> {
    if !matrix.iter().all(|v| v.is_finite()) {
        return None;
    }
    let (q, r) = matrix.qr().unpack();
    (r.diagonal().iter().all(|v| v.is_finite() && *v > 0.0) && q.iter().all(|v| v.is_finite()))
        .then_some((q, r))
}

/// Factor a finite full-rank 4D homogeneous transform into orthogonal and
/// upper-triangular factors with nonnegative diagonal scales.
#[must_use]
pub fn affine_qr4(matrix: Matrix4<f64>) -> Option<(Matrix4<f64>, Matrix4<f64>)> {
    if !matrix.iter().all(|v| v.is_finite()) {
        return None;
    }
    let (q, r) = matrix.qr().unpack();
    (r.diagonal().iter().all(|v| v.is_finite() && *v > 0.0) && q.iter().all(|v| v.is_finite()))
        .then_some((q, r))
}

/// Complex roots of `x³ + a x² + b x + c`, obtained from a
/// scaled companion matrix. Callers select real roots according to their
/// numerical tolerance.
#[must_use]
pub fn monic_cubic_roots([a, b, c]: [f64; 3]) -> Option<[nalgebra::Complex<f64>; 3]> {
    let scale = a.abs().max(b.abs().sqrt()).max(c.abs().cbrt());
    if ![a, b, c, scale].iter().all(|v| v.is_finite()) {
        return None;
    }
    if scale == 0.0 {
        return Some([nalgebra::Complex::new(0.0, 0.0); 3]);
    }
    let matrix = Matrix3::new(
        0.0,
        0.0,
        -(c / scale / scale / scale),
        1.0,
        0.0,
        -(b / scale / scale),
        0.0,
        1.0,
        -a / scale,
    );
    let roots = Schur::try_new(matrix, f64::EPSILON, 256)?.complex_eigenvalues();
    Some(roots.map(|v| v * scale).into())
}

/// Roots of `x⁴ + a x³ + b x² + c x + d`, obtained from a scaled companion matrix.
#[must_use]
pub fn monic_quartic_roots([a, b, c, d]: [f64; 4]) -> Option<[nalgebra::Complex<f64>; 4]> {
    let scale = a
        .abs()
        .max(b.abs().sqrt())
        .max(c.abs().cbrt())
        .max(d.abs().sqrt().sqrt());
    if ![a, b, c, d, scale].iter().all(|v| v.is_finite()) {
        return None;
    }
    if scale == 0.0 {
        return Some([nalgebra::Complex::new(0.0, 0.0); 4]);
    }
    let matrix = Matrix4::new(
        0.0,
        0.0,
        0.0,
        -(d / scale / scale / scale / scale),
        1.0,
        0.0,
        0.0,
        -(c / scale / scale / scale),
        0.0,
        1.0,
        0.0,
        -(b / scale / scale),
        0.0,
        0.0,
        1.0,
        -a / scale,
    );
    let roots = Schur::try_new(matrix, f64::EPSILON, 512)?.complex_eigenvalues();
    Some(roots.map(|v| v * scale).into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qr_preserves_affine_factors_and_reflection() {
        let transform = Matrix4::new(
            -2.0, 1.0, 0.0, 3.0, 0.0, 3.0, 2.0, 4.0, 0.0, 0.0, 4.0, 5.0, 0.0, 0.0, 0.0, 1.0,
        );
        let Some((q, r)) = affine_qr4(transform) else {
            panic!("full rank affine transform")
        };
        assert!((q.transpose() * q - Matrix4::identity()).norm() < 1e-12);
        assert!((q * r - transform).norm() < 1e-12);
        assert!(q.determinant() < 0.0);
        assert!((r.diagonal() - nalgebra::Vector4::new(2.0, 3.0, 4.0, 1.0)).norm() < 1e-12);
    }

    #[test]
    fn companion_roots_scale_across_magnitudes() {
        for scale in [1e-9_f64, 1.0, 1e9] {
            let Some(mut roots) =
                monic_cubic_roots([-6.0 * scale, 11.0 * scale.powi(2), -6.0 * scale.powi(3)])
            else {
                panic!("three distinct real roots")
            };
            roots.sort_by(|a, b| a.re.total_cmp(&b.re));
            for (root, expected) in roots.into_iter().zip([1.0, 2.0, 3.0]) {
                assert!((root.re / scale - expected).abs() < 1e-12);
                assert!(root.im.abs() / scale < 1e-12);
            }
            let Some(mut roots) =
                monic_quartic_roots([0.0, -5.0 * scale.powi(2), 0.0, 4.0 * scale.powi(4)])
            else {
                panic!("four distinct real roots")
            };
            roots.sort_by(|a, b| a.re.total_cmp(&b.re));
            for (root, expected) in roots.into_iter().zip([-2.0, -1.0, 1.0, 2.0]) {
                assert!((root.re / scale - expected).abs() < 1e-12);
                assert!(root.im.abs() / scale < 1e-12);
            }
        }
    }

    #[test]
    fn zero_polynomials_have_zero_roots() {
        assert_eq!(
            monic_cubic_roots([0.0; 3]),
            Some([nalgebra::Complex::new(0.0, 0.0); 3])
        );
        assert_eq!(
            monic_quartic_roots([0.0; 4]),
            Some([nalgebra::Complex::new(0.0, 0.0); 4])
        );
    }
}
