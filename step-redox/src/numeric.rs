use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use num_traits::ToPrimitive;

/// Internal numerical consistency bound, measured in millimetres.
///
/// This is for evaluating geometry the kernel has *already* constructed and
/// for proving exact geometric identities. It is not a supplier-model
/// feature-healing threshold.
pub const NUMERICAL_GEOMETRY_TOL_MM: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.numerical_geometry_mm;

/// Default absolute deviation budget for recognizing an *intended* common
/// feature in imported supplier CAD, measured in millimetres (10 micrometres).
///
/// The budget is not a permission to connect arbitrarily close unrelated
/// bodies, collapse genuine small features, or accept a non-manifold result.
/// A healing operation must still build a coherent common curve/surface and
/// validate its output topology. Keep numerical consistency tests separate.
pub const DEFAULT_FEATURE_HEAL_TOL_MM: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.feature_heal_mm;

/// Can two measurements be treated as the same *intended* feature?
///
/// Always compare independent source measurements against this threshold;
/// do not use it in place of exact topology identity or kernel validation.
#[must_use]
pub fn near_imported_feature_mm(residual_mm: f64) -> bool {
    residual_mm.is_finite() && (0.0..=DEFAULT_FEATURE_HEAL_TOL_MM).contains(&residual_mm)
}

const MAX_EXACT_F64_INTEGER: u64 = 1_u64 << 53;
const TWO_POW_32: f64 = 4_294_967_296.0;

/// Convert an `i64` to `f64` only when the integer is exactly representable.
///
/// STEP permits integer-valued numeric parameters wider than the exact integer
/// range of IEEE-754 binary64. Geometry recovery must not silently round them.
pub fn exact_i64_to_f64(value: i64) -> Option<f64> {
    let magnitude = value.unsigned_abs();
    if magnitude > MAX_EXACT_F64_INTEGER {
        return None;
    }

    let high = u32::try_from(magnitude >> 32).ok()?;
    let low = u32::try_from(magnitude & u64::from(u32::MAX)).ok()?;
    let magnitude_f64 = f64::from(high).mul_add(TWO_POW_32, f64::from(low));
    Some(if value.is_negative() {
        -magnitude_f64
    } else {
        magnitude_f64
    })
}

/// Round a finite `f64` to an `i64`, rejecting values outside the integer range.
///
/// This is for geometric quantization where rounding is part of the algorithm.
/// It avoids Rust's saturating float-to-integer `as` semantics hiding overflow.
pub fn rounded_f64_to_i64(value: f64) -> Option<i64> {
    value
        .is_finite()
        .then(|| value.round())
        .and_then(|value| value.to_i64())
}

/// Floor a finite `f64` to an `i64`, rejecting values outside the integer range.
///
/// This is for lattice-coordinate decomposition where flooring is part of the
/// algorithm and out-of-range coordinates must fail closed.
pub fn floored_f64_to_i64(value: f64) -> Option<i64> {
    value
        .is_finite()
        .then(|| value.floor())
        .and_then(|value| value.to_i64())
}

/// Convert a `usize` to `f64` only when the integer is exactly representable.
pub fn exact_usize_to_f64(value: usize) -> Option<f64> {
    let value = u64::try_from(value).ok()?;
    if value > MAX_EXACT_F64_INTEGER {
        return None;
    }

    let high = u32::try_from(value >> 32).ok()?;
    let low = u32::try_from(value & u64::from(u32::MAX)).ok()?;
    Some(f64::from(high).mul_add(TWO_POW_32, f64::from(low)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn imported_feature_tolerance_accepts_micron_scale_fits_without_weakening_numerics() {
        assert!(near_imported_feature_mm(0.006_701_7));
        assert!(near_imported_feature_mm(0.01));
        assert!(!near_imported_feature_mm(0.010_001));
        assert!(!near_imported_feature_mm(0.05));
        assert!(!near_imported_feature_mm(-1.0));
        assert!(!near_imported_feature_mm(f64::NAN));
        assert!(!near_imported_feature_mm(f64::INFINITY));
    }

    #[test]
    fn i64_conversion_preserves_exact_binary64_integer_range() {
        assert_eq!(exact_i64_to_f64(0), Some(0.0));
        assert_eq!(exact_i64_to_f64(-42), Some(-42.0));
        assert_eq!(
            exact_i64_to_f64(9_007_199_254_740_992),
            Some(9_007_199_254_740_992.0)
        );
        assert_eq!(
            exact_i64_to_f64(-9_007_199_254_740_992),
            Some(-9_007_199_254_740_992.0)
        );
        assert_eq!(exact_i64_to_f64(9_007_199_254_740_993), None);
        assert_eq!(exact_i64_to_f64(-9_007_199_254_740_993), None);
        assert_eq!(exact_i64_to_f64(i64::MIN), None);
        assert_eq!(exact_i64_to_f64(i64::MAX), None);
    }

    #[test]
    fn float_quantization_rejects_nonfinite_and_out_of_range_values() {
        assert_eq!(rounded_f64_to_i64(2.6), Some(3));
        assert_eq!(rounded_f64_to_i64(-2.6), Some(-3));
        assert_eq!(floored_f64_to_i64(2.9), Some(2));
        assert_eq!(floored_f64_to_i64(-2.1), Some(-3));
        assert_eq!(rounded_f64_to_i64(f64::NAN), None);
        assert_eq!(floored_f64_to_i64(f64::INFINITY), None);
        assert_eq!(rounded_f64_to_i64(1.0e30), None);
        assert_eq!(floored_f64_to_i64(-1.0e30), None);
    }

    #[test]
    fn usize_conversion_preserves_exact_binary64_integer_range() {
        assert_eq!(exact_usize_to_f64(0), Some(0.0));
        assert_eq!(exact_usize_to_f64(42), Some(42.0));
        if let (Ok(last_exact), Ok(first_inexact)) = (
            usize::try_from(9_007_199_254_740_992_u64),
            usize::try_from(9_007_199_254_740_993_u64),
        ) {
            assert_eq!(
                exact_usize_to_f64(last_exact),
                Some(9_007_199_254_740_992.0)
            );
            assert_eq!(exact_usize_to_f64(first_inexact), None);
        }
    }
}
