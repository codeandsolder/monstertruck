#![expect(
    clippy::many_single_char_names,
    reason = "numerical geometry algorithms use compact paired parameter notation"
)]
#![expect(
    clippy::similar_names,
    reason = "numerical geometry algorithms use compact paired parameter notation"
)]

//! Scalar-generic algorithm helpers.
//!
//! These mirror the `f64`-bound algorithms in [`crate::algo`] but work with
//! any scalar satisfying `cgmath::BaseFloat`.

use num_traits::Float;

#[inline]
pub(super) fn scalar_from_usize<S: Float>(mut value: usize) -> S {
    let mut result = S::zero();
    let mut bit = S::one();
    while value != 0 {
        if value & 1 == 1 {
            result = result + bit;
        }
        value >>= 1;
        if value != 0 {
            bit = bit + bit;
        }
    }
    result
}

/// Curve algorithms.
pub mod curve;
/// Surface algorithms.
pub mod surface;
