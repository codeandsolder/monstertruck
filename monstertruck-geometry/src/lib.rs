//! Geometric primitives: knot vectors, B-splines, NURBS, and T-splines.
//!
//! # Examples
//!
//! ```
//! use monstertruck_geometry::prelude::*;
//!
//! // A quadratic Bezier curve from (0,0,0) through (1,1,0) to (2,0,0)
//! let knot_vec = KnotVector::bezier_knot(2);
//! let ctrl_pts = vec![
//!     Point3::new(0.0, 0.0, 0.0),
//!     Point3::new(1.0, 1.0, 0.0),
//!     Point3::new(2.0, 0.0, 0.0),
//! ];
//! let curve = BsplineCurve::new_unchecked(knot_vec, ctrl_pts);
//!
//! let mid = curve.evaluate(0.5);    // Point3 at parameter t=0.5
//! let tan = curve.derivative(0.5);  // tangent vector
//! let (t0, t1) = curve.range_tuple(); // (0.0, 1.0)
//! ```

#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

use monstertruck_core::bounding_box::Bounded;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, ops::Bound};

const INCLUDE_CURVE_TRIALS: usize = 100;
const PRESEARCH_DIVISION: usize = 50;

/// Converts a discrete index/count to the geometry scalar domain without scattering lossy casts.
///
/// Like Rust's integer-to-float conversion, values above f64's exact-integer range are rounded.
/// Keeping the conversion in one helper makes that approximation boundary explicit and auditable.
pub(crate) fn usize_to_f64(value: usize) -> f64 {
    // `usize` always has a finite `f64` representation; values above 2^53 are rounded
    // exactly as by Rust's primitive integer-to-float conversion.
    value.to_f64().unwrap_or(f64::INFINITY)
}

/// Number of uniform segments needed so every segment is strictly shorter than `max_step`.
///
/// Returns `None` for non-finite/negative spans, non-positive steps, or a count that cannot be
/// represented by `usize`. Keeping the float-to-integer conversion here makes every closed-form
/// parameter divider use the same checked boundary.
fn strict_segment_count(span: f64, max_step: f64) -> Option<usize> {
    if !span.is_finite() || span < 0.0 || !max_step.is_finite() || max_step <= 0.0 {
        return None;
    }
    (span / max_step).floor().to_usize()?.checked_add(1)
}

/// Converts a non-negative continuous segment estimate to an integer count.
///
/// A zero estimate still means one segment: callers interpolate using `i / segments`, so allowing
/// zero here would turn a zero-width parameter range into `0 / 0` and manufacture NaNs.
fn ceil_segment_count(estimate: f64) -> Option<usize> {
    if !estimate.is_finite() || estimate < 0.0 {
        return None;
    }
    Some(estimate.ceil().to_usize()?.max(1))
}

fn pow_usize(mut base: f64, mut exponent: usize) -> f64 {
    let mut result = 1.0;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result *= base;
        }
        exponent >>= 1;
        if exponent != 0 {
            base *= base;
        }
    }
    result
}

/// Closed-form surface divisions share the adaptive divider's measured work ceiling.
fn closed_form_surface_segment_counts(
    u_span: f64,
    v_span: f64,
    max_step: f64,
) -> Option<(usize, usize)> {
    let u_segments = strict_segment_count(u_span, max_step)?;
    let v_segments = strict_segment_count(v_span, max_step)?;
    let cells = u_segments.checked_mul(v_segments)?;
    let max_cells =
        usize::try_from(monstertruck_traits::algo::surface::MAX_PARAMETER_DIVISION_CELLS).ok()?;
    (cells <= max_cells).then_some((u_segments, v_segments))
}

/// re-export `monstertruck_core`
pub mod base {
    pub use monstertruck_core::{
        assert_near, assert_near2, bounding_box::BoundingBox, cgmath64::*, hash, hash::HashGen,
        prop_assert_near, prop_assert_near2, tolerance::*,
    };
    pub use monstertruck_traits::*;
}
/// NURBS and B-spline curves, surfaces, and knot vectors.
pub mod nurbs;

/// Error types for geometry operations.
pub mod errors;

/// Concrete geometric primitives: [`Plane`](crate::specifieds::Plane), [`Sphere`](crate::specifieds::Sphere), [`Line`](crate::specifieds::Line), etc.
pub mod specifieds;

/// Composite geometry: revolved curves, intersection curves, processor wrappers.
pub mod decorators;

/// T-Spline and T-NURCC surface types.
pub mod t_spline;

mod analytic_surface;
/// [`DeterministicContentHash`](monstertruck_core::DeterministicContentHash) impls for geometry types.
mod content_hash_impls;

/// Trait for extracting an exact polynomial B-spline surface representation.
mod bspline_conversion;
mod parameter_boundary;

/// re-export all modules.
pub mod prelude {
    use crate::{
        analytic_surface, base, bspline_conversion, decorators, errors, nurbs, parameter_boundary,
        specifieds, t_spline,
    };
    pub use analytic_surface::{
        AnalyticSurfaceKind, HomogeneousExtrusionSurface, SphericalRevolutionSurface,
        SurfaceParameterAxis, TryIntoAnalyticSurfaceKind,
    };
    pub use base::*;
    pub use bspline_conversion::{
        HomogeneousSurfaceConversion, SupportsExactPatchDomains, SurfaceFrameAxes,
        SurfaceParameterRectangle, TryIntoBsplineSurface, TryIntoHomogeneousBsplineCurve,
        TryIntoHomogeneousBsplineSurface,
    };
    pub use decorators::*;
    pub use errors::*;
    pub use nurbs::*;
    pub use parameter_boundary::BoundaryCurve2D;
    pub use specifieds::*;
    pub use t_spline::*;
}

#[cfg(test)]
mod subdivision_count_tests {
    use super::*;

    #[test]
    fn zero_estimate_still_has_one_segment() {
        assert_eq!(ceil_segment_count(0.0), Some(1));
    }

    #[test]
    fn invalid_or_unrepresentable_segment_counts_are_refused() {
        assert_eq!(strict_segment_count(f64::INFINITY, 1.0), None);
        assert_eq!(strict_segment_count(1.0, 0.0), None);
        assert_eq!(ceil_segment_count(f64::INFINITY), None);
        assert_eq!(ceil_segment_count(-1.0), None);
    }

    #[test]
    fn usize_power_handles_orders_beyond_i32() {
        let huge_even = i32::MAX as usize + 1;
        let huge_odd = huge_even + 1;
        assert_eq!(pow_usize(-1.0, huge_even), 1.0);
        assert_eq!(pow_usize(-1.0, huge_odd), -1.0);
    }
}
