//! User-facing geometry accuracy contracts.
//!
//! This is deliberately separate from [`crate::tolerance`]. Numerical kernel
//! tolerances keep floating-point algorithms robust; a geometry policy says
//! whether an operation is allowed to change represented geometry at all.

use crate::tolerance::KERNEL_EPSILON;
use std::{error::Error, fmt};

/// Accuracy contract for geometry-changing operations.
///
/// Monstertruck's kernel coordinates do not currently carry physical units, so
/// [`GeometryPolicy::Standard`] is expressed in the caller's current model
/// units. A physical preset such as "0.01 mm" therefore belongs at a unit-aware
/// adapter boundary, not in this unitless core type.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeometryPolicy {
    /// Admit only operations that establish geometric equivalence.
    ///
    /// Numerical thresholds may still be used to evaluate or prove geometry;
    /// they are not an approximation budget. An algorithm that only has an
    /// approximate implementation must refuse this policy.
    Equivalent,
    /// Admit bounded geometric changes under a global tolerance in current model units.
    Standard {
        /// Global geometric tolerance budget in current model units.
        tolerance: f64,
    },
}

impl GeometryPolicy {
    /// Constructs the Standard modeling-tolerance policy in current model units.
    #[must_use]
    pub const fn standard(tolerance: f64) -> Self { Self::Standard { tolerance } }

    /// Returns the modeling tolerance, or [`None`] for [`Self::Equivalent`].
    #[must_use]
    pub const fn tolerance(self) -> Option<f64> {
        match self {
            Self::Equivalent => None,
            Self::Standard { tolerance } => Some(tolerance),
        }
    }

    /// Validates the policy against the current kernel's numerical resolution.
    ///
    /// # Errors
    ///
    /// [`GeometryPolicyError::InvalidStandardTolerance`] is returned when a
    /// Standard tolerance is non-finite or smaller than the kernel can
    /// robustly resolve.
    pub fn validate(self) -> Result<Self, GeometryPolicyError> {
        match self {
            Self::Equivalent => Ok(self),
            Self::Standard { tolerance }
                if tolerance.is_finite() && tolerance >= KERNEL_EPSILON =>
            {
                Ok(self)
            }
            Self::Standard { tolerance } => {
                Err(GeometryPolicyError::InvalidStandardTolerance { tolerance })
            }
        }
    }
}

/// Invalid user-facing geometry policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GeometryPolicyError {
    /// Standard's tolerance was non-finite or below kernel numerical resolution.
    InvalidStandardTolerance {
        /// Rejected modeling tolerance in current model units.
        tolerance: f64,
    },
}

impl fmt::Display for GeometryPolicyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStandardTolerance { tolerance } => write!(
                formatter,
                "Standard tolerance must be finite and at least KERNEL_EPSILON ({KERNEL_EPSILON}); got {tolerance}",
            ),
        }
    }
}

impl Error for GeometryPolicyError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equivalent_has_no_error_budget() {
        assert_eq!(GeometryPolicy::Equivalent.tolerance(), None);
        assert_eq!(
            GeometryPolicy::Equivalent.validate(),
            Ok(GeometryPolicy::Equivalent)
        );
    }

    #[test]
    fn standard_uses_model_space_budget() {
        let policy = GeometryPolicy::standard(0.01);
        assert_eq!(policy.tolerance(), Some(0.01));
        assert_eq!(policy.validate(), Ok(policy));
    }

    #[test]
    fn standard_rejects_unrepresentable_budgets() {
        for tolerance in [f64::NAN, f64::INFINITY, -1.0, 0.0, KERNEL_EPSILON * 0.5] {
            assert!(matches!(
                GeometryPolicy::standard(tolerance).validate(),
                Err(GeometryPolicyError::InvalidStandardTolerance { .. })
            ));
        }
    }
}
