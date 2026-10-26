use monstertruck_geometry::prelude::*;
use std::num::NonZeroUsize;

/// Profile shape for fillet operations.
#[derive(Debug, Clone, Default)]
pub enum FilletProfile {
    /// Circular arc cross-section (traditional fillet).
    #[default]
    Round,
    /// Flat ruled surface (chamfer/bevel) using the ordinary symmetric radius.
    Chamfer,
    /// Flat ruled surface with independent setbacks from the two carrier surfaces.
    ///
    /// `first_setback_mm` applies to the first face passed to the single-edge
    /// fillet operation and `second_setback_mm` to the second. This variant is
    /// intentionally single-edge only until wire ordering has equally explicit
    /// per-carrier semantics.
    AsymmetricChamfer {
        /// Setback measured along the first carrier surface normal construction.
        first_setback_mm: f64,
        /// Setback measured along the second carrier surface normal construction.
        second_setback_mm: f64,
    },
    /// V-shaped ridge: two straight segments meeting at the transit point.
    Ridge,
    /// User-provided 2D profile curve. Domain `[0, 1]`, maps `(0, 0)` -> `contact0`,
    /// `(1, 0)` -> `contact1`, y-axis = displacement toward transit.
    Custom(Box<BsplineCurve<Point2>>),
}

impl FilletProfile {
    /// Whether this profile consumes [`FilletOptions::radius`].
    ///
    /// Asymmetric chamfers are fully defined by their two carrier setbacks; an
    /// unrelated radius specification must not affect validation or geometry.
    pub(crate) const fn uses_radius(&self) -> bool {
        !matches!(self, Self::AsymmetricChamfer { .. })
    }
}

#[cfg(test)]
mod profile_contract_tests {
    use super::*;

    #[test]
    fn asymmetric_chamfer_does_not_consume_radius() {
        assert!(
            !FilletProfile::AsymmetricChamfer {
                first_setback_mm: 0.2,
                second_setback_mm: 0.3,
            }
            .uses_radius()
        );
        assert!(FilletProfile::Round.uses_radius());
        assert!(FilletProfile::Chamfer.uses_radius());
        assert!(FilletProfile::Ridge.uses_radius());
    }
}

/// Radius specification for fillet operations.
pub enum RadiusSpec {
    /// Constant radius along the entire edge/wire.
    Constant(f64),
    /// Variable radius as a function of normalized parameter `t` in `[0, 1]`.
    ///
    /// Supported for single-edge and wire fillets.
    /// For closed wires, endpoint continuity is required:
    /// `f(0.0)` must be near `f(1.0)`.
    Variable(Box<dyn Fn(f64) -> f64>),
    /// Per-edge radius. Length must match the edge count passed to [`fillet_edges`](super::fillet_edges).
    PerEdge(Vec<f64>),
}

impl std::fmt::Debug for RadiusSpec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Constant(r) => f.debug_tuple("Constant").field(r).finish(),
            Self::Variable(_) => f.debug_tuple("Variable").field(&"<fn>").finish(),
            Self::PerEdge(v) => f.debug_tuple("PerEdge").field(v).finish(),
        }
    }
}

/// Options for fillet operations.
#[derive(Debug)]
pub struct FilletOptions {
    /// Radius specification.
    pub radius: RadiusSpec,
    /// Number of divisions for the rolling ball algorithm. Default: 5.
    pub divisions: NonZeroUsize,
    /// Profile shape. Default: [`FilletProfile::Round`].
    pub profile: FilletProfile,
}

impl Default for FilletOptions {
    fn default() -> Self {
        Self {
            radius: RadiusSpec::Constant(0.1),
            divisions: NonZeroUsize::new(5).unwrap_or(NonZeroUsize::MIN),
            profile: FilletProfile::default(),
        }
    }
}

impl FilletOptions {
    /// Creates options with a constant radius.
    pub fn constant(radius: f64) -> Self {
        Self {
            radius: RadiusSpec::Constant(radius),
            ..Default::default()
        }
    }

    /// Creates options with a variable radius function.
    pub fn variable(radius: impl Fn(f64) -> f64 + 'static) -> Self {
        Self {
            radius: RadiusSpec::Variable(Box::new(radius)),
            ..Default::default()
        }
    }

    /// Sets the fillet radius specification.
    pub fn with_radius(mut self, radius: RadiusSpec) -> Self {
        self.radius = radius;
        self
    }

    /// Sets the division count used by the fillet algorithm.
    pub fn with_division(mut self, division: NonZeroUsize) -> Self {
        self.divisions = division;
        self
    }

    /// Sets the fillet profile.
    pub fn with_profile(mut self, profile: FilletProfile) -> Self {
        self.profile = profile;
        self
    }
}
