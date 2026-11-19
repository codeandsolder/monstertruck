use monstertruck_geometry::prelude::*;
use std::num::NonZeroUsize;

use super::error::FilletError;

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

fn valid_radius_value(value: f64) -> bool { value.is_finite() && value > 0.0 }

/// Runtime check for values returned by arbitrary radius functions between
/// pre-validation samples. Geometry helpers accept a plain f64 callback, so
/// record invalid evaluations and convert them to hard errors afterward.
#[derive(Default)]
pub struct RuntimeRadiusValidation(std::cell::Cell<bool>);

impl RuntimeRadiusValidation {
    pub fn observe(&self, radius: f64) -> f64 {
        if !valid_radius_value(radius) {
            self.0.set(true);
        }
        radius
    }

    pub const fn ensure_valid(&self) -> std::result::Result<(), FilletError> {
        if self.0.get() {
            Err(FilletError::InvalidRadius)
        } else {
            Ok(())
        }
    }
}

impl RadiusSpec {
    /// Validate the defined input domain before topology mutation and return
    /// an upper bound for the selected radii. Variable functions are arbitrary,
    /// so their representative domain samples are supplemented by checks at
    /// actual geometry evaluation sites.
    pub(crate) fn validated_upper_bound(&self) -> std::result::Result<f64, FilletError> {
        const VARIABLE_RADIUS_SAMPLES: usize = 32;
        let valid = |value: f64| {
            valid_radius_value(value)
                .then_some(value)
                .ok_or(FilletError::InvalidRadius)
        };
        match self {
            Self::Constant(radius) => valid(*radius),
            Self::PerEdge(radii) => radii.iter().copied().try_fold(0.0f64, |max, radius| {
                valid(radius).map(|radius| max.max(radius))
            }),
            Self::Variable(function) => (0..=VARIABLE_RADIUS_SAMPLES)
                .map(|i| {
                    function(crate::count_as_f64(i) / crate::count_as_f64(VARIABLE_RADIUS_SAMPLES))
                })
                .try_fold(0.0f64, |max, radius| {
                    valid(radius).map(|radius| max.max(radius))
                }),
        }
    }
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
    /// Validate the radius consumed by the selected profile. Asymmetric
    /// chamfer setbacks determine their own geometry and ignore the radius.
    pub(crate) fn validated_radius_upper_bound(&self) -> std::result::Result<f64, FilletError> {
        if self.profile.uses_radius() {
            self.radius.validated_upper_bound()
        } else {
            Ok(0.0)
        }
    }

    /// Creates options with a constant radius.
    #[must_use]
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
    #[must_use]
    pub fn with_radius(mut self, radius: RadiusSpec) -> Self {
        self.radius = radius;
        self
    }

    /// Sets the division count used by the fillet algorithm.
    #[must_use]
    pub const fn with_division(mut self, division: NonZeroUsize) -> Self {
        self.divisions = division;
        self
    }

    /// Sets the fillet profile.
    #[must_use]
    pub fn with_profile(mut self, profile: FilletProfile) -> Self {
        self.profile = profile;
        self
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
