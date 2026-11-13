use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

/// Dimension-aware tolerances for the classic Boolean backend.
///
/// `mesh_seed` is a discretization choice and may be deliberately coarser than
/// the requested Boolean numerical tolerance. It must never be reused as proof of
/// coplanarity, equality, containment, or topology identity.
#[derive(Clone, Copy, Debug)]
pub(super) struct BooleanTolerance {
    pub(super) length: f64,
    pub(super) area_floor: f64,
    pub(super) angular: f64,
    pub(super) classifier_offset: f64,
    pub(super) mesh_seed: f64,
}

impl BooleanTolerance {
    pub(super) fn new(requested: f64) -> Self {
        // The caller's positive length tolerance is authoritative: a larger
        // global kernel resolution must not silently widen geometric identity.
        assert!(
            requested.is_finite() && requested > 0.0,
            "Boolean length tolerance must be finite and positive"
        );
        let length = requested;
        Self {
            length,
            area_floor: length * length,
            angular: DEFAULT_NUMERICAL_TOLERANCES.angular,
            classifier_offset: 8.0 * length,
            mesh_seed: requested.max(DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor),
        }
    }

    pub(super) fn area_for_span(self, span: f64) -> f64 {
        16.0f64.mul_add(
            self.area_floor,
            8.0 * self.length * span.abs().max(self.length),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn semantic_tolerances_do_not_inherit_mesh_floor() {
        let tol = BooleanTolerance::new(1.0e-6);
        assert_eq!(tol.length, 1.0e-6);
        assert_eq!(tol.classifier_offset, 8.0e-6);
        assert_eq!(tol.mesh_seed, 1.0e-3);
    }

    #[test]
    fn area_tolerance_has_area_units_and_scales_with_extent() {
        let tol = BooleanTolerance::new(1.0e-5);
        assert!(tol.area_for_span(10.0) > tol.area_for_span(1.0));
        assert!(tol.area_for_span(1.0) > tol.area_floor);
    }
}
