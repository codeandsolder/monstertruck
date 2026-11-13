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
        // EXPERIMENT ONLY: isolate mesh-discretization granularity from
        // geometric admission tolerances and curve adaptive chord sampling.
        static EXPERIMENT_MESH_FLOOR: std::sync::OnceLock<f64> = std::sync::OnceLock::new();
        let mesh_floor = *EXPERIMENT_MESH_FLOOR.get_or_init(|| {
            std::env::var("MT_EXPERIMENT_BOOLEAN_MESH_FLOOR_MM")
                .ok()
                .map(|text| {
                    let value: f64 = text
                        .parse()
                        .expect("MT_EXPERIMENT_BOOLEAN_MESH_FLOOR_MM must be a number");
                    assert!(
                        value.is_finite() && value > 0.0,
                        "MT_EXPERIMENT_BOOLEAN_MESH_FLOOR_MM must be finite and positive"
                    );
                    value
                })
                .unwrap_or(0.0)
        });
        Self {
            length,
            area_floor: length * length,
            angular: DEFAULT_NUMERICAL_TOLERANCES.angular,
            classifier_offset: DEFAULT_NUMERICAL_TOLERANCES.boolean_classifier_offset_factor
                * length,
            mesh_seed: requested
                .max(DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor)
                .max(mesh_floor),
        }
    }

    /// Maximum allowed separation between independently corrected points on
    /// the two exact SSI carrier surfaces, in model-coordinate length units.
    /// Not a mesh discretization allowance.
    pub(super) fn surface_support_residual(self) -> f64 {
        DEFAULT_NUMERICAL_TOLERANCES.boolean_surface_projection_factor * self.length
    }

    /// Maximum endpoint/branch discrepancy allowed when validating an
    /// intersection fragment against exact carrier geometry.
    pub(super) fn fragment_residual(self) -> f64 {
        DEFAULT_NUMERICAL_TOLERANCES.boolean_fragment_validation_factor * self.length
    }

    /// Exact curve/surface endpoint admission, independent of the seed mesh.
    pub(super) fn boundary_residual(self) -> f64 {
        DEFAULT_NUMERICAL_TOLERANCES.boolean_boundary_search_factor * self.length
    }

    /// Cross-face endpoint stitching radius. This is a numerical
    /// canonicalization tolerance, never the permitted modeling-error budget.
    pub(super) fn cross_face_stitch_residual(self) -> f64 {
        DEFAULT_NUMERICAL_TOLERANCES.cross_face_stitch_factor * self.length
    }

    pub(super) fn area_for_span(self, span: f64) -> f64 {
        DEFAULT_NUMERICAL_TOLERANCES
            .boolean_area_floor_factor
            .mul_add(
                self.area_floor,
                DEFAULT_NUMERICAL_TOLERANCES.boolean_area_span_factor * self.length * span.abs(),
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
    fn sub_kernel_epsilon_length_is_not_widened_to_kernel_epsilon() {
        let requested = 1.0e-8;
        let tolerance = BooleanTolerance::new(requested);
        assert_eq!(tolerance.length, requested);
        assert_eq!(tolerance.area_floor, requested * requested);
        assert_eq!(
            tolerance.classifier_offset,
            DEFAULT_NUMERICAL_TOLERANCES.boolean_classifier_offset_factor * requested
        );
        assert_eq!(
            tolerance.mesh_seed, DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor,
            "only the seed mesh gets the independent, coarser sampling floor"
        );
    }

    #[test]
    #[should_panic(expected = "Boolean length tolerance must be finite and positive")]
    fn zero_length_tolerance_is_not_silently_widened() { let _ = BooleanTolerance::new(0.0); }

    #[test]
    fn exact_ssi_residuals_scale_with_numerical_length_not_global_epsilon() {
        let tighter = BooleanTolerance::new(1.0e-8);
        let base = BooleanTolerance::new(1.0e-6);
        assert_eq!(tighter.surface_support_residual(), 4.0e-8);
        assert_eq!(tighter.fragment_residual(), 4.0e-8);
        assert_eq!(base.surface_support_residual(), 4.0e-6);
        assert_eq!(base.fragment_residual(), 4.0e-6);
        assert_eq!(tighter.mesh_seed, base.mesh_seed);
    }

    #[test]
    fn boundary_and_cross_face_precision_scales_with_operation() {
        let precise = BooleanTolerance::new(1.0e-8);
        assert_eq!(precise.boundary_residual(), 8.0e-8);
        assert_eq!(precise.cross_face_stitch_residual(), 8.0e-8);
        let base = BooleanTolerance::new(1.0e-6);
        assert_eq!(base.boundary_residual(), 8.0e-6);
        assert_eq!(base.cross_face_stitch_residual(), 8.0e-6);
        assert_eq!(precise.mesh_seed, base.mesh_seed);
    }

    #[test]
    fn area_tolerance_has_area_units_and_scales_with_extent() {
        let tol = BooleanTolerance::new(1.0e-5);
        assert!(tol.area_for_span(10.0) > tol.area_for_span(1.0));
        assert!(tol.area_for_span(1.0) > tol.area_floor);
    }
}
