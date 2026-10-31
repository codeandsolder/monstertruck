use crate::cgmath64::{Zero, cgmath};
use cgmath::AbsDiffEq;
use std::fmt::Debug;

/// Central numerical tolerances used by the kernel implementation.
///
/// These are robustness and representation thresholds, not a modeling error
/// budget. Operations that are allowed to change geometry must take an
/// explicit caller tolerance instead of silently widening these values.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NumericalTolerances {
    /// Scalar epsilon for the established `f64` kernel.
    pub kernel_f64: f64,
    /// Scalar epsilon for explicitly `f32` geometry.
    pub kernel_f32: f32,
    /// Threshold for structural/algebraic zero tests that must be much tighter
    /// than ordinary geometric coincidence.
    pub structural: f64,
    /// Dimensionless angular/direction threshold.
    pub angular: f64,
    /// Relative parameter-space endpoint tolerance used while preserving trims.
    pub trim_parameter_relative: f64,
    /// Relative parameter-space tolerance used while splitting curves.
    pub split_parameter_relative: f64,
    /// Absolute floor paired with `split_parameter_relative`.
    pub split_parameter_absolute: f64,
    /// Tight parameter-space projection threshold.
    pub parameter_projection: f64,
    /// Relative padding for finite parameter ranges derived from external data.
    pub parameter_range_relative: f64,
    /// Relative tolerance for periodic-parameter branch comparisons.
    pub periodic_relative: f64,
    /// Relative threshold for ordinary rigid/similarity-transform recognition.
    pub transform_similarity: f64,
    /// Stricter similarity threshold used by proof-bearing same-domain logic.
    pub same_domain_similarity: f64,
    /// Homogeneous extrusion identity threshold.
    pub analytic_extrusion_identity: f64,
    /// Planar analytic-carrier identity threshold.
    pub analytic_plane_identity: f64,
    /// Spherical/revolution analytic-recognition threshold.
    pub analytic_spherical_revolution: f64,
    /// Relative threshold for rejecting nearly parallel surface normals.
    pub transversality_relative: f64,
    /// Coarsest numerical mesh seed used by algorithms that need a seed mesh.
    pub mesh_seed_floor: f64,
    /// Relative expansion used to construct finite half-space clipping bounds.
    pub half_space_relative_margin: f64,
    /// UV closure threshold for already-sampled tessellation boundaries.
    pub mesh_uv_closure: f64,
    /// Spatial distance used by tessellation to rescue samples at analytic poles.
    pub mesh_singularity_distance: f64,
    /// Squared UV-distance threshold inherited by closed-face healing.
    pub healing_boundary_closure_squared: f64,
    /// Parameter-space probe step used to recover a normal at a singular NURBS point.
    pub normal_probe_step: f64,
    /// Residual accepted by derivative consistency assertions.
    pub derivative_residual: f64,
    /// Multiplier from requested tessellation tolerance to surface-span segment length.
    pub surface_span_factor: f64,
    /// Default closure threshold for text-outline reconstruction.
    pub text_closure: f64,
    /// Central finite-difference step for fallback derivative evaluation.
    pub finite_difference_step: f64,
    /// Minimum sampling tolerance as a multiple of the kernel epsilon in healing.
    pub healing_sampling_floor_factor: f64,
    /// Multiplier used for best-effort STEP missing-edge recovery.
    pub recovery_geometry_factor: f64,
    /// Maximum continuation jump relative to the guide step in missing-edge recovery.
    pub recovery_continuation_jump_factor: f64,
    /// Kernel-epsilon floor used to choose stable radial evidence in same-domain proofs.
    pub same_domain_direction_floor_factor: f64,
    /// Kernel-precision envelope for cross-face stitching.
    pub cross_face_stitch_factor: f64,
    /// Projection padding around Boolean weld candidates.
    pub boolean_weld_projection_factor: f64,
    /// Residual envelope for Boolean surface projection evidence.
    pub boolean_surface_projection_factor: f64,
    /// Endpoint envelope for recovered Boolean iso-curves.
    pub boolean_curve_endpoint_factor: f64,
    /// Kernel-epsilon cell width multiplier for Boolean polyline hashing.
    pub boolean_polyline_grid_factor: f64,
    /// Kernel-epsilon floor multiplier for finite half-space bounds.
    pub boolean_half_space_kernel_factor: f64,
    /// Kernel-epsilon envelope for stitched Boolean intersection fragments.
    pub boolean_fragment_validation_factor: f64,
    /// Kernel-epsilon envelope for Boolean boundary parameter search.
    pub boolean_boundary_search_factor: f64,
    /// Kernel-epsilon cell width multiplier used by VTK point hashing.
    pub mesh_point_hash_factor: f64,
    /// Offset multiplier for STL coordinate quantization.
    pub stl_quantization_offset_factor: f64,
    /// Step multiplier for STL coordinate quantization.
    pub stl_quantization_step_factor: f64,
}

/// Default numerical policy for the current kernel implementation.
pub const DEFAULT_NUMERICAL_TOLERANCES: NumericalTolerances = NumericalTolerances {
    kernel_f64: 1.0e-6,
    kernel_f32: 1.0e-4,
    structural: 1.0e-12,
    angular: 1.0e-6,
    trim_parameter_relative: 1.0e-6,
    split_parameter_relative: 1.0e-8,
    split_parameter_absolute: 1.0e-10,
    parameter_projection: 1.0e-9,
    parameter_range_relative: 1.0e-9,
    periodic_relative: 1.0e-6,
    transform_similarity: 1.0e-9,
    same_domain_similarity: 1.0e-10,
    analytic_extrusion_identity: 1.0e-8,
    analytic_plane_identity: 1.0e-6,
    analytic_spherical_revolution: 1.0e-7,
    transversality_relative: 1.0e-8,
    mesh_seed_floor: 1.0e-3,
    half_space_relative_margin: 1.0e-3,
    mesh_uv_closure: 1.0e-3,
    mesh_singularity_distance: 1.0e-3,
    healing_boundary_closure_squared: 1.0e-3,
    normal_probe_step: 1.0e-4,
    derivative_residual: 1.0e-4,
    surface_span_factor: 8.0,
    text_closure: 1.0e-7,
    finite_difference_step: 1.0e-6,
    healing_sampling_floor_factor: 10.0,
    recovery_geometry_factor: 100.0,
    recovery_continuation_jump_factor: 8.0,
    same_domain_direction_floor_factor: 16.0,
    cross_face_stitch_factor: 8.0,
    boolean_weld_projection_factor: 2.0,
    boolean_surface_projection_factor: 4.0,
    boolean_curve_endpoint_factor: 8.0,
    boolean_polyline_grid_factor: 2.0,
    boolean_half_space_kernel_factor: 100.0,
    boolean_fragment_validation_factor: 4.0,
    boolean_boundary_search_factor: 8.0,
    mesh_point_hash_factor: 50.0,
    stl_quantization_offset_factor: 0.25,
    stl_quantization_step_factor: 0.5,
};

/// Numerical epsilon for ordinary `f64` kernel predicates.
pub const KERNEL_EPSILON: f64 = DEFAULT_NUMERICAL_TOLERANCES.kernel_f64;

/// Squared `f64` kernel epsilon.
pub const KERNEL_EPSILON_SQUARED: f64 = KERNEL_EPSILON * KERNEL_EPSILON;

/// Numerical epsilon for explicitly `f32` kernel predicates.
pub const KERNEL_EPSILON_F32: f32 = DEFAULT_NUMERICAL_TOLERANCES.kernel_f32;

/// Defines a tolerance in the whole package
pub trait Tolerance: AbsDiffEq<Epsilon = f64> + Debug {
    /// The "distance" is less than `KERNEL_EPSILON`.
    fn near(&self, other: &Self) -> bool { self.abs_diff_eq(other, KERNEL_EPSILON) }

    /// The "distance" is less than [`KERNEL_EPSILON_SQUARED`].
    fn near2(&self, other: &Self) -> bool { self.abs_diff_eq(other, KERNEL_EPSILON_SQUARED) }
}

impl<T: AbsDiffEq<Epsilon = f64> + Debug> Tolerance for T {}

/// Asserts that `left.near(&right)` (using `Tolerance`).
#[macro_export]
macro_rules! assert_near {
    ($left: expr, $right: expr $(,)?) => {{
        let (left, right) = ($left, $right);
        assert!(
            $crate::tolerance::Tolerance::near(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}",
        )
    }};
    ($left: expr, $right: expr, $($arg: tt)+) => {{
        let (left, right) = ($left, $right);
        assert!(
            $crate::tolerance::Tolerance::near(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}: {}",
            format_args!($($arg)+),
        )
    }};
}

/// Similar to `assert_near!`, but returns a test failure instead of panicking if the condition fails.
#[macro_export]
macro_rules! prop_assert_near {
    ($left: expr, $right: expr $(,)?) => {{
        let (left, right) = ($left, $right);
        prop_assert!(
            $crate::tolerance::Tolerance::near(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}",
        )
    }};
    ($left: expr, $right: expr, $($arg: tt)+) => {{
        let (left, right) = ($left, $right);
        prop_assert!(
            $crate::tolerance::Tolerance::near(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?}, right: {right:?}: {}",
            format_args!($($arg)+),
        )
    }};
}

const _: () = {
    assert!(KERNEL_EPSILON_SQUARED == KERNEL_EPSILON * KERNEL_EPSILON);
    assert!(DEFAULT_NUMERICAL_TOLERANCES.structural < KERNEL_EPSILON);
    assert!((KERNEL_EPSILON_F32 as f64) > KERNEL_EPSILON);
    assert!(DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor >= KERNEL_EPSILON);
    assert!(DEFAULT_NUMERICAL_TOLERANCES.normal_probe_step >= KERNEL_EPSILON);
};

#[test]
#[should_panic(expected = "assertion failed: `left` is near `right`")]
fn assert_near_without_msg() { assert_near!(1.0, 2.0) }

#[test]
#[should_panic(expected = "assertion failed: `left` is near `right`")]
fn assert_near_with_msg() { assert_near!(1.0, 2.0, "{}", "test OK") }

/// Asserts that `left.near2(&right)` (using `Tolerance`).
#[macro_export]
macro_rules! assert_near2 {
    ($left: expr, $right: expr $(,)?) => {{
        let (left, right) = ($left, $right);
        assert!(
            $crate::tolerance::Tolerance::near2(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}",
        )
    }};
    ($left: expr, $right: expr, $($arg: tt)+) => {{
        let (left, right) = ($left, $right);
        assert!(
            $crate::tolerance::Tolerance::near2(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}: {}",
            format_args!($($arg)+),
        )
    }};
}

/// Similar to `assert_near2!`, but returns a test failure instead of panicking if the condition fails.
#[macro_export]
macro_rules! prop_assert_near2 {
    ($left: expr, $right: expr $(,)?) => {{
        let (left, right) = ($left, $right);
        prop_assert!(
            $crate::tolerance::Tolerance::near2(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}",
        )
    }};
    ($left: expr, $right: expr, $($arg: tt)+) => {
        let (left, right) = ($left, $right);
        prop_assert!(
            $crate::tolerance::Tolerance::near2(&left, &right),
            "assertion failed: `left` is near `right`\nleft: {left:?},\nright: {right:?}: {}",
            format_args!($($arg)+),
        )
    };
}

#[test]
#[should_panic(expected = "assertion failed: `left` is near `right`")]
fn assert_near2_without_msg() { assert_near2!(1.0, 2.0) }

#[test]
#[should_panic(expected = "assertion failed: `left` is near `right`")]
fn assert_near2_with_msg() { assert_near2!(1.0, 2.0, "{}", "test OK") }

/// The structs defined the origin. `f64`, `Vector`, and so on.
pub trait Origin: Tolerance + Zero {
    /// near origin
    #[inline]
    fn so_small(&self) -> bool { self.near(&Self::zero()) }

    /// near origin in square order
    #[inline]
    fn so_small2(&self) -> bool { self.near2(&Self::zero()) }
}

impl<T: Tolerance + Zero> Origin for T {}
