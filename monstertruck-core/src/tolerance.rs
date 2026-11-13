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
    /// Default numerical length supplied to classic Boolean operations.
    /// Independent of the generic kernel zero threshold and of the maximum
    /// geometric-change budget selected by a caller.
    pub boolean_default_numeric_length: f64,
    /// Absolute convergence step in native curve parameter coordinates; not a spatial error budget.
    pub newton_curve_parameter_step: f64,
    /// Absolute convergence step in native surface UV coordinates; not a millimetre distance.
    pub newton_surface_parameter_step: f64,
    /// Absolute convergence step in mixed surface-UV plus curve parameter coordinates.
    pub newton_surface_curve_parameter_step: f64,
    /// Absolute convergence step in four surface UV coordinates during SSI.
    pub newton_intersection_parameter_step: f64,
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
    /// Unitless squared-normalization residual for analytic-curve support axes.
    pub analytic_curve_axis_unit: f64,
    /// Dimensionless leader segment extrapolation allowance when reconstructing SSI.
    pub boolean_leader_parameter_slack: f64,
    /// Dimensionless progress monotonicity guard along a leader polyline.
    pub boolean_leader_progress_slack: f64,
    /// Dimensionless overlap-parameter comparison between SSI fragments.
    pub boolean_fragment_overlap_parameter_slack: f64,
    /// Relative span padding for planar polygon overlap.
    pub boolean_planar_relative_span: f64,
    /// Dimensionless UV-loop area threshold relative to the loop's UV bounding box.
    pub boolean_uv_loop_area_relative: f64,
    /// UV-boundary contact distance relative to the combined UV extent.
    pub boolean_uv_boundary_contact_relative: f64,
    /// Axis-component cutoff for choosing a well-conditioned cylinder chart basis.
    pub boolean_cylinder_reference_axis_component_limit: f64,
    /// Fractional endpoint exclusion used by planar overlap sampling.
    pub boolean_planar_end_fraction: f64,
    /// Dimensionless ordered probe ratios for planar overlap disambiguation.
    pub boolean_planar_primary_probe_ratios: [f64; 4],
    /// Dimensionless ordered refinement ratios for planar overlap disambiguation.
    pub boolean_planar_secondary_probe_ratios: [f64; 2],
    /// Dimensionless allowed deviation of aligned planar carrier normals.
    pub boolean_planar_normal_cosine_slack: f64,
    /// Dimensionless allowed deviation of planar trim tangent alignment.
    pub boolean_planar_tangent_cosine_slack: f64,
    /// Native UV-coordinate tolerance used for tessellation projection diagnostics.
    pub meshing_debug_uv_edge_tolerance: f64,
    /// Native UV-coordinate threshold for STEP trim seam probe.
    pub step_uv_clamp_probe: f64,
    /// Millimetre residual for STEP trim curve reconstruction.
    pub step_trim_residual_mm: f64,
    /// STEP file uncertainty metadata (current model units).
    pub step_output_distance_accuracy: f64,
    /// Non-geometric STEP number formatting threshold, dimensionless.
    pub step_output_scientific_notation_threshold: f64,
    /// Scale multiplying Boolean numeric tolerance when joining near-fragments.
    pub boolean_fragment_stitch_scale: f64,
    /// Scale multiplying Boolean numerical tolerance in coplanar matching.
    pub boolean_coplanar_scale: f64,
    /// Classifier witness offset relative to explicit Boolean length tolerance.
    pub boolean_classifier_offset_factor: f64,
    /// Area comparison's absolute squared-length coefficient.
    pub boolean_area_floor_factor: f64,
    /// Area comparison's length-times-span coefficient.
    pub boolean_area_span_factor: f64,
    /// Minimum cosine between adjacent SSI fragment tangents before joining.
    pub boolean_fragment_forward_alignment_cosine: f64,
    /// Minimum absolute cosine between recovered cylinder normal and radial direction.
    pub boolean_cylinder_normal_alignment_cosine: f64,
    /// Ordered dimensionless inset fractions for finding a strict interior witness in thin planar regions.
    pub boolean_interior_witness_insets: [f64; 10],
    /// Bounded sample fractions for bidirectional whole-curve coincidence evidence.
    pub boolean_coincidence_fractions: [f64; 5],
    /// Bounded sample fractions for trimmed subedge containment and same-endpoint correspondence evidence.
    pub boolean_containment_fractions: [f64; 9],
    /// Offset multiplier for STL coordinate quantization.
    pub stl_quantization_offset_factor: f64,
    /// Step multiplier for STL coordinate quantization.
    pub stl_quantization_step_factor: f64,
}

/// Default numerical policy for the current kernel implementation.
pub const DEFAULT_NUMERICAL_TOLERANCES: NumericalTolerances = NumericalTolerances {
    kernel_f64: 1.0e-6,
    boolean_default_numeric_length: 1.0e-6,
    newton_curve_parameter_step: 1.0e-6,
    newton_surface_parameter_step: 1.0e-6,
    newton_surface_curve_parameter_step: 1.0e-6,
    newton_intersection_parameter_step: 1.0e-6,
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
    analytic_curve_axis_unit: 1.0e-8,
    boolean_leader_parameter_slack: 1.0e-8,
    boolean_leader_progress_slack: 1.0e-8,
    boolean_fragment_overlap_parameter_slack: 1.0e-5,
    boolean_planar_relative_span: 1.0e-10,
    boolean_uv_loop_area_relative: 1.0e-6,
    boolean_uv_boundary_contact_relative: 1.0e-6,
    boolean_cylinder_reference_axis_component_limit: 0.8,
    boolean_planar_end_fraction: 1.0e-6,
    boolean_planar_primary_probe_ratios: [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5],
    boolean_planar_secondary_probe_ratios: [1.0e-7, 1.0e-6],
    boolean_planar_normal_cosine_slack: 1.0e-10,
    boolean_planar_tangent_cosine_slack: 1.0e-8,
    meshing_debug_uv_edge_tolerance: 1.0e-9,
    step_uv_clamp_probe: 1.0e-3,
    step_trim_residual_mm: 1.0e-8,
    step_output_distance_accuracy: 1.0e-6,
    step_output_scientific_notation_threshold: 1.0e-2,
    boolean_fragment_stitch_scale: 4.0,
    boolean_coplanar_scale: 2.0,
    boolean_classifier_offset_factor: 8.0,
    boolean_area_floor_factor: 16.0,
    boolean_area_span_factor: 8.0,
    boolean_fragment_forward_alignment_cosine: 0.8,
    boolean_cylinder_normal_alignment_cosine: 0.9,
    boolean_interior_witness_insets: [
        0.1, 0.03, 0.01, 0.003, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12,
    ],
    boolean_coincidence_fractions: [0.0, 0.25, 0.5, 0.75, 1.0],
    boolean_containment_fractions: [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0],
    stl_quantization_offset_factor: 0.25,
    stl_quantization_step_factor: 0.5,
};

/// STEP-specific accuracy policy.
///
/// Unlike the unitless kernel thresholds, all distance fields here are in
/// millimetres after STEP unit normalization. Unitless and angular values are
/// explicitly labeled. These are measured source admissibility/recognition
/// thresholds; none grants permission to collapse unrelated topology.
///
/// Deliberately separate fields for different semantic decisions even where
/// the current numbers happen to coincide: future tuning must not silently
/// couple unrelated proof and recovery operations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepRecoveryTolerances {
    /// Numerical consistency of a reconstructed STEP support, not permission to heal a feature.
    /// Unit: mm.
    pub numerical_geometry_mm: f64,
    /// Maximum source-feature reconciliation distance, subject to separate topology checks.
    /// Unit: mm.
    pub feature_heal_mm: f64,
    /// Default permitted STEP CAD modeling-error budget.
    /// Unit: mm.
    pub cad_project_standard_mm: f64,
    /// Matching user-selected edge endpoints against their literal coordinates.
    /// Unit: mm.
    pub cad_edge_selector_mm: f64,
    /// Detecting revolution full turns and zero angular sweeps.
    /// Unit: radians.
    pub cad_axis_angle_rad: f64,
    /// Detecting vertices and edges along a revolution axis.
    /// Unit: mm.
    pub cad_axis_locus_mm: f64,
    /// Joining successive CAD profile segments.
    /// Unit: mm.
    pub cad_profile_join_mm: f64,
    /// Rejecting vanishing profile arc sweep.
    /// Unit: radians.
    pub cad_profile_angle_rad: f64,
    /// Rejecting zero or nonunit CAD edge selector normals.
    /// Unit: dimensionless.
    pub cad_normal_unit: f64,
    /// Surface recovery/recognition residual against input coordinates.
    /// Unit: mm.
    pub surface_geometry_mm: f64,
    /// Homogeneous NURBS weight comparison for recovery.
    /// Unit: dimensionless.
    pub surface_weight_relative: f64,
    /// Knot/parameter endpoint comparison for surface recognition.
    /// Unit: native surface parameter.
    pub surface_parameter: f64,
    /// Geometric residual in extrusion, loft, revolution, and sweep reconstruction.
    /// Unit: mm.
    pub solid_geometry_mm: f64,
    /// Direction cosines/orthogonality for extrusion and loft recognition.
    /// Unit: dimensionless.
    pub solid_direction: f64,
    /// Direction cosines in swept-profile recognition.
    /// Unit: dimensionless.
    pub sweep_direction: f64,
    /// Axis/orientation agreement for revolution recognition.
    /// Unit: dimensionless.
    pub revolution_direction: f64,
    /// Loft support normal agreement.
    /// Unit: dimensionless.
    pub loft_direction: f64,
    /// Coarse axis grouping to identify dominant extrusions.
    /// Unit: dimensionless.
    pub dominant_extrusion_direction: f64,
    /// Direction agreement for edge-treatment recognition.
    /// Unit: dimensionless.
    pub edge_treatment_direction: f64,
    /// Direction agreement in inverse fillet/chamfer identification.
    /// Unit: dimensionless.
    pub inverse_edge_treatment_direction: f64,
    /// Direction agreement in reduced boundary representation.
    /// Unit: dimensionless.
    pub reduced_brep_direction: f64,
    /// Thickness and overlap matching for sheet-metal recovery.
    /// Unit: mm.
    pub sheet_geometry_mm: f64,
    /// Direction agreement for formed sheet faces and bends.
    /// Unit: dimensionless.
    pub sheet_direction: f64,
    /// Tighter formed-sheet normal and axis agreement.
    /// Unit: dimensionless.
    pub formed_sheet_direction: f64,
    /// Degenerate-angle threshold in extrusion inference.
    /// Unit: radians.
    pub solid_extrusion_angle_rad: f64,
    /// Spatial residual for periodic chains and periodic bodies.
    /// Unit: mm.
    pub periodic_geometry_mm: f64,
    /// Coordinate matching while resizing periodic features.
    /// Unit: mm.
    pub periodic_resize_mm: f64,
    /// Pattern center and lattice matching during feature recognition.
    /// Unit: mm.
    pub pattern_geometry_mm: f64,
    /// Pattern basis orientation and angular consistency.
    /// Unit: dimensionless.
    pub pattern_orientation: f64,
    /// Quantization of feature partition vertices.
    /// Unit: mm.
    pub partition_point_quant_mm: f64,
    /// Quantization of partition-normal directions.
    /// Unit: dimensionless.
    pub partition_normal_quant: f64,
    /// Position matching when interning equivalent supports.
    /// Unit: mm.
    pub geometric_intern_position_mm: f64,
    /// Scalar offset matching when interning equivalent supports.
    /// Unit: mm.
    pub geometric_intern_scalar_mm: f64,
    /// Normal and direction matching when interning supports.
    /// Unit: dimensionless.
    pub geometric_intern_direction: f64,
    /// Position quantization for independent boundary proof.
    /// Unit: mm.
    pub boundary_point_quant_mm: f64,
    /// Scalar offsets quantization in boundary proof.
    /// Unit: mm.
    pub boundary_scalar_quant_mm: f64,
    /// Axis/normal quantization in independent boundary proof.
    /// Unit: dimensionless.
    pub boundary_direction_quant: f64,
    /// Factor converting boundary scalar quantum to admitted support locus residual.
    /// Unit: multiplier.
    pub boundary_locus_scale: f64,
    /// Factor converting boundary direction quantum to proof normal residual.
    /// Unit: multiplier.
    pub boundary_angular_scale: f64,
    /// Line segment angular collinearity comparison.
    /// Unit: dimensionless.
    pub line_collinear: f64,
    /// Ordering slack in line recovery.
    /// Unit: dimensionless.
    pub line_monotonic: f64,
    /// Direction discretization for line recovery.
    /// Unit: dimensionless.
    pub line_direction_quant: f64,
    /// Endpoint direction consistency in line recovery.
    /// Unit: dimensionless.
    pub line_endpoint_direction: f64,
    /// Squared chord length for degenerate recovered lines.
    /// Unit: mm².
    pub line_degenerate_chord_squared: f64,
    /// Very small vector accepted for source BREP direction normalization.
    /// Unit: dimensionless.
    pub source_brep_direction: f64,
    /// Geometric agreement of repeated source curves.
    /// Unit: mm.
    pub curve_replica_geometry_mm: f64,
    /// Pure translation residual for exact curve replication.
    /// Unit: mm.
    pub curve_replica_transform_mm: f64,
    /// Pure translation residual for exact surface replication.
    /// Unit: mm.
    pub surface_replica_transform_mm: f64,
    /// Side and orientation checks in planar feature recognition.
    /// Unit: dimensionless.
    pub planar_feature_side: f64,
    /// Planar feature candidate residual.
    /// Unit: mm.
    pub planar_feature_fit_mm: f64,
    /// Upper source-data uncertainty for revolution reconstruction.
    /// Unit: mm.
    pub revolution_max_source_uncertainty_mm: f64,
    /// Axis grouping criterion in compound CAD recovery.
    /// Unit: dimensionless.
    pub axis_group_dot: f64,
    /// Station clustering in cap/loft CAD reconstruction.
    /// Unit: mm.
    pub station_position_mm: f64,
    /// Frame orthonormality for reconstructed source cap.
    /// Unit: dimensionless.
    pub source_cap_frame: f64,
    /// Orthonormality for reconstructed source sweep profiles.
    /// Unit: dimensionless.
    pub source_sweep_frame: f64,
    /// Strict original sweep support fitting error.
    /// Unit: mm.
    pub source_sweep_strict_residual_mm: f64,
    /// Larger source sweep support fitting error.
    /// Unit: mm.
    pub source_sweep_loose_residual_mm: f64,
    /// Tie breaker between alternative pattern basis fits.
    /// Unit: dimensionless cost.
    pub pattern_basis_cost_delta: f64,
    /// Minimum position quantization for feature pattern candidates.
    /// Unit: mm.
    pub pattern_quantization_floor: f64,
    /// Reject nearly degenerate two-dimensional pattern bases.
    /// Unit: mm⁴.
    pub pattern_determinant_squared: f64,
    /// Denominator floor for pattern coordinates.
    /// Unit: dimensionless.
    pub pattern_lattice_normal_floor: f64,
    /// Lattice fractional-coordinate quantization floor.
    /// Unit: dimensionless.
    pub pattern_fraction_quant_floor: f64,
    /// Minimum origin matching scale for pattern transforms.
    /// Unit: mm.
    pub pattern_origin_floor_mm: f64,
    /// Tolerance for pattern filling-ratio ranking ties.
    /// Unit: dimensionless.
    pub pattern_comparison_residual: f64,
    /// Minimum vector norm accepted for pattern direction normalization.
    /// Unit: dimensionless.
    pub pattern_normalization_min: f64,
    /// Near-zero geometric chord threshold for source support extraction.
    /// Unit: mm.
    pub source_straightness_length_mm: f64,
    /// Acceptable axis norm deviation in source feature frames.
    /// Unit: dimensionless.
    pub source_axis_unit_cosine: f64,
    /// Axis-normal agreement for source features.
    /// Unit: dimensionless.
    pub source_axis_direction_cosine: f64,
    /// Minimal accepted bend angle and supplement for source sweeps.
    /// Unit: radians.
    pub source_sweep_angle_rad: f64,
    /// Tangent and profile orientation consistency in source sweep recovery.
    /// Unit: dimensionless.
    pub source_sweep_tangent_cosine: f64,
    /// Progress and joined segment angular residual in sweep recovery.
    /// Unit: radians.
    pub source_sweep_parameter_slack: f64,
    /// Parameter-comparison slack during source sweep recovery.
    /// Unit: native curve parameter.
    pub source_sweep_curve_parameter: f64,
    /// Knot equality and source support weight matching.
    /// Unit: native knot/weight parameter.
    pub source_sweep_knot_residual: f64,
    /// Reject near-zero B-spline homogeneous support weight.
    /// Unit: native homogeneous weight.
    pub source_sweep_small_weight: f64,
    /// Minimal valid sheet bend angle.
    /// Unit: radians.
    pub source_sheet_angle_rad: f64,
    /// Sheet face normal alignment and near-parallel classification.
    /// Unit: dimensionless.
    pub source_sheet_normal_cosine: f64,
    /// Sheet seam and tangent alignment.
    /// Unit: dimensionless.
    pub source_sheet_tangent_cosine: f64,
    /// Sheet bend target-angle fitting.
    /// Unit: radians.
    pub source_sheet_bend_angle_rad: f64,
    /// Squared size below which source direction/area evidence is degenerate.
    /// Unit: mm².
    pub source_identity_zero_norm_squared: f64,
    /// Source coordinate frame orthonormality threshold.
    /// Unit: dimensionless.
    pub source_identity_axis_unit: f64,
    /// Zero-motion matching in source rigid transformations.
    /// Unit: mm.
    pub source_identity_transform: f64,
    /// Reciprocal position quantization step for source identity keys.
    /// Unit: mm^-1.
    pub source_identity_quantization_inv_mm: f64,
    /// Inverse normal quantization step for source identity keys.
    /// Unit: dimensionless.
    pub source_identity_normal_quantization: f64,
    /// Point-on-original-source support recovery check.
    /// Unit: mm.
    pub source_feature_incidence_mm: f64,
    /// Analytic surface candidate fitting in boundary proof.
    /// Unit: mm.
    pub source_analytic_fit_mm: f64,
    /// Finite parameter sampling fraction for boundary tangent proof.
    /// Unit: dimensionless.
    pub source_boundary_tangent_fraction: f64,
    /// Minimum valid source meridian-angle comparison.
    /// Unit: radians.
    pub source_revolution_angle_floor_rad: f64,
    /// Lower bound for parameter-space arc geometry slop.
    /// Unit: native curve parameter.
    pub source_revolution_param_floor: f64,
    /// Orthonormality of source CAD frames and transformed vectors.
    /// Unit: dimensionless.
    pub source_cad_frame_orthogonality: f64,
    /// Allowed angular mismatch in reverse fillet arc recovery.
    /// Unit: radians.
    pub source_rounding_arc_residual_rad: f64,
    /// Allowed angular excess beyond a full rounding sweep.
    /// Unit: radians.
    pub source_rounding_arc_max_slack_rad: f64,
    /// Bound on source rounding parameter interval comparisons.
    /// Unit: native curve parameter.
    pub source_rounding_parameter_slack: f64,
    /// Allowed slack when decomposing corner arcs into quarter turns.
    /// Unit: radians.
    pub source_rounding_quarter_turn_slack_rad: f64,
    /// Straight source curve fit against expected 3D locus.
    /// Unit: mm.
    pub source_line_residual_mm: f64,
    /// Point matching while reconstructing source cap boundary.
    /// Unit: mm.
    pub source_cap_point_match_mm: f64,
    /// Direction agreement in source cap and circle sweep frames.
    /// Unit: dimensionless.
    pub source_cap_direction_cosine: f64,
    /// Source cap and sphere common-support matching.
    /// Unit: mm.
    pub source_cap_rounding_residual_mm: f64,
    /// Relative radius offset bound in source support interning.
    /// Unit: dimensionless.
    pub source_intern_radius_scale_relative: f64,
    /// Direction quantization for periodic source feature chains.
    /// Unit: dimensionless.
    pub source_periodic_direction_quant: f64,
    /// Reject near-degenerate pattern-recovery support vectors.
    /// Unit: mm².
    pub pattern_inverse_norm_squared: f64,
    /// Maximum acceptable fitted sweep support residual.
    /// Unit: mm.
    pub source_sweep_absolute_residual_mm: f64,
    /// Scale-dependent sweep support comparison.
    /// Unit: dimensionless.
    pub source_sweep_relative_compare: f64,
    /// Scale-dependent sheet geometry comparison.
    /// Unit: dimensionless.
    pub source_sheet_relative_compare: f64,
    /// Inverse coordinate quantization for periodic graph weld keys.
    /// Unit: model-coordinate^-1.
    pub source_periodic_weld_quantization_inv: f64,
    /// Minimum nondegenerate reference area for periodic graph editing.
    /// Unit: mm².
    pub source_periodic_area_min: f64,
    /// Tiny allowed source-fitted sweep tolerance budget overrun.
    /// Unit: mm.
    pub source_sweep_tolerance_slack_mm: f64,
    /// Position check for cap-slot boundary recovery.
    /// Unit: mm.
    pub source_slot_boundary_mm: f64,
    /// Arc sweep endpoint slop in CAD feature reconstruction.
    /// Unit: radians.
    pub source_cad_arc_angle_slack_rad: f64,
}

/// Current STEP reconstruction tolerances; every field has one owning definition.
pub const DEFAULT_STEP_RECOVERY_TOLERANCES: StepRecoveryTolerances = StepRecoveryTolerances {
    numerical_geometry_mm: 1.0e-7,
    feature_heal_mm: 0.01,
    cad_project_standard_mm: 0.01,
    cad_edge_selector_mm: 1.0e-7,
    cad_axis_angle_rad: 1.0e-10,
    cad_axis_locus_mm: 1.0e-9,
    cad_profile_join_mm: 1.0e-8,
    cad_profile_angle_rad: 1.0e-12,
    cad_normal_unit: 1.0e-12,
    surface_geometry_mm: 1.0e-5,
    surface_weight_relative: 1.0e-12,
    surface_parameter: 1.0e-12,
    solid_geometry_mm: 1.0e-7,
    solid_direction: 1.0e-10,
    sweep_direction: 1.0e-9,
    revolution_direction: 1.0e-9,
    loft_direction: 1.0e-10,
    dominant_extrusion_direction: 1.0e-8,
    edge_treatment_direction: 1.0e-8,
    inverse_edge_treatment_direction: 1.0e-10,
    reduced_brep_direction: 1.0e-10,
    sheet_geometry_mm: 1.0e-5,
    sheet_direction: 1.0e-9,
    formed_sheet_direction: 1.0e-10,
    solid_extrusion_angle_rad: 1.0e-10,
    periodic_geometry_mm: 1.0e-5,
    periodic_resize_mm: 1.0e-7,
    pattern_geometry_mm: 1.0e-7,
    pattern_orientation: 1.0e-10,
    partition_point_quant_mm: 1.0e-8,
    partition_normal_quant: 1.0e-10,
    geometric_intern_position_mm: 1.0e-5,
    geometric_intern_scalar_mm: 1.0e-5,
    geometric_intern_direction: 1.0e-10,
    boundary_point_quant_mm: 1.0e-7,
    boundary_scalar_quant_mm: 1.0e-7,
    boundary_direction_quant: 1.0e-10,
    boundary_locus_scale: 4.0,
    boundary_angular_scale: 4.0,
    line_collinear: 1.0e-12,
    line_monotonic: 1.0e-12,
    line_direction_quant: 1.0e-12,
    line_endpoint_direction: 2.0e-13,
    line_degenerate_chord_squared: 1.0e-24,
    source_brep_direction: 1.0e-15,
    curve_replica_geometry_mm: 1.0e-5,
    curve_replica_transform_mm: 1.0e-12,
    surface_replica_transform_mm: 1.0e-12,
    planar_feature_side: 1.0e-8,
    planar_feature_fit_mm: 1.0e-7,
    revolution_max_source_uncertainty_mm: 1.0e-5,
    axis_group_dot: 1.0e-10,
    station_position_mm: 1.0e-7,
    source_cap_frame: 1.0e-10,
    source_sweep_frame: 1.0e-8,
    source_sweep_strict_residual_mm: 1.0e-6,
    source_sweep_loose_residual_mm: 1.0e-5,
    pattern_basis_cost_delta: 1.0e-12,
    pattern_quantization_floor: 1.0e-9,
    pattern_determinant_squared: 1.0e-18,
    pattern_lattice_normal_floor: 1.0e-12,
    pattern_fraction_quant_floor: 1.0e-10,
    pattern_origin_floor_mm: 1.0e-9,
    pattern_comparison_residual: 1.0e-12,
    pattern_normalization_min: 1.0e-15,
    source_straightness_length_mm: 1.0e-12,
    source_axis_unit_cosine: 1.0e-9,
    source_axis_direction_cosine: 1.0e-8,
    source_sweep_angle_rad: 1.0e-8,
    source_sweep_tangent_cosine: 1.0e-7,
    source_sweep_parameter_slack: 1.0e-7,
    source_sweep_curve_parameter: 1.0e-10,
    source_sweep_knot_residual: 1.0e-12,
    source_sweep_small_weight: 1.0e-14,
    source_sheet_angle_rad: 1.0e-8,
    source_sheet_normal_cosine: 1.0e-7,
    source_sheet_tangent_cosine: 1.0e-8,
    source_sheet_bend_angle_rad: 1.0e-7,
    source_identity_zero_norm_squared: 1.0e-24,
    source_identity_axis_unit: 1.0e-9,
    source_identity_transform: 1.0e-12,
    source_identity_quantization_inv_mm: 1.0e9,
    source_identity_normal_quantization: 1.0e5,
    source_feature_incidence_mm: 1.0e-7,
    source_analytic_fit_mm: 1.0e-7,
    source_boundary_tangent_fraction: 1.0e-5,
    source_revolution_angle_floor_rad: 1.0e-12,
    source_revolution_param_floor: 1.0e-12,
    source_cad_frame_orthogonality: 1.0e-10,
    source_rounding_arc_residual_rad: 1.0e-6,
    source_rounding_arc_max_slack_rad: 1.0e-7,
    source_rounding_parameter_slack: 1.0e-9,
    source_rounding_quarter_turn_slack_rad: 1.0e-10,
    source_line_residual_mm: 1.0e-7,
    source_cap_point_match_mm: 1.0e-8,
    source_cap_direction_cosine: 1.0e-8,
    source_cap_rounding_residual_mm: 1.0e-12,
    source_intern_radius_scale_relative: 1.0e-10,
    source_periodic_direction_quant: 1.0e-10,
    pattern_inverse_norm_squared: 1.0e-24,
    source_sweep_absolute_residual_mm: 1.0e-5,
    source_sweep_relative_compare: 1.0e-9,
    source_sheet_relative_compare: 1.0e-9,
    source_periodic_weld_quantization_inv: 1.0e12,
    source_periodic_area_min: 1.0e-12,
    source_sweep_tolerance_slack_mm: 1.0e-15,
    source_slot_boundary_mm: 1.0e-7,
    source_cad_arc_angle_slack_rad: 1.0e-9,
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
