# Geometry tolerance policy -- authoritative inventory and reconstruction targets

**Scope:** Monstertruck kernel, Boolean/intersection/meshing, STEP reader, and
`step-redox` feature reconstruction. The **only definition of numerical
values** in this inventory is `monstertruck-core/src/tolerance.rs`.
The following tables are a generated *description* of those definitions,
not a second set of constants. If a field changes, refresh these tables.
`KERNEL_EPSILON` is a computed alias of `kernel_f64`, not an independent value.

**Status (2026-11-13):** Centralization is being migrated incrementally.
This catalog distinguishes all the formerly ad-hoc named numerical constants
listed below. An inline literal is not automatically a tolerance: mathematical
constants (0, 1, 2, PI, powers for format conversion) belong in algorithms,
while comparisons that admit/reject geometry must reference a named policy
field. Do not conflate default *numeric resolution* with a user's maximum
*permitted geometric change*.

## Units and meanings

- **Millimetres (mm):** only in the STEP adapter after source units are converted.
  The Monstertruck core is unitless; its 'spatial' thresholds have units of
  the coordinate system of the supplied geometry. If supplied STEP geometry is
  in mm then its spatial values are mm.
- **UV / curve parameter units:** not necessarily mm. A universal 1e-6
  parameter step has no defensible direct geometric meaning; relate it to
  the local Jacobian/condition number, curvature, parameter span and requested
  spatial residual when building the next operation-local policy.
- **Dimensionless (direction cosine / normalized residual):** a dot-product
  difference of 1e-8 corresponds to about 0.008 degrees for unit vectors near
  parallel; an angular difference in radians is a different quantity.
- **Squared distance:** must compare with squared lengths; never pass a
  linear threshold directly to a squared predicate.
- **Maximum modeling error budget:** a limit on permitted source geometry
  changes, not an instruction to make every point equal within that distance.

## Kernel inventory: current active defaults

| Field | Current | Meaning / units |
|---|---:|---|
| `kernel_f64` | `1.0e-6` | Scalar epsilon for the established `f64` kernel. |
| `boolean_default_numeric_length` | `1.0e-6` | Default numerical length supplied to classic Boolean operations. Independent of the generic kernel zero threshold and of the maximum geometric-change budget selected by a caller. |
| `newton_curve_parameter_step` | `1.0e-6` | Absolute convergence step in native curve parameter coordinates; not a spatial error budget. |
| `newton_surface_parameter_step` | `1.0e-6` | Absolute convergence step in native surface UV coordinates; not a millimetre distance. |
| `newton_surface_curve_parameter_step` | `1.0e-6` | Absolute convergence step in mixed surface-UV plus curve parameter coordinates. |
| `newton_intersection_parameter_step` | `1.0e-6` | Absolute convergence step in four surface UV coordinates during SSI. |
| `kernel_f32` | `1.0e-4` | Scalar epsilon for explicitly `f32` geometry. |
| `structural` | `1.0e-12` | Threshold for structural/algebraic zero tests that must be much tighter than ordinary geometric coincidence. |
| `angular` | `1.0e-6` | Dimensionless angular/direction threshold. |
| `trim_parameter_relative` | `1.0e-6` | Relative parameter-space endpoint tolerance used while preserving trims. |
| `trim_parameter_roundoff_ulps` | `8.0` | Native trim cut boundary ulps, independent of the world-space kernel distance. |
| `healing_projection_parameter_relative` | `1.0e-6` | Native UV span-relative projection overshoot accepted during healed trim reconstruction. |
| `healing_projection_parameter_roundoff_ulps` | `8.0` | Absolute-parameter-magnitude floating-point ulps during healed trim UV normalization. |
| `split_parameter_relative` | `1.0e-8` | Relative parameter-space tolerance used while splitting curves. |
| `split_parameter_roundoff_ulps` | `8.0` | Float-representation roundoff measured in each native edge/trim parameter domain; no unrelated fixed floor. |
| `parameter_projection` | `1.0e-9` | Tight parameter-space projection threshold. |
| `parameter_range_relative` | `1.0e-9` | Relative padding for finite parameter ranges derived from external data. |
| `periodic_relative` | `1.0e-6` | Relative tolerance for periodic-parameter branch comparisons. |
| `transform_similarity` | `1.0e-9` | Relative threshold for ordinary rigid/similarity-transform recognition. |
| `same_domain_similarity` | `1.0e-10` | Stricter similarity threshold used by proof-bearing same-domain logic. |
| `analytic_extrusion_identity` | `1.0e-8` | Homogeneous extrusion identity threshold. |
| `analytic_plane_identity` | `1.0e-6` | Planar analytic-carrier identity threshold. |
| `analytic_spherical_revolution` | `1.0e-7` | Spherical/revolution analytic-recognition threshold. |
| `transversality_relative` | `1.0e-8` | Relative threshold for rejecting nearly parallel surface normals. |
| `mesh_seed_floor` | `1.0e-3` | Coarsest numerical mesh seed used by algorithms that need a seed mesh. |
| `half_space_relative_margin` | `1.0e-3` | Relative expansion used to construct finite half-space clipping bounds. |
| `mesh_open_boundary_direction_relative` | `1.0e-6` | Native UV span-relative directional discrimination when closing an open meshing contour. |
| `mesh_projection_parameter_relative` | `1.0e-6` | Native parameter interval relative rounding allowance for numerical boundary projection. |
| `mesh_boundary_wrap_relative` | `1.0e-6` | Native UV span-relative seam offset used for bounded/periodic parameter normalization. |
| `mesh_boundary_roundoff_ulps` | `8.0` | Absolute native-parameter ulps for representable UV boundary rounding. |
| `mesh_open_pair_connector_score_relative` | `1.0e-6` | Dimensionless margin comparing two sums of squared UV endpoint-connector distances. |
| `mesh_periodic_singularity_tangent_length` | `1.0e-6` | World-length equivalent of per-radian periodic tangent magnitude when detecting mesh poles and validating transverse probe changes. |
| `mesh_periodic_singularity_tangent_relative` | `1.0e-6` | Dimensionless relative reduction of periodic tangent magnitude required at a recovered singularity candidate. |
| `mesh_periodic_section_spatial_residual` | `1.0e-6` | World-space residual when testing whether a sampled periodic loop lies on one fixed geometric section. |
| `mesh_periodic_span_relative` | `1.0e-6` | Relative native-period threshold when qualifying full-turn spans and compacting wrapped loops. |
| `mesh_periodic_seam_endpoint_relative` | `1.0e-6` | Relative native-period endpoint exclusion when constructing meshing seam loops. |
| `mesh_periodic_period_match_relative` | `1.0e-6` | Dimensionless relative compatibility threshold between paired meshing periods. |
| `mesh_uv_closure` | `1.0e-3` | UV closure threshold for already-sampled tessellation boundaries. |
| `mesh_singularity_distance` | `1.0e-3` | Spatial distance used by tessellation to rescue samples at analytic poles. |
| `healing_boundary_uv_relative` | `1.0e-6` | Per-native-axis span-relative closed-loop endpoint discrepancy in healing's 2D chart. |
| `healing_boundary_uv_roundoff_ulps` | `8.0` | Representable-coordinate ulps for UV loop endpoint closure, independent of 3D distance. |
| `normal_probe_step` | `1.0e-4` | Parameter-space probe step used to recover a normal at a singular NURBS point. |
| `derivative_residual` | `1.0e-4` | Residual accepted by derivative consistency assertions. |
| `surface_span_factor` | `8.0` | Multiplier from requested tessellation tolerance to surface-span segment length. |
| `text_closure` | `1.0e-7` | Default closure threshold for text-outline reconstruction. |
| `finite_difference_step` | `1.0e-6` | Central finite-difference step for fallback derivative evaluation. |
| `healing_sampling_floor_factor` | `10.0` | Minimum sampling tolerance as a multiple of the kernel epsilon in healing. |
| `recovery_geometry_factor` | `100.0` | Multiplier used for best-effort STEP missing-edge recovery. |
| `recovery_continuation_jump_factor` | `8.0` | Maximum continuation jump relative to the guide step in missing-edge recovery. |
| `same_domain_direction_floor_factor` | `16.0` | Kernel-epsilon floor used to choose stable radial evidence in same-domain proofs. |
| `cross_face_stitch_factor` | `8.0` | Kernel-precision envelope for cross-face stitching. |
| `boolean_weld_projection_factor` | `2.0` | Projection padding around Boolean weld candidates. |
| `boolean_surface_projection_factor` | `4.0` | Residual envelope for Boolean surface projection evidence. |
| `boolean_curve_endpoint_factor` | `8.0` | Endpoint envelope for recovered Boolean iso-curves. |
| `boolean_polyline_grid_factor` | `2.0` | Kernel-epsilon cell width multiplier for Boolean polyline hashing. |
| `boolean_half_space_kernel_factor` | `100.0` | Kernel-epsilon floor multiplier for finite half-space bounds. |
| `boolean_fragment_validation_factor` | `4.0` | Kernel-epsilon envelope for stitched Boolean intersection fragments. |
| `boolean_boundary_search_factor` | `8.0` | Kernel-epsilon envelope for Boolean boundary parameter search. |
| `mesh_point_hash_factor` | `50.0` | Kernel-epsilon cell width multiplier used by VTK point hashing. |
| `analytic_curve_axis_unit` | `1.0e-8` | Unitless squared-normalization residual for analytic-curve support axes. |
| `boolean_leader_parameter_slack` | `1.0e-8` | Dimensionless leader segment extrapolation allowance when reconstructing SSI. |
| `boolean_leader_progress_slack` | `1.0e-8` | Dimensionless progress monotonicity guard along a leader polyline. |
| `boolean_fragment_overlap_parameter_slack` | `1.0e-5` | Dimensionless overlap-parameter comparison between SSI fragments. |
| `boolean_planar_relative_span` | `1.0e-10` | Relative span padding for planar polygon overlap. |
| `boolean_iso_parameter_relative` | `1.0e-6` | Native spline/iso parameter domain-relative interval slack; not a spatial length. |
| `boolean_iso_parameter_roundoff_ulps` | `8.0` | Machine ulps at absolute parameter magnitude added to relative slack. |
| `boolean_iso_validation_segments` | `32` | Subinterval count for independent iso candidate surface-incidence sampling. |
| `boolean_inverse_search_trials` | `100` | Newton iteration budget in Boolean conversion inverse surface parameter search. |
| `boolean_uv_loop_area_relative` | `1.0e-6` | Signed polygon area relative to its UV bounding-box area; dimensionless, not a world-space length. |
| `boolean_uv_boundary_contact_relative` | `1.0e-6` | UV contact distance relative to combined UV bounding extent; dimensionless, not a world-space length. |
| `boolean_ruled_structural_length` | `1.0e-10` | Model-coordinate structural residual for matched rational ruled extrusions; not a topology-weld allowance. |
| `boolean_ruled_direction` | `1.0e-10` | Dimensionless maximum sine between ruled sweep axes; also constrained by world-space deviation over complete sweep. |
| `boolean_ruled_weight_relative` | `1.0e-10` | Relative homogeneous profile control weights and within-row sweep-weight agreement. |
| `boolean_ruled_knot_relative` | `1.0e-10` | Normalized profile/sweep knot-position agreement under positive affine reparameterization. |
| `boolean_cylinder_structural_length` | `2.0e-9` | Model-coordinate length for structural rational-cylinder recognition; not Boolean welding tolerance. |
| `boolean_cylinder_axis_direction` | `2.0e-9` | Dimensionless cross-axis alignment for certified circular-cylinder carriers. |
| `boolean_cylinder_homogeneous_weight` | `2.0e-9` | Absolute homogeneous weight agreement for linear rational-cylinder sweep rows. |
| `boolean_cylinder_knot_relative` | `2.0e-9` | Native circular-cylinder profile knot span-relative interior split/dedup slack. |
| `boolean_cylinder_polynomial_relative` | `2.0e-9` | Relative residual of each homogeneous quadratic circle polynomial coefficient against its contributing-term magnitude. |
| `boolean_cylinder_reference_axis_component_limit` | `0.8` | Unitless switch for selecting a transverse cylinder-chart reference axis with a well-conditioned cross product. |
| `boolean_planar_end_fraction` | `1.0e-6` | Fractional endpoint exclusion used by planar overlap sampling. |
| `boolean_planar_primary_probe_ratios` | `[1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]` | Dimensionless ordered probe ratios for planar overlap disambiguation. |
| `boolean_planar_secondary_probe_ratios` | `[1.0e-7, 1.0e-6]` | Dimensionless ordered refinement ratios for planar overlap disambiguation. |
| `boolean_planar_normal_cosine_slack` | `1.0e-10` | Dimensionless allowed deviation of aligned planar carrier normals. |
| `boolean_planar_tangent_cosine_slack` | `1.0e-8` | Dimensionless allowed deviation of planar trim tangent alignment. |
| `meshing_debug_uv_edge_tolerance` | `1.0e-9` | Native UV-coordinate tolerance used for tessellation projection diagnostics. |
| `step_torus_inverse_sine_roundoff_ulps` | `128.0` | Unitless ulps for STEP torus polar inverse asin argument clamping; no world-size multiplier. |
| `step_trim_knot_parameter_relative` | `1.0e-6` | Maximum native STEP knot-domain-relative overshoot admitted for numeric endpoint clamping. |
| `step_trim_knot_parameter_roundoff_ulps` | `8.0` | Absolute-parameter-magnitude roundoff added to the STEP native-domain clamp. |
| `step_uv_clamp_probe` | `1.0e-3` | Native UV-coordinate threshold for STEP trim seam probe. |
| `step_trim_residual_mm` | `1.0e-8` | Millimetre residual for STEP trim curve reconstruction. |
| `step_output_distance_accuracy` | `1.0e-6` | STEP file uncertainty metadata (current model units). |
| `step_output_scientific_notation_threshold` | `1.0e-2` | Non-geometric STEP number formatting threshold, dimensionless. |
| `boolean_fragment_stitch_scale` | `4.0` | Scale multiplying Boolean numeric tolerance when joining near-fragments. |
| `boolean_coplanar_scale` | `2.0` | Scale multiplying Boolean numerical tolerance in coplanar matching. |
| `boolean_classifier_offset_factor` | `8.0` | Classifier witness offset relative to explicit Boolean length tolerance. |
| `boolean_area_floor_factor` | `16.0` | Area comparison's absolute squared-length coefficient. |
| `boolean_area_span_factor` | `8.0` | Area comparison's length-times-span coefficient. |
| `boolean_fragment_forward_alignment_cosine` | `0.8` | Minimum cosine between adjacent SSI fragment tangents before joining. |
| `boolean_cylinder_normal_alignment_cosine` | `0.9` | Minimum absolute cosine between recovered cylinder normal and radial direction. |
| `boolean_interior_witness_insets` | `[ 0.1, 0.03, 0.01, 0.003, 1.0e-3, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12, ]` | Ordered dimensionless inset fractions for finding a strict interior witness in thin planar regions. |
| `boolean_coincidence_fractions` | `[0.0, 0.25, 0.5, 0.75, 1.0]` | Bounded sample fractions for bidirectional whole-curve coincidence evidence. |
| `boolean_containment_fractions` | `[0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]` | Bounded sample fractions for trimmed subedge containment and same-endpoint correspondence evidence. |
| `stl_quantization_offset_factor` | `0.25` | Offset multiplier for STL coordinate quantization. |
| `stl_quantization_step_factor` | `0.5` | Step multiplier for STL coordinate quantization. |

## STEP reconstruction inventory: current active defaults

| Field | Current | Meaning / units |
|---|---:|---|
| `numerical_geometry_mm` | `1.0e-7` | Numerical consistency of a reconstructed STEP support, not permission to heal a feature. Unit: mm. |
| `feature_heal_mm` | `0.01` | Maximum source-feature reconciliation distance, subject to separate topology checks. Unit: mm. |
| `cad_project_standard_mm` | `0.01` | Default permitted STEP CAD modeling-error budget. Unit: mm. |
| `cad_edge_selector_mm` | `1.0e-7` | Matching user-selected edge endpoints against their literal coordinates. Unit: mm. |
| `cad_axis_angle_rad` | `1.0e-10` | Detecting revolution full turns and zero angular sweeps. Unit: radians. |
| `cad_axis_locus_mm` | `1.0e-9` | Detecting vertices and edges along a revolution axis. Unit: mm. |
| `cad_profile_join_mm` | `1.0e-8` | Joining successive CAD profile segments. Unit: mm. |
| `cad_profile_angle_rad` | `1.0e-12` | Rejecting vanishing profile arc sweep. Unit: radians. |
| `cad_normal_unit` | `1.0e-12` | Rejecting zero or nonunit CAD edge selector normals. Unit: dimensionless. |
| `surface_geometry_mm` | `1.0e-5` | Surface recovery/recognition residual against input coordinates. Unit: mm. |
| `surface_weight_relative` | `1.0e-12` | Homogeneous NURBS weight comparison for recovery. Unit: dimensionless. |
| `surface_parameter` | `1.0e-12` | Knot/parameter endpoint comparison for surface recognition. Unit: native surface parameter. |
| `solid_geometry_mm` | `1.0e-7` | Geometric residual in extrusion, loft, revolution, and sweep reconstruction. Unit: mm. |
| `solid_direction` | `1.0e-10` | Direction cosines/orthogonality for extrusion and loft recognition. Unit: dimensionless. |
| `sweep_direction` | `1.0e-9` | Direction cosines in swept-profile recognition. Unit: dimensionless. |
| `revolution_direction` | `1.0e-9` | Axis/orientation agreement for revolution recognition. Unit: dimensionless. |
| `loft_direction` | `1.0e-10` | Loft support normal agreement. Unit: dimensionless. |
| `dominant_extrusion_direction` | `1.0e-8` | Coarse axis grouping to identify dominant extrusions. Unit: dimensionless. |
| `edge_treatment_direction` | `1.0e-8` | Direction agreement for edge-treatment recognition. Unit: dimensionless. |
| `inverse_edge_treatment_direction` | `1.0e-10` | Direction agreement in inverse fillet/chamfer identification. Unit: dimensionless. |
| `reduced_brep_direction` | `1.0e-10` | Direction agreement in reduced boundary representation. Unit: dimensionless. |
| `sheet_geometry_mm` | `1.0e-5` | Thickness and overlap matching for sheet-metal recovery. Unit: mm. |
| `sheet_direction` | `1.0e-9` | Direction agreement for formed sheet faces and bends. Unit: dimensionless. |
| `formed_sheet_direction` | `1.0e-10` | Tighter formed-sheet normal and axis agreement. Unit: dimensionless. |
| `solid_extrusion_angle_rad` | `1.0e-10` | Degenerate-angle threshold in extrusion inference. Unit: radians. |
| `periodic_geometry_mm` | `1.0e-5` | Spatial residual for periodic chains and periodic bodies. Unit: mm. |
| `periodic_resize_mm` | `1.0e-7` | Coordinate matching while resizing periodic features. Unit: mm. |
| `pattern_geometry_mm` | `1.0e-7` | Pattern center and lattice matching during feature recognition. Unit: mm. |
| `pattern_orientation` | `1.0e-10` | Pattern basis orientation and angular consistency. Unit: dimensionless. |
| `partition_point_quant_mm` | `1.0e-8` | Quantization of feature partition vertices. Unit: mm. |
| `partition_normal_quant` | `1.0e-10` | Quantization of partition-normal directions. Unit: dimensionless. |
| `geometric_intern_position_mm` | `1.0e-5` | Position matching when interning equivalent supports. Unit: mm. |
| `geometric_intern_scalar_mm` | `1.0e-5` | Scalar offset matching when interning equivalent supports. Unit: mm. |
| `geometric_intern_direction` | `1.0e-10` | Normal and direction matching when interning supports. Unit: dimensionless. |
| `boundary_point_quant_mm` | `1.0e-7` | Position quantization for independent boundary proof. Unit: mm. |
| `boundary_scalar_quant_mm` | `1.0e-7` | Scalar offsets quantization in boundary proof. Unit: mm. |
| `boundary_direction_quant` | `1.0e-10` | Axis/normal quantization in independent boundary proof. Unit: dimensionless. |
| `boundary_locus_scale` | `4.0` | Factor converting boundary scalar quantum to admitted support locus residual. Unit: multiplier. |
| `boundary_angular_scale` | `4.0` | Factor converting boundary direction quantum to proof normal residual. Unit: multiplier. |
| `line_collinear` | `1.0e-12` | Line segment angular collinearity comparison. Unit: dimensionless. |
| `line_monotonic` | `1.0e-12` | Ordering slack in line recovery. Unit: dimensionless. |
| `line_direction_quant` | `1.0e-12` | Direction discretization for line recovery. Unit: dimensionless. |
| `line_endpoint_direction` | `2.0e-13` | Endpoint direction consistency in line recovery. Unit: dimensionless. |
| `line_degenerate_chord_squared` | `1.0e-24` | Squared chord length for degenerate recovered lines. Unit: mm². |
| `source_brep_direction` | `1.0e-15` | Very small vector accepted for source BREP direction normalization. Unit: dimensionless. |
| `curve_replica_geometry_mm` | `1.0e-5` | Geometric agreement of repeated source curves. Unit: mm. |
| `curve_replica_transform_mm` | `1.0e-12` | Pure translation residual for exact curve replication. Unit: mm. |
| `surface_replica_transform_mm` | `1.0e-12` | Pure translation residual for exact surface replication. Unit: mm. |
| `planar_feature_side` | `1.0e-8` | Side and orientation checks in planar feature recognition. Unit: dimensionless. |
| `planar_feature_fit_mm` | `1.0e-7` | Planar feature candidate residual. Unit: mm. |
| `revolution_max_source_uncertainty_mm` | `1.0e-5` | Upper source-data uncertainty for revolution reconstruction. Unit: mm. |
| `axis_group_dot` | `1.0e-10` | Axis grouping criterion in compound CAD recovery. Unit: dimensionless. |
| `station_position_mm` | `1.0e-7` | Station clustering in cap/loft CAD reconstruction. Unit: mm. |
| `source_cap_frame` | `1.0e-10` | Frame orthonormality for reconstructed source cap. Unit: dimensionless. |
| `source_sweep_frame` | `1.0e-8` | Orthonormality for reconstructed source sweep profiles. Unit: dimensionless. |
| `source_sweep_strict_residual_mm` | `1.0e-6` | Strict original sweep support fitting error. Unit: mm. |
| `source_sweep_loose_residual_mm` | `1.0e-5` | Larger source sweep support fitting error. Unit: mm. |
| `pattern_basis_cost_delta` | `1.0e-12` | Tie breaker between alternative pattern basis fits. Unit: dimensionless cost. |
| `pattern_quantization_floor` | `1.0e-9` | Minimum position quantization for feature pattern candidates. Unit: mm. |
| `pattern_determinant_squared` | `1.0e-18` | Reject nearly degenerate two-dimensional pattern bases. Unit: mm⁴. |
| `pattern_lattice_normal_floor` | `1.0e-12` | Denominator floor for pattern coordinates. Unit: dimensionless. |
| `pattern_fraction_quant_floor` | `1.0e-10` | Lattice fractional-coordinate quantization floor. Unit: dimensionless. |
| `pattern_origin_floor_mm` | `1.0e-9` | Minimum origin matching scale for pattern transforms. Unit: mm. |
| `pattern_comparison_residual` | `1.0e-12` | Tolerance for pattern filling-ratio ranking ties. Unit: dimensionless. |
| `pattern_normalization_min` | `1.0e-15` | Minimum vector norm accepted for pattern direction normalization. Unit: dimensionless. |
| `source_straightness_length_mm` | `1.0e-12` | Near-zero geometric chord threshold for source support extraction. Unit: mm. |
| `source_axis_unit_cosine` | `1.0e-9` | Acceptable axis norm deviation in source feature frames. Unit: dimensionless. |
| `source_axis_direction_cosine` | `1.0e-8` | Axis-normal agreement for source features. Unit: dimensionless. |
| `source_sweep_angle_rad` | `1.0e-8` | Minimal accepted bend angle and supplement for source sweeps. Unit: radians. |
| `source_sweep_tangent_cosine` | `1.0e-7` | Tangent and profile orientation consistency in source sweep recovery. Unit: dimensionless. |
| `source_sweep_parameter_slack` | `1.0e-7` | Progress and joined segment angular residual in sweep recovery. Unit: radians. |
| `source_sweep_curve_parameter` | `1.0e-10` | Parameter-comparison slack during source sweep recovery. Unit: native curve parameter. |
| `source_sweep_knot_residual` | `1.0e-12` | Knot equality and source support weight matching. Unit: native knot/weight parameter. |
| `source_sweep_small_weight` | `1.0e-14` | Reject near-zero B-spline homogeneous support weight. Unit: native homogeneous weight. |
| `source_sheet_angle_rad` | `1.0e-8` | Minimal valid sheet bend angle. Unit: radians. |
| `source_sheet_normal_cosine` | `1.0e-7` | Sheet face normal alignment and near-parallel classification. Unit: dimensionless. |
| `source_sheet_tangent_cosine` | `1.0e-8` | Sheet seam and tangent alignment. Unit: dimensionless. |
| `source_sheet_bend_angle_rad` | `1.0e-7` | Sheet bend target-angle fitting. Unit: radians. |
| `source_identity_zero_norm_squared` | `1.0e-24` | Squared size below which source direction/area evidence is degenerate. Unit: mm². |
| `source_identity_axis_unit` | `1.0e-9` | Source coordinate frame orthonormality threshold. Unit: dimensionless. |
| `source_identity_transform` | `1.0e-12` | Zero-motion matching in source rigid transformations. Unit: mm. |
| `source_identity_quantization_inv_mm` | `1.0e9` | Reciprocal position quantization step for source identity keys. Unit: mm^-1. |
| `source_identity_normal_quantization` | `1.0e5` | Inverse normal quantization step for source identity keys. Unit: dimensionless. |
| `source_feature_incidence_mm` | `1.0e-7` | Point-on-original-source support recovery check. Unit: mm. |
| `source_analytic_fit_mm` | `1.0e-7` | Analytic surface candidate fitting in boundary proof. Unit: mm. |
| `source_boundary_tangent_fraction` | `1.0e-5` | Finite parameter sampling fraction for boundary tangent proof. Unit: dimensionless. |
| `source_revolution_angle_floor_rad` | `1.0e-12` | Minimum valid source meridian-angle comparison. Unit: radians. |
| `source_revolution_param_floor` | `1.0e-12` | Lower bound for parameter-space arc geometry slop. Unit: native curve parameter. |
| `source_cad_frame_orthogonality` | `1.0e-10` | Orthonormality of source CAD frames and transformed vectors. Unit: dimensionless. |
| `source_rounding_arc_residual_rad` | `1.0e-6` | Allowed angular mismatch in reverse fillet arc recovery. Unit: radians. |
| `source_rounding_arc_max_slack_rad` | `1.0e-7` | Allowed angular excess beyond a full rounding sweep. Unit: radians. |
| `source_rounding_parameter_slack` | `1.0e-9` | Bound on source rounding parameter interval comparisons. Unit: native curve parameter. |
| `source_rounding_quarter_turn_slack_rad` | `1.0e-10` | Allowed slack when decomposing corner arcs into quarter turns. Unit: radians. |
| `source_line_residual_mm` | `1.0e-7` | Straight source curve fit against expected 3D locus. Unit: mm. |
| `source_cap_point_match_mm` | `1.0e-8` | Point matching while reconstructing source cap boundary. Unit: mm. |
| `source_cap_direction_cosine` | `1.0e-8` | Direction agreement in source cap and circle sweep frames. Unit: dimensionless. |
| `source_cap_rounding_residual_mm` | `1.0e-12` | Source cap and sphere common-support matching. Unit: mm. |
| `source_intern_radius_scale_relative` | `1.0e-10` | Relative radius offset bound in source support interning. Unit: dimensionless. |
| `source_periodic_direction_quant` | `1.0e-10` | Direction quantization for periodic source feature chains. Unit: dimensionless. |
| `pattern_inverse_norm_squared` | `1.0e-24` | Reject near-degenerate pattern-recovery support vectors. Unit: mm². |
| `source_sweep_absolute_residual_mm` | `1.0e-5` | Maximum acceptable fitted sweep support residual. Unit: mm. |
| `source_sweep_relative_compare` | `1.0e-9` | Scale-dependent sweep support comparison. Unit: dimensionless. |
| `source_sheet_relative_compare` | `1.0e-9` | Scale-dependent sheet geometry comparison. Unit: dimensionless. |
| `source_periodic_weld_quantization_inv` | `1.0e12` | Inverse coordinate quantization for periodic graph weld keys. Unit: model-coordinate^-1. |
| `source_periodic_area_min` | `1.0e-12` | Minimum nondegenerate reference area for periodic graph editing. Unit: mm². |
| `source_sweep_tolerance_slack_mm` | `1.0e-15` | Tiny allowed source-fitted sweep tolerance budget overrun. Unit: mm. |
| `source_slot_boundary_mm` | `1.0e-7` | Position check for cap-slot boundary recovery. Unit: mm. |
| `source_cad_arc_angle_slack_rad` | `1.0e-9` | Arc sweep endpoint slop in CAD feature reconstruction. Unit: radians. |

## Candidate ranges for our source-STEP reconstruction target

These are **engineering starting hypotheses, not changes already made**.
They are deliberately *operation-specific* and must be checked against
source STEP examples and geometric difference, not just face count.

| Family / job | Candidate target | Why / required check |
|---|---|---|
| Structural algebraic identity (exact knot structure, homogeneous proportionality) | f64 relative error ~1e-12--1e-10, scaled by coefficient magnitude and conditioning | Reject incompatible supports without trying to 'heal' them; separate source representation quantization from geometric equivalence |
| STEP source fitting / recognition of intended features | 1e-6--1e-4 mm for clean CAD, up to ~0.01 mm only when source discrepancy is measured and topology can still be disambiguated | Prefer residuals measured from original STEP entities and a scale-dependent bound; a fixed 0.01 mm across models can collapse smaller details |
| CAD reconstruction max geometry deviation | User-selected bound, often ~0.005--0.02 mm for supplier CAD recovery; for high-fidelity editing use ~0.001 mm | Output checks must measure Hausdorff-like boundary distance, continuity, orientation and signed volume, not merely face count |
| Surface/curve sampling for *seeding* intersections | Curvature/feature-dependent chord deviation ~0.005--0.05 mm, with a minimum-feature guard; subdivide adaptively | Tolerance-only benchmark reached a plateau at 0.01 mm sampling; coarse samples must never be used alone to prove incidence |
| Exact surface-intersection residual | Conditioned spatial residual, typically aim 1e-7--1e-5 mm on mm-normalized source; escalate precision only where required | Need local Jacobian, floating-point scale, tangent near-singular handling; strict 4e-6 mm global cutoff is not a proof of branch identity |
| Newton parameter convergence | No single meaningful universal value. Stop based on spatial residual + scaled parameter update, e.g. normalized relative 1e-8--1e-6 and conditioning-aware acceptance | Current 1e-6 is in *parameter coordinates*, can mean vastly different mm on different surfaces; log residual and step separately |
| Curve/face/edge coincidence identity | Prefer shared provenance and exact support/parameter relations; use 1--10 µm spatial closeness as a **candidate filter only**, tighten if features are closer | Never merge solely on coordinate proximity; direction, branch topology, orientation and connected trim graph must agree |
| Seam stitching and trim endpoint matching | Start ~1e-6--1e-4 mm for clean STEP; otherwise use measured source seam uncertainty and feature-scale separation | Verify no crossing, loop closure and two-sided support incidence; avoid hardcoded fixed multiples of kernel epsilon |
| Analytic plane/extrusion/revolution recognition | Algebraic certificate where possible, otherwise compare normalized residual at ~1e-10--1e-8 relative to shape scale | Recognition enables huge solver savings; don't reject perfect analytic supports due to world-origin offset |
| Direction/normal matching | Dot-product residual ~1e-10--1e-7 depending conditioning and source noise | Unit-vector normalization and orientation matters; never confuse angular radians with cosine residual |
| STEP export reported accuracy | Encode actual conservative shape accuracy (not a fixed 1e-6 claiming 1nm precision) | Metadata must track measured input/output residual; no constant can certify the final result |
| Failure thresholds (conditioning, max iteration/cell counts) | Typed resource/conditioning policy, not a geometric-error tolerance | Limit explosion deterministically and report why evaluation was refused; never silently accept an inaccurately sampled result |

## Call graph and known non-propagation

`step-redox::GeometryPolicy::Standard { tolerance_mm }` represents a **maximum
permitted model change**, and `step-redox::MonstertruckKernel::evaluate_boolean`
converts it to `monstertruck-solid::GeometryPolicy::Standard`.

`monstertruck-solid/src/transversal/integrate/mod.rs::classic_boolean_tolerance`
currently validates that policy and then **deliberately returns
`KERNEL_EPSILON`**, which is independent of the user-supplied Standard
budget. This is a concrete gap in the future accuracy-budget propagation
architecture, not an unnoticed constant: there is a dedicated regression test
for this policy choice. Do not silently reinterpret a maximum modeling-error
budget as unconditional permission to weld unrelated edges.

**Correction after direct test validation:** the old
`BooleanTolerance::new(requested)` conversion `requested.abs().max(KERNEL_EPSILON)`
was **not needed** to satisfy the ten existing planar tests at 1e-8 and 1e-7.
Those tests failed only because a previous assert incorrectly rejected
legitimate sub-epsilon requests. Removing the clamp outright and preserving
the positive finite caller-supplied length yielded **136/136 passing
`monstertruck-solid` tests**, including all ten. The math/topology tests are
valid for their exact geometry fixtures; there is no evidence that imposing
a 1e-6 lower bound is required for the algorithms.

The constructor now requires a finite *positive* numeric tolerance, retains
its exact value, and separately computes the **coarser seed mesh tolerance**
using the named `mesh_seed_floor` (which is intentionally only a sampling
choice). Coincident-edge reuse similarly accepts any valid positive tolerance
instead of refusing all values smaller than `KERNEL_EPSILON`. This is
distinct from the `GeometryPolicy::Standard` maximum modeling-change budget.

The public low-level classic numerical validator **also** previously rejected
positive requests below `KERNEL_EPSILON`. Following the successful no-clamp
tests, this validator now accepts any finite positive request as well.
There is no numerical-resolution floor on geometric coincidence or topology
admission. Separate finite positive validation remains mandatory, while the
coarser mesh seeding tolerance is explicitly marked as sampling-only.

The Boolean-specific parameters for classifier offsets, area comparison,
fragment continuity, cylinder normal alignment, and witness/curve sampling
schedules are now named fields in `DEFAULT_NUMERICAL_TOLERANCES`. The two
small-fraction schedules are **sample evidence**, not proof of analytic
curve identity. They should be replaced with genuine geometric certificates
where possible; changing sampling density alone cannot establish topology.

The public classic Boolean entry point now selects a **separately named**
`boolean_default_numeric_length` instead of coupling its numerical setting
to generic `KERNEL_EPSILON`. This keeps the existing 1e-6 behavior until
tuning and preserves the distinct semantics of `GeometryPolicy::Standard`
(upper bound on modeling error). The low-level classic numerical interface
accepts smaller explicit finite positive values.

The classic Boolean path then uses `BooleanTolerance::new` to derive
`length`, `classifier_offset`, and `mesh_seed` from its numerical
argument, and many surface/curve Newton and intersection helper paths still
use other kernel policy values. Grouped defaults are no replacement for
passing semantically correct tolerances through those APIs.

## Migration contract (no hidden defaults)

1. Each operation chooses **named** numerical convergence, sampling, proof,
   and permitted-geometry-change budgets; no caller should construct those
   from unexplained in-code constants.
2. The chosen policy is propagated explicitly. Missing, NaN, negative, and
   incompatible units are **errors**, not silently corrected with
   `.abs()`, `.max(KERNEL_EPSILON)`, or a fallback default.
3. Kernel-global defaults are presets and initial API values, **not** a
   substitute for an operation-specific tolerance context.
4. Numerical methods must not claim success solely because successive
   solver parameters are close: test residuals and conditioning.
5. Record error distribution, feature disambiguation, topology validity,
   solver counts, and wall time when the values are tuned later.

## Evidence already collected

- [Boolean tolerance sweep PR #103](https://github.com/codeandsolder/monstertruck/pull/103)
  changed curve chord sampling and *requested* modeling budget, but **not**
  the fixed classic Boolean numerical tolerance.
- The sampling-only sweep at 1nm/1um/10um/100um reduced NURBS searches from
  951,774 to 854,934 to 440,053 to 439,358, with stage4 single-solid
  47-face topology still present; no full shape-difference measurement.
- Changing the *requested* 10nm modeling budget to 1um and 10um did nothing
  to Newton counts **because** `classic_boolean_tolerance` ignores that
  budget by design; those runs do *not* prove effective numeric tolerance
  sensitivity is absent.
- For deeper provenance and stage4/5 topology details see the project
  [Notion technical record](https://app.notion.com/p/410f74be902081cb8cb5c176b3be5f40).

## Not the same thing: literal numeric constants

Do **not** move PI, zero/one identity tests, polynomial coefficients,
unit-conversion exponents, scalar precision constants (`f64::EPSILON`),
or formula-derived subdivisions into this registry merely to eliminate
digits from source. They are definitions of mathematics or the file format,
not tunable acceptance thresholds.

## SSI precision propagation (2026-11-13)

The marching surface/surface intersection lane now receives the complete
operation's `BooleanTolerance` from `create_loops_stores` through
`intersection_curves`, `build_intersection_curve`, curved fragment correction,
and overlapping-fragment validation. Previously, these methods independently
multiplied global `KERNEL_EPSILON`, which meant sub-epsilon operation requests
were effectively ignored by critical proof predicates.

- **Exact carrier-support residual:**
  `boolean_surface_projection_factor * BooleanTolerance.length`.
- **Validated fragment residual:**
  `boolean_fragment_validation_factor * BooleanTolerance.length`.
- **Candidate displacement / broad-phase seed:** `BooleanTolerance.mesh_seed`
  (explicitly a coarser sampling budget; never a coincidence proof).
- **Redundant-strand removal, closure, and local continuation checks:**
  the operation's `BooleanTolerance.length`, rather than global epsilon.
- **Linear fragment validation:** has separate support-surface and fragment
  thresholds even though their present default multipliers both happen to be 4.

The numerical values remain unchanged for the existing default `1e-6`
Boolean operation. A regression checks that a point `5e-7` away from an
exact planar carrier is *rejected* at a `1e-8` request but can be accepted
at `1e-6`; passing a policy value now changes actual SSI proof behavior.
`cargo +1.99.1 test -p monstertruck-solid --lib` reports 142 passing tests,
with strict Clippy and formatting passing.

**Outstanding:** LoopStore's trim-boundary lookup, endpoint canonicalization,
`curve_surface_projection`, and cross-face stitching still use independent
global tolerances. Those functions need the same operation-local policy,
including distinct UV/parameter search and world-coordinate acceptance,
before numerical tuning can be considered end-to-end. Do not conflate
`GeometryPolicy::Standard` (maximum error permitted) with the actual
numerical tolerance spent by topology classification.

## Classic Boolean boundary and cross-face tolerance propagation (2026-11-13)

The LoopsStore now receives the operation's typed BooleanTolerance through
BranchContext and BoundaryEndpoint. The previous CROSS_FACE_STITCH_TOLERANCE
constant, evaluated once from global epsilon, has been removed.

- Polygon boundary insertion uses BooleanTolerance.boundary_residual()
  rather than a fixed global 8 * KERNEL_EPSILON.
- Exact boundary supported-vertex candidates keep the mesh-seed broad-phase
  search radius separate from the exact carrier-support residual
  (BooleanTolerance.surface_support_residual()) and candidate ambiguity
  (BooleanTolerance.boundary_residual()).
- Curve/surface tangent endpoint correction requires independently projected
  exact carriers to agree within BooleanTolerance.length; predictor
  displacement may extend to BooleanTolerance.boundary_residual().
- Cross-face endpoint canonicalization uses
  BooleanTolerance.cross_face_stitch_residual() for unique candidate search;
  collapse of short SSI edges and deduplication of candidate vertices uses
  BooleanTolerance.length, never mesh seed.
- Coplanar-imprint endpoint insertion now passes its explicit operation-local
  length through to the polygon insertion path.

Existing default 1e-6 operations retain 4e-6 carrier-support, 8e-6 boundary,
and 8e-6 cross-face stitch thresholds. New geometric regressions exercise a
real 5e-7 separation and verify that a 1e-8 numerical request rejects both
cross-face stitching and an off-carrier tangent projection.

Still to fix: the general search_boundary_parameter inverse search calls
SearchParameter (which has its own numerical precision) and may accept its
result as exact without independently checking the returned spatial residual.
Its coarse nearest-search allowance must stay distinct from exact admission
and from parameter-space endpoint classification.

## Exact boundary inverse search must certify spatial error

The legacy Line::search_parameter (and other curve implementations) may use
kernel-global Point::near to report a nominally exact inverse. LoopsStore's
search_boundary_parameter used to accept any such answer without checking the
actual distance between curve.evaluate(parameter) and the query point.
Consequently, an operation requesting 1e-8 precision could accept a 5e-7
lateral discrepancy as exact solely because the curve's own inverse search
uses 1e-6. The outer Boolean policy was not reaching this admission decision.

The exact-incidence path now checks Euclidean evaluation residual against the
operation's own BooleanTolerance.length, independent of whether the inverse
solver said the parameter was exact. A separate nearest-search pass retains
its explicit mesh discovery radius and returns the curve-projected coordinate;
that broad-phase allowance is not labeled proof of exact incidence. The
hinted exact-search path performs the same independent residual check before
considering the separately budgeted nearest candidate. A regression queries
a line at 5e-7 lateral offset: no exact hit at a 1e-8 request, yet its nearest
projection remains discoverable when explicitly supplying a coarse search
radius. Source values are unchanged for default operations.

Follow-up: topological endpoint coalescing needs independent second-carrier
validation even for a candidate discovered by a wide mesh seed. The
operation-local accuracy budget must be used to certify the actual locus,
not merely to find the candidate.

## Opposite-carrier incidence and operation-local curve equivalence (2026-11-13)

A boundary vertex discovered within the coarse mesh search radius is **not**
necessarily on the other Boolean carrier. `add_geom_vertex` now checks the
actual canonical exact endpoint through the other surface's nearest-parameter
search followed by an independent forward-evaluated spatial residual. It must
fit `BooleanTolerance.surface_support_residual()` **before** `set_point` or
`change_vertex` mutates topology. The same admission function is used in
`supported_boundary_vertex`, and that check now precedes candidate grouping;
an unsupported neighboring vertex cannot inherit the group's valid status.

Tests cover (a) rejecting a unique nearby boundary that lies 5e-7 off the
opposite planar carrier with a requested 1e-8 numerical length, for both
hinted and unhinted searches, and verifying the original edge and vertex are
unchanged; and (b) accepting an approximate seed when its *exact* canonical
endpoint really does lie on both carriers.

Coplanar imprinting separately validates the forward residual of every
`SearchParameter` result against the operation length before counting its edge
as incident. This prevents a global `Point::near` result from splitting a
parallel, genuinely distinct edge 5e-7 away at requested 1e-8 tolerance.

Curve-equivalence evidence is also now parameterized by the explicit typed
`BooleanTolerance`. The bounded sample/projection tests for whole-curve
coincidence, trimmed-subedge containment and endpoint-aligned matching use
spatial distances against **operation.length**, not their own implicit global
`near()` thresholds. The operation policy is threaded through LoopStore
boundary candidate ambiguity, repeated-edge reuse, coplanar imprinting, and
edge insertion. The final region-weld helper explicitly selects the current
Boolean numerical default, rather than receiving a hidden tolerance from
these generic validators. It still needs a fully operation-local policy at
its upstream interface in a later pass.

A red-then-green regression independently demonstrated the preexisting bug:
two parallel edges 5e-7 apart were incorrectly treated as **one** nearest
boundary candidate with operation length 1e-8. With explicit local coincidence
tolerance the two candidates remain ambiguous, and a separate curve-sampling
regression confirms the same precision distinction. The changed helpers are
validation evidence, not mathematical proof that arbitrary oscillating curves
cannot diverge between samples; structural/provenance-based equality remains
the longer-term target.

## Endpoint classification must use geometry, not parameter epsilons

The old `ParameterKind::try_new(parameter, (front, back))` used the generic
`f64::near` epsilon to decide `Front`, `Back` or `Inner`. That is dimensionally
wrong: curve parameters are native (often normalized, angles, knots or arbitrary
scales) and their delta has no fixed relationship to model-coordinate length.
For a 100000-unit line, parameter delta 1e-8 is 0.001 model units; for a
1e-6-unit line, parameter delta 0.01 is only 1e-8 model units.

The new classification explicitly receives the bounded curve and the active
`BooleanTolerance`. It requires finite ordered parameter bounds and an in-range
parameter; exact endpoint parameters map to their explicit ends. For an
interior parameter, evaluate the point and compare its **3D distances to both
trim endpoints** against the operation's spatial numerical length. If within
budget of precisely one endpoint, canonicalize to it; if within neither,
retain the interior parameter; if within both, return ambiguity rather than
inventing precedence. An extrapolated inverse root may classify as the corresponding endpoint only when its forward-evaluated model-space location lies within the operation length tolerance of that endpoint and not the other endpoint. There is no raw-parameter epsilon.
The exact-trim search invoked by coplanar imprint additionally independently
verifies `curve.evaluate(parameter).distance(query)` against the same length.

This is an endpoint *classification*, not a proof that the source curve and
the sampled SSI leader coincide between samples. For self-intersections and
closed periodic curves, topology/provenance may still be necessary to resolve
ambiguous repeated positions; a proximity-only rule must not pick one side.
Validation should include long and short line parameterizations, ambiguous
sub-feature edges, out-of-domain inverse roots, sub-epsilon coincident trim
proof, full Boolean regression, and actual source-derived STEP reconstruction.

## Spatial index and operation-local Boolean weld (2026-11-13)

The previous classic Boolean weld used global 1e-6 kernel epsilon as
both its three-dimensional hash cell size and the actual coincidence
test. It divided world coordinates by this epsilon, then converted the
result to signed 64-bit cell indices. Values outside the integer range
silently saturated, and nonfinite values were assigned an arbitrary cell
instead of signaling malformed input. This was not a sound spatial index.

The new weld context is constructed with the same explicit
BooleanTolerance as the Boolean operation. Vertex lookup uses a
deterministically ordered BTreeMap of original floating-point X
coordinates for candidate search. It then checks full three-dimensional
Euclidean distance against exactly the requested length; the search
structure itself does not define geometric equality and does not
quantize or saturate any coordinate. Nonfinite points reject the weld
explicitly. This supports arbitrarily large finite model coordinates.

Both source planar noding and final Boolean region welding now receive
the same operation tolerance, threaded through their actual callers.
Exact/nearest projection candidates, curve identity checks, split point
acceptance, and consistency of split analytic curve endpoints likewise
use operation numerical length instead of unrelated global epsilon.
The explicitly named boolean_weld_projection_factor, currently 2,
pads only the spatial candidate bounding box and is not an identity
or model-edit allowance. The old parameter-span cutoff at f64::EPSILON
is replaced by finite ascending-domain validation, because raw
curve-parameter deltas have no fixed spatial meaning.

Tests check preserving separate vertices 5e-7 apart at a requested
1e-8 while allowing merge at an explicit 1e-6, near-coincidence at
5e-9, curve coincidence at distinct requests, finite X coordinates
whose prior hash indices would saturate i64, and explicit rejection of
nonfinite coordinates. Previous regression coverage of exact line
splits and planar union is preserved.

Performance: BTreeMap candidate lookup is logarithmic in distinct X
coordinates plus the number of candidates in the requested spatial
interval. There is no extra fixed numerical grid tolerance to tune.
Any large-source performance effect must be measured before claiming
an improvement or regression.

## Operation-local planar-chart admission (2026-11-13)

The Boolean planar overlay has two separate tasks: decide whether a projected
chart crossing represents a new **exact** edge node, and choose a convenient
chart for a known analytic carrier. Neither may silently widen a caller's
numeric coincidence tolerance to the global `KERNEL_EPSILON`.

- `needs_exact_noding` now receives the operation's numerical length and
  ignores an apparent UV T-junction only when its original world-space
  endpoints actually agree within that length. Its scale-relative UV
  crossing tests retain their separately named chart-space threshold; the
  UV ratio is not a spatial modeling error allowance.
- `CylinderChart::new` qualifies the first world-space vertex of a certified
  circular carrier using the same requested length. The rational-cylinder
  certificate itself remains independent of that qualification.
- `classification::sample_points_on_face` rejects degenerate tessellation
  triangles against the operation-local `BooleanTolerance.area_floor` in
  squared world-coordinate units, not the square of the fixed kernel
  epsilon. A rejected mesh still takes the preexisting trimmed-region
  interior-sampling path.

Focused regressions distinguish 500 nm-separated endpoints and cylinder
radial discrepancies at requested 10 nm and 1 um numerical lengths. These
tests check numerical identity only, not permission to modify the model by
that amount. Remaining `divide_face` UV-area and UV boundary-contact
thresholds need separate Jacobian/conditioning analysis: their parameter
units cannot be replaced by a raw world-space tolerance.


Verification on the Rust 1.99.1 Asus worktree: 158/158 `monstertruck-solid`
unit tests pass and all-target strict Clippy passes. The original STEP
`mixed_complete_constructive_ir_replays_without_intermediate_step_import`
replay completed stage 4 (void) in 95.434 s, producing the asserted valid
47-face solid. Stage 5 still fails at the previously known three-face oriented
open shell. Solver totals are unchanged from previous verified branch tip:
948,784 3D Newton calls, 1,004,854 4D Newton calls, and 951,220 NURBS
presearch queries. This change does not claim a measured speedup; the
95.434 s wall time is confounded by host load.

## UV degeneracy and chart-basis thresholds (2026-11-13)

Two previously unrelated uses of `KERNEL_EPSILON` in
`classic/divide_face.rs` are now explicitly named in the shared numerical
policy. The unsigned/signed UV polygon-area degeneracy gate is controlled
by `boolean_uv_loop_area_relative` multiplied by the loop's UV bounding-box
area. The disjoint-hole boundary-contact gate is controlled separately by
`boolean_uv_boundary_contact_relative` multiplied by the combined UV
extent, and squared before comparing with squared point-to-segment distance.
Both quantities are **parameterization-relative**, not lengths in mm;
neither is derived from the Boolean operation's world-space length.
Likewise the cylinder chart's transverse basis switch uses the dimensionless
`boolean_cylinder_reference_axis_component_limit`, separate from its
operation-local geometric radial admission.

The initial values `1e-6`, `1e-6` and `0.8` respectively preserve the
existing exact numerical thresholds; this pass makes their independent
roles inspectable, not experimentally validated or automatically suitable
for every possible UV mapping. The existing scale-invariance,
touching-boundary, and nondegenerate UV loop regression tests remain
unchanged. Future parameterization-conditioning experiments must consider
nonuniform scaling/shearing, not confuse source-space tolerance with UV
contact. The full 158-test solid suite and strict all-target solid Clippy
remain green after this source-equivalent extraction.
## Exact Boolean iso-section conversion (2026-11-13)

The classic Boolean's conversion of a marching SSI curve to an analytically
certified iso-parameter NURBS previously used the kernel-global `1e-6` as
**both** a 3D spatial length and an absolute UV/knot parameter interval.
Even with an explicit Boolean numerical length of `1e-8`, a recovered
candidate could pass checks permitting `4e-6` surface residual and `8e-6`
endpoint residual. No requested modeling-change budget authorized this.

Now the `BooleanTolerance` received by the Boolean operation travels
through `convert_regions`, `convert_alt_shell`, `altshell_to_shell`,
`exact_iso_leader` and all exact iso candidate, surface, and endpoint
proof helpers. The weld's `certified_alt_curve_box` also uses the same
operation length when recognizing an exact iso carrier.

**Spatial evidence:** `surface_support_residual()` (4× caller's Boolean
numeric length) checks candidate/world carrier incidence. Endpoint
agreement uses `boundary_residual()` (8× caller's numeric length).
Both are world-coordinate numerical admission, **not permission to
alter original source geometry**.

**Parameter evidence:** Exact NURBS trimming and the U/V iso-coordinate
classification use `iso_parameter_slack((lo,hi))` -- the *relative native
spline parameter interval* (current policy 1e-6 of span) plus eight
machine ulps at the absolute parameter magnitude. Finite parameters
marginally beyond a domain are clamped before exact spline cutting;
nonfinite values, degenerate intervals and actual excursions are
rejected. This avoids dropping a valid 1e-9-wide knot domain simply
because a global 1e-6 length threshold was confused with its parameter
units. The fixed 32-interval candidate validation and 100 Newton
inverse-parameter iterations now have separately named policies.

Additional local tests prove scale-relative native-knot trimming and
operation-length independence. With the previous experiment-only
curve sampling controls, the real 47-face stage4 CAD subtraction passed
in 3.166s; stage5 still has its independent unknown-face classifier
failure. The ordinary production-default precision and independently
measured source-to-result geometry remain separate validation gates.

## Structural rational-cylinder certification threshold separation (2026-11-13)

The circular-cylinder structural certificate previously obtained one
`same_domain_similarity * 20.0` scalar and used it simultaneously for
**coordinate length**, **unit-vector alignment**, **homogeneous weights**,
**native profile knot locations**, and **squared-length polynomial
coefficients**. A cylinder-chart basis also used a raw `0.8` axis
component. Those are not interchangeable tolerances.

The structural recognition thresholds now live in separate named numerical
policy entries (table above). The cylinder-chart basis switch shares the
existing named `boolean_cylinder_reference_axis_component_limit`, since
it is the same angular-conditioning decision. The profile knot interval
slack is scaled to its native parameter range and includes named ulps
for representability at large offsets; it is not a world-space length.
The rational-circle polynomial identity is tested coefficient-by-coefficient
with a dimensionless residual relative to the maximum total absolute
coefficient contribution in the polynomial (squared length times squared
homogeneous weight). One consistent scale accounts for cancellation
and floating-point error in coefficients that vanish analytically. This removes the arbitrary `1.0 + radius2`
scale and avoids giving tiny circles a much larger certificate
tolerance merely because someone assumed 1-unit circles.

The default recognition values equal the previous implicit 2e-9
where meaningful. This is structural algebraic certificate precision,
not the requested Boolean numerical length and not a source-model
editing allowance. All existing near-cylinder refusal regressions
stay enabled, plus a native-knot-reparameterization regression. As a
separate check, the source 47-face constructive stage4 replay must
remain valid, without assuming any stage5 fix.

## Rational common ruled-extrusion identity (2026-11-13)

`classic/loops_store/ruled.rs::common_ruled_extrusion_strip` used one
`same_domain_similarity` scalar for unrelated measurements: native U/V
knot distances, world-coordinate sweep and offset lengths, direction
cross products (including a cross product of **two unnormalized sweep
vectors**, whose units are model-length squared), dimensionless
homogeneous weight ratios, and overlap length.

The four separate policy entries above specify their respective
physical dimensions. Profile and sweep-domain knot vectors are compared
as `(knot - first)/(last - first)` in their own **positive native
parameter spans**, rather than testing absolute knot values: equal
control nets and affine-reparameterized knots describe the same
B-spline support. A malformed/zero-range/nonfinite knot vector
cannot certify identity. Linear swept V knots must still have clamped
end multiplicities, judged relative to the V domain.

Normalized sweep-axis cross products enforce the dimensionless angular
limit, but a very small angular mismatch over a very long sweep can
become geometrically large. A separate perpendicular displacement
bound of `boolean_ruled_structural_length` is therefore checked over
the entire sweep length, along with origin offsets, row vectors,
profile control positions, and strictly positive axial overlap.
Rational weights are positive, finite and compared as dimensionless
ratios to their respective first control weights; a nonfinite
Cartesian projected control point is rejected.

New focused tests verify same-support equality despite affine
rescaling/translation of both U and V native knot domains; reject
non-affine changes to internal knots; and reject a tiny angle whose
long-sweep lateral deviation exceeds the structural spatial bound.
These are algebraic carrier-recognition tolerances, never the
user's Boolean topological coincidence length or permitted source
modeling change. Experimental sampling floors and numeric overrides
remain restricted to the separate benchmark worktree.


## Real STEP knot-domain clamp vs placeholder-axis handling (2026-11-13)

`monstertruck-io/src/step/load/step_geometry/geom_impls/parameter_trims.rs`
previously snapped any STEP UV parameter sitting within an absolute
`1e-6` of the limits of a *genuine knot-bounded axis*, regardless of
the knot span. That wrongly equated the global spatial kernel epsilon
with a native parameter distance. The already documented "placeholder
axis is not a bounded domain" distinction is retained unchanged.

`clamp_near_range` now admits only a finite positive native knot
interval and computes a separately named relative limit proportional
to that actual interval, plus roundoff measured in machine ulps of
the absolute parameter magnitude. All other values remain unchanged,
including excursions too far beyond the real domain and invalid ranges.
For the common [0,1] parameter domain, the prior 1e-6 behavior is
retained; for a 1e-9 knot interval, an overshoot of 1e-11 no longer
gets silently snapped to the boundary. This is a **numerical UV
roundoff rule**, not a maximum geometric displacement or source-feature
healing tolerance. No source STEP model geometry is intentionally
changed by this update.

Existing fixture tests distinguish genuine knot boundaries from
unbounded analytic placeholder axes. Additional direct regressions
check tiny finite knot domains, legitimate representable endpoint
drift, unit-domain compatibility, invalid ranges and nonfinite
values. The historical measurement table above describes the earlier
fixed-kernel-epsilon implementation and must not be interpreted as
post-change metrics.


## UV contour closure and seam normalization (2026-11-13)

The tessellator's open contour closure had three unrelated
`KERNEL_EPSILON` usages in UV-space: choosing whether an open boundary
travels along the U or V axis; nudging a candidate parameter before
periodic range normalization; and breaking ties between alternative
pairing directions using a *sum of squared UV distances*. The latter
compared a model-coordinate length epsilon with a squared UV quantity.

The numeric policy now independently defines UV-domain-relative
direction and seam-wrap slacks, absolute-parameter-magnitude ulps,
and a dimensionless tie-break fraction for the squared-UV connector
scores. For finite positive parameter intervals, direction and
wrapping scale with the native U or V span, instead of silently
treating one UV parameter unit as one model-coordinate millimetre.
Invalid/degenerate parameter intervals cannot authorize a seam wrap
or corner-closure axis decision.

The paired-contour direction choice now compares sums of squared
UV lengths with a dimensionless relative score margin, so a uniform
rescaling of UV coordinates does not flip the chosen connectivity
merely because a fixed 1e-6 was added in the wrong dimension.

Focused tests check 1e-18 to 1e18 squared-UV scale invariance and
small versus unit native-parameter spans. The full meshing library
tests and strict all-target clippy are the verification gate. This
changes only *sampling and triangulation/contour normalization*:
it does not widen exact Boolean topology identity or source CAD
modeling-error permissions. Other absolute thresholds in the
periodic/singularity and mesh/STEP domains remain to be analyzed.


## Meshing native-period comparisons (2026-11-13)

The second periodic-meshing pass removed global `KERNEL_EPSILON`
from full-period loop classification, periodic jump compaction,
seam-window selection, first/last point exclusion on manufactured
seams, and period equality when pairing loops. The original
`1e-6` allowance was an *absolute world-coordinate* length
silently added to native UV angles or periodic parameters.

These decisions now use named relative values in
`DEFAULT_NUMERICAL_TOLERANCES`: span qualification
(`mesh_periodic_span_relative`), manufactured seam endpoint
exclusion (`mesh_periodic_seam_endpoint_relative`), and
paired-period agreement (`mesh_periodic_period_match_relative`).
The earlier `mesh_boundary_wrap_relative` governs seam-window
placement separately. The shared helper
`periodic_parameter_slack(period, native_coordinates, relative)`
returns `period * relative + ulps * EPSILON * max(period, |coordinates|)`,
rejecting nonfinite/invalid periodic inputs rather than comparing
invalid values as valid periods. As a result, a model with a 1e-9
native periodic parameter range cannot accept a 1e-6 *parameter*
drift merely because 1e-6 is the kernel's world-distance epsilon.

Three focused regressions verify period scaling from 1e-9 to 1e9,
seam-window selection at a tiny period, representational ulps at
a large coordinate offset, and refusal of invalid inputs.
Validation: `monstertruck-meshing --lib` 43 passed, 1 ignored;
strict all-target Clippy PASS. The experiment-only accelerated
source CAD stage4 still reaches the expected valid 47-face solid
in 3.067 seconds; stage5 retains the unrelated unknown-face
classifier failure. No default sampling floor or requested Boolean
numerical length was changed.

Explicitly out of scope in this commit: the distinct
`singular_periodic_other_parameter` derivative/radius
conditioning, the other-axis criterion previously derived from
`periodic_relative`, and STEP same-domain parameterization.
Those have separate dimensions and need separate tests.


## STEP same-domain inverse torus parameterization (2026-11-13)

The STEP torus analytic inverse constructed
`sine = relative.z / torus.small_radius()`. This is a
dimensionless ratio whose admissible inverse-trigonometric interval
is exactly [-1,1], with a *dimensionless numerical roundoff* envelope.
The previous code multiplied a fixed `128 * EPSILON` by a
world-coordinate scale (large torus radius, small radius, or point
magnitude), arbitrarily widening the sine domain for large models.
The outer geometric support proof is distinct from this inverse
trig-domain gate and must not govern its admissibility.

`torus_sine_with_roundoff` now validates finite ratios against
`1 + step_torus_inverse_sine_roundoff_ulps * EPSILON * max(1, |sine|)`
and clamps only representable overrun, independent of the torus's
world size. The small radius must be finite and positive before
division. Subsequent radial support and evaluated world-space
carrier checks retain their existing spatial predicates.

The STEP same-domain periodic unwrap previously refused valid
periods smaller than the *world-coordinate* `KERNEL_EPSILON`;
a period is now considered valid if finite and positive. New
focused tests cover torus radii 1e-6 through 1e9, near-±1
roundoff admission and genuine out-of-domain rejection, nonfinite
ratios and a native parameter period 1e-9. The full STEP IO lib
suite passed 84 tests (9 ignored), and strict all-target IO Clippy
also passed. Other same-domain matching/world-space acceptance
tolerances remain independent of this fix.


## Periodic closed-section evidence in model space (2026-11-13)

The periodic ring detector in
`monstertruck-meshing/src/tessellation/triangulation/boundary/periodic.rs`
previously declared a loop to be at a constant nonperiodic coordinate
when `other_max - other_min <= periodic_relative (=1e-6)`.
That applied an unrelated fixed UV difference to a potentially
world-length axial parameter or an angular parameter and was **not**
invariant under surface reparameterization or model scale.

`sampled_constant_periodic_section` now holds the periodic coordinate
at the first sampled value, evaluates the source surface at each
sample's nonperiodic coordinate, and tests the resulting **world-space
distance** to a fixed-section reference under the separately named
`mesh_periodic_section_spatial_residual = 1e-6`. This tests
only the *sampled contour* (appropriate for a meshing decision), not
an algebraic certificate that the entire carrier is a constant
section. It is not reused for topology welding or exact geometric
admission.

Tests use physically equivalent revolutionary cylinder profiles with
0.01, 1 and 1000 model-unit axial sweeps: a 0.25e-6 model-unit
cross-section drift is accepted, whereas a 4e-6 model-unit drift is
rejected for each scale regardless of the completely different
native UV parameter offsets. The preexisting cone/ring meshing
regressions are preserved and the full meshing library and strict
Clippy pass. The real source 47-face stage4 replay is separately
required before production cherry-pick to identify any unexpected
meshing/topology effects.

The still-open singularity inverse
`singular_periodic_other_parameter` contains a distinct model
derivative/parameter slope vs spatial precision decision and is
not modified by this patch.


## Periodic singularity recovery — separate native UV, first derivatives and model length (2026-11-13)

The previous `singular_periodic_other_parameter` algorithm
computed `|dS/dperiodic_parameter|`, a quantity with units
**model length / periodic UV**, called it `radius` and compared it
directly to world-distance `KERNEL_EPSILON`. Its finite difference
`(|sample_derivative|-|origin_derivative|) / delta_other` has
units **model length / (periodic UV × other UV)** and was tested
against that same world distance. The candidate's derivative residual
used a third implicit allowance inherited from generic
`periodic_relative`. Its probe step dropped any native range narrower
than `KERNEL_EPSILON`, falling back to a completely arbitrary `1 UV`
regardless of model size or nonperiodic parameter scaling.

The recovery path now explicitly obtains the finite, positive
native periodic-axis period and converts the derivative magnitude
to a **characteristic physical periodic tangent length per radian**
with `|dS/dparameter| × native_period / TAU`. For circular
revolutions it equals the actual circle radius; for a generic periodic
surface it is a local tangent-scale proxy, **not** an algebraic radius
or singularity certificate. Two distinct entries in the central
numerical policy control world-scale tangent smallness
(`mesh_periodic_singularity_tangent_length`) and dimensionless
relative reduction at a candidate
(`mesh_periodic_singularity_tangent_relative`).

Transverse parameter probes use either a finite positive native
domain span or an inferred **UV step corresponding to one periodic
tangent-length of physical travel**, from the ratio to the
nonperiodic surface derivative. Both probe distances remain in
native *other-axis* UV units. No absolute 1-UV fallback is needed.
The finite difference rejects insufficient *physical tangent-length
change*, rather than comparing a derivative slope to a model
distance. Candidate parameterization is computed from that ratio
with the proper UV dimension; nonfinite or invalid values refuse.

Finally, a derivative vanishing at **one** periodic coordinate
isn't enough evidence that a ring collapses. Both near-zero
direct and extrapolated pole candidates now require sampled
3D agreement of the source surface at quarter-turn, half-turn,
three-quarter-turn and full-turn positions, measured against the
existing independently named
`mesh_periodic_section_spatial_residual` world-length
criterion. This is a **sampling check appropriate to meshing**, not a
continuous proof of the entire surface. Empty curves and invalid
periodic-axis indices return `None` rather than indexing blindly.
No source CAD geometry or Boolean admissibility criterion changes.

Focused regressions exercise geometrically identical infinite
conical support with native axial parameter scales
`1e-9, 1e-4, 1 and 1e3`, recovering the identical 3D apex
for each, as well as a nondegenerate constant-radius cylinder
that cannot have an apex and invalid/empty inputs. The existing
full-period cone-cap/ring tests remain enabled. The source CAD
constructive replay's stage4 was again valid at 47 faces with the
**separately experimental** faster sampling settings; stage5 is
still the pre-existing unknown-face classifier failure.


## Mesh-boundary projection parameter normalization and pole recovery (2026-11-13)

In `monstertruck-meshing/src/tessellation/triangulation/boundary/projection.rs`,
native U/V parameter boundary snapping used a fixed *model-distance*
`KERNEL_EPSILON` irrespective of the span of the interval being
normalized. The period branch also judged a native UV span
`span.so_small()` using that spatial epsilon, collapsing every
valid periodic parameter domain narrower than 1e-6 down to its
lower endpoint. A purported 3D pole check evaluated the surface
one unit above the previous V parameter, regardless of its actual
period, and only supported bounded-U/periodic-V charts.

The projector now shares `periodic.rs::native_boundary_slack`,
with a separate `mesh_projection_parameter_relative` policy
and the existing absolute-UV-magnitude floating-point ulp policy.
Finite positive native ranges receive relative/ulp endpoint rounding;
a small but positive UV period is never confused with a degenerate
zero-width domain. Invalid/nonpositive/nonfinite periodic values
are rejected before computing normalized coordinates.

`PolyBoundaryPiece::recover_periodic_pole` now handles either
UV orientation, requiring an actual positive finite period on
the other axis and a finite bounded-axis endpoint whose 3D
surface position lies within the dedicated
`mesh_singularity_distance` model-space mesh rescue threshold.
It checks that **the whole periodic orbit collapses** using
`sampled_periodic_collapse` at 1/4, 1/2, 3/4 and one complete
native turn, rather than testing only a guessed +1 UV displacement.
This is a meshing rescue, not an exact support identity proof.
The periodic collapse helper explicitly rejects invalid periods
and nonfinite seed parameters before evaluation.

New tests cover native UV intervals of 1e-9, the original
unit-interval rounding behavior, invalid periodic values,
both bounded-U/periodic-V and bounded-V/periodic-U cone poles,
a non-polar cylinder, far-away points, and invalid sampled periods.
Existing meshing tests and the accelerated source CAD replay
remain required prior to promotion.

The reported [0,1] plane ranges in the core
`Plane::parameter_range` are synthetic, not intrinsic
geometric boundaries (see that implementation's existing
documentation). Simply scaling the numeric clamp does not by
itself distinguish such artificial ranges: that requires
explicit parameter-domain provenance in the surface interface,
and is a separate semantic issue to resolve rather than
heuristically assuming every reported bound is physical.

## Intrinsic surface parameter bounds versus synthetic display/meshing domains (2026-11-13)

A distinct semantic problem remained after the numerical native-UV
clamp fix: core `Plane::parameter_range` deliberately reports a
fake [0,1] by [0,1] square, despite the mathematical plane
having no finite parameter limits. The artificial square is still
needed by older bounded-surface consumers, so simply replacing
`parameter_range` with unbounded values would change unrelated
API/meshing behavior and expose existing assumptions. However,
the numerical projector previously took `try_range_tuple` as
authoritative and would incorrectly snap valid planar UV coordinates
just past u=1 or below v=0.

`monstertruck_traits::ParametricSurface` now explicitly reports
axis-wise provenance with
`parameter_bounds_are_intrinsic() -> (bool, bool)`.
The default is (true,true) for ordinary bounded mathematical
supports, and the core Plane overrides to (false,false).
`&S`, `Box<S>`, `Processor<S,T>` (including swapped UV
orientation), and `OffsetSurface<S,N>` forward the underlying
bound truth instead of silently turning an artificial source
range into a trusted range. Both the enum and tuple-struct
`ParametricSurface` and `ParametricSurface3D` procedural
derives forward the method for nested STEP surface enums and
other wrappers. The method is a semantic contract about each
surface's *reported bounds*, not a retroactive change to the
values returned by `parameter_range`.

The meshing boundary projector uses reported finite bounds
**only when their corresponding intrinsic bit is true**.
Its parameter normalization still uses the explicit native
interval-relative plus ulp threshold for genuine knot-defined
surfaces. For a Plane, no clamping against the fabricated
[0,1] square occurs; sample coordinates u>1 and v<0 are retained.
Tests exercise bare Plane, transformed Processor<Plane>,
a genuine bilinear B-spline knot domain, and the two nested
STEP enum layers via the same trait method, verifying that
the metadata survives procedural derive forwarding.

This corrects projection semantics without deleting the
synthetic square still needed by bounded-surface renderers.
It does not yet fix **other** consumers which incorrectly
assume that every reported range is intrinsic, including the
healing/analytic SSI sites called out in
`monstertruck-geometry/src/specifieds/plane.rs`.
Those should migrate to the same explicit provenance contract,
with independent fixture tests for each operation rather than
forcing all consumers into a global all-or-nothing range change.


## Healing trim reconstruction: intrinsic domain provenance and native parameter precision (2026-11-13)

`monstertruck-healing/src/lib.rs` previously duplicated a
`normalize_uv` closure in the sampled-trim and linear-trim
regeneration routes. Both used `try_range_tuple()` as
unconditionally intrinsic, turning an unbounded Plane's
synthetic [0,1] × [0,1] meshing rectangle into a **hard-clamped
geometric domain**. Both also called
`normalize_parameter_axis`, which judged native periodic
intervals via spatial `span.so_small()` and collapsed every
positive parameter span smaller than the 1e-6 model-space
kernel epsilon.

Now both recovery routes share
`normalize_surface_parameter(surface, uv, previous)`.
It checks `surface.parameter_bounds_are_intrinsic()`
separately for U and V, retains unbounded/synthetic-plane UV
values exactly, and uses separate finite native parameter
range-relative plus absolute-magnitude ulp thresholds from
`healing_projection_parameter_relative` and
`healing_projection_parameter_roundoff_ulps` on real
bounded axes only. Unlike the old hard `value.clamp`, values
substantially outside an actual finite support **fail
normalization** rather than silently teleporting a trim onto a
different 3D position. True zero-width ranges are treated as
degenerate, not any interval less than world-space epsilon.
Invalid periods/nonfinite hints are rejected; periodic
unwrap under a hint uses the nearest integer multiple of
the actual native period, without the previous arbitrary
five-turn search limit.

The separate `oriented_trim_segment` also had four
mixed-unit sites: checking a native trim interval with
`trim_len.so_small()`, then defining its trim split
tolerances as `max(trim_len*trim_parameter_relative,
KERNEL_EPSILON)`, including each subsequent cut and
inner segment. These world-distance floors swallowed
entire tiny but valid native trim intervals. The
`trim_parameter_slack((t0,t1))` helper now computes
`native_positive_span * trim_parameter_relative
 + trim_parameter_roundoff_ulps * EPSILON *
 max(|t0|,|t1|)`, checking positive finite spans.
Before cutting, source fractions must be finite and
inside the original trim range to within that same
native parameter slack; out-of-range cut values are
not silently clamped. The ***geometric*** distance
verification between the preserved 3D trim and segment
still uses the caller's world-space `tol`, separately.
Topology identity and source modeling tolerances were
not enlarged.

Regressions cover:
- A Plane with actual UV (3.9,-2.1) and an inverted
  Processor<Plane>, both surviving normalization despite
  synthetic [0,1] reported ranges;
- A true two-direction bounded B-spline surface,
  snapping a 5e-7 native overshoot but refusing a
  0.02-unit arbitrary projection displacement;
- Tiny native periodic span 1e-9; unwrapping more than
  two turns; valid degenerate/nondegenerate intervals,
  invalid finite/nonfinite parameters;
- Scale-invariant native trim cut/slack behavior from
  1e-10 to 1e8; and an actual 3D line backed by
  a knot-bounded B-spline curve with **1e-9 native
  trim span** but **1.0 world-unit physical length**:
  cutting out the middle quarter-to-three-quarter
  segment must preserve the legitimate 3D geometry.

Validation: `monstertruck-healing --lib` 27 passed, 1
ignored; strict `monstertruck-healing --all-targets`
Clippy passed on Rust 1.99.1 with sccache. Existing
experimental source constructive replay still reaches
the expected valid 47-face stage4, and stage5 retains
the known shell1 unknown-face classification issue.

**Remaining:** `same_surface` in the healing path
still samples synthetic parameter rectangle corners
as a heuristic, which is not a support-equality proof.
Other analytic SSI domain/clip logic may live outside
the published monstertruck crates; inspect its actual
code and data before changing another consumer.


## Pass-through edge and trim split units; preserve full wire topology (2026-11-13)

`monstertruck-healing/src/split_pass_through_edges.rs`
previously applied `split_parameter_absolute=1e-10` as a
fixed floor to both native edge and trim parameter spaces. It
used `.so_small()` (a model-space world-distance epsilon)
on `edge_len` and `trim_scale`. A valid 1e-12 native
edge-parameter span was assigned an arbitrary half-way split
`alpha=0.5` even when the actual fractional position was
0.2; positive 1e-9 trim spans were discarded. Moreover
`split_trim_curve` deduplicated the resulting *trim-native*
coordinates against the **edge-native** `param_eps`, even
though the two parameter scales can be arbitrarily different.

The new `split_parameter_slack((first,last))` takes the
actual positive finite native interval and computes
`span * split_parameter_relative
 + split_parameter_roundoff_ulps * EPSILON *
 max(|first|, |last|)`. The central policy removes
`split_parameter_absolute`, replacing it with the ulp-based
roundoff factor. `split_parameters` computes edge-native
slack; `split_trim_curve` computes **distinct** edge-domain
and trim-domain slacks for determining valid interior splits
and deduplicating mapped trims. It never invents a 0.5
fraction when the native edge interval is small.

The separate *topological cardinality* invariant is now
explicit: if one source 3D edge is split into N new edge
uses, there must be exactly N valid, individually split
UV trim curves before preserving their geometry. The old
`zip(edge_indices, trim_segments)` silently discarded
edge uses if the trim side returned too few segments,
leaving the face wire shorter and potentially invalid.
`split_trim_curve` now returns `Option<Vec<T>>` and
refuses invalid/nonfinite/out-of-range/merged cuts.
`replace_trimmed_edge_uses` retains **all** N edge uses
when the trim mapping fails, setting each
`trim_curve=None` so subsequent trim regeneration may
recover their geometry, rather than quietly removing the
topological edges. A cardinality check precedes any zip.

Regression coverage: genuine 1e-9-native-span B-spline
trim mapped from a [0,1] edge domain splits exactly at
one-quarter/three-quarter; edge-native span 1e-12
maps the actual 0.2 fraction into a [0,1] trim,
not midpoint 0.5; an unrepresentable zero-width trim
mapping preserves three edge indices with missing trims;
and native slacks vary by correct parameter scale for
spans 1e-12 through 1e8, distinguishing coordinate
representation ulps from geometric acceptance.

No world-space modeling tolerance, SSI threshold, or
experimental mesh refinement setting is changed by
this patch. This behavior should be submitted as an
independent upstream correctness fix with a compact
failing-before/passing-after reproduction and a
topology-cardinality assertion; separate upstream
issue tracking is maintained in the project Notion
hardening bug catalog.


## Closed face loop classification: UV parameter evidence **and** spatial evidence (2026-11-13)

`monstertruck-healing/src/split_closed_faces/mod.rs::boundary_type`
previously tested `first.distance2(last) <
healing_boundary_closure_squared` with one fixed
`1e-3` applied to **UV distance squared**, despite the
native chart dimensions being arbitrary. This
threshold could falsely mark a physically open wire
as closed when the UV gap was small in absolute
numeric coordinates (e.g. 0.02 in a unit chart),
or reject a physically closed wire represented in a
chart with 1e9 native units despite equal real 3D
geometry. It also treated zero-area contours as
`Negative`, potentially causing a degenerate loop
to enter hole-handling logic.

The two existing callers, `non_closed_wires_in_param_divisor`
and `assort_boundary`, already have both the actual
`ParametricSurface3D` and the operation's world-space
`tol`; they now pass both to `boundary_type`.
The corrected algorithm checks closure **separately on U
and V**, each against `native_axis_span *
healing_boundary_uv_relative + representation-ulps *
max(abs(min), abs(max))` (the two new central policy
fields). It also evaluates the actual 3D surface at
both endpoints and demands finite squared distance
`<= tol*tol`. This prevents a tiny UV ratio from
masking a physically large gap on an ill-scaled
surface, without conflating parameter and spatial
units. Invalid/nonfinite inputs, insufficient
vertices and zero/nonfinite oriented area classify
as `NotClosed` rather than a genuine
Positive/Negative loop. A genuine negative-oriented
finite-area wire remains Negative.

Regressions use bilinear B-spline surfaces with
**identical physical unit-square geometry** in
native U/V charts parametrized independently as
1e-9/1e9, 1e9/1e-9, and 1/1. Physically nearly
closed contours with 0.5e-6 fractional UV endpoint
drift are consistently recognized as Positive or
Negative according to orientation, while 1e-4
fractional open gaps are refused. A second
B-spline has 1e9 world-length scale along one axis
but ordinary [0,1] UV; even a 0.5e-6 fractional UV
gap is refused because its **world-space gap**
is huge. Degenerate and invalid sample loops,
negative/NaN/infinite tolerances, and the canonical
STEP healing fixtures are tested separately.

The original source face topology and native
periodicity are not inferred from coincident 3D
points alone. A contour must close in the
**unwrapped UV chart** *as well as* spatially,
so a periodic full-turn seam is not silently
mistaken for an oriented planar UV loop.
