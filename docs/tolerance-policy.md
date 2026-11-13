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
| `split_parameter_relative` | `1.0e-8` | Relative parameter-space tolerance used while splitting curves. |
| `split_parameter_absolute` | `1.0e-10` | Absolute floor paired with `split_parameter_relative`. |
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
| `mesh_uv_closure` | `1.0e-3` | UV closure threshold for already-sampled tessellation boundaries. |
| `mesh_singularity_distance` | `1.0e-3` | Spatial distance used by tessellation to rescue samples at analytic poles. |
| `healing_boundary_closure_squared` | `1.0e-3` | Squared UV-distance threshold inherited by closed-face healing. |
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
| `boolean_planar_end_fraction` | `1.0e-6` | Fractional endpoint exclusion used by planar overlap sampling. |
| `boolean_planar_primary_probe_ratios` | `[1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5]` | Dimensionless ordered probe ratios for planar overlap disambiguation. |
| `boolean_planar_secondary_probe_ratios` | `[1.0e-7, 1.0e-6]` | Dimensionless ordered refinement ratios for planar overlap disambiguation. |
| `boolean_planar_normal_cosine_slack` | `1.0e-10` | Dimensionless allowed deviation of aligned planar carrier normals. |
| `boolean_planar_tangent_cosine_slack` | `1.0e-8` | Dimensionless allowed deviation of planar trim tangent alignment. |
| `meshing_debug_uv_edge_tolerance` | `1.0e-9` | Native UV-coordinate tolerance used for tessellation projection diagnostics. |
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
