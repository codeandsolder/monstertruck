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
| `boolean_planar_primary_probe_ratios` | `[1.0e-8` | Dimensionless ordered probe ratios for planar overlap disambiguation. |
| `boolean_planar_secondary_probe_ratios` | `[1.0e-7` | Dimensionless ordered refinement ratios for planar overlap disambiguation. |
| `boolean_planar_normal_cosine_slack` | `1.0e-10` | Dimensionless allowed deviation of aligned planar carrier normals. |
| `boolean_planar_tangent_cosine_slack` | `1.0e-8` | Dimensionless allowed deviation of planar trim tangent alignment. |
| `meshing_debug_uv_edge_tolerance` | `1.0e-9` | Native UV-coordinate tolerance used for tessellation projection diagnostics. |
| `step_uv_clamp_probe` | `1.0e-3` | Native UV-coordinate threshold for STEP trim seam probe. |
| `step_trim_residual_mm` | `1.0e-8` | Millimetre residual for STEP trim curve reconstruction. |
| `step_output_distance_accuracy` | `1.0e-6` | STEP file uncertainty metadata (current model units). |
| `step_output_scientific_notation_threshold` | `1.0e-2` | Non-geometric STEP number formatting threshold, dimensionless. |
| `boolean_fragment_stitch_scale` | `4.0` | Scale multiplying Boolean numeric tolerance when joining near-fragments. |
| `boolean_coplanar_scale` | `2.0` | Scale multiplying Boolean numerical tolerance in coplanar matching. |
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
