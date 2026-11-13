# Effective Boolean tolerance sweep -- 2026-11-13

**Scope:** source-derived six-stage `4a6356` CAD IR fixture on Asus, Rust 1.99.1, active PR104 parent `e3e682d9`. This is an **experiment branch**. Do not merge runtime `MT_EXPERIMENT_*` environment overrides into production without an explicit accuracy/conditioning policy.

## Independent controls

- `MT_EXPERIMENT_CURVE_SAMPLE_FLOOR_MM` raises *adaptive modeled-curve chord sampling* floor. It is not a topology-equality tolerance.
- `MT_EXPERIMENT_CLASSIC_NUMERIC_LENGTH` changes the **actual** numerical length used by classic Boolean coincidence, endpoint admission and some dependent work. This is independent from `GeometryPolicy::Standard` maximum modeling error; increasing numerical length spends an unproven portion of that budget.
- `MT_EXPERIMENT_BOOLEAN_MESH_FLOOR_MM` raises only the mesh discretization floor in `BooleanTolerance::new`. The original seed floor is `0.001` model units.
- `MT_EXPERIMENT_EXPORT_STAGES_ONLY` is a test-diagnostic flag for exporting CAD stages without evaluating a separately failing upper-round cutter beforehand.

All lengths are in model units (mm for this STEP example). The original test still fails at stage 5; stage 4 has an independently asserted valid single 47-face solid. No stage6 or full source recovery is claimed.

## Measured matrix (Asus, single real source, wall time load-sensitive)

| Curve floor | Numeric length | Mesh floor | Stage2 s | Stage3 s | Stage4 s | NURBS seeds | 3D Newton | 4D Newton | Stage5 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| default | 1e-6 | 1e-3 | 13.850 | 15.934 | 95.434 | 951220 | 948784 | 1004854 | oriented open shell |
| 0.01 | 1e-6 | 1e-3 | 9.186 | 12.899 | 38.460 | 439519 | 437440 | 492060 | oriented open shell |
| 0.01 | 1e-5 | 1e-3 | 8.939 | 11.886 | 25.998 | 245869 | 243758 | 298168 | oriented open shell (different loop) |
| 0.01 | 1e-4 | 1e-3 | 8.898 | 11.825 | 23.646 | 158577 | 156297 | 211923 | interface-construction failure |
| 0.01 | 1e-6 | 0.01 | 2.291 | 3.049 | 23.968 | 375933 | 374868 | 427218 | unknown-face classification failure |
| 0.01 | 1e-5 | 0.01 | 2.328 | 3.264 | 21.379 | 206843 | 205774 | 258244 | unknown-face classification failure |

For each row stage4 passed its existing 47-face single-solid assertions. The last two rows have a different later failure; that matters even though stage4 passed. The large improvement from default to curve floor 0.01 is driven by fewer expensive Newton searches, not merely wall-clock variation. Remaining floor and geometry-budget effects are not yet mapped to a general production policy.

## Independent OpenCascade evidence

With `MIXED_CONSTRUCTIVE_IR_OUTPUT_DIR` and `MT_EXPERIMENT_EXPORT_STAGES_ONLY=1`, exported the stage-4 STEP solids for the unmodified-default setting and the fastest combined 0.01/1e-5/0.01 setting.

- Both: 1 solid, 47 faces, valid OCC BRep.
- Original volume: 7.392843574930893 mm^3; fast candidate volume: 7.392842432812021 mm^3.
- OpenCascade `BRepAlgoAPI_Cut` original-minus-candidate and candidate-minus-original: both report **0.0 mm^3** residual and valid results.
- Signed/integrated volume discrepancy is ~1.142118872e-6 mm^3 (relative ~1.55e-7); therefore results are OCC-Boolean-equivalent at its effective resolution but **not certified bitwise or analytic-trim identical**.
- Stage5 still fails; its error changed, so loosening tolerances is **not** a stage5 repair.

Source-only validation is deliberately not inferred from matching face count, exported STEP BRep validity, or empty OCC differences alone. Full original-oriented trims, continuous deviation bound and stage5/stage6 remain unverified.

## Reproduction

Run (with one variable set at a time for independent axes) from this worktree:

```bash
MT_NEWTON_PROFILE=1 MONSTERTRUCK_PRESEARCH_TRACE=1 \
  MT_EXPERIMENT_CURVE_SAMPLE_FLOOR_MM=0.01 \
  MT_EXPERIMENT_CLASSIC_NUMERIC_LENGTH=1e-5 \
  MT_EXPERIMENT_BOOLEAN_MESH_FLOOR_MM=0.01 \
  cargo +1.99.1 test -p step-redox --features cad-kernel-monstertruck --lib \
    mixed_complete_constructive_ir_replays_without_intermediate_step_import -- --nocapture
```

The original branch `refactor/tolerance-policy-audit-20261113` is unmodified. These knobs and the stage-only export guard exist only on this isolated experimental branch. Run logs are ignored under `experiments/effective-tolerance-sweep-20261113/*.log`; exported comparison STEP files are kept in the dedicated `/home/jan/build/monstertruck-effective-tolerance-sweep-20261113/geometry/` folder on Asus.

## Engineering interpretation

The answer is **three separate tolerances**, not one magic value. A practical speed experiment uses 0.01-mm adaptive curve sampling; numerical length 1e-5 and mesh floor 0.01 give smaller additional savings but change the failing downstream topology. The measured fastest stage4 is still ~21 s, far from subsecond. Further major improvements must eliminate repeated exact SSI/NURBS reprojections and reuse solved UV/3D traces across phase boundaries rather than continuing to loosen geometric identity.

## 2026-11-13 continuation: expensive SSI division and imprint reuse

This section is **experimental and uses the same real source fixture**, the same 47-face stage-4 assertion, default Rust 1.99.1 `cargo test` debug/test profile, and normal sccache. The optional `MT_SSI_WORKLOAD=1` instrumentation is cherry-picked from existing diagnostics commit `224efafe`; its counters are *entire-test* including the known failing stage5, not stage4-only.

The prior `MT_EXPERIMENT_CURVE_SAMPLE_FLOOR_MM` affects `monstertruck_modeling::Curve::parameter_division`, but **did not affect** `monstertruck_geometry::SurfaceCurve::parameter_division` / `IntersectionCurve::parameter_division`. Those lower-level adaptive divisions continued using a much tighter numeric length; each exact evaluation performs NURBS search + 4D Newton unless at an endpoint. A separate experiment-only `MT_EXPERIMENT_SSI_CURVE_SAMPLE_FLOOR_MM` now raises their adaptive chord tolerance, keeping the actual Boolean geometric admission length unchanged.

| SSI curve sampling floor (mm) | Stage 4 (s) | Full replay SSI division samples | Full replay exact double projections | Full replay NURBS presearches | Stage4 |
| --- | ---: | ---: | ---: | ---: | --- |
| Unmodified (effective numeric 1e-5) | 21.379 | 47,616 | 208,513 | 206,843 | 47-face PASS |
| 0.01 | 3.479 | 3,855 | 52,828 | 26,756 | 47-face PASS |
| 0.03 | 3.207 | 2,484 | 48,067 | 21,347 | 47-face PASS |
| 0.1 | 3.080 | 2,132 | 46,849 | 19,967 | 47-face PASS |
| 0.3 | 3.089 | 2,009 | 46,422 | 19,514 | 47-face PASS |
| 1.0 | 3.022 | ~2,000 | ~46,000 | 19,227 | 47-face PASS |

All rows also use `MT_EXPERIMENT_CURVE_SAMPLE_FLOOR_MM=0.01`, actual numeric length `MT_EXPERIMENT_CLASSIC_NUMERIC_LENGTH=1e-5`, and mesh floor `MT_EXPERIMENT_BOOLEAN_MESH_FLOOR_MM=0.01`. Like the preceding test, stage5 **still fails** (unknown-face classification). Timings vary by load and instrumentation; the major speedup is supported independently by **~4x fewer double projections and ~8x fewer NURBS presearch queries**.

Independent OCC analysis: exported 0.03-mm and 0.1-mm SSI-floor stage4 STEP solids, both `BRepCheck_Analyzer` valid one-solid 47-face, and both match the preceding 21-second fast candidate by **exactly the same OCC-integrated volume 7.392842432812021 mm³** and two directed OCCT `BRepAlgoAPI_Cut` residual volumes **0.0 / 0.0 mm³**. These are geometry-equivalence checks at OCC resolution, not an exact oriented-trim or analytic deviation certificate. Importantly **the SSI sampling floor is NOT currently a production-safe setting**: `parameter_division` is also invoked by some sampled curve-coincidence checks, so raising its floor changes coverage of those admissibility samples even when numerical thresholds are unchanged. Need separate discretization-only API before merging a non-experimental default.

Stack instrumentation at 10k-call intervals identified a high-cost surviving call path: `loops_store::reuse_imprinted_coincident_edge` inside `imprint_edges_on_faces`. This path evaluates and reciprocally reprojects exact curve samples; additional work remains even after making SSI sampling as coarse as 1 mm, so further tolerance loosening has negligible benefit (3.08 s → 3.02 s over 10x change).

Structural work now applied in the same experiment branch: precompute the exact `edge.oriented_curve().evaluate(midpoint)` **once per unique source edge** outside the inner target-face loop in `imprint_edges_on_faces` while preserving the precise same candidate samples for each face. On the same 0.1 mm SSI-floor fixture, the repeated call count improved from 19,967 to **19,037 NURBS presearches**, SSI evaluations from 60,421 to **58,083**, double projections from 46,849 to **44,511**. Stage4 remained a 47-face valid solid; its measured 3.817 s is higher than the earlier 3.080 s due to test contention, so **no wall-time speedup is claimed** for that isolated micro-optimization.

Next high-payoff structural change: stop recomputing complete two-direction SSI curve-coincidence evidence on every repeated face/imprint pairing. Cache **verifiable** geometric relation results for identical exact edge/trim identities within one Boolean operation, or establish certified identical carrier-and-interval proofs; preserve complete validation for unrelated same-endpoint thin lenses. More targeted CPU profiling may help: `perf record` on Asus is blocked by `kernel.perf_event_paranoid=2` for user-space monitoring; the correct root fix would be to change the sysctl for an authorized profiling session rather than fake a profile.

Even the fastest tested under 3.1 s stage4 is **not under 1 s** in the current unoptimized Rust test configuration. No production solver defaults are changed. All new environment controls are local to this experimental branch, not PR104. The original PR104 worktree remains untouched.
