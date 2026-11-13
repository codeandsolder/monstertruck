# Reproducible Boolean tolerance / cost sweep (2026-11-13)

**Experiment only.** Intended to separate mesh/curve approximation density
from accepted geometric error and solver precision. No defaults are changed.

- `MT_EXPERIMENT_CURVE_SAMPLE_FLOOR_MM`: optional minimum adaptive chord
  tolerance in `monstertruck_modeling::Curve::parameter_division` (unit: mm
  for this fixture). This changes how many samples are generated, **not**
  the caller's original exact-incidence/weld acceptance threshold.
- `MT_EXPERIMENT_BOOLEAN_TOL_MM`: optional modeling tolerance in the
  single source-derived constructive-replay test, default `1e-5` mm.
  This controls the Boolean operation's actual requested tolerance and
  is intentionally independent of sampling.
- `MT_PROFILE_CURVE_DIVISION=1`: trace per-division effective versus
  requested tolerances, sample count and time. Extra output belongs in
  this local experiment folder.
- `MT_NEWTON_PROFILE=1`, `MONSTERTRUCK_PRESEARCH_TRACE=1`: log
  solver counts and timings. These are already present in the CPU branch.

Fixture: `step-redox` test
`mixed_complete_constructive_ir_replays_without_intermediate_step_import`.
Stages 1-4 produce 17/31/37/47 faces, then stage 5 fails due to
an **independent existing** oriented shell defect. The test therefore
normally has a nonzero exit code; evaluate whether the stage-4 47-face
solid was produced and the reported numerical work. Face count and
manifold status alone do not prove geometric equivalence.

## Complete sampling-only sweep (Asus, same source fixture)

| Curve sampling floor | Stage4 wall | NURBS searches | 3D Newton calls | 4D Newton calls | Stage4 |
| --- | ---: | ---: | ---: | ---: | --- |
| No override | 70.921s | 951,774 | 949,338 | 1,005,172 | 47-face PASS |
| 0.001mm (1µm) | 62.966s | 854,934 | 852,498 | 908,332 | 47-face PASS |
| 0.01mm (10µm) | 39.556s | 440,053 | 437,974 | 492,369 | 47-face PASS |
| 0.1mm (100µm) | 38.349s | 439,358 | 437,302 | 491,560 | 47-face PASS |

The original stage5 invalid oriented shell still blocks the remainder of
the replay; face count and manifoldness do NOT establish geometric fidelity.

The 0.01mm -> 0.1mm result has essentially **no further effect** on the
NURBS workload: only 695 searches saved by a further 10x coarsening.
The modeled-curve division trace reported 2,498 divisions producing
6,655 points at 0.1mm, versus roughly 439,000 NURBS searches. This is
not a simple over-tessellation problem. The remaining workload is driven by
repeated nonlinear SSI reprojection of curves elsewhere in the Boolean
pipeline, plus other topology operations.

## Modeling-tolerance axis

On Acer (do not cross-compare wall times with Asus), raising the requested
modeling tolerance from 1e-5mm to 0.001mm (100x) left the total Newton
counts **identical** (949,338 3D and 1,005,172 4D; 951,774 seeds).
The model's outward-facing tolerance is not driving most of those repeated
internal calculations. Increasing the requested modeling tolerance again to **0.01mm (1000x baseline)**
on the same Acer host also left the entire call counts **identical**,
with the expected 47-face stage4 solid and the unchanged stage5 error.
At 0.001mm modeling tolerance stage4 was 78.352s; at 0.01mm
it was 75.680s, variation without an actual reduction in work.

All these measurements use the unoptimized Rust test profile and substantial
host-load variation. Compare solver counts to distinguish actual work removal
from timing noise. Do not infer that any coarse sampling option is ready for
production without a geometric-error comparison.
