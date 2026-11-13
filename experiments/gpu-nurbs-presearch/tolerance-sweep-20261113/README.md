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

Original optimized branch (`perf/ssi-reuse-division-samples-20261113`):
stage4 ~70.921s, 951,774 NURBS seeds, 949,338 3D Newton, 1,005,172
4D Newton. With a 0.001mm sampling floor: stage4 ~62.966s, 854,934
NURBS seeds, 852,498 3D Newton, 908,332 4D Newton, same 47-face
stage4 and same stage5 failure. Wall times fluctuate with host load.

All commands use Rust 1.99.1 and the system's normal sccache-integrated
Cargo wrapper. Do not treat microbench timings as production-optimized
performance. Do not merge experimental environment knobs as permanent
public API; the eventual solution needs a coherent error-budget policy.
