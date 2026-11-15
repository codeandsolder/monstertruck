# STEP CAD corpus and OpenCascade CI tooling

This directory is CI-only Python tooling, not a geometry implementation or
runtime dependency. It uses uv with the local pyproject.toml and uv.lock.

- harness.py and semantic_probe.py validate the committed STEP fixtures.
- occt_* utilities provide independent OpenCascade geometry comparisons and
  corpus diagnostics. They do not participate in constructive recovery.
- Core STEP parsing, feature detection, geometry and topology belong to Rust.

From the repository root, call uv run --frozen --project .github/ci/step-redox
python .github/ci/step-redox/harness.py, passing fixture manifest, executable
paths and report directories. Run full corpus validation on release/manual CI.
Do not relax repository-wide strict Python policy to accommodate this tooling.


## Source-independent Boolean occupancy checks

`occt_boolean_occupancy.py` verifies a Boolean result using **three independently
imported STEP shapes**: its left operand, right operand, and candidate result.
It checks the actual material set operation, not surface/face counts, volume
integrals, or the candidate kernel's own point classifier. Supports difference,
intersection, union, and XOR, including multi-solid operands. A point on an
OCCT-classified boundary is explicitly inconclusive, never coerced to inside or
outside. At least one unambiguous witness is required to pass.

```sh
uv run --frozen --project .github/ci/step-redox python \
  .github/ci/step-redox/occt_boolean_occupancy.py \
  --operation difference --left host.step --right cutter.step \
  --result result.step --point 0.1 0.2 0.3 \
  --report occupancy-report.json
```

Use repeated `--point X Y Z`, `--points-json points.json`, or a deterministic
`--box XMIN XMAX YMIN YMAX ZMIN ZMAX --samples N --seed 42` to sample the
actual solid material. Exit status 1 means a violated Boolean witness; 2 means
no conclusive witness. This is a *finite independent correctness regression*,
not a proof of full occupied-volume equivalence; for that use source-matched
carrier/trim and independent volume certificates in addition.
