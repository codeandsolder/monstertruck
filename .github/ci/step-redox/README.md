# STEP CAD corpus and OpenCascade CI tooling

This directory is CI-only Python tooling, not a geometry implementation or
runtime dependency. It uses uv with the local pyproject.toml and uv.lock.

- harness.py validates committed STEP geometry fixtures. The Rust step-validate-semantics binary independently validates CLI JSON expectations.
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

## Tight absolute paired-volume diagnostics

For two independently imported single-solid STEP files, including differently
subdivided reconstructed faces:

```sh
uv run --frozen --project .github/ci/step-redox python \
  .github/ci/step-redox/occt_compare_solids.py source.step reconstructed.step \
  --single-solid
```

The single-solid and source-entity comparison modes report `paired_absolute_volume`
in addition to the existing independent validity, directed cut and adaptive
volume results. OCCT's whole-solid routines choose integration references from
vertex occurrence averages. Face splitting can change this reference and
therefore the measured volume even for matching boundaries. The paired diagnostic
integrates signed exact surface faces (including inner shells) using the same
four geometry-derived locations for both operands: the joint bounding-box
center and the three upper axis face centers. It uses Gauss–Kronrod with spline
spans at relative quadrature targets 1e-12 and 1e-13, summing face contributions
with compensated summation. Every paired delta and individual quadrature drift
must meet the existing absolute volume tolerance, 1e-10 mm³. No best-reference
selection, tolerance expansion, tessellation, or omitted invalid faces is used.

Individual volume spread across references remains visible: source and candidate
can share a trim-closure residual even when their paired volumes agree. This is
an independent numerical diagnostic, not an exact boundary certificate and not
an override for the existing `geometry_equal` decision. A passing volume test
alone cannot detect equal-volume changes to shape or placement.

Focused analytic/subdivision/cavity/negative tests:

```sh
cd .github/ci/step-redox
uv run --frozen python -m unittest -v test_occt_paired_volume
```
