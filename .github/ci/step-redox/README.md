# Legacy native OpenCascade diagnostic tools

The release fixture and semantic gates use Rust and an independent embedded
OCCT-WASM kernel. See the [validation guide](../../../step-redox/validation/README.md)
and [validator CLI](../../../step-redox/tools/occt-wasm-validator/README.md).

The remaining native OCP scripts are diagnostic and differential-oracle
utilities, not release dependencies or geometry implementations. They use
uv with the local pyproject.toml and uv.lock. The native harness is retained
during qualification of the migrated release gate. Independent Boolean
occupancy is provided by the Rust validator's occupancy subcommand.

Do not relax repository-wide strict Python policy for these tools.

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
