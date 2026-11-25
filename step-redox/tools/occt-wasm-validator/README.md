# Independent STEP validator using OCCT-WASM

The standalone Cargo workspace consumes the independent OCCT kernel directly
from [occt-wasm-truck](https://github.com/codeandsolder/occt-wasm-truck).
The fork contains and tests the embedded WASI kernel; consumers do not build
OCCT or apply source patches.

Until the first upstream-synchronized release passes its release-age gate,
the dependency is pinned to the validated fork revision. Subsequent releases
use stable SemVer tags and Renovate.

Build from the Monstertruck workspace root with the repository toolchain:

```sh
cargo build --manifest-path step-redox/tools/occt-wasm-validator/Cargo.toml --locked --release
```

Use the resulting `step-redox-occt-wasm-validator` executable:

```sh
step-redox-occt-wasm-validator analyze source.step
step-redox-occt-wasm-validator analyze-batch --mesh source1.step source2.step
step-redox-occt-wasm-validator compare-metrics reference.json candidate.json --thresholds thresholds.json
step-redox-occt-wasm-validator compare reference.step candidate.step --thresholds thresholds.json --out comparison
step-redox-occt-wasm-validator mesh-hash source.step
step-redox-occt-wasm-validator paired-volume reference.step candidate.step --absolute-tolerance 1e-10
```

The geometry report contains actual STEP transfer counts, repeated topology
occurrences, signed default BRepGProp mass properties, loose bounds and the
per-solid inventory. Batch analysis compiles the WASM module once and uses
isolated instances for each input.

Thresholds and case selection are described in
[the validation manifest documentation](../../validation/README.md).
The full comparison uses the same geometry tolerances and alignment modes
as the manifest. It produces eight reference, candidate and difference
views, a contact sheet, and per-view image metrics. The Rust software rasterizer
does not guarantee bit-identical Pillow pixels. Strict mesh hashing uses
double-precision world-space triangles and nearest-even quantization at
1e-9 mm. The digest report includes the original triangle count and SHA-256.
Enable the same measurement in analyze or analyze-batch with --mesh.
Exact digests are kernel-version-specific: native and WASM builds may produce
different coordinates at identical triangle counts. Compare these digests only
within the same kernel build and tessellation settings.
Paired volume integrates signed faces at two quadrature targets
and four common references.

For independent Boolean occupancy:

```sh
step-redox-occt-wasm-validator occupancy --operation difference \
  --left host.step --right cutter.step --result result.step \
  --point 0.1 0.2 0.3 --report occupancy.json
```

Use repeated `--point X Y Z`, a `--points-json` file, or `--box XMIN XMAX
YMIN YMAX ZMIN ZMAX --samples N --seed 42`. Operations are difference,
intersection, union and XOR. Boundary witnesses are inconclusive;
unknown classification is separately reported. Exit 1 means an actual
violation and exit 2 means insufficient conclusive evidence or unknown
classification. The report records all witness coordinates for exact replay.
The seeded sampler uses ChaCha12 rather than Python's Mersenne Twister.
Finite occupancy evidence is not a proof of complete geometric equivalence.

The release runner acquires SHA-256-verified fixtures, expands generator
arguments without a shell, reuses one compiled OCCT module for all cases,
and writes per-case and aggregate reports. See the validation README for
the complete invocation. Use `fetch-fixtures` for acquisition without
kernel initialization and repeat `--case NAME` to select a subset.

For explicitly bounded execution, see
`step-redox/tools/validate-occt-wasm.sh`.

For per-solid Boolean volume diagnostics, use:

```sh
step-redox-occt-wasm-validator boolean-residuals reference.step candidate.step \
  --absolute-tolerance 1e-9
```

This reports both directed cut volumes and the missing volume relative to
the common intersection for each matched solid occurrence. Matching uses
the same deterministic geometry metric correspondence as the release
validator, without truncating uneven inventories. The report identifies
both occurrence indices. These are diagnostics, not an exact boundary
certificate or a proof of globally optimal material correspondence.

## Exact single-solid diagnostic comparison

Compare the same individual solid across independently exported STEP files,
including reconstructed solids with different face subdivisions:

```sh
step-redox-occt-wasm-validator solid-diagnostics source.step candidate.step
```

For a source containing repeated or overlapping solids, select the exact
Part 21 textual entity number (the `#123` label, **not** OCCT's model rank):

```sh
step-redox-occt-wasm-validator solid-diagnostics full-source.step recovered.step \
  --reference-entity-id 123
```

`--candidate-entity-id` selects a source entity in the other operand
independently. Without these options, each entire file must transfer to
exactly one solid. Entity selection uses OCCT's parsed STEP model, rejecting
missing, ambiguous, or nontransferable identifiers.

The JSON includes signed default BRep properties, topology, center/bounds,
native whole-solid adaptive Gauss–Kronrod volumes and relative error
estimates, both directed Boolean cut volumes, and the shared-reference
paired absolute face-integral diagnostic. `geometry_equal` requires both
valid solids, consistent adaptive volumes, and both material cuts below
the existing absolute/relative policy. The stricter paired-volume result is
reported **separately**, not used to override that decision. Face-count
equality is deliberately not required. Exit status 0 indicates equality,
1 a measured mismatch, and an execution error fails with an error message.
The tolerances and adaptive integration targets retain the previous native
oracle's definitions; the new implementation uses the independently built
WASM kernel without native Python/OCP.

For diagnosing a nonzero paired-volume residual, add `--face-flux` to
`solid-diagnostics` or `paired-volume`. The optional detail is part of the
paired-volume report and reuses exactly its integration terms. This includes the ordered, oriented shell/face occurrences,
carrier types, bounds, UV bounds, surface areas, and signed integrals and
quadrature error estimates at the four shared reference points for both
quadrature tolerances. The report
also sums the signed contributions by carrier type using the same accurate
summation as paired volume. Carrier families can change under exact
reconstruction, so a family difference alone is not a geometric error or a
face correspondence. The optional report does not change acceptance policy.

For spherical faces, the optional geodesic diagnostic reports an independent
analytic signed area and volume flux for a convex single-wire polygon with
minor great-circle edges. Native circles are checked analytically; finite
nonperiodic rational B-spline trims are refined into homogeneous Bézier spans,
then the full span's plane, squared radius and angular progression are checked
using Bernstein coefficients. The bounds are computed in floating point,
not interval arithmetic. Unsupported or ambiguous patches include an explicit
reason. This result is diagnostic only and does not replace the quadrature
volume or change the acceptance gate.

When at least one spherical face can be independently evaluated, the report
also includes `analytic_substitution`: the mixed analytic and quadrature
volume at each reference, its direct signed-term difference, quadrature
drift, reference spread and substituted/unsupported face counts. The
original quadrature report and exit status remain unchanged.

The spatial spherical boundary diagnostic additionally integrates the actual
3-D trim rather than the imported UV pcurve. Native circle trims are
represented as exact rational quadratic spans; finite clamped rational
B-splines retain their homogeneous Bézier spans. Full-domain Bernstein
squared-radius bounds, directed wire closure by native vertex identity, and
two independent certified nonsingular gauges are required. Tiny gaps between
trim endpoints and successive spans are closed by explicit radially projected
short chords; their total normalized length is reported separately. Solid angle is recovered modulo the full
sphere using the actual face orientation, and vector area follows from
Stokes' boundary integral. Two quadrature tolerances, their estimated errors
and gauge disagreement are reported separately. Multiple wires, small-circle
trims and complementary spherical regions are supported; ambiguous empty/full
regions and unsupported curves remain explicit rejections. The
`spatial_substitution` diagnostic combines these fluxes with the original
quadrature for other faces, without changing the raw gate or its exit status.
Its separate `observed_within_tolerance` diagnostic requires every spherical
face to be covered and the measured paired difference, coarse/fine drift,
and shared-reference spread (plus the summed estimated spatial quadrature
error) to fit the requested absolute tolerance. This is not a rigorous
interval enclosure or an alternative process exit criterion.
