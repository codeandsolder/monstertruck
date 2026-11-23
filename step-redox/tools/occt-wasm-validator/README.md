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
