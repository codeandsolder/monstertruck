# STEP validation harness

This directory is the regression suite for step-redox.

The fixture manifest records public download URLs, SHA-256 checksums, and source
notes. Raw upstream STEP files are not committed; compressed fixture bundles
provide reproducible offline acquisition. The Rust validator verifies cached
files and bundled decompression results before accepting them, and uses the
recorded download URLs only when the bundle and verified cache are unavailable.

Each case runs step-redox or the semantic count-generator CLI and validates
the result with several independent signals. The same manifest therefore covers
ordinary optimizer roundtrips and cross-model generation. The strongest
pair now exercises the same recovered header grammar in both directions:
48 contacts → 72 contacts and 72 contacts → 48 contacts, each compared with an
independently published sibling model.

- OCCT import and BRepCheck validity
- solid, shell, face, edge, and vertex counts
- bounding-box center and extents
- compound-level volume, surface area, total edge length
- per-solid topology, position, extents, volume, and area
- summed per-solid material volume/area for multi-solid assemblies
- center of mass and inertia tensor
- deterministic software renders from +X, -X, +Y, -Y, +Z, -Z, and two isometric views
- render mean absolute error, changed-pixel fraction, and silhouette XOR

For mapped multi-solid assemblies, OpenCascade's one-shot compound volume can differ by several ppm even when every imported solid matches independently to numerical noise. The harness therefore keeps the compound integral visible but uses the summed per-solid volume as the stricter occupied-material invariant.

The independent embedded OCCT-WASM kernel performs tessellation. A Rust
software rasterizer generates the reference, candidate and difference images
without GPU or display-server dependencies. See the
[validator documentation](../tools/occt-wasm-validator/README.md).

Thresholds are explicit per test case in fixtures.json. Exact topology can be
required independently from geometric tolerances. A case comparing two
independent vendor models may use bbox-center alignment to permit a known rigid
translation without permitting scaling or shape changes.

Local invocation from the workspace root:

    cargo build --release -p step-redox --features cad-kernel-monstertruck \
      --bin step-redox --bin step-count-resize --bin step-chain-resize --bin step-validate-semantics
    cargo build --release --locked \
      --manifest-path step-redox/tools/occt-wasm-validator/Cargo.toml
    step-redox/tools/occt-wasm-validator/target/release/step-redox-occt-wasm-validator run \
      --manifest step-redox/validation/fixtures.json \
      --cache step-redox/.cache/step-redox-validation \
      --out step-redox/validation/out \
      --step-redox target/release/step-redox \
      --step-count-resize target/release/step-count-resize \
      --step-chain-resize target/release/step-chain-resize

The executable paths are explicit: use your actual exported artifact paths
when building with a managed Cargo wrapper. The runner skips disabled
manifest cases; repeat `--case NAME` for a targeted subset.

After fixture acquisition, run the independent semantic probe:

    target/release/step-validate-semantics \
      --manifest step-redox/validation/fixtures.json \
      --cache step-redox/.cache/step-redox-validation \
      --out step-redox/validation/out/semantic-probes \
      --step-redox target/release/step-redox

It checks the cached fixture hashes independently before testing periodic-chain
topology, CAD proof provenance, instance patterns, and formed-sheet metadata.

The release workflow retains the aggregate report.json, summary.md, per-case
reports and reference/candidate/difference images as an artifact.

Fixture policy:

1. Prefer public first-party manufacturer CAD URLs when practical.
2. Otherwise use stable public first-party distribution/model endpoints.
3. Never commit the downloaded STEP files.
4. Pin every fixture by SHA-256.
5. Keep enough provenance in source_note to locate/review the source later.

## Constructive recovery coverage survey

From the unified Monstertruck workspace root, run a batch of STEP documents in
one kernel-backed process:

    cargo +stable run -p step-redox --bin step-recovery-survey \
      --features cad-kernel-monstertruck -- \
      --output /tmp/step-recovery-survey.json \
      step-redox/validation/fixtures/native_staged_rectangular_loft.step \
      step-redox/validation/fixtures/native_mixed_treatment_4a6356.step \
      step-redox/validation/fixtures/production_radial_slot.step

Without `--output`, the command writes JSON to stdout. Each input produces a
record (even if reading/parsing fails) with a SHA-256 for successfully read
bytes, elapsed time, body-level constructive/semantic/unrecovered counts,
recovered operation families, and original source IDs for unrecovered solids.
Aggregate counters distinguish fully and partly constructive documents from
semantic-only, wholly unrecovered, empty, and failed inputs.

This survey measures **automatic constructive IR selection**, not an independent
source-exact equality oracle. A recovered radial-slot body can be editable
while its exported B-rep differs from the original manufacturer geometry.
Track independent oriented analytic-boundary proofs, OCCT differential checks,
and explicit source-healing experiments separately. The checked-in `native_*`
fixtures contain a mixture of real-source captures and generated regressions;
never interpret their success rate as an unbiased JLC/catalog success rate.
For real manufacturer surveys, preserve input provenance and compare the
reported SHA-256 to the pinned hashes in `validation/fixtures.json`.

When `--output` is provided, a sibling `*.cases.jsonl` checkpoint journal is
also written. Each complete case is synchronized to disk before the next source is
attempted. Fresh runs refuse to overwrite either an existing report or an
existing journal. To resume an interrupted run, repeat the same command with
`--resume`. The checkpoint header binds the ordered absolute input paths and
the SHA-256 of the executable to the run; completed source files are rehashed
before their checkpoint results are reused. Changed source bytes or an
updated executable require a new output path. Only a torn, unterminated final
journal record is discarded; complete corrupt records fail closed. A file lock
prevents concurrent writers to the same journal. The aggregate report is
published atomically without replacing an existing file after every input has
a completed result. The journal format is identified by the v2 schema header.

## Source spline-to-carrier incidence audit

For source documents that fail independent boundary qualification, inspect the
actual authored spline edges against their own literal plane, cylinder and
cone carriers:

    cargo run --release -p step-redox --features cad-kernel-monstertruck \
      --bin step-spline-edge-audit -- input.step SOLID_ID --tolerance-mm 1e-5

The report resolves the solid's authored length unit, includes its declared
representation uncertainty when available, and brackets the maximum
whole-edge distance using positive-weight Bernstein control bounds. Complete
open edge coverage requires clamped continuous knots and unique endpoint
visits; a valid curved edge need not have a globally monotone chord projection.
Unsupported trim domains and unexamined face/curve families remain explicit.
The caller's tolerance classifies the reported bounds only: it does not
change the constructive recovery gate, source geometry or source uncertainty.

For explicit source-distinct repairs of a proven defective crossing-cylinder
spline, preserve the original carrier surfaces and edge/coedge/vertex
identities and construct a bounded finite 3D leader for the analytic branch:

    cargo run --release -p step-redox --features cad-kernel-monstertruck \
      --bin step-heal-source-edges -- original.step repaired.step \
      --endpoint-tolerance-mm 1e-8 --target-gap-mm 1e-9 \
      cylinders --repair 0:171:1049 --repair 0:171:1063 \
      --leader-tolerance-mm 2e-10 --max-leader-segments 8192

The command refuses existing destinations. The manifest records the full
source-to-analytic separation bracket, Hermite approximation and serialized
control-hull bounds, and all original and new entity identities. The new
document is independently checked as a new source; it is not admitted as
equivalent to the original defective source. The same CLI's `hyperbolas`
subcommand accepts the five-field DATA_SECTION:SOLID:EDGE:PLANE_FACE:CONE_FACE
repair specification and the two common error budgets.
