# OCCT-WASM validation bindings

This patch targets the exact published occt-wasm 7.2.2 source revision
`eb94d1ad3b5f57ad90a671386f4ff804b34c3bd5`, OCCT submodule
`c7f60ffe85f28c67358d1771454044f0521f2fdd`. It is not a patch to the
geometry algorithm: it expands the independent reference kernel's public
validation facade.

- Topological occurrence enumeration/counting uses TopExp_Explorer, not
  TopExp::MapShapes (which deduplicates occurrences).
- Optimal bounding boxes expose both triangulation and shape-tolerance flags.
- Signed face integrals expose the actual OCCT BRepGProp_VinertGK result
  at a caller-chosen reference and epsilon, splitting spline integration spans.
  The flattened result contains alternating signed mass and relative error.
  Valid, closed, single-solid input is required; failed/nonfinite quadrature
  is an error, not a zero contribution.
- Point classification preserves OUT/IN/ON/UNKNOWN as 0/1/2/3 respectively.
- Double-precision triangle coordinates preserve the native mesh hash quantum.
- Binary-safe STEP import does not require UTF-8 transcoding.
- Strict kernel integration tests exercise all seven new bindings.

The upstream binding generator owns the C++, WASI and Rust glue. Run the
dedicated pinned-builder workflow to regenerate the embedded Brotli WASM
and test the actual kernel. A Cargo patch must select a rebuilt crate, not
just regenerated Rust stubs against the original WASM bytes.

The GitHub installation cannot create an upstream fork (403). Until a fork
is available, keep this small upstream patch here; the workflow produces a
tested packaged crate suitable for an isolated path-based Cargo patch.
The active geometry oracle remains unchanged until numerical parity passes.

The validator's new commands (with the rebuilt crate selected) are:

```sh
step-redox-occt-wasm-validator analyze source.step
step-redox-occt-wasm-validator paired-volume source.step candidate.step --absolute-tolerance 1e-10
```

Paired volume prints the complete shared-reference report and returns nonzero
if either the volume delta or the quadrature drift exceeds tolerance.
The four reference locations are the midpoint and three positive-axis box
faces. Each reference is evaluated at 1e-12 and 1e-13. Signed per-face terms
are summed by a non-overlapping expansion with final nearest-even rounding.
The strict mesh hash retains the native oracle's 1e-9 mm (one picometer)
quantization, winding/order independence, and SHA-256 byte packing.

Run hostile or malformed fixture inputs under the existing runtime-limits
wrapper with the explicitly exported executable, e.g.
`validate-occt-wasm.sh /path/to/validator 180 6291456 analyze source.step`.
This keeps runtime limits separate from compilation and avoids reusing an
arbitrarily stale executable in a local target directory.
