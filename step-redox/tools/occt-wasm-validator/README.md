# Independent STEP validator using OCCT-WASM

This separate Cargo workspace consumes the patched OCCT-WASM crate directly
from [occt-wasm-truck](https://github.com/codeandsolder/occt-wasm-truck). The
fork contains and tests the actual embedded WASI kernel; consumers do not build
OCCT or apply source patches.

The draft uses the exact validated fork revision while the initial stable
release is quarantined. Use a stable fork SemVer tag once it exists
(`v7.3.0`, then subsequent upstream crate releases via Renovate).

```sh
cargo build --manifest-path step-redox/tools/occt-wasm-validator/Cargo.toml --locked --release
step-redox-occt-wasm-validator analyze source.step
step-redox-occt-wasm-validator mesh-hash source.step
step-redox-occt-wasm-validator paired-volume reference.step candidate.step --absolute-tolerance 1e-10
```

The analyze JSON follows the native OCP harness schema, including actual
STEP transfer root counts, repeated topology occurrences, signed default
BRepGProp mass properties, loose bounds and the per-solid inventory.
Mesh hashing uses native double world-space triangles and nearest-even
quantization at 1e-9 mm. Paired volume uses independent signed face
Gauss-Kronrod integration with two tolerances and four shared references.

For bounded execution, use `step-redox/tools/validate-occt-wasm.sh` with the
explicit built validator path, timeout and virtual-memory limit. The native
oracle remains available while differential parity is established.
