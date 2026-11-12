# Optional real-kernel CUDA FP64 seed bridge

This is an **experimental** live integration into the existing
`NurbsSurface::presearch_separable` cached-sample path. It is not the
production default and does not accelerate Newton refinement.

With `cuda-presearch` compiled and `MONSTERTRUCK_PRESEARCH_CUDA=1` set, the
kernel retains each 51×51 grid in VRAM, sends each query to a persistent
CUDA stream, and receives an FP64 minimum-seed index. It **always** computes
and compares the exact existing CPU FP64 answer before returning a seed to
the CAD solver, because GPU arithmetic can differ on numerical tie boundaries.
Any CUDA error falls back to the CPU result. This is a correctness-first
instrumented path, not yet a faster GPU algorithm.

On Asus (RTX 2060 Max-Q, CUDA 13.4):

```sh
LD_LIBRARY_PATH=/opt/cuda/lib64 MONSTERTRUCK_PRESEARCH_CUDA=1 \
  MONSTERTRUCK_PRESEARCH_TRACE=1 cargo +1.99.1 test \
  -p monstertruck-geometry --features cuda-presearch \
  --lib gpu_grid_matches_ordered_fp64_distance_scan -- --nocapture
```

The CUDA kernel is built into Cargo's ephemeral `OUT_DIR` only when the
feature is enabled; plain Rust builds never invoke NVCC. GPU query telemetry
is included in `NURBS_PRESEARCH_METRICS` (query count, mismatches, cumulative
time). Neither output binaries nor raw success payloads are persisted.

The **bottleneck now is about 1.16 million iterative surface-projection calls,
not the seed scan**. Real GPU acceleration requires batching the Newton
system solves across independent projections while retaining validated exact
incidence and deterministic topology. The current synchronous per-query
bridge demonstrates that the hardware can be accessed from live CAD but
is not the batch-refinement solution.
