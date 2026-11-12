# NURBS presearch on GPU? Measurement before integration

**Conclusion (2026-11-13):** A CUDA f64 nearest-grid scan has useful throughput
only for batches of dozens or more. A cached **f64 KD-tree on the CPU beats
CUDA f64**, even for thousands of queries, without changing the sequential
inverse-projection API. A CUDA f32 scan is much faster, but silently changes
nearest seeds on near-tie queries. Do **not** wire an unchecked f32 result into
CAD geometry or replace existing numerical precision criteria.

## Reproduce on an NVIDIA GPU

On Asus: RTX 2060 Max-Q (SM 7.5, 6 GB), CUDA 13.4 and Ryzen 9 4900HS.

```sh
cd experiments/gpu-nurbs-presearch
/opt/cuda/bin/nvcc -std=c++17 -O3 -arch=sm_75 --fmad=false \
  -Xcompiler=-O3,-march=native,-ffp-contract=off,-fopenmp \
  -o benchmark benchmark.cu
OMP_WAIT_POLICY=ACTIVE GOMP_SPINCOUNT=1000000 \
  OMP_PROC_BIND=TRUE OMP_PLACES=cores ./benchmark
```

Benchmark code intentionally exercises the *same f64 cached-grid argmin*
with strict first-index tie resolution (original row-major grid index).
There are 2,601 precomputed grid samples; the GPU retains the grid in VRAM,
uses pinned host buffers for batches of queries/winners, and does an f64
tree reduction within each CTA. The f64 CPU KD-tree is built once per grid
and then reused. Its AABB pruning deliberately searches extra candidates
around rounding boundaries, so the original f64 argmin always wins.
It tests 8,192 points including near-tie stress cases. All GPU f64 and
CPU KD-tree results must agree with the brute-force CPU baseline.

The **f32 GPU experiment is explicitly non-equivalent**: it translates
coordinates to the local object origin before narrowing to f32 and uses
the same grid/query points. The latter half of the queries are close to
Voronoi tie boundaries, where seed changes are expected. Their mismatches
are measured but not permitted as CAD answers.

## Representative measurements

**All times are median wall-clock per batch, including host-to-device
copy, kernel launch, device-to-host copy, and synchronization for GPU.**
The RTX 2060 Max-Q is weak at f64. Asus was simultaneously running other
Rust builds (load average roughly 13--20), so CPU samples are noisy.
CPU parallel runs use persistent active OpenMP workers; the default passive
sleep configuration was dominated by scheduler delays and is not comparable.
Do not treat these as controlled production end-to-end speedups.

A representative run after adding the f64 KD-tree:

| Batch | Brute CPU1 | KD CPU1 | KD CPU8 | GPU f64 E2E | GPU f32 E2E |
|---:|---:|---:|---:|---:|---:|
| 1 | 6.32 us | 0.23 us | — | 34.80 us | 23.47 us |
| 32 | 176.75 us | 13.28 us | 3.81 us | 51.03 us | 23.39 us |
| 128 | 745.85 us | 56.23 us | 8.70 us | 101.76 us | 30.69 us |
| 512 | 3.17 ms | 0.27 ms | 0.028 ms | 0.33 ms | 0.043 ms |
| 8,192 | 47.74 ms | 4.18 ms | 0.80 ms | 4.70 ms | 0.40 ms |

KD-tree build: **1.01 ms once per grid** (2,601 points), not per query.
At batch 8,192: GPU f64 **0 incorrect**, CPU f64 KD-tree **0 incorrect**,
GPU local-origin f32 **1,885 differing winners** (mostly adversarial near-ties).
The isolated GPU f64 launch+sync costs ~25--28 us even without transfers.
The GPU f32 speed is an *upper bound for a heuristic candidate stage*, not
a verified answer at CAD numerical precision.

## What is and is not measured

The 51x51 grid positions are a synthetic smooth circular-arc extrusion
at roughly x=67 mm, y=5 mm, z=2 mm, **not a trace of 4a6356**.
This experiments only with selecting the nearest precomputed sample.
It does **not** evaluate the NURBS basis, compute Newton refinement,
verify exact source boundaries, or model Boolean face dependency scheduling.

The production kernel currently calls inverse search synchronously,
so getting the GPU's larger-batch throughput would require restructuring
calls by surface and query group. Further optimization should first
replace the full cached-grid scan with a generic f64 spatial index and
retain the current strict first-minimum and cache invalidation semantics.
Any GPU f32 candidate generator would require sound error-bound
certification of the nearest sample and a full f64 fallback on ambiguity.
