# `monstertruck-healing`

<!-- cargo-rdme start -->

Shape healing for solids imported from other CAD systems.

The compressed-shell repair passes that turn a loaded B-rep into something a
boolean kernel can work with: closed-edge and closed-face splitting, seam
stripping, pass-through-edge dedup, shell orientation normalization and opt-in
conservative naked-edge welding.

These passes are POST-CSG and kernel-independent -- nothing here references a
boolean backend, which is why the crate sits below both the published
`monstertruck-solid` marching kernel and any external SSI boolean backend
rather than inside either.

<!-- cargo-rdme end -->

> Forked from [`truck-shapeops`](https://crates.io/crates/truck-shapeops) v0.4.0 by [ricosjp](https://github.com/ricosjp/truck).

## Naked-edge welding

`snap_naked_edges` is an explicit, tolerance-bounded operation on a
`CompressedShell`; it is **not** part of the automatic STEP healing or Boolean
pipeline. It identifies singly used edge pairs, bounds their geometric
separation, and reports the remapping of old vertices and edges and any
candidates it declined to join. A failed repair leaves the shell unchanged.

The geometric proof assumes that the curve type supplies an exact homogeneous
B-spline representation. It does **not** independently certify that the
surviving shared edge lies on both incident face surfaces or reconstruct
2D trims. Do not use it as a replacement for exact trim/SSI validation or
silently apply it to provenance-carrying STEP geometry.

## License

Apache License 2.0
