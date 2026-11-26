//! Topologically global edge repair and simplification for inverse meshing.
//!
//! Split source edges at existing geometric T-junction vertices before building
//! edge incidence. After patching, remove straight-through vertices only if
//! their two incident *unique* boundary edges have exactly the same pair of
//! adjacent patches; this guarantees both sides of a seam receive the same cut.

use super::{EdgeIncidence, EdgeKey, InputPolygon, MeshRecoveryError, PlanarPatch, edge_key};
use crate::{InnerSpace, Point3, Vector3};
use rustc_hash::FxHashMap as HashMap;
use std::collections::HashSet;

// Error-free difference (high component, subtraction tail).
fn two_diff(a: f64, b: f64) -> [f64; 2] {
    let hi = a - b;
    let virtual_b = a - hi;
    let virtual_a = hi + virtual_b;
    [hi, (a - virtual_a) + (virtual_b - b)]
}

// Grow a nonoverlapping floating-point expansion by one scalar, using
// error-free addition. At most sixteen components can arise from the
// eight two-component products in a determinant of two differences.
fn add_expansion_term(expansion: &mut [f64; 16], len: &mut usize, mut carry: f64) {
    if carry == 0.0 {
        return;
    }
    let mut write = 0;
    for read in 0..*len {
        let component = expansion[read];
        let sum = carry + component;
        let virtual_b = sum - carry;
        let virtual_a = sum - virtual_b;
        let error = (carry - virtual_a) + (component - virtual_b);
        if error != 0.0 {
            expansion[write] = error;
            write += 1;
        }
        carry = sum;
    }
    if carry != 0.0 {
        expansion[write] = carry;
        write += 1;
    }
    *len = write;
}

fn difference_of_products(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> f64 {
    if [a[1], b[1], c[1], d[1]].into_iter().all(|tail| tail == 0.0) {
        // Compensated difference of products. FMA retains the low bits
        // lost by rounding the second product before subtraction.
        let product = c[0] * d[0];
        return a[0].mul_add(b[0], -product) + (-c[0]).mul_add(d[0], product);
    }
    let mut expansion = [0.0; 16];
    let mut len = 0;
    for (left, right, sign) in [(a, b, 1.0), (c, d, -1.0)] {
        for x in left {
            for y in right {
                let product = x * y;
                let error = x.mul_add(y, -product);
                add_expansion_term(&mut expansion, &mut len, sign * error);
                add_expansion_term(&mut expansion, &mut len, sign * product);
            }
        }
    }
    expansion[..len].iter().sum()
}

// Adaptive orientation guarantees a determinant's sign, not its relative
// magnitude. Normals of skinny triangles require the latter too. Compute
// a compensated cross product of the original coordinate differences.
pub(super) fn oriented_triangle_area(a: Point3, b: Point3, c: Point3) -> Vector3 {
    let first: [[f64; 2]; 3] = std::array::from_fn(|i| two_diff(b[i], a[i]));
    let second: [[f64; 2]; 3] = std::array::from_fn(|i| two_diff(c[i], a[i]));
    Vector3::new(
        difference_of_products(first[1], second[2], first[2], second[1]),
        difference_of_products(first[2], second[0], first[0], second[2]),
        difference_of_products(first[0], second[1], first[1], second[0]),
    )
}

fn projected_orientation(a: Point3, b: Point3, c: Point3, axis: usize) -> f64 {
    use robust::{Coord, orient2d};
    let project = |p: Point3| Coord {
        x: p[(axis + 1) % 3],
        y: p[(axis + 2) % 3],
    };
    orient2d(project(a), project(b), project(c))
}

/// The three projected 2D orientation tests certify collinearity in 3D
/// using adaptive exact predicates for the *stored f64 input coordinates*.
/// A fitted near-line is intentionally not a shared topological edge.
pub(super) fn on_line(a: Point3, b: Point3, p: Point3) -> bool {
    (0..3).all(|axis| projected_orientation(a, b, p, axis) == 0.0)
}

/// Returns the coordinate along the directed segment only for a strictly
/// interior vertex exactly collinear with the source f64 endpoints.
pub(super) fn on_segment(a: Point3, b: Point3, p: Point3) -> Option<f64> {
    if !on_line(a, b, p) {
        return None;
    }
    let delta = b - a;
    let length_squared = delta.magnitude2();
    if length_squared <= 0.0 || !length_squared.is_finite() {
        return None;
    }
    let along = (p - a).dot(delta) / length_squared;
    (0.0..1.0).contains(&along).then_some(along)
}

/// Input mesh polygons may share a geometric seam with different subdivisions.
/// Insert existing vertices on longer edges, preserving winding and geometry.
///
/// Three pre-sorted axes provide exact segment bounding-box range queries.
/// Select the sparsest nonzero coordinate interval, then reject candidates
/// outside the other two exact coordinate bounds before robust predicates.
/// This preserves exact T-junction detection while avoiding a quadratic scan
/// of unrelated vertices on densely tessellated CAD surfaces.
pub(super) fn split_source_t_junctions<'a>(
    positions: &[Point3],
    polygons: impl IntoIterator<Item = &'a mut InputPolygon>,
) {
    let polygons = polygons.into_iter().collect::<Vec<_>>();
    let mut used = vec![false; positions.len()];
    let mut spokes = vec![Vec::new(); positions.len()];
    for polygon in &polygons {
        for (i, &a) in polygon.vertices.iter().enumerate() {
            let b = polygon.vertices[(i + 1) % polygon.vertices.len()];
            used[a] = true;
            used[b] = true;
            spokes[a].push(b);
            spokes[b].push(a);
        }
    }
    for around in &mut spokes {
        around.sort_unstable();
        around.dedup();
    }
    let point_ids: Vec<usize> = used
        .iter()
        .enumerate()
        .filter_map(|(id, &used)| used.then_some(id))
        .collect();
    let sorted: [Vec<usize>; 3] = std::array::from_fn(|axis| {
        let mut indices = point_ids.clone();
        indices.sort_unstable_by(|&a, &b| {
            positions[a][axis]
                .total_cmp(&positions[b][axis])
                .then(a.cmp(&b))
        });
        indices
    });

    for polygon in polygons {
        let mut split = Vec::with_capacity(polygon.vertices.len());
        for i in 0..polygon.vertices.len() {
            let a = polygon.vertices[i];
            let b = polygon.vertices[(i + 1) % polygon.vertices.len()];
            split.push(a);

            let vector = positions[b] - positions[a];
            let nearest = (0..3)
                .filter(|&axis| vector[axis] != 0.0)
                .map(|axis| {
                    let first = positions[a][axis].min(positions[b][axis]);
                    let last = positions[a][axis].max(positions[b][axis]);
                    let candidates = &sorted[axis];
                    let begin = candidates.partition_point(|&id| positions[id][axis] <= first);
                    let end = candidates.partition_point(|&id| positions[id][axis] < last);
                    (axis, begin, end)
                })
                .min_by_key(|&(axis, begin, end)| (end - begin, axis));
            let Some((axis, begin, end)) = nearest else {
                continue;
            };
            let mut inner = sorted[axis][begin..end]
                .iter()
                .filter_map(|&id| {
                    let point = positions[id];
                    if (0..3).any(|coordinate| {
                        point[coordinate] < positions[a][coordinate].min(positions[b][coordinate])
                            || point[coordinate]
                                > positions[a][coordinate].max(positions[b][coordinate])
                    }) {
                        return None;
                    }
                    let parameter = on_segment(positions[a], positions[b], point)?;
                    spokes[id]
                        .iter()
                        .any(|&neighbor| on_line(positions[a], positions[b], positions[neighbor]))
                        .then_some((parameter, id))
                })
                .collect::<Vec<_>>();
            inner.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            split.extend(inner.into_iter().map(|(_, id)| id));
        }
        polygon.vertices = split;
    }
}

/// An undirected source seam is identified by the *pair of patches it joins*.
/// An exterior edge is considered to join the patch to the outside.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Seam {
    left: usize,
    right: Option<usize>,
}

fn seam_for_edge(edge: EdgeKey, incidence: &EdgeIncidence, membership: &[usize]) -> Seam {
    let uses = &incidence[&edge];
    let a = membership[uses[0].0];
    let b = uses.get(1).map(|use_| membership[use_.0]);
    match b {
        Some(other) if other < a => Seam {
            left: other,
            right: Some(a),
        },
        _ => Seam { left: a, right: b },
    }
}

/// Global degree-2 candidates are removable only on a straight seam carrying
/// the same two adjacent patches. A patch corner, a hole branch, or the end of
/// a different patch automatically becomes an anchor.
fn removable_boundary_vertices(
    patches: &[PlanarPatch],
    positions: &[Point3],
    incidence: &EdgeIncidence,
    membership: &[usize],
) -> HashSet<usize> {
    let mut unique_edges = HashMap::<EdgeKey, Seam>::default();
    for patch in patches {
        for cycle in &patch.boundary_loops {
            for i in 0..cycle.len() {
                let key = edge_key(cycle[i], cycle[(i + 1) % cycle.len()]).0;
                unique_edges
                    .entry(key)
                    .or_insert_with(|| seam_for_edge(key, incidence, membership));
            }
        }
    }
    let mut incident = HashMap::<usize, Vec<(usize, Seam)>>::default();
    for (&[a, b], &seam) in &unique_edges {
        incident.entry(a).or_default().push((b, seam));
        incident.entry(b).or_default().push((a, seam));
    }
    incident
        .into_iter()
        .filter_map(|(v, neighbors)| {
            let [(left, a), (right, b)] = neighbors.as_slice() else {
                return None;
            };
            if a != b || left == right {
                return None;
            }
            on_segment(positions[*left], positions[*right], positions[v]).map(|_| v)
        })
        .collect()
}

/// Verify each proposed collapsed chain against the final chord, not merely
/// against its local pair of segments. Revoke *the whole chain* if it bends,
/// so a run of many tiny departures cannot accumulate an unacceptable drift.
fn revoke_nonstraight_chains(
    patches: &[PlanarPatch],
    positions: &[Point3],
    removable: &mut HashSet<usize>,
) {
    loop {
        let mut revoked = HashSet::new();
        for patch in patches {
            for loop_ in &patch.boundary_loops {
                let anchors = (0..loop_.len())
                    .filter(|&i| !removable.contains(&loop_[i]))
                    .collect::<Vec<_>>();
                if anchors.len() < 3 {
                    revoked.extend(loop_.iter().copied());
                    continue;
                }
                for i in 0..anchors.len() {
                    let from = anchors[i];
                    let to = anchors[(i + 1) % anchors.len()];
                    let mut cursor = (from + 1) % loop_.len();
                    let mut skipped = Vec::new();
                    let mut valid = true;
                    while cursor != to {
                        let v = loop_[cursor];
                        skipped.push(v);
                        if on_segment(positions[loop_[from]], positions[loop_[to]], positions[v])
                            .is_none()
                        {
                            valid = false;
                        }
                        cursor = (cursor + 1) % loop_.len();
                    }
                    if !valid {
                        revoked.extend(skipped);
                    }
                }
            }
        }
        let previous = removable.len();
        removable.retain(|v| !revoked.contains(v));
        if removable.len() == previous {
            break;
        }
    }
}

/// One globally agreed set of removable seam vertices is applied to *every*
/// face loop. Source-position indices remain stable for evidence and diagnostics.
pub(super) fn simplify_shared_seams(
    patches: &mut [PlanarPatch],
    positions: &[Point3],
    incidence: &EdgeIncidence,
    membership: &[usize],
) -> Result<(), MeshRecoveryError> {
    let mut removable = removable_boundary_vertices(patches, positions, incidence, membership);
    revoke_nonstraight_chains(patches, positions, &mut removable);
    for (patch_id, patch) in patches.iter_mut().enumerate() {
        for loop_ in &mut patch.boundary_loops {
            loop_.retain(|v| !removable.contains(v));
            if loop_.len() < 3 {
                return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
            }
        }
    }
    Ok(())
}
