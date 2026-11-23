//! Topologically global edge repair and simplification for inverse meshing.
//!
//! Split source edges at existing geometric T-junction vertices before building
//! edge incidence. After patching, remove straight-through vertices only if
//! their two incident *unique* boundary edges have exactly the same pair of
//! adjacent patches; this guarantees both sides of a seam receive the same cut.

use super::{EdgeIncidence, EdgeKey, InputPolygon, MeshRecoveryError, PlanarPatch, edge_key};
use crate::{InnerSpace, Point3};
use rustc_hash::FxHashMap as HashMap;
use std::collections::HashSet;

/// The three projected 2D orientation tests certify collinearity in 3D
/// using adaptive exact predicates for the *stored f64 input coordinates*.
/// A fitted near-line is intentionally not a shared topological edge.
fn on_line(a: Point3, b: Point3, p: Point3) -> bool {
    use robust::{Coord, orient2d};
    let xy = |v: Point3| Coord { x: v.x, y: v.y };
    let xz = |v: Point3| Coord { x: v.x, y: v.z };
    let yz = |v: Point3| Coord { x: v.y, y: v.z };
    orient2d(xy(a), xy(b), xy(p)) == 0.0
        && orient2d(xz(a), xz(b), xz(p)) == 0.0
        && orient2d(yz(a), yz(b), yz(p)) == 0.0
}

/// Returns the coordinate along the directed segment only for a strictly
/// interior vertex exactly collinear with the source f64 endpoints.
fn on_segment(a: Point3, b: Point3, p: Point3) -> Option<f64> {
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
/// Sorting by the longest coordinate axis allows range queries of candidate
/// points without scanning all mesh positions for each source edge.
pub(super) fn split_source_t_junctions(positions: &[Point3], polygons: &mut [InputPolygon]) {
    let mut used = vec![false; positions.len()];
    let mut spokes = vec![Vec::new(); positions.len()];
    for polygon in polygons.iter() {
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

    for polygon in polygons.iter_mut() {
        let mut split = Vec::with_capacity(polygon.vertices.len());
        for i in 0..polygon.vertices.len() {
            let a = polygon.vertices[i];
            let b = polygon.vertices[(i + 1) % polygon.vertices.len()];
            split.push(a);

            let vector = positions[b] - positions[a];
            let axis = (0..3)
                .max_by(|&x, &y| vector[x].abs().total_cmp(&vector[y].abs()))
                .unwrap_or(0);
            let first = positions[a][axis].min(positions[b][axis]);
            let last = positions[a][axis].max(positions[b][axis]);
            let candidates = &sorted[axis];
            let begin = candidates.partition_point(|&p| positions[p][axis] <= first);
            let end = candidates.partition_point(|&p| positions[p][axis] < last);
            let mut inner = candidates[begin..end]
                .iter()
                .filter(|&&id| {
                    spokes[id]
                        .iter()
                        .any(|&neighbor| on_line(positions[a], positions[b], positions[neighbor]))
                })
                .filter_map(|&id| {
                    on_segment(positions[a], positions[b], positions[id])
                        .map(|position| (position, id))
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
