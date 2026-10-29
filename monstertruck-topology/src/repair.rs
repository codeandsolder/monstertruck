//! Canonicalization and repair passes for compressed topology.
//!
//! These operations deliberately mutate malformed/intermediate compressed B-reps
//! into topology that can be validated and extracted. Serialization itself lives
//! in [`crate::compress`]; repair policy does not.

use crate::compress::{
    CompressedEdge, CompressedEdgeIndex, CompressedEdgeUse, CompressedFace, CompressedShell,
    CompressedTrimmedFace, CompressedTrimmedShell,
};
use rustc_hash::FxHashMap as HashMap;

trait EdgeUseLike: Clone {
    fn index(&self) -> usize;
    fn orientation(&self) -> bool;
}

impl EdgeUseLike for CompressedEdgeIndex {
    fn index(&self) -> usize { self.index }
    fn orientation(&self) -> bool { self.orientation }
}

impl<T: Clone> EdgeUseLike for CompressedEdgeUse<T> {
    fn index(&self) -> usize { self.index }
    fn orientation(&self) -> bool { self.orientation }
}

fn split_non_simple_boundaries<C, E: EdgeUseLike>(
    edges: &[CompressedEdge<C>],
    boundaries: &mut Vec<Vec<E>>,
) {
    let wire_start_vertices = |wire: &[E]| -> Vec<usize> {
        wire.iter()
            .map(|edge_use| {
                let (v0, v1) = edges[edge_use.index()].vertices;
                if edge_use.orientation() { v0 } else { v1 }
            })
            .collect()
    };
    let mut wire_index = 0usize;
    while wire_index < boundaries.len() {
        let wire = &boundaries[wire_index];
        let starts = wire_start_vertices(wire);
        let mut seen: HashMap<usize, usize> = HashMap::default();
        let revisit = starts
            .iter()
            .enumerate()
            .find_map(|(j, &vertex)| seen.insert(vertex, j).map(|i| (i, j)));
        if let Some((i, j)) = revisit {
            let sub_loop = wire[i..j].to_vec();
            let outer = wire[..i]
                .iter()
                .chain(wire[j..].iter())
                .cloned()
                .collect::<Vec<_>>();
            if sub_loop.len() >= 2 && !outer.is_empty() {
                // Only an exact retrace of the SAME compressed edge is a spike.
                // Distinct geometric/topological edges may connect the same two
                // vertices (a legitimate digon), so endpoint pairs are not an
                // identity key here.
                let mut traversal_balance: HashMap<usize, isize> = HashMap::default();
                sub_loop.iter().for_each(|edge_use| {
                    *traversal_balance.entry(edge_use.index()).or_insert(0) +=
                        if edge_use.orientation() { 1 } else { -1 };
                });
                let is_spike = traversal_balance.values().all(|balance| *balance == 0);
                boundaries[wire_index] = outer;
                if !is_spike {
                    boundaries.push(sub_loop);
                }
                continue;
            }
        }
        wire_index += 1;
    }

    boundaries.retain(|wire| {
        let mut edge_balance: HashMap<usize, isize> = HashMap::default();
        for edge_use in wire {
            *edge_balance.entry(edge_use.index()).or_insert(0) +=
                if edge_use.orientation() { 1 } else { -1 };
        }
        !edge_balance.values().all(|balance| *balance == 0)
    });
}

/// Split face boundary walks that revisit a vertex into simple sub-loops.
///
/// Zero-area retraces whose edge traversals cancel pairwise are discarded. A
/// complete boundary that is only the exact same compressed edge traversed in
/// opposite directions is also discarded; distinct edges joining the same two
/// vertices are preserved because they may form a legitimate digon.
pub fn split_non_simple_compressed_wires<P, C, S>(shell: &mut CompressedShell<P, C, S>) {
    let edges = &shell.edges;
    shell.faces.iter_mut().for_each(|face| {
        split_non_simple_boundaries(edges, &mut face.boundaries);
    });
    shell.faces.retain(|face| !face.boundaries.is_empty());
}

/// Trim-preserving counterpart of [`split_non_simple_compressed_wires`].
///
/// Every [`CompressedEdgeUse`] is moved intact between loops, so its exact
/// face-local trim curve stays attached to the same oriented edge use.
pub fn split_non_simple_compressed_trimmed_wires<P, C, S, T: Clone>(
    shell: &mut CompressedTrimmedShell<P, C, S, T>,
) {
    let edges = &shell.edges;
    shell.faces.iter_mut().for_each(|face| {
        split_non_simple_boundaries(edges, &mut face.boundaries);
    });
    shell.faces.retain(|face| !face.boundaries.is_empty());
}

/// Split singular compressed vertices into one topological vertex per local edge fan.
///
/// This mirrors [`crate::Shell::singular_vertices`]: for each vertex, boundary
/// corners connect the two incident edges in the vertex link. If that link has
/// multiple connected components, the point is duplicated and each extra edge
/// fan is retargeted to its own vertex. Geometry is unchanged; only topological
/// identity is separated. This is what a pinched face needs after a revisiting
/// boundary has been split into simple loops: keep the loops on the same face,
/// but stop making otherwise-disjoint loops share one topology vertex.
///
/// Returns the number of duplicate vertices created.
pub fn split_singular_compressed_trimmed_vertices<P: Clone, C, S, T>(
    shell: &mut CompressedTrimmedShell<P, C, S, T>,
) -> usize {
    use std::collections::{BTreeMap, BTreeSet};

    let original_vertex_count = shell.vertices.len();
    let mut links = (0..original_vertex_count)
        .map(|_| BTreeMap::<usize, BTreeSet<usize>>::new())
        .collect::<Vec<_>>();

    for face in &shell.faces {
        for wire in &face.boundaries {
            if wire.is_empty() {
                continue;
            }
            for index in 0..wire.len() {
                let current = &wire[index];
                let next = &wire[(index + 1) % wire.len()];
                let (cv0, cv1) = shell.edges[current.index].vertices;
                let (nv0, nv1) = shell.edges[next.index].vertices;
                let back = if current.orientation { cv1 } else { cv0 };
                let front = if next.orientation { nv0 } else { nv1 };
                if back != front || back >= original_vertex_count {
                    continue;
                }
                links[back]
                    .entry(current.index)
                    .or_default()
                    .insert(next.index);
                links[back]
                    .entry(next.index)
                    .or_default()
                    .insert(current.index);
            }
        }
    }

    let mut created = 0usize;
    for (vertex, adjacency) in links.iter().enumerate().take(original_vertex_count) {
        if adjacency.len() <= 1 {
            continue;
        }
        let mut unseen = adjacency.keys().copied().collect::<BTreeSet<_>>();
        let mut components = Vec::<Vec<usize>>::new();
        while let Some(seed) = unseen.pop_first() {
            let mut component = vec![seed];
            let mut stack = vec![seed];
            while let Some(edge) = stack.pop() {
                if let Some(neighbors) = adjacency.get(&edge) {
                    for &neighbor in neighbors {
                        if unseen.remove(&neighbor) {
                            component.push(neighbor);
                            stack.push(neighbor);
                        }
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        if components.len() <= 1 {
            continue;
        }
        components.sort_by_key(|component| component[0]);
        for component in components.into_iter().skip(1) {
            let duplicate = shell.vertices.len();
            shell.vertices.push(shell.vertices[vertex].clone());
            created += 1;
            for edge_index in component {
                let edge = &mut shell.edges[edge_index];
                if edge.vertices.0 == vertex {
                    edge.vertices.0 = duplicate;
                }
                if edge.vertices.1 == vertex {
                    edge.vertices.1 = duplicate;
                }
            }
        }
    }
    created
}

/// Diagnostics from a transactional same-domain face-unification pass.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SameDomainUnification {
    /// Candidate adjacent face pairs whose supports were reported equivalent.
    pub proposals: usize,
    /// Face-pair merges that passed trim reconstruction and complete shell validation.
    pub accepted: usize,
    /// Candidates whose surviving directed boundary could not be traced unambiguously.
    pub rejected_ambiguous_boundary: usize,
    /// Candidates for which at least one surviving 3D edge could not be represented
    /// on the chosen reference surface.
    pub rejected_missing_trim: usize,
    /// Candidates whose rebuilt shell failed face, closure, connectivity, or
    /// manifold validation.
    pub rejected_topology: usize,
}

/// Oriented relationship between two parameterized surfaces occupying the same support.
///
/// This is deliberately stronger than a boolean same-support predicate. A face
/// whose carrier parametrization is reversed is mergeable only when its face
/// orientation is reversed as well; otherwise the patches have opposite physical normals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceDomainRelation {
    /// Surface parameterizations induce the same normal orientation.
    SameOrientation,
    /// Surface parameterizations induce opposite normal orientations.
    ReversedOrientation,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SameDomainMergeFailure {
    AmbiguousBoundary,
    MissingTrim,
}

#[inline]
fn oriented_edge_vertices<C, T>(
    edges: &[CompressedEdge<C>],
    edge_use: &CompressedEdgeUse<T>,
) -> (usize, usize) {
    let (front, back) = edges[edge_use.index].vertices;
    if edge_use.orientation {
        (front, back)
    } else {
        (back, front)
    }
}

fn merged_same_domain_face<C, S: Clone, T, F>(
    edges: &[CompressedEdge<C>],
    edge_face_uses: &HashMap<usize, Vec<(usize, bool)>>,
    face0: (usize, &CompressedTrimmedFace<S, T>),
    face1: (usize, &CompressedTrimmedFace<S, T>),
    reverse_face1: bool,
    retrim: &F,
) -> Result<CompressedTrimmedFace<S, T>, SameDomainMergeFailure>
where
    F: Fn(&C, &S) -> Option<T>,
{
    use std::collections::{BTreeMap, BTreeSet};

    let (face_index0, face0) = face0;
    let (face_index1, face1) = face1;
    let uses0 = face0
        .boundaries
        .iter()
        .flatten()
        .map(|edge_use| (edge_use.index, edge_use.orientation))
        .collect::<Vec<_>>();
    let uses1 = face1
        .boundaries
        .iter()
        .flatten()
        .map(|edge_use| (edge_use.index, edge_use.orientation ^ reverse_face1))
        .collect::<Vec<_>>();
    let indices0 = uses0
        .iter()
        .map(|&(index, _)| index)
        .collect::<BTreeSet<_>>();
    let indices1 = uses1
        .iter()
        .map(|&(index, _)| index)
        .collect::<BTreeSet<_>>();
    let shared = indices0
        .intersection(&indices1)
        .copied()
        .collect::<BTreeSet<_>>();
    if shared.is_empty() {
        return Err(SameDomainMergeFailure::AmbiguousBoundary);
    }

    // A removable internal edge must be used exactly once by each candidate
    // face, in opposite directions, and by no third face. Anything else is a
    // seam/non-manifold/keep-edge situation rather than a safe domain partition.
    for &edge_index in &shared {
        let local0 = uses0
            .iter()
            .filter(|&&(index, _)| index == edge_index)
            .collect::<Vec<_>>();
        let local1 = uses1
            .iter()
            .filter(|&&(index, _)| index == edge_index)
            .collect::<Vec<_>>();
        let Some(global) = edge_face_uses.get(&edge_index) else {
            return Err(SameDomainMergeFailure::AmbiguousBoundary);
        };
        if local0.len() != 1
            || local1.len() != 1
            || global.len() != 2
            || global[0].0 == global[1].0
            || local0[0].1 == local1[0].1
            || !global
                .iter()
                .all(|(face, _)| *face == face_index0 || *face == face_index1)
        {
            return Err(SameDomainMergeFailure::AmbiguousBoundary);
        }
    }

    let reference = &face0.surface;
    let mut surviving = Vec::<CompressedEdgeUse<T>>::new();
    for (face, reverse) in [(face0, false), (face1, reverse_face1)] {
        for edge_use in face.boundaries.iter().flatten() {
            if shared.contains(&edge_use.index) {
                continue;
            }
            let trim_curve = retrim(&edges[edge_use.index].curve, reference)
                .ok_or(SameDomainMergeFailure::MissingTrim)?;
            surviving.push(CompressedEdgeUse {
                index: edge_use.index,
                orientation: edge_use.orientation ^ reverse,
                trim_curve: Some(trim_curve),
            });
        }
    }
    if surviving.is_empty() {
        return Err(SameDomainMergeFailure::AmbiguousBoundary);
    }

    // A valid set of disjoint closed loops has one outgoing and one incoming
    // surviving use at every boundary vertex. Requiring that invariant before
    // tracing deliberately rejects pinches, T-junctions, and ambiguous hole
    // ownership instead of guessing a loop decomposition.
    let mut by_start = BTreeMap::<usize, CompressedEdgeUse<T>>::new();
    let mut incoming = BTreeMap::<usize, usize>::new();
    for edge_use in surviving {
        let (start, end) = oriented_edge_vertices(edges, &edge_use);
        if by_start.insert(start, edge_use).is_some() {
            return Err(SameDomainMergeFailure::AmbiguousBoundary);
        }
        *incoming.entry(end).or_default() += 1;
    }
    if incoming.values().any(|count| *count != 1)
        || incoming.len() != by_start.len()
        || incoming.keys().any(|vertex| !by_start.contains_key(vertex))
    {
        return Err(SameDomainMergeFailure::AmbiguousBoundary);
    }

    let mut boundaries = Vec::new();
    while let Some((&start, _)) = by_start.first_key_value() {
        let mut vertex = start;
        let mut wire = Vec::new();
        loop {
            let edge_use = by_start
                .remove(&vertex)
                .ok_or(SameDomainMergeFailure::AmbiguousBoundary)?;
            let (_, end) = oriented_edge_vertices(edges, &edge_use);
            wire.push(edge_use);
            if end == start {
                break;
            }
            vertex = end;
            if wire.len() > edges.len() {
                return Err(SameDomainMergeFailure::AmbiguousBoundary);
            }
        }
        boundaries.push(wire);
    }

    Ok(CompressedTrimmedFace {
        boundaries,
        orientation: face0.orientation,
        surface: face0.surface.clone(),
    })
}

/// Transactionally unify adjacent trimmed faces that occupy the same geometric support.
///
/// `same_domain` answers whether the supports coincide and whether their
/// parameterizations induce the same or reversed normal orientation. For every
/// proposed merge,
/// all surviving face-local trims are supplied on the first face's reference
/// surface through `retrim`; callers may rebuild them from 3-D geometry or
/// transport exact pcurves through a proven chart map. The candidate shell is
/// then fully extracted and checked as
/// a solid boundary before the mutation is committed. Thus a bad equivalence
/// predicate, an ambiguous boundary, or a failed reparameterization can only
/// reject a candidate -- it cannot leave the shell partially modified.
///
/// The input shell itself must already be an extractable closed manifold. This
/// function is a normalizer, not a substitute for structural healing.
pub fn unify_same_domain_compressed_trimmed_faces<P, C, S, T, D, R>(
    shell: &mut CompressedTrimmedShell<P, C, S, T>,
    same_domain: D,
    retrim: R,
) -> crate::Result<SameDomainUnification>
where
    P: Clone,
    C: Clone,
    S: Clone,
    T: Clone,
    D: Fn(&S, &S) -> Option<SurfaceDomainRelation>,
    R: Fn(&C, &S) -> Option<T>,
{
    use crate::trimmed::TrimmedShell;
    use std::collections::BTreeSet;

    let original: TrimmedShell<P, C, S, T> = shell.clone().try_into()?;
    original.check_solid_boundary()?;

    let mut report = SameDomainUnification::default();
    loop {
        let mut edge_face_uses = HashMap::<usize, Vec<(usize, bool)>>::default();
        for (face_index, face) in shell.faces.iter().enumerate() {
            for edge_use in face.boundaries.iter().flatten() {
                edge_face_uses
                    .entry(edge_use.index)
                    .or_default()
                    .push((face_index, edge_use.orientation));
            }
        }
        let mut pairs = BTreeSet::<(usize, usize)>::new();
        for uses in edge_face_uses.values() {
            if uses.len() == 2 && uses[0].0 != uses[1].0 {
                let (a, b) = (uses[0].0.min(uses[1].0), uses[0].0.max(uses[1].0));
                pairs.insert((a, b));
            }
        }

        let mut committed = false;
        for (face_index0, face_index1) in pairs {
            let face0 = &shell.faces[face_index0];
            let face1 = &shell.faces[face_index1];
            let Some(relation) = same_domain(&face0.surface, &face1.surface) else {
                continue;
            };
            let reverse_face1 = matches!(relation, SurfaceDomainRelation::ReversedOrientation);
            let orientation_compatible = if reverse_face1 {
                face0.orientation != face1.orientation
            } else {
                face0.orientation == face1.orientation
            };
            if !orientation_compatible {
                continue;
            }
            report.proposals += 1;
            let merged = match merged_same_domain_face(
                &shell.edges,
                &edge_face_uses,
                (face_index0, face0),
                (face_index1, face1),
                reverse_face1,
                &retrim,
            ) {
                Ok(merged) => merged,
                Err(SameDomainMergeFailure::AmbiguousBoundary) => {
                    report.rejected_ambiguous_boundary += 1;
                    continue;
                }
                Err(SameDomainMergeFailure::MissingTrim) => {
                    report.rejected_missing_trim += 1;
                    continue;
                }
            };

            let mut candidate = shell.clone();
            candidate.faces[face_index0] = merged;
            candidate.faces.remove(face_index1);
            let valid = TrimmedShell::try_from(candidate.clone())
                .and_then(|runtime| runtime.check_solid_boundary().map(|_| runtime));
            if valid.is_err() {
                report.rejected_topology += 1;
                continue;
            }
            *shell = candidate;
            report.accepted += 1;
            committed = true;
            break;
        }
        if !committed {
            break;
        }
    }
    Ok(report)
}

/// Split a compressed face whose individually simple boundary wires share vertices.
///
/// Wires in a vertex-sharing component cannot coexist in one topology face, so
/// each is emitted as a separate face on the same surface. Mutually disjoint
/// boundaries, including ordinary outer-loop-plus-hole faces, stay together.
pub fn split_pinched_compressed_faces<P, C, S: Clone>(shell: &mut CompressedShell<P, C, S>) {
    const fn component_root(components: &mut [usize], index: usize) -> usize {
        let mut root = index;
        while components[root] != root {
            root = components[root];
        }
        let mut cursor = index;
        while components[cursor] != root {
            let next = components[cursor];
            components[cursor] = root;
            cursor = next;
        }
        root
    }

    let faces = std::mem::take(&mut shell.faces);
    let mut new_faces = Vec::with_capacity(faces.len());
    for face in faces {
        if face.boundaries.len() <= 1 {
            new_faces.push(face);
            continue;
        }
        let vertex_sets = face
            .boundaries
            .iter()
            .map(|wire| {
                wire.iter()
                    .flat_map(|edge_use| {
                        let vertices: [usize; 2] = shell.edges[edge_use.index].vertices.into();
                        vertices
                    })
                    .collect::<std::collections::BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let wire_count = face.boundaries.len();
        let mut components = (0..wire_count).collect::<Vec<_>>();
        for i in 0..wire_count {
            for j in (i + 1)..wire_count {
                if !vertex_sets[i].is_disjoint(&vertex_sets[j]) {
                    let (a, b) = (
                        component_root(&mut components, i),
                        component_root(&mut components, j),
                    );
                    if a != b {
                        components[b] = a;
                    }
                }
            }
        }
        let mut component_sizes: HashMap<usize, usize> = HashMap::default();
        for index in 0..wire_count {
            *component_sizes
                .entry(component_root(&mut components, index))
                .or_insert(0) += 1;
        }
        if component_sizes.values().all(|size| *size == 1) {
            new_faces.push(face);
            continue;
        }
        let CompressedFace {
            boundaries,
            orientation,
            surface,
        } = face;
        let mut kept_boundaries = Vec::new();
        for (index, wire) in boundaries.into_iter().enumerate() {
            if component_sizes[&component_root(&mut components, index)] == 1 {
                kept_boundaries.push(wire);
            } else {
                new_faces.push(CompressedFace {
                    boundaries: vec![wire],
                    orientation,
                    surface: surface.clone(),
                });
            }
        }
        if !kept_boundaries.is_empty() {
            new_faces.push(CompressedFace {
                boundaries: kept_boundaries,
                orientation,
                surface,
            });
        }
    }
    shell.faces = new_faces;
}

#[cfg(test)]
mod trimmed_tests {
    use super::*;
    use crate::trimmed::TrimmedShell;
    use std::collections::BTreeSet;

    fn edge(vertices: (usize, usize)) -> CompressedEdge<()> {
        CompressedEdge {
            vertices,
            curve: (),
        }
    }

    fn edge_use(index: usize, trim: u32) -> CompressedEdgeUse<u32> {
        CompressedEdgeUse {
            index,
            orientation: true,
            trim_curve: Some(trim),
        }
    }

    #[test]
    fn trimmed_non_simple_split_keeps_trims_and_disjoins_vertex_fans() {
        // One boundary walks two triangular lobes through vertex 0.
        let mut shell = CompressedTrimmedShell {
            vertices: vec![(); 5],
            edges: vec![
                edge((0, 1)),
                edge((1, 2)),
                edge((2, 0)),
                edge((0, 3)),
                edge((3, 4)),
                edge((4, 0)),
            ],
            faces: vec![CompressedTrimmedFace {
                boundaries: vec![(0..6).map(|i| edge_use(i, 100 + i as u32)).collect()],
                orientation: true,
                surface: (),
            }],
        };

        split_non_simple_compressed_trimmed_wires(&mut shell);
        assert_eq!(shell.faces[0].boundaries.len(), 2);
        assert_eq!(split_singular_compressed_trimmed_vertices(&mut shell), 1);
        assert_eq!(shell.vertices.len(), 6);

        let boundaries = &shell.faces[0].boundaries;
        let vertices = boundaries
            .iter()
            .map(|wire| {
                wire.iter()
                    .flat_map(|edge_use| {
                        let edge = &shell.edges[edge_use.index];
                        [edge.vertices.0, edge.vertices.1]
                    })
                    .collect::<BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        assert!(vertices[0].is_disjoint(&vertices[1]));

        let mut trims = boundaries
            .iter()
            .flatten()
            .map(|edge_use| edge_use.trim_curve.unwrap())
            .collect::<Vec<_>>();
        trims.sort_unstable();
        assert_eq!(trims, vec![100, 101, 102, 103, 104, 105]);
        let _: TrimmedShell<(), (), (), u32> =
            shell.try_into().expect("repaired face must extract");
    }

    #[test]
    fn same_domain_unification_removes_only_the_internal_partition() {
        fn use_(index: usize, orientation: bool) -> CompressedEdgeUse<u32> {
            CompressedEdgeUse {
                index,
                orientation,
                trim_curve: Some(index as u32),
            }
        }
        fn face(surface: u8, uses: &[(usize, bool)]) -> CompressedTrimmedFace<u8, u32> {
            CompressedTrimmedFace {
                boundaries: vec![
                    uses.iter()
                        .map(|&(index, orientation)| use_(index, orientation))
                        .collect(),
                ],
                orientation: true,
                surface,
            }
        }

        // Closed cube topology whose top quad is exported as two coplanar
        // triangles separated by edge 12 (4 -> 6). Every other face uses a
        // distinct support id, so exactly one same-domain merge is available.
        let edge_vertices = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0), // bottom
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4), // top
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7), // verticals
            (4, 6), // top partition
        ];
        let mut shell = CompressedTrimmedShell {
            vertices: vec![(); 8],
            edges: edge_vertices
                .into_iter()
                .enumerate()
                .map(|(index, vertices)| CompressedEdge {
                    vertices,
                    curve: index as u32,
                })
                .collect(),
            faces: vec![
                face(0, &[(3, false), (2, false), (1, false), (0, false)]),
                face(1, &[(0, true), (9, true), (4, false), (8, false)]),
                face(2, &[(1, true), (10, true), (5, false), (9, false)]),
                face(3, &[(2, true), (11, true), (6, false), (10, false)]),
                face(4, &[(3, true), (8, true), (7, false), (11, false)]),
                face(9, &[(4, true), (5, true), (12, false)]),
                face(9, &[(12, true), (6, true), (7, true)]),
            ],
        };

        let before: TrimmedShell<(), u32, u8, u32> =
            shell.clone().try_into().expect("fixture must extract");
        before
            .check_solid_boundary()
            .expect("fixture must be a solid boundary");

        let report = unify_same_domain_compressed_trimmed_faces(
            &mut shell,
            |a, b| (a == b).then_some(SurfaceDomainRelation::SameOrientation),
            |curve, _surface| Some(*curve + 100),
        )
        .expect("normalizer requires and must preserve a valid shell");
        assert_eq!(report.accepted, 1);
        assert_eq!(report.rejected_topology, 0);
        assert_eq!(report.rejected_missing_trim, 0);
        assert_eq!(shell.faces.len(), 6);

        let top = shell.faces.iter().find(|face| face.surface == 9).unwrap();
        assert_eq!(top.boundaries.len(), 1);
        let top_edges = top.boundaries[0]
            .iter()
            .map(|edge_use| edge_use.index)
            .collect::<BTreeSet<_>>();
        assert_eq!(top_edges, BTreeSet::from([4, 5, 6, 7]));
        assert!(!top_edges.contains(&12));
        assert!(
            top.boundaries[0]
                .iter()
                .all(|edge_use| edge_use.trim_curve == Some(edge_use.index as u32 + 100))
        );

        let after: TrimmedShell<(), u32, u8, u32> =
            shell.try_into().expect("normalized shell must extract");
        after
            .check_solid_boundary()
            .expect("normalized shell must stay a solid boundary");
    }
}
