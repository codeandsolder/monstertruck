//! Canonicalization and repair passes for compressed topology.
//!
//! These operations deliberately mutate malformed/intermediate compressed B-reps
//! into topology that can be validated and extracted. Serialization itself lives
//! in [`crate::compress`]; repair policy does not.

use crate::compress::{CompressedEdgeIndex, CompressedFace, CompressedShell};
use rustc_hash::FxHashMap as HashMap;

/// Split face boundary walks that revisit a vertex into simple sub-loops.
///
/// Zero-area retraces whose edge traversals cancel pairwise are discarded. A
/// complete boundary that is only the exact same compressed edge traversed in
/// opposite directions is also discarded; distinct edges joining the same two
/// vertices are preserved because they may form a legitimate digon.
pub fn split_non_simple_compressed_wires<P, C, S>(shell: &mut CompressedShell<P, C, S>) {
    let edges = &shell.edges;
    let wire_start_vertices = |wire: &[CompressedEdgeIndex]| -> Vec<usize> {
        wire.iter()
            .map(|edge_use| {
                let (v0, v1) = edges[edge_use.index].vertices;
                if edge_use.orientation { v0 } else { v1 }
            })
            .collect()
    };
    shell.faces.iter_mut().for_each(|face| {
        let mut wire_index = 0usize;
        while wire_index < face.boundaries.len() {
            let wire = &face.boundaries[wire_index];
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
                    .copied()
                    .collect::<Vec<_>>();
                if sub_loop.len() >= 2 && !outer.is_empty() {
                    let mut traversal_balance: HashMap<(usize, usize), isize> = HashMap::default();
                    for edge_use in &sub_loop {
                        let (v0, v1) = edges[edge_use.index].vertices;
                        let (start, end) = if edge_use.orientation {
                            (v0, v1)
                        } else {
                            (v1, v0)
                        };
                        let key = (start.min(end), start.max(end));
                        *traversal_balance.entry(key).or_insert(0) +=
                            if start <= end { 1 } else { -1 };
                    }
                    let is_spike = traversal_balance.values().all(|balance| *balance == 0);
                    face.boundaries[wire_index] = outer;
                    if !is_spike {
                        face.boundaries.push(sub_loop);
                    }
                    continue;
                }
            }
            wire_index += 1;
        }

        face.boundaries.retain(|wire| {
            let mut edge_balance: HashMap<usize, isize> = HashMap::default();
            for edge_use in wire {
                *edge_balance.entry(edge_use.index).or_insert(0) +=
                    if edge_use.orientation { 1 } else { -1 };
            }
            !edge_balance.values().all(|balance| *balance == 0)
        });
    });
    shell.faces.retain(|face| !face.boundaries.is_empty());
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
