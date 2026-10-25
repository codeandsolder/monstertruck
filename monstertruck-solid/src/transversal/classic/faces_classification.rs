//! Boolean arrangement-region classification.
//!
//! Divider output is allowed to contain non-simple intermediate boundary walks.
//! Those are arrangement regions, not valid B-rep [`Face`]s. Keep them as an
//! explicit private type until ownership is known; only the final selected
//! regions are welded/repaired into valid topology.

use super::loops_store::ShapesOpStatus;
use monstertruck_topology::{compress::*, repair, *};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

#[derive(Clone, Debug)]
pub(super) struct BooleanRegion<P, C, S> {
    boundaries: Vec<Wire<P, C>>,
    surface: S,
    orientation: bool,
    key: usize,
    source_face: usize,
}

impl<P, C, S> BooleanRegion<P, C, S> {
    pub(super) fn new(
        boundaries: Vec<Wire<P, C>>,
        surface: S,
        orientation: bool,
        key: usize,
        source_face: usize,
    ) -> Self {
        Self {
            boundaries,
            surface,
            orientation,
            key,
            source_face,
        }
    }

    pub(super) const fn key(&self) -> usize { self.key }
    pub(super) const fn source_face(&self) -> usize { self.source_face }
    pub(super) const fn absolute_boundaries(&self) -> &Vec<Wire<P, C>> { &self.boundaries }
    pub(super) const fn orientation(&self) -> bool { self.orientation }
    pub(super) fn invert(&mut self) { self.orientation = !self.orientation; }

    pub(super) fn surface(&self) -> S
    where S: Clone {
        self.surface.clone()
    }

    pub(super) fn try_face(&self) -> Option<Face<P, C, S>>
    where
        P: Clone,
        C: Clone,
        S: Clone, {
        let mut face = Face::try_new(self.boundaries.clone(), self.surface.clone()).ok()?;
        if !self.orientation {
            face.invert();
        }
        Some(face)
    }
}

fn cancel_retraced_edge_uses<P: Clone, C: Clone>(
    boundaries: &[Wire<P, C>],
) -> Option<(Vec<Wire<P, C>>, bool)> {
    let uses = boundaries
        .iter()
        .flat_map(|wire| wire.iter().cloned())
        .collect::<Vec<_>>();
    let mut by_edge = HashMap::<EdgeId<C>, (Vec<usize>, Vec<usize>)>::default();
    for (index, edge) in uses.iter().enumerate() {
        let entry = by_edge.entry(edge.id()).or_default();
        if edge.orientation() {
            entry.0.push(index);
        } else {
            entry.1.push(index);
        }
    }
    let mut cancelled = vec![false; uses.len()];
    for (forward, reverse) in by_edge.values() {
        for (&lhs, &rhs) in forward.iter().zip(reverse) {
            cancelled[lhs] = true;
            cancelled[rhs] = true;
        }
    }
    if !cancelled.iter().any(|cancelled| *cancelled) {
        return Some((boundaries.to_vec(), false));
    }
    let mut remaining = uses
        .into_iter()
        .enumerate()
        .filter_map(|(index, edge)| (!cancelled[index]).then_some(edge))
        .collect::<Vec<_>>();
    let mut wires = Vec::new();
    while !remaining.is_empty() {
        let first = remaining.remove(0);
        let start = first.front().id();
        let mut chain = vec![first];
        loop {
            let back = chain.last()?.back().id();
            if back == start {
                break;
            }
            let mut next = None;
            for (index, edge) in remaining.iter().enumerate() {
                if edge.front().id() != back {
                    continue;
                }
                if next.is_some() {
                    return None;
                }
                next = Some(index);
            }
            chain.push(remaining.remove(next?));
        }
        wires.push(Wire::from(chain));
    }
    Some((wires, true))
}

fn materialize_region_set<P: Clone, C: Clone, S: Clone>(
    regions: &[BooleanRegion<P, C, S>],
) -> Option<Shell<P, C, S>> {
    let direct = regions
        .iter()
        .map(BooleanRegion::try_face)
        .collect::<Option<Vec<_>>>();
    if let Some(faces) = direct {
        return Some(Shell::from(faces));
    }

    let mut vertices = Vec::<P>::new();
    let mut vertex_map = HashMap::<VertexId<P>, usize>::default();
    let mut edges = Vec::<CompressedEdge<C>>::new();
    let mut edge_objects = Vec::<Edge<P, C>>::new();
    let mut edge_map = HashMap::<EdgeId<C>, usize>::default();
    let mut faces = Vec::<CompressedFace<S>>::with_capacity(regions.len());

    for region in regions {
        let mut vertex_index = |vertex: &Vertex<P>| {
            if let Some(index) = vertex_map.get(&vertex.id()).copied() {
                index
            } else {
                let index = vertices.len();
                vertices.push(vertex.point());
                vertex_map.insert(vertex.id(), index);
                index
            }
        };
        let boundaries = region
            .absolute_boundaries()
            .iter()
            .map(|wire| {
                wire.iter()
                    .map(|edge| {
                        let index = if let Some(index) = edge_map.get(&edge.id()).copied() {
                            index
                        } else {
                            let index = edges.len();
                            edges.push(CompressedEdge {
                                vertices: (
                                    vertex_index(edge.absolute_front()),
                                    vertex_index(edge.absolute_back()),
                                ),
                                curve: edge.curve(),
                            });
                            edge_objects.push(edge.absolute_clone());
                            edge_map.insert(edge.id(), index);
                            index
                        };
                        CompressedEdgeIndex {
                            index,
                            orientation: edge.orientation(),
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        faces.push(CompressedFace {
            boundaries,
            orientation: region.orientation(),
            surface: region.surface(),
        });
    }

    let mut compressed = CompressedShell {
        vertices,
        edges,
        faces,
        vertex_stable_ids: None,
        edge_stable_ids: None,
        face_stable_ids: None,
    };
    repair::split_non_simple_compressed_wires(&mut compressed);
    repair::split_pinched_compressed_faces(&mut compressed);

    let faces = compressed
        .faces
        .into_iter()
        .map(|compressed_face| {
            let boundaries = compressed_face
                .boundaries
                .into_iter()
                .map(|wire| {
                    wire.into_iter()
                        .map(|edge_use| {
                            let edge = edge_objects.get(edge_use.index)?.clone();
                            Some(if edge_use.orientation {
                                edge
                            } else {
                                edge.inverse()
                            })
                        })
                        .collect::<Option<Vec<_>>>()
                        .map(Wire::from)
                })
                .collect::<Option<Vec<_>>>()?;
            let mut face = Face::try_new(boundaries, compressed_face.surface).ok()?;
            if !compressed_face.orientation {
                face.invert();
            }
            Some(face)
        })
        .collect::<Option<Vec<_>>>()?;
    Some(Shell::from(faces))
}

pub(super) fn materialize_regions<P: Clone, C: Clone, S: Clone>(
    regions: &[BooleanRegion<P, C, S>],
) -> Option<Shell<P, C, S>> {
    let mut changed = false;
    let mut normalized = Vec::with_capacity(regions.len());
    for region in regions {
        let Some((boundaries, did_change)) =
            cancel_retraced_edge_uses(region.absolute_boundaries())
        else {
            // A real retrace can still leave several valid ways to stitch the
            // surviving directed edges. Do not guess locally: keep the original
            // arrangement and let the final global repair see all selected regions.
            return materialize_region_set(regions);
        };
        changed |= did_change;
        if boundaries.is_empty() {
            continue;
        }
        if did_change {
            normalized.push(BooleanRegion::new(
                boundaries,
                region.surface(),
                region.orientation(),
                region.key(),
                region.source_face(),
            ));
        } else {
            normalized.push(region.clone());
        }
    }
    if !changed {
        return materialize_region_set(regions);
    }
    materialize_region_set(&normalized)
}

#[cfg(test)]
mod materialization_tests {
    use super::*;

    #[test]
    fn distinct_edge_digon_is_not_treated_as_retrace() {
        let front = Vertex::new(());
        let back = Vertex::new(());
        let edge0 = Edge::new(&front, &back, ());
        let edge1 = Edge::new(&back, &front, ());
        let region = BooleanRegion::new(vec![Wire::from(vec![edge0, edge1])], (), true, 0, 0);

        let materialized = materialize_regions(&[region]);
        assert!(materialized.is_some());
        let shell = materialized.unwrap_or_default();
        assert_eq!(shell.len(), 1);
        assert_eq!(shell[0].absolute_boundaries()[0].len(), 2);
    }
}

#[derive(Clone, Debug)]
pub(super) struct FacesClassification<P, C, S> {
    regions: Vec<BooleanRegion<P, C, S>>,
    status: Vec<ShapesOpStatus>,
}

impl<P, C, S> Default for FacesClassification<P, C, S> {
    fn default() -> Self {
        Self {
            regions: Vec::new(),
            status: Vec::new(),
        }
    }
}

impl<P, C, S> FacesClassification<P, C, S> {
    pub(super) fn push_region(
        &mut self,
        boundaries: Vec<Wire<P, C>>,
        surface: S,
        orientation: bool,
        status: ShapesOpStatus,
        source_face: usize,
    ) {
        let key = self.regions.len();
        self.regions.push(BooleanRegion::new(
            boundaries,
            surface,
            orientation,
            key,
            source_face,
        ));
        self.status.push(status);
    }

    pub(super) fn push(&mut self, face: Face<P, C, S>, status: ShapesOpStatus, source_face: usize)
    where
        P: Clone,
        C: Clone,
        S: Clone, {
        self.push_region(
            face.absolute_boundaries().clone(),
            face.surface(),
            face.orientation(),
            status,
            source_face,
        );
    }

    pub(super) fn and_or_unknown(&self) -> [Vec<BooleanRegion<P, C, S>>; 3]
    where
        P: Clone,
        C: Clone,
        S: Clone, {
        let [mut and, mut or, mut unknown] = [Vec::new(), Vec::new(), Vec::new()];
        for (region, status) in self.regions.iter().zip(&self.status) {
            match status {
                ShapesOpStatus::And => and.push(region.clone()),
                ShapesOpStatus::Or => or.push(region.clone()),
                ShapesOpStatus::Unknown => unknown.push(region.clone()),
            }
        }
        [and, or, unknown]
    }

    pub(super) fn all_regions(&self) -> Vec<BooleanRegion<P, C, S>>
    where
        P: Clone,
        C: Clone,
        S: Clone, {
        self.regions.clone()
    }

    pub(super) fn integrate_by_component(&mut self) {
        let boundary_edges_for = |wanted: ShapesOpStatus| {
            let mut counts = HashMap::<EdgeId<C>, usize>::default();
            for (region, status) in self.regions.iter().zip(&self.status) {
                if *status != wanted {
                    continue;
                }
                for edge in region.boundaries.iter().flat_map(|wire| wire.iter()) {
                    *counts.entry(edge.id()).or_insert(0) += 1;
                }
            }
            counts
                .into_iter()
                .filter_map(|(edge, count)| (count == 1).then_some(edge))
                .collect::<HashSet<_>>()
        };
        let and_boundary = boundary_edges_for(ShapesOpStatus::And);
        let or_boundary = boundary_edges_for(ShapesOpStatus::Or);

        let unknown_indices = self
            .status
            .iter()
            .enumerate()
            .filter_map(|(index, status)| (*status == ShapesOpStatus::Unknown).then_some(index))
            .collect::<Vec<_>>();
        let mut edge_to_unknown = HashMap::<EdgeId<C>, Vec<usize>>::default();
        for &index in &unknown_indices {
            for edge in self.regions[index]
                .boundaries
                .iter()
                .flat_map(|wire| wire.iter())
            {
                edge_to_unknown.entry(edge.id()).or_default().push(index);
            }
        }

        let mut unseen = unknown_indices.into_iter().collect::<HashSet<_>>();
        while let Some(&start) = unseen.iter().next() {
            let mut stack = vec![start];
            let mut component = Vec::new();
            unseen.remove(&start);
            while let Some(index) = stack.pop() {
                component.push(index);
                for edge in self.regions[index]
                    .boundaries
                    .iter()
                    .flat_map(|wire| wire.iter())
                {
                    if let Some(neighbors) = edge_to_unknown.get(&edge.id()) {
                        for &neighbor in neighbors {
                            if unseen.remove(&neighbor) {
                                stack.push(neighbor);
                            }
                        }
                    }
                }
            }

            let mut component_edges = HashMap::<EdgeId<C>, usize>::default();
            for &index in &component {
                for edge in self.regions[index]
                    .boundaries
                    .iter()
                    .flat_map(|wire| wire.iter())
                {
                    *component_edges.entry(edge.id()).or_insert(0) += 1;
                }
            }
            let component_boundary = component_edges
                .into_iter()
                .filter_map(|(edge, count)| (count == 1).then_some(edge))
                .collect::<HashSet<_>>();
            let touches_and = component_boundary
                .iter()
                .any(|edge| and_boundary.contains(edge));
            let touches_or = component_boundary
                .iter()
                .any(|edge| or_boundary.contains(edge));
            let status = match (touches_and, touches_or) {
                (true, false) => ShapesOpStatus::And,
                (false, true) => ShapesOpStatus::Or,
                _ => continue,
            };
            for index in component {
                if let Some(current) = self.status.get_mut(index) {
                    *current = status;
                }
            }
        }
    }
}
