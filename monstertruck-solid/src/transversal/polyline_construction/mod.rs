//! Assemble unordered interference segments into ordered polylines.
//!
//! Truck-derived marching-SSI support: [`monstertruck_meshing`]'s
//! mesh-mesh interference extraction returns an unordered soup of
//! `(Point3, Point3)` segments; this stitches them into connected
//! [`PolylineCurve`] chains (open or closed) by walking a tolerance-keyed
//! adjacency graph. Resurrected from the 0.3.2 published crate (the last
//! release that shipped a self-contained boolean backend).

use monstertruck_core::{cgmath64::*, tolerance::*};
use monstertruck_meshing::prelude::PolylineCurve;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::collections::VecDeque;

pub fn construct_polylines(lines: &[(Point3, Point3)]) -> Vec<PolylineCurve<Point3>> {
    let mut graph: Graph = lines.iter().collect();
    let mut res = Vec::new();
    while !graph.is_empty() {
        let (mut idx, node) = graph.first_node();
        let mut wire: VecDeque<_> = vec![node.coord].into();
        while let Some((idx0, pt)) = graph.next_node(idx) {
            idx = idx0;
            wire.push_back(pt);
        }
        let mut idx = PointIndex::from(wire[0]);
        while let Some((idx0, pt)) = graph.next_node(idx) {
            idx = idx0;
            wire.push_front(pt);
        }
        res.push(PolylineCurve(wire.into()));
    }
    res
}

/// Merge separate interference chains whose endpoints agree within the caller's
/// geometric tolerance. Only endpoints of distinct chains are considered, and
/// every endpoint must have at most one possible mate. Ambiguous junctions are
/// rejected rather than guessed.
pub(super) fn stitch_nearby_polylines(
    mut polylines: Vec<PolylineCurve<Point3>>,
    tolerance: f64,
) -> Option<Vec<PolylineCurve<Point3>>> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let tolerance2 = tolerance * tolerance;

    loop {
        #[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
        struct Endpoint {
            chain: usize,
            front: bool,
        }
        #[derive(Clone, Copy, Debug)]
        struct Candidate {
            a: Endpoint,
            b: Endpoint,
            distance2: f64,
        }

        let mut candidates = Vec::<Candidate>::new();
        for first in 0..polylines.len() {
            let first_points = &polylines[first].0;
            if first_points.is_empty() {
                return None;
            }
            for (second, second_chain) in polylines.iter().enumerate().skip(first + 1) {
                let second_points = &second_chain.0;
                if second_points.is_empty() {
                    return None;
                }
                for first_front in [true, false] {
                    let a = if first_front {
                        first_points[0]
                    } else {
                        *first_points.last()?
                    };
                    for second_front in [true, false] {
                        let b = if second_front {
                            second_points[0]
                        } else {
                            *second_points.last()?
                        };
                        let distance2 = a.distance2(b);
                        if distance2 <= tolerance2 {
                            candidates.push(Candidate {
                                a: Endpoint {
                                    chain: first,
                                    front: first_front,
                                },
                                b: Endpoint {
                                    chain: second,
                                    front: second_front,
                                },
                                distance2,
                            });
                        }
                    }
                }
            }
        }
        if candidates.is_empty() {
            return Some(polylines);
        }

        let mut endpoint_counts = HashMap::<Endpoint, usize>::default();
        for candidate in &candidates {
            *endpoint_counts.entry(candidate.a).or_default() += 1;
            *endpoint_counts.entry(candidate.b).or_default() += 1;
        }
        if endpoint_counts.values().any(|&count| count != 1) {
            return None;
        }

        candidates.sort_by(|left, right| {
            left.distance2
                .total_cmp(&right.distance2)
                .then_with(|| left.a.cmp(&right.a))
                .then_with(|| left.b.cmp(&right.b))
        });
        let candidate = candidates[0];
        let (low, high, low_front, high_front) = if candidate.a.chain < candidate.b.chain {
            (
                candidate.a.chain,
                candidate.b.chain,
                candidate.a.front,
                candidate.b.front,
            )
        } else {
            (
                candidate.b.chain,
                candidate.a.chain,
                candidate.b.front,
                candidate.a.front,
            )
        };

        let mut right = polylines.remove(high).0;
        let mut left = polylines.remove(low).0;
        if low_front {
            left.reverse();
        }
        if !high_front {
            right.reverse();
        }
        let left_end = *left.last()?;
        let right_start = *right.first()?;
        let midpoint = left_end.midpoint(right_start);
        *left.last_mut()? = midpoint;
        *right.first_mut()? = midpoint;
        left.extend(right.into_iter().skip(1));
        polylines.push(PolylineCurve(left));
    }
}

#[derive(Clone, Debug, Copy, Hash, PartialEq, Eq)]
struct PointIndex([i64; 3]);

impl From<Point3> for PointIndex {
    #[inline(always)]
    fn from(pt: Point3) -> PointIndex {
        let idx = pt.add_element_wise(TOLERANCE) / (2.0 * TOLERANCE);
        // SAFETY: point coordinates are finite, so the cast to `i64` always succeeds.
        PointIndex(idx.cast::<i64>().unwrap().into())
    }
}

struct Node {
    coord: Point3,
    // Interference lines describe a geometric set. Adjacent triangle pairs can
    // report the same segment more than once, so parallel edges between the
    // same tolerance-quantized endpoints must collapse here.
    adjacency: HashSet<PointIndex>,
}

impl Node {
    #[inline(always)]
    fn new(coord: Point3, adjacency: HashSet<PointIndex>) -> Node {
        Node { coord, adjacency }
    }

    fn pop_one_adjacency(&mut self) -> PointIndex {
        // SAFETY: nodes are removed from the graph when their adjacency set becomes empty.
        let idx = *self.adjacency.iter().next().unwrap();
        self.adjacency.remove(&idx);
        idx
    }
}

struct Graph(HashMap<PointIndex, Node>);

impl std::ops::Deref for Graph {
    type Target = HashMap<PointIndex, Node>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for Graph {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Graph {
    fn add_half_edge(&mut self, pt0: Point3, pt1: Point3) {
        let idx0 = pt0.into();
        let idx1 = pt1.into();
        if let Some(node) = self.get_mut(&idx0) {
            node.adjacency.insert(idx1);
        } else {
            self.insert(idx0, Node::new(pt0, HashSet::from_iter([idx1])));
        }
    }

    fn add_edge(&mut self, line: (Point3, Point3)) {
        if !line.0.near(&line.1) {
            self.add_half_edge(line.0, line.1);
            self.add_half_edge(line.1, line.0);
        }
    }

    #[inline(always)]
    fn first_node(&self) -> (PointIndex, &Node) {
        // SAFETY: only called inside `while !graph.is_empty()`.
        let (idx, node) = self.iter().next().unwrap();
        (*idx, node)
    }

    fn next_node(&mut self, idx: PointIndex) -> Option<(PointIndex, Point3)> {
        let node = self.get_mut(&idx)?;
        let idx0 = node.pop_one_adjacency();
        if node.adjacency.is_empty() {
            self.remove(&idx);
        }
        let node = self.get_mut(&idx0)?;
        node.adjacency.remove(&idx);
        let pt = node.coord;
        if node.adjacency.is_empty() {
            self.remove(&idx0);
        }
        Some((idx0, pt))
    }
}

impl<'a> FromIterator<&'a (Point3, Point3)> for Graph {
    fn from_iter<I: IntoIterator<Item = &'a (Point3, Point3)>>(iter: I) -> Graph {
        let mut res = Graph(HashMap::default());
        iter.into_iter().for_each(|line| res.add_edge(*line));
        res
    }
}

#[cfg(test)]
mod tests;
