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

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct ChainEndpoint {
    pub(super) chain: usize,
    pub(super) front: bool,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct ChainJoin {
    pub(super) first: ChainEndpoint,
    pub(super) second: ChainEndpoint,
    pub(super) distance2: f64,
}

pub(super) fn endpoint_point(polyline: &PolylineCurve<Point3>, front: bool) -> Option<Point3> {
    if front {
        polyline.first().copied()
    } else {
        polyline.last().copied()
    }
}

/// Enumerate endpoint pairs accepted by a caller-supplied geometric validation.
///
/// This owns only chain topology: endpoint orientation, deterministic ordering,
/// and distance metadata. Whether a pair is actually joinable remains the
/// caller's policy (simple proximity, SSI carrier validation, etc.).
pub(super) fn enumerate_endpoint_joins(
    polylines: &[PolylineCurve<Point3>],
    mut accepts: impl FnMut(&PolylineCurve<Point3>, bool, &PolylineCurve<Point3>, bool) -> bool,
) -> Option<Vec<ChainJoin>> {
    let mut joins = Vec::new();
    for first in 0..polylines.len() {
        if polylines[first].is_empty() {
            return None;
        }
        for second in first + 1..polylines.len() {
            if polylines[second].is_empty() {
                return None;
            }
            for first_front in [true, false] {
                let a = endpoint_point(&polylines[first], first_front)?;
                for second_front in [true, false] {
                    if !accepts(
                        &polylines[first],
                        first_front,
                        &polylines[second],
                        second_front,
                    ) {
                        continue;
                    }
                    let b = endpoint_point(&polylines[second], second_front)?;
                    joins.push(ChainJoin {
                        first: ChainEndpoint {
                            chain: first,
                            front: first_front,
                        },
                        second: ChainEndpoint {
                            chain: second,
                            front: second_front,
                        },
                        distance2: a.distance2(b),
                    });
                }
            }
        }
    }
    Some(joins)
}

/// Apply one selected endpoint join using the canonical orientation/snap rule.
pub(super) fn merge_endpoint_join(
    polylines: &mut Vec<PolylineCurve<Point3>>,
    join: ChainJoin,
) -> Option<()> {
    let (low, high, low_front, high_front) = if join.first.chain < join.second.chain {
        (
            join.first.chain,
            join.second.chain,
            join.first.front,
            join.second.front,
        )
    } else {
        (
            join.second.chain,
            join.first.chain,
            join.second.front,
            join.first.front,
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
    Some(())
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
        let mut candidates =
            enumerate_endpoint_joins(&polylines, |first, first_front, second, second_front| {
                let Some(a) = endpoint_point(first, first_front) else {
                    return false;
                };
                let Some(b) = endpoint_point(second, second_front) else {
                    return false;
                };
                a.distance2(b) <= tolerance2
            })?;
        if candidates.is_empty() {
            return Some(polylines);
        }

        let mut endpoint_counts = HashMap::<ChainEndpoint, usize>::default();
        for candidate in &candidates {
            *endpoint_counts.entry(candidate.first).or_default() += 1;
            *endpoint_counts.entry(candidate.second).or_default() += 1;
        }
        if endpoint_counts.values().any(|&count| count != 1) {
            return None;
        }

        candidates.sort_by(|left, right| {
            left.distance2
                .total_cmp(&right.distance2)
                .then_with(|| left.first.cmp(&right.first))
                .then_with(|| left.second.cmp(&right.second))
        });
        merge_endpoint_join(&mut polylines, candidates[0])?;
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
    fn new(coord: Point3, adjacency: HashSet<PointIndex>) -> Node { Node { coord, adjacency } }

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
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl std::ops::DerefMut for Graph {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
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
