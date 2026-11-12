//! Assemble unordered interference segments into ordered polylines.
//!
//! Truck-derived marching-SSI support: [`monstertruck_meshing`]'s
//! mesh-mesh interference extraction returns an unordered soup of
//! `(Point3, Point3)` segments; this stitches them into connected
//! [`PolylineCurve`] chains (open or closed) by walking a tolerance-keyed
//! adjacency graph. Resurrected from the 0.3.2 published crate (the last
//! release that shipped a self-contained boolean backend).

use monstertruck_core::{
    cgmath64::{EuclideanSpace, InnerSpace, MetricSpace, Point3},
    tolerance::{DEFAULT_NUMERICAL_TOLERANCES, KERNEL_EPSILON, Tolerance},
};
use monstertruck_meshing::prelude::PolylineCurve;
use num_traits::ToPrimitive;
use rustc_hash::FxHashMap as HashMap;
use std::collections::{BTreeMap, BTreeSet};

/// Decompose an unordered segment set into maximal non-branching chains.
///
/// Branches are geometric junctions, not permission to choose an arbitrary
/// continuation. Every unique undirected edge is emitted exactly once; output
/// order, orientation and node representatives are independent of input order.
pub fn construct_polylines(lines: &[(Point3, Point3)]) -> Vec<PolylineCurve<Point3>> {
    let mut graph: Graph = lines.iter().collect();
    // Preserve the *original* degree classification. Removing the first arm
    // of a junction must not turn the remaining two arms into a false through
    // path, or results would depend on edge-removal order.
    let terminals = graph
        .iter()
        .filter(|(_, node)| node.adjacency.len() != 2)
        .map(|(index, _)| *index)
        .collect::<BTreeSet<_>>();
    let mut polylines = Vec::new();
    for terminal in &terminals {
        while graph.has_edges(*terminal) {
            polylines.push(graph.walk_chain(*terminal, &terminals));
        }
    }
    // The only remaining edges belong to components where every vertex had
    // degree two: closed cycles with a canonical starting point and direction.
    while let Some((start, _)) = graph.first_node() {
        polylines.push(graph.walk_chain(start, &terminals));
    }
    polylines
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

/// Merge a certified, sampled *curved* gap. Unlike an endpoint midpoint
/// snap, this retains the exact interior SSI points supplied by the corrector.
/// The caller must certify the entire bridge on both analytic carriers.
pub(super) fn merge_endpoint_join_with_bridge(
    polylines: &mut Vec<PolylineCurve<Point3>>,
    join: ChainJoin,
    bridge: &[Point3],
    tolerance: f64,
) -> Option<()> {
    if bridge.len() < 2 || join.first.chain == join.second.chain {
        return None;
    }
    let first = polylines.get(join.first.chain)?;
    let second = polylines.get(join.second.chain)?;
    let a = endpoint_point(first, join.first.front)?;
    let b = endpoint_point(second, join.second.front)?;
    if a.distance(*bridge.first()?) > tolerance || b.distance(*bridge.last()?) > tolerance {
        return None;
    }
    let (low, high, low_front, high_front, reverse_bridge) = if join.first.chain < join.second.chain
    {
        (
            join.first.chain,
            join.second.chain,
            join.first.front,
            join.second.front,
            false,
        )
    } else {
        (
            join.second.chain,
            join.first.chain,
            join.second.front,
            join.first.front,
            true,
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
    let mut middle = bridge[1..bridge.len() - 1].to_vec();
    if reverse_bridge {
        middle.reverse();
    }
    left.extend(middle);
    left.extend(right);
    polylines.push(PolylineCurve(left));
    Some(())
}

/// Remove short duplicate strands produced when separate mesh triangles
/// describe the same intersection path with a redundant overlapping arm.
///
/// An endpoint junction alone is not sufficient: every vertex and segment
/// midpoint of the candidate must lie on the longer piecewise-linear path
/// within kernel precision, and the projected positions must progress
/// monotonically along that path. A real branching curve with a distinct
/// locus remains untouched and will be rejected by ambiguous stitching.
pub(super) fn remove_redundant_interference_strands(
    polylines: Vec<PolylineCurve<Point3>>,
    tolerance: f64,
) -> Vec<PolylineCurve<Point3>> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return polylines;
    }
    let lengths = polylines
        .iter()
        .map(|poly| {
            poly.windows(2)
                .map(|segment| segment[0].distance(segment[1]))
                .sum::<f64>()
        })
        .collect::<Vec<_>>();

    let project = |poly: &PolylineCurve<Point3>, point: Point3| {
        let mut cumulative = 0.0;
        let mut closest = None::<(f64, f64)>;
        for segment in poly.windows(2) {
            let direction = segment[1] - segment[0];
            let squared = direction.magnitude2();
            if squared <= 0.0 || !squared.is_finite() {
                continue;
            }
            let length = squared.sqrt();
            let fraction = (point - segment[0]).dot(direction) / squared;
            let fraction = fraction.clamp(0.0, 1.0);
            let foot = segment[0] + direction * fraction;
            let distance2 = foot.distance2(point);
            let progress = cumulative + length * fraction;
            if closest.is_none_or(|(best, _)| distance2 < best) {
                closest = Some((distance2, progress));
            }
            cumulative += length;
        }
        closest
    };
    let contained = |candidate: &PolylineCurve<Point3>, container: &PolylineCurve<Point3>| {
        let Some(start) = candidate.first().copied() else {
            return false;
        };
        let Some(end) = candidate.last().copied() else {
            return false;
        };
        let Some(container_start) = container.first().copied() else {
            return false;
        };
        let Some(container_end) = container.last().copied() else {
            return false;
        };
        // The duplicated arm must begin/end at a shared graph junction.
        if [start, end].into_iter().all(|point| {
            point.distance2(container_start) > tolerance * tolerance
                && point.distance2(container_end) > tolerance * tolerance
        }) {
            return false;
        }
        let mut steps = Vec::new();
        for segment in candidate.windows(2) {
            if steps.is_empty() {
                steps.push(segment[0]);
            }
            steps.push(segment[0].midpoint(segment[1]));
            steps.push(segment[1]);
        }
        let Some(progress) = steps
            .into_iter()
            .map(|point| project(container, point))
            .collect::<Option<Vec<_>>>()
        else {
            return false;
        };
        if progress
            .iter()
            .any(|(distance2, _)| *distance2 > tolerance * tolerance)
        {
            return false;
        }
        let increasing = progress
            .windows(2)
            .all(|pair| pair[1].1 + tolerance >= pair[0].1);
        let decreasing = progress
            .windows(2)
            .all(|pair| pair[0].1 + tolerance >= pair[1].1);
        increasing || decreasing
    };
    let mut keep = vec![true; polylines.len()];
    for (index, candidate) in polylines.iter().enumerate() {
        if !lengths[index].is_finite() || lengths[index] <= 0.0 {
            continue;
        }
        for (other, container) in polylines.iter().enumerate() {
            if index == other || lengths[other] <= lengths[index] + tolerance {
                continue;
            }
            if contained(candidate, container) {
                keep[index] = false;
                break;
            }
        }
    }
    polylines
        .into_iter()
        .zip(keep)
        .filter_map(|(poly, retain)| retain.then_some(poly))
        .collect()
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

#[derive(Clone, Debug, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum CoordinateIndex {
    Bucket(i64),
    Exact(u64),
}

impl CoordinateIndex {
    fn from_coordinate(value: f64) -> Option<Self> {
        if !value.is_finite() {
            return None;
        }
        let grid = DEFAULT_NUMERICAL_TOLERANCES.boolean_polyline_grid_factor * KERNEL_EPSILON;
        let scaled = (value + KERNEL_EPSILON) / grid;
        Some(scaled.to_i64().map_or_else(
            || {
                // Beyond the i64 grid domain, adjacent representable f64 values are
                // already farther apart than the kernel tolerance. Exact bits are
                // therefore a lossless key for this coordinate while the remaining
                // coordinates can still use normal tolerance buckets.
                Self::Exact(value.to_bits())
            },
            Self::Bucket,
        ))
    }
}

#[derive(Clone, Debug, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct PointIndex([CoordinateIndex; 3]);

impl PointIndex {
    fn from_point(pt: Point3) -> Option<Self> {
        Some(Self([
            CoordinateIndex::from_coordinate(pt.x)?,
            CoordinateIndex::from_coordinate(pt.y)?,
            CoordinateIndex::from_coordinate(pt.z)?,
        ]))
    }
}

struct Node {
    coord: Point3,
    // Interference lines describe a geometric set. Adjacent triangle pairs can
    // report the same segment more than once, so parallel edges between the
    // same tolerance-quantized endpoints must collapse here.
    adjacency: BTreeSet<PointIndex>,
}

impl Node {
    const fn new(coord: Point3, adjacency: BTreeSet<PointIndex>) -> Self {
        Self { coord, adjacency }
    }
}

struct Graph(BTreeMap<PointIndex, Node>);

fn point_cmp(left: Point3, right: Point3) -> std::cmp::Ordering {
    left.x
        .total_cmp(&right.x)
        .then_with(|| left.y.total_cmp(&right.y))
        .then_with(|| left.z.total_cmp(&right.z))
}

impl std::ops::Deref for Graph {
    type Target = BTreeMap<PointIndex, Node>;
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl std::ops::DerefMut for Graph {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl Graph {
    fn add_half_edge(&mut self, idx0: PointIndex, pt0: Point3, idx1: PointIndex) {
        if let Some(node) = self.get_mut(&idx0) {
            // Several slightly different endpoints can quantize to one node.
            // Choose its representative deterministically, not "first seen".
            if point_cmp(pt0, node.coord).is_lt() {
                node.coord = pt0;
            }
            node.adjacency.insert(idx1);
        } else {
            self.insert(idx0, Node::new(pt0, BTreeSet::from([idx1])));
        }
    }

    fn add_edge(&mut self, line: (Point3, Point3)) {
        if line.0.near(&line.1) {
            return;
        }
        let (Some(idx0), Some(idx1)) = (
            PointIndex::from_point(line.0),
            PointIndex::from_point(line.1),
        ) else {
            return;
        };
        if idx0 == idx1 {
            // A nonzero segment can still collapse into a single grid node.
            return;
        }
        self.add_half_edge(idx0, line.0, idx1);
        self.add_half_edge(idx1, line.1, idx0);
    }

    fn has_edges(&self, idx: PointIndex) -> bool {
        self.get(&idx)
            .is_some_and(|node| !node.adjacency.is_empty())
    }

    fn first_node(&self) -> Option<(PointIndex, &Node)> {
        self.first_key_value().map(|(idx, node)| (*idx, node))
    }

    fn next_node(&mut self, idx: PointIndex) -> Option<(PointIndex, Point3)> {
        let next = *self.get(&idx)?.adjacency.first()?;
        let node = self.get_mut(&idx)?;
        node.adjacency.remove(&next);
        if node.adjacency.is_empty() {
            self.remove(&idx);
        }
        let node = self.get_mut(&next)?;
        node.adjacency.remove(&idx);
        let point = node.coord;
        if node.adjacency.is_empty() {
            self.remove(&next);
        }
        Some((next, point))
    }

    fn walk_chain(
        &mut self,
        start: PointIndex,
        terminals: &BTreeSet<PointIndex>,
    ) -> PolylineCurve<Point3> {
        let mut points = vec![self[&start].coord];
        let mut current = start;
        while let Some((next, point)) = self.next_node(current) {
            points.push(point);
            if next == start || terminals.contains(&next) {
                break;
            }
            current = next;
        }
        PolylineCurve(points)
    }
}

impl<'a> FromIterator<&'a (Point3, Point3)> for Graph {
    fn from_iter<I: IntoIterator<Item = &'a (Point3, Point3)>>(iter: I) -> Self {
        let mut res = Self(BTreeMap::new());
        iter.into_iter().for_each(|line| res.add_edge(*line));
        res
    }
}

#[cfg(test)]
mod tests;
