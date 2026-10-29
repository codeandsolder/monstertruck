//! Intersection-loop store for the classic Boolean backend.
//!
//! This began as the published 0.3.2 loop store, but now also owns exact-boundary
//! endpoint resolution, coplanar interface matching, and partition stitching.
//! Geometric equality checks here are deliberately separated from broad-phase
//! signatures: sampled or summary evidence may select candidates, never prove
//! trimmed-region identity by itself.

#![allow(
    clippy::many_single_char_names,
    reason = "loop-intersection formulas follow compact geometric notation"
)]

use super::{
    curve_equivalence::{
        CurveCoincidence, validate_endpoint_aligned_coincidence, validate_sampled_coincidence,
        validate_sampled_containment,
    },
    intersection_curve,
    tolerance::BooleanTolerance,
};
mod planar;

use monstertruck_core::{cgmath64::*, tolerance::TOLERANCE};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::{Vertex, *};
use planar::{
    CROSS_FACE_STITCH_TOLERANCE, FullPlanarInterfaceMatch, imprint_edges_on_faces,
    planar_polygons_same_trimmed_face, shells_full_planar_interface_pair,
    shells_have_partial_planar_overlap, shells_have_positive_aabb_overlap,
};
pub(in crate::transversal::classic) use planar::{
    direct_planar_face_match, single_face_matches_planar_partition,
};
use rustc_hash::FxHashMap as HashMap;

type PolylineCurve = monstertruck_meshing::prelude::PolylineCurve<Point3>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ShapesOpStatus {
    Unknown,
    And,
    Or,
}

impl ShapesOpStatus {
    fn not(self) -> Self {
        match self {
            Self::Unknown => Self::Unknown,
            Self::And => Self::Or,
            Self::Or => Self::And,
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct BoundaryWire<P, C> {
    wire: Wire<P, C>,
    status: ShapesOpStatus,
}

impl<P, C> BoundaryWire<P, C> {
    #[inline(always)]
    pub(super) fn new(wire: Wire<P, C>, status: ShapesOpStatus) -> Self { Self { wire, status } }
    #[inline(always)]
    pub(super) fn status(&self) -> ShapesOpStatus { self.status }
    #[inline(always)]
    pub(super) fn inverse(&self) -> Self {
        Self {
            wire: self.wire.inverse(),
            status: self.status.not(),
        }
    }
}

impl ShapesOpStatus {
    fn from_is_curve<C, S0, S1>(curve: &IntersectionCurve<C, S0, S1>) -> Option<ShapesOpStatus>
    where
        C: ParametricCurve3D + BoundedCurve,
        S0: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
        S1: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>, {
        let (t0, t1) = curve.range_tuple();
        let t = (t0 + t1) / 2.0;
        let (_, pt0, pt1) = curve.search_triple(t, 100)?;
        let der = curve.leader().derivative(t);
        let normal0 = curve.surface0().normal(pt0[0], pt0[1]);
        let normal1 = curve.surface1().normal(pt1[0], pt1[1]);
        match normal0.cross(der).dot(normal1) > 0.0 {
            true => Some(ShapesOpStatus::Or),
            false => Some(ShapesOpStatus::And),
        }
    }
}

impl<P, C> std::ops::Deref for BoundaryWire<P, C> {
    type Target = Wire<P, C>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target { &self.wire }
}

impl<P, C> std::ops::DerefMut for BoundaryWire<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.wire }
}

#[derive(Clone, Debug)]
pub(super) struct Loops<P, C>(Vec<BoundaryWire<P, C>>);
#[derive(Clone, Debug)]
pub(super) struct LoopsStore<P, C>(Vec<Loops<P, C>>);

impl<P, C> std::ops::Deref for Loops<P, C> {
    type Target = Vec<BoundaryWire<P, C>>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P, C> std::ops::DerefMut for Loops<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P, C> std::ops::Deref for LoopsStore<P, C> {
    type Target = Vec<Loops<P, C>>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P, C> std::ops::DerefMut for LoopsStore<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P, C> FromIterator<BoundaryWire<P, C>> for Loops<P, C> {
    #[inline(always)]
    fn from_iter<I: IntoIterator<Item = BoundaryWire<P, C>>>(iter: I) -> Self {
        Self(Vec::from_iter(iter))
    }
}

impl<'a, P, C, S> From<&'a Face<P, C, S>> for Loops<P, C> {
    #[inline(always)]
    fn from(face: &'a Face<P, C, S>) -> Loops<P, C> {
        face.absolute_boundaries()
            .iter()
            .map(|wire| BoundaryWire::new(wire.clone(), ShapesOpStatus::Unknown))
            .collect()
    }
}

impl<'a, P: 'a, C: 'a, S: 'a> FromIterator<&'a Face<P, C, S>> for LoopsStore<P, C> {
    fn from_iter<I: IntoIterator<Item = &'a Face<P, C, S>>>(iter: I) -> Self {
        Self(iter.into_iter().map(Loops::from).collect())
    }
}

impl<'a, P, C> IntoIterator for &'a LoopsStore<P, C> {
    type Item = <&'a Vec<Loops<P, C>> as IntoIterator>::Item;
    type IntoIter = <&'a Vec<Loops<P, C>> as IntoIterator>::IntoIter;
    fn into_iter(self) -> Self::IntoIter { self.0.iter() }
}

#[derive(Clone, Debug, Copy, PartialEq)]
enum ParameterKind {
    Front,
    Back,
    Inner(f64),
}

impl ParameterKind {
    fn try_new(t: f64, (t0, t1): (f64, f64)) -> Option<ParameterKind> {
        if t0.near(&t) {
            Some(ParameterKind::Front)
        } else if t1.near(&t) {
            Some(ParameterKind::Back)
        } else if t0 < t && t < t1 {
            Some(ParameterKind::Inner(t))
        } else {
            None
        }
    }
}

impl<P: Copy, C: Clone> Loops<P, C> {
    fn search_parameter(&self, pt: P) -> Option<(usize, usize, ParameterKind)>
    where C: BoundedCurve<Point = P> + SearchParameter<CurveParameter, Point = P> {
        self.iter()
            .enumerate()
            .flat_map(move |(i, wire)| wire.iter().enumerate().map(move |(j, edge)| (i, j, edge)))
            .find_map(|(i, j, edge)| {
                let curve = edge.curve();
                curve.search_parameter(pt, None, 1).and_then(|t| {
                    let kind = ParameterKind::try_new(t, curve.range_tuple())?;
                    Some((i, j, kind))
                })
            })
    }

    fn change_vertex(
        &mut self,
        old_vertex: &Vertex<P>,
        new_vertex: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) {
        self.iter_mut().for_each(|wire| {
            wire.retain_mut(|edge| {
                let edge_id = edge.id();
                let absolute_front = edge.absolute_front();
                let absolute_back = edge.absolute_back();
                let replaces_front = absolute_front == old_vertex;
                let replaces_back = absolute_back == old_vertex;
                if !replaces_front && !replaces_back {
                    return true;
                }

                // Vertex coalescing can make an already-split boundary segment
                // degenerate. This occurs when coplanar imprinting discovers
                // that two independently introduced topology vertices are the
                // same exact boundary point. The correct merged topology omits
                // that zero-length edge; constructing a self-edge would be both
                // invalid and, through `Edge::new`, a panic.
                if (replaces_front && absolute_back == new_vertex)
                    || (replaces_back && absolute_front == new_vertex)
                {
                    if edge.count() == 1 {
                        emap.remove(&edge_id);
                    }
                    return false;
                }

                let mut new_edge = emap
                    .entry(edge_id)
                    .or_insert_with(|| {
                        if replaces_front {
                            Edge::new_unchecked(new_vertex, absolute_back, edge.curve())
                        } else {
                            Edge::new_unchecked(absolute_front, new_vertex, edge.curve())
                        }
                    })
                    .clone();
                if !edge.orientation() {
                    new_edge.invert();
                }
                // Remove the edge from the HashMap when it is no longer there because Id reassignment will occur.
                if edge.count() == 1 {
                    emap.remove(&edge_id);
                }
                *edge = new_edge;
                true
            });
        });
    }

    fn swap_edge_into_wire(&mut self, edge_id: EdgeId<C>, new_wire: &Wire<P, C>) {
        self.iter_mut().for_each(|wire| {
            let mut iter = wire.iter().enumerate();
            if let Some((idx, edge)) = iter.find(|(_, edge)| edge.id() == edge_id) {
                if edge.orientation() {
                    wire.swap_edge_into_wire(idx, new_wire.clone());
                } else {
                    wire.swap_edge_into_wire(idx, new_wire.inverse());
                }
            }
        });
    }

    #[inline(always)]
    fn add_independent_loop(&mut self, r#loop: BoundaryWire<P, C>) {
        self.push(r#loop.inverse());
        self.push(r#loop);
    }

    fn reuse_coincident_edge(&mut self, edge0: &Edge<P, C>) -> bool
    where
        P: Tolerance,
        C: ParametricCurve<Point = P>
            + BoundedCurve<Point = P>
            + SearchNearestParameter<CurveParameter, Point = P>
            + Invertible, {
        let curves_sampled_match = |edge: &Edge<P, C>| {
            validate_sampled_coincidence::<P, _>(&edge.oriented_curve(), &edge0.oriented_curve())
                == CurveCoincidence::SampledValidation
        };
        let coincident = self.iter().flat_map(|wire| wire.iter()).find_map(|edge| {
            let same_direction = edge.front().point().near(&edge0.front().point())
                && edge.back().point().near(&edge0.back().point());
            let opposite_direction = edge.front().point().near(&edge0.back().point())
                && edge.back().point().near(&edge0.front().point());
            (same_direction || opposite_direction)
                .then(|| curves_sampled_match(edge))
                .filter(|matched| *matched)
                .map(|_| (edge.front().clone(), edge.back().clone(), same_direction))
        });
        let Some((old_front, old_back, same_direction)) = coincident else {
            return false;
        };
        let (new_front, new_back) = if same_direction {
            (edge0.front().clone(), edge0.back().clone())
        } else {
            (edge0.back().clone(), edge0.front().clone())
        };
        let mut emap = HashMap::default();
        if old_front.id() != new_front.id() {
            self.change_vertex(&old_front, &new_front, &mut emap);
        }
        if old_back.id() != new_back.id() {
            self.change_vertex(&old_back, &new_back, &mut emap);
        }
        self.iter_mut().for_each(|wire| {
            wire.iter_mut().for_each(|edge| {
                let same_direction = edge.front().point().near(&edge0.front().point())
                    && edge.back().point().near(&edge0.back().point());
                let opposite_direction = edge.front().point().near(&edge0.back().point())
                    && edge.back().point().near(&edge0.front().point());
                if same_direction {
                    *edge = edge0.clone();
                } else if opposite_direction {
                    *edge = edge0.inverse();
                }
            });
        });
        true
    }

    fn add_edge(
        &mut self,
        edge0: Edge<P, C>,
        status: ShapesOpStatus,
        reuse_coincident: bool,
    ) -> [Option<(usize, usize)>; 2]
    where
        P: Tolerance,
        C: ParametricCurve<Point = P>
            + BoundedCurve<Point = P>
            + SearchNearestParameter<CurveParameter, Point = P>
            + Invertible,
    {
        if reuse_coincident && self.reuse_coincident_edge(&edge0) {
            return [None, None];
        }

        let matching_occurrences = |vertex: &Vertex<P>| {
            self.iter()
                .enumerate()
                .flat_map(|(wire_index, wire)| {
                    wire.iter()
                        .enumerate()
                        .filter_map(move |(edge_index, edge)| {
                            (edge.front() == vertex).then_some((wire_index, edge_index))
                        })
                })
                .collect::<Vec<_>>()
        };
        let a_occurrences = matching_occurrences(edge0.back());
        let b_occurrences = matching_occurrences(edge0.front());

        // A divider endpoint can legitimately occur in several sibling wires after
        // an earlier divider has split a face. Choosing the first occurrence of
        // each endpoint independently can therefore select two different siblings
        // even when both endpoints coexist on one sibling, incorrectly merging the
        // sibling regions into one non-simple walk. Prefer the first deterministic
        // same-wire pairing; preserve the historical cross-wire behavior only when
        // no such pairing exists.
        let same_wire_pair = a_occurrences.iter().find_map(|&a| {
            b_occurrences
                .iter()
                .find(|&&(wire_index, _)| wire_index == a.0)
                .copied()
                .map(|b| (a, b))
        });
        let (a, b) = match same_wire_pair {
            Some((a, b)) => (Some(a), Some(b)),
            None => (
                a_occurrences.first().copied(),
                b_occurrences.first().copied(),
            ),
        };
        if let Some((wire_index0, edge_index0)) = a {
            self[wire_index0].rotate_left(edge_index0);
            self[wire_index0].push_front(edge0.clone());
            self[wire_index0].push_back(edge0.inverse());
        }
        match (a, b) {
            (Some((wire_index0, edge_index0)), Some((wire_index1, edge_index1))) => {
                if wire_index0 == wire_index1 {
                    let len = self[wire_index0].len() - 2;
                    let edge_index1 = (len + edge_index1 - edge_index0) % len + 1;
                    let new_wire = self[wire_index0].split_off(edge_index1);
                    self[wire_index0].status = status;
                    self.push(BoundaryWire::new(new_wire, status.not()));
                } else {
                    let mut new_wire0 = self[wire_index1].clone();
                    let mut new_wire1 = new_wire0.split_off(edge_index1);
                    new_wire0.append(&mut self[wire_index0]);
                    new_wire0.append(&mut new_wire1);
                    self[wire_index0] = new_wire0;
                    self.swap_remove(wire_index1);
                }
            }
            (None, Some((wire_index1, edge_index1))) => {
                self[wire_index1].rotate_left(edge_index1);
                self[wire_index1].push_front(edge0.inverse());
                self[wire_index1].push_back(edge0);
            }
            (None, None) => self.push(BoundaryWire::new(
                vec![edge0.inverse(), edge0].into(),
                ShapesOpStatus::Unknown,
            )),
            _ => {}
        }
        [a, b]
    }
}

impl<C> Loops<Point3, C> {
    fn has_partitioned_coincident_path(&self, edge0: &Edge<Point3, C>) -> bool
    where C: Clone
            + ParametricCurve<Point = Point3>
            + BoundedCurve<Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3> {
        let reference = edge0.curve();
        let reference_range = reference.range_tuple();
        let start = edge0.front().id();
        let goal = edge0.back().id();

        let mut edge_indices = HashMap::default();
        let mut edges = Vec::new();
        let mut adjacency = HashMap::<VertexId<Point3>, Vec<usize>>::default();
        for edge in self.iter().flat_map(|wire| wire.iter()) {
            let edge_id = edge.id();
            if edge_indices.contains_key(&edge_id) {
                continue;
            }
            let edge = edge.absolute_clone();
            let index = edges.len();
            let front = edge.absolute_front().id();
            let back = edge.absolute_back().id();
            edge_indices.insert(edge_id, index);
            edges.push(edge);
            adjacency.entry(front).or_default().push(index);
            adjacency.entry(back).or_default().push(index);
        }

        if !adjacency.contains_key(&start) || !adjacency.contains_key(&goal) {
            return false;
        }

        // Geometric coincidence is expensive for intersection curves. Explore the
        // topological graph first and prove only edges reachable from the incoming
        // edge's start. This is equivalent to filtering every edge up front, while
        // avoiding inverse-mapping unrelated wires and branches.
        let mut coincident = vec![None; edges.len()];
        let mut frontier = vec![start];
        let mut seen = HashMap::default();
        seen.insert(start, ());
        while let Some(current) = frontier.pop() {
            let Some(incident) = adjacency.get(&current) else {
                continue;
            };
            for &index in incident {
                let edge = &edges[index];
                let front = edge.absolute_front().id();
                let back = edge.absolute_back().id();
                // Leave the single-edge case to `add_edge`, which also
                // canonicalizes the shared EdgeId across the two operands.
                if (front == start && back == goal) || (front == goal && back == start) {
                    continue;
                }
                let next = if front == current {
                    back
                } else if back == current {
                    front
                } else {
                    continue;
                };
                if next != goal && seen.contains_key(&next) {
                    continue;
                }

                let is_coincident = *coincident[index].get_or_insert_with(|| {
                    let curve = edge.curve();
                    let (t0, t1) = curve.range_tuple();
                    let forward = front == current;
                    let parameters = if forward {
                        [
                            (t1 - t0).mul_add(0.25, t0),
                            f64::midpoint(t0, t1),
                            (t1 - t0).mul_add(0.75, t0),
                        ]
                    } else {
                        [
                            (t1 - t0).mul_add(0.75, t0),
                            f64::midpoint(t0, t1),
                            (t1 - t0).mul_add(0.25, t0),
                        ]
                    };
                    let other_endpoint = if forward {
                        edge.absolute_back().point()
                    } else {
                        edge.absolute_front().point()
                    };
                    let points = parameters
                        .into_iter()
                        .map(|parameter| curve.evaluate(parameter))
                        .chain(std::iter::once(other_endpoint));
                    let mut hint = None;
                    points.into_iter().all(|point| {
                        let validated_parameter = |search_hint| {
                            let parameter =
                                reference.search_nearest_parameter(point, search_hint, 100)?;
                            ParameterKind::try_new(parameter, reference_range)?;
                            reference
                                .evaluate(parameter)
                                .near(&point)
                                .then_some(parameter)
                        };
                        let parameter = validated_parameter(hint)
                            .or_else(|| hint.and_then(|_| validated_parameter(None)));
                        let Some(parameter) = parameter else {
                            return false;
                        };
                        hint = Some(parameter);
                        true
                    })
                });
                if !is_coincident {
                    continue;
                }
                if next == goal {
                    return true;
                }
                if seen.insert(next, ()).is_none() {
                    frontier.push(next);
                }
            }
        }
        false
    }
}

impl<P: Copy + Tolerance, C: Clone> LoopsStore<P, C> {
    #[inline(always)]
    fn change_vertex(
        &mut self,
        old_vertex: &Vertex<P>,
        new_vertex: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) {
        self.iter_mut()
            .for_each(|loops| loops.change_vertex(old_vertex, new_vertex, emap));
    }

    #[inline(always)]
    fn swap_edge_into_wire(&mut self, edge_id: EdgeId<C>, new_wire: &Wire<P, C>) {
        self.iter_mut()
            .for_each(|loops| loops.swap_edge_into_wire(edge_id, new_wire))
    }
}

impl<C> LoopsStore<Point3, C> {
    fn search_boundary_parameter(
        &self,
        loops_index: usize,
        point: Point3,
        nearest_tolerance: f64,
    ) -> Option<(usize, usize, ParameterKind, Option<Point3>)>
    where
        C: Clone
            + BoundedCurve<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let loops = &self[loops_index];
        let mut candidate = None;
        let mut candidate_edge: Option<Edge<Point3, C>> = None;
        let mut candidate_endpoint: Option<VertexId<Point3>> = None;

        let mut accept_candidate = |wire_index: usize,
                                    edge_index: usize,
                                    edge: &Edge<Point3, C>,
                                    kind: ParameterKind,
                                    projected: Option<Point3>| {
            let endpoint = match kind {
                ParameterKind::Front => Some(edge.absolute_front().id()),
                ParameterKind::Back => Some(edge.absolute_back().id()),
                ParameterKind::Inner(_) => None,
            };
            if let Some(existing) = candidate_edge.as_ref() {
                let same_edge = existing.id() == edge.id();
                let same_endpoint = candidate_endpoint.is_some() && candidate_endpoint == endpoint;
                let same_endpoints = (existing
                    .absolute_front()
                    .point()
                    .near(&edge.absolute_front().point())
                    && existing
                        .absolute_back()
                        .point()
                        .near(&edge.absolute_back().point()))
                    || (existing
                        .absolute_front()
                        .point()
                        .near(&edge.absolute_back().point())
                        && existing
                            .absolute_back()
                            .point()
                            .near(&edge.absolute_front().point()));
                let same_curve_copy = same_endpoints
                    && (validate_sampled_coincidence::<Point3, _>(
                        &existing.curve(),
                        &edge.curve(),
                    ) == CurveCoincidence::SampledValidation
                        || validate_endpoint_aligned_coincidence::<Point3, _>(
                            &existing.curve(),
                            &edge.curve(),
                        ) == CurveCoincidence::SampledValidation);
                if same_edge || same_endpoint || same_curve_copy {
                    return true;
                }
                return false;
            }

            let canonical_point = match kind {
                ParameterKind::Front => edge.absolute_front().point(),
                ParameterKind::Back => edge.absolute_back().point(),
                ParameterKind::Inner(_) => projected.unwrap_or(point),
            };
            candidate = Some((
                wire_index,
                edge_index,
                kind,
                projected.map(|_| canonical_point),
            ));
            candidate_edge = Some(edge.clone());
            candidate_endpoint = endpoint;
            true
        };

        // Prefer exact incidence, but still reject multiple distinct topological
        // loci instead of silently taking the first edge whose inverse map succeeds.
        let mut found_exact = false;
        for (wire_index, wire) in loops.iter().enumerate() {
            for (edge_index, edge) in wire.iter().enumerate() {
                let curve = edge.curve();
                let Some(parameter) = curve.search_parameter(point, None, 1) else {
                    continue;
                };
                let Some(kind) = ParameterKind::try_new(parameter, curve.range_tuple()) else {
                    continue;
                };
                found_exact = true;
                if !accept_candidate(wire_index, edge_index, edge, kind, None) {
                    return None;
                }
            }
        }
        if found_exact {
            return candidate;
        }

        // SSI endpoints originate from tessellated carrier faces. Candidate
        // discovery may therefore use the caller's tessellation-error budget,
        // but the selected point is projected back onto exact boundary geometry
        // before topology is mutated. Distinct nearby loci remain ambiguous and
        // fail closed; only repeated uses/shared endpoints or sampled-validated
        // copies of the same complete edge are deduplicated.
        for (wire_index, wire) in loops.iter().enumerate() {
            for (edge_index, edge) in wire.iter().enumerate() {
                let curve = edge.curve();
                let Some(parameter) = curve.search_nearest_parameter(point, None, 10) else {
                    continue;
                };
                let projected = curve.evaluate(parameter);
                if projected.distance(point) > nearest_tolerance {
                    continue;
                }
                let Some(kind) = ParameterKind::try_new(parameter, curve.range_tuple()) else {
                    continue;
                };
                let canonical_point = match kind {
                    ParameterKind::Front => edge.absolute_front().point(),
                    ParameterKind::Back => edge.absolute_back().point(),
                    ParameterKind::Inner(_) => projected,
                };
                if !accept_candidate(wire_index, edge_index, edge, kind, Some(canonical_point)) {
                    return None;
                }
            }
        }
        candidate
    }

    fn add_imprint_vertex(
        &mut self,
        loops_index: usize,
        vertex: &Vertex<Point3>,
        emap: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
    ) -> Option<()>
    where
        C: Clone
            + Cut<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let point = vertex.point();
        let mut endpoint_vertices = Vec::<Vertex<Point3>>::new();
        let mut interior = Vec::<(EdgeId<C>, Edge<Point3, C>, f64)>::new();

        for wire in self[loops_index].iter() {
            for edge in wire.iter() {
                let curve = edge.curve();
                let Some(parameter) = curve.search_parameter(point, None, 1) else {
                    continue;
                };
                let Some(kind) = ParameterKind::try_new(parameter, curve.range_tuple()) else {
                    continue;
                };
                match kind {
                    ParameterKind::Front => {
                        let endpoint = edge.absolute_front().clone();
                        if !endpoint_vertices
                            .iter()
                            .any(|known| known.id() == endpoint.id())
                        {
                            endpoint_vertices.push(endpoint);
                        }
                    }
                    ParameterKind::Back => {
                        let endpoint = edge.absolute_back().clone();
                        if !endpoint_vertices
                            .iter()
                            .any(|known| known.id() == endpoint.id())
                        {
                            endpoint_vertices.push(endpoint);
                        }
                    }
                    ParameterKind::Inner(parameter) => {
                        if interior.iter().any(|(edge_id, _, _)| *edge_id == edge.id()) {
                            continue;
                        }
                        if let Some((_, reference, _)) = interior.first() {
                            let same_endpoints = (reference
                                .absolute_front()
                                .point()
                                .near(&edge.absolute_front().point())
                                && reference
                                    .absolute_back()
                                    .point()
                                    .near(&edge.absolute_back().point()))
                                || (reference
                                    .absolute_front()
                                    .point()
                                    .near(&edge.absolute_back().point())
                                    && reference
                                        .absolute_back()
                                        .point()
                                        .near(&edge.absolute_front().point()));
                            let reference_curve = reference.curve();
                            let edge_curve = edge.curve();
                            let same_curve = validate_sampled_coincidence::<Point3, _>(
                                &reference_curve,
                                &edge_curve,
                            ) == CurveCoincidence::SampledValidation
                                || validate_sampled_containment::<Point3, _>(
                                    &reference_curve,
                                    &edge_curve,
                                ) == CurveCoincidence::SampledValidation
                                || validate_sampled_containment::<Point3, _>(
                                    &edge_curve,
                                    &reference_curve,
                                ) == CurveCoincidence::SampledValidation
                                || (same_endpoints
                                    && validate_endpoint_aligned_coincidence::<Point3, _>(
                                        &reference_curve,
                                        &edge_curve,
                                    ) == CurveCoincidence::SampledValidation);
                            if !same_curve {
                                return None;
                            }
                        }
                        interior.push((edge.id(), edge.absolute_clone(), parameter));
                    }
                }
            }
        }

        if endpoint_vertices.is_empty() && interior.is_empty() {
            self.add_polygon_vertex(loops_index, vertex, emap)?;
            return Some(());
        }

        // A divided face can contain several topological representations of one
        // boundary locus, including a full edge beside an already-partitioned
        // sub-edge. Every proven-coincident copy must receive the split vertex or
        // the face set becomes unsynchronized. The scan above validates the
        // carrier relationship before any topology is mutated.
        for (edge_id, edge, parameter) in interior {
            let (edge0, edge1) = edge.cut_with_parameter(vertex, parameter)?;
            let new_wire: Wire<_, _> = vec![edge0, edge1].into();
            self.swap_edge_into_wire(edge_id, &new_wire);
        }

        for old_vertex in endpoint_vertices {
            if old_vertex.id() != vertex.id() {
                self.change_vertex(&old_vertex, vertex, emap);
            }
        }
        Some(())
    }

    fn add_polygon_vertex(
        &mut self,
        loops_index: usize,
        vertex: &Vertex<Point3>,
        emap: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
    ) -> Option<(usize, usize, ParameterKind)>
    where
        C: Clone
            + Cut<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let (wire_index, edge_index, kind, projected) =
            self.search_boundary_parameter(loops_index, vertex.point(), 8.0 * TOLERANCE)?;
        if let Some(projected) = projected {
            vertex.set_point(projected);
        }
        match kind {
            ParameterKind::Front => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_front()
                    .clone();
                self.change_vertex(&old_vertex, vertex, emap);
            }
            ParameterKind::Back => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_back()
                    .clone();
                self.change_vertex(&old_vertex, vertex, emap);
            }
            ParameterKind::Inner(parameter) => {
                let edge = self[loops_index][wire_index][edge_index].absolute_clone();
                let edge_id = edge.id();
                let (edge0, edge1) = edge.cut_with_parameter(vertex, parameter)?;
                let new_wire: Wire<_, _> = vec![edge0, edge1].into();
                self.swap_edge_into_wire(edge_id, &new_wire);
            }
        }
        Some((wire_index, edge_index, kind))
    }

    fn resolve_boundary_hint(
        &self,
        loops_index: usize,
        (wire_index, edge_index): (usize, usize),
        point: Point3,
        nearest_tolerance: f64,
    ) -> Option<(usize, usize, ParameterKind, Option<Point3>)>
    where
        C: BoundedCurve<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let edge = self.get(loops_index)?.get(wire_index)?.get(edge_index)?;
        let curve = edge.curve();
        let classify = |parameter: f64| {
            let projected = curve.evaluate(parameter);
            if projected.distance(point) > nearest_tolerance {
                return None;
            }
            let kind = ParameterKind::try_new(parameter, curve.range_tuple())?;
            let canonical_point = match kind {
                ParameterKind::Front => edge.absolute_front().point(),
                ParameterKind::Back => edge.absolute_back().point(),
                ParameterKind::Inner(_) => projected,
            };
            Some((wire_index, edge_index, kind, Some(canonical_point)))
        };

        curve
            .search_parameter(point, None, 1)
            .and_then(classify)
            .or_else(|| {
                curve
                    .search_nearest_parameter(point, None, 10)
                    .and_then(classify)
            })
    }

    fn add_geom_vertex<S>(
        &mut self,
        loops_index: usize,
        boundary_hint: Option<(usize, usize)>,
        vertex: &Vertex<Point3>,
        another_surface: &S,
        nearest_tolerance: f64,
        emap: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
    ) -> Option<()>
    where
        C: BoundedCurve<Point = Point3>
            + Cut<Point = Point3, Vector = Vector3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
        S: ParametricSurface3D
            + SearchParameter<SurfaceParameter, Point = Point3>
            + SearchNearestParameter<SurfaceParameter, Point = Point3>,
    {
        // Polygon and exact loop stores can diverge after earlier endpoint
        // canonicalization. A polygon boundary index is therefore only a
        // candidate hint: accept it when the corresponding exact edge still
        // resolves the endpoint within the tessellation budget, otherwise fall
        // back to the global ambiguity-safe exact-boundary search.
        let hinted = boundary_hint.and_then(|hint| {
            self.resolve_boundary_hint(loops_index, hint, vertex.point(), nearest_tolerance)
        });
        let (wire_index, edge_index, kind, projected) = hinted.or_else(|| {
            self.search_boundary_parameter(loops_index, vertex.point(), nearest_tolerance)
        })?;
        if let Some(projected) = projected {
            vertex.set_point(projected);
        }
        let mut coalesce_endpoint = |this: &mut Self, old_vertex: Vertex<Point3>| {
            vertex.set_point(old_vertex.point());
            this.change_vertex(&old_vertex, vertex, emap);
        };
        match kind {
            ParameterKind::Front => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_front()
                    .clone();
                coalesce_endpoint(self, old_vertex);
            }
            ParameterKind::Back => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_back()
                    .clone();
                coalesce_endpoint(self, old_vertex);
            }
            ParameterKind::Inner(_) => {
                let curve = self[loops_index][wire_index][edge_index].curve();
                // Polyline and exact curves can disagree about whether a recovered SSI
                // point is just inside a trim or exactly on its endpoint. Classify the
                // solved exact parameter before deciding whether topology needs a cut.
                let (projected, parameter, _) = curve_surface_projection(
                    &curve,
                    None,
                    another_surface,
                    None,
                    vertex.point(),
                    100,
                )?;
                match ParameterKind::try_new(parameter, curve.range_tuple())? {
                    ParameterKind::Front => {
                        let old_vertex = self[loops_index][wire_index][edge_index]
                            .absolute_front()
                            .clone();
                        coalesce_endpoint(self, old_vertex);
                    }
                    ParameterKind::Back => {
                        let old_vertex = self[loops_index][wire_index][edge_index]
                            .absolute_back()
                            .clone();
                        coalesce_endpoint(self, old_vertex);
                    }
                    ParameterKind::Inner(parameter) => {
                        vertex.set_point(projected);
                        let edge = self[loops_index][wire_index][edge_index].absolute_clone();
                        let edge_id = edge.id();
                        let (edge0, edge1) = edge.cut_with_parameter(vertex, parameter)?;
                        let new_wire: Wire<_, _> = vec![edge0, edge1].into();
                        self.swap_edge_into_wire(edge_id, &new_wire);
                    }
                }
            }
        }
        Some(())
    }
}

fn curve_surface_projection<C, S>(
    curve: &C,
    curve_hint: Option<f64>,
    surface: &S,
    surface_hint: Option<(f64, f64)>,
    point: Point3,
    trials: usize,
) -> Option<(Point3, f64, Point2)>
where
    C: ParametricCurve3D
        + SearchParameter<CurveParameter, Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3>,
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if trials == 0 {
        return None;
    }
    let t = curve
        .search_parameter(point, curve_hint, 10)
        .or_else(|| curve.search_nearest_parameter(point, curve_hint, 10))?;
    let pt0 = curve.evaluate(t);
    let (u, v) = surface
        .search_parameter(point, surface_hint, 10)
        .or_else(|| surface.search_nearest_parameter(point, surface_hint, 10))?;
    let pt1 = surface.evaluate(u, v);

    // SSI endpoints come from tessellated carrier surfaces, so a seed can be a
    // few kernel tolerances away from the exact curve/surface intersection even
    // when the two independent exact projections already agree. Accept that
    // validation before taking the Newton step below: near tangencies make
    // `curve.derivative(t).dot(surface.normal(u, v))` singular and an otherwise
    // solved endpoint can be thrown arbitrarily far away. The 8x envelope is
    // the same fixed kernel-precision bound used to canonicalize adjacent SSI
    // endpoints; it deliberately does not grow with the caller's modeling
    // tolerance.
    let seed_tolerance = 8.0 * TOLERANCE;
    if pt0.near(&pt1)
        && point.distance(pt0) <= seed_tolerance
        && point.distance(pt1) <= seed_tolerance
    {
        // Keep the point exactly on the boundary curve. `pt0.near(pt1)` proves
        // that it is also on the other surface to kernel precision, while `t`
        // remains the exact split parameter for the edge below.
        return Some((pt0, t, Point2::new(u, v)));
    }

    let l = curve.derivative(t);
    let n = surface.normal(u, v);
    let denominator = l.dot(n);
    if !denominator.is_finite() || denominator == 0.0 {
        return None;
    }
    let t0 = (pt1 - pt0).dot(n) / denominator;
    curve_surface_projection(
        curve,
        Some(t),
        surface,
        Some((u, v)),
        pt0 + t0 * l,
        trials - 1,
    )
}

fn unique_nearby_loop_point<C>(
    lhs: &Loops<Point3, C>,
    rhs: &Loops<Point3, C>,
    point: Point3,
    tolerance: f64,
) -> Option<Point3> {
    let mut candidate = None;
    for nearby in lhs
        .iter()
        .chain(rhs.iter())
        .flat_map(|wire| wire.vertex_iter())
        .map(|vertex| vertex.point())
        .filter(|nearby| (*nearby - point).magnitude() <= tolerance)
    {
        match candidate {
            None => candidate = Some(nearby),
            Some(existing) if existing.near(&nearby) => {}
            Some(_) => return None,
        }
    }
    candidate
}

fn create_independent_loop<P, C, D>(mut poly_curve0: C) -> Wire<P, D>
where
    C: Cut<Point = P>,
    D: From<C>, {
    let (t0, t1) = poly_curve0.range_tuple();
    let t = (t0 + t1) / 2.0;
    let poly_curve1 = poly_curve0.cut(t);
    let v0 = Vertex::new(poly_curve0.front());
    let v1 = Vertex::new(poly_curve1.front());
    let edge0 = Edge::new_unchecked(&v0, &v1, poly_curve0.into());
    let edge1 = Edge::new_unchecked(&v1, &v0, poly_curve1.into());
    wire![edge0, edge1]
}

struct OperandLoopStores<C> {
    geom: LoopsStore<Point3, C>,
    poly: LoopsStore<Point3, PolylineCurve>,
}

struct BoundaryEndpoint<'a, S> {
    face_index: usize,
    poly_vertex: &'a Vertex<Point3>,
    geom_vertex: &'a Vertex<Point3>,
    other_surface: &'a S,
    nearest_tolerance: f64,
}

impl<C> OperandLoopStores<C> {
    fn add_boundary_endpoint<S>(
        &mut self,
        endpoint: BoundaryEndpoint<'_, S>,
        poly_edge_map: &mut HashMap<EdgeId<PolylineCurve>, Edge<Point3, PolylineCurve>>,
        geom_edge_map: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
    ) -> Option<bool>
    where
        C: BoundedCurve<Point = Point3>
            + Cut<Point = Point3, Vector = Vector3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>
            + Clone,
        S: ParametricSurface3D
            + SearchParameter<SurfaceParameter, Point = Point3>
            + SearchNearestParameter<SurfaceParameter, Point = Point3>,
    {
        let BoundaryEndpoint {
            face_index,
            poly_vertex,
            geom_vertex,
            other_surface,
            nearest_tolerance,
        } = endpoint;
        let Some((wire_index, edge_index, _)) =
            self.poly
                .add_polygon_vertex(face_index, poly_vertex, poly_edge_map)
        else {
            return Some(false);
        };
        self.geom.add_geom_vertex(
            face_index,
            Some((wire_index, edge_index)),
            geom_vertex,
            other_surface,
            nearest_tolerance,
            geom_edge_map,
        )?;
        Some(true)
    }

    fn add_intersection_edge(
        &mut self,
        face_index: usize,
        poly_edge: Edge<Point3, PolylineCurve>,
        geom_edge: Edge<Point3, C>,
        status: ShapesOpStatus,
        reuse_coincident: bool,
    ) where
        C: ParametricCurve<Point = Point3>
            + BoundedCurve<Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>
            + Invertible
            + Clone,
    {
        self.poly[face_index].add_edge(poly_edge, status, reuse_coincident);
        self.geom[face_index].add_edge(geom_edge, status, reuse_coincident);
    }
}

fn loop_stores_are_synchronized<C>(
    geom: &LoopsStore<Point3, C>,
    poly: &LoopsStore<Point3, PolylineCurve>,
) -> bool {
    geom.len() == poly.len()
        && geom
            .iter()
            .zip(poly.iter())
            .all(|(geom_loops, poly_loops)| {
                geom_loops.len() == poly_loops.len()
                    && geom_loops
                        .iter()
                        .zip(poly_loops.iter())
                        .all(|(geom_wire, poly_wire)| {
                            geom_wire.len() == poly_wire.len()
                                && geom_wire.iter().zip(poly_wire.iter()).all(
                                    |(geom_edge, poly_edge)| {
                                        geom_edge.front().point().near(&poly_edge.front().point())
                                            && geom_edge
                                                .back()
                                                .point()
                                                .near(&poly_edge.back().point())
                                    },
                                )
                        })
            })
}

pub(super) struct InterfaceLoops<C> {
    pub(super) geom: [LoopsStore<Point3, C>; 2],
    pub(super) sampled: [Option<LoopsStore<Point3, PolylineCurve>>; 2],
    pub(super) coplanar_overlap: bool,
    pub(super) full_interface_adjacency: Option<FullPlanarInterfaceMatch>,
}

pub(super) fn create_loops_stores<C, S>(
    geom_shell0: &Shell<Point3, C, S>,
    poly_shell0: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    geom_shell1: &Shell<Point3, C, S>,
    poly_shell1: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
    imprint_coplanar: bool,
) -> Option<InterfaceLoops<C>>
where
    C: SearchNearestParameter<CurveParameter, Point = Point3>
        + SearchParameter<CurveParameter, Point = Point3>
        + Cut<Point = Point3, Vector = Vector3>
        + From<IntersectionCurve<PolylineCurve, S, S>>
        + Invertible,
    S: ParametricSurface3D
        + Clone
        + Invertible
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let mut operand0 = OperandLoopStores {
        geom: geom_shell0.face_iter().collect(),
        poly: poly_shell0.face_iter().collect(),
    };
    let mut operand1 = OperandLoopStores {
        geom: geom_shell1.face_iter().collect(),
        poly: poly_shell1.face_iter().collect(),
    };
    let store0_len = operand0.geom.len();
    let store1_len = operand1.geom.len();
    let boolean_tolerance = BooleanTolerance::new(tol);
    let coplanar_tol = 2.0 * boolean_tolerance.length;
    let boundary_search_tolerance = boolean_tolerance.mesh_seed;
    let positive_aabb_overlap =
        shells_have_positive_aabb_overlap(poly_shell0, poly_shell1, coplanar_tol)?;
    let full_coplanar_interface = if imprint_coplanar {
        shells_full_planar_interface_pair(poly_shell0, poly_shell1, coplanar_tol)?
    } else {
        None
    };
    let full_interface_adjacency = if positive_aabb_overlap {
        None
    } else {
        full_coplanar_interface.clone()
    };
    let partial_coplanar_overlap = imprint_coplanar
        && positive_aabb_overlap
        && shells_have_partial_planar_overlap(poly_shell0, poly_shell1, coplanar_tol)?;
    let coplanar_overlap = full_coplanar_interface.is_some() || partial_coplanar_overlap;
    (0..store0_len)
        .flat_map(move |i| (0..store1_len).map(move |j| (i, j)))
        .try_for_each(|(face_index0, face_index1)| {
            // Solids that meet on one complete trimmed planar face with
            // opposite outward normals need no SSI at all. Generic SSI on
            // their coincident continuation surfaces (e.g. coaxial cylinders)
            // is degenerate; downstream coplanar classification removes the
            // shared internal face and topology welding stitches the rim.
            if full_interface_adjacency.is_some() {
                return Some(());
            }
            let ori0 = geom_shell0[face_index0].orientation();
            let ori1 = geom_shell1[face_index1].orientation();
            let surface0 = geom_shell0[face_index0].surface();
            let surface1 = geom_shell1[face_index1].surface();
            let polygon0 = poly_shell0[face_index0].surface()?;
            let polygon1 = poly_shell1[face_index1].surface()?;
            // Exactly coincident trimmed planar faces are already a complete
            // interface. Running generic SSI on coincident surfaces is both
            // unnecessary and degenerate; ownership is resolved by the
            // coplanar face classifier after division.
            if imprint_coplanar
                && planar_polygons_same_trimmed_face(&polygon0, &polygon1, coplanar_tol)
            {
                return Some(());
            }
            let curves = intersection_curve::intersection_curves(
                surface0.clone(),
                &polygon0,
                surface1.clone(),
                &polygon1,
                tol,
            )?;
            curves
                .into_iter()
                .try_for_each(|(mut polyline, mut intersection_curve)| {
                    let status = ShapesOpStatus::from_is_curve(&intersection_curve)?;
                    let (status0, status1) = match (ori0, ori1) {
                        (true, true) => (status, status.not()),
                        (true, false) => (status.not(), status.not()),
                        (false, true) => (status, status),
                        (false, false) => (status.not(), status),
                    };
                    if polyline.front().near(&polyline.back()) {
                        let poly_wire = create_independent_loop(polyline);
                        operand0.poly[face_index0]
                            .add_independent_loop(BoundaryWire::new(poly_wire.clone(), status0));
                        operand1.poly[face_index1]
                            .add_independent_loop(BoundaryWire::new(poly_wire, status1));
                        let geom_wire = create_independent_loop(intersection_curve);
                        operand0.geom[face_index0]
                            .add_independent_loop(BoundaryWire::new(geom_wire.clone(), status0));
                        operand1.geom[face_index1]
                            .add_independent_loop(BoundaryWire::new(geom_wire, status1));
                    } else {
                        // Adjacent face-pairs solve SSI independently. Their shared
                        // endpoint can differ by a few kernel tolerances even though both
                        // segments belong to one topological intersection. Snap only to a
                        // unique already-known endpoint on either face; ambiguous nearby
                        // topology remains fail-closed. Topological identity must not scale with
                        // the caller's modeling tolerance; use a fixed kernel-precision envelope.
                        let stitch_tolerance = CROSS_FACE_STITCH_TOLERANCE;
                        let mut front = polyline.front();
                        let mut back = polyline.back();
                        if let Some(point) = unique_nearby_loop_point(
                            &operand0.geom[face_index0],
                            &operand1.geom[face_index1],
                            front,
                            stitch_tolerance,
                        ) {
                            front = point;
                        }
                        if let Some(point) = unique_nearby_loop_point(
                            &operand0.geom[face_index0],
                            &operand1.geom[face_index1],
                            back,
                            stitch_tolerance,
                        ) {
                            back = point;
                        }
                        // Cross-face endpoint stitching can collapse a numerically
                        // tiny SSI fragment onto one already-known intersection
                        // vertex. Such a fragment is a point contact after
                        // canonicalization, not a one-dimensional trim. Discard it
                        // before splitting boundary edges or constructing topology.
                        if front.near(&back) {
                            return Some(());
                        }
                        *polyline.first_mut()? = front;
                        *polyline.last_mut()? = back;
                        *intersection_curve.leader_mut().first_mut()? = front;
                        *intersection_curve.leader_mut().last_mut()? = back;

                        let pv0 = Vertex::new(front);
                        let pv1 = Vertex::new(back);
                        let gv0 = Vertex::new(front);
                        let gv1 = Vertex::new(back);
                        let mut pemap0 = HashMap::default();
                        let mut pemap1 = HashMap::default();
                        let mut gemap0 = HashMap::default();
                        let mut gemap1 = HashMap::default();
                        operand0.add_boundary_endpoint(
                            BoundaryEndpoint {
                                face_index: face_index0,
                                poly_vertex: &pv0,
                                geom_vertex: &gv0,
                                other_surface: &surface1,
                                nearest_tolerance: boundary_search_tolerance,
                            },
                            &mut pemap0,
                            &mut gemap0,
                        )?;
                        operand0.add_boundary_endpoint(
                            BoundaryEndpoint {
                                face_index: face_index0,
                                poly_vertex: &pv1,
                                geom_vertex: &gv1,
                                other_surface: &surface1,
                                nearest_tolerance: boundary_search_tolerance,
                            },
                            &mut pemap1,
                            &mut gemap1,
                        )?;
                        operand1.add_boundary_endpoint(
                            BoundaryEndpoint {
                                face_index: face_index1,
                                poly_vertex: &pv0,
                                geom_vertex: &gv0,
                                other_surface: &surface0,
                                nearest_tolerance: boundary_search_tolerance,
                            },
                            &mut pemap0,
                            &mut gemap0,
                        )?;
                        operand1.add_boundary_endpoint(
                            BoundaryEndpoint {
                                face_index: face_index1,
                                poly_vertex: &pv1,
                                geom_vertex: &gv1,
                                other_surface: &surface0,
                                nearest_tolerance: boundary_search_tolerance,
                            },
                            &mut pemap1,
                            &mut gemap1,
                        )?;
                        *intersection_curve.leader_mut().first_mut()? = gv0.point();
                        *intersection_curve.leader_mut().last_mut()? = gv1.point();

                        let pedge = Edge::new_unchecked(&pv0, &pv1, polyline);
                        let gedge = Edge::new_unchecked(&gv0, &gv1, intersection_curve.into());
                        operand0.add_intersection_edge(
                            face_index0,
                            pedge.clone(),
                            gedge.clone(),
                            status0,
                            partial_coplanar_overlap,
                        );
                        operand1.add_intersection_edge(
                            face_index1,
                            pedge,
                            gedge,
                            status1,
                            partial_coplanar_overlap,
                        );
                    }
                    Some(())
                })
        })?;

    // Mesh/mesh interference intentionally ignores coplanar contact. By this
    // point ordinary SSI has already split the boundary edges at transverse
    // intersections, so imprint those existing edge segments onto any face of
    // the opposite operand that contains them. This supplies the missing trim
    // lines for coplanar overlap without manufacturing approximate SSI curves.
    if partial_coplanar_overlap {
        let source0 = operand0.geom.clone();
        let source1 = operand1.geom.clone();
        imprint_edges_on_faces(&source1, poly_shell0, &mut operand0.geom, tol)?;
        imprint_edges_on_faces(&source0, poly_shell1, &mut operand1.geom, tol)?;
    }

    let sampled0 = (!partial_coplanar_overlap
        && loop_stores_are_synchronized(&operand0.geom, &operand0.poly))
    .then_some(operand0.poly);
    let sampled1 = (!partial_coplanar_overlap
        && loop_stores_are_synchronized(&operand1.geom, &operand1.poly))
    .then_some(operand1.poly);
    Some(InterfaceLoops {
        geom: [operand0.geom, operand1.geom],
        sampled: [sampled0, sampled1],
        coplanar_overlap,
        full_interface_adjacency,
    })
}

#[cfg(test)]
mod projection_tests {
    use super::*;

    #[test]
    fn synchronized_sampling_evidence_rejects_structural_drift() {
        let g0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let g1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let geom_edge = Edge::new_unchecked(&g0, &g1, Line(g0.point(), g1.point()));
        let geom = LoopsStore(vec![Loops(vec![BoundaryWire::new(
            vec![geom_edge].into(),
            ShapesOpStatus::Unknown,
        )])]);

        let p0 = Vertex::new(g0.point());
        let p1 = Vertex::new(g1.point());
        let poly_edge = Edge::new_unchecked(
            &p0,
            &p1,
            PolylineCurve(vec![p0.point(), Point3::new(0.5, 0.0, 2.0e-4), p1.point()]),
        );
        let poly = LoopsStore(vec![Loops(vec![BoundaryWire::new(
            vec![poly_edge].into(),
            ShapesOpStatus::Unknown,
        )])]);
        assert!(loop_stores_are_synchronized(&geom, &poly));

        let stale_back = Vertex::new(Point3::new(1.0, 0.01, 0.0));
        let stale_edge = Edge::new_unchecked(
            &p0,
            &stale_back,
            PolylineCurve(vec![p0.point(), stale_back.point()]),
        );
        let stale = LoopsStore(vec![Loops(vec![BoundaryWire::new(
            vec![stale_edge].into(),
            ShapesOpStatus::Unknown,
        )])]);
        assert!(!loop_stores_are_synchronized(&geom, &stale));
    }

    #[test]
    fn change_vertex_removes_edge_collapsed_by_coalescing() {
        let v0 = Vertex::new(0usize);
        let v1 = Vertex::new(0usize);
        let v2 = Vertex::new(1usize);
        let edge01 = Edge::new_unchecked(&v0, &v1, ());
        let edge12 = Edge::new_unchecked(&v1, &v2, ());
        let edge20 = Edge::new_unchecked(&v2, &v0, ());
        let wire = BoundaryWire::new(vec![edge01, edge12, edge20].into(), ShapesOpStatus::Unknown);
        let mut loops = Loops(vec![wire]);
        let mut emap = HashMap::default();

        loops.change_vertex(&v1, &v0, &mut emap);

        assert_eq!(loops[0].len(), 2);
        assert!(loops[0].iter().all(|edge| edge.front() != edge.back()));
        assert!(loops[0].vertex_iter().all(|vertex| vertex != v1));
    }

    #[test]
    fn loop_boundary_search_accepts_unique_bounded_nearest_projection() {
        let v0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let v1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * TOLERANCE, 0.0);

        let (_, _, kind, projected) = store
            .search_boundary_parameter(0, point, 8.0 * TOLERANCE)
            .expect("unique bounded nearest projection should recover the boundary parameter");

        assert!(matches!(kind, ParameterKind::Inner(_)));
        assert_eq!(projected, Some(Point3::new(0.5, 0.0, 0.0)));
    }

    #[test]
    fn loop_boundary_vertex_insertion_snaps_recovered_projection_before_cut() {
        let v0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let v1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let mut store = LoopsStore(vec![loops]);
        let vertex = Vertex::new(Point3::new(0.5, 1.125 * TOLERANCE, 0.0));
        let mut emap = HashMap::default();

        let (_, _, kind) = store
            .add_polygon_vertex(0, &vertex, &mut emap)
            .expect("recovered boundary projection should split the exact edge");

        assert!(matches!(kind, ParameterKind::Inner(_)));
        assert_eq!(vertex.point(), Point3::new(0.5, 0.0, 0.0));
        assert_eq!(store[0][0].len(), 2);
    }

    #[test]
    fn exact_geometry_resolves_boundary_endpoint_independently() {
        let front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let old_back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&front, &old_back, Line(front.point(), old_back.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let mut store = LoopsStore(vec![loops]);
        let recovered = Vertex::new(old_back.point());
        let mut emap = HashMap::default();

        let inserted = store.add_geom_vertex(
            0,
            None,
            &recovered,
            &Plane::xy(),
            8.0 * TOLERANCE,
            &mut emap,
        );
        assert!(
            inserted.is_some(),
            "exact endpoint must be reclassified before topology mutation"
        );

        assert_eq!(store[0][0].len(), 1);
        assert!(store[0][0].vertex_iter().any(|vertex| vertex == recovered));
        assert!(store[0][0].vertex_iter().all(|vertex| vertex != old_back));
    }

    #[test]
    fn loop_boundary_search_uses_explicit_tessellation_candidate_budget() {
        let v0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let v1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 3.0e-4, 0.0);

        assert!(
            store
                .search_boundary_parameter(0, point, 8.0 * TOLERANCE)
                .is_none()
        );
        let recovered = store.search_boundary_parameter(0, point, 1.0e-3);
        assert!(
            recovered.is_some(),
            "mesh-error budget may discover a unique exact boundary candidate"
        );
        let (_, _, kind, projected) = recovered.unwrap_or((0, 0, ParameterKind::Front, None));
        assert!(matches!(kind, ParameterKind::Inner(_)));
        assert_eq!(projected, Some(Point3::new(0.5, 0.0, 0.0)));
    }

    #[test]
    fn loop_boundary_search_deduplicates_distinct_vertices_at_same_geometric_endpoint() {
        let a0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let a1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let b0 = Vertex::new(a0.point());
        let b1 = Vertex::new(a1.point());
        let edge_a = Edge::new_unchecked(&a0, &a1, Line(a0.point(), a1.point()));
        let edge_b = Edge::new_unchecked(&b0, &b1, Line(b0.point(), b1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge_a, edge_b].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.0, 3.0e-4, 0.0);

        let recovered = store.search_boundary_parameter(0, point, 1.0e-3);
        assert!(
            recovered.is_some(),
            "coincident endpoint copies are one geometric boundary locus"
        );
        let (_, _, kind, projected) = recovered.unwrap_or((0, 0, ParameterKind::Front, None));
        assert!(matches!(kind, ParameterKind::Front));
        assert_eq!(projected, Some(a0.point()));
    }

    #[test]
    fn loop_boundary_search_rejects_distinct_exact_edges_at_same_coordinate() {
        let horizontal_front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let horizontal_back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let vertical_front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let vertical_back = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let horizontal = Edge::new_unchecked(
            &horizontal_front,
            &horizontal_back,
            Line(horizontal_front.point(), horizontal_back.point()),
        );
        let vertical = Edge::new_unchecked(
            &vertical_front,
            &vertical_back,
            Line(vertical_front.point(), vertical_back.point()),
        );
        let loops = Loops(vec![
            BoundaryWire::new(vec![horizontal].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![vertical].into(), ShapesOpStatus::Unknown),
        ]);
        let store = LoopsStore(vec![loops]);

        assert!(
            store
                .search_boundary_parameter(0, Point3::new(0.0, 0.0, 0.0), 1.0e-3)
                .is_none()
        );
    }

    #[test]
    fn imprint_splits_all_proven_copies_of_one_boundary_edge() {
        let a0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let a1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let b0 = Vertex::new(a0.point());
        let b1 = Vertex::new(a1.point());
        let edge_a = Edge::new_unchecked(&a0, &a1, Line(a0.point(), a1.point()));
        let edge_b = Edge::new_unchecked(&b0, &b1, Line(b0.point(), b1.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![edge_a].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge_b].into(), ShapesOpStatus::Unknown),
        ]);
        let mut store = LoopsStore(vec![loops]);
        let split = Vertex::new(Point3::new(0.5, 0.0, 0.0));
        let mut emap = HashMap::default();

        assert!(store.add_imprint_vertex(0, &split, &mut emap).is_some());
        assert_eq!(store[0][0].len(), 2);
        assert_eq!(store[0][1].len(), 2);
        assert!(
            store[0]
                .iter()
                .all(|wire| wire.vertex_iter().any(|v| v == split))
        );
    }

    #[test]
    fn imprint_splits_full_edge_and_partitioned_subedge_on_same_locus() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let d = Vertex::new(Point3::new(3.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let c = Vertex::new(Point3::new(2.0, 0.0, 0.0));
        let full = Edge::new_unchecked(&a, &d, Line(a.point(), d.point()));
        let sub = Edge::new_unchecked(&b, &c, Line(b.point(), c.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![full].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![sub].into(), ShapesOpStatus::Unknown),
        ]);
        let mut store = LoopsStore(vec![loops]);
        let split = Vertex::new(Point3::new(1.5, 0.0, 0.0));
        let mut emap = HashMap::default();

        assert!(store.add_imprint_vertex(0, &split, &mut emap).is_some());
        assert_eq!(store[0][0].len(), 2);
        assert_eq!(store[0][1].len(), 2);
        assert!(
            store[0]
                .iter()
                .all(|wire| wire.vertex_iter().any(|v| v == split))
        );
    }

    #[test]
    fn imprint_rejects_distinct_crossing_boundary_edges() {
        let h0 = Vertex::new(Point3::new(-1.0, 0.0, 0.0));
        let h1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let v0 = Vertex::new(Point3::new(0.0, -1.0, 0.0));
        let v1 = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let horizontal = Edge::new_unchecked(&h0, &h1, Line(h0.point(), h1.point()));
        let vertical = Edge::new_unchecked(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![horizontal].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![vertical].into(), ShapesOpStatus::Unknown),
        ]);
        let mut store = LoopsStore(vec![loops]);
        let split = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let mut emap = HashMap::default();

        assert!(store.add_imprint_vertex(0, &split, &mut emap).is_none());
        assert_eq!(store[0][0].len(), 1);
        assert_eq!(store[0][1].len(), 1);
    }

    #[test]
    fn loop_boundary_search_deduplicates_inverse_use_of_same_edge() {
        let v0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let v1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge.clone(), edge.inverse()].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * TOLERANCE, 0.0);

        assert!(
            store
                .search_boundary_parameter(0, point, 8.0 * TOLERANCE)
                .is_some()
        );
    }

    #[test]
    fn loop_boundary_search_accepts_distinct_edges_at_one_shared_vertex() {
        let left = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let shared = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let right = Vertex::new(Point3::new(2.0, 0.0, 0.0));
        let edge_left = Edge::new_unchecked(&left, &shared, Line(left.point(), shared.point()));
        let edge_right = Edge::new_unchecked(&shared, &right, Line(shared.point(), right.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge_left, edge_right].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(1.0, 1.125 * TOLERANCE, 0.0);

        let (_, _, kind, projected) = store
            .search_boundary_parameter(0, point, 8.0 * TOLERANCE)
            .expect("incident edges sharing one vertex are one boundary locus");

        assert!(matches!(kind, ParameterKind::Front | ParameterKind::Back));
        assert_eq!(projected, Some(shared.point()));
    }

    #[test]
    fn loop_boundary_search_rejects_ambiguous_bounded_nearest_projection() {
        let a0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let a1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let b0 = Vertex::new(Point3::new(0.0, 2.25 * TOLERANCE, 0.0));
        let b1 = Vertex::new(Point3::new(1.0, 2.25 * TOLERANCE, 0.0));
        let edge_a = Edge::new_unchecked(&a0, &a1, Line(a0.point(), a1.point()));
        let edge_b = Edge::new_unchecked(&b0, &b1, Line(b0.point(), b1.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![edge_a].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge_b].into(), ShapesOpStatus::Unknown),
        ]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * TOLERANCE, 0.0);

        assert!(
            store
                .search_boundary_parameter(0, point, 8.0 * TOLERANCE)
                .is_none()
        );
    }

    #[test]
    fn cross_face_endpoint_stitch_rejects_model_scale_nearby_topology() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&a, &b, Line(a.point(), b.point()));
        let lhs = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let rhs = Loops(Vec::new());

        let kernel_near = Point3::new(0.0, 0.0, 4.0 * TOLERANCE);
        assert_eq!(
            unique_nearby_loop_point(&lhs, &rhs, kernel_near, CROSS_FACE_STITCH_TOLERANCE),
            Some(a.point())
        );

        // The high-Z clip uses a modeling tolerance of 0.05. Scaling the old
        // 8x stitch envelope by that tolerance made z=0.8 and z=1.0 topology
        // close enough to be incorrectly canonicalized as one endpoint.
        let model_scale_distinct = Point3::new(0.0, 0.0, 0.2);
        assert!(
            unique_nearby_loop_point(
                &lhs,
                &rhs,
                model_scale_distinct,
                CROSS_FACE_STITCH_TOLERANCE
            )
            .is_none()
        );
    }

    #[test]
    fn planar_partition_ignores_edge_touching_neighbor() {
        fn rectangle(x0: f64, x1: f64) -> PolygonMesh {
            PolygonMesh::new(
                StandardAttributes {
                    positions: vec![
                        Point3::new(x0, 0.0, 0.0),
                        Point3::new(x1, 0.0, 0.0),
                        Point3::new(x1, 1.0, 0.0),
                        Point3::new(x0, 1.0, 0.0),
                    ],
                    ..Default::default()
                },
                Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
            )
            .expect("fixture mesh is valid")
        }

        let target = rectangle(0.0, 2.0);
        let fragments: Shell<Point3, PolylineCurve, Option<PolygonMesh>> = Shell::from(vec![
            Face::new_unchecked(Vec::new(), Some(rectangle(0.0, 1.0))),
            Face::new_unchecked(Vec::new(), Some(rectangle(1.0, 2.0))),
            // This face shares x=2 with the target but contributes no target area.
            Face::new_unchecked(Vec::new(), Some(rectangle(2.0, 3.0))),
        ]);

        assert!(
            single_face_matches_planar_partition(&target, &fragments, 1.0e-6)
                .flatten()
                .is_some()
        );
    }

    #[test]
    fn add_edge_prefers_same_sibling_wire_for_shared_endpoint() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(0.0, -1.0, 0.0));
        let c = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let d = Vertex::new(Point3::new(1.0, 1.0, 0.0));
        let e = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let f = Vertex::new(Point3::new(1.0, -1.0, 0.0));
        let line = |front: &Vertex<Point3>, back: &Vertex<Point3>| {
            Edge::new_unchecked(front, back, Line(front.point(), back.point()))
        };

        // C belongs to both sibling regions, while B belongs only to the lower
        // region. The incoming B->C divider must therefore split the lower
        // sibling rather than merging the two siblings through C.
        let upper = BoundaryWire::new(
            vec![line(&c, &d), line(&d, &e), line(&e, &a), line(&a, &c)].into(),
            ShapesOpStatus::Or,
        );
        let lower = BoundaryWire::new(
            vec![line(&a, &b), line(&b, &f), line(&f, &c), line(&c, &a)].into(),
            ShapesOpStatus::And,
        );
        let mut loops = Loops(vec![upper, lower]);
        let incoming = line(&b, &c);

        let _ = loops.add_edge(incoming, ShapesOpStatus::And, false);

        // The operation must split the lower sibling into the two triangles
        // A-B-C and B-F-C without merging the upper A-C-D-E region.
        assert_eq!(loops.len(), 3, "one sibling region should split into two");
        assert!(
            loops
                .iter()
                .all(|wire| wire.is_closed() && wire.is_simple())
        );

        let has = |wire: &BoundaryWire<_, _>, vertex: &Vertex<Point3>| {
            wire.vertex_iter().any(|item| item == *vertex)
        };
        let upper = loops
            .iter()
            .find(|wire| has(wire, &d) || has(wire, &e))
            .expect("upper sibling must survive");
        assert_eq!(upper.status(), ShapesOpStatus::Or);
        assert!(has(upper, &a) && has(upper, &c));
        assert!(!has(upper, &b) && !has(upper, &f));

        let lower_pieces = loops
            .iter()
            .filter(|wire| has(wire, &b) && has(wire, &c))
            .collect::<Vec<_>>();
        assert_eq!(lower_pieces.len(), 2);
        assert!(lower_pieces.iter().any(|wire| has(wire, &a)));
        assert!(lower_pieces.iter().any(|wire| has(wire, &f)));
    }

    #[test]
    fn partitioned_coincident_path_is_already_imprinted() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let c = Vertex::new(Point3::new(2.0, 0.0, 0.0));
        let d = Vertex::new(Point3::new(3.0, 0.0, 0.0));
        let e = Vertex::new(Point3::new(3.0, 1.0, 0.0));
        let f = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let wire = BoundaryWire::new(
            vec![
                Edge::new_unchecked(&a, &b, Line(a.point(), b.point())),
                Edge::new_unchecked(&b, &c, Line(b.point(), c.point())),
                Edge::new_unchecked(&c, &d, Line(c.point(), d.point())),
                Edge::new_unchecked(&d, &e, Line(d.point(), e.point())),
                Edge::new_unchecked(&e, &f, Line(e.point(), f.point())),
                Edge::new_unchecked(&f, &a, Line(f.point(), a.point())),
            ]
            .into(),
            ShapesOpStatus::Or,
        );
        let loops = Loops(vec![wire]);
        let incoming = Edge::new_unchecked(&a, &d, Line(a.point(), d.point()));

        assert!(loops.has_partitioned_coincident_path(&incoming));
    }

    #[test]
    fn single_coincident_edge_still_uses_edge_id_canonicalization_path() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let d = Vertex::new(Point3::new(3.0, 0.0, 0.0));
        let e = Vertex::new(Point3::new(3.0, 1.0, 0.0));
        let f = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let wire = BoundaryWire::new(
            vec![
                Edge::new_unchecked(&a, &d, Line(a.point(), d.point())),
                Edge::new_unchecked(&d, &e, Line(d.point(), e.point())),
                Edge::new_unchecked(&e, &f, Line(e.point(), f.point())),
                Edge::new_unchecked(&f, &a, Line(f.point(), a.point())),
            ]
            .into(),
            ShapesOpStatus::Or,
        );
        let loops = Loops(vec![wire]);
        let incoming = Edge::new_unchecked(&a, &d, Line(a.point(), d.point()));

        assert!(!loops.has_partitioned_coincident_path(&incoming));
    }

    #[test]
    fn curve_surface_projection_accepts_agreeing_tangent_projections() {
        let line = Line(Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0));
        let sphere = Sphere::new(Point3::origin(), 1.0);
        let seed = Point3::new(1.0 + 1.5 * TOLERANCE, 0.0, 0.0);

        let (point, parameter, _) = curve_surface_projection(&line, None, &sphere, None, seed, 100)
            .expect("agreeing exact projections should resolve a tangent SSI seed");

        assert!(point.near(&Point3::new(1.0, 0.0, 0.0)));
        assert!((parameter - 0.5).abs() <= TOLERANCE);
        assert!(sphere.include(point));
    }
}
