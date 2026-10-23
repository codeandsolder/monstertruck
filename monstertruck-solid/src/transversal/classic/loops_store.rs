//! Intersection-loop store for the classic Boolean backend.
//!
//! This began as the published 0.3.2 loop store, but now also owns exact-boundary
//! endpoint resolution, coplanar interface matching, and partition stitching.
//! Geometric equality checks here are deliberately separated from broad-phase
//! signatures: sampled or summary evidence may select candidates, never prove
//! trimmed-region identity by itself.

#![allow(clippy::many_single_char_names)]

use super::{
    curve_equivalence::{CurveCoincidence, validate_sampled_coincidence},
    intersection_curve,
    tolerance::BooleanTolerance,
};
use monstertruck_core::{cgmath64::*, tolerance::TOLERANCE};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::{Vertex, *};
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
                            Edge::new(new_vertex, absolute_back, edge.curve())
                        } else {
                            Edge::new(absolute_front, new_vertex, edge.curve())
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
        if reuse_coincident {
            let curves_sampled_match = |edge: &Edge<P, C>| {
                validate_sampled_coincidence::<P, _>(
                    &edge.oriented_curve(),
                    &edge0.oriented_curve(),
                ) == CurveCoincidence::SampledValidation
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
            if let Some((old_front, old_back, same_direction)) = coincident {
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
                return [None, None];
            }
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
        let point_on_reference = |point: Point3| {
            reference
                .search_nearest_parameter(point, None, 100)
                .and_then(|parameter| {
                    ParameterKind::try_new(parameter, reference_range)?;
                    Some(reference.evaluate(parameter).near(&point))
                })
                .unwrap_or(false)
        };

        let mut edge_ids = Vec::new();
        let mut candidates = Vec::new();
        for edge in self.iter().flat_map(|wire| wire.iter()) {
            let edge_id = edge.id();
            if edge_ids.contains(&edge_id) {
                continue;
            }
            edge_ids.push(edge_id);
            let curve = edge.curve();
            let (t0, t1) = curve.range_tuple();
            let samples = [
                edge.absolute_front().point(),
                curve.evaluate((t1 - t0).mul_add(0.25, t0)),
                curve.evaluate(f64::midpoint(t0, t1)),
                curve.evaluate((t1 - t0).mul_add(0.75, t0)),
                edge.absolute_back().point(),
            ];
            if samples.into_iter().all(point_on_reference) {
                candidates.push(edge.absolute_clone());
            }
        }

        let start = edge0.front().id();
        let goal = edge0.back().id();
        let mut frontier = vec![start];
        let mut seen = vec![start];
        while let Some(current) = frontier.pop() {
            for edge in &candidates {
                let front = edge.absolute_front().id();
                let back = edge.absolute_back().id();
                // Leave the single-edge case to `add_edge`, which also
                // canonicalizes the shared EdgeId across the two operands.
                if (front == start && back == goal) || (front == goal && back == start) {
                    continue;
                }
                let next = if front == current {
                    Some(back)
                } else if back == current {
                    Some(front)
                } else {
                    None
                };
                let Some(next) = next else { continue };
                if next == goal {
                    return true;
                }
                if !seen.contains(&next) {
                    seen.push(next);
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
    ) -> Option<(usize, usize, ParameterKind, Option<Point3>)>
    where
        C: Clone
            + BoundedCurve<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let loops = &self[loops_index];
        let exact = loops
            .iter()
            .enumerate()
            .flat_map(|(wire_index, wire)| {
                wire.iter()
                    .enumerate()
                    .map(move |(edge_index, edge)| (wire_index, edge_index, edge))
            })
            .find_map(|(wire_index, edge_index, edge)| {
                let curve = edge.curve();
                let parameter = curve.search_parameter(point, None, 1)?;
                let kind = ParameterKind::try_new(parameter, curve.range_tuple())?;
                Some((wire_index, edge_index, kind, None))
            });
        if exact.is_some() {
            return exact;
        }

        // SSI endpoints originate from tessellated carrier faces. The cleaned
        // exact-surface solution can land a few kernel tolerances away from an
        // exact topological boundary. Admit a bounded nearest projection only
        // when the nearby topology proves a unique locus. Multiple occurrences
        // of one edge, or distinct edges meeting at the same VertexId, are one
        // locus; unrelated nearby edges remain ambiguous and fail closed.
        let tolerance = 8.0 * TOLERANCE;
        let mut candidate = None;
        let mut candidate_identity: Option<(EdgeId<C>, Option<VertexId<Point3>>)> = None;
        for (wire_index, wire) in loops.iter().enumerate() {
            for (edge_index, edge) in wire.iter().enumerate() {
                let curve = edge.curve();
                let Some(parameter) = curve.search_nearest_parameter(point, None, 10) else {
                    continue;
                };
                let projected = curve.evaluate(parameter);
                if projected.distance(point) > tolerance {
                    continue;
                }
                let Some(kind) = ParameterKind::try_new(parameter, curve.range_tuple()) else {
                    continue;
                };
                let endpoint_id = match kind {
                    ParameterKind::Front => Some(edge.absolute_front().id()),
                    ParameterKind::Back => Some(edge.absolute_back().id()),
                    ParameterKind::Inner(_) => None,
                };
                let canonical_point = match kind {
                    ParameterKind::Front => edge.absolute_front().point(),
                    ParameterKind::Back => edge.absolute_back().point(),
                    ParameterKind::Inner(_) => projected,
                };
                if let Some((existing_edge, existing_endpoint)) = candidate_identity {
                    let same_edge = existing_edge == edge.id();
                    let same_endpoint =
                        existing_endpoint.is_some() && existing_endpoint == endpoint_id;
                    if same_edge || same_endpoint {
                        continue;
                    }
                    return None;
                }
                candidate_identity = Some((edge.id(), endpoint_id));
                candidate = Some((wire_index, edge_index, kind, Some(canonical_point)));
            }
        }
        candidate
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
            self.search_boundary_parameter(loops_index, vertex.point())?;
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

    fn add_geom_vertex<S>(
        &mut self,
        (loops_index, wire_index, edge_index): (usize, usize, usize),
        vertex: &Vertex<Point3>,
        kind: ParameterKind,
        another_surface: &S,
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
        let (wire_index, edge_index, kind, projected) = if self
            .get(loops_index)
            .and_then(|loops| loops.get(wire_index))
            .and_then(|wire| wire.get(edge_index))
            .is_some()
        {
            (wire_index, edge_index, kind, None)
        } else {
            self.search_boundary_parameter(loops_index, vertex.point())?
        };
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

fn planar_polygon_plane(polygon: &PolygonMesh, tol: f64) -> Option<(Point3, Vector3)> {
    let tolerance = BooleanTolerance::new(tol);
    let positions = polygon.positions();
    let origin = *positions.first()?;
    let mut normal = None;
    'outer: for i in 1..positions.len() {
        for j in i + 1..positions.len() {
            let candidate = (positions[i] - origin).cross(positions[j] - origin);
            if candidate.magnitude() > tolerance.area_floor {
                normal = Some(candidate.normalize());
                break 'outer;
            }
        }
    }
    let normal = normal?;
    positions
        .iter()
        .all(|point| (*point - origin).dot(normal).abs() <= tolerance.length)
        .then_some((origin, normal))
}

fn planar_polygons_coplanar(lhs: &PolygonMesh, rhs: &PolygonMesh, tol: f64) -> bool {
    let Some((lhs_origin, lhs_normal)) = planar_polygon_plane(lhs, tol) else {
        return false;
    };
    let Some((rhs_origin, rhs_normal)) = planar_polygon_plane(rhs, tol) else {
        return false;
    };
    let tolerance = BooleanTolerance::new(tol);
    lhs_normal.cross(rhs_normal).magnitude() <= tolerance.angular
        && (rhs_origin - lhs_origin).dot(lhs_normal).abs() <= tolerance.length
}

#[derive(Clone, Copy, Debug)]
struct PlanarPolygonSignature {
    area: f64,
    centroid: Point3,
    bounds: [f64; 4],
}

fn planar_basis(normal: Vector3) -> Option<(Vector3, Vector3)> {
    let reference = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
        Vector3::unit_x()
    } else if normal.y.abs() <= normal.z.abs() {
        Vector3::unit_y()
    } else {
        Vector3::unit_z()
    };
    let u = normal.cross(reference);
    (!u.magnitude2().so_small()).then(|| {
        let u = u.normalize();
        (u, normal.cross(u).normalize())
    })
}

fn planar_polygon_set_signature<'a>(
    polygons: impl IntoIterator<Item = &'a PolygonMesh>,
    origin: Point3,
    normal: Vector3,
) -> Option<PlanarPolygonSignature> {
    let (u, v) = planar_basis(normal)?;
    let mut area = 0.0;
    let mut weighted_centroid = Vector3::new(0.0, 0.0, 0.0);
    let mut min_u = f64::INFINITY;
    let mut max_u = f64::NEG_INFINITY;
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    let mut any = false;

    for polygon in polygons {
        any = true;
        for triangle in polygon.faces().triangle_iter() {
            let a = polygon.positions()[triangle[0].pos];
            let b = polygon.positions()[triangle[1].pos];
            let c = polygon.positions()[triangle[2].pos];
            let triangle_area = 0.5 * (b - a).cross(c - a).magnitude();
            if !triangle_area.is_finite() || triangle_area.so_small() {
                continue;
            }
            area += triangle_area;
            weighted_centroid += (a.to_vec() + b.to_vec() + c.to_vec()) * (triangle_area / 3.0);
        }
        for point in polygon.positions() {
            let delta = *point - origin;
            let pu = delta.dot(u);
            let pv = delta.dot(v);
            min_u = min_u.min(pu);
            max_u = max_u.max(pu);
            min_v = min_v.min(pv);
            max_v = max_v.max(pv);
        }
    }
    if !any || !area.is_finite() || area.so_small() {
        return None;
    }
    let centroid = Point3::from_vec(weighted_centroid / area);
    [min_u, max_u, min_v, max_v]
        .iter()
        .all(|value| value.is_finite())
        .then_some(PlanarPolygonSignature {
            area,
            centroid,
            bounds: [min_u, max_u, min_v, max_v],
        })
}

fn planar_polygon_signature(
    polygon: &PolygonMesh,
    origin: Point3,
    normal: Vector3,
) -> Option<PlanarPolygonSignature> {
    planar_polygon_set_signature(std::iter::once(polygon), origin, normal)
}

fn planar_signatures_match(
    lhs: PlanarPolygonSignature,
    rhs: PlanarPolygonSignature,
    tol: f64,
) -> bool {
    let tolerance = BooleanTolerance::new(tol);
    let span = (lhs.bounds[1] - lhs.bounds[0])
        .abs()
        .max((lhs.bounds[3] - lhs.bounds[2]).abs())
        .max((rhs.bounds[1] - rhs.bounds[0]).abs())
        .max((rhs.bounds[3] - rhs.bounds[2]).abs());
    let area_tol = tolerance.area_for_span(span);
    (lhs.area - rhs.area).abs() <= area_tol
        && lhs.centroid.distance(rhs.centroid) <= 2.0 * tolerance.length
        && lhs
            .bounds
            .iter()
            .zip(rhs.bounds)
            .all(|(lhs, rhs)| (*lhs - rhs).abs() <= 2.0 * tolerance.length)
}

fn cross2(lhs: Vector2, rhs: Vector2) -> f64 { lhs.x * rhs.y - lhs.y * rhs.x }

fn projected_triangle(
    polygon: &PolygonMesh,
    triangle: [StandardVertex; 3],
    origin: Point3,
    u: Vector3,
    v: Vector3,
) -> [Point2; 3] {
    triangle.map(|vertex| {
        let delta = polygon.positions()[vertex.pos] - origin;
        Point2::new(delta.dot(u), delta.dot(v))
    })
}

fn polygon_area2(points: &[Point2]) -> f64 {
    if points.len() < 3 {
        return 0.0;
    }
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .take(points.len())
        .map(|(lhs, rhs)| lhs.x * rhs.y - lhs.y * rhs.x)
        .sum()
}

fn clip_convex_polygon(subject: &[Point2], clip: [Point2; 3], tol: f64) -> Vec<Point2> {
    let mut clip = clip;
    if polygon_area2(&clip) < 0.0 {
        clip.swap(1, 2);
    }
    let mut output = subject.to_vec();
    for edge_index in 0..3 {
        if output.is_empty() {
            break;
        }
        let a = clip[edge_index];
        let b = clip[(edge_index + 1) % 3];
        let edge = b - a;
        let cross_tol = tol * edge.magnitude().max(tol);
        let input = std::mem::take(&mut output);
        let Some(&mut_previous) = input.last() else {
            break;
        };
        let mut previous = mut_previous;
        let mut previous_distance = cross2(edge, previous - a);
        let mut previous_inside = previous_distance >= -cross_tol;
        for current in input {
            let current_distance = cross2(edge, current - a);
            let current_inside = current_distance >= -cross_tol;
            if current_inside != previous_inside {
                let denominator = previous_distance - current_distance;
                if denominator.abs() > f64::EPSILON {
                    let ratio = previous_distance / denominator;
                    output.push(previous + (current - previous) * ratio);
                }
            }
            if current_inside {
                output.push(current);
            }
            previous = current;
            previous_distance = current_distance;
            previous_inside = current_inside;
        }
    }
    output
}

fn triangle_intersection_area(lhs: [Point2; 3], rhs: [Point2; 3], tol: f64) -> f64 {
    let clipped = clip_convex_polygon(&lhs, rhs, tol);
    0.5 * polygon_area2(&clipped).abs()
}

fn planar_polygon_overlap_area(
    lhs: &PolygonMesh,
    rhs: &PolygonMesh,
    origin: Point3,
    normal: Vector3,
    tol: f64,
) -> Option<f64> {
    let (u, v) = planar_basis(normal)?;
    let lhs_triangles = lhs
        .faces()
        .triangle_iter()
        .map(|triangle| projected_triangle(lhs, triangle, origin, u, v))
        .collect::<Vec<_>>();
    let rhs_triangles = rhs
        .faces()
        .triangle_iter()
        .map(|triangle| projected_triangle(rhs, triangle, origin, u, v))
        .collect::<Vec<_>>();
    Some(
        lhs_triangles
            .iter()
            .flat_map(|lhs| {
                rhs_triangles
                    .iter()
                    .map(move |rhs| triangle_intersection_area(*lhs, *rhs, tol))
            })
            .sum(),
    )
}

fn planar_area_tolerance(signature: PlanarPolygonSignature, tol: f64) -> f64 {
    let span = (signature.bounds[1] - signature.bounds[0])
        .abs()
        .max((signature.bounds[3] - signature.bounds[2]).abs());
    BooleanTolerance::new(tol).area_for_span(span)
}

fn planar_polygon_set_has_interior_overlap(
    polygons: &[PolygonMesh],
    origin: Point3,
    normal: Vector3,
    tol: f64,
) -> Option<bool> {
    for left in 0..polygons.len() {
        let left_signature = planar_polygon_signature(&polygons[left], origin, normal)?;
        for right in left + 1..polygons.len() {
            let right_signature = planar_polygon_signature(&polygons[right], origin, normal)?;
            if !projected_bounds_overlap(left_signature.bounds, right_signature.bounds, tol) {
                continue;
            }
            let overlap = planar_polygon_overlap_area(
                &polygons[left],
                &polygons[right],
                origin,
                normal,
                tol,
            )?;
            let area_tol = planar_area_tolerance(left_signature, tol)
                .max(planar_area_tolerance(right_signature, tol));
            if overlap > area_tol {
                return Some(true);
            }
        }
    }
    Some(false)
}

fn planar_polygons_same_trimmed_face(lhs: &PolygonMesh, rhs: &PolygonMesh, tol: f64) -> bool {
    let Some((lhs_origin, lhs_normal)) = planar_polygon_plane(lhs, tol) else {
        return false;
    };
    let Some((rhs_origin, rhs_normal)) = planar_polygon_plane(rhs, tol) else {
        return false;
    };
    let tolerance = BooleanTolerance::new(tol);
    if lhs_normal.cross(rhs_normal).magnitude() > tolerance.angular
        || (rhs_origin - lhs_origin).dot(lhs_normal).abs() > tolerance.length
    {
        return false;
    }
    let Some(lhs_signature) = planar_polygon_signature(lhs, lhs_origin, lhs_normal) else {
        return false;
    };
    let Some(rhs_signature) = planar_polygon_signature(rhs, lhs_origin, lhs_normal) else {
        return false;
    };
    if !planar_signatures_match(lhs_signature, rhs_signature, tol) {
        return false;
    }
    let Some(overlap) = planar_polygon_overlap_area(lhs, rhs, lhs_origin, lhs_normal, tol) else {
        return false;
    };
    let area_tol =
        planar_area_tolerance(lhs_signature, tol).max(planar_area_tolerance(rhs_signature, tol));
    (lhs_signature.area - overlap).abs() <= area_tol
        && (rhs_signature.area - overlap).abs() <= area_tol
}

fn projected_bounds_overlap(lhs: [f64; 4], rhs: [f64; 4], tol: f64) -> bool {
    lhs[1].min(rhs[1]) - lhs[0].max(rhs[0]) > tol && lhs[3].min(rhs[3]) - lhs[2].max(rhs[2]) > tol
}

/// Match one trimmed planar face against a partition of the same physical
/// interface on the other shell. Revolving a radial cap, for example, may
/// represent one disk as two semicircular faces while a sweep uses one disk.
pub(super) fn single_face_matches_planar_partition(
    target: &PolygonMesh,
    fragments: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<Option<usize>> {
    let Some((origin, normal)) = planar_polygon_plane(target, tol) else {
        return Some(None);
    };
    let Some(target_signature) = planar_polygon_signature(target, origin, normal) else {
        return Some(None);
    };
    let mut matched = Vec::<PolygonMesh>::new();
    let mut representative = None;
    for (index, face) in fragments.iter().enumerate() {
        let polygon = face.surface()?;
        if !planar_polygons_coplanar(target, &polygon, tol) {
            continue;
        }
        let signature = planar_polygon_signature(&polygon, origin, normal)?;
        if !projected_bounds_overlap(target_signature.bounds, signature.bounds, tol) {
            continue;
        }
        matched.push(polygon);
        representative.get_or_insert(index);
    }
    if matched.is_empty() {
        return Some(None);
    }
    let fragment_signature = planar_polygon_set_signature(matched.iter(), origin, normal)?;
    if !planar_signatures_match(target_signature, fragment_signature, tol)
        || planar_polygon_set_has_interior_overlap(&matched, origin, normal, tol)?
    {
        return Some(None);
    }
    let overlap = matched.iter().try_fold(0.0, |sum, fragment| {
        planar_polygon_overlap_area(target, fragment, origin, normal, tol).map(|area| sum + area)
    })?;
    let area_tol = planar_area_tolerance(target_signature, tol)
        .max(planar_area_tolerance(fragment_signature, tol));
    let Some(representative) = representative else {
        return Some(None);
    };
    Some(
        ((target_signature.area - overlap).abs() <= area_tol
            && (fragment_signature.area - overlap).abs() <= area_tol)
            .then_some(representative),
    )
}

fn shells_full_planar_interface_pair(
    lhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    rhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<Option<(usize, usize)>> {
    for (left_index, left) in lhs.iter().enumerate() {
        let left_polygon = left.surface()?;
        if let Some(right_index) = single_face_matches_planar_partition(&left_polygon, rhs, tol)? {
            return Some(Some((left_index, right_index)));
        }
    }
    for (right_index, right) in rhs.iter().enumerate() {
        let right_polygon = right.surface()?;
        if let Some(left_index) = single_face_matches_planar_partition(&right_polygon, lhs, tol)? {
            return Some(Some((left_index, right_index)));
        }
    }
    Some(None)
}

fn planar_polygons_partially_overlap(lhs: &PolygonMesh, rhs: &PolygonMesh, tol: f64) -> bool {
    if !planar_polygons_coplanar(lhs, rhs, tol) || planar_polygons_same_trimmed_face(lhs, rhs, tol)
    {
        return false;
    }
    lhs.collide_with_neighborhood_of(rhs.positions(), tol)
        || rhs.collide_with_neighborhood_of(lhs.positions(), tol)
}

fn shell_polygon_bounds(
    shell: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
) -> Option<(Point3, Point3)> {
    let mut points = shell
        .iter()
        .filter_map(|face| face.surface())
        .flat_map(|polygon| polygon.positions().clone());
    let first = points.next()?;
    Some(points.fold((first, first), |(mut min, mut max), point| {
        min.x = min.x.min(point.x);
        min.y = min.y.min(point.y);
        min.z = min.z.min(point.z);
        max.x = max.x.max(point.x);
        max.y = max.y.max(point.y);
        max.z = max.z.max(point.z);
        (min, max)
    }))
}

fn shells_have_positive_aabb_overlap(
    lhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    rhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<bool> {
    let (lhs_min, lhs_max) = shell_polygon_bounds(lhs)?;
    let (rhs_min, rhs_max) = shell_polygon_bounds(rhs)?;
    Some(
        lhs_max.x.min(rhs_max.x) - lhs_min.x.max(rhs_min.x) > tol
            && lhs_max.y.min(rhs_max.y) - lhs_min.y.max(rhs_min.y) > tol
            && lhs_max.z.min(rhs_max.z) - lhs_min.z.max(rhs_min.z) > tol,
    )
}

fn shells_have_full_planar_interface(
    lhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    rhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<bool> {
    Some(shells_full_planar_interface_pair(lhs, rhs, tol)?.is_some())
}

fn shells_have_partial_planar_overlap(
    lhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    rhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<bool> {
    if !shells_have_positive_aabb_overlap(lhs, rhs, tol)? {
        return Some(false);
    }
    for left in lhs.iter() {
        let left_polygon = left.surface()?;
        for right in rhs.iter() {
            let right_polygon = right.surface()?;
            if planar_polygons_partially_overlap(&left_polygon, &right_polygon, tol) {
                return Some(true);
            }
        }
    }
    Some(false)
}

fn imprint_edges_on_faces<C>(
    source: &LoopsStore<Point3, C>,
    target_poly_shell: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    target: &mut LoopsStore<Point3, C>,
    tol: f64,
) -> Option<()>
where
    C: ParametricCurve3D
        + BoundedCurve<Point = Point3>
        + Cut<Point = Point3, Vector = Vector3>
        + SearchParameter<CurveParameter, Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3>
        + Invertible
        + Clone,
{
    let mut unique = Vec::<Edge<Point3, C>>::new();
    source
        .iter()
        .flat_map(|loops| loops.iter())
        .flat_map(|wire| wire.iter())
        .for_each(|edge| {
            if !unique.iter().any(|candidate| candidate.id() == edge.id()) {
                unique.push(edge.clone());
            }
        });

    let proximity = 2.0 * BooleanTolerance::new(tol).length;
    for face_index in 0..target.len() {
        let polygon = target_poly_shell[face_index].surface()?;
        for edge in &unique {
            let curve = edge.oriented_curve();
            let (t0, t1) = curve.range_tuple();
            let samples = [
                edge.front().point(),
                curve.subs((t0 + t1) * 0.5),
                edge.back().point(),
            ];
            if !polygon.neighborhood_include(&samples, proximity) {
                continue;
            }
            if target[face_index]
                .search_parameter(edge.front().point())
                .is_none()
                || target[face_index]
                    .search_parameter(edge.back().point())
                    .is_none()
            {
                continue;
            }

            let mut emap = HashMap::default();
            target.add_polygon_vertex(face_index, edge.front(), &mut emap)?;
            target.add_polygon_vertex(face_index, edge.back(), &mut emap)?;
            if target[face_index].has_partitioned_coincident_path(edge) {
                continue;
            }
            target[face_index].add_edge(edge.clone(), ShapesOpStatus::And, true);
        }
    }
    Some(())
}

const CROSS_FACE_STITCH_TOLERANCE: f64 = 8.0 * TOLERANCE;

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
    let edge0 = Edge::new(&v0, &v1, poly_curve0.into());
    let edge1 = Edge::new(&v1, &v0, poly_curve1.into());
    wire![edge0, edge1]
}

#[allow(dead_code)]
pub(super) struct LoopsStoreQuadruple<C> {
    pub(super) geom_loops_store0: LoopsStore<Point3, C>,
    pub(super) poly_loops_store0: LoopsStore<Point3, PolylineCurve>,
    pub(super) geom_loops_store1: LoopsStore<Point3, C>,
    pub(super) poly_loops_store1: LoopsStore<Point3, PolylineCurve>,
    pub(super) coplanar_overlap: bool,
}

pub(super) fn create_loops_stores<C, S>(
    geom_shell0: &Shell<Point3, C, S>,
    poly_shell0: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    geom_shell1: &Shell<Point3, C, S>,
    poly_shell1: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
    imprint_coplanar: bool,
) -> Option<LoopsStoreQuadruple<C>>
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
    let mut geom_loops_store0: LoopsStore<_, _> = geom_shell0.face_iter().collect();
    let mut poly_loops_store0: LoopsStore<_, _> = poly_shell0.face_iter().collect();
    let mut geom_loops_store1: LoopsStore<_, _> = geom_shell1.face_iter().collect();
    let mut poly_loops_store1: LoopsStore<_, _> = poly_shell1.face_iter().collect();
    let store0_len = geom_loops_store0.len();
    let store1_len = geom_loops_store1.len();
    let coplanar_tol = 2.0 * BooleanTolerance::new(tol).length;
    let positive_aabb_overlap =
        shells_have_positive_aabb_overlap(poly_shell0, poly_shell1, coplanar_tol)?;
    let full_coplanar_interface = imprint_coplanar
        && shells_have_full_planar_interface(poly_shell0, poly_shell1, coplanar_tol)?;
    let full_interface_adjacency = full_coplanar_interface && !positive_aabb_overlap;
    let partial_coplanar_overlap = imprint_coplanar
        && positive_aabb_overlap
        && shells_have_partial_planar_overlap(poly_shell0, poly_shell1, coplanar_tol)?;
    let coplanar_overlap = full_coplanar_interface || partial_coplanar_overlap;
    (0..store0_len)
        .flat_map(move |i| (0..store1_len).map(move |j| (i, j)))
        .try_for_each(|(face_index0, face_index1)| {
            // Solids that meet on one complete trimmed planar face with
            // opposite outward normals need no SSI at all. Generic SSI on
            // their coincident continuation surfaces (e.g. coaxial cylinders)
            // is degenerate; downstream coplanar classification removes the
            // shared internal face and topology welding stitches the rim.
            if full_interface_adjacency {
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
                        poly_loops_store0[face_index0]
                            .add_independent_loop(BoundaryWire::new(poly_wire.clone(), status0));
                        poly_loops_store1[face_index1]
                            .add_independent_loop(BoundaryWire::new(poly_wire, status1));
                        let geom_wire = create_independent_loop(intersection_curve);
                        geom_loops_store0[face_index0]
                            .add_independent_loop(BoundaryWire::new(geom_wire.clone(), status0));
                        geom_loops_store1[face_index1]
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
                            &geom_loops_store0[face_index0],
                            &geom_loops_store1[face_index1],
                            front,
                            stitch_tolerance,
                        ) {
                            front = point;
                        }
                        if let Some(point) = unique_nearby_loop_point(
                            &geom_loops_store0[face_index0],
                            &geom_loops_store1[face_index1],
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
                        *polyline.first_mut().unwrap() = front;
                        *polyline.last_mut().unwrap() = back;
                        *intersection_curve.leader_mut().first_mut().unwrap() = front;
                        *intersection_curve.leader_mut().last_mut().unwrap() = back;

                        let pv0 = Vertex::new(front);
                        let pv1 = Vertex::new(back);
                        let gv0 = Vertex::new(front);
                        let gv1 = Vertex::new(back);
                        let mut pemap0 = HashMap::default();
                        let mut pemap1 = HashMap::default();
                        let mut gemap0 = HashMap::default();
                        let mut gemap1 = HashMap::default();
                        let idx00 =
                            poly_loops_store0.add_polygon_vertex(face_index0, &pv0, &mut pemap0);
                        if let Some((wire_index, edge_index, kind)) = idx00 {
                            geom_loops_store0.add_geom_vertex(
                                (face_index0, wire_index, edge_index),
                                &gv0,
                                kind,
                                &surface1,
                                &mut gemap0,
                            )?;
                            let polyline = intersection_curve.leader_mut();
                            *polyline.first_mut().unwrap() = gv0.point();
                        }
                        let idx01 =
                            poly_loops_store0.add_polygon_vertex(face_index0, &pv1, &mut pemap1);
                        if let Some((wire_index, edge_index, kind)) = idx01 {
                            geom_loops_store0.add_geom_vertex(
                                (face_index0, wire_index, edge_index),
                                &gv1,
                                kind,
                                &surface1,
                                &mut gemap1,
                            )?;
                            let polyline = intersection_curve.leader_mut();
                            *polyline.last_mut().unwrap() = gv1.point();
                        }
                        let idx10 =
                            poly_loops_store1.add_polygon_vertex(face_index1, &pv0, &mut pemap0);
                        if let Some((wire_index, edge_index, kind)) = idx10 {
                            geom_loops_store1.add_geom_vertex(
                                (face_index1, wire_index, edge_index),
                                &gv0,
                                kind,
                                &surface0,
                                &mut gemap0,
                            )?;
                            let polyline = intersection_curve.leader_mut();
                            *polyline.first_mut().unwrap() = gv0.point();
                        }
                        let idx11 =
                            poly_loops_store1.add_polygon_vertex(face_index1, &pv1, &mut pemap1);
                        if let Some((wire_index, edge_index, kind)) = idx11 {
                            geom_loops_store1.add_geom_vertex(
                                (face_index1, wire_index, edge_index),
                                &gv1,
                                kind,
                                &surface0,
                                &mut gemap1,
                            )?;
                            let polyline = intersection_curve.leader_mut();
                            *polyline.last_mut().unwrap() = gv1.point();
                        }
                        let pedge = Edge::new(&pv0, &pv1, polyline);
                        let gedge = Edge::new(&gv0, &gv1, intersection_curve.into());
                        poly_loops_store0[face_index0].add_edge(
                            pedge.clone(),
                            status0,
                            partial_coplanar_overlap,
                        );
                        geom_loops_store0[face_index0].add_edge(
                            gedge.clone(),
                            status0,
                            partial_coplanar_overlap,
                        );
                        poly_loops_store1[face_index1].add_edge(
                            pedge,
                            status1,
                            partial_coplanar_overlap,
                        );
                        geom_loops_store1[face_index1].add_edge(
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
        let source0 = geom_loops_store0.clone();
        let source1 = geom_loops_store1.clone();
        imprint_edges_on_faces(&source1, poly_shell0, &mut geom_loops_store0, tol)?;
        imprint_edges_on_faces(&source0, poly_shell1, &mut geom_loops_store1, tol)?;
    }

    Some(LoopsStoreQuadruple {
        geom_loops_store0,
        poly_loops_store0,
        geom_loops_store1,
        poly_loops_store1,
        coplanar_overlap,
    })
}

#[cfg(test)]
mod projection_tests {
    use super::*;

    #[test]
    fn change_vertex_removes_edge_collapsed_by_coalescing() {
        let v0 = Vertex::new(0usize);
        let v1 = Vertex::new(0usize);
        let v2 = Vertex::new(1usize);
        let edge01 = Edge::new(&v0, &v1, ());
        let edge12 = Edge::new(&v1, &v2, ());
        let edge20 = Edge::new(&v2, &v0, ());
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
        let edge = Edge::new(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * TOLERANCE, 0.0);

        let (_, _, kind, projected) = store
            .search_boundary_parameter(0, point)
            .expect("unique bounded nearest projection should recover the boundary parameter");

        assert!(matches!(kind, ParameterKind::Inner(_)));
        assert_eq!(projected, Some(Point3::new(0.5, 0.0, 0.0)));
    }

    #[test]
    fn loop_boundary_vertex_insertion_snaps_recovered_projection_before_cut() {
        let v0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let v1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new(&v0, &v1, Line(v0.point(), v1.point()));
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
    fn exact_geometry_prefers_existing_polygon_boundary_hint() {
        let first_front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let first_back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let hinted_front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let hinted_back = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let first_edge = Edge::new(
            &first_front,
            &first_back,
            Line(first_front.point(), first_back.point()),
        );
        let hinted_edge = Edge::new(
            &hinted_front,
            &hinted_back,
            Line(hinted_front.point(), hinted_back.point()),
        );
        let loops = Loops(vec![
            BoundaryWire::new(vec![first_edge].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![hinted_edge].into(), ShapesOpStatus::Unknown),
        ]);
        let mut store = LoopsStore(vec![loops]);
        let recovered = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let mut emap = HashMap::default();

        store
            .add_geom_vertex(
                (0, 1, 0),
                &recovered,
                ParameterKind::Front,
                &Plane::xy(),
                &mut emap,
            )
            .expect("existing polygon hint should disambiguate coincident exact boundaries");

        assert!(
            store[0][0]
                .vertex_iter()
                .any(|vertex| vertex == first_front)
        );
        assert!(store[0][1].vertex_iter().any(|vertex| vertex == recovered));
        assert!(
            store[0][1]
                .vertex_iter()
                .all(|vertex| vertex != hinted_front)
        );
    }

    #[test]
    fn exact_geometry_falls_back_when_polygon_hint_is_stale() {
        let front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let old_back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new(&front, &old_back, Line(front.point(), old_back.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let mut store = LoopsStore(vec![loops]);
        let recovered = Vertex::new(old_back.point());
        let mut emap = HashMap::default();

        store
            .add_geom_vertex(
                (0, 99, 99),
                &recovered,
                ParameterKind::Inner(0.999),
                &Plane::xy(),
                &mut emap,
            )
            .expect("stale polygon hint should fall back to exact boundary search");

        assert_eq!(store[0][0].len(), 1);
        assert!(store[0][0].vertex_iter().any(|vertex| vertex == recovered));
        assert!(store[0][0].vertex_iter().all(|vertex| vertex != old_back));
    }

    #[test]
    fn exact_geometry_reclassifies_inner_hint_at_boundary_endpoint() {
        let front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let old_back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new(&front, &old_back, Line(front.point(), old_back.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let mut store = LoopsStore(vec![loops]);
        let recovered = Vertex::new(old_back.point());
        let mut emap = HashMap::default();

        store
            .add_geom_vertex(
                (0, 0, 0),
                &recovered,
                ParameterKind::Inner(0.999),
                &Plane::xy(),
                &mut emap,
            )
            .expect("exact endpoint must be reclassified before topology mutation");

        assert_eq!(store[0][0].len(), 1);
        assert!(store[0][0].vertex_iter().any(|vertex| vertex == recovered));
        assert!(store[0][0].vertex_iter().all(|vertex| vertex != old_back));
    }

    #[test]
    fn loop_boundary_search_deduplicates_inverse_use_of_same_edge() {
        let v0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let v1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new(&v0, &v1, Line(v0.point(), v1.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge.clone(), edge.inverse()].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * TOLERANCE, 0.0);

        assert!(store.search_boundary_parameter(0, point).is_some());
    }

    #[test]
    fn loop_boundary_search_accepts_distinct_edges_at_one_shared_vertex() {
        let left = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let shared = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let right = Vertex::new(Point3::new(2.0, 0.0, 0.0));
        let edge_left = Edge::new(&left, &shared, Line(left.point(), shared.point()));
        let edge_right = Edge::new(&shared, &right, Line(shared.point(), right.point()));
        let loops = Loops(vec![BoundaryWire::new(
            vec![edge_left, edge_right].into(),
            ShapesOpStatus::Unknown,
        )]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(1.0, 1.125 * TOLERANCE, 0.0);

        let (_, _, kind, projected) = store
            .search_boundary_parameter(0, point)
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
        let edge_a = Edge::new(&a0, &a1, Line(a0.point(), a1.point()));
        let edge_b = Edge::new(&b0, &b1, Line(b0.point(), b1.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![edge_a].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge_b].into(), ShapesOpStatus::Unknown),
        ]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * TOLERANCE, 0.0);

        assert!(store.search_boundary_parameter(0, point).is_none());
    }

    #[test]
    fn cross_face_endpoint_stitch_rejects_model_scale_nearby_topology() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new(&a, &b, Line(a.point(), b.point()));
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
            Edge::new(front, back, Line(front.point(), back.point()))
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

        let result = loops.add_edge(incoming, ShapesOpStatus::And, false);

        assert_eq!(result[0].map(|x| x.0), Some(1));
        assert_eq!(result[1].map(|x| x.0), Some(1));
        assert_eq!(loops.len(), 3);
        assert!(
            loops
                .iter()
                .all(|wire| wire.is_closed() && wire.is_simple())
        );
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
                Edge::new(&a, &b, Line(a.point(), b.point())),
                Edge::new(&b, &c, Line(b.point(), c.point())),
                Edge::new(&c, &d, Line(c.point(), d.point())),
                Edge::new(&d, &e, Line(d.point(), e.point())),
                Edge::new(&e, &f, Line(e.point(), f.point())),
                Edge::new(&f, &a, Line(f.point(), a.point())),
            ]
            .into(),
            ShapesOpStatus::Or,
        );
        let loops = Loops(vec![wire]);
        let incoming = Edge::new(&a, &d, Line(a.point(), d.point()));

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
                Edge::new(&a, &d, Line(a.point(), d.point())),
                Edge::new(&d, &e, Line(d.point(), e.point())),
                Edge::new(&e, &f, Line(e.point(), f.point())),
                Edge::new(&f, &a, Line(f.point(), a.point())),
            ]
            .into(),
            ShapesOpStatus::Or,
        );
        let loops = Loops(vec![wire]);
        let incoming = Edge::new(&a, &d, Line(a.point(), d.point()));

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

#[cfg(test)]
mod planar_region_identity_tests {
    use super::*;

    fn diamond() -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(2.0, 1.0, 0.0),
                    Point3::new(1.0, 2.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
        )
    }

    fn disjoint_corner_triangles() -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(2.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(2.0, 2.0, 0.0),
                    Point3::new(0.0, 2.0, 0.0),
                    Point3::new(2.0, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [3, 4, 5]]),
        )
    }

    #[test]
    fn area_centroid_and_bounds_are_not_region_identity() {
        assert!(!planar_polygons_same_trimmed_face(
            &diamond(),
            &disjoint_corner_triangles(),
            1.0e-8
        ));
    }

    #[test]
    fn duplicate_overlapping_fragments_do_not_fake_partition_coverage() {
        let target = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(2.0, 0.0, 0.0),
                    Point3::new(2.0, 2.0, 0.0),
                    Point3::new(0.0, 2.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
        );
        let diamond = diamond();
        let fragments: Shell<Point3, PolylineCurve, Option<PolygonMesh>> = Shell::from(vec![
            Face::new_unchecked(Vec::new(), Some(diamond.clone())),
            Face::new_unchecked(Vec::new(), Some(diamond)),
        ]);

        assert_eq!(
            single_face_matches_planar_partition(&target, &fragments, 1.0e-8),
            Some(None)
        );
    }
}
