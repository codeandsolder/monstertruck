//! Classic (0.3.2) boolean loops store.
//!
//! Ported verbatim from the published 0.3.2 crate's `transversal::loops_store`,
//! adjusted only for the graft: the marching intersection-curve backend now
//! yields a raw [`IntersectionCurve`] leader (the 0.3.2 parameter wrapper was
//! dropped), so the `.into()` that unwrapped the wrapper is gone.

#![allow(clippy::many_single_char_names)]

use super::intersection_curve;
use monstertruck_core::cgmath64::*;
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
    pub(super) fn new(wire: Wire<P, C>, status: ShapesOpStatus) -> Self {
        Self { wire, status }
    }
    #[inline(always)]
    pub(super) fn status(&self) -> ShapesOpStatus {
        self.status
    }
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
        S1: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
    {
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
    fn deref(&self) -> &Self::Target {
        &self.wire
    }
}

impl<P, C> std::ops::DerefMut for BoundaryWire<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.wire
    }
}

#[derive(Clone, Debug)]
pub(super) struct Loops<P, C>(Vec<BoundaryWire<P, C>>);
#[derive(Clone, Debug)]
pub(super) struct LoopsStore<P, C>(Vec<Loops<P, C>>);

impl<P, C> std::ops::Deref for Loops<P, C> {
    type Target = Vec<BoundaryWire<P, C>>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<P, C> std::ops::DerefMut for Loops<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<P, C> std::ops::Deref for LoopsStore<P, C> {
    type Target = Vec<Loops<P, C>>;
    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<P, C> std::ops::DerefMut for LoopsStore<P, C> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
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
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
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
    where
        C: BoundedCurve<Point = P> + SearchParameter<CurveParameter, Point = P>,
    {
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
        self.iter_mut()
            .flat_map(|wire| wire.iter_mut())
            .for_each(|edge| {
                let mut new_edge = if edge.absolute_front() == old_vertex {
                    emap.entry(edge.id()).or_insert_with(|| {
                        Edge::new(new_vertex, edge.absolute_back(), edge.curve())
                    })
                } else if edge.absolute_back() == old_vertex {
                    emap.entry(edge.id()).or_insert_with(|| {
                        Edge::new(edge.absolute_front(), new_vertex, edge.curve())
                    })
                } else {
                    return;
                }
                .clone();
                if !edge.orientation() {
                    new_edge.invert();
                }
                // Remove the edge from the HashMap when it is no longer there because Id reassignment will occur.
                if edge.count() == 1 {
                    emap.remove(&edge.id());
                }
                *edge = new_edge;
            })
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
            let curves_match = |edge: &Edge<P, C>| {
                let lhs = edge.oriented_curve();
                let rhs = edge0.oriented_curve();
                let sample_on = |curve: &C, other: &C| {
                    let (t0, t1) = curve.range_tuple();
                    [0.25, 0.5, 0.75].into_iter().all(|fraction| {
                        let point = curve.evaluate(t0 + (t1 - t0) * fraction);
                        other
                            .search_nearest_parameter(point, None, 100)
                            .map(|parameter| other.evaluate(parameter).near(&point))
                            .unwrap_or(false)
                    })
                };
                let midpoint = |curve: &C| {
                    let (t0, t1) = curve.range_tuple();
                    curve.evaluate((t0 + t1) * 0.5)
                };
                midpoint(&lhs).near(&midpoint(&rhs))
                    || (sample_on(&lhs, &rhs) && sample_on(&rhs, &lhs))
            };
            let coincident = self.iter().flat_map(|wire| wire.iter()).find_map(|edge| {
                let same_direction = edge.front().point().near(&edge0.front().point())
                    && edge.back().point().near(&edge0.back().point());
                let opposite_direction = edge.front().point().near(&edge0.back().point())
                    && edge.back().point().near(&edge0.front().point());
                (same_direction || opposite_direction)
                    .then(|| curves_match(edge))
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

        let a = self.iter().enumerate().find_map(|(i, wire)| {
            wire.iter().enumerate().find_map(|(j, edge)| {
                if edge.front() == edge0.back() {
                    Some((i, j))
                } else {
                    None
                }
            })
        });
        let b = self.iter().enumerate().find_map(|(i, wire)| {
            wire.iter().enumerate().find_map(|(j, edge)| {
                if edge.front() == edge0.front() {
                    Some((i, j))
                } else {
                    None
                }
            })
        });
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

    fn add_polygon_vertex(
        &mut self,
        loops_index: usize,
        v: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) -> Option<(usize, usize, ParameterKind)>
    where
        C: Cut<Point = P> + SearchParameter<CurveParameter, Point = P>,
    {
        let pt = v.point();
        let (wire_index, edge_index, kind) = self[loops_index].search_parameter(pt)?;
        match kind {
            ParameterKind::Front => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_front()
                    .clone();
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Back => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_back()
                    .clone();
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Inner(t) => {
                let edge = self[loops_index][wire_index][edge_index].absolute_clone();
                let edge_id = edge.id();
                let (edge0, edge1) = edge.cut_with_parameter(v, t)?;
                let new_wire: Wire<_, _> = vec![edge0, edge1].into();
                self.swap_edge_into_wire(edge_id, &new_wire);
            }
        }
        Some((wire_index, edge_index, kind))
    }
}

impl<C> LoopsStore<Point3, C> {
    fn add_geom_vertex<S>(
        &mut self,
        (loops_index, wire_index, edge_index): (usize, usize, usize),
        v: &Vertex<Point3>,
        kind: ParameterKind,
        another_surface: &S,
        emap: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
    ) -> Option<()>
    where
        C: Cut<Point = Point3, Vector = Vector3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
        S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
    {
        match kind {
            ParameterKind::Front => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_front()
                    .clone();
                v.set_point(old_vertex.point());
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Back => {
                let old_vertex = self[loops_index][wire_index][edge_index]
                    .absolute_back()
                    .clone();
                v.set_point(old_vertex.point());
                self.change_vertex(&old_vertex, v, emap);
            }
            ParameterKind::Inner(_) => {
                let curve = self[loops_index][wire_index][edge_index].curve();
                let (pt, t, _) =
                    curve_surface_projection(&curve, None, another_surface, None, v.point(), 100)?;
                v.set_point(pt);
                let edge = self[loops_index][wire_index][edge_index].absolute_clone();
                let edge_id = edge.id();
                let (edge0, edge1) = edge.cut_with_parameter(v, t)?;
                let new_wire: Wire<_, _> = vec![edge0, edge1].into();
                self.swap_edge_into_wire(edge_id, &new_wire);
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
    C: ParametricCurve3D + SearchNearestParameter<CurveParameter, Point = Point3>,
    S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if trials == 0 {
        return None;
    }
    let t = curve.search_nearest_parameter(point, curve_hint, 10)?;
    let pt0 = curve.evaluate(t);
    let (u, v) = surface.search_nearest_parameter(point, surface_hint, 10)?;
    let pt1 = surface.evaluate(u, v);
    if point.near(&pt0) && point.near(&pt1) && pt0.near(&pt1) {
        Some((point, t, Point2::new(u, v)))
    } else {
        let l = curve.derivative(t);
        let n = surface.normal(u, v);
        let t0 = (pt1 - pt0).dot(n) / l.dot(n);
        curve_surface_projection(
            curve,
            Some(t),
            surface,
            Some((u, v)),
            pt0 + t0 * l,
            trials - 1,
        )
    }
}

fn planar_polygon_plane(polygon: &PolygonMesh, tol: f64) -> Option<(Point3, Vector3)> {
    let positions = polygon.positions();
    let origin = *positions.first()?;
    let mut normal = None;
    'outer: for i in 1..positions.len() {
        for j in i + 1..positions.len() {
            let candidate = (positions[i] - origin).cross(positions[j] - origin);
            if candidate.magnitude() > tol {
                normal = Some(candidate.normalize());
                break 'outer;
            }
        }
    }
    let normal = normal?;
    positions
        .iter()
        .all(|point| (*point - origin).dot(normal).abs() <= tol)
        .then_some((origin, normal))
}

fn planar_polygons_coplanar(lhs: &PolygonMesh, rhs: &PolygonMesh, tol: f64) -> bool {
    let Some((lhs_origin, lhs_normal)) = planar_polygon_plane(lhs, tol) else {
        return false;
    };
    let Some((rhs_origin, rhs_normal)) = planar_polygon_plane(rhs, tol) else {
        return false;
    };
    lhs_normal.cross(rhs_normal).magnitude() <= tol
        && (rhs_origin - lhs_origin).dot(lhs_normal).abs() <= tol
}

fn planar_polygons_same_trimmed_face(lhs: &PolygonMesh, rhs: &PolygonMesh, tol: f64) -> bool {
    planar_polygons_coplanar(lhs, rhs, tol)
        && lhs.neighborhood_include(rhs.positions(), tol)
        && rhs.neighborhood_include(lhs.positions(), tol)
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
    for left in lhs.iter() {
        let left_polygon = left.surface()?;
        for right in rhs.iter() {
            let right_polygon = right.surface()?;
            if planar_polygons_same_trimmed_face(&left_polygon, &right_polygon, tol) {
                return Some(true);
            }
        }
    }
    Some(false)
}

fn shells_have_opposed_full_planar_interface<C, S>(
    geom_lhs: &Shell<Point3, C, S>,
    poly_lhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    geom_rhs: &Shell<Point3, C, S>,
    poly_rhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<bool>
where
    C: Clone,
    S: ParametricSurface3D
        + Clone
        + Invertible
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    for (left_index, left) in poly_lhs.iter().enumerate() {
        let left_polygon = left.surface()?;
        for (right_index, right) in poly_rhs.iter().enumerate() {
            let right_polygon = right.surface()?;
            if !planar_polygons_same_trimmed_face(&left_polygon, &right_polygon, tol) {
                continue;
            }
            let point = *left_polygon.positions().first()?;
            let left_surface = geom_lhs[left_index].oriented_surface();
            let right_surface = geom_rhs[right_index].oriented_surface();
            let left_parameter = left_surface
                .search_parameter(point, None, 100)
                .or_else(|| left_surface.search_nearest_parameter(point, None, 100))?;
            let right_parameter = right_surface
                .search_parameter(point, None, 100)
                .or_else(|| right_surface.search_nearest_parameter(point, None, 100))?;
            let left_normal = left_surface.normal(left_parameter.0, left_parameter.1);
            let right_normal = right_surface.normal(right_parameter.0, right_parameter.1);
            let left_norm = left_normal.magnitude();
            let right_norm = right_normal.magnitude();
            if left_norm <= tol || right_norm <= tol {
                continue;
            }
            let cosine = left_normal.dot(right_normal) / (left_norm * right_norm);
            if cosine <= -1.0 + 1.0e-6 {
                return Some(true);
            }
        }
    }
    Some(false)
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

    let proximity = tol.max(1.0e-4) * 2.0;
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
            target[face_index].add_edge(edge.clone(), ShapesOpStatus::And, true);
        }
    }
    Some(())
}

fn create_independent_loop<P, C, D>(mut poly_curve0: C) -> Wire<P, D>
where
    C: Cut<Point = P>,
    D: From<C>,
{
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
    let coplanar_tol = tol.max(1.0e-4) * 2.0;
    let positive_aabb_overlap =
        shells_have_positive_aabb_overlap(poly_shell0, poly_shell1, coplanar_tol)?;
    let full_coplanar_interface = imprint_coplanar
        && shells_have_full_planar_interface(poly_shell0, poly_shell1, coplanar_tol)?;
    let full_interface_adjacency = full_coplanar_interface
        && shells_have_opposed_full_planar_interface(
            geom_shell0,
            poly_shell0,
            geom_shell1,
            poly_shell1,
            coplanar_tol,
        )?;
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
            )?;
            curves.into_iter()
            .try_for_each(|(polyline, mut intersection_curve)| {
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
                    let pv0 = Vertex::new(polyline.front());
                    let pv1 = Vertex::new(polyline.back());
                    let gv0 = Vertex::new(polyline.front());
                    let gv1 = Vertex::new(polyline.back());
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
