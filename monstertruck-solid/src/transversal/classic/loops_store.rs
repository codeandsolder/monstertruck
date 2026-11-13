//! Intersection-loop store for the classic Boolean backend.
//!
//! This began as the published 0.3.2 loop store, but now also owns exact-boundary
//! endpoint resolution, coplanar interface matching, and partition stitching.
//! Geometric equality checks here are deliberately separated from broad-phase
//! signatures: sampled or summary evidence may select candidates, never prove
//! trimmed-region identity by itself.

#![expect(
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
mod ruled;
pub(super) use ruled::{
    certified_circular_cylinder, common_circular_extrusion_cylinder, common_ruled_extrusion_strip,
};

use monstertruck_core::{
    cgmath64::{EuclideanSpace, InnerSpace, MetricSpace, Point2, Point3, Vector2, Vector3},
    tolerance::DEFAULT_NUMERICAL_TOLERANCES,
};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::{Edge, EdgeId, Face, Shell, Vertex, VertexId, Wire, wire};
use planar::{
    FullPlanarInterfaceMatch, imprint_edges_on_faces, planar_polygons_same_trimmed_face,
    shells_full_planar_interface_pair, shells_have_partial_planar_overlap,
    shells_have_positive_aabb_overlap,
};
pub(in crate::transversal::classic) use planar::{
    direct_planar_face_match, single_face_matches_planar_partition,
};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

type PolylineCurve = monstertruck_meshing::prelude::PolylineCurve<Point3>;
type PathGraph<C> = (Vec<Edge<Point3, C>>, HashMap<VertexId<Point3>, Vec<usize>>);
type BranchBoundaryVertices = (
    Vertex<Point3>,
    Vertex<Point3>,
    Vertex<Point3>,
    Vertex<Point3>,
);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ShapesOpStatus {
    Unknown,
    And,
    Or,
}

impl ShapesOpStatus {
    const fn not(self) -> Self {
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
    pub(super) const fn new(wire: Wire<P, C>, status: ShapesOpStatus) -> Self {
        Self { wire, status }
    }
    pub(super) const fn status(&self) -> ShapesOpStatus { self.status }
    pub(super) fn inverse(&self) -> Self {
        Self {
            wire: self.wire.inverse(),
            status: self.status.not(),
        }
    }
}

impl ShapesOpStatus {
    fn from_is_curve<C, S0, S1>(curve: &IntersectionCurve<C, S0, S1>) -> Option<Self>
    where
        C: ParametricCurve3D + BoundedCurve,
        S0: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
        S1: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>, {
        let (t0, t1) = curve.range_tuple();
        let t = f64::midpoint(t0, t1);
        let (_, pt0, pt1) = curve.search_triple(t, 100)?;
        let der = curve.leader().derivative(t);
        let normal0 = curve.surface0().normal(pt0[0], pt0[1]);
        let normal1 = curve.surface1().normal(pt1[0], pt1[1]);
        if normal0.cross(der).dot(normal1) > 0.0 {
            Some(Self::Or)
        } else {
            Some(Self::And)
        }
    }
}

/// Give each regular SSI branch one geometrically defined tangent orientation.
///
/// Marching mesh segments are undirected. Region labels derived from the
/// ordered carrier-surface normals must not depend on the edge order of a
/// tessellation adjacency graph. Normalize the *same* raw and exact projected
/// curves together, and reject degenerate intersections with no consistent
/// oriented tangent rather than inventing a face-region label.
fn orient_intersection_branch<S>(
    polyline: &mut PolylineCurve,
    intersection_curve: &mut IntersectionCurve<PolylineCurve, S, S>,
) -> Option<ShapesOpStatus>
where
    S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3> + Clone,
{
    // Read the side once: reversing the leader flips its oriented sign.
    // Re-solving the nonlinear carrier intersection after inversion can fail
    // even when the original curve was valid (e.g. a nearly axial edge).
    if ShapesOpStatus::from_is_curve(intersection_curve)? == ShapesOpStatus::And {
        polyline.invert();
        intersection_curve.invert();
    }
    Some(ShapesOpStatus::Or)
}

impl<P, C> std::ops::Deref for BoundaryWire<P, C> {
    type Target = Wire<P, C>;
    fn deref(&self) -> &Self::Target { &self.wire }
}

impl<P, C> std::ops::DerefMut for BoundaryWire<P, C> {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.wire }
}

#[derive(Clone, Debug)]
pub(super) struct Loops<P, C>(Vec<BoundaryWire<P, C>>);
#[derive(Clone, Debug)]
pub(super) struct LoopsStore<P, C>(Vec<Loops<P, C>>);

impl<P, C> std::ops::Deref for Loops<P, C> {
    type Target = Vec<BoundaryWire<P, C>>;
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P, C> std::ops::DerefMut for Loops<P, C> {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P, C> std::ops::Deref for LoopsStore<P, C> {
    type Target = Vec<Loops<P, C>>;
    fn deref(&self) -> &Self::Target { &self.0 }
}

impl<P, C> std::ops::DerefMut for LoopsStore<P, C> {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}

impl<P, C> FromIterator<BoundaryWire<P, C>> for Loops<P, C> {
    fn from_iter<I: IntoIterator<Item = BoundaryWire<P, C>>>(iter: I) -> Self {
        Self(Vec::from_iter(iter))
    }
}

impl<'a, P, C, S> From<&'a Face<P, C, S>> for Loops<P, C> {
    fn from(face: &'a Face<P, C, S>) -> Self {
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
    fn try_new(t: f64, (t0, t1): (f64, f64)) -> Option<Self> {
        if t0.near(&t) {
            Some(Self::Front)
        } else if t1.near(&t) {
            Some(Self::Back)
        } else if t0 < t && t < t1 {
            Some(Self::Inner(t))
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

    fn add_independent_loop(&mut self, r#loop: BoundaryWire<P, C>) {
        self.push(r#loop.inverse());
        self.push(r#loop);
    }

    fn reuse_coincident_edge(&mut self, edge0: &Edge<P, C>, tolerances: BooleanTolerance) -> bool
    where
        P: Tolerance + MetricSpace<Metric = f64>,
        C: ParametricCurve<Point = P>
            + BoundedCurve<Point = P>
            + SearchNearestParameter<CurveParameter, Point = P>
            + Invertible, {
        let curves_sampled_match = |edge: &Edge<P, C>| {
            validate_sampled_coincidence::<P, _>(
                &edge.oriented_curve(),
                &edge0.oriented_curve(),
                tolerances,
            ) == CurveCoincidence::SampledValidation
        };
        let coincident = self.iter().flat_map(|wire| wire.iter()).find_map(|edge| {
            let same_direction = edge.front().point().distance(edge0.front().point())
                <= tolerances.length
                && edge.back().point().distance(edge0.back().point()) <= tolerances.length;
            let opposite_direction = edge.front().point().distance(edge0.back().point())
                <= tolerances.length
                && edge.back().point().distance(edge0.front().point()) <= tolerances.length;
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
                let same_direction = edge.front().point().distance(edge0.front().point())
                    <= tolerances.length
                    && edge.back().point().distance(edge0.back().point()) <= tolerances.length;
                let opposite_direction = edge.front().point().distance(edge0.back().point())
                    <= tolerances.length
                    && edge.back().point().distance(edge0.front().point()) <= tolerances.length;
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
        tolerances: BooleanTolerance,
    ) -> [Option<(usize, usize)>; 2]
    where
        P: Tolerance + MetricSpace<Metric = f64>,
        C: ParametricCurve<Point = P>
            + BoundedCurve<Point = P>
            + SearchNearestParameter<CurveParameter, Point = P>
            + Invertible,
    {
        if reuse_coincident && self.reuse_coincident_edge(&edge0, tolerances) {
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

fn partition_path_graph<C: Clone>(loops: &Loops<Point3, C>) -> PathGraph<C> {
    let mut edge_index_by_id = HashMap::default();
    let mut unique_edges = Vec::new();
    let mut adjacency = HashMap::<VertexId<Point3>, Vec<usize>>::default();
    for edge in loops.iter().flat_map(|wire| wire.iter()) {
        let edge_id = edge.id();
        if edge_index_by_id.contains_key(&edge_id) {
            continue;
        }
        let edge = edge.absolute_clone();
        let index = unique_edges.len();
        let front = edge.absolute_front().id();
        let back = edge.absolute_back().id();
        edge_index_by_id.insert(edge_id, index);
        unique_edges.push(edge);
        adjacency.entry(front).or_default().push(index);
        adjacency.entry(back).or_default().push(index);
    }
    (unique_edges, adjacency)
}

fn same_boundary_candidate<C>(
    existing: &Edge<Point3, C>,
    existing_endpoint: Option<VertexId<Point3>>,
    edge: &Edge<Point3, C>,
    endpoint: Option<VertexId<Point3>>,
    tolerances: BooleanTolerance,
) -> bool
where
    C: Clone
        + BoundedCurve<Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3>,
{
    if existing.id() == edge.id() || existing_endpoint.is_some() && existing_endpoint == endpoint {
        return true;
    }
    let same_endpoints = (existing
        .absolute_front()
        .point()
        .distance(edge.absolute_front().point())
        <= tolerances.length
        && existing
            .absolute_back()
            .point()
            .distance(edge.absolute_back().point())
            <= tolerances.length)
        || (existing
            .absolute_front()
            .point()
            .distance(edge.absolute_back().point())
            <= tolerances.length
            && existing
                .absolute_back()
                .point()
                .distance(edge.absolute_front().point())
                <= tolerances.length);
    same_endpoints
        && (validate_sampled_coincidence::<Point3, _>(&existing.curve(), &edge.curve(), tolerances)
            == CurveCoincidence::SampledValidation
            || validate_endpoint_aligned_coincidence::<Point3, _>(
                &existing.curve(),
                &edge.curve(),
                tolerances,
            ) == CurveCoincidence::SampledValidation)
}

impl<C> Loops<Point3, C> {
    /// Canonicalize an edge already present in this face's trim arrangement.
    ///
    /// Independent STEP carriers may represent the same physical edge with
    /// slightly different parameterizations and sub-kernel geometric noise.
    /// The imprint stage is explicitly allowed the caller's Boolean length
    /// tolerance, unlike the kernel-global exact incidence methods. Require
    /// matching oriented endpoints AND agreement of the whole sampled path;
    /// otherwise a genuine thin lens sharing its endpoints must not collapse.
    fn reuse_imprinted_coincident_edge(
        &mut self,
        incoming: &Edge<Point3, C>,
        length_tol: f64,
    ) -> bool
    where
        C: ParametricCurve<Point = Point3>
            + BoundedCurve<Point = Point3>
            + ParameterDivision1D<Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>
            + Invertible
            + Clone,
    {
        // Reject an invalid tolerance; silently clamping it would change
        // which distinct trimmed edges can be welded into one identity.
        if !length_tol.is_finite() || length_tol <= 0.0 {
            return false;
        }
        let same_endpoints = |edge: &Edge<Point3, C>| {
            if edge.front().point().distance(incoming.front().point()) <= length_tol
                && edge.back().point().distance(incoming.back().point()) <= length_tol
            {
                Some(true)
            } else if edge.front().point().distance(incoming.back().point()) <= length_tol
                && edge.back().point().distance(incoming.front().point()) <= length_tol
            {
                Some(false)
            } else {
                None
            }
        };
        let same_path = |edge: &Edge<Point3, C>, _same_direction: bool| {
            // Trimmed curves can have unrelated parameterizations (especially
            // after STEP import and independent SSI reconstruction). Equality
            // at a normalized parameter fraction is NOT geometric equality.
            // Require both directed paths to project back onto the other
            // exact curve within the caller's actual length tolerance.
            let a = edge.curve();
            let b = incoming.curve();
            let samples_on = |curve: &C, other: &C| {
                let (t0, t1) = curve.range_tuple();
                // The adaptive division has *already evaluated* every
                // returned point on the exact curve. Do not throw those
                // results away and redo millions of expensive 4D SSI solves
                // merely to re-check their exact same parameter.
                let (parameters, points) = curve.parameter_division((t0, t1), length_tol);
                if parameters.len() != points.len() {
                    return false;
                }
                let mut previous = None;
                let mut verify = |point: Point3| {
                    let nearest = |hint| {
                        let parameter = other.search_nearest_parameter(point, hint, 100)?;
                        (other.evaluate(parameter).distance(point) <= length_tol)
                            .then_some(parameter)
                    };
                    let parameter = nearest(previous).or_else(|| nearest(None));
                    if let Some(parameter) = parameter {
                        previous = Some(parameter);
                        true
                    } else {
                        false
                    }
                };
                let dense = parameters
                    .into_iter()
                    .zip(points)
                    .map(|(parameter, point)| {
                        let fraction = (parameter - t0) / (t1 - t0);
                        // Original validation normalized and reconstructed every
                        // parameter. Reuse the stored exact evaluation when those
                        // floating-point operations return exactly the same bits.
                        // Otherwise retain the original re-evaluation to preserve
                        // all existing strict acceptance behavior, even at seams.
                        let mapped = (t1 - t0).mul_add(fraction, t0);
                        (
                            fraction,
                            if mapped.to_bits() == parameter.to_bits() {
                                point
                            } else {
                                curve.evaluate(mapped)
                            },
                        )
                    });
                let coarse = (0..=16).map(|i| {
                    let fraction = f64::from(i) / 16.0;
                    (fraction, curve.evaluate((t1 - t0).mul_add(fraction, t0)))
                });
                dense.chain(coarse).all(|(fraction, point)| {
                    fraction.is_finite() && (0.0..=1.0).contains(&fraction) && verify(point)
                })
            };
            samples_on(&a, &b) && samples_on(&b, &a)
        };
        let matching = self.iter().flat_map(|wire| wire.iter()).find_map(|edge| {
            let direction = same_endpoints(edge)?;
            same_path(edge, direction)
                .then(|| (edge.front().clone(), edge.back().clone(), direction))
        });
        let Some((old_front, old_back, direction)) = matching else {
            return false;
        };
        let (new_front, new_back) = if direction {
            (incoming.front().clone(), incoming.back().clone())
        } else {
            (incoming.back().clone(), incoming.front().clone())
        };
        let mut map = HashMap::default();
        if old_front.id() != new_front.id() {
            self.change_vertex(&old_front, &new_front, &mut map);
        }
        if old_back.id() != new_back.id() {
            self.change_vertex(&old_back, &new_back, &mut map);
        }

        // Only canonicalize edges individually proven to coincide. Another
        // contour can share both vertices without describing the same trim.
        self.iter_mut().for_each(|wire| {
            wire.iter_mut().for_each(|edge| {
                if let Some(direction) = same_endpoints(edge)
                    && same_path(edge, direction)
                {
                    *edge = if direction {
                        incoming.clone()
                    } else {
                        incoming.inverse()
                    };
                }
            });
        });
        true
    }

    fn has_partitioned_coincident_path(&self, edge0: &Edge<Point3, C>) -> bool
    where C: Clone
            + ParametricCurve<Point = Point3>
            + BoundedCurve<Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3> {
        let reference = edge0.curve();
        let reference_range = reference.range_tuple();
        let start = edge0.front().id();
        let goal = edge0.back().id();

        let (unique_edges, adjacency) = partition_path_graph(self);

        if !adjacency.contains_key(&start) || !adjacency.contains_key(&goal) {
            return false;
        }

        // Geometric coincidence is expensive for intersection curves. Explore the
        // topological graph first and prove only edges reachable from the incoming
        // edge's start. This is equivalent to filtering every edge up front, while
        // avoiding inverse-mapping unrelated wires and branches.
        let mut coincident = vec![None; unique_edges.len()];
        let mut frontier = vec![start];
        let mut seen = HashSet::default();
        seen.insert(start);
        while let Some(current) = frontier.pop() {
            let Some(incident) = adjacency.get(&current) else {
                continue;
            };
            for &index in incident {
                let edge = &unique_edges[index];
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
                if next != goal && seen.contains(&next) {
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
                if seen.insert(next) {
                    frontier.push(next);
                }
            }
        }
        false
    }
}

impl<P: Copy + Tolerance, C: Clone> LoopsStore<P, C> {
    fn change_vertex(
        &mut self,
        old_vertex: &Vertex<P>,
        new_vertex: &Vertex<P>,
        emap: &mut HashMap<EdgeId<C>, Edge<P, C>>,
    ) {
        self.iter_mut()
            .for_each(|loops| loops.change_vertex(old_vertex, new_vertex, emap));
    }

    fn swap_edge_into_wire(&mut self, edge_id: EdgeId<C>, new_wire: &Wire<P, C>) {
        self.iter_mut()
            .for_each(|loops| loops.swap_edge_into_wire(edge_id, new_wire));
    }
}

/// Prove that a candidate canonical boundary point lies on the other exact
/// Boolean carrier, not merely within the coarse tessellation search radius.
///
/// Inverse search may find an approximate UV for an off-surface point. Only its
/// forward-evaluated spatial residual is evidence of incidence.
fn on_other_carrier<S>(point: Point3, other_surface: &S, tolerances: BooleanTolerance) -> bool
where S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3> {
    other_surface
        .search_nearest_parameter(point, None, 25)
        .is_some_and(|(u, v)| {
            other_surface.evaluate(u, v).distance(point) <= tolerances.surface_support_residual()
        })
}

impl<C> LoopsStore<Point3, C> {
    fn search_boundary_parameter(
        &self,
        loops_index: usize,
        point: Point3,
        nearest_tolerance: f64,
        tolerances: BooleanTolerance,
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
                return same_boundary_candidate(
                    existing,
                    candidate_endpoint,
                    edge,
                    endpoint,
                    tolerances,
                );
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
                // SearchParameter can use its own kernel epsilon. A returned
                // parameter is only exact at the operation's requested spatial
                // precision if its forward evaluation agrees independently.
                if curve.evaluate(parameter).distance(point) > tolerances.length {
                    continue;
                }
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
        tolerances: BooleanTolerance,
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
                // As with exact boundary inverse search, the curve's own
                // SearchParameter tolerance is not authoritative. Two trims
                // separated by less than the global epsilon must not be
                // imprinted onto each other at a tighter operation tolerance.
                if curve.evaluate(parameter).distance(point) > tolerances.length {
                    continue;
                }
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
                                .distance(edge.absolute_front().point())
                                <= tolerances.length
                                && reference
                                    .absolute_back()
                                    .point()
                                    .distance(edge.absolute_back().point())
                                    <= tolerances.length)
                                || (reference
                                    .absolute_front()
                                    .point()
                                    .distance(edge.absolute_back().point())
                                    <= tolerances.length
                                    && reference
                                        .absolute_back()
                                        .point()
                                        .distance(edge.absolute_front().point())
                                        <= tolerances.length);
                            let reference_curve = reference.curve();
                            let edge_curve = edge.curve();
                            let same_curve = validate_sampled_coincidence::<Point3, _>(
                                &reference_curve,
                                &edge_curve,
                                tolerances,
                            ) == CurveCoincidence::SampledValidation
                                || validate_sampled_containment::<Point3, _>(
                                    &reference_curve,
                                    &edge_curve,
                                    tolerances,
                                ) == CurveCoincidence::SampledValidation
                                || validate_sampled_containment::<Point3, _>(
                                    &edge_curve,
                                    &reference_curve,
                                    tolerances,
                                ) == CurveCoincidence::SampledValidation
                                || (same_endpoints
                                    && validate_endpoint_aligned_coincidence::<Point3, _>(
                                        &reference_curve,
                                        &edge_curve,
                                        tolerances,
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
            self.add_polygon_vertex(loops_index, vertex, tolerances, emap)?;
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
        tolerances: BooleanTolerance,
        emap: &mut HashMap<EdgeId<C>, Edge<Point3, C>>,
    ) -> Option<(usize, usize, ParameterKind)>
    where
        C: Clone
            + Cut<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let (wire_index, edge_index, kind, projected) = self.search_boundary_parameter(
            loops_index,
            vertex.point(),
            tolerances.boundary_residual(),
            tolerances,
        )?;
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
        tolerances: BooleanTolerance,
    ) -> Option<(usize, usize, ParameterKind, Option<Point3>)>
    where
        C: BoundedCurve<Point = Point3>
            + SearchParameter<CurveParameter, Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>,
    {
        let edge = self.get(loops_index)?.get(wire_index)?.get(edge_index)?;
        let curve = edge.curve();
        let classify = |parameter: f64, residual_limit: f64| {
            let projected = curve.evaluate(parameter);
            if projected.distance(point) > residual_limit {
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
            .and_then(|t| classify(t, tolerances.length))
            .or_else(|| {
                curve
                    .search_nearest_parameter(point, None, 10)
                    .and_then(|t| classify(t, nearest_tolerance))
            })
    }

    /// Resolve a sampled SSI endpoint to an *existing* exact face vertex.
    ///
    /// After prior SSI splits, the polygonal wire index need no longer identify
    /// the same exact wire/edge. Search the entire face's existing vertices,
    /// selecting a unique nearest geometric locus within the tessellation
    /// seed error, and independently verify that locus on the intersecting
    /// surface before changing any topology. Distinct nearly equidistant
    /// locations are ambiguous, even if their edge inverse solvers disagree.
    fn supported_boundary_vertex<S>(
        &self,
        face_index: usize,
        point: Point3,
        nearest_tolerance: f64,
        tolerances: BooleanTolerance,
        other_surface: &S,
    ) -> Option<Vec<Vertex<Point3>>>
    where
        S: ParametricSurface3D + SearchNearestParameter<SurfaceParameter, Point = Point3>,
    {
        let ambiguity_tol = tolerances.boundary_residual();
        let mut candidates: Vec<(Point3, f64, Vec<Vertex<Point3>>)> = Vec::new();
        for vertex in self
            .get(face_index)?
            .iter()
            .flat_map(|wire| wire.vertex_iter())
        {
            let location = vertex.point();
            let distance = location.distance(point);
            if distance > nearest_tolerance {
                continue;
            }
            if !on_other_carrier(location, other_surface, tolerances) {
                continue;
            }
            if let Some((_, _, vertices)) = candidates
                .iter_mut()
                .find(|(existing, _, _)| existing.distance(location) <= tolerances.length)
            {
                if !vertices.iter().any(|known| known.id() == vertex.id()) {
                    vertices.push(vertex.clone());
                }
                continue;
            }
            candidates.push((location, distance, vec![vertex.clone()]));
        }
        candidates.sort_by(|lhs, rhs| lhs.1.total_cmp(&rhs.1));
        let (_, closest_distance, selected) = candidates.first()?;
        // Do not choose by index among two distinguishable but numerically
        // indistinguishable exact vertices.
        if candidates
            .get(1)
            .is_some_and(|(_, next_distance, _)| next_distance - closest_distance <= ambiguity_tol)
        {
            return None;
        }
        Some(selected.clone())
    }

    fn add_geom_vertex<S>(
        &mut self,
        loops_index: usize,
        boundary_hint: Option<(usize, usize)>,
        sampled_on_boundary_vertex: bool,
        vertex: &Vertex<Point3>,
        another_surface: &S,
        nearest_tolerance: f64,
        tolerances: BooleanTolerance,
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
            self.resolve_boundary_hint(
                loops_index,
                hint,
                vertex.point(),
                nearest_tolerance,
                tolerances,
            )
        });
        let resolved = hinted.or_else(|| {
            self.search_boundary_parameter(
                loops_index,
                vertex.point(),
                nearest_tolerance,
                tolerances,
            )
        });
        if resolved.is_none() && sampled_on_boundary_vertex {
            if let Some(vertices) = self.supported_boundary_vertex(
                loops_index,
                vertex.point(),
                nearest_tolerance,
                tolerances,
                another_surface,
            ) {
                let canonical = vertices[0].point();
                vertex.set_point(canonical);
                for old in vertices {
                    self.change_vertex(&old, vertex, emap);
                }
                return Some(());
            }
        }
        let (wire_index, edge_index, kind, projected) = resolved?;
        // Nearest-curve search is only discovery. Before coalescing an
        // *existing* exact endpoint, certify its location against the
        // intersecting exact surface, independently of the mesh error budget.
        // Perform this before set_point/change_vertex mutates either store.
        let wire_edge = &self[loops_index][wire_index][edge_index];
        let endpoint = match kind {
            ParameterKind::Front => Some(wire_edge.absolute_front().point()),
            ParameterKind::Back => Some(wire_edge.absolute_back().point()),
            ParameterKind::Inner(_) => None,
        };
        if endpoint.is_some_and(|point| !on_other_carrier(point, another_surface, tolerances)) {
            return None;
        }
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
                    tolerances,
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
    tolerances: BooleanTolerance,
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

    // SSI seeds originate from tessellated carriers. The *seed* must lie
    // within the operation's boundary residual of both projections, while the
    // independently projected curve/surface points must themselves agree
    // within the operation's numerical length (rather than the independently
    // coarser SSI support-projection envelope). Near tangencies can
    // make a Newton update singular even when the two exact projections
    // already agree; accept that verified case before attempting a step.
    // The Standard modeling-error budget never widens these thresholds.
    let seed_tolerance = tolerances.boundary_residual();
    if pt0.distance(pt1) <= tolerances.length
        && point.distance(pt0) <= seed_tolerance
        && point.distance(pt1) <= seed_tolerance
    {
        // Keep the point exactly on the boundary curve. Independent carrier
        // evaluations agreed to the operation's numerical length, and `t`
        // remains the exact split parameter of this trim edge.
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
        tolerances,
        trials - 1,
    )
}

fn unique_nearby_loop_point<C>(
    lhs: &Loops<Point3, C>,
    rhs: &Loops<Point3, C>,
    point: Point3,
    tolerances: BooleanTolerance,
) -> Option<Point3> {
    let mut candidate = None;
    let tolerance = tolerances.cross_face_stitch_residual();
    for nearby in lhs
        .iter()
        .chain(rhs.iter())
        .flat_map(|wire| wire.vertex_iter())
        .map(|vertex| vertex.point())
        .filter(|nearby| (*nearby - point).magnitude() <= tolerance)
    {
        match candidate {
            None => candidate = Some(nearby),
            Some(existing) if existing.distance(nearby) <= tolerances.length => {}
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
    let t = f64::midpoint(t0, t1);
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

struct BranchContext<'a, C, S> {
    operand0: &'a mut OperandLoopStores<C>,
    operand1: &'a mut OperandLoopStores<C>,
    face_indices: [usize; 2],
    surfaces: [&'a S; 2],
    nearest_tolerance: f64,
    tolerances: BooleanTolerance,
    reuse_coincident: bool,
}

#[derive(Clone, Copy)]
struct BoundaryEndpoint<'a, S> {
    face_index: usize,
    poly_vertex: &'a Vertex<Point3>,
    geom_vertex: &'a Vertex<Point3>,
    other_surface: &'a S,
    nearest_tolerance: f64,
    tolerances: BooleanTolerance,
}

impl<C> OperandLoopStores<C> {
    fn add_boundary_endpoint<S>(
        &mut self,
        endpoint: &BoundaryEndpoint<'_, S>,
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
            tolerances,
        } = *endpoint;
        let Some((wire_index, edge_index, polygon_kind)) =
            self.poly
                .add_polygon_vertex(face_index, poly_vertex, tolerances, poly_edge_map)
        else {
            return Some(false);
        };
        self.geom.add_geom_vertex(
            face_index,
            Some((wire_index, edge_index)),
            matches!(polygon_kind, ParameterKind::Front | ParameterKind::Back),
            geom_vertex,
            other_surface,
            nearest_tolerance,
            tolerances,
            geom_edge_map,
        )?;
        // Both stores must share the same canonical endpoint coordinate.
        // A split of an exact curve is not represented by the polygon vertex
        // retaining its old mesh-only position.
        poly_vertex.set_point(geom_vertex.point());
        Some(true)
    }

    fn add_intersection_edge(
        &mut self,
        face_index: usize,
        poly_edge: Edge<Point3, PolylineCurve>,
        geom_edge: Edge<Point3, C>,
        status: ShapesOpStatus,
        reuse_coincident: bool,
        tolerances: BooleanTolerance,
    ) where
        C: ParametricCurve<Point = Point3>
            + BoundedCurve<Point = Point3>
            + SearchNearestParameter<CurveParameter, Point = Point3>
            + Invertible
            + Clone,
    {
        self.poly[face_index].add_edge(poly_edge, status, reuse_coincident, tolerances);
        self.geom[face_index].add_edge(geom_edge, status, reuse_coincident, tolerances);
    }
}

fn add_branch_boundary_vertices<C, S>(
    context: &mut BranchContext<'_, C, S>,
    front: Point3,
    back: Point3,
) -> Option<BranchBoundaryVertices>
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
    let pv0 = Vertex::new(front);
    let pv1 = Vertex::new(back);
    let gv0 = Vertex::new(front);
    let gv1 = Vertex::new(back);
    let mut pemap0 = HashMap::default();
    let mut pemap1 = HashMap::default();
    let mut gemap0 = HashMap::default();
    let mut gemap1 = HashMap::default();
    for (operand, face_index, other_surface) in [
        (
            &mut *context.operand0,
            context.face_indices[0],
            context.surfaces[1],
        ),
        (
            &mut *context.operand1,
            context.face_indices[1],
            context.surfaces[0],
        ),
    ] {
        operand.add_boundary_endpoint(
            &BoundaryEndpoint {
                face_index,
                poly_vertex: &pv0,
                geom_vertex: &gv0,
                other_surface,
                nearest_tolerance: context.nearest_tolerance,
                tolerances: context.tolerances,
            },
            &mut pemap0,
            &mut gemap0,
        )?;
        operand.add_boundary_endpoint(
            &BoundaryEndpoint {
                face_index,
                poly_vertex: &pv1,
                geom_vertex: &gv1,
                other_surface,
                nearest_tolerance: context.nearest_tolerance,
                tolerances: context.tolerances,
            },
            &mut pemap1,
            &mut gemap1,
        )?;
    }
    Some((pv0, pv1, gv0, gv1))
}

fn add_open_intersection_branch<C, S>(
    context: &mut BranchContext<'_, C, S>,
    mut polyline: PolylineCurve,
    mut intersection_curve: IntersectionCurve<PolylineCurve, S, S>,
    statuses: [ShapesOpStatus; 2],
) -> Option<()>
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
    let mut front = polyline.front();
    let mut back = polyline.back();
    if let Some(point) = unique_nearby_loop_point(
        &context.operand0.geom[context.face_indices[0]],
        &context.operand1.geom[context.face_indices[1]],
        front,
        context.tolerances,
    ) {
        front = point;
    }
    if let Some(point) = unique_nearby_loop_point(
        &context.operand0.geom[context.face_indices[0]],
        &context.operand1.geom[context.face_indices[1]],
        back,
        context.tolerances,
    ) {
        back = point;
    }
    // Canonicalization can collapse a numerically tiny SSI fragment to a point.
    if front.distance(back) <= context.tolerances.length {
        return Some(());
    }
    *polyline.first_mut()? = front;
    *polyline.last_mut()? = back;
    *intersection_curve.leader_mut().first_mut()? = front;
    *intersection_curve.leader_mut().last_mut()? = back;

    let (pv0, pv1, gv0, gv1) = add_branch_boundary_vertices(context, front, back)?;
    *intersection_curve.leader_mut().first_mut()? = gv0.point();
    *intersection_curve.leader_mut().last_mut()? = gv1.point();

    let pedge = Edge::new_unchecked(&pv0, &pv1, polyline);
    let gedge = Edge::new_unchecked(&gv0, &gv1, intersection_curve.into());
    context.operand0.add_intersection_edge(
        context.face_indices[0],
        pedge.clone(),
        gedge.clone(),
        statuses[0],
        context.reuse_coincident,
        context.tolerances,
    );
    context.operand1.add_intersection_edge(
        context.face_indices[1],
        pedge,
        gedge,
        statuses[1],
        context.reuse_coincident,
        context.tolerances,
    );
    Some(())
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
}

struct CoplanarInterfaceAnalysis {
    partial_overlap: bool,
    any_overlap: bool,
    full_interface_adjacency: Option<FullPlanarInterfaceMatch>,
}

fn analyze_coplanar_interfaces(
    poly_shell0: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    poly_shell1: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tolerance: f64,
    imprint_coplanar: bool,
) -> Option<CoplanarInterfaceAnalysis> {
    let positive_aabb_overlap =
        shells_have_positive_aabb_overlap(poly_shell0, poly_shell1, tolerance)?;
    let full_interface = if imprint_coplanar {
        shells_full_planar_interface_pair(poly_shell0, poly_shell1, tolerance).ok()?
    } else {
        None
    };
    let full_interface_adjacency = if positive_aabb_overlap {
        None
    } else {
        full_interface
            .clone()
            .filter(|matched| !matched.lhs_faces.is_empty() && !matched.rhs_faces.is_empty())
    };
    // Partial planar interfaces can be exactly zero-volume contacts.
    // In that case a full-interface shortcut does not apply, but the
    // face-local positive-area test still must be considered.
    let partial_overlap = imprint_coplanar
        && full_interface_adjacency.is_none()
        && shells_have_partial_planar_overlap(poly_shell0, poly_shell1, tolerance)?;
    Some(CoplanarInterfaceAnalysis {
        partial_overlap,
        any_overlap: full_interface.is_some() || partial_overlap,
        full_interface_adjacency,
    })
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
        + ParameterDivision1D<Point = Point3>
        + Cut<Point = Point3, Vector = Vector3>
        + From<IntersectionCurve<PolylineCurve, S, S>>
        + Invertible,
    S: ParametricSurface3D
        + Clone
        + Invertible
        + TryIntoHomogeneousBsplineSurface
        + TryIntoAnalyticSurfaceKind
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
    let coplanar_tol =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_coplanar_scale * boolean_tolerance.length;
    let boundary_search_tolerance = boolean_tolerance.mesh_seed;
    let coplanar =
        analyze_coplanar_interfaces(poly_shell0, poly_shell1, coplanar_tol, imprint_coplanar)?;
    // Proof-first common-carrier discovery: identical rational profile control
    // nets with collinear affine sweep rows describe the same ruled support.
    // Cached once per input face; sampling never authorizes surface identity.
    let net0 = geom_shell0
        .iter()
        .map(|f| f.surface().try_into_homogeneous_bspline_surface())
        .collect::<Vec<_>>();
    let net1 = geom_shell1
        .iter()
        .map(|f| f.surface().try_into_homogeneous_bspline_surface())
        .collect::<Vec<_>>();
    // Full carrier recognition/control-net certification is done once per
    // face, not once for every O(F0*F1) candidate pair.
    let planes0 = geom_shell0
        .iter()
        .zip(net0.iter())
        .map(|(face, net)| intersection_curve::certified_plane(&face.surface(), net.as_ref()))
        .collect::<Vec<_>>();
    let planes1 = geom_shell1
        .iter()
        .zip(net1.iter())
        .map(|(face, net)| intersection_curve::certified_plane(&face.surface(), net.as_ref()))
        .collect::<Vec<_>>();
    let mut shared_ruled_pairs = HashSet::default();
    for (i, left) in net0.iter().enumerate() {
        for (j, right) in net1.iter().enumerate() {
            if let (Some(a), Some(b)) = (left, right)
                && (ruled::common_ruled_extrusion_strip(a, b)
                    || ruled::common_circular_extrusion_cylinder(a, b))
            {
                shared_ruled_pairs.insert((i, j));
            }
        }
    }
    let has_shared_ruled = !shared_ruled_pairs.is_empty();
    (0..store0_len)
        .flat_map(move |i| (0..store1_len).map(move |j| (i, j)))
        .try_for_each(|(face_index0, face_index1)| {
            // Solids that meet on one complete trimmed planar face with
            // opposite outward normals need no SSI at all. Generic SSI on
            // their coincident continuation surfaces (e.g. coaxial cylinders)
            // is degenerate; downstream coplanar classification removes the
            // shared internal face and topology welding stitches the rim.
            if coplanar.full_interface_adjacency.is_some() {
                return Some(());
            }
            let ori0 = geom_shell0[face_index0].orientation();
            let ori1 = geom_shell1[face_index1].orientation();
            let surface0 = geom_shell0[face_index0].surface();
            let surface1 = geom_shell1[face_index1].surface();
            let polygon0 = poly_shell0[face_index0].surface()?;
            let polygon1 = poly_shell1[face_index1].surface()?;
            // A common ruled carrier has a two-dimensional overlap: a
            // transverse SSI curve is undefined. Partition its real trimmed
            // boundaries via exact edge imprinting below instead.
            if shared_ruled_pairs.contains(&(face_index0, face_index1)) {
                return Some(());
            }
            // Exactly coincident trimmed planar faces are already a complete
            // interface. Running generic SSI on coincident surfaces is both
            // unnecessary and degenerate; ownership is resolved by the
            // coplanar face classifier after division.
            if planar_polygons_same_trimmed_face(&polygon0, &polygon1, coplanar_tol) {
                // The two complete trims already bound the same planar region.
                // A transverse surface/surface curve is not defined on that
                // two-dimensional coincidence, for difference as for union.
                return Some(());
            }
            let curves = intersection_curve::intersection_curves(
                surface0.clone(),
                &polygon0,
                surface1.clone(),
                &polygon1,
                boolean_tolerance,
                (planes0[face_index0], planes1[face_index1]),
            )?;
            curves
                .into_iter()
                .try_for_each(|(mut polyline, mut intersection_curve)| {
                    let status =
                        orient_intersection_branch(&mut polyline, &mut intersection_curve)?;
                    let (status0, status1) = match (ori0, ori1) {
                        (true, true) => (status, status.not()),
                        (true, false) => (status.not(), status.not()),
                        (false, true) => (status, status),
                        (false, false) => (status.not(), status),
                    };
                    if polyline.front().distance(polyline.back()) <= boolean_tolerance.length {
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
                        let mut context = BranchContext {
                            operand0: &mut operand0,
                            operand1: &mut operand1,
                            face_indices: [face_index0, face_index1],
                            surfaces: [&surface0, &surface1],
                            nearest_tolerance: boundary_search_tolerance,
                            tolerances: boolean_tolerance,
                            reuse_coincident: coplanar.partial_overlap,
                        };
                        add_open_intersection_branch(
                            &mut context,
                            polyline,
                            intersection_curve,
                            [status0, status1],
                        )?;
                    }
                    Some(())
                })
        })?;

    // Mesh/mesh interference intentionally ignores coplanar contact. By this
    // point ordinary SSI has already split the boundary edges at transverse
    // intersections, so imprint those existing edge segments onto any face of
    // the opposite operand that contains them. This supplies the missing trim
    // lines for coplanar overlap without manufacturing approximate SSI curves.
    if coplanar.partial_overlap || has_shared_ruled {
        let source0 = operand0.geom.clone();
        let source1 = operand1.geom.clone();
        imprint_edges_on_faces(&source1, poly_shell0, &mut operand0.geom, tol)?;
        imprint_edges_on_faces(&source0, poly_shell1, &mut operand1.geom, tol)?;
    }

    let sampled0 = (!(coplanar.partial_overlap || has_shared_ruled)
        && loop_stores_are_synchronized(&operand0.geom, &operand0.poly))
    .then_some(operand0.poly);
    let sampled1 = (!(coplanar.partial_overlap || has_shared_ruled)
        && loop_stores_are_synchronized(&operand1.geom, &operand1.poly))
    .then_some(operand1.poly);
    Some(InterfaceLoops {
        geom: [operand0.geom, operand1.geom],
        sampled: [sampled0, sampled1],
        coplanar_overlap: coplanar.any_overlap,
    })
}

#[cfg(test)]
mod projection_tests {
    use super::*;
    use monstertruck_core::tolerance::KERNEL_EPSILON;

    #[test]
    fn ssi_carrier_normal_orientation_is_independent_of_polyline_direction() {
        let start = Point3::new(0.0, -1.0, 0.0);
        let end = Point3::new(0.0, 1.0, 0.0);
        let mut normalized = Vec::new();
        for points in [vec![start, end], vec![end, start]] {
            let mut polyline = PolylineCurve(points);
            let mut curve = IntersectionCurve::new(Plane::xy(), Plane::yz(), polyline.clone());
            assert_eq!(
                orient_intersection_branch(&mut polyline, &mut curve),
                Some(ShapesOpStatus::Or)
            );
            assert_eq!(curve.leader().0, polyline.0);
            normalized.push(polyline.0);
        }
        assert_eq!(normalized[0], normalized[1]);
    }

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
        let point = Point3::new(0.5, 1.125 * KERNEL_EPSILON, 0.0);

        let (_, _, kind, projected) = store
            .search_boundary_parameter(
                0,
                point,
                8.0 * KERNEL_EPSILON,
                BooleanTolerance::new(KERNEL_EPSILON),
            )
            .expect("unique bounded nearest projection should recover the boundary parameter");

        assert!(matches!(kind, ParameterKind::Inner(_)));
        assert_eq!(projected, Some(Point3::new(0.5, 0.0, 0.0)));
    }

    #[test]
    fn boundary_exact_inverse_requires_requested_spatial_precision() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let z = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&a, &z, Line(a.point(), z.point()));
        let store = LoopsStore(vec![Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )])]);
        let query = Point3::new(0.5, 5.0e-7, 0.0);
        let precision = BooleanTolerance::new(1.0e-8);
        // The line inverse itself uses the generic kernel epsilon and
        // regards this point as on the line. The Boolean MUST reject that
        // classification at the tighter operation-specific tolerance.
        assert!(
            store
                .search_boundary_parameter(0, query, precision.boundary_residual(), precision)
                .is_none()
        );
        // A separate mesh candidate search is allowed to find the line:
        // candidate proximity is not an exact-incidence certificate.
        let (_, _, kind, projected) = store
            .search_boundary_parameter(0, query, 1.0e-3, precision)
            .expect("coarse candidate search should return the exact projection");
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
        let vertex = Vertex::new(Point3::new(0.5, 1.125 * KERNEL_EPSILON, 0.0));
        let mut emap = HashMap::default();

        let (_, _, kind) = store
            .add_polygon_vertex(0, &vertex, BooleanTolerance::new(KERNEL_EPSILON), &mut emap)
            .expect("recovered boundary projection should split the exact edge");

        assert!(matches!(kind, ParameterKind::Inner(_)));
        assert_eq!(vertex.point(), Point3::new(0.5, 0.0, 0.0));
        assert_eq!(store[0][0].len(), 2);
    }

    #[test]
    fn coarse_endpoint_candidates_must_belong_to_other_exact_carrier() {
        // The input SSI seed lies on z=0 but the unique nearby boundary edge
        // is actually on z=5e-7. Global Point::near and mesh search would
        // have incorrectly welded that boundary into the plane.
        let gap = 5.0e-7;
        for hint in [None, Some((0, 0))] {
            let old_front = Vertex::new(Point3::new(0.0, 0.0, gap));
            let old_back = Vertex::new(Point3::new(1.0, 0.0, gap));
            let edge = Edge::new_unchecked(
                &old_front,
                &old_back,
                Line(old_front.point(), old_back.point()),
            );
            let mut store = LoopsStore(vec![Loops(vec![BoundaryWire::new(
                vec![edge].into(),
                ShapesOpStatus::Unknown,
            )])]);
            let seed = Vertex::new(Point3::new(0.0, 0.0, 0.0));
            let mut edge_map = HashMap::default();
            assert!(
                store
                    .add_geom_vertex(
                        0,
                        hint,
                        false,
                        &seed,
                        &Plane::xy(),
                        1.0e-3,
                        BooleanTolerance::new(1.0e-8),
                        &mut edge_map,
                    )
                    .is_none(),
                "neither a hinted nor a global coarse candidate proves opposite-carrier incidence"
            );
            assert_eq!(store[0][0].len(), 1);
            assert_eq!(store[0][0][0].absolute_front().id(), old_front.id());
            assert_eq!(old_front.point(), Point3::new(0.0, 0.0, gap));
            assert_eq!(seed.point(), Point3::new(0.0, 0.0, 0.0));
            assert!(edge_map.is_empty());
        }
    }

    #[test]
    fn coarse_endpoint_search_may_snap_to_proven_exact_carrier() {
        // A mesh seed may be off the carrier. Its exact boundary candidate,
        // however, lies on the other surface and is safe to canonicalize.
        let old_front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let old_back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(
            &old_front,
            &old_back,
            Line(old_front.point(), old_back.point()),
        );
        let mut store = LoopsStore(vec![Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )])]);
        let seed = Vertex::new(Point3::new(0.0, 0.0, 5.0e-7));
        let mut edge_map = HashMap::default();
        assert!(
            store
                .add_geom_vertex(
                    0,
                    Some((0, 0)),
                    false,
                    &seed,
                    &Plane::xy(),
                    1.0e-3,
                    BooleanTolerance::new(1.0e-8),
                    &mut edge_map,
                )
                .is_some()
        );
        assert_eq!(seed.point(), old_front.point());
        assert_eq!(store[0][0][0].absolute_front().id(), seed.id());
        assert_eq!(store[0][0].len(), 1);
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
            false,
            &recovered,
            &Plane::xy(),
            8.0 * KERNEL_EPSILON,
            BooleanTolerance::new(KERNEL_EPSILON),
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
                .search_boundary_parameter(
                    0,
                    point,
                    8.0 * KERNEL_EPSILON,
                    BooleanTolerance::new(KERNEL_EPSILON)
                )
                .is_none()
        );
        let recovered = store.search_boundary_parameter(
            0,
            point,
            1.0e-3,
            BooleanTolerance::new(KERNEL_EPSILON),
        );
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

        let recovered = store.search_boundary_parameter(
            0,
            point,
            1.0e-3,
            BooleanTolerance::new(KERNEL_EPSILON),
        );
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
                .search_boundary_parameter(
                    0,
                    Point3::new(0.0, 0.0, 0.0),
                    1.0e-3,
                    BooleanTolerance::new(KERNEL_EPSILON)
                )
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

        assert!(
            store
                .add_imprint_vertex(0, &split, BooleanTolerance::new(KERNEL_EPSILON), &mut emap)
                .is_some()
        );
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

        assert!(
            store
                .add_imprint_vertex(0, &split, BooleanTolerance::new(KERNEL_EPSILON), &mut emap)
                .is_some()
        );
        assert_eq!(store[0][0].len(), 2);
        assert_eq!(store[0][1].len(), 2);
        assert!(
            store[0]
                .iter()
                .all(|wire| wire.vertex_iter().any(|v| v == split))
        );
    }

    #[test]
    fn imprint_does_not_split_nearby_distinct_edges_with_sub_epsilon_tolerance() {
        let gap = 5.0e-7;
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let c = Vertex::new(Point3::new(0.0, 0.0, gap));
        let d = Vertex::new(Point3::new(1.0, 0.0, gap));
        let edge0 = Edge::new_unchecked(&a, &b, Line(a.point(), b.point()));
        let edge1 = Edge::new_unchecked(&c, &d, Line(c.point(), d.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![edge0].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge1].into(), ShapesOpStatus::Unknown),
        ]);
        let mut store = LoopsStore(vec![loops]);
        let split = Vertex::new(Point3::new(0.5, 0.0, 0.0));
        let mut edge_map = HashMap::default();

        assert!(
            store
                .add_imprint_vertex(0, &split, BooleanTolerance::new(1.0e-8), &mut edge_map)
                .is_some()
        );
        assert_eq!(store[0][0].len(), 2);
        assert_eq!(
            store[0][1].len(),
            1,
            "the parallel edge 5e-7 away must not be split or welded"
        );
        assert_eq!(store[0][1][0].absolute_front().id(), c.id());
        assert_eq!(store[0][1][0].absolute_back().id(), d.id());
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

        assert!(
            store
                .add_imprint_vertex(0, &split, BooleanTolerance::new(KERNEL_EPSILON), &mut emap)
                .is_none()
        );
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
        let point = Point3::new(0.5, 1.125 * KERNEL_EPSILON, 0.0);

        assert!(
            store
                .search_boundary_parameter(
                    0,
                    point,
                    8.0 * KERNEL_EPSILON,
                    BooleanTolerance::new(KERNEL_EPSILON)
                )
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
        let point = Point3::new(1.0, 1.125 * KERNEL_EPSILON, 0.0);

        let (_, _, kind, projected) = store
            .search_boundary_parameter(
                0,
                point,
                8.0 * KERNEL_EPSILON,
                BooleanTolerance::new(KERNEL_EPSILON),
            )
            .expect("incident edges sharing one vertex are one boundary locus");

        assert!(matches!(kind, ParameterKind::Front | ParameterKind::Back));
        assert_eq!(projected, Some(shared.point()));
    }

    #[test]
    fn loop_boundary_search_rejects_ambiguous_bounded_nearest_projection() {
        let a0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let a1 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let b0 = Vertex::new(Point3::new(0.0, 2.25 * KERNEL_EPSILON, 0.0));
        let b1 = Vertex::new(Point3::new(1.0, 2.25 * KERNEL_EPSILON, 0.0));
        let edge_a = Edge::new_unchecked(&a0, &a1, Line(a0.point(), a1.point()));
        let edge_b = Edge::new_unchecked(&b0, &b1, Line(b0.point(), b1.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![edge_a].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge_b].into(), ShapesOpStatus::Unknown),
        ]);
        let store = LoopsStore(vec![loops]);
        let point = Point3::new(0.5, 1.125 * KERNEL_EPSILON, 0.0);

        assert!(
            store
                .search_boundary_parameter(
                    0,
                    point,
                    8.0 * KERNEL_EPSILON,
                    BooleanTolerance::new(KERNEL_EPSILON)
                )
                .is_none()
        );
    }

    #[test]
    fn nearby_distinct_boundary_edges_cannot_pass_global_coincidence_epsilon() {
        let gap = 5.0e-7;
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let c = Vertex::new(Point3::new(0.0, 0.0, gap));
        let d = Vertex::new(Point3::new(1.0, 0.0, gap));
        let edge0 = Edge::new_unchecked(&a, &b, Line(a.point(), b.point()));
        let edge1 = Edge::new_unchecked(&c, &d, Line(c.point(), d.point()));
        let store = LoopsStore(vec![Loops(vec![
            BoundaryWire::new(vec![edge0].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge1].into(), ShapesOpStatus::Unknown),
        ])]);
        let query = Point3::new(0.5, 0.0, gap / 2.0);
        let tolerances = BooleanTolerance::new(1.0e-8);
        // Both are plausible nearest-search candidates. Neither is exact,
        // and they do not constitute one topological locus at 1e-8.
        assert!(
            store
                .search_boundary_parameter(0, query, 1.0e-3, tolerances)
                .is_none(),
            "two globally-near but operation-distinct trim edges are ambiguous"
        );
    }

    #[test]
    fn sampled_endpoint_selects_unique_exact_vertex_on_opposite_surface() {
        let near = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let far = Vertex::new(Point3::new(1.0, 0.0003, 0.0));
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(2.0, 0.0003, 0.0));
        let edge0 = Edge::new_unchecked(&a, &near, Line(a.point(), near.point()));
        let edge1 = Edge::new_unchecked(&far, &b, Line(far.point(), b.point()));
        let loops = Loops(vec![
            BoundaryWire::new(vec![edge0].into(), ShapesOpStatus::Unknown),
            BoundaryWire::new(vec![edge1].into(), ShapesOpStatus::Unknown),
        ]);
        let store = LoopsStore(vec![loops]);
        let selected = store
            .supported_boundary_vertex(
                0,
                Point3::new(1.0, 0.00002, 0.0),
                0.001,
                BooleanTolerance::new(KERNEL_EPSILON),
                &Plane::xy(),
            )
            .expect("one supported exact vertex is decisively nearer");
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].id(), near.id());
    }

    #[test]
    fn sampled_endpoint_rejects_ambiguous_or_unsupported_exact_vertices() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(0.0, 0.000012, 0.0));
        let far0 = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let far1 = Vertex::new(Point3::new(1.0, 0.000012, 0.0));
        let loops = Loops(vec![
            BoundaryWire::new(
                vec![Edge::new_unchecked(
                    &a,
                    &far0,
                    Line(a.point(), far0.point()),
                )]
                .into(),
                ShapesOpStatus::Unknown,
            ),
            BoundaryWire::new(
                vec![Edge::new_unchecked(
                    &b,
                    &far1,
                    Line(b.point(), far1.point()),
                )]
                .into(),
                ShapesOpStatus::Unknown,
            ),
        ]);
        let store = LoopsStore(vec![loops]);
        assert!(
            store
                .supported_boundary_vertex(
                    0,
                    Point3::new(0.0, 0.000006, 0.0),
                    0.001,
                    BooleanTolerance::new(KERNEL_EPSILON),
                    &Plane::xy()
                )
                .is_none(),
            "two similarly close distinct vertices are ambiguous"
        );

        let lifted = LoopsStore(vec![Loops(vec![BoundaryWire::new(
            vec![Edge::new_unchecked(
                &Vertex::new(Point3::new(0.0, 0.0, 0.02)),
                &Vertex::new(Point3::new(1.0, 0.0, 0.02)),
                Line(Point3::new(0.0, 0.0, 0.02), Point3::new(1.0, 0.0, 0.02)),
            )]
            .into(),
            ShapesOpStatus::Unknown,
        )])]);
        assert!(
            lifted
                .supported_boundary_vertex(
                    0,
                    Point3::new(0.0, 0.0, 0.02),
                    0.001,
                    BooleanTolerance::new(KERNEL_EPSILON),
                    &Plane::xy()
                )
                .is_none(),
            "vertex off opposite carrier cannot be reused"
        );
    }

    #[test]
    fn cross_face_stitch_respects_sub_epsilon_numerical_length() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let edge = Edge::new_unchecked(&a, &b, Line(a.point(), b.point()));
        let lhs = Loops(vec![BoundaryWire::new(
            vec![edge].into(),
            ShapesOpStatus::Unknown,
        )]);
        let rhs = Loops(Vec::new());
        let near = Point3::new(0.0, 0.0, 5.0e-7);
        assert_eq!(
            unique_nearby_loop_point(&lhs, &rhs, near, BooleanTolerance::new(1.0e-8)),
            None,
            "a tight Boolean may not canonicalize different endpoints"
        );
        assert_eq!(
            unique_nearby_loop_point(&lhs, &rhs, near, BooleanTolerance::new(1.0e-6)),
            Some(a.point()),
            "a larger numerical budget admits the same candidate"
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

        let kernel_near = Point3::new(0.0, 0.0, 4.0 * KERNEL_EPSILON);
        assert_eq!(
            unique_nearby_loop_point(
                &lhs,
                &rhs,
                kernel_near,
                BooleanTolerance::new(KERNEL_EPSILON)
            ),
            Some(a.point())
        );

        // A loose user-facing modeling budget must never scale this kernel-level
        // stitch envelope: z=0.8 and z=1.0 topology remain distinct endpoints.
        let model_scale_distinct = Point3::new(0.0, 0.0, 0.2);
        assert!(
            unique_nearby_loop_point(
                &lhs,
                &rhs,
                model_scale_distinct,
                BooleanTolerance::new(KERNEL_EPSILON)
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
                .ok()
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

        let _ = loops.add_edge(
            incoming,
            ShapesOpStatus::And,
            false,
            BooleanTolerance::new(KERNEL_EPSILON),
        );

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
    fn imprint_reuses_wiggled_same_edge_within_explicit_boolean_tolerance() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let c = Vertex::new(Point3::new(1.0, 1.0, 0.0));
        let d = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let smooth = Edge::new_unchecked(
            &a,
            &b,
            PolylineCurve(vec![
                a.point(),
                Point3::new(0.25, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.0),
                Point3::new(0.75, 0.0, 0.0),
                b.point(),
            ]),
        );
        let line = |front: &Vertex<Point3>, back: &Vertex<Point3>| {
            Edge::new_unchecked(
                front,
                back,
                PolylineCurve(vec![front.point(), back.point()]),
            )
        };
        let mut loops = Loops(vec![BoundaryWire::new(
            vec![smooth.clone(), line(&b, &c), line(&c, &d), line(&d, &a)].into(),
            ShapesOpStatus::Or,
        )]);
        let a2 = Vertex::new(Point3::new(2.0e-7, 1.0e-7, 0.0));
        let b2 = Vertex::new(Point3::new(1.0 + 3.0e-7, 2.0e-7, 0.0));
        let incoming = Edge::new_unchecked(
            &a2,
            &b2,
            PolylineCurve(vec![
                a2.point(),
                Point3::new(0.25, -4.0e-7, 0.0),
                Point3::new(0.5, 2.0e-6, 0.0),
                Point3::new(0.75, -5.0e-7, 0.0),
                b2.point(),
            ]),
        );
        assert!(
            loops.reuse_imprinted_coincident_edge(&incoming, 1.0e-5),
            "a tolerance-supported coincident source edge must be reused"
        );
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].len(), 4);
        assert!(
            loops[0].iter().any(|edge| edge.id() == incoming.id()),
            "the existing target uses must share the canonical incoming EdgeId"
        );
        assert!(
            !Loops(vec![BoundaryWire::new(
                vec![smooth, line(&b, &c), line(&c, &d), line(&d, &a)].into(),
                ShapesOpStatus::Or,
            )])
            .reuse_imprinted_coincident_edge(&incoming, 1.0e-6),
            "the same noise is distinguishable under a tighter Boolean tolerance"
        );
    }

    #[test]
    fn imprint_coincidence_does_not_require_matching_curve_parameters() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let existing = Edge::new_unchecked(
            &a,
            &b,
            PolylineCurve(vec![
                a.point(),
                Point3::new(0.25, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.0),
                Point3::new(0.75, 0.0, 0.0),
                b.point(),
            ]),
        );
        let incoming = Edge::new_unchecked(
            &a,
            &b,
            PolylineCurve(vec![
                a.point(),
                Point3::new(0.02, 0.0, 0.0),
                Point3::new(0.11, 0.0, 0.0),
                Point3::new(0.42, 0.0, 0.0),
                Point3::new(0.93, 0.0, 0.0),
                b.point(),
            ]),
        );
        let mut loops = Loops(vec![BoundaryWire::new(
            vec![existing].into(),
            ShapesOpStatus::Unknown,
        )]);
        assert!(
            loops.reuse_imprinted_coincident_edge(&incoming, 1.0e-6),
            "spatially identical curves with unrelated knot parameters must share an imprint edge"
        );
        assert_eq!(loops[0][0].id(), incoming.id());
    }

    #[test]
    fn imprint_keeps_real_thin_lenses_and_matching_endpoints_with_distinct_paths() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let start = Edge::new_unchecked(&a, &b, PolylineCurve(vec![a.point(), b.point()]));
        let candidate = |offsets: &[f64]| {
            Edge::new_unchecked(
                &a,
                &b,
                PolylineCurve(
                    offsets
                        .iter()
                        .enumerate()
                        .map(|(index, &y)| {
                            Point3::new(index as f64 / (offsets.len() - 1) as f64, y, 0.0)
                        })
                        .collect(),
                ),
            )
        };
        let mut existing = Loops(vec![BoundaryWire::new(
            vec![start].into(),
            ShapesOpStatus::Unknown,
        )]);
        assert!(
            !existing.reuse_imprinted_coincident_edge(&candidate(&[0.0, 5.0e-5, 0.0]), 1.0e-5,),
            "a nonzero thin lens beyond requested tolerance cannot collapse"
        );
        assert!(
            !existing.reuse_imprinted_coincident_edge(
                &candidate(&[0.0, 4.0e-5, 0.0, -4.0e-5, 0.0]),
                1.0e-5,
            ),
            "matching endpoints and midpoint must not hide a distinct path"
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
    fn tangent_boundary_projection_rejects_distinct_sub_epsilon_supports() {
        // This line and plane are parallel and genuinely 5e-7 apart.
        // Their projections coincide only if the caller explicitly permits
        // that much exact support error. A singular tangent update cannot
        // serve as an alternative proof of incidence.
        let z = 5.0e-7;
        let line = Line(Point3::new(0.0, 0.0, z), Point3::new(1.0, 0.0, z));
        let seed = Point3::new(0.5, 0.0, z);
        assert!(
            curve_surface_projection(
                &line,
                None,
                &Plane::xy(),
                None,
                seed,
                BooleanTolerance::new(1.0e-8),
                100,
            )
            .is_none(),
            "the exact supports are distinct at 1e-8 numerical length"
        );
        assert!(
            curve_surface_projection(
                &line,
                None,
                &Plane::xy(),
                None,
                seed,
                BooleanTolerance::new(1.0e-6),
                100,
            )
            .is_some(),
            "the same projected pair may pass a deliberately wider budget"
        );
    }

    #[test]
    fn curve_surface_projection_accepts_agreeing_tangent_projections() {
        let line = Line(Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0));
        let sphere = Sphere::new(Point3::origin(), 1.0);
        let seed = Point3::new(1.0 + 1.5 * KERNEL_EPSILON, 0.0, 0.0);

        let (point, parameter, _) = curve_surface_projection(
            &line,
            None,
            &sphere,
            None,
            seed,
            BooleanTolerance::new(KERNEL_EPSILON),
            100,
        )
        .expect("agreeing exact projections should resolve a tangent SSI seed");

        assert!(point.near(&Point3::new(1.0, 0.0, 0.0)));
        assert!((parameter - 0.5).abs() <= KERNEL_EPSILON);
        assert!(sphere.include(point));
    }
}
