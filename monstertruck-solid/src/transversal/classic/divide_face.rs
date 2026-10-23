//! Face division for the classic Boolean backend.
//!
//! The implementation descends from the published 0.3.2 kernel, but split regions
//! are now arranged by geometric nesting rather than status/shared-vertex cases.
//! Intermediate regions remain `FaceFragment`s until compressed-topology repair and
//! extraction establish valid [`Face`] invariants.

#![allow(clippy::many_single_char_names)]

use super::faces_classification::FacesClassification;
use super::{loops_store::*, tolerance::BooleanTolerance};
use monstertruck_meshing::prelude::*;
use monstertruck_topology::{compress::*, repair, *};
use rustc_hash::FxHashMap as HashMap;
use std::ops::Deref;

fn project_to_parameter<S>(
    surface: &S,
    point: Point3,
    hint: Option<Point2>,
    tol: f64,
) -> Option<Point2>
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let hint = hint.map(|uv| (uv.x, uv.y));
    if let Some(parameter) = surface.search_parameter(point, hint, 100) {
        return Some(parameter.into());
    }

    let accept_nearest = |parameter: (f64, f64)| {
        let projected = surface.evaluate(parameter.0, parameter.1);
        (projected.distance2(point).sqrt() <= tol).then(|| parameter.into())
    };

    if let Some(parameter) = surface.search_nearest_parameter(point, hint, 100)
        && let Some(parameter) = accept_nearest(parameter)
    {
        return Some(parameter);
    }

    if hint.is_some() {
        if let Some(parameter) = surface.search_parameter(point, None, 100) {
            return Some(parameter.into());
        }
        if let Some(parameter) = surface.search_nearest_parameter(point, None, 100)
            && let Some(parameter) = accept_nearest(parameter)
        {
            return Some(parameter);
        }
    }

    None
}

fn create_parameter_boundary<C, S>(
    face: &Face<Point3, C, S>,
    wire: &Wire<Point3, C>,
    polys: &mut HashMap<EdgeId<C>, PolylineCurve<Point3>>,
    tol: f64,
) -> Option<PolylineCurve<Point2>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let surface = face.surface();
    let pt = wire.front_vertex().unwrap().point();
    let p = project_to_parameter(&surface, pt, None, tol)?;
    let vec = wire.edge_iter().try_fold(vec![p], |mut vec, edge| {
        let poly = polys.entry(edge.id()).or_insert_with(|| {
            let curve = edge.curve();
            let div = curve.parameter_division(curve.range_tuple(), tol).1;
            PolylineCurve(div)
        });
        let mut p = *vec.last().unwrap();
        let closure = |q: &Point3| -> Option<Point2> {
            p = project_to_parameter(&surface, *q, Some(p), tol)?;
            Some(p)
        };
        let add: Option<Vec<Point2>> = match edge.orientation() {
            true => poly.iter().skip(1).map(closure).collect(),
            false => poly.iter().rev().skip(1).map(closure).collect(),
        };
        vec.append(&mut add?);
        Some(vec)
    })?;
    Some(PolylineCurve(vec))
}

#[derive(Clone, Debug)]
struct WireChunk<'a, C> {
    poly: PolylineCurve<Point2>,
    wire: &'a BoundaryWire<Point3, C>,
    reverse_for_face: bool,
}

type FaceWithShapesOpStatus<C, S> = (Face<Point3, C, S>, ShapesOpStatus);
#[derive(Clone, Debug)]
struct FaceFragment<C, S> {
    wires: Vec<Wire<Point3, C>>,
    surface: S,
    orientation: bool,
}

impl<C: Clone, S: Clone> FaceFragment<C, S> {
    fn canonicalize(self) -> Option<Vec<Face<Point3, C, S>>> {
        let mut vertices = Vec::<Point3>::new();
        let mut vertex_map = HashMap::<VertexId<Point3>, usize>::default();
        let mut edges = Vec::<CompressedEdge<C>>::new();
        let mut edge_map = HashMap::<EdgeId<C>, usize>::default();

        let mut vertex_index = |vertex: &Vertex<Point3>| {
            if let Some(index) = vertex_map.get(&vertex.id()).copied() {
                index
            } else {
                let index = vertices.len();
                vertices.push(vertex.point());
                vertex_map.insert(vertex.id(), index);
                index
            }
        };

        let boundaries = self
            .wires
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

        let mut compressed = CompressedShell {
            vertices,
            edges,
            faces: vec![CompressedFace {
                boundaries,
                orientation: self.orientation,
                surface: self.surface,
            }],
            vertex_stable_ids: None,
            edge_stable_ids: None,
            face_stable_ids: None,
        };
        repair::split_non_simple_compressed_wires(&mut compressed);
        repair::split_pinched_compressed_faces(&mut compressed);
        Shell::extract(compressed)
            .ok()
            .map(IntoIterator::into_iter)
            .map(Iterator::collect)
    }
}

fn strict_interior_point(poly: &PolylineCurve<Point2>, tol: f64) -> Option<Point2> {
    let area = poly.area();
    if area.abs() <= polyline_area_tolerance(poly, tol) {
        return None;
    }
    let mut domain = poly.clone();
    if area < 0.0 {
        domain.invert();
    }
    let points = &poly.0;
    if points.len() < 2 {
        return None;
    }
    for index in 0..points.len() {
        let a = points[index];
        let b = points[(index + 1) % points.len()];
        let edge = b - a;
        let length = edge.magnitude();
        if length <= tol {
            continue;
        }
        let left = Vector2::new(-edge.y, edge.x) / length;
        let inward = if area > 0.0 { left } else { -left };
        let midpoint = a.midpoint(b);
        for scale in [4.0, 16.0, 64.0] {
            let candidate = midpoint + inward * (tol * scale).min(length * 0.125);
            if domain.include(candidate) {
                return Some(candidate);
            }
        }
    }
    None
}

fn same_edge_set<C>(lhs: &BoundaryWire<Point3, C>, rhs: &BoundaryWire<Point3, C>) -> bool {
    lhs.len() == rhs.len()
        && lhs
            .edge_iter()
            .all(|edge| rhs.edge_iter().any(|other| edge.id() == other.id()))
}

fn polyline_area_tolerance(poly: &PolylineCurve<Point2>, tol: f64) -> f64 {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for point in poly.iter() {
        min_x = min_x.min(point.x);
        max_x = max_x.max(point.x);
        min_y = min_y.min(point.y);
        max_y = max_y.max(point.y);
    }
    let span = (max_x - min_x).abs().max((max_y - min_y).abs());
    BooleanTolerance::new(tol).area_for_span(span)
}

fn divide_one_face<C, S>(
    face: &Face<Point3, C, S>,
    loops: &Loops<Point3, C>,
    tol: f64,
) -> Option<Vec<FaceWithShapesOpStatus<C, S>>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let (mut pre_faces, mut negative_wires) = (Vec::new(), Vec::new());
    let mut map = HashMap::default();
    loops.iter().try_for_each(|wire| {
        let poly = create_parameter_boundary(face, wire, &mut map, tol)?;
        let area = poly.area();
        if area.abs() < polyline_area_tolerance(&poly, tol) {
            return Some(());
        }
        match area > 0.0 {
            true => pre_faces.push(vec![WireChunk {
                poly,
                wire,
                reverse_for_face: false,
            }]),
            false => negative_wires.push(WireChunk {
                poly,
                wire,
                reverse_for_face: false,
            }),
        }
        Some(())
    })?;
    negative_wires.into_iter().try_for_each(|mut chunk| {
        let witness = strict_interior_point(&chunk.poly, tol)?;
        let chunk_area_abs = chunk.poly.area().abs();

        // Build the nesting relation from geometry, not Boolean status or shared
        // vertices. A negative loop is a hole of the smallest strictly larger
        // positive region containing an interior witness. Boundary contact alone
        // never establishes containment.
        let parent = pre_faces
            .iter()
            .enumerate()
            .filter(|(_, pre_face)| !pre_face.is_empty())
            .filter(|(_, pre_face)| {
                let outer = &pre_face[0];
                let area_tol = polyline_area_tolerance(&outer.poly, tol)
                    .max(polyline_area_tolerance(&chunk.poly, tol));
                outer.poly.area().abs() > chunk_area_abs + area_tol && outer.poly.include(witness)
            })
            .min_by(|(_, lhs), (_, rhs)| {
                lhs[0]
                    .poly
                    .area()
                    .abs()
                    .total_cmp(&rhs[0].poly.area().abs())
            })
            .map(|(index, _)| index);

        if let Some(index) = parent {
            pre_faces[index].push(chunk);
            return Some(());
        }

        // With no strictly larger parent, an equal inverse boundary cancels the
        // region. This is intentionally checked after nesting so a duplicated
        // inner cut becomes "outer-with-hole + inner disk" rather than deleting
        // the inner disk.
        let equal = pre_faces.iter().position(|pre_face| {
            if pre_face.is_empty() {
                return false;
            }
            let outer = &pre_face[0];
            let area_tol = polyline_area_tolerance(&outer.poly, tol)
                .max(polyline_area_tolerance(&chunk.poly, tol));
            (outer.poly.area().abs() - chunk_area_abs).abs() <= area_tol
                && same_edge_set(outer.wire, chunk.wire)
        });
        if let Some(index) = equal {
            pre_faces[index].clear();
            return Some(());
        }

        // A clockwise region that is not contained by another positive region
        // is its own sibling cell. Flip it into the canonical positive
        // orientation; no status-specific special case is required.
        chunk.poly.invert();
        chunk.reverse_for_face = true;
        pre_faces.push(vec![chunk]);
        Some(())
    })?;
    let vec: Vec<_> = pre_faces
        .into_iter()
        .filter(|pre_face| !pre_face.is_empty())
        .map(|pre_face| {
            let surface = face.surface();
            let op = pre_face
                .iter()
                .find(|chunk| chunk.wire.status() != ShapesOpStatus::Unknown);
            let status = match op {
                Some(chunk) => chunk.wire.status(),
                None => ShapesOpStatus::Unknown,
            };
            let wires: Vec<Wire<Point3, C>> = pre_face
                .into_iter()
                .map(|chunk| {
                    let wire = chunk.wire.deref().clone();
                    if chunk.reverse_for_face {
                        wire.inverse()
                    } else {
                        wire
                    }
                })
                .collect();
            let fragment = FaceFragment {
                wires,
                surface,
                orientation: face.orientation(),
            };
            fragment.canonicalize().map(|faces| {
                faces
                    .into_iter()
                    .map(|face| (face, status))
                    .collect::<Vec<_>>()
            })
        })
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    Some(vec)
}

pub(super) fn divide_faces<C, S>(
    shell: &Shell<Point3, C, S>,
    loops_store: &LoopsStore<Point3, C>,
    tol: f64,
) -> Option<FacesClassification<Point3, C, S>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let mut res = FacesClassification::<Point3, C, S>::default();
    shell
        .iter()
        .zip(loops_store)
        .try_for_each(|(face, loops)| {
            if loops
                .iter()
                .all(|wire| wire.status() == ShapesOpStatus::Unknown)
            {
                res.push(face.clone(), ShapesOpStatus::Unknown);
            } else {
                let vec = divide_one_face(face, loops, tol)?;
                vec.into_iter()
                    .for_each(|(face, status)| res.push(face, status));
            }
            Some(())
        })?;
    Some(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::{BsplineCurve, KnotVector, Plane};

    #[test]
    fn approximate_surface_projection_is_bounded_by_boolean_tolerance() {
        let plane = Plane::xy();
        let point = Point3::new(0.25, 0.75, 0.01);

        let parameter = project_to_parameter(&plane, point, None, 0.02)
            .expect("a nearest projection inside the Boolean tolerance must be accepted");
        assert!(parameter.near(&Point2::new(0.25, 0.75)));

        assert!(
            project_to_parameter(&plane, point, None, 0.005).is_none(),
            "a nearest projection outside the Boolean tolerance must be rejected"
        );
    }

    fn test_line(
        front: &Vertex<Point3>,
        back: &Vertex<Point3>,
    ) -> Edge<Point3, BsplineCurve<Point3>> {
        Edge::new(
            front,
            back,
            BsplineCurve::new(
                KnotVector::bezier_knot(1),
                vec![front.point(), back.point()],
            ),
        )
    }

    #[test]
    fn touching_complementary_loops_are_sibling_faces() {
        let l0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let a = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let r0 = Vertex::new(Point3::new(2.0, 0.0, 0.0));
        let l1 = Vertex::new(Point3::new(0.0, 2.0, 0.0));
        let b = Vertex::new(Point3::new(1.0, 2.0, 0.0));
        let r1 = Vertex::new(Point3::new(2.0, 2.0, 0.0));

        let outer: Wire<_, _> = vec![
            test_line(&l0, &r0),
            test_line(&r0, &r1),
            test_line(&r1, &l1),
            test_line(&l1, &l0),
        ]
        .into();
        let face = Face::debug_new(vec![outer], Plane::xy()).expect("valid source face");

        let left: Wire<_, _> = vec![
            test_line(&l0, &a),
            test_line(&a, &b),
            test_line(&b, &l1),
            test_line(&l1, &l0),
        ]
        .into();
        let right_positive: Wire<_, _> = vec![
            test_line(&a, &r0),
            test_line(&r0, &r1),
            test_line(&r1, &b),
            test_line(&b, &a),
        ]
        .into();
        let loops: Loops<_, _> = vec![
            BoundaryWire::new(left, ShapesOpStatus::And),
            BoundaryWire::new(right_positive.inverse(), ShapesOpStatus::Or),
        ]
        .into_iter()
        .collect();

        let divided = divide_one_face(&face, &loops, 0.01).expect("split must succeed");
        assert_eq!(divided.len(), 2);
        assert!(divided.iter().all(|(face, _)| face.boundaries().len() == 1));
        assert!(
            divided
                .iter()
                .any(|(_, status)| *status == ShapesOpStatus::And)
        );
        assert!(
            divided
                .iter()
                .any(|(_, status)| *status == ShapesOpStatus::Or)
        );
    }

    #[test]
    fn inverse_copy_of_closed_cut_still_partitions_outer_face_and_inner_disk() {
        let o0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let o1 = Vertex::new(Point3::new(4.0, 0.0, 0.0));
        let o2 = Vertex::new(Point3::new(4.0, 4.0, 0.0));
        let o3 = Vertex::new(Point3::new(0.0, 4.0, 0.0));
        let h0 = Vertex::new(Point3::new(1.0, 1.0, 0.0));
        let h1 = Vertex::new(Point3::new(3.0, 1.0, 0.0));
        let h2 = Vertex::new(Point3::new(3.0, 3.0, 0.0));
        let h3 = Vertex::new(Point3::new(1.0, 3.0, 0.0));

        let outer: Wire<_, _> = vec![
            test_line(&o0, &o1),
            test_line(&o1, &o2),
            test_line(&o2, &o3),
            test_line(&o3, &o0),
        ]
        .into();
        let face = Face::debug_new(vec![outer.clone()], Plane::xy()).expect("valid source face");
        let inner: Wire<_, _> = vec![
            test_line(&h0, &h1),
            test_line(&h1, &h2),
            test_line(&h2, &h3),
            test_line(&h3, &h0),
        ]
        .into();
        let loops: Loops<_, _> = vec![
            BoundaryWire::new(outer, ShapesOpStatus::Unknown),
            BoundaryWire::new(inner.clone(), ShapesOpStatus::Or),
            BoundaryWire::new(inner.inverse(), ShapesOpStatus::And),
        ]
        .into_iter()
        .collect();

        let divided = divide_one_face(&face, &loops, 0.01).expect("closed cut must partition");
        assert_eq!(divided.len(), 2);
        assert!(divided.iter().any(|(face, _)| face.boundaries().len() == 2));
        assert!(divided.iter().any(|(face, _)| face.boundaries().len() == 1));
    }

    #[test]
    fn disjoint_negative_loop_remains_a_hole() {
        let o0 = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let o1 = Vertex::new(Point3::new(4.0, 0.0, 0.0));
        let o2 = Vertex::new(Point3::new(4.0, 4.0, 0.0));
        let o3 = Vertex::new(Point3::new(0.0, 4.0, 0.0));
        let h0 = Vertex::new(Point3::new(1.0, 1.0, 0.0));
        let h1 = Vertex::new(Point3::new(3.0, 1.0, 0.0));
        let h2 = Vertex::new(Point3::new(3.0, 3.0, 0.0));
        let h3 = Vertex::new(Point3::new(1.0, 3.0, 0.0));

        let outer: Wire<_, _> = vec![
            test_line(&o0, &o1),
            test_line(&o1, &o2),
            test_line(&o2, &o3),
            test_line(&o3, &o0),
        ]
        .into();
        let face = Face::debug_new(vec![outer.clone()], Plane::xy()).expect("valid source face");
        let hole_positive: Wire<_, _> = vec![
            test_line(&h0, &h1),
            test_line(&h1, &h2),
            test_line(&h2, &h3),
            test_line(&h3, &h0),
        ]
        .into();
        let loops: Loops<_, _> = vec![
            BoundaryWire::new(outer, ShapesOpStatus::And),
            BoundaryWire::new(hole_positive.inverse(), ShapesOpStatus::Unknown),
        ]
        .into_iter()
        .collect();

        let divided = divide_one_face(&face, &loops, 0.01).expect("hole nesting must succeed");
        assert_eq!(divided.len(), 1);
        assert_eq!(divided[0].0.boundaries().len(), 2);
        assert_eq!(divided[0].1, ShapesOpStatus::And);
    }

    #[test]
    fn exact_surface_projection_does_not_depend_on_boolean_tolerance() {
        let plane = Plane::xy();
        let point = Point3::new(0.25, 0.75, 0.0);

        let parameter = project_to_parameter(&plane, point, None, f64::EPSILON)
            .expect("an exact surface point must use the exact inverse map");
        assert!(parameter.near(&Point2::new(0.25, 0.75)));
    }
}
