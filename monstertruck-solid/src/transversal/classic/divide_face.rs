//! Classic (0.3.2) face division.
//!
//! Ported verbatim from the published 0.3.2 crate's `transversal::divide_face`:
//! projects each loop's edges into the face parameter domain, splits the face
//! into positively oriented pre-faces with their contained holes, and tags each
//! with its `and`/`or`/`unknown` status.

#![allow(clippy::many_single_char_names)]

use super::faces_classification::FacesClassification;
use super::loops_store::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::*;
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
        if area.abs() < tol {
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
        let sibling = pre_faces.iter().position(|pre_face| {
            if pre_face.is_empty() {
                return false;
            }
            let outer = &pre_face[0];
            let complementary = matches!(
                (outer.wire.status(), chunk.wire.status()),
                (ShapesOpStatus::And, ShapesOpStatus::Or)
                    | (ShapesOpStatus::Or, ShapesOpStatus::And)
            );
            let same_edge_set = outer.wire.len() == chunk.wire.len()
                && outer
                    .wire
                    .edge_iter()
                    .all(|edge| chunk.wire.edge_iter().any(|other| edge.id() == other.id()));
            complementary
                && !same_edge_set
                && outer
                    .wire
                    .vertex_iter()
                    .any(|vertex| chunk.wire.vertex_iter().any(|other| vertex == other))
        });
        if sibling.is_some() {
            chunk.poly.invert();
            chunk.reverse_for_face = true;
            pre_faces.push(vec![chunk]);
            return Some(());
        }

        let pt = chunk.poly.front();
        let chunk_area_abs = chunk.poly.area().abs();
        // A negative loop may be nested inside several positive loops. The
        // equal-area opposite-status positive loop is the sibling region on
        // the other side of the same divider, not its parent. Prefer the
        // smallest *strictly larger* enclosing loop; fall back to an equal
        // enclosing loop only when no larger parent exists, preserving the
        // whole-face cancellation case.
        let idx = pre_faces
            .iter()
            .enumerate()
            .filter(|(_, pre_face)| {
                !pre_face.is_empty()
                    && pre_face[0].poly.include(pt)
                    && pre_face[0].poly.area().abs() + tol >= chunk_area_abs
            })
            .min_by(|(_, lhs), (_, rhs)| {
                let lhs_area = lhs[0].poly.area().abs();
                let rhs_area = rhs[0].poly.area().abs();
                let lhs_equal = lhs_area <= chunk_area_abs + tol;
                let rhs_equal = rhs_area <= chunk_area_abs + tol;
                lhs_equal
                    .cmp(&rhs_equal)
                    .then_with(|| lhs_area.total_cmp(&rhs_area))
            })
            .map(|(index, _)| index);
        if let Some(i) = idx {
            let outer_area = pre_faces[i][0].poly.area();
            let chunk_area = chunk.poly.area();
            // If the sum of areas is zero, the face is canceled.
            // This happens when an intersection loop exactly matches the face boundary.
            if (outer_area + chunk_area).abs() < tol {
                pre_faces[i].clear();
            } else {
                pre_faces[i].push(chunk);
            }
        }
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
            let mut new_face = Face::debug_new(wires, surface).ok()?;
            if !face.orientation() {
                new_face.invert();
            }
            Some((new_face, status))
        })
        .collect::<Option<Vec<_>>>()?;
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
