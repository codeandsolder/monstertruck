//! Face division for the classic Boolean backend.
//!
//! The implementation descends from the published 0.3.2 kernel, but split regions
//! are now arranged by geometric nesting rather than status/shared-vertex cases.
//! Intermediate regions remain `FaceFragment`s until compressed-topology repair and
//! extraction establish valid [`Face`] invariants.

#![allow(clippy::many_single_char_names)]

use super::faces_classification::FacesClassification;
use super::loops_store::*;
use super::tolerance::BooleanTolerance;
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

fn create_parameter_boundary_on_surface<C, S>(
    surface: &S,
    wire: &Wire<Point3, C>,
    sampled_wire: Option<&Wire<Point3, PolylineCurve<Point3>>>,
    polys: &mut HashMap<EdgeId<C>, PolylineCurve<Point3>>,
    tol: f64,
) -> Option<PolylineCurve<Point2>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let sampling_tol = BooleanTolerance::new(tol).mesh_seed;
    let projection_tol = if sampled_wire.is_some() {
        sampling_tol
    } else {
        tol
    };
    let pt = wire.front_vertex()?.point();
    let p = project_to_parameter(surface, pt, None, projection_tol)?;
    let mut vec = vec![p];
    if let Some(sampled) = sampled_wire {
        for edge in sampled.iter() {
            let poly = edge.curve();
            let mut previous = vec.last().copied()?;
            let mut project = |q: &Point3| -> Option<Point2> {
                previous = project_to_parameter(surface, *q, Some(previous), sampling_tol)?;
                Some(previous)
            };
            let mut add = match edge.orientation() {
                true => poly
                    .iter()
                    .skip(1)
                    .map(&mut project)
                    .collect::<Option<Vec<_>>>()?,
                false => poly
                    .iter()
                    .rev()
                    .skip(1)
                    .map(&mut project)
                    .collect::<Option<Vec<_>>>()?,
            };
            vec.append(&mut add);
        }
        return Some(PolylineCurve(vec));
    }

    for edge in wire.iter() {
        let poly = polys.entry(edge.id()).or_insert_with(|| {
            let curve = edge.curve();
            PolylineCurve(
                curve
                    .parameter_division(curve.range_tuple(), sampling_tol)
                    .1,
            )
        });
        let mut previous = vec.last().copied()?;
        let mut project = |q: &Point3| -> Option<Point2> {
            previous = project_to_parameter(surface, *q, Some(previous), tol)?;
            Some(previous)
        };
        let mut add = match edge.orientation() {
            true => poly
                .iter()
                .skip(1)
                .map(&mut project)
                .collect::<Option<Vec<_>>>()?,
            false => poly
                .iter()
                .rev()
                .skip(1)
                .map(&mut project)
                .collect::<Option<Vec<_>>>()?,
        };
        vec.append(&mut add);
    }
    Some(PolylineCurve(vec))
}

fn create_parameter_boundary<C, S>(
    face: &Face<Point3, C, S>,
    wire: &Wire<Point3, C>,
    sampled_wire: Option<&Wire<Point3, PolylineCurve<Point3>>>,
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
    create_parameter_boundary_on_surface(&face.surface(), wire, sampled_wire, polys, tol)
}

#[derive(Clone, Debug)]
struct WireChunk<'a, C> {
    poly: PolylineCurve<Point2>,
    wire: &'a BoundaryWire<Point3, C>,
    reverse_for_face: bool,
}

#[derive(Clone, Debug)]
struct FaceFragment<C, S> {
    wires: Vec<Wire<Point3, C>>,
    surface: S,
    orientation: bool,
}

#[cfg(test)]
impl<C: Clone, S> FaceFragment<C, S> {
    fn boundaries(&self) -> Vec<Wire<Point3, C>> {
        if self.orientation {
            self.wires.clone()
        } else {
            self.wires.iter().map(Wire::inverse).collect()
        }
    }
}

fn polyline_interior_candidates(poly: &PolylineCurve<Point2>) -> Vec<Point2> {
    let area = polyline_signed_area(poly);
    if area.abs() <= polyline_area_epsilon(poly) {
        return Vec::new();
    }
    let mut domain = poly.clone();
    if area < 0.0 {
        domain.invert();
    }
    let points = &poly.0;
    if points.len() < 2 {
        return Vec::new();
    }
    let (mut min_x, mut max_x, mut min_y, mut max_y) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for point in points {
        min_x = min_x.min(point.x);
        max_x = max_x.max(point.x);
        min_y = min_y.min(point.y);
        max_y = max_y.max(point.y);
    }
    let scale = (max_x - min_x).abs().max((max_y - min_y).abs());
    if !scale.is_finite() || scale <= f64::EPSILON {
        return Vec::new();
    }
    let mut candidates = Vec::new();
    for index in 0..points.len() {
        let a = points[index];
        let b = points[(index + 1) % points.len()];
        let edge = b - a;
        let length = edge.magnitude();
        if length <= f64::EPSILON {
            continue;
        }
        let left = Vector2::new(-edge.y, edge.x) / length;
        let inward = if area > 0.0 { left } else { -left };
        let midpoint = a.midpoint(b);
        for relative in [1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2] {
            let distance = (scale * relative).min(length * 0.125);
            let candidate = midpoint + inward * distance;
            if domain.include(candidate) {
                candidates.push(candidate);
                break;
            }
        }
    }
    candidates
}

fn strict_interior_point(poly: &PolylineCurve<Point2>, _tol: f64) -> Option<Point2> {
    polyline_interior_candidates(poly).into_iter().next()
}

fn split_non_simple_wire<P: Clone, C: Clone>(wire: &Wire<P, C>) -> Vec<Wire<P, C>> {
    let mut pending = vec![wire.iter().cloned().collect::<Vec<_>>()];
    let mut simple = Vec::new();
    while let Some(edges) = pending.pop() {
        if edges.len() < 2 {
            continue;
        }
        let starts = edges
            .iter()
            .map(|edge| edge.front().id())
            .collect::<Vec<_>>();
        let mut seen = HashMap::default();
        let revisit = starts
            .iter()
            .enumerate()
            .find_map(|(j, vertex)| seen.insert(*vertex, j).map(|i| (i, j)));
        let Some((i, j)) = revisit else {
            simple.push(Wire::from(edges));
            continue;
        };
        let sub = edges[i..j].to_vec();
        let outer = edges[..i]
            .iter()
            .chain(edges[j..].iter())
            .cloned()
            .collect::<Vec<_>>();
        if sub.len() >= 2 {
            pending.push(sub);
        }
        if outer.len() >= 2 {
            pending.push(outer);
        }
    }
    simple
}

pub(super) fn sample_points_on_region<C, S>(
    boundaries: &[Wire<Point3, C>],
    surface: &S,
    tol: f64,
) -> Option<Vec<Point3>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3> + Clone,
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let mut polys = HashMap::default();
    let parameter_loops = boundaries
        .iter()
        .flat_map(split_non_simple_wire)
        .filter_map(|wire| {
            create_parameter_boundary_on_surface(surface, &wire, None, &mut polys, tol)
        })
        .filter(|poly| polyline_signed_area(poly).abs() > polyline_area_epsilon(poly))
        .collect::<Vec<_>>();
    if parameter_loops.is_empty() {
        return None;
    }
    let winding_at = |point: Point2| {
        parameter_loops.iter().fold(0isize, |winding, poly| {
            if poly.include(point) {
                winding
                    + if polyline_signed_area(poly) > 0.0 {
                        1
                    } else {
                        -1
                    }
            } else {
                winding
            }
        })
    };
    let mut points = Vec::new();
    for poly in parameter_loops
        .iter()
        .filter(|poly| polyline_signed_area(poly) > 0.0)
    {
        for uv in polyline_interior_candidates(poly) {
            if winding_at(uv) > 0 {
                points.push(surface.evaluate(uv.x, uv.y));
                break;
            }
        }
    }
    (!points.is_empty()).then_some(points)
}

fn same_edge_set<C>(lhs: &BoundaryWire<Point3, C>, rhs: &BoundaryWire<Point3, C>) -> bool {
    lhs.len() == rhs.len()
        && lhs
            .edge_iter()
            .all(|edge| rhs.edge_iter().any(|other| edge.id() == other.id()))
}

fn polyline_signed_area(poly: &PolylineCurve<Point2>) -> f64 {
    let Some(origin) = poly.first().copied() else {
        return 0.0;
    };
    if poly.len() < 3 {
        return 0.0;
    }
    (1..poly.len() - 1)
        .map(|index| {
            let lhs = poly[index] - origin;
            let rhs = poly[index + 1] - origin;
            lhs.x * rhs.y - lhs.y * rhs.x
        })
        .sum::<f64>()
        * 0.5
}

/// Numerical degeneracy threshold in the face's own parameter domain.
///
/// Boolean `tol` is a world-space length. UV coordinates can be arbitrarily
/// scaled or sheared by an equivalent surface parameterization, so feeding that
/// length directly into a UV-area threshold is not meaningful. Instead compare
/// signed area with the loop's own UV bounding-box area, yielding a
/// dimensionless relative test.
fn polyline_area_epsilon(poly: &PolylineCurve<Point2>) -> f64 {
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
    let width = (max_x - min_x).abs();
    let height = (max_y - min_y).abs();
    let bbox_area = width * height;
    if bbox_area.is_finite() {
        TOLERANCE * bbox_area
    } else {
        f64::INFINITY
    }
}

fn divide_one_face<C, S>(
    face: &Face<Point3, C, S>,
    loops: &Loops<Point3, C>,
    sampled_loops: Option<&Loops<Point3, PolylineCurve<Point3>>>,
    tol: f64,
) -> Option<Vec<(FaceFragment<C, S>, ShapesOpStatus)>>
where
    C: BoundedCurve<Point = Point3> + ParameterDivision1D<Point = Point3>,
    S: Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let (mut pre_faces, mut negative_wires) = (Vec::new(), Vec::new());
    let mut map = HashMap::default();
    loops
        .iter()
        .enumerate()
        .try_for_each(|(wire_index, wire)| {
            let sampled_wire = sampled_loops
                .and_then(|sampled| sampled.get(wire_index))
                .map(Deref::deref);
            let poly = create_parameter_boundary(face, wire, sampled_wire, &mut map, tol)?;
            let area = polyline_signed_area(&poly);
            if area.abs() <= polyline_area_epsilon(&poly) {
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
        let chunk_area_abs = polyline_signed_area(&chunk.poly).abs();

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
                let area_tol =
                    polyline_area_epsilon(&outer.poly).max(polyline_area_epsilon(&chunk.poly));
                polyline_signed_area(&outer.poly).abs() > chunk_area_abs + area_tol
                    && outer.poly.include(witness)
            })
            .min_by(|(_, lhs), (_, rhs)| {
                polyline_signed_area(&lhs[0].poly)
                    .abs()
                    .total_cmp(&polyline_signed_area(&rhs[0].poly).abs())
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
            let area_tol =
                polyline_area_epsilon(&outer.poly).max(polyline_area_epsilon(&chunk.poly));
            (polyline_signed_area(&outer.poly).abs() - chunk_area_abs).abs() <= area_tol
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
    Some(
        pre_faces
            .into_iter()
            .filter(|pre_face| !pre_face.is_empty())
            .map(|pre_face| {
                let surface = face.surface();
                let status = pre_face
                    .iter()
                    .find(|chunk| chunk.wire.status() != ShapesOpStatus::Unknown)
                    .map_or(ShapesOpStatus::Unknown, |chunk| chunk.wire.status());
                let wires = pre_face
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
                (
                    FaceFragment {
                        wires,
                        surface,
                        orientation: face.orientation(),
                    },
                    status,
                )
            })
            .collect(),
    )
}

#[cfg(test)]
mod uv_area_tests {
    use super::*;

    #[test]
    fn real_thin_uv_region_is_not_degenerate() {
        let poly = PolylineCurve(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0e-9),
            Point2::new(0.0, 1.0e-9),
        ]);
        let area = polyline_signed_area(&poly).abs();
        assert!(area > polyline_area_epsilon(&poly));
    }

    #[test]
    fn collinear_uv_loop_is_degenerate() {
        let poly = PolylineCurve(vec![
            Point2::new(-2.0, -2.0),
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 2.0),
        ]);
        assert!(polyline_signed_area(&poly).abs() <= polyline_area_epsilon(&poly));
    }

    #[test]
    fn uv_degeneracy_ratio_is_scale_invariant() {
        let unit = PolylineCurve(vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 0.2),
            Point2::new(0.0, 0.2),
        ]);
        let scaled = PolylineCurve(
            unit.iter()
                .map(|point| Point2::new(point.x * 1.0e6, point.y * 1.0e6))
                .collect(),
        );
        let unit_ratio = polyline_signed_area(&unit).abs() / polyline_area_epsilon(&unit);
        let scaled_ratio = polyline_signed_area(&scaled).abs() / polyline_area_epsilon(&scaled);
        assert!((unit_ratio / scaled_ratio - 1.0).abs() < TOLERANCE);
    }
}

pub(super) fn divide_faces<C, S>(
    shell: &Shell<Point3, C, S>,
    loops_store: &LoopsStore<Point3, C>,
    sampled_store: Option<&LoopsStore<Point3, PolylineCurve<Point3>>>,
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
        .enumerate()
        .try_for_each(|(source_face, (face, loops))| {
            let sampled_loops = sampled_store.and_then(|sampled| sampled.get(source_face));
            if loops
                .iter()
                .all(|wire| wire.status() == ShapesOpStatus::Unknown)
            {
                res.push(face.clone(), ShapesOpStatus::Unknown, source_face);
            } else {
                divide_one_face(face, loops, sampled_loops, tol)?
                    .into_iter()
                    .for_each(|(fragment, status)| {
                        res.push_region(
                            fragment.wires,
                            fragment.surface,
                            fragment.orientation,
                            status,
                            source_face,
                        );
                    });
            }
            Some(())
        })?;
    Some(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::{BsplineCurve, KnotVector, Line, Plane};

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
            BsplineCurve::new_unchecked(
                KnotVector::bezier_knot(1),
                vec![front.point(), back.point()],
            ),
        )
    }

    #[test]
    fn synchronized_mesh_samples_use_mesh_projection_budget() {
        let front = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let back = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let exact: Wire<_, _> =
            vec![Edge::new(&front, &back, Line(front.point(), back.point()))].into();
        let sampled_front = Vertex::new(front.point());
        let sampled_back = Vertex::new(back.point());
        let sampled_curve = PolylineCurve(vec![
            front.point(),
            Point3::new(0.5, 0.0, 5.0e-4),
            back.point(),
        ]);
        let sampled: Wire<_, _> =
            vec![Edge::new(&sampled_front, &sampled_back, sampled_curve)].into();
        let mut polys = HashMap::default();

        let boundary = create_parameter_boundary_on_surface(
            &Plane::xy(),
            &exact,
            Some(&sampled),
            &mut polys,
            1.0e-6,
        );
        assert!(
            boundary.is_some(),
            "trusted mesh evidence must project within the mesh budget"
        );
        if let Some(boundary) = boundary {
            assert_eq!(boundary.len(), 3);
            assert!(boundary[1].near(&Point2::new(0.5, 0.0)));
        }
    }

    #[test]
    fn exact_fallback_does_not_inherit_mesh_projection_budget() {
        let plane = Plane::xy();
        let point = Point3::new(0.5, 0.0, 5.0e-4);
        assert!(project_to_parameter(&plane, point, None, 1.0e-6).is_none());
        assert!(
            project_to_parameter(&plane, point, None, BooleanTolerance::new(1.0e-6).mesh_seed)
                .is_some()
        );
    }

    #[test]
    fn non_simple_region_sampling_respects_touching_negative_notch() {
        let a = Vertex::new(Point3::new(0.0, 0.0, 0.0));
        let b = Vertex::new(Point3::new(4.0, 0.0, 0.0));
        let c = Vertex::new(Point3::new(4.0, 4.0, 0.0));
        let d = Vertex::new(Point3::new(0.0, 4.0, 0.0));
        let e = Vertex::new(Point3::new(1.0, 0.0, 0.0));
        let f = Vertex::new(Point3::new(0.0, 1.0, 0.0));
        let boundary: Wire<_, _> = vec![
            test_line(&a, &b),
            test_line(&b, &c),
            test_line(&c, &d),
            test_line(&d, &a),
            test_line(&a, &f),
            test_line(&f, &e),
            test_line(&e, &a),
        ]
        .into();

        let samples =
            sample_points_on_region(&[boundary], &Plane::xy(), 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(samples.iter().all(|point| {
            point.x >= 0.0
                && point.x <= 4.0
                && point.y >= 0.0
                && point.y <= 4.0
                && !(point.x > 0.0 && point.y > 0.0 && point.x + point.y < 1.0)
        }));
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

        let divided = divide_one_face(&face, &loops, None, 0.01).expect("split must succeed");
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

        let divided =
            divide_one_face(&face, &loops, None, 0.01).expect("closed cut must partition");
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

        let divided =
            divide_one_face(&face, &loops, None, 0.01).expect("hole nesting must succeed");
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
