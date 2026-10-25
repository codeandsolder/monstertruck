//! Self-contained classic Boolean pipeline, descended from the published 0.3.2 backend.
//!
//! The historical backend remains the algorithmic base, but this subtree now has
//! explicit staged analysis/classification, typed failure boundaries, canonical
//! topology repair, and stricter geometric evidence semantics. It is compiled only
//! when the marching SSI is the active backend
//! (default features); a build with an external SSI backend never touches it
//! and runs the current `integrate` pipeline unchanged. The public entry points in
//! `transversal::integrate` cfg-dispatch `and`/`or`/`difference`/
//! `symmetric_difference` here.
//!
//! Graft-boundary adaptations vs. the 0.3.2 source (all documented inline):
//! - The Alternative intersection-curve arm is a raw
//!   `IntersectionCurve<PolylineCurve<Point3>, S, S>` (0.3.2 wrapped it with
//!   parameter polylines the pipeline never consumed).
//! - `altshell_to_shell` converts that arm back to `C` via
//!   `SurfaceCurve::with_boundaries(..).into()` (the conversion the current
//!   `ShapeOpsCurve` bound provides), not `IntersectionCurve::new(..).into()`.
//! - Internal stages return typed `ShapeOpsError` variants; only an assembled
//!   but topologically invalid candidate may select a documented fallback path.

mod curve_equivalence;
mod divide_face;
mod faces_classification;
mod intersection_curve;
mod loops_store;
mod tolerance;

use super::integrate::{
    InvalidOutputShellConditionData, ShapeOpsCurve, ShapeOpsError, ShapeOpsSurface,
};
use crate::alternative::Alternative;
use curve_equivalence::{
    CurveCoincidence, validate_endpoint_aligned_coincidence, validate_sampled_coincidence,
};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use monstertruck_topology::*;
use std::collections::{HashMap, VecDeque};
use tolerance::BooleanTolerance;

type ClassicResult<T> = std::result::Result<T, ShapeOpsError>;

type AltCurve<C, S> = Alternative<C, IntersectionCurve<PolylineCurve<Point3>, S, S>>;
type AltCurveShell<C, S> = Shell<Point3, AltCurve<C, S>, S>;
type AltRegion<C, S> = faces_classification::BooleanRegion<Point3, AltCurve<C, S>, S>;

fn weld_coincident_regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
) -> Option<AltCurveShell<C, S>> {
    type AltEdge<C, S> = Edge<Point3, AltCurve<C, S>>;

    type VertexCell = [i64; 3];
    type VertexKey = VertexId<Point3>;

    let mut vertex_cells = HashMap::<VertexCell, Vec<Vertex<Point3>>>::new();
    let mut edges = HashMap::<(VertexKey, VertexKey), Vec<AltEdge<C, S>>>::new();
    let vertex_cell = |point: Point3| -> VertexCell {
        [
            (point.x / TOLERANCE).floor() as i64,
            (point.y / TOLERANCE).floor() as i64,
            (point.z / TOLERANCE).floor() as i64,
        ]
    };
    let canonical_vertex = |point: Point3, cells: &mut HashMap<VertexCell, Vec<Vertex<Point3>>>| {
        let key = vertex_cell(point);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor = [
                        key[0].saturating_add(dx),
                        key[1].saturating_add(dy),
                        key[2].saturating_add(dz),
                    ];
                    if let Some(vertex) = cells
                        .get(&neighbor)
                        .and_then(|bucket| bucket.iter().find(|vertex| vertex.point().near(&point)))
                    {
                        return vertex.clone();
                    }
                }
            }
        }
        let vertex = Vertex::new(point);
        cells.entry(key).or_default().push(vertex.clone());
        vertex
    };
    let curves_sampled_coincident = |lhs: &AltEdge<C, S>, rhs: &AltEdge<C, S>| {
        let lhs = lhs.oriented_curve();
        let rhs = rhs.oriented_curve();
        validate_sampled_coincidence::<Point3, _>(&lhs, &rhs) == CurveCoincidence::SampledValidation
            || validate_endpoint_aligned_coincidence::<Point3, _>(&lhs, &rhs)
                == CurveCoincidence::SampledValidation
    };

    type SourceEdgeId<C, S> = EdgeId<AltCurve<C, S>>;
    type SourceVertexId = VertexId<Point3>;
    let mut source_edges_by_id =
        HashMap::<SourceEdgeId<C, S>, (SourceVertexId, SourceVertexId, AltCurve<C, S>)>::new();
    let mut source_vertex_points = HashMap::<SourceVertexId, Point3>::new();
    for edge in regions
        .iter()
        .flat_map(|region| region.absolute_boundaries())
        .flat_map(|wire| wire.iter())
    {
        let absolute = edge.absolute_clone();
        let front = absolute.absolute_front();
        let back = absolute.absolute_back();
        source_vertex_points
            .entry(front.id())
            .or_insert(front.point());
        source_vertex_points
            .entry(back.id())
            .or_insert(back.point());
        source_edges_by_id
            .entry(absolute.id())
            .or_insert_with(|| (front.id(), back.id(), absolute.oriented_curve()));
    }
    let source_edges = source_edges_by_id.into_values().collect::<Vec<_>>();

    for point in source_vertex_points.values().copied() {
        canonical_vertex(point, &mut vertex_cells);
    }

    let canonical_edge = |rebuilt: AltEdge<C, S>,
                          edges: &mut HashMap<(VertexKey, VertexKey), Vec<AltEdge<C, S>>>|
     -> AltEdge<C, S> {
        let forward = (rebuilt.front().id(), rebuilt.back().id());
        let reverse = (rebuilt.back().id(), rebuilt.front().id());
        if let Some(canonical) = edges
            .get(&forward)
            .into_iter()
            .chain(edges.get(&reverse))
            .flatten()
            .find(|candidate| curves_sampled_coincident(candidate, &rebuilt))
        {
            if canonical.front() == rebuilt.front() {
                canonical.clone()
            } else {
                canonical.inverse()
            }
        } else {
            edges.entry(forward).or_default().push(rebuilt.clone());
            rebuilt
        }
    };

    let rebuilt_regions = regions
        .iter()
        .map(|region| {
            let wires = region
                .absolute_boundaries()
                .iter()
                .map(|wire| {
                    let rebuilt = wire
                        .iter()
                        .flat_map(|edge| {
                            let front = canonical_vertex(edge.front().point(), &mut vertex_cells);
                            let back = canonical_vertex(edge.back().point(), &mut vertex_cells);
                            if front == back {
                                return Vec::new();
                            }

                            let curve = edge.oriented_curve();
                            let (t0, t1) = curve.range_tuple();
                            let span = t1 - t0;
                            let parameter_on_curve = |point: Point3| {
                                curve
                                    .search_parameter(point, None, 100)
                                    .or_else(|| curve.search_nearest_parameter(point, None, 100))
                                    .filter(|&parameter| curve.evaluate(parameter).near(&point))
                            };
                            // A split point is only useful when an existing source edge proves
                            // that the reference edge is already partitioned there. Scan those
                            // source edges directly instead of testing every canonical vertex and,
                            // for each matching vertex, scanning every source edge again. Besides
                            // being equivalent evidence, this removes an O(E * V * E) inverse-map
                            // search from large Boolean results.
                            let mut split_vertices = Vec::new();
                            if span.abs() > f64::EPSILON {
                                let projected_vertices = source_vertex_points
                                    .iter()
                                    .filter_map(|(id, point)| {
                                        parameter_on_curve(*point)
                                            .map(|parameter| (*id, (*point, parameter)))
                                    })
                                    .collect::<HashMap<_, _>>();
                                for (candidate_front, candidate_back, candidate_curve) in
                                    &source_edges
                                {
                                    let Some(&(front_point, front_parameter)) =
                                        projected_vertices.get(candidate_front)
                                    else {
                                        continue;
                                    };
                                    let Some(&(back_point, back_parameter)) =
                                        projected_vertices.get(candidate_back)
                                    else {
                                        continue;
                                    };
                                    let (candidate_t0, candidate_t1) =
                                        candidate_curve.range_tuple();
                                    let midpoint = candidate_curve
                                        .evaluate((candidate_t0 + candidate_t1) * 0.5);
                                    if parameter_on_curve(midpoint).is_none() {
                                        continue;
                                    }

                                    for (point, parameter) in [
                                        (front_point, front_parameter),
                                        (back_point, back_parameter),
                                    ] {
                                        let progress = (parameter - t0) / span;
                                        if !(progress > 0.0 && progress < 1.0) {
                                            continue;
                                        }
                                        let vertex = canonical_vertex(point, &mut vertex_cells);
                                        if vertex != front
                                            && vertex != back
                                            && !split_vertices
                                                .iter()
                                                .any(|(_, _, existing)| existing == &vertex)
                                        {
                                            split_vertices.push((progress, parameter, vertex));
                                        }
                                    }
                                }
                            }
                            split_vertices.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));
                            split_vertices.dedup_by(|lhs, rhs| lhs.2 == rhs.2);

                            let mut pieces = Vec::with_capacity(split_vertices.len() + 1);
                            let mut remaining = curve;
                            let mut segment_front = front;
                            for (_, parameter, split_vertex) in split_vertices {
                                let later = remaining.cut(parameter);
                                remaining
                                    .snap_endpoints(segment_front.point(), split_vertex.point());
                                let segment = Edge::new(&segment_front, &split_vertex, remaining);
                                pieces.push(canonical_edge(segment, &mut edges));
                                remaining = later;
                                segment_front = split_vertex;
                            }
                            remaining.snap_endpoints(segment_front.point(), back.point());
                            let segment = Edge::new(&segment_front, &back, remaining);
                            pieces.push(canonical_edge(segment, &mut edges));
                            pieces
                        })
                        .collect::<Vec<_>>();
                    Wire::from(rebuilt)
                })
                .collect::<Vec<_>>();
            faces_classification::BooleanRegion::new(
                wires,
                region.surface(),
                region.orientation(),
                region.key(),
                region.source_face(),
            )
        })
        .collect::<Vec<_>>();
    faces_classification::materialize_regions(&rebuilt_regions)
}

#[cfg(test)]
fn weld_coincident_topology<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &AltCurveShell<C, S>,
) -> AltCurveShell<C, S> {
    let regions = shell
        .iter()
        .enumerate()
        .map(|(key, face)| {
            faces_classification::BooleanRegion::new(
                face.absolute_boundaries().clone(),
                face.surface(),
                face.orientation(),
                key,
                key,
            )
        })
        .collect::<Vec<_>>();
    weld_coincident_regions(&regions).unwrap_or_default()
}

#[cfg(test)]
mod weld_tests {
    use super::*;
    use monstertruck_modeling::{Curve, Surface};

    fn alt_line(front: Point3, back: Point3) -> AltCurve<Curve, Surface> {
        Alternative::FirstType(Curve::from(Line(front, back)))
    }

    #[test]
    fn weld_drops_face_collapsed_to_exact_edge_retrace() {
        let p0 = Point3::origin();
        let p1 = Point3::new(TOLERANCE * 0.5, 0.0, 0.0);
        let p2 = Point3::new(1.0, 1.0, 0.0);
        let v0 = Vertex::new(p0);
        let v1 = Vertex::new(p1);
        let v2 = Vertex::new(p2);
        let wire: Wire<Point3, AltCurve<Curve, Surface>> = vec![
            Edge::new(&v0, &v1, alt_line(p0, p1)),
            Edge::new(&v1, &v2, alt_line(p1, p2)),
            Edge::new(&v2, &v0, alt_line(p2, p0)),
        ]
        .into();
        let face = Face::new_unchecked(
            vec![wire],
            Surface::Plane(Plane::new(
                Point3::origin(),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            )),
        );
        let shell: AltCurveShell<Curve, Surface> = Shell::from(vec![face]);

        let welded = weld_coincident_topology(&shell);

        assert!(
            welded.is_empty(),
            "a face reduced to one edge traversed both ways has zero area"
        );
    }

    #[test]
    fn weld_materializes_touching_negative_notch_as_simple_boundary() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(4.0, 0.0, 0.0);
        let c = Point3::new(4.0, 4.0, 0.0);
        let d = Point3::new(0.0, 4.0, 0.0);
        let e = Point3::new(1.0, 0.0, 0.0);
        let f = Point3::new(0.0, 1.0, 0.0);
        let va = Vertex::new(a);
        let vb = Vertex::new(b);
        let vc = Vertex::new(c);
        let vd = Vertex::new(d);
        let ve = Vertex::new(e);
        let vf = Vertex::new(f);
        let wire: Wire<Point3, AltCurve<Curve, Surface>> = vec![
            Edge::new(&va, &vb, alt_line(a, b)),
            Edge::new(&vb, &vc, alt_line(b, c)),
            Edge::new(&vc, &vd, alt_line(c, d)),
            Edge::new(&vd, &va, alt_line(d, a)),
            Edge::new(&va, &vf, alt_line(a, f)),
            Edge::new(&vf, &ve, alt_line(f, e)),
            Edge::new(&ve, &va, alt_line(e, a)),
        ]
        .into();
        let surface = Surface::Plane(Plane::new(
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ));
        let shell: AltCurveShell<Curve, Surface> =
            Shell::from(vec![Face::new_unchecked(vec![wire], surface)]);

        let welded = weld_coincident_topology(&shell);

        assert_eq!(welded.len(), 1);
        assert_eq!(welded[0].absolute_boundaries().len(), 1);
        let boundary = &welded[0].absolute_boundaries()[0];
        assert!(boundary.is_simple());
        assert_eq!(boundary.len(), 5);
        assert!(
            boundary
                .vertex_iter()
                .all(|vertex| !vertex.point().near(&a))
        );
    }

    #[test]
    fn weld_partitions_coincident_edges_at_shared_vertices() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(2.0, 0.0, 0.0);
        let p3 = Point3::new(0.0, 1.0, 0.0);
        let v0 = Vertex::new(p0);
        let v1 = Vertex::new(p1);
        let v2 = Vertex::new(p2);
        let v3 = Vertex::new(p3);

        let whole: Wire<Point3, AltCurve<Curve, Surface>> = vec![
            Edge::new(&v0, &v2, alt_line(p0, p2)),
            Edge::new(&v2, &v3, alt_line(p2, p3)),
            Edge::new(&v3, &v0, alt_line(p3, p0)),
        ]
        .into();
        let partitioned: Wire<Point3, AltCurve<Curve, Surface>> = vec![
            Edge::new(&v0, &v1, alt_line(p0, p1)),
            Edge::new(&v1, &v2, alt_line(p1, p2)),
            Edge::new(&v2, &v3, alt_line(p2, p3)),
            Edge::new(&v3, &v0, alt_line(p3, p0)),
        ]
        .into();
        let surface = Surface::Plane(Plane::new(
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ));
        let shell: AltCurveShell<Curve, Surface> = Shell::from(vec![
            Face::new_unchecked(vec![whole], surface.clone()),
            Face::new_unchecked(vec![partitioned], surface),
        ]);

        let welded = weld_coincident_topology(&shell);
        let whole = &welded[0].absolute_boundaries()[0];
        let partitioned = &welded[1].absolute_boundaries()[0];

        assert_eq!(whole.len(), 4);
        assert_eq!(partitioned.len(), 4);
        assert!(whole.iter().all(|edge| {
            partitioned
                .iter()
                .any(|candidate| candidate.id() == edge.id())
        }));
    }
}

/// Convert a polyline leader to an exactly equivalent degree-1 B-spline.
///
/// `PolylineCurve` uses integer segment parameters, so scale the uniform knot
/// vector to the same `[0, segment_count]` domain. This preserves both the
/// guide geometry and parameterization exactly.
fn polyline_leader_bspline(polyline: &PolylineCurve<Point3>) -> Option<BsplineCurve<Point3>> {
    let segment_count = polyline.len().checked_sub(1)?;
    if segment_count == 0 {
        return None;
    }
    let mut knots = KnotVector::uniform_knot(1, segment_count);
    knots.transform(segment_count as f64, 0.0);
    Some(BsplineCurve::new(knots, polyline.iter().copied().collect()))
}

#[derive(Clone, Copy)]
enum IsoParameter {
    U,
    V,
}

fn trim_homogeneous_curve(
    mut curve: NurbsCurve<Vector4>,
    start: f64,
    end: f64,
) -> Option<NurbsCurve<Vector4>> {
    if !start.is_finite() || !end.is_finite() || (start - end).abs() <= TOLERANCE {
        return None;
    }
    let reverse = start > end;
    let (lo, hi) = if reverse { (end, start) } else { (start, end) };
    let (range0, range1) = curve.range_tuple();
    if lo < range0 - TOLERANCE || hi > range1 + TOLERANCE {
        return None;
    }
    if lo > range0 + TOLERANCE {
        curve = curve.cut(lo);
    }
    let (_, current_end) = curve.range_tuple();
    if hi < current_end - TOLERANCE {
        let _tail = curve.cut(hi);
    }
    if reverse {
        curve.invert();
    }
    Some(curve)
}

fn homogeneous_iso_curve<S: ShapeOpsSurface>(
    surface: &S,
    iso: IsoParameter,
    constant: f64,
    start: f64,
    end: f64,
) -> Option<NurbsCurve<Vector4>> {
    let mut homogeneous = surface.try_into_homogeneous_bspline_surface()?;
    let (u_range, v_range) = homogeneous.range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range, v_range);
    let curve = match iso {
        IsoParameter::U => {
            if constant < u0 - TOLERANCE || constant > u1 + TOLERANCE {
                return None;
            }
            let row = if (constant - u0).abs() <= TOLERANCE {
                0
            } else if (constant - u1).abs() <= TOLERANCE {
                homogeneous.control_points().len().checked_sub(1)?
            } else {
                let _right = homogeneous.cut_u(constant);
                homogeneous.control_points().len().checked_sub(1)?
            };
            NurbsCurve::new(homogeneous.curve_v(row))
        }
        IsoParameter::V => {
            if constant < v0 - TOLERANCE || constant > v1 + TOLERANCE {
                return None;
            }
            let column = if (constant - v0).abs() <= TOLERANCE {
                0
            } else if (constant - v1).abs() <= TOLERANCE {
                homogeneous.control_points().first()?.len().checked_sub(1)?
            } else {
                let _right = homogeneous.cut_v(constant);
                homogeneous.control_points().first()?.len().checked_sub(1)?
            };
            NurbsCurve::new(homogeneous.curve_u(column))
        }
    };
    trim_homogeneous_curve(curve, start, end)
}

fn point_on_surface<S: ShapeOpsSurface>(surface: &S, point: Point3) -> bool {
    surface
        .search_parameter(point, None, 100)
        .or_else(|| surface.search_nearest_parameter(point, None, 100))
        .is_some_and(|(u, v)| surface.evaluate(u, v).distance(point) <= 4.0 * TOLERANCE)
}

fn prove_iso_candidate<S: ShapeOpsSurface>(
    candidate: &NurbsCurve<Vector4>,
    surface0: &S,
    surface1: &S,
    leader: &PolylineCurve<Point3>,
) -> bool {
    let (t0, t1) = candidate.range_tuple();
    if candidate.subs(t0).distance(leader.front()) > 8.0 * TOLERANCE
        || candidate.subs(t1).distance(leader.back()) > 8.0 * TOLERANCE
    {
        return false;
    }
    (0..=32).all(|index| {
        let ratio = index as f64 / 32.0;
        let point = candidate.subs(t0 * (1.0 - ratio) + t1 * ratio);
        point_on_surface(surface0, point) && point_on_surface(surface1, point)
    })
}

fn leader_matches_iso_locus<S: ShapeOpsSurface>(
    surface: &S,
    params: &[(f64, f64)],
    leader: &PolylineCurve<Point3>,
    iso: IsoParameter,
    constant: f64,
) -> bool {
    params.len() == leader.len()
        && constant.is_finite()
        && params
            .iter()
            .copied()
            .zip(leader.iter().copied())
            .all(|((u, v), point)| {
                let (iso_u, iso_v) = match iso {
                    IsoParameter::U => (constant, v),
                    IsoParameter::V => (u, constant),
                };
                surface.evaluate(iso_u, iso_v).distance(point) <= 4.0 * TOLERANCE
            })
}

fn exact_iso_leader<S: ShapeOpsSurface>(
    intersection: &IntersectionCurve<PolylineCurve<Point3>, S, S>,
) -> Option<NurbsCurve<Vector4>> {
    let leader = intersection.leader();
    if leader.len() < 2 {
        return None;
    }
    for surface in [intersection.surface0(), intersection.surface1()] {
        let Some(params) = leader
            .iter()
            .copied()
            .map(|point| {
                surface
                    .search_parameter(point, None, 100)
                    .or_else(|| surface.search_nearest_parameter(point, None, 100))
            })
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };
        for iso in [IsoParameter::U, IsoParameter::V] {
            let values = params
                .iter()
                .map(|&(u, v)| match iso {
                    IsoParameter::U => u,
                    IsoParameter::V => v,
                })
                .collect::<Vec<_>>();
            let constant = values.iter().sum::<f64>() / values.len() as f64;
            // Inverse parameter search on tightly trimmed B-spline surfaces can
            // wander by hundreds of ulps even for an exact iso-parameter
            // intersection. Prove constancy in 3D instead: every marching
            // sample must lie on the exact iso locus within the Boolean validation
            // tolerance. The candidate is independently checked against both
            // carrier surfaces below before it can replace the marching leader.
            if !leader_matches_iso_locus(surface, &params, leader, iso, constant) {
                continue;
            }
            let start = match iso {
                IsoParameter::U => params.first()?.1,
                IsoParameter::V => params.first()?.0,
            };
            let end = match iso {
                IsoParameter::U => params.last()?.1,
                IsoParameter::V => params.last()?.0,
            };
            let Some(candidate) = homogeneous_iso_curve(surface, iso, constant, start, end) else {
                continue;
            };
            if prove_iso_candidate(
                &candidate,
                intersection.surface0(),
                intersection.surface1(),
                leader,
            ) {
                return Some(candidate);
            }
        }
    }
    None
}

/// Convert an Alternative-curve shell back to the target curve type.
///
/// Surface intersection curves already carry a cleaned marching leader. When
/// that leader is provably an iso-parameter curve of either exact carrier,
/// replace the mesh-derived polyline by the exact homogeneous section before
/// serialization. Otherwise keep the leader exactly as a degree-1 B-spline.
/// `SurfaceCurve` still carries both supporting surfaces in either case.
fn altshell_to_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    altshell: &AltCurveShell<C, S>,
) -> Option<Shell<Point3, C, S>> {
    altshell.try_mapped(
        |p| Some(*p),
        |c| match c {
            Alternative::FirstType(c) => Some(c.clone()),
            Alternative::SecondType(ic) => {
                let surface0 = ic.surface0().clone();
                let surface1 = ic.surface1().clone();
                let boundary0: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                let boundary1: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                if let Some(nurbs) = exact_iso_leader(ic) {
                    Some(
                        SurfaceCurve::with_boundaries(
                            surface0, surface1, nurbs, boundary0, boundary1,
                        )
                        .into(),
                    )
                } else {
                    let bspline = polyline_leader_bspline(ic.leader())?;
                    Some(
                        SurfaceCurve::with_boundaries(
                            surface0, surface1, bspline, boundary0, boundary1,
                        )
                        .into(),
                    )
                }
            }
        },
        |s| Some(s.clone()),
    )
}

/// Classify a tessellated shell point using the signed winding convention shared
/// with [`IncludingPointInDomain`], with several deterministic ray directions to
/// avoid a single degenerate ray deciding the result.
fn classify_inside_with_polyshell(
    poly_shell: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    point: Point3,
) -> Option<bool> {
    let offsets = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.613, -0.271, 0.149),
        Vector3::new(-0.347, 0.509, -0.221),
        Vector3::new(0.193, 0.401, 0.577),
        Vector3::new(-0.433, -0.127, 0.389),
    ];
    let (inside_votes, outside_votes) = offsets
        .into_iter()
        .map(|offset| hash::take_one_unit(point + offset))
        .try_fold((0usize, 0usize), |(inside, outside), direction| {
            let count = poly_shell.iter().try_fold(0isize, |count, face| {
                let polygon = face.surface()?;
                Some(count + polygon.signed_crossing_faces(point, direction))
            })?;
            if count >= 1 {
                Some((inside + 1, outside))
            } else {
                Some((inside, outside + 1))
            }
        })?;
    Some(inside_votes > outside_votes)
}

fn sample_points_on_face<C, S>(face: &Face<Point3, C, S>, tol: f64) -> Option<Vec<Point3>>
where
    C: PolylineableCurve + Clone,
    S: MeshableSurface + Clone, {
    let shell = Shell::from(vec![face.clone()]);
    let triangulated = shell.triangulation(BooleanTolerance::new(tol).mesh_seed);
    let polygon = triangulated.first()?.surface()?;
    let positions = polygon.positions();
    let mut candidates = polygon
        .faces()
        .triangle_iter()
        .filter_map(|triangle| {
            let a = positions[triangle[0].pos];
            let b = positions[triangle[1].pos];
            let c = positions[triangle[2].pos];
            let area2 = (b - a).cross(c - a).magnitude();
            (area2.is_finite() && area2 > TOLERANCE * TOLERANCE).then_some((
                area2,
                Point3::from_vec((a.to_vec() + b.to_vec() + c.to_vec()) / 3.0),
            ))
        })
        .collect::<Vec<_>>();
    if candidates.is_empty() {
        return None;
    }

    const MAX_CLASSIFICATION_SAMPLES: usize = 9;
    candidates.sort_by(|lhs, rhs| rhs.0.total_cmp(&lhs.0));
    let mut samples = vec![candidates.swap_remove(0).1];
    while samples.len() < MAX_CLASSIFICATION_SAMPLES && !candidates.is_empty() {
        let next = candidates
            .iter()
            .enumerate()
            .map(|(index, (_, point))| {
                let distance2 = samples
                    .iter()
                    .map(|sample| sample.distance2(*point))
                    .fold(f64::INFINITY, f64::min);
                (index, distance2)
            })
            .max_by(|lhs, rhs| lhs.1.total_cmp(&rhs.1))
            .map(|(index, _)| index)?;
        samples.push(candidates.swap_remove(next).1);
    }
    Some(samples)
}

fn classify_unknown_face<C, S>(
    poly_shell: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    face: &Face<Point3, C, S>,
    tol: f64,
) -> Option<bool>
where
    C: PolylineableCurve + Clone,
    S: MeshableSurface + Clone,
{
    let points = sample_points_on_face(face, tol)?;
    let (inside, outside) =
        points
            .into_iter()
            .try_fold((0usize, 0usize), |(inside, outside), point| {
                if classify_inside_with_polyshell(poly_shell, point)? {
                    Some((inside + 1, outside))
                } else {
                    Some((inside, outside + 1))
                }
            })?;
    Some(inside >= outside)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FaceRelation {
    Inside,
    Outside,
    CoincidentSameFacing,
}

fn classify_face_relation<C, S>(
    poly_shell: &Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>,
    face: &Face<Point3, C, S>,
    tol: f64,
) -> Option<FaceRelation>
where
    C: PolylineableCurve + Clone,
    S: MeshableSurface
        + Clone
        + Invertible
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let points = sample_points_on_face(face, tol)?;
    let surface = face.oriented_surface();
    let offset = BooleanTolerance::new(tol).classifier_offset;
    let (inside, outside, coincident, total) = points.into_iter().try_fold(
        (0usize, 0usize, 0usize, 0usize),
        |(inside, outside, coincident, total), point| {
            let parameter = surface
                .search_parameter(point, None, 100)
                .or_else(|| surface.search_nearest_parameter(point, None, 100));
            let Some((u, v)) = parameter else {
                let is_inside = classify_inside_with_polyshell(poly_shell, point)?;
                return Some((
                    inside + usize::from(is_inside),
                    outside + usize::from(!is_inside),
                    coincident,
                    total + 1,
                ));
            };
            let normal = surface.normal(u, v);
            let outward_inside =
                classify_inside_with_polyshell(poly_shell, point + normal * offset)?;
            let inward_inside =
                classify_inside_with_polyshell(poly_shell, point - normal * offset)?;
            if outward_inside {
                Some((inside + 1, outside, coincident, total + 1))
            } else if inward_inside {
                Some((inside, outside, coincident + 1, total + 1))
            } else {
                Some((inside, outside + 1, coincident, total + 1))
            }
        },
    )?;
    let relation = if total != 0 && coincident > inside && coincident > outside {
        FaceRelation::CoincidentSameFacing
    } else if inside > outside {
        FaceRelation::Inside
    } else {
        FaceRelation::Outside
    };
    Some(relation)
}

fn planar_partition_same_facing<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    region: &AltRegion<C, S>,
    candidate: &PolygonMesh,
    other_regions: &[AltRegion<C, S>],
    other_poly: &PolyShell,
    tol: f64,
) -> Option<bool> {
    let representative = match loops_store::direct_planar_face_match(candidate, other_poly, tol)? {
        Some(index) => index,
        None => {
            loops_store::single_face_matches_planar_partition(candidate, other_poly, tol)??
                .representative
        }
    };
    let other = other_regions.get(representative)?;
    let point = other
        .absolute_boundaries()
        .first()?
        .vertex_iter()
        .next()?
        .point();
    let mut lhs = region.surface();
    if !region.orientation() {
        lhs.invert();
    }
    let mut rhs = other.surface();
    if !other.orientation() {
        rhs.invert();
    }
    let (lhs_u, lhs_v) = lhs
        .search_parameter(point, None, 100)
        .or_else(|| lhs.search_nearest_parameter(point, None, 100))?;
    let (rhs_u, rhs_v) = rhs
        .search_parameter(point, None, 100)
        .or_else(|| rhs.search_nearest_parameter(point, None, 100))?;
    Some(lhs.normal(lhs_u, lhs_v).dot(rhs.normal(rhs_u, rhs_v)) > 0.0)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairMode {
    AndOr,
    Difference,
}

type PolyShell = Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>>;

#[derive(Clone, Debug)]
struct FullInterfaceAdjacency {
    shell0: Vec<usize>,
    shell1: Vec<usize>,
}

impl FullInterfaceAdjacency {
    fn contains(&self, shell_index: usize, source_face: usize) -> bool {
        match shell_index {
            0 => self.shell0.contains(&source_face),
            1 => self.shell1.contains(&source_face),
            _ => false,
        }
    }
}

#[derive(Clone, Debug)]
struct PairAnalysis<C, S> {
    tolerance: f64,
    poly_shell0: PolyShell,
    poly_shell1: PolyShell,
    cls0: faces_classification::FacesClassification<Point3, AltCurve<C, S>, S>,
    cls1: faces_classification::FacesClassification<Point3, AltCurve<C, S>, S>,
    coplanar_overlap: bool,
    full_interface_adjacency: Option<FullInterfaceAdjacency>,
}

fn analyze_pair<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<PairAnalysis<C, S>> {
    if tol < TOLERANCE {
        return Err(ShapeOpsError::InvalidTolerance);
    }
    let seed_tol = BooleanTolerance::new(tol).mesh_seed;
    let poly_shell0 = shell0.triangulation(seed_tol);
    let poly_shell1 = shell1.triangulation(seed_tol);
    let altshell0: AltCurveShell<C, S> =
        shell0.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let altshell1: AltCurveShell<C, S> =
        shell1.mapped(|x| *x, |c| Alternative::FirstType(c.clone()), Clone::clone);
    let loops_store::InterfaceLoops {
        geom: [loops_store0, loops_store1],
        coplanar_overlap,
        full_interface_adjacency,
    } = loops_store::create_loops_stores(
        &altshell0,
        &poly_shell0,
        &altshell1,
        &poly_shell1,
        tol,
        mode == PairMode::AndOr,
    )
    .ok_or(ShapeOpsError::InterfaceAnalysisFailed)?;
    let full_interface_adjacency = full_interface_adjacency.map(|matched| FullInterfaceAdjacency {
        shell0: matched.lhs_faces,
        shell1: matched.rhs_faces,
    });
    let cls0 = divide_face::divide_faces(&altshell0, &loops_store0, tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 0 })?;
    let cls1 = divide_face::divide_faces(&altshell1, &loops_store1, tol)
        .ok_or(ShapeOpsError::DivideFacesFailed { shell_index: 1 })?;
    Ok(PairAnalysis {
        tolerance: tol,
        poly_shell0,
        poly_shell1,
        cls0,
        cls1,
        coplanar_overlap,
        full_interface_adjacency,
    })
}

fn classify_unknown<C, S>(
    poly_shell: &PolyShell,
    region: &faces_classification::BooleanRegion<Point3, C, S>,
    tol: f64,
    shell_index: usize,
) -> ClassicResult<bool>
where
    C: PolylineableCurve + ParameterDivision1D<Point = Point3> + Clone,
    S: MeshableSurface
        + Clone
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if let Some(face) = region.try_face() {
        return classify_unknown_face(poly_shell, &face, tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index });
    }
    let points =
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    let (inside, outside) = points
        .into_iter()
        .try_fold((0usize, 0usize), |(inside, outside), point| {
            let is_inside = classify_inside_with_polyshell(poly_shell, point)?;
            Some((
                inside + usize::from(is_inside),
                outside + usize::from(!is_inside),
            ))
        })
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    Ok(inside >= outside)
}

fn classify_relation<C, S>(
    poly_shell: &PolyShell,
    region: &faces_classification::BooleanRegion<Point3, C, S>,
    tol: f64,
    shell_index: usize,
) -> ClassicResult<FaceRelation>
where
    C: PolylineableCurve + ParameterDivision1D<Point = Point3> + Clone,
    S: MeshableSurface
        + Clone
        + Invertible
        + ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    if let Some(face) = region.try_face() {
        return classify_face_relation(poly_shell, &face, tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index });
    }
    let points =
        divide_face::sample_points_on_region(region.absolute_boundaries(), &region.surface(), tol)
            .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    let mut surface = region.surface();
    if !region.orientation() {
        surface.invert();
    }
    let offset = BooleanTolerance::new(tol).classifier_offset;
    let (inside, outside, coincident, total) = points
        .into_iter()
        .try_fold(
            (0usize, 0usize, 0usize, 0usize),
            |(inside, outside, coincident, total), point| {
                let parameter = surface
                    .search_parameter(point, None, 100)
                    .or_else(|| surface.search_nearest_parameter(point, None, 100));
                let Some((u, v)) = parameter else {
                    let is_inside = classify_inside_with_polyshell(poly_shell, point)?;
                    return Some((
                        inside + usize::from(is_inside),
                        outside + usize::from(!is_inside),
                        coincident,
                        total + 1,
                    ));
                };
                let normal = surface.normal(u, v);
                let outward_inside =
                    classify_inside_with_polyshell(poly_shell, point + normal * offset)?;
                let inward_inside =
                    classify_inside_with_polyshell(poly_shell, point - normal * offset)?;
                if outward_inside {
                    Some((inside + 1, outside, coincident, total + 1))
                } else if inward_inside {
                    Some((inside, outside, coincident + 1, total + 1))
                } else {
                    Some((inside, outside + 1, coincident, total + 1))
                }
            },
        )
        .ok_or(ShapeOpsError::UnknownClassificationFailed { shell_index })?;
    Ok(
        if total != 0 && coincident > inside && coincident > outside {
            FaceRelation::CoincidentSameFacing
        } else if inside > outside {
            FaceRelation::Inside
        } else {
            FaceRelation::Outside
        },
    )
}

fn convert_alt_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell: &AltCurveShell<C, S>,
    operation: &'static str,
) -> ClassicResult<Shell<Point3, C, S>> {
    altshell_to_shell(shell).ok_or(ShapeOpsError::AltShellConversionFailed { operation })
}

fn weld_and_convert<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
    operation: &'static str,
) -> ClassicResult<Shell<Point3, C, S>> {
    let shell = weld_coincident_regions(regions)
        .ok_or(ShapeOpsError::RegionMaterializationFailed { operation })?;
    convert_alt_shell(&shell, operation)
}

fn shell_is_valid<P, C, S>(shell: &Shell<P, C, S>) -> bool {
    shell.is_empty()
        || shell
            .connected_components()
            .iter()
            .all(|component| component.check_solid_boundary().is_ok())
}

fn process_difference_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();

    let [_inside0, mut difference0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        if !classify_unknown(&analysis.poly_shell1, &region, tol, 0)? {
            difference0.push(region);
        }
    }
    let [mut difference1, _outside1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if classify_unknown(&analysis.poly_shell0, &region, tol, 1)? {
            difference1.push(region);
        }
    }
    for region in &mut difference1 {
        region.invert();
    }
    difference0.append(&mut difference1);
    let mut label_shell = weld_and_convert(&difference0, "difference")?;
    let _ = orient_regular_shell(&mut label_shell);
    if label_shell
        .connected_components()
        .iter()
        .all(|component| component.check_solid_boundary().is_ok())
    {
        return Ok([label_shell, Shell::default()]);
    }

    let [and0, or0, unknown0] = analysis.cls0.and_or_unknown();
    let [and1, or1, unknown1] = analysis.cls1.and_or_unknown();
    let mut fallback = Vec::new();
    for region in and0.into_iter().chain(or0).chain(unknown0) {
        if !classify_unknown(&analysis.poly_shell1, &region, tol, 0)? {
            fallback.push(region);
        }
    }
    for mut region in and1.into_iter().chain(or1).chain(unknown1) {
        if classify_unknown(&analysis.poly_shell0, &region, tol, 1)? {
            region.invert();
            fallback.push(region);
        }
    }
    Ok([weld_and_convert(&fallback, "difference")?, Shell::default()])
}

fn process_non_coplanar_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();
    let [mut and0, mut or0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        if classify_unknown(&analysis.poly_shell1, &region, tol, 0)? {
            and0.push(region);
        } else {
            or0.push(region);
        }
    }
    let [mut and1, mut or1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if classify_unknown(&analysis.poly_shell0, &region, tol, 1)? {
            and1.push(region);
        } else {
            or1.push(region);
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    Ok([
        weld_and_convert(&and0, "and")?,
        weld_and_convert(&or0, "or")?,
    ])
}

fn partition_evidence<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: &[AltRegion<C, S>],
    mesh_tol: f64,
) -> (Vec<AltRegion<C, S>>, PolyShell, HashMap<usize, PolygonMesh>) {
    let mut valid_regions = Vec::new();
    let mut poly_faces = Vec::new();
    let mut by_key = HashMap::new();
    for region in regions {
        let Some(face) = region.try_face() else {
            continue;
        };
        let triangulated = Shell::from(vec![face]).triangulation(mesh_tol);
        let Some(poly_face) = triangulated.first() else {
            continue;
        };
        let Some(polygon) = poly_face.surface() else {
            continue;
        };
        by_key.insert(region.key(), polygon.clone());
        valid_regions.push(region.clone());
        poly_faces.push(poly_face.clone());
    }
    (valid_regions, Shell::from(poly_faces), by_key)
}

fn process_coplanar_analysis<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut analysis: PairAnalysis<C, S>,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let tol = analysis.tolerance;
    analysis.cls0.integrate_by_component();
    analysis.cls1.integrate_by_component();

    let [mut and0, mut or0, unknown0] = analysis.cls0.and_or_unknown();
    for region in unknown0 {
        if analysis
            .full_interface_adjacency
            .as_ref()
            .is_some_and(|interface| interface.contains(0, region.source_face()))
        {
            continue;
        }
        match classify_relation(&analysis.poly_shell1, &region, tol, 0)? {
            FaceRelation::Inside => and0.push(region),
            FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(region),
        }
    }
    let [mut and1, mut or1, unknown1] = analysis.cls1.and_or_unknown();
    for region in unknown1 {
        if analysis
            .full_interface_adjacency
            .as_ref()
            .is_some_and(|interface| interface.contains(1, region.source_face()))
        {
            continue;
        }
        match classify_relation(&analysis.poly_shell0, &region, tol, 1)? {
            FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(region),
            FaceRelation::Outside => or1.push(region),
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    let label_and_shell = weld_and_convert(&and0, "and")?;
    let label_or_shell = weld_and_convert(&or0, "or")?;
    if !label_or_shell.is_empty()
        && shell_is_valid(&label_and_shell)
        && shell_is_valid(&label_or_shell)
    {
        return Ok([label_and_shell, label_or_shell]);
    }

    let all0 = analysis.cls0.all_regions();
    let all1 = analysis.cls1.all_regions();
    let partition_tolerance = BooleanTolerance::new(tol);
    let partition_geom_tol = partition_tolerance.length;
    let partition_mesh_tol = partition_tolerance.mesh_seed;
    let (valid0, poly0, candidate0) = partition_evidence(&all0, partition_mesh_tol);
    let (valid1, poly1, candidate1) = partition_evidence(&all1, partition_mesh_tol);

    let [and0_labeled, or0_labeled, unknown0] = analysis.cls0.and_or_unknown();
    let mut and0 = Vec::new();
    let mut or0 = Vec::new();
    for region in and0_labeled.into_iter().chain(or0_labeled).chain(unknown0) {
        let coincident = candidate0.get(&region.key()).is_some_and(|candidate| {
            planar_partition_same_facing(&region, candidate, &valid1, &poly1, partition_geom_tol)
                == Some(true)
        });
        let relation = if coincident {
            FaceRelation::CoincidentSameFacing
        } else {
            classify_relation(&analysis.poly_shell1, &region, tol, 0)?
        };
        match relation {
            FaceRelation::Inside => and0.push(region),
            FaceRelation::Outside | FaceRelation::CoincidentSameFacing => or0.push(region),
        }
    }

    let [and1_labeled, or1_labeled, unknown1] = analysis.cls1.and_or_unknown();
    let mut and1 = Vec::new();
    let mut or1 = Vec::new();
    for region in and1_labeled.into_iter().chain(or1_labeled).chain(unknown1) {
        let coincident = candidate1.get(&region.key()).is_some_and(|candidate| {
            planar_partition_same_facing(&region, candidate, &valid0, &poly0, partition_geom_tol)
                == Some(true)
        });
        let relation = if coincident {
            FaceRelation::CoincidentSameFacing
        } else {
            classify_relation(&analysis.poly_shell0, &region, tol, 1)?
        };
        match relation {
            FaceRelation::Inside | FaceRelation::CoincidentSameFacing => and1.push(region),
            FaceRelation::Outside => or1.push(region),
        }
    }
    and0.append(&mut and1);
    or0.append(&mut or1);
    Ok([
        weld_and_convert(&and0, "and")?,
        weld_and_convert(&or0, "or")?,
    ])
}

fn process_one_pair_of_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    shell0: &Shell<Point3, C, S>,
    shell1: &Shell<Point3, C, S>,
    tol: f64,
    mode: PairMode,
) -> ClassicResult<[Shell<Point3, C, S>; 2]> {
    let analysis = analyze_pair(shell0, shell1, tol, mode)?;
    match mode {
        PairMode::Difference => process_difference_analysis(analysis),
        PairMode::AndOr if analysis.coplanar_overlap => process_coplanar_analysis(analysis),
        PairMode::AndOr => process_non_coplanar_analysis(analysis),
    }
}

/// Re-orient a regular manifold shell from shared-edge constraints.
fn orient_regular_shell<P, C, S>(shell: &mut Shell<P, C, S>) -> bool {
    if shell.shell_condition() != ShellCondition::Regular {
        return true;
    }

    let mut adjacency = vec![Vec::<(usize, bool)>::new(); shell.len()];
    for left in 0..shell.len() {
        for right in left + 1..shell.len() {
            let mut relation = None;
            for left_edge in shell[left].edge_iter() {
                for right_edge in shell[right].edge_iter() {
                    if left_edge.id() != right_edge.id() {
                        continue;
                    }
                    let must_differ = left_edge.orientation() == right_edge.orientation();
                    match relation {
                        Some(existing) if existing != must_differ => return false,
                        Some(_) => {}
                        None => relation = Some(must_differ),
                    }
                }
            }
            if let Some(must_differ) = relation {
                adjacency[left].push((right, must_differ));
                adjacency[right].push((left, must_differ));
            }
        }
    }

    let mut flips = vec![None; shell.len()];
    for root in 0..shell.len() {
        if flips[root].is_some() {
            continue;
        }
        flips[root] = Some(false);
        let mut queue = VecDeque::from([root]);
        while let Some(face) = queue.pop_front() {
            let Some(flip) = flips[face] else {
                return false;
            };
            for &(neighbor, must_differ) in &adjacency[face] {
                let required = flip ^ must_differ;
                match flips[neighbor] {
                    Some(existing) if existing != required => return false,
                    Some(_) => {}
                    None => {
                        flips[neighbor] = Some(required);
                        queue.push_back(neighbor);
                    }
                }
            }
        }
    }

    for (face, flip) in shell.iter_mut().zip(flips) {
        if flip == Some(true) {
            face.invert();
        }
    }
    shell.shell_condition() != ShellCondition::Regular
}

fn finalize<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    operation: &'static str,
    mut shell: Shell<Point3, C, S>,
) -> ClassicResult<Solid<Point3, C, S>> {
    let _ = orient_regular_shell(&mut shell);
    let boundaries = shell.connected_components();
    for (shell_index, boundary) in boundaries.iter().enumerate() {
        if boundary.check_solid_boundary().is_ok() {
            continue;
        }
        let open_boundaries = boundary.extract_boundaries();
        let first_boundary = open_boundaries.first();
        let singular = boundary.singular_vertices();
        let mut edge_uses = Vec::<(EdgeId<C>, Edge<Point3, C>, usize)>::new();
        for edge in boundary.edge_iter() {
            if let Some((_, _, uses)) = edge_uses.iter_mut().find(|(id, _, _)| *id == edge.id()) {
                *uses += 1;
            } else {
                edge_uses.push((edge.id(), edge.clone(), 1));
            }
        }
        let overused = edge_uses
            .iter()
            .filter(|(_, _, uses)| *uses > 2)
            .collect::<Vec<_>>();
        let first_overused = overused.first();
        return Err(ShapeOpsError::InvalidOutputShellCondition(Box::new(
            InvalidOutputShellConditionData {
                operation,
                shell_index,
                face_count: boundary.len(),
                empty: boundary.is_empty(),
                connected: boundary.is_connected(),
                condition: boundary.shell_condition(),
                boundary_loops: open_boundaries.len(),
                first_boundary_len: first_boundary.map(|wire| wire.len()),
                first_boundary_front: first_boundary
                    .and_then(Wire::front_vertex)
                    .map(Vertex::point),
                first_boundary_back: first_boundary
                    .and_then(Wire::back_vertex)
                    .map(Vertex::point),
                overused_edges: overused.len(),
                first_overused_uses: first_overused.map(|(_, _, uses)| *uses),
                first_overused_front: first_overused.map(|(_, edge, _)| edge.front().point()),
                first_overused_back: first_overused.map(|(_, edge, _)| edge.back().point()),
                singular_vertices: singular.len(),
                first_singular: singular.first().map(Vertex::point),
            },
        )));
    }
    Solid::try_new(boundaries)
        .map_err(|source| ShapeOpsError::InvalidOutputShell { operation, source })
}

/// AND operation between two solids (classic backend).
fn and_shells<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "and";
    let pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
    };
    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next().unwrap();
    let shell1 = iter1.next().unwrap();
    let [mut and_shell, _] = pair(shell0, shell1)?;
    for shell in iter0 {
        let [res, _] = pair(&and_shell, shell)?;
        and_shell = res;
    }
    for shell in iter1 {
        let [res, _] = pair(&and_shell, shell)?;
        and_shell = res;
    }
    finalize(operation, and_shell)
}

/// AND operation between two solids (classic backend).
pub(crate) fn and<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    and_shells(solid0, solid1, tol)
}

/// OR operation between two solids (classic backend).
pub(crate) fn or<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "or";
    let pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
    };
    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next().unwrap();
    let shell1 = iter1.next().unwrap();
    let [_, mut or_shell] = pair(shell0, shell1)?;
    for shell in iter0 {
        let [_, res] = pair(&or_shell, shell)?;
        or_shell = res;
    }
    for shell in iter1 {
        let [_, res] = pair(&or_shell, shell)?;
        or_shell = res;
    }
    finalize(operation, or_shell)
}

/// Difference: the region inside `solid0` but outside `solid1` (classic backend).
pub(crate) fn difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let operation = "difference";
    let difference_pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::Difference)
    };
    let and_pair = |a: &Shell<Point3, C, S>, b: &Shell<Point3, C, S>| {
        process_one_pair_of_shells(a, b, tol, PairMode::AndOr)
    };

    let mut iter0 = solid0.boundaries().iter();
    let mut iter1 = solid1.boundaries().iter();
    let shell0 = iter0.next().unwrap();
    let shell1 = iter1.next().unwrap();
    let [mut difference_shell, _] = difference_pair(shell0, shell1)?;

    for shell in iter0 {
        let [res, _] = and_pair(&difference_shell, shell)?;
        difference_shell = res;
    }
    for shell in iter1 {
        let [res, _] = difference_pair(&difference_shell, shell)?;
        difference_shell = res;
    }
    finalize(operation, difference_shell)
}

/// Symmetric difference (XOR): the region inside exactly one solid (classic backend).
pub(crate) fn symmetric_difference<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    solid0: &Solid<Point3, C, S>,
    solid1: &Solid<Point3, C, S>,
    tol: f64,
) -> ClassicResult<Solid<Point3, C, S>> {
    let d0 = difference(solid0, solid1, tol)?;
    let d1 = difference(solid1, solid0, tol)?;
    or(&d0, &d1, tol)
}

#[cfg(test)]
mod structural_classification_tests {
    use super::*;
    use monstertruck_modeling::Curve;

    fn tetrahedron_mesh(reverse: bool) -> PolygonMesh {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let mut faces = vec![[0usize, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        if reverse {
            faces.iter_mut().for_each(|face| face.swap(1, 2));
        }
        PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )
    }

    #[test]
    fn signed_winding_does_not_turn_negative_shell_into_inside() {
        let reversed = tetrahedron_mesh(true);
        let shell: Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>> =
            Shell::from(vec![Face::new_unchecked(Vec::new(), Some(reversed))]);
        assert_eq!(
            classify_inside_with_polyshell(&shell, Point3::new(0.1, 0.1, 0.1)),
            Some(false)
        );
    }

    #[test]
    fn face_classification_samples_stay_inside_concave_trim() {
        let points = [
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (2.0, 3.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
        ];
        let vertices = points
            .iter()
            .map(|&(x, y)| Vertex::new(Point3::new(x, y, 0.0)))
            .collect::<Vec<_>>();
        let edges = (0..vertices.len())
            .map(|index| {
                let next = (index + 1) % vertices.len();
                Edge::new(
                    &vertices[index],
                    &vertices[next],
                    Curve::from(Line(vertices[index].point(), vertices[next].point())),
                )
            })
            .collect::<Vec<_>>();
        let face = Face::new_unchecked(
            vec![Wire::from(edges)],
            monstertruck_modeling::Surface::Plane(Plane::xy()),
        );
        let boundary = PolylineCurve(
            points
                .iter()
                .map(|&(x, y)| Point2::new(x, y))
                .chain(std::iter::once(Point2::new(points[0].0, points[0].1)))
                .collect(),
        );
        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(samples.len() <= 9);
        assert!(
            samples
                .iter()
                .all(|point| boundary.include(Point2::new(point.x, point.y)))
        );
    }

    #[test]
    fn face_classification_sampling_keeps_regions_thinner_than_mesh_seed() {
        let points = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0e-4), (0.0, 1.0e-4)];
        let vertices = points
            .iter()
            .map(|&(x, y)| Vertex::new(Point3::new(x, y, 0.0)))
            .collect::<Vec<_>>();
        let edges = (0..vertices.len())
            .map(|index| {
                let next = (index + 1) % vertices.len();
                Edge::new(
                    &vertices[index],
                    &vertices[next],
                    Curve::from(Line(vertices[index].point(), vertices[next].point())),
                )
            })
            .collect::<Vec<_>>();
        let face = Face::new_unchecked(
            vec![Wire::from(edges)],
            monstertruck_modeling::Surface::Plane(Plane::xy()),
        );

        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(
            samples
                .iter()
                .all(|point| point.x > 0.0 && point.x < 1.0 && point.y > 0.0 && point.y < 1.0e-4)
        );
    }

    fn scaled_tetrahedron_mesh(origin: Point3, scale: f64, reverse: bool) -> PolygonMesh {
        let positions = vec![
            origin,
            origin + Vector3::new(scale, 0.0, 0.0),
            origin + Vector3::new(0.0, scale, 0.0),
            origin + Vector3::new(0.0, 0.0, scale),
        ];
        let mut faces = vec![[0usize, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        if reverse {
            faces.iter_mut().for_each(|face| face.swap(1, 2));
        }
        PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )
    }

    #[test]
    fn signed_winding_respects_cavity_orientation() {
        let outer = scaled_tetrahedron_mesh(Point3::origin(), 4.0, false);
        let cavity = scaled_tetrahedron_mesh(Point3::new(0.5, 0.5, 0.5), 1.0, true);
        let shell: Shell<Point3, PolylineCurve<Point3>, Option<PolygonMesh>> = Shell::from(vec![
            Face::new_unchecked(Vec::new(), Some(outer)),
            Face::new_unchecked(Vec::new(), Some(cavity)),
        ]);

        assert_eq!(
            classify_inside_with_polyshell(&shell, Point3::new(0.25, 0.25, 0.25)),
            Some(true)
        );
        assert_eq!(
            classify_inside_with_polyshell(&shell, Point3::new(0.6, 0.6, 0.6)),
            Some(false)
        );
    }

    #[test]
    fn face_classification_samples_avoid_holes() {
        fn wire(points: &[(f64, f64)]) -> Wire<Point3, Curve> {
            let vertices = points
                .iter()
                .map(|&(x, y)| Vertex::new(Point3::new(x, y, 0.0)))
                .collect::<Vec<_>>();
            (0..vertices.len())
                .map(|index| {
                    let next = (index + 1) % vertices.len();
                    Edge::new(
                        &vertices[index],
                        &vertices[next],
                        Curve::from(Line(vertices[index].point(), vertices[next].point())),
                    )
                })
                .collect::<Vec<_>>()
                .into()
        }

        let outer = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
        let hole = [(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)];
        let face = Face::new_unchecked(
            vec![wire(&outer), wire(&hole)],
            monstertruck_modeling::Surface::Plane(Plane::xy()),
        );

        let samples = sample_points_on_face(&face, 1.0e-6).unwrap_or_default();
        assert!(!samples.is_empty());
        assert!(samples.iter().all(|point| {
            let x = point.x;
            let y = point.y;
            x <= 1.0 + TOLERANCE
                || x >= 3.0 - TOLERANCE
                || y <= 1.0 + TOLERANCE
                || y >= 3.0 - TOLERANCE
        }));
    }
}
