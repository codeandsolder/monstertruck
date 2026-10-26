//! Canonicalize selected Boolean regions into shared manifold topology.

use super::curve_equivalence::{
    CurveCoincidence, validate_endpoint_aligned_coincidence, validate_sampled_coincidence,
};
use super::*;

fn certified_homogeneous_curve_box(curve: &NurbsCurve<Vector4>) -> Option<BoundingBox<Point3>> {
    let mut bounds = BoundingBox::new();
    let mut any = false;
    for control in curve.control_points() {
        let weight = control.w;
        if !(weight > 0.0 && weight.is_finite()) {
            return None;
        }
        let point = Point3::new(control.x / weight, control.y / weight, control.z / weight);
        if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
            return None;
        }
        bounds.push(point);
        any = true;
    }
    any.then_some(bounds)
}

fn certified_alt_curve_box<C, S>(curve: &AltCurve<C, S>) -> Option<BoundingBox<Point3>>
where
    C: BoundedCurve<Point = Point3>,
    S: ShapeOpsSurface, {
    match curve {
        Alternative::FirstType(curve) => curve.certified_spatial_bounds(),
        Alternative::SecondType(curve) => conversion::exact_iso_leader(curve)
            .and_then(|exact| certified_homogeneous_curve_box(&exact)),
    }
}

fn point_inside_padded_box(bounds: BoundingBox<Point3>, point: Point3, pad: f64) -> bool {
    let min = bounds.min();
    let max = bounds.max();
    (0..3).all(|axis| point[axis] >= min[axis] - pad && point[axis] <= max[axis] + pad)
}

pub(super) fn weld_coincident_regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
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
        validate_endpoint_aligned_coincidence::<Point3, _>(&lhs, &rhs)
            == CurveCoincidence::SampledValidation
            || validate_sampled_coincidence::<Point3, _>(&lhs, &rhs)
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
    let source_midpoint = |curve: &AltCurve<C, S>| {
        let (t0, t1) = curve.range_tuple();
        curve.evaluate(f64::midpoint(t0, t1))
    };
    let source_edges = source_edges_by_id
        .into_values()
        .map(|(front, back, curve)| {
            let midpoint = source_midpoint(&curve);
            (front, back, midpoint, curve)
        })
        .collect::<Vec<_>>();

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
                            // A certified box may only reject impossible source evidence.
                            // Every surviving point still goes through the historical inverse-map
                            // and geometric-nearness proof below; `None` preserves the old scan.
                            let curve_bounds = certified_alt_curve_box(&curve);
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
                                let projection_pad = 2.0 * TOLERANCE;
                                let projected_vertices = source_vertex_points
                                    .iter()
                                    .filter(|(_, point)| {
                                        curve_bounds.is_none_or(|bounds| {
                                            point_inside_padded_box(bounds, **point, projection_pad)
                                        })
                                    })
                                    .filter_map(|(id, point)| {
                                        parameter_on_curve(*point)
                                            .map(|parameter| (*id, (*point, parameter)))
                                    })
                                    .collect::<HashMap<_, _>>();
                                for (
                                    candidate_front,
                                    candidate_back,
                                    candidate_midpoint,
                                    _candidate_curve,
                                ) in &source_edges
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
                                    if curve_bounds.is_some_and(|bounds| {
                                        !point_inside_padded_box(
                                            bounds,
                                            *candidate_midpoint,
                                            projection_pad,
                                        )
                                    }) {
                                        continue;
                                    }
                                    if parameter_on_curve(*candidate_midpoint).is_none() {
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
        assert!(
            boundary
                .vertex_iter()
                .all(|vertex| !vertex.point().near(&a))
        );
        for expected in [b, c, d, e, f] {
            assert!(
                boundary
                    .vertex_iter()
                    .any(|vertex| vertex.point().near(&expected)),
                "materialized notch boundary must preserve geometric corner {expected:?}",
            );
        }
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

        assert!(
            whole.vertex_iter().any(|vertex| vertex.point().near(&p1))
                && partitioned
                    .vertex_iter()
                    .any(|vertex| vertex.point().near(&p1)),
            "the shared geometric vertex must partition both coincident boundaries",
        );
        assert!(whole.iter().all(|edge| {
            partitioned
                .iter()
                .any(|candidate| candidate.id() == edge.id())
        }));
        assert!(
            partitioned
                .iter()
                .all(|edge| { whole.iter().any(|candidate| candidate.id() == edge.id()) })
        );
    }
}
