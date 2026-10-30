//! Planar-region evidence and coplanar interface imprinting.

use super::*;

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

pub(super) fn planar_polygons_same_trimmed_face(
    lhs: &PolygonMesh,
    rhs: &PolygonMesh,
    tol: f64,
) -> bool {
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub(in crate::transversal::classic) struct PlanarPartitionMatch {
    pub(in crate::transversal::classic) representative: usize,
    pub(in crate::transversal::classic) members: Vec<usize>,
}

pub(in crate::transversal::classic) fn direct_planar_face_match(
    target: &PolygonMesh,
    fragments: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<Option<usize>> {
    for (index, face) in fragments.iter().enumerate() {
        let polygon = face.surface()?;
        if planar_polygons_same_trimmed_face(target, &polygon, tol) {
            return Some(Some(index));
        }
    }
    Some(None)
}

/// Match one trimmed planar face against a partition of the same physical
/// interface on the other shell. Revolving a radial cap, for example, may
/// represent one disk as two semicircular faces while a sweep uses one disk.
pub(in crate::transversal::classic) fn single_face_matches_planar_partition(
    target: &PolygonMesh,
    fragments: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<Option<PlanarPartitionMatch>> {
    let Some((origin, normal)) = planar_polygon_plane(target, tol) else {
        return Some(None);
    };
    let Some(target_signature) = planar_polygon_signature(target, origin, normal) else {
        return Some(None);
    };
    let mut matched = Vec::<PolygonMesh>::new();
    let mut matched_indices = Vec::<usize>::new();
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
        matched_indices.push(index);
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
            .then_some(PlanarPartitionMatch {
                representative,
                members: matched_indices,
            }),
    )
}

#[derive(Clone, Debug)]
pub(in crate::transversal::classic) struct FullPlanarInterfaceMatch {
    pub(in crate::transversal::classic) lhs_faces: Vec<usize>,
    pub(in crate::transversal::classic) rhs_faces: Vec<usize>,
}

pub(super) fn shells_full_planar_interface_pair(
    lhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    rhs: &Shell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tol: f64,
) -> Option<Option<FullPlanarInterfaceMatch>> {
    for (left_index, left) in lhs.iter().enumerate() {
        let left_polygon = left.surface()?;
        if let Some(right) = single_face_matches_planar_partition(&left_polygon, rhs, tol)? {
            return Some(Some(FullPlanarInterfaceMatch {
                lhs_faces: vec![left_index],
                rhs_faces: right.members,
            }));
        }
    }
    for (right_index, right) in rhs.iter().enumerate() {
        let right_polygon = right.surface()?;
        if let Some(left) = single_face_matches_planar_partition(&right_polygon, lhs, tol)? {
            return Some(Some(FullPlanarInterfaceMatch {
                lhs_faces: left.members,
                rhs_faces: vec![right_index],
            }));
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

pub(super) fn shells_have_positive_aabb_overlap(
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

pub(super) fn shells_have_partial_planar_overlap(
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

pub(super) fn imprint_edges_on_faces<C>(
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
                curve.evaluate((t0 + t1) * 0.5),
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

            if target[face_index].reuse_coincident_edge(edge) {
                continue;
            }
            let mut emap = HashMap::default();
            target.add_imprint_vertex(face_index, edge.front(), &mut emap)?;
            target.add_imprint_vertex(face_index, edge.back(), &mut emap)?;
            if target[face_index].has_partitioned_coincident_path(edge) {
                continue;
            }
            target[face_index].add_edge(edge.clone(), ShapesOpStatus::And, true);
        }
    }
    Some(())
}

pub(super) const CROSS_FACE_STITCH_TOLERANCE: f64 = 8.0 * TOLERANCE;

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
        .expect("fixture mesh is valid")
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
        .expect("fixture mesh is valid")
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
        )
        .expect("fixture mesh is valid");
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
