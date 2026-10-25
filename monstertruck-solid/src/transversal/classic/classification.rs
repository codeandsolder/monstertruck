//! Geometric inside/outside and coincident-face classification.

use super::*;

/// Classify a tessellated shell point using the signed winding convention shared
/// with [`IncludingPointInDomain`], with several deterministic ray directions to
/// avoid a single degenerate ray deciding the result.
pub(super) fn classify_inside_with_polyshell(
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
                // Tessellation stores each polygon in the underlying surface's
                // orientation and represents topological reversal on the Face.
                // Apply that orientation to the signed crossing contribution;
                // cloning and inverting the whole polygon per ray is unnecessary.
                let polygon = face.surface()?;
                let crossing = polygon.signed_crossing_faces(point, direction);
                Some(
                    count
                        + if face.orientation() {
                            crossing
                        } else {
                            -crossing
                        },
                )
            })?;
            if count >= 1 {
                Some((inside + 1, outside))
            } else {
                Some((inside, outside + 1))
            }
        })?;
    Some(inside_votes > outside_votes)
}

pub(super) fn sample_points_on_face<C, S>(
    face: &Face<Point3, C, S>,
    tol: f64,
) -> Option<Vec<Point3>>
where
    C: PolylineableCurve + Clone,
    S: MeshableSurface + Clone,
{
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

pub(super) fn classify_unknown_face<C, S>(
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
pub(super) enum FaceRelation {
    Inside,
    Outside,
    CoincidentSameFacing,
}

pub(super) fn classify_face_relation<C, S>(
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

pub(super) fn planar_partition_same_facing<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
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
