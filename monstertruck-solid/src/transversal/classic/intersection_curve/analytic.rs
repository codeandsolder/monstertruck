//! Direct intersection of *certified planar* carriers.
//!
//! A transverse pair has one exact supporting line. Intersect that line with
//! each triangulated trim separately, merge its occupied intervals, and take
//! their intersection. Curved boundaries remain tessellated seeds; no mesh
//! sample is promoted to a proof of curved-edge incidence.
use monstertruck_core::{
    cgmath64::{EuclideanSpace, InnerSpace, Point3, Vector3},
    tolerance::DEFAULT_NUMERICAL_TOLERANCES,
};
use monstertruck_geometry::prelude::{
    AnalyticSurfaceKind, BsplineSurface, Plane, TryIntoAnalyticSurfaceKind, Vector4,
};
use monstertruck_meshing::prelude::{PolygonMesh, PolylineCurve};

const COORD_EPS: f64 = DEFAULT_NUMERICAL_TOLERANCES.structural;
const MIN_TRANSVERSAL_SINE: f64 = DEFAULT_NUMERICAL_TOLERANCES.angular;

type Interval = (f64, f64);

pub(in crate::transversal::classic) fn certified_plane<S: TryIntoAnalyticSurfaceKind>(
    surface: &S,
    net: Option<&BsplineSurface<Vector4>>,
) -> Option<Plane> {
    let AnalyticSurfaceKind::Plane(plane) = surface.try_into_analytic_surface_kind()? else {
        return None;
    };
    // The analytic kind recognizer may tolerate departures from a plane up to
    // modeling/kernel epsilon. This fast path requires a much stricter carrier
    // identity: every homogeneous control lies on the plane, which certifies
    // the full rational B-spline support (assuming positive weights).
    let net = net?;
    let normal = plane.normal();
    let origin = plane.origin();
    net.control_points()
        .iter()
        .flatten()
        .all(|point| {
            point.w.is_finite()
                && point.w > 0.0
                && (normal.dot(Vector3::new(point.x, point.y, point.z) / point.w - origin.to_vec()))
                    .abs()
                    <= COORD_EPS * origin.to_vec().magnitude().max(1.0)
        })
        .then_some(plane)
}

fn merge_intervals(mut intervals: Vec<Interval>, eps: f64) -> Option<Vec<Interval>> {
    if intervals
        .iter()
        .any(|(a, b)| !a.is_finite() || !b.is_finite() || a > b)
    {
        return None;
    }
    intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut merged: Vec<Interval> = Vec::new();
    for (low, high) in intervals {
        if let Some(last) = merged.last_mut()
            && low <= last.1 + eps
        {
            last.1 = last.1.max(high);
            continue;
        }
        merged.push((low, high));
    }
    Some(merged)
}

fn trimmed_line_intervals(
    polygon: &PolygonMesh,
    cutting_plane: &Plane,
    line_origin: Point3,
    line_direction: Vector3,
    eps: f64,
) -> Option<Vec<Interval>> {
    let origin = cutting_plane.origin();
    let normal = cutting_plane.normal();
    let mut intervals = Vec::new();
    let positions = polygon.positions();
    for triangle in polygon.faces().triangle_iter() {
        let points = triangle.map(|vertex| positions[vertex.pos]);
        let distances = points.map(|point| (point - origin).dot(normal));
        if distances.iter().any(|distance| !distance.is_finite()) {
            return None;
        }
        // A triangle entirely on the other plane is a degeneracy (or a
        // numerically indistinguishable near-parallel pair): fall back.
        if distances.iter().all(|distance| distance.abs() <= eps) {
            return None;
        }
        let mut parameters = Vec::with_capacity(3);
        for index in 0..3 {
            let next = (index + 1) % 3;
            let a = distances[index];
            let b = distances[next];
            if a.abs() <= eps {
                parameters.push((points[index] - line_origin).dot(line_direction));
            }
            if (a < -eps && b > eps) || (a > eps && b < -eps) {
                let ratio = a / (a - b);
                let crossing = points[index] + (points[next] - points[index]) * ratio;
                parameters.push((crossing - line_origin).dot(line_direction));
            }
        }
        if parameters.len() < 2 {
            continue; // vertex-only touch; not a one-dimensional trim
        }
        parameters.sort_by(f64::total_cmp);
        let low = parameters[0];
        let high = *parameters.last()?;
        if high - low > eps {
            intervals.push((low, high));
        }
    }
    merge_intervals(intervals, eps)
}

/// An analytical branch returns `Some` including an empty set (disjoint
/// supports/trims). `None` means the pair requires generic numerical SSI.
pub(super) fn planar_intersection(
    plane0: Plane,
    polygon0: &PolygonMesh,
    plane1: Plane,
    polygon1: &PolygonMesh,
) -> Option<Vec<PolylineCurve<Point3>>> {
    let normal0 = plane0.normal();
    let normal1 = plane1.normal();
    let axis = normal0.cross(normal1);
    let axis2 = axis.magnitude2();
    if !axis2.is_finite() {
        return None;
    }
    if axis2 <= MIN_TRANSVERSAL_SINE * MIN_TRANSVERSAL_SINE {
        // A nearly-parallel pair can still intersect very far from the
        // supplied origins. The coplanar classifier and generic SSI own the
        // ambiguous region; never turn an angular tolerance into disjointness.
        return None;
    }
    let direction = axis / axis2.sqrt();
    let delta = (plane1.origin() - plane0.origin()).dot(normal1);
    // Lies on both planes: normal1 . (axis x normal0) = |axis|^2.
    let line_origin = plane0.origin() + axis.cross(normal0) * (delta / axis2);
    let scale = line_origin.to_vec().magnitude().max(1.0);
    let eps = COORD_EPS * scale;
    let a = trimmed_line_intervals(polygon0, &plane1, line_origin, direction, eps)?;
    let b = trimmed_line_intervals(polygon1, &plane0, line_origin, direction, eps)?;
    let mut result = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let lo = a[i].0.max(b[j].0);
        let hi = a[i].1.min(b[j].1);
        if hi - lo > eps {
            result.push(PolylineCurve(vec![
                line_origin + direction * lo,
                line_origin + direction * hi,
            ]));
        }
        if a[i].1 < b[j].1 {
            i += 1;
        } else {
            j += 1;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::TryIntoHomogeneousBsplineSurface;
    use monstertruck_meshing::prelude::{Faces, StandardAttributes};

    fn square_xy(x0: f64, x1: f64, y0: f64, y1: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(x0, y0, 0.0),
                    Point3::new(x1, y0, 0.0),
                    Point3::new(x1, y1, 0.0),
                    Point3::new(x0, y1, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 2, 3]]),
        )
        .expect("valid square")
    }

    fn square_xz(x0: f64, x1: f64, z0: f64, z1: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(x0, 0.0, z0),
                    Point3::new(x1, 0.0, z0),
                    Point3::new(x1, 0.0, z1),
                    Point3::new(x0, 0.0, z1),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 2, 3]]),
        )
        .expect("valid square")
    }

    fn split_xy(gap: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(-2.0, -1.0, 0.0),
                    Point3::new(-gap, -1.0, 0.0),
                    Point3::new(-gap, 1.0, 0.0),
                    Point3::new(-2.0, 1.0, 0.0),
                    Point3::new(gap, -1.0, 0.0),
                    Point3::new(2.0, -1.0, 0.0),
                    Point3::new(2.0, 1.0, 0.0),
                    Point3::new(gap, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]]),
        )
        .expect("valid disconnected trims")
    }

    fn segments_x(result: &[PolylineCurve<Point3>]) -> Vec<(f64, f64)> {
        let mut segments = result
            .iter()
            .map(|line| {
                assert_eq!(line.len(), 2);
                for point in line.iter() {
                    assert!(point.y.abs() < 1.0e-10);
                    assert!(point.z.abs() < 1.0e-10);
                }
                let a = line[0].x;
                let b = line[1].x;
                (a.min(b), a.max(b))
            })
            .collect::<Vec<_>>();
        segments.sort_by(|a, b| a.0.total_cmp(&b.0));
        segments
    }

    #[test]
    fn only_certainly_planar_homogeneous_nets_are_admitted() {
        use monstertruck_geometry::prelude::{BsplineSurface, KnotVector, NurbsSurface, Vector4};
        let plane = Plane::xy();
        let plane_net = plane.try_into_homogeneous_bspline_surface();
        assert!(certified_plane(&plane, plane_net.as_ref()).is_some());
        let mut almost_plane = BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(1), KnotVector::bezier_knot(1)),
            vec![
                vec![
                    Vector4::new(0.0, 0.0, 0.0, 1.0),
                    Vector4::new(0.0, 1.0, 0.0, 1.0),
                ],
                vec![
                    Vector4::new(1.0, 0.0, 0.0, 1.0),
                    Vector4::new(1.0, 1.0, 0.0, 1.0),
                ],
            ],
        );
        let accepted = NurbsSurface::new(almost_plane.clone());
        assert!(certified_plane(&accepted, Some(&almost_plane)).is_some());
        almost_plane.control_points_mut().last().unwrap().z = 1.0e-7;
        let perturbed = NurbsSurface::new(almost_plane.clone());
        assert!(certified_plane(&perturbed, Some(&almost_plane)).is_none());
    }

    #[test]
    fn transverse_planar_trims_intersect_on_analytic_support_line() {
        let xy = Plane::xy();
        let zx = Plane::zx();
        let a = square_xy(-2.0, 3.0, -1.0, 1.0);
        let b = square_xz(-1.0, 1.0, -1.0, 1.0);
        let forward = planar_intersection(xy, &a, zx, &b).expect("analytic plane pair");
        assert_eq!(segments_x(&forward), [(-1.0, 1.0)]);
        let reverse = planar_intersection(zx, &b, xy, &a).expect("swapped planes");
        assert_eq!(segments_x(&reverse), segments_x(&forward));
    }

    #[test]
    fn analytic_line_preserves_disconnected_trim_gaps_and_touches() {
        let xy = Plane::xy();
        let zx = Plane::zx();
        let a = split_xy(0.01);
        let b = square_xz(-3.0, 3.0, -1.0, 1.0);
        let pieces = planar_intersection(xy, &a, zx, &b).expect("analytic plane pair");
        assert_eq!(segments_x(&pieces), [(-2.0, -0.01), (0.01, 2.0)]);

        let far = square_xz(3.0, 4.0, -1.0, 1.0);
        assert!(planar_intersection(xy, &a, zx, &far).unwrap().is_empty());
        let touching = square_xz(2.0, 3.0, -1.0, 1.0);
        assert!(
            planar_intersection(xy, &a, zx, &touching)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn parallel_or_ill_conditioned_planes_use_existing_classifier() {
        let xy = Plane::xy();
        let translated = Plane::new(
            Point3::new(0.0, 0.0, 0.001),
            Point3::new(1.0, 0.0, 0.001),
            Point3::new(0.0, 1.0, 0.001),
        );
        let mesh = square_xy(-1.0, 1.0, -1.0, 1.0);
        assert!(planar_intersection(xy, &mesh, xy, &mesh).is_none());
        assert!(planar_intersection(xy, &mesh, translated, &mesh).is_none());

        let tilted = Plane::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 1.0e-8),
        );
        assert!(planar_intersection(xy, &mesh, tilted, &mesh).is_none());
    }
}
