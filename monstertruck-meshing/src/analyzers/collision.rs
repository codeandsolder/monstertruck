use super::{
    BoundingBox, EuclideanSpace, InnerSpace, MetricSpace, Origin, Point3, PolygonMesh,
    StandardVertex, TOLERANCE2, Vector3, hash,
};
use array_macro::array;

/// Find collisions between two polygon meshes and extract interference lines.
///
/// # Details
///
/// The algorithm for mesh extraction is as follows.
/// All the triangles are projected toward the axis to form intervals.
/// By looking at the overlap of these intervals, we can narrow down the interfering triangles.
/// The overlap of the intervals is obtained by sorting the endpoints.
///
/// # Remarks
///
/// We see that the surfaces that are in contact are not interfering with each other.
/// Therefore, polygon meshes included in one plane are not interfering.
pub trait Collision {
    /// If `self` and `other` collide, then returns only one interference line.
    /// Otherwise, returns `None`.
    fn collide_with(&self, other: &PolygonMesh) -> Option<(Point3, Point3)>;
    /// Extract all interference lines between `self` and `other`.
    /// # Remarks
    /// The results is not arranged so that included lines make continuous maximal polyline curve.
    fn extract_interference(&self, other: &PolygonMesh) -> Vec<(Point3, Point3)>;
}

impl Collision for PolygonMesh {
    #[inline]
    fn collide_with(&self, other: &PolygonMesh) -> Option<(Point3, Point3)> {
        are_colliding(self, other)
    }
    #[inline]
    fn extract_interference(&self, other: &PolygonMesh) -> Vec<(Point3, Point3)> {
        collision(self, other)
    }
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum EndPointType {
    Front,
    Back,
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
struct EndPoint {
    entity: f64,
    r#type: EndPointType,
    segnum: usize,
    index: usize,
}

impl EndPoint {
    #[inline]
    const fn new(entity: f64, r#type: EndPointType, segnum: usize, index: usize) -> Self {
        Self {
            entity,
            r#type,
            segnum,
            index,
        }
    }
    #[inline]
    fn from_seg(seg: (f64, f64), segnum: usize, index: usize) -> Vec<Self> {
        vec![
            Self::new(seg.0, EndPointType::Front, segnum, index),
            Self::new(seg.1, EndPointType::Back, segnum, index),
        ]
    }
}

fn tri_to_seg(tri: [Point3; 3], unit: Vector3) -> (f64, f64) {
    let a = tri[0].to_vec().dot(unit);
    let b = tri[1].to_vec().dot(unit);
    let c = tri[2].to_vec().dot(unit);
    (f64::min(f64::min(a, b), c), f64::max(f64::max(a, b), c))
}

fn sorted_endpoints<I, J>(first_triangles: I, second_triangles: J, unit: Vector3) -> Vec<EndPoint>
where
    I: IntoIterator<Item = [Point3; 3]>,
    J: IntoIterator<Item = [Point3; 3]>, {
    let mut res: Vec<EndPoint> = first_triangles
        .into_iter()
        .enumerate()
        .filter(|(_, tri)| !(tri[1] - tri[0]).cross(tri[2] - tri[0]).so_small())
        .flat_map(|(i, tri)| EndPoint::from_seg(tri_to_seg(tri, unit), 0, i))
        .chain(
            second_triangles
                .into_iter()
                .enumerate()
                .filter(|(_, tri)| !(tri[1] - tri[0]).cross(tri[2] - tri[0]).so_small())
                .flat_map(|(i, tri)| EndPoint::from_seg(tri_to_seg(tri, unit), 1, i)),
        )
        .collect();
    res.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Greater));
    res
}

fn colliding_segment_pairs(sort_endpoints: Vec<EndPoint>) -> impl Iterator<Item = (usize, usize)> {
    let mut current = [Vec::<usize>::new(), Vec::<usize>::new()];
    sort_endpoints
        .into_iter()
        .filter_map(
            move |EndPoint {
                      r#type,
                      segnum,
                      index,
                      ..
                  }| match r#type {
                EndPointType::Front => {
                    current[segnum].push(index);
                    Some(
                        current[1 - segnum]
                            .clone()
                            .into_iter()
                            .map(move |i| if segnum == 0 { (index, i) } else { (i, index) }),
                    )
                }
                EndPointType::Back => {
                    let i = current[segnum].iter().position(|idx| *idx == index)?;
                    current[segnum].swap_remove(i);
                    None
                }
            },
        )
        .flatten()
}

fn disjoint_bdbs(first_triangle: [Point3; 3], second_triangle: [Point3; 3]) -> bool {
    let bdb0: BoundingBox<Point3> = first_triangle.iter().collect();
    let bdb1: BoundingBox<Point3> = second_triangle.iter().collect();
    bdb0.max()[0] < bdb1.min()[0]
        || bdb1.max()[0] < bdb0.min()[0]
        || bdb0.max()[1] < bdb1.min()[1]
        || bdb1.max()[1] < bdb0.min()[1]
        || bdb0.max()[2] < bdb1.min()[2]
        || bdb1.max()[2] < bdb0.min()[2]
}

fn collide_seg_triangle(seg: [Point3; 2], tri: [Point3; 3]) -> Option<Point3> {
    let ab = tri[1] - tri[0];
    let bc = tri[2] - tri[1];
    let ca = tri[0] - tri[2];
    let nor = ab.cross(ca);
    if nor.so_small() {
        return None;
    }
    let ap = seg[0] - tri[0];
    let aq = seg[1] - tri[0];
    let front_plane_distance = ap.dot(nor);
    let back_plane_distance = aq.dot(nor);
    if front_plane_distance * back_plane_distance > 0.0 {
        return None;
    }
    let h = seg[0]
        + front_plane_distance / (front_plane_distance - back_plane_distance) * (seg[1] - seg[0]);
    if f64::signum(ab.cross(nor).dot(h - tri[0]) + TOLERANCE2)
        + f64::signum(bc.cross(nor).dot(h - tri[1]) + TOLERANCE2)
        + f64::signum(ca.cross(nor).dot(h - tri[2]) + TOLERANCE2)
        >= 2.0
    {
        Some(h)
    } else {
        None
    }
}

fn collide_triangles(
    first_triangle: [Point3; 3],
    second_triangle: [Point3; 3],
) -> Option<(Point3, Point3)> {
    let mut tuple = (None, None);
    for point in [
        collide_seg_triangle([first_triangle[0], first_triangle[1]], second_triangle),
        collide_seg_triangle([first_triangle[1], first_triangle[2]], second_triangle),
        collide_seg_triangle([first_triangle[2], first_triangle[0]], second_triangle),
        collide_seg_triangle([second_triangle[0], second_triangle[1]], first_triangle),
        collide_seg_triangle([second_triangle[1], second_triangle[2]], first_triangle),
        collide_seg_triangle([second_triangle[2], second_triangle[0]], first_triangle),
    ] {
        match tuple {
            (None, _) => tuple.0 = point,
            (Some(_), None) => tuple.1 = point,
            (Some(ref mut first), Some(ref mut second)) => {
                if let Some(point) = point {
                    let first_distance = point.distance2(*first);
                    let second_distance = point.distance2(*second);
                    let endpoint_distance = first.distance2(*second);
                    if endpoint_distance < first_distance {
                        *second = point;
                    } else if endpoint_distance < second_distance {
                        *first = point;
                    }
                }
            }
        }
    }
    match tuple {
        (Some(first), Some(second)) => Some((first, second)),
        _ => None,
    }
}

const fn make_pos_tri(poly: &PolygonMesh, face: [StandardVertex; 3]) -> [Point3; 3] {
    array![i => poly.positions()[face[i].pos]; 3]
}

fn collision(poly0: &PolygonMesh, poly1: &PolygonMesh) -> Vec<(Point3, Point3)> {
    let unit = if poly0.positions().is_empty() {
        return Vec::new();
    } else {
        hash::take_one_unit(poly0.positions()[0])
    };
    let first_faces = poly0.faces().triangle_iter().collect::<Vec<_>>();
    let second_faces = poly1.faces().triangle_iter().collect::<Vec<_>>();
    let first_triangles = first_faces.iter().map(|face| make_pos_tri(poly0, *face));
    let second_triangles = second_faces.iter().map(|face| make_pos_tri(poly1, *face));
    colliding_segment_pairs(sorted_endpoints(first_triangles, second_triangles, unit))
        .filter_map(|(idx0, idx1)| {
            let first_triangle = make_pos_tri(poly0, first_faces[idx0]);
            let second_triangle = make_pos_tri(poly1, second_faces[idx1]);
            if disjoint_bdbs(first_triangle, second_triangle) {
                None
            } else {
                collide_triangles(first_triangle, second_triangle)
            }
        })
        .collect()
}

fn are_colliding(poly0: &PolygonMesh, poly1: &PolygonMesh) -> Option<(Point3, Point3)> {
    let unit = if poly0.positions().is_empty() {
        return None;
    } else {
        hash::take_one_unit(poly0.positions()[0])
    };
    let first_faces = poly0.faces().triangle_iter().collect::<Vec<_>>();
    let second_faces = poly1.faces().triangle_iter().collect::<Vec<_>>();
    let first_triangles = first_faces.iter().map(|face| make_pos_tri(poly0, *face));
    let second_triangles = second_faces.iter().map(|face| make_pos_tri(poly1, *face));
    colliding_segment_pairs(sorted_endpoints(first_triangles, second_triangles, unit)).find_map(
        |(idx0, idx1)| {
            let first_triangle = make_pos_tri(poly0, first_faces[idx0]);
            let second_triangle = make_pos_tri(poly1, second_faces[idx1]);
            if disjoint_bdbs(first_triangle, second_triangle) {
                None
            } else {
                collide_triangles(first_triangle, second_triangle)
            }
        },
    )
}

#[test]
fn collide_triangles_test() {
    let first_triangle = [
        Point3::origin(),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let second_triangle = [
        Point3::new(0.0, 0.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
    ];
    assert!(collide_triangles(first_triangle, second_triangle).is_some());
    let first_triangle = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let second_triangle = [
        Point3::new(0.0, 0.0, 0.5),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
    ];
    assert!(collide_triangles(first_triangle, second_triangle).is_none());
}
