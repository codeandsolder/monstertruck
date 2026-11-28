use super::{
    AffinePoint, InnerProduct, Origin, Point3, PolygonMesh, Vector3, distance2_point_triangle,
    hash, nonpositive_tolerance,
};
use array_macro::array;

use super::super::{
    active_sweep::ActiveSweep, compare_sweep_coordinate, triangle_projection_interval,
};

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum EndPointType {
    Front,
    Middle,
    Back,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct EndPoint {
    entity: f64,
    r#type: EndPointType,
    index: usize,
}

impl EndPoint {
    #[inline]
    const fn new(entity: f64, r#type: EndPointType, index: usize) -> Self {
        Self {
            entity,
            r#type,
            index,
        }
    }
    #[inline]
    const fn from_seg(seg: (f64, f64), index: usize) -> [Self; 2] {
        [
            Self::new(seg.0, EndPointType::Front, index),
            Self::new(seg.1, EndPointType::Back, index),
        ]
    }
}

fn sorted_endpoints<'a, I, J>(iter0: I, iter1: J, unit: Vector3, tol: f64) -> Vec<EndPoint>
where
    I: IntoIterator<Item = [Point3; 3]>,
    J: IntoIterator<Item = &'a Point3>, {
    let mut res: Vec<EndPoint> = iter0
        .into_iter()
        .enumerate()
        .filter(|(_, tri)| !(tri[1] - tri[0]).cross(&(tri[2] - tri[0])).so_small())
        .flat_map(|(i, tri)| EndPoint::from_seg(triangle_projection_interval(tri, unit, tol), i))
        .chain(iter1.into_iter().enumerate().map(|(i, point)| {
            EndPoint::new(point.position_vector().dot(&unit), EndPointType::Middle, i)
        }))
        .collect();
    res.sort_unstable_by(|a, b| {
        compare_sweep_coordinate(a.entity, b.entity)
            .then_with(|| a.r#type.cmp(&b.r#type))
            .then_with(|| a.index.cmp(&b.index))
    });
    res
}

fn sorted_endpoints_by_polymesh_points(
    polygon: &PolygonMesh,
    points: &[Point3],
    tol: f64,
) -> Vec<EndPoint> {
    let unit = if !polygon.positions().is_empty() {
        hash::take_one_unit(polygon.positions()[0])
    } else if !points.is_empty() {
        hash::take_one_unit(points[0])
    } else {
        return Vec::new();
    };
    sorted_endpoints(
        polygon
            .faces()
            .triangle_iter()
            .map(|tri| array![i => polygon.positions()[tri[i].pos]; 3]),
        points.iter(),
        unit,
        tol,
    )
}

pub fn pointcloud_in_polygon_neighborhood(
    polygon: &PolygonMesh,
    points: &[Point3],
    tol: f64,
) -> bool {
    nonpositive_tolerance!(tol, 0.0);
    let triangles = polygon.faces().triangle_iter().collect::<Vec<_>>();
    let mut current = ActiveSweep::new(triangles.len());
    sorted_endpoints_by_polymesh_points(polygon, points, tol)
        .into_iter()
        .all(move |EndPoint { r#type, index, .. }| match r#type {
            EndPointType::Front => current.insert(index),
            EndPointType::Back => current.remove(index),
            EndPointType::Middle => current.indices().iter().any(|i| {
                let tri = triangles[*i];
                let tri = array![i => polygon.positions()[tri[i].pos]; 3];
                distance2_point_triangle(points[index], tri) < tol * tol
            }),
        })
}
