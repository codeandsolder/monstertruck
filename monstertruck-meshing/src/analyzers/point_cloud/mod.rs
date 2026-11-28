use crate::common::point_segment_distance2;
use array_macro::array;

use super::{
    AffinePoint, BoundingBox, InnerProduct, KERNEL_EPSILON, MetricSpace, Origin, Point3,
    PolygonMesh, Vector3, hash, nonpositive_tolerance,
};
mod hashed_point_cloud;
use hashed_point_cloud::HashedPointCloud;
mod projection_index;
mod sort_end_points;
pub use projection_index::ProjectionNeighborhood;

/// Investigates positional relation between polygon mesh and point cloud.
pub trait WithPointCloud {
    /// Whether all faces of the polygon mesh `self` has intersection with the neighborhood of `point_cloud`.
    /// # Arguments
    /// - `tol`: the radius of the neighborhoods of points in point cloud.
    /// # Panics
    /// `tol` must be greater than or equal to `KERNEL_EPSILON`.
    /// # Examples
    /// ```
    /// use monstertruck_meshing::prelude::*;
    /// let mesh = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.0, 0.0),
    ///             Point3::new(1.0, 0.0, 0.0),
    ///             Point3::new(0.0, 1.0, 0.0),
    ///             Point3::new(0.0, 0.0, 2.0),
    ///             Point3::new(1.0, 0.0, 2.0),
    ///             Point3::new(0.0, 1.0, 2.0),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(vec![[0, 1, 2], [3, 4, 5]]),
    /// ).expect("valid example mesh");
    ///
    /// let mut point_cloud = vec![Point3::new(0.25, 0.25, 0.0)];
    /// assert!(!mesh.is_clung_to_by(&point_cloud, 0.001));
    /// point_cloud.push(Point3::new(0.25, 0.25, 2.0));
    /// assert!(mesh.is_clung_to_by(&point_cloud, 0.001));
    /// ```
    fn is_clung_to_by(&self, point_cloud: &[Point3], tol: f64) -> bool;
    /// Whether the neighborhood of the polygon mesh `self` includes `point_cloud`.
    /// # Panics
    /// `tol` must be greater than or equal to `0.0`.
    fn neighborhood_include(&self, point_cloud: &[Point3], tol: f64) -> bool;
    /// Whether the polygon mesh `self` and `point_cloud` collides.
    /// # Panics
    /// `tol` must be greater than or equal to `0.0`.
    fn collide_with_neighborhood_of(&self, point_cloud: &[Point3], tol: f64) -> bool;
}

impl WithPointCloud for PolygonMesh {
    #[inline]
    fn is_clung_to_by(&self, point_cloud: &[Point3], tol: f64) -> bool {
        nonpositive_tolerance!(tol);
        let Some(point_cloud) = HashedPointCloud::from_points(point_cloud, tol * 2.0) else {
            return false;
        };
        // Covering every face is equivalent to the maximum of their nearest
        // distances being below tolerance, but permits early rejection.
        let mut triangles = self.faces().triangle_iter().peekable();
        triangles.peek().is_some()
            && triangles.all(|tri| {
                let triangle = array![i => self.positions()[tri[i].pos]; 3];
                point_cloud.distance2(&triangle) < tol * tol
            })
    }
    #[inline]
    fn collide_with_neighborhood_of(&self, point_cloud: &[Point3], tol: f64) -> bool {
        nonpositive_tolerance!(tol, 0.0);
        if tol == 0.0 {
            return false;
        }
        let Some(point_cloud) = HashedPointCloud::from_points(point_cloud, tol * 2.0) else {
            return false;
        };
        point_cloud.is_colliding(&self, tol)
    }
    #[inline]
    fn neighborhood_include(&self, point_cloud: &[Point3], tol: f64) -> bool {
        sort_end_points::pointcloud_in_polygon_neighborhood(self, point_cloud, tol)
    }
}

// https://iquilezles.org/www/articles/distfunctions/distfunctions.htm
fn distance2_point_triangle(point: Point3, triangle: [Point3; 3]) -> f64 {
    let ab = triangle[1] - triangle[0];
    let ap = point - triangle[0];
    let bc = triangle[2] - triangle[1];
    let bp = point - triangle[1];
    let ca = triangle[0] - triangle[2];
    let cp = point - triangle[2];
    let nor = ab.cross(&ca);

    let coef = f64::signum(ab.cross(&nor).dot(&ap))
        + f64::signum(bc.cross(&nor).dot(&bp))
        + f64::signum(ca.cross(&nor).dot(&cp));
    if coef < 2.0 || nor.norm().so_small() {
        let a = point_segment_distance2(point, triangle[0], triangle[1]);
        let b = point_segment_distance2(point, triangle[1], triangle[2]);
        let c = point_segment_distance2(point, triangle[2], triangle[0]);
        f64::min(f64::min(a, b), c)
    } else {
        nor.dot(&ap) * nor.dot(&ap) / nor.norm_squared()
    }
}

#[test]
fn degenerate_triangles_have_finite_point_and_segment_distances() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let point = Point3::new(3.0, 4.0, 5.0);
    assert_eq!(
        distance2_point_triangle(point, [origin, origin, origin]),
        50.0
    );

    // A collapsed triangular face is a segment, irrespective of which
    // vertex was repeated. Its distance can be checked analytically.
    let end = Point3::new(2.0, 0.0, 0.0);
    let near = Point3::new(1.0, 1.0, 0.0);
    for triangle in [
        [origin, origin, end],
        [origin, end, end],
        [origin, end, Point3::new(1.0, 0.0, 0.0)],
    ] {
        assert_eq!(distance2_point_triangle(near, triangle), 1.0);
    }
    assert_eq!(
        distance2_point_triangle(Point3::new(3.0, 4.0, 0.0), [origin, origin, end]),
        17.0
    );
}
