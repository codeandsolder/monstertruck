#[cfg(any(feature = "analyzers", feature = "tessellation"))]
use crate::MetricSpace;
use crate::{EuclideanSpace, Faces, InnerSpace, Point3, Vector3, Vertex, Zero};

// The closest point on a zero-length segment is its sole endpoint. Dividing
// by the squared length would otherwise produce NaN for collapsed mesh edges.
#[cfg(any(feature = "analyzers", feature = "tessellation"))]
pub fn point_segment_distance2(point: Point3, front: Point3, back: Point3) -> f64 {
    let segment = back - front;
    let length2 = segment.dot(segment);
    if length2 == 0.0 {
        point.distance2(front)
    } else {
        let t = ((point - front).dot(segment) / length2).clamp(0.0, 1.0);
        point.distance2(front + segment * t)
    }
}

mod adjacency;
mod face_normal;
pub use adjacency::Adjacency;
pub use face_normal::FaceNormal;
