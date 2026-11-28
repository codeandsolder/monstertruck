#[cfg(any(feature = "analyzers", feature = "tessellation"))]
use crate::MetricSpace;
use crate::{Faces, Point3, Vector3, Vertex, Zero};

// The closest point on a zero-length segment is its sole endpoint. Dividing
// by the squared length would otherwise produce NaN for collapsed mesh edges.
#[cfg(any(feature = "analyzers", feature = "tessellation"))]
pub fn point_segment_distance2(point: Point3, front: Point3, back: Point3) -> f64 {
    let segment = back - front;
    let length2 = segment.dot(&segment);
    if length2 == 0.0 {
        point.distance2(front)
    } else {
        let t = ((point - front).dot(&segment) / length2).clamp(0.0, 1.0);
        point.distance2(front + segment * t)
    }
}

mod adjacency;
mod face_normal;
pub use adjacency::Adjacency;
pub use face_normal::FaceNormal;

/// Truncate coordinates to signed grid indices, rejecting nonfinite or overflowing values.
#[cfg(any(feature = "filters", all(feature = "vtk", not(target_arch = "wasm32"))))]
pub fn checked_grid_coordinates<const N: usize>(coordinates: [f64; N]) -> Option<[i64; N]> {
    use num_traits::ToPrimitive as _;
    let mut result = [0; N];
    for (out, value) in result.iter_mut().zip(coordinates) {
        *out = value.to_i64()?;
    }
    Some(result)
}

#[cfg(all(test, feature = "filters"))]
#[test]
fn grid_cast_rejects_nonfinite_and_overflow() {
    assert_eq!(checked_grid_coordinates([1.9, -1.9]), Some([1, -1]));
    assert_eq!(checked_grid_coordinates([f64::NAN]), None);
    assert_eq!(checked_grid_coordinates([f64::INFINITY]), None);
    assert_eq!(checked_grid_coordinates([2.0_f64.powi(63)]), None);
}
