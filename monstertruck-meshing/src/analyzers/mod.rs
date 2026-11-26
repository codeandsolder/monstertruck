use crate::{
    Adjacency, BoundingBox, EuclideanSpace, FaceNormal, Faces, InnerSpace, KERNEL_EPSILON,
    KERNEL_EPSILON_SQUARED, Matrix3, MetricSpace, Origin, Point3, PolygonMesh, PolylineCurve,
    SquareMatrix, StandardAttributes, StandardVertex, Vector3, Vector4, Vertex, hash,
    nonpositive_tolerance,
};
use std::cmp::Ordering;

// Preserve signed-zero equality while making NaN comparisons deterministic.
const fn sweep_coordinate_key(value: f64) -> u64 {
    const SIGN: u64 = 1 << 63;
    let bits = if value == 0.0 { 0 } else { value.to_bits() };
    if bits & SIGN == 0 { bits ^ SIGN } else { !bits }
}

fn compare_sweep_coordinate(left: f64, right: f64) -> Ordering {
    left.partial_cmp(&right)
        .unwrap_or_else(|| left.total_cmp(&right))
}

// Shared triangle sweep bounds keep collision and point-cloud broad phases aligned.
fn triangle_projection_interval(tri: [Point3; 3], unit: Vector3, tol: f64) -> (f64, f64) {
    let a = tri[0].to_vec().dot(unit);
    let b = tri[1].to_vec().dot(unit);
    let c = tri[2].to_vec().dot(unit);
    (
        f64::min(f64::min(a, b), c) - tol,
        f64::max(f64::max(a, b), c) + tol,
    )
}

#[test]
fn sweep_coordinate_order_is_total_and_preserves_signed_zero() {
    let values = [
        f64::NEG_INFINITY,
        -1.0,
        -0.0,
        0.0,
        1.0,
        f64::INFINITY,
        f64::NAN,
        f64::from_bits(0xfff8_0000_0000_0001),
    ];
    assert_eq!(compare_sweep_coordinate(-0.0, 0.0), Ordering::Equal);
    for &a in &values {
        for &b in &values {
            assert_eq!(
                sweep_coordinate_key(a).cmp(&sweep_coordinate_key(b)),
                compare_sweep_coordinate(a, b)
            );
            assert_eq!(
                compare_sweep_coordinate(a, b),
                compare_sweep_coordinate(b, a).reverse()
            );
            for &c in &values {
                if compare_sweep_coordinate(a, b) != Ordering::Greater
                    && compare_sweep_coordinate(b, c) != Ordering::Greater
                {
                    assert_ne!(compare_sweep_coordinate(a, c), Ordering::Greater);
                }
            }
        }
    }
}

#[test]
fn sortable_sweep_key_matches_total_float_order_for_arbitrary_bits() {
    let mut state = 0xd4a3_01cf_9412_8abd_u64;
    let mut values = Vec::with_capacity(4096);
    for _ in 0..4096 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        values.push(f64::from_bits(state));
    }
    let mut original = values.clone();
    original.sort_by(|a, b| compare_sweep_coordinate(*a, *b));
    values.sort_by_key(|&value| sweep_coordinate_key(value));
    assert!(
        values
            .iter()
            .zip(original)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

mod active_sweep;
mod collision;
mod in_out_judge;
mod point_cloud;
mod splitting;
mod topology;
mod volume;

pub use collision::Collision;
pub use in_out_judge::{IncludingPointInDomain, SignedRayIndex};
pub use monstertruck_topology::shell::ShellCondition;
pub use point_cloud::{ProjectionNeighborhood, WithPointCloud};
pub use splitting::ExperimentalSplitters;
pub use splitting::Splitting;
pub use topology::Topology;
pub use volume::CalcVolume;
