use super::{EuclideanSpace, InnerSpace, Point3, Vector3, Vertex, Zero};

#[derive(Clone, Copy, Debug)]
pub struct FaceNormal {
    /// Only consumed by the `filters` feature's normal-clustering code.
    /// Tagged `dead_code` so `cargo check --no-default-features --features analyzers`
    /// stays clean under the workspace's `-D warnings`.
    #[allow(
        dead_code,
        reason = "face identity is retained with the normal for diagnostics even when a consumer only reads the vector"
    )]
    pub face_id: usize,
    pub normal: Vector3,
}

impl FaceNormal {
    pub fn new(positions: &[Point3], face: &[Vertex], face_id: usize) -> Self {
        let center = face
            .iter()
            .fold(Vector3::zero(), |sum, v| sum + positions[v.pos].to_vec())
            / crate::count_as_f64(face.len());
        let normal = face
            .windows(2)
            .chain(std::iter::once([face[face.len() - 1], face[0]].as_ref()))
            .fold(Vector3::zero(), |sum, v| {
                let vec0 = positions[v[0].pos].to_vec() - center;
                let vec1 = positions[v[1].pos].to_vec() - center;
                sum + vec0.cross(vec1)
            })
            .normalize();
        Self { face_id, normal }
    }
}
