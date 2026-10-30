use super::{EuclideanSpace, InnerSpace, Point3, Vector3, Vertex, Zero};
use monstertruck_mesh::Origin;

#[derive(Clone, Copy, Debug)]
pub struct FaceNormal {
    pub face_id: usize,
    pub normal: Vector3,
}

impl FaceNormal {
    pub fn new(positions: &[Point3], face: &[Vertex], face_id: usize) -> Option<Self> {
        if face.len() < 3 || face.iter().any(|vertex| vertex.pos >= positions.len()) {
            return None;
        }
        let center = face.iter().fold(Vector3::zero(), |sum, vertex| {
            sum + positions[vertex.pos].to_vec()
        }) / crate::count_as_f64(face.len());
        let first = face[0];
        let last = face[face.len() - 1];
        let normal = face
            .windows(2)
            .chain(std::iter::once([last, first].as_ref()))
            .fold(Vector3::zero(), |sum, edge| {
                let vec0 = positions[edge[0].pos].to_vec() - center;
                let vec1 = positions[edge[1].pos].to_vec() - center;
                sum + vec0.cross(vec1)
            });
        if normal.magnitude2().so_small() || !normal.magnitude2().is_finite() {
            return None;
        }
        Some(Self {
            face_id,
            normal: normal.normalize(),
        })
    }
}
