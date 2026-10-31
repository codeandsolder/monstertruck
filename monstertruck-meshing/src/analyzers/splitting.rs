use super::{
    Adjacency, FaceNormal, Faces, InnerSpace, MetricSpace, Point3, PolygonMesh, StandardAttributes,
    Vector3, Vertex, nonpositive_tolerance,
};
use std::f64::consts::PI;

#[cfg(test)]
use crate::EuclideanSpace;

/// Splitting the faces into several clusters.
pub trait Splitting {
    /// Creates a sub mesh by the face indices.
    /// # Examples
    /// ```
    /// use monstertruck_mesh::*;
    /// use monstertruck_meshing::analyzers::*;
    ///
    /// // cube
    /// let mesh = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.0, 0.0),
    ///             Point3::new(1.0, 0.0, 0.0),
    ///             Point3::new(1.0, 1.0, 0.0),
    ///             Point3::new(0.0, 1.0, 0.0),
    ///             Point3::new(0.0, 0.0, 1.0),
    ///             Point3::new(1.0, 0.0, 1.0),
    ///             Point3::new(1.0, 1.0, 1.0),
    ///             Point3::new(0.0, 1.0, 1.0),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(&[
    ///         &[3, 2, 1, 0], &[0, 1, 5, 4], &[1, 2, 6, 5],
    ///         &[2, 3, 7, 6], &[3, 0, 4, 7], &[4, 5, 6, 7],
    ///     ]),
    /// ).expect("fixture mesh is valid");
    ///
    /// let submesh = mesh.create_mesh_by_face_indices(&[0, 3, 5]);
    ///
    /// // the same attributes vector
    /// assert_eq!(mesh.positions(), submesh.positions());
    ///
    /// let faces0 = Faces::from_iter(&[
    ///     &[3, 2, 1, 0], &[2, 3, 7, 6], &[4, 5, 6, 7],
    /// ]);
    /// assert_eq!(submesh.faces().len(), faces0.len());
    /// assert_eq!(submesh.faces()[0], faces0[0]);
    /// assert_eq!(submesh.faces()[1], faces0[1]);
    /// assert_eq!(submesh.faces()[2], faces0[2]);
    /// ```
    fn create_mesh_by_face_indices(&self, indices: &[usize]) -> PolygonMesh;
    /// Extracts polygons such that there exists its normal is the same as its face normal.
    /// # Arguments
    /// - `tol`: tolerance to be regarded as the same normal as the face normal
    /// # Returns
    /// - The first polygon consists the faces included in planes.
    /// - The second polygon is the extracted remainder.
    /// # Panics
    /// `tol` must be greater than or equal to `0.0`.
    /// # Examples
    /// ```
    /// use monstertruck_mesh::*;
    /// use monstertruck_meshing::analyzers::*;
    /// let mesh = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.5, 0.0),
    ///             Point3::new(0.0, 0.5, 1.0),
    ///             Point3::new(1.0, 0.5, 1.0),
    ///             Point3::new(1.0, 0.5, 0.0),
    ///             Point3::new(0.0, 0.0, 2.0),
    ///             Point3::new(1.0, 0.0, 2.0),
    ///         ],
    ///         normals: vec![
    ///             Vector3::new(0.0, 1.0, 0.0),
    ///             // displaced normals for smooth rendering
    ///             Vector3::new(0.0, 1.0, 1.0).normalize(),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(&[
    ///         &[(0, None, Some(0)), (1, None, Some(0)), (2, None, Some(0)), (3, None, Some(0))],
    ///         &[(2, None, Some(0)), (1, None, Some(0)), (4, None, Some(1)), (5, None, Some(1))],
    ///     ]),
    /// ).expect("fixture mesh is valid");
    /// let (plane, remained) = mesh.extract_planes(KERNEL_EPSILON); // KERNEL_EPSILON == 1.0e-6
    /// assert_eq!(plane.len(), 1); assert_eq!(plane[0], 0);
    /// assert_eq!(remained.len(), 1); assert_eq!(remained[0], 1);
    /// ```
    fn extract_planes(&self, tol: f64) -> (Vec<usize>, Vec<usize>);
    /// Splits into the components.
    /// # Details
    /// Two polygons are considered to be in the same component if they share an edge
    /// whose vertices has the same positions (and normals if `use_normal == true`).
    /// # Examples
    /// ```
    /// use monstertruck_mesh::*;
    /// use monstertruck_meshing::{analyzers::*, filters::*};
    ///
    /// // cube consisting tri_faces
    /// let mut mesh = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.0, 0.0),
    ///             Point3::new(1.0, 0.0, 0.0),
    ///             Point3::new(1.0, 1.0, 0.0),
    ///             Point3::new(0.0, 1.0, 0.0),
    ///             Point3::new(0.0, 0.0, 1.0),
    ///             Point3::new(1.0, 0.0, 1.0),
    ///             Point3::new(1.0, 1.0, 1.0),
    ///             Point3::new(0.0, 1.0, 1.0),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(&[
    ///         &[3, 2, 0], &[1, 0, 2], &[0, 1, 4], &[5, 4, 1],
    ///         &[1, 2, 5], &[6, 5, 2], &[2, 3, 6], &[7, 6, 3],
    ///         &[3, 0, 7], &[4, 7, 0], &[4, 5, 7], &[6, 7, 5],
    ///     ]),
    /// ).expect("fixture mesh is valid");
    ///
    /// // sign up normals
    /// mesh.add_naive_normals(true).put_together_same_attrs(KERNEL_EPSILON);
    ///
    /// // into components with normals
    /// let components = mesh.components(true);
    /// // The number of components is six because the mesh is a cube.
    /// assert_eq!(components.len(), 6);
    ///
    /// // into components without normals
    /// let components = mesh.components(false);
    /// // The number of components is one because the mesh is a connected cube.
    /// assert_eq!(components.len(), 1);
    /// ```
    fn components(&self, use_normal: bool) -> Vec<Vec<usize>>;
}

impl Splitting for PolygonMesh {
    fn create_mesh_by_face_indices(&self, indices: &[usize]) -> PolygonMesh {
        let positions = self.positions().to_vec();
        let uv_coords = self.uv_coords().to_vec();
        let normals = self.normals().to_vec();
        let faces: Faces = indices.iter().map(|i| &self.faces()[*i]).collect();
        Self::new_unchecked(
            StandardAttributes {
                positions,
                uv_coords,
                normals,
            },
            faces,
        )
    }

    fn extract_planes(&self, tol: f64) -> (Vec<usize>, Vec<usize>) {
        nonpositive_tolerance!(tol, 0.0);
        self.faces_into_two_clusters(|face: &[Vertex]| {
            is_in_the_plane(self.positions(), self.normals(), face, tol * tol)
        })
    }

    fn components(&self, use_normal: bool) -> Vec<Vec<usize>> {
        let face_adjacency = self.faces().face_adjacency(use_normal);
        get_components(&face_adjacency)
    }
}

#[doc(hidden)]
pub trait ExperimentalSplitters {
    fn faces_into_two_clusters<F: Fn(&[Vertex]) -> bool>(
        &self,
        func: F,
    ) -> (Vec<usize>, Vec<usize>);
    fn clustering_faces_by_gcurvature(
        &self,
        threshold: f64,
        preferred_upper: bool,
    ) -> (Vec<usize>, Vec<usize>);
    fn get_gcurve(&self) -> Vec<f64>;
}

/// splitting the global faces
impl ExperimentalSplitters for PolygonMesh {
    /// separate all faces into two clusters, one with `func` returning true and the other with false.
    /// # Returns
    /// (the vector of all faces `f` which `func(f) == true`, the one of the other faces)
    fn faces_into_two_clusters<F: Fn(&[Vertex]) -> bool>(
        &self,
        func: F,
    ) -> (Vec<usize>, Vec<usize>) {
        let mut true_faces = Vec::new();
        let mut false_faces = Vec::new();
        for (i, face) in self.face_iter().enumerate() {
            if func(face) {
                true_faces.push(i);
            } else {
                false_faces.push(i);
            }
        }
        (true_faces, false_faces)
    }

    fn clustering_faces_by_gcurvature(
        &self,
        threshold: f64,
        preferred_upper: bool,
    ) -> (Vec<usize>, Vec<usize>) {
        let gcurve = self.get_gcurve();
        self.faces_into_two_clusters(|face: &[Vertex]| {
            is_signed_up_upper(face, &gcurve, preferred_upper, threshold)
        })
    }

    fn get_gcurve(&self) -> Vec<f64> {
        let positions = self.positions();
        let mut angles = vec![0.0; positions.len()];
        let mut weights = vec![0.0; positions.len()];
        for face in self.tri_faces() {
            angles[face[0].pos] += get_angle(positions, face, 0, 1, 2);
            angles[face[1].pos] += get_angle(positions, face, 1, 2, 0);
            angles[face[2].pos] += get_angle(positions, face, 2, 0, 1);
            add_weights(&mut weights, positions, face);
        }
        for face in self.quad_faces() {
            angles[face[0].pos] += get_angle(positions, face, 0, 1, 3);
            angles[face[1].pos] += get_angle(positions, face, 1, 2, 0);
            angles[face[2].pos] += get_angle(positions, face, 2, 3, 1);
            angles[face[3].pos] += get_angle(positions, face, 3, 0, 1);
            add_weights(&mut weights, positions, face);
        }
        for face in self.other_faces() {
            let n = face.len() - 1;
            angles[face[0].pos] += get_angle(positions, face, 0, 1, n);
            for i in 1..n {
                angles[face[i].pos] += get_angle(positions, face, i, i + 1, i - 1);
            }
            angles[face[n].pos] += get_angle(positions, face, n, 0, n - 1);
            add_weights(&mut weights, positions, face);
        }

        angles
            .into_iter()
            .zip(weights)
            .map(|(ang, weight)| PI.mul_add(2.0, -ang) / weight)
            .collect()
    }
}

/// divide the graph to the connected components.
/// # Arguments
/// * adjacency - the adjacency matrix
/// # Return
/// * the list of the indices of faces contained in each components
fn get_components(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut unchecked = vec![true; adjacency.len()];
    let mut components = Vec::new();
    loop {
        let Some(first) = unchecked.iter().position(|x| *x) else {
            return components;
        };
        let mut stack = vec![first];
        let mut component = vec![first];
        unchecked[first] = false;
        while let Some(cursor) = stack.pop() {
            for i in &adjacency[cursor] {
                if unchecked[*i] {
                    unchecked[*i] = false;
                    component.push(*i);
                    stack.push(*i);
                }
            }
        }
        components.push(component);
    }
}

fn is_in_the_plane(positions: &[Point3], normals: &[Vector3], face: &[Vertex], tol2: f64) -> bool {
    let Some(normal) = FaceNormal::new(positions, face, 0).map(|face| face.normal) else {
        return false;
    };
    face.iter().any(|vertex| {
        vertex
            .nor
            .and_then(|idx| normals.get(idx))
            .is_some_and(|candidate| normal.distance2(*candidate) < tol2)
    })
}

fn is_signed_up_upper(
    face: &[Vertex],
    gcurve: &[f64],
    preferred_upper: bool,
    threshold: f64,
) -> bool {
    if preferred_upper {
        face.as_ref().iter().any(|v| gcurve[v.pos] > threshold)
    } else {
        face.as_ref().iter().all(|v| gcurve[v.pos] > threshold)
    }
}

fn get_angle(positions: &[Point3], face: &[Vertex], idx0: usize, idx1: usize, idx2: usize) -> f64 {
    let vec0 = positions[face[idx1].pos] - positions[face[idx0].pos];
    let vec1 = positions[face[idx2].pos] - positions[face[idx0].pos];
    vec0.angle(vec1).0
}

fn add_weights(weights: &mut [f64], positions: &[Point3], face: &[Vertex]) {
    let area = (2..face.len()).fold(0.0, |sum, i| {
        let vec0 = positions[face[i - 1].pos] - positions[face[0].pos];
        let vec1 = positions[face[i].pos] - positions[face[0].pos];
        sum + (vec0.cross(vec1)).magnitude() / 2.0
    }) / (crate::count_as_f64(face.len()));
    for v in face {
        weights[v.pos] += area;
    }
}

#[test]
fn into_components_test() {
    let faces = vec![
        [(0, None, Some(0)), (1, None, Some(1)), (2, None, Some(2))].as_ref(),
        &[(0, None, Some(0)), (5, None, Some(5)), (1, None, Some(1))],
        &[(0, None, Some(7)), (5, None, Some(5)), (6, None, Some(6))],
        &[(5, None, Some(5)), (6, None, Some(6)), (7, None, Some(7))],
        &[
            (1, None, Some(1)),
            (4, None, Some(4)),
            (3, None, Some(3)),
            (2, None, Some(2)),
        ],
    ]
    .into_iter()
    .collect();
    let positions = vec![Point3::origin(); 8];
    let normals = vec![Vector3::unit_x(); 8];
    let mesh = PolygonMesh::new(
        StandardAttributes {
            positions,
            normals,
            ..Default::default()
        },
        faces,
    )
    .expect("fixture mesh is valid");
    let comp = mesh.components(true);
    assert_eq!(comp.len(), 2);
    assert_eq!(comp[0], vec![0, 1, 4]);
    assert_eq!(comp[1], vec![2, 3]);
}
