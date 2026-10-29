use crate::{Attributes, Faces, PolygonMesh, errors};
use rustc_hash::FxHashMap;
use std::fmt::Debug;
use std::hash::Hash;

impl<V: Copy + Hash + Debug + Eq, A: Attributes<V>> PolygonMesh<V, A> {
    /// Contract attributes and expand polygon.
    ///
    /// # Examples
    /// ```
    /// use monstertruck_mesh::*;
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let polygon = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.0, 0.0),
    ///             Point3::new(1.0, 0.0, 0.0),
    ///             Point3::new(0.0, 1.0, 0.0),
    ///             Point3::new(1.0, 1.0, 0.0),
    ///         ],
    ///         normals: vec![
    ///             Vector3::new(0.0, 0.0, 1.0),
    ///             Vector3::new(0.0, 0.0, -1.0),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(&[
    ///         &[(0, None, Some(0)), (1, None, Some(0)), (2, None, Some(0))],
    ///         &[(3, None, Some(1)), (1, None, Some(1)), (2, None, Some(1))],
    ///     ])
    /// )?;
    /// let expands = polygon
    ///     .expands(|attr| (attr.position, attr.normal.unwrap()))?;
    /// assert_eq!(
    ///     expands,
    ///     PolygonMesh::<usize, Vec<(Point3, Vector3)>>::new(
    ///         vec![
    ///            (Point3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
    ///            (Point3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
    ///            (Point3::new(0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, 1.0)),
    ///            (Point3::new(1.0, 1.0, 0.0), Vector3::new(0.0, 0.0, -1.0)),
    ///            (Point3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.0, -1.0)),
    ///            (Point3::new(0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, -1.0)),
    ///         ],
    ///         Faces::from_iter(&[[0, 1, 2], [3, 4, 5]]),
    ///     )?
    /// );
    /// # Ok(())
    /// # }
    /// ```
    /// # Errors
    /// Returns an out-of-range error if a face references an attribute that is not present.
    pub fn expands<T: Copy>(
        &self,
        contraction: impl Fn(A::Output) -> T,
    ) -> Result<PolygonMesh<usize, Vec<T>>, errors::Error<V>> {
        let mut attributes = Vec::<T>::new();
        let mut vertex_map = FxHashMap::<V, usize>::default();
        let mut faces = Faces::<usize>::default();
        for face in self.face_iter() {
            let mut expanded = Vec::with_capacity(face.len());
            for &vertex in face {
                let idx = if let Some(idx) = vertex_map.get(&vertex).copied() {
                    idx
                } else {
                    let attribute = self
                        .attributes
                        .get(vertex)
                        .ok_or(errors::Error::OutOfRange(vertex))?;
                    let idx = attributes.len();
                    attributes.push(contraction(attribute));
                    vertex_map.insert(vertex, idx);
                    idx
                };
                expanded.push(idx);
            }
            faces.push(expanded);
        }
        Ok(PolygonMesh::new_unchecked(attributes, faces))
    }
}

#[cfg(test)]
mod totality_tests {
    use super::*;
    use crate::{Point3, StandardAttributes, StandardVertex};

    #[test]
    fn expansion_reports_invalid_source_indices() {
        let mesh = PolygonMesh::new_unchecked(
            StandardAttributes {
                positions: vec![Point3::new(0.0, 0.0, 0.0)],
                ..Default::default()
            },
            Faces::from_iter([[
                StandardVertex::from(0usize),
                StandardVertex::from(1usize),
                StandardVertex::from(0usize),
            ]]),
        );
        let result = mesh.expands(|attribute| attribute);
        assert!(matches!(result, Err(errors::Error::OutOfRange(vertex)) if vertex.pos == 1));
    }
}
