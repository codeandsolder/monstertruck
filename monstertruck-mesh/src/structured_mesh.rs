use crate::{
    Deserialize, Faces, Point3, PolygonMesh, StandardAttributes, StandardVertex, StructuredMesh,
    Vector2, Vector3, errors,
};
use errors::Error;
type Result<T> = std::result::Result<T, Error>;

impl StandardVertex {
    #[inline]
    const fn tuple(x: usize, uv: bool, nor: bool) -> Self {
        let pos = x;
        let uv = if uv { Some(x) } else { None };
        let nor = if nor { Some(x) } else { None };
        Self { pos, uv, nor }
    }
}

impl StructuredMesh {
    /// Creates a structured polygon without `uv_division` and `normal`.
    #[inline]
    /// # Errors
    /// Returns an irregular-array error if the position grid is not rectangular.
    pub fn from_positions(positions: Vec<Vec<Point3>>) -> Result<Self> {
        check_matrix_regularity(&positions)?;
        Ok(Self::from_positions_unchecked(positions))
    }

    /// Creates a structured polygon without `uv_division` and `normal`.
    #[inline]
    #[must_use]
    pub const fn from_positions_unchecked(positions: Vec<Vec<Point3>>) -> Self {
        Self {
            positions,
            uv_division: None,
            normals: None,
        }
    }
    /// Creates a structured polygon without normals.
    #[inline]
    /// # Errors
    /// Returns an error if the position grid and UV divisions have incompatible lengths or the divisions are unsorted.
    pub fn from_positions_and_uvs(
        positions: Vec<Vec<Point3>>,
        (u_div, v_div): (Vec<f64>, Vec<f64>),
    ) -> Result<Self> {
        check_matrix_vectors_compatibility(&positions, &u_div, &v_div)?;
        check_vectors_regularity(&u_div, &v_div)?;
        Ok(Self::from_positions_and_uvs_unchecked(
            positions,
            (u_div, v_div),
        ))
    }
    /// Creates a structured polygon without normals.
    #[inline]
    #[must_use]
    pub const fn from_positions_and_uvs_unchecked(
        positions: Vec<Vec<Point3>>,
        uv_divisions: (Vec<f64>, Vec<f64>),
    ) -> Self {
        Self {
            positions,
            uv_division: Some(uv_divisions),
            normals: None,
        }
    }
    /// Creates a structured polygon without uv divisions.
    #[inline]
    /// # Errors
    /// Returns an error if either grid is irregular or the position and normal grids have incompatible shapes.
    pub fn from_positions_and_normals(
        positions: Vec<Vec<Point3>>,
        normals: Vec<Vec<Vector3>>,
    ) -> Result<Self> {
        check_matrix_regularity(&positions)?;
        check_matrices_compatibility(&positions, &normals)?;
        Ok(Self::from_positions_and_normals_unchecked(
            positions, normals,
        ))
    }
    /// Creates a structured polygon without uv divisions.
    #[inline]
    #[must_use]
    pub const fn from_positions_and_normals_unchecked(
        positions: Vec<Vec<Point3>>,
        normals: Vec<Vec<Vector3>>,
    ) -> Self {
        Self {
            positions,
            uv_division: None,
            normals: Some(normals),
        }
    }

    /// Creates new structured mesh.
    /// Checks whether the size of vectors are compatible before creation.
    #[inline]
    /// # Errors
    /// Returns an error if the position or normal grids and UV divisions are incompatible or the divisions are unsorted.
    pub fn new(
        positions: Vec<Vec<Point3>>,
        (u_div, v_div): (Vec<f64>, Vec<f64>),
        normals: Vec<Vec<Vector3>>,
    ) -> Result<Self> {
        check_matrix_vectors_compatibility(&positions, &u_div, &v_div)?;
        check_matrix_vectors_compatibility(&normals, &u_div, &v_div)?;
        check_vectors_regularity(&u_div, &v_div)?;
        Ok(Self::new_unchecked(positions, (u_div, v_div), normals))
    }

    /// Creates new structured mesh.
    /// Does not check whether the size of vectors are compatible before creation.
    #[inline]
    #[must_use]
    pub const fn new_unchecked(
        positions: Vec<Vec<Point3>>,
        uv_division: (Vec<f64>, Vec<f64>),
        normals: Vec<Vec<Vector3>>,
    ) -> Self {
        Self {
            positions,
            uv_division: Some(uv_division),
            normals: Some(normals),
        }
    }

    /// Returns the matrix of all positions.
    #[inline]
    #[must_use]
    pub const fn positions(&self) -> &Vec<Vec<Point3>> { &self.positions }

    /// Returns the vector of the mutable references to the rows of the positions matrix.
    #[inline]
    pub fn positions_mut(&mut self) -> Vec<&mut [Point3]> {
        self.positions
            .iter_mut()
            .map(std::convert::AsMut::as_mut)
            .collect()
    }

    /// Returns the divisions of uv coordinates.
    #[inline]
    #[must_use]
    pub fn uv_division(&self) -> Option<(&Vec<f64>, &Vec<f64>)> {
        self.uv_division.as_ref().map(|tuple| (&tuple.0, &tuple.1))
    }

    /// Returns the mutable slice of uv coordinates division.
    #[inline]
    pub fn uv_division_mut(&mut self) -> Option<(&mut [f64], &mut [f64])> {
        self.uv_division
            .as_mut()
            .map(|tuple| (tuple.0.as_mut(), tuple.1.as_mut()))
    }

    /// Returns the matrix of all normals.
    #[inline]
    #[must_use]
    pub const fn normals(&self) -> Option<&Vec<Vec<Vector3>>> { self.normals.as_ref() }

    /// Returns the vector of the mutable references to the rows of the normals matrix.
    #[inline]
    pub fn normals_mut(&mut self) -> Option<Vec<&mut [Vector3]>> {
        self.normals.as_mut().map(|normals| {
            normals
                .iter_mut()
                .map(std::convert::AsMut::as_mut)
                .collect()
        })
    }

    /// Creates new polygon mesh by destructing `self`.
    #[inline]
    pub fn destruct(self) -> PolygonMesh {
        let Self {
            positions,
            uv_division,
            normals,
        } = self;
        let m = positions.len();
        let n = positions.first().map_or(0, Vec::len);
        let positions = positions.into_iter().flatten().collect();
        let uv_coords = uv_division.map_or_else(Vec::new, move |(udiv, vdiv)| {
            udiv.into_iter()
                .flat_map(|u| vdiv.iter().map(move |v| Vector2::new(u, *v)))
                .collect()
        });
        let normals = normals.map_or_else(Vec::new, |n| n.into_iter().flatten().collect());
        let uv = !uv_coords.is_empty();
        let nor = !normals.is_empty();
        let quad_faces: Vec<_> = itertools::iproduct!(1..m, 1..n)
            .map(move |(i, j)| {
                [
                    StandardVertex::tuple((i - 1) * n + j - 1, uv, nor),
                    StandardVertex::tuple(i * n + j - 1, uv, nor),
                    StandardVertex::tuple(i * n + j, uv, nor),
                    StandardVertex::tuple((i - 1) * n + j, uv, nor),
                ]
            })
            .collect();
        let faces = Faces {
            quad_faces,
            ..Default::default()
        };
        PolygonMesh {
            attributes: StandardAttributes {
                positions,
                uv_coords,
                normals,
            },
            faces,
        }
    }
}

impl<'de> Deserialize<'de> for StructuredMesh {
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        #[derive(Deserialize)]
        pub struct StructuredMesh_ {
            positions: Vec<Vec<Point3>>,
            uv_division: Option<(Vec<f64>, Vec<f64>)>,
            normals: Option<Vec<Vec<Vector3>>>,
        }
        let StructuredMesh_ {
            positions,
            uv_division,
            normals,
        } = StructuredMesh_::deserialize(deserializer)?;
        match (uv_division, normals) {
            (Some(uv_division), Some(normals)) => Self::new(positions, uv_division, normals),
            (Some(uv_division), None) => Self::from_positions_and_uvs(positions, uv_division),
            (None, Some(normals)) => Self::from_positions_and_normals(positions, normals),
            (None, None) => Self::from_positions(positions),
        }
        .map_err(serde::de::Error::custom)
    }
}

#[inline]
fn check_matrix_regularity<T>(matrix: &[Vec<T>]) -> Result<()> {
    let Some(first) = matrix.first() else {
        return Ok(());
    };
    for arr in matrix {
        if arr.len() != first.len() {
            return Err(Error::IrregularArray);
        }
    }
    Ok(())
}

#[inline]
fn check_matrices_compatibility<S, T>(matrix0: &[Vec<S>], matrix1: &[Vec<T>]) -> Result<()> {
    if matrix0.len() != matrix1.len() {
        return Err(Error::DifferentLengthArrays);
    }
    let width0 = matrix0.first().map_or(0, Vec::len);
    let width1 = matrix1.first().map_or(0, Vec::len);
    if width0 != width1 {
        return Err(Error::DifferentLengthArrays);
    }
    if matrix0.iter().any(|row| row.len() != width0)
        || matrix1.iter().any(|row| row.len() != width1)
    {
        return Err(Error::DifferentLengthArrays);
    }
    Ok(())
}

#[inline]
fn check_vectors_regularity(vec0: &[f64], vec1: &[f64]) -> Result<()> {
    for i in 1..vec0.len() {
        if vec0[i - 1] > vec0[i] {
            return Err(Error::UnsortedDivision);
        }
    }
    for i in 1..vec1.len() {
        if vec1[i - 1] > vec1[i] {
            return Err(Error::UnsortedDivision);
        }
    }
    Ok(())
}

#[inline]
fn check_matrix_vectors_compatibility<T>(
    matrix: &[Vec<T>],
    vec0: &[f64],
    vec1: &[f64],
) -> Result<()> {
    if matrix.len() != vec0.len() {
        return Err(Error::DifferentLengthArrays);
    }
    for arr in matrix {
        if arr.len() != vec1.len() {
            return Err(Error::DifferentLengthArrays);
        }
    }
    Ok(())
}

#[cfg(test)]
mod totality_tests {
    use super::*;

    #[test]
    fn empty_structured_mesh_destructs_to_empty_polygon() {
        let mesh =
            StructuredMesh::from_positions(Vec::new()).expect("empty structured mesh is valid");
        let polygon = mesh.destruct();
        assert!(polygon.positions().is_empty());
        assert!(polygon.faces().is_empty());
    }
}
