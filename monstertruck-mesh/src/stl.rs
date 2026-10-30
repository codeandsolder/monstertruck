use crate::{
    Deserialize, Faces, InnerSpace, Point3, PolygonMesh, Serialize, StandardAttributes,
    StandardVertex, TOLERANCE, Vector3, array, errors, faces,
};
use bytemuck::{Pod, Zeroable};
use num_traits::ToPrimitive;
use rustc_hash::FxHashMap as HashMap;
use std::io::{BufRead, BufReader, Lines, Read, Write};

const FACESIZE: usize = size_of::<StlFace>();
const CHUNKSIZE: usize = FACESIZE + 2;

type Vertex = StandardVertex;
type Result<T> = std::result::Result<T, errors::Error>;

fn invalid_data(message: impl Into<String>) -> errors::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, message.into()).into()
}

fn stl_f32(value: f64) -> f32 {
    value.to_f32().unwrap_or_else(|| {
        if value.is_nan() {
            f32::NAN
        } else if value.is_sign_negative() {
            f32::NEG_INFINITY
        } else {
            f32::INFINITY
        }
    })
}

fn quantize_coordinate(value: f32) -> i64 {
    let scaled = TOLERANCE.mul_add(0.25, f64::from(value)) / (TOLERANCE * 0.5);
    scaled.to_i64().unwrap_or_else(|| {
        if scaled.is_nan() {
            0
        } else if scaled.is_sign_negative() {
            i64::MIN
        } else {
            i64::MAX
        }
    })
}

fn dequantize_coordinate(value: i64) -> f64 {
    value.to_f64().unwrap_or_else(|| {
        if value.is_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        }
    }) * TOLERANCE
        * 0.5
}

fn syntax_error() -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, "syntax error")
}

/// STL naive mesh.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize, Pod, Zeroable)]
pub struct StlFace {
    /// normal vector
    pub normal: [f32; 3],
    /// the positions of vertices
    pub vertices: [[f32; 3]; 3],
}

impl StlFace {
    #[inline]
    fn is_empty(&self) -> bool { self == &Self::default() }
}

/// STL reading iterator.
#[derive(Debug)]
pub enum StlReader<R: Read> {
    #[doc(hidden)]
    Ascii(Lines<BufReader<R>>),
    #[doc(hidden)]
    Binary(R, usize),
}

/// STL type.
#[derive(Copy, Clone, Debug, Default, Serialize, Deserialize)]
pub enum StlType {
    /// Determine STL type automatically.
    ///
    /// # Reading
    /// If the first 5 bytes are..
    /// - "solid" => ascii format
    /// - otherwise => binary format
    ///
    /// # Writing
    /// Always binary format.
    #[default]
    Automatic,
    /// ASCII format.
    Ascii,
    /// Binary format.
    Binary,
}

impl<R: Read> StlReader<R> {
    #[inline]
    fn text_reader(reader: R) -> Self { Self::Ascii(BufReader::new(reader).lines()) }
    fn binary_reader(mut reader: R, header_judge: bool) -> Result<Self> {
        let mut header = [0; 5];
        reader.read_exact(&mut header)?;
        if header_judge && &header == b"solid" {
            return Ok(Self::text_reader(reader));
        }
        let mut header = [0; 75];
        reader.read_exact(&mut header)?;
        let mut length_bytes = [0; 4];
        reader.read_exact(&mut length_bytes)?;
        let length = usize::try_from(u32::from_le_bytes(length_bytes))
            .map_err(|_| invalid_data("binary STL face count does not fit usize"))?;
        Ok(Self::Binary(reader, length))
    }
    /// Creates new STL reader.
    #[inline]
    /// # Errors
    /// Returns an error if the STL header or count cannot be read or represented on this target.
    pub fn new(reader: R, stl_type: StlType) -> Result<Self> {
        match stl_type {
            StlType::Automatic => Self::binary_reader(reader, true),
            StlType::Binary => Self::binary_reader(reader, false),
            StlType::Ascii => Ok(Self::text_reader(reader)),
        }
    }
    /// Returns the STL type.
    #[inline]
    pub const fn stl_type(&self) -> StlType {
        match self {
            Self::Ascii(_) => StlType::Ascii,
            Self::Binary(_, _) => StlType::Binary,
        }
    }
}

impl<R: Read> Iterator for StlReader<R> {
    type Item = Result<StlFace>;
    fn next(&mut self) -> Option<Self::Item> {
        let res = match self {
            Self::Binary(reader, length) => {
                if *length == 0 {
                    Ok(None)
                } else {
                    *length -= 1;
                    binary_one_read(reader)
                }
            }
            Self::Ascii(lines) => ascii_one_read(lines),
        };
        match res {
            Ok(Some(got)) => Some(Ok(got)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

fn ascii_one_read<R: BufRead>(lines: &mut Lines<R>) -> Result<Option<StlFace>> {
    let mut face = StlFace::default();
    let mut num_ver = 0;
    loop {
        let line = if let Some(got) = lines.next() {
            got?
        } else {
            if face.is_empty() {
                return Ok(None);
            }
            return Err(syntax_error().into());
        };
        let line = line.trim();
        if line.len() < 8 {
            continue;
        }
        if line.starts_with("facet") {
            let args: Vec<_> = line.split_whitespace().collect();
            face.normal = array![i => args[i + 2].parse::<f32>()?; 3];
        } else if line.starts_with("vertex") {
            let args: Vec<_> = line.split_whitespace().collect();
            if num_ver > 2 {
                return Err(syntax_error().into());
            }
            face.vertices[num_ver] = array![i => args[i + 1].parse::<f32>()?; 3];
            num_ver += 1;
        } else if line.starts_with("endfacet") {
            if num_ver != 3 {
                return Err(syntax_error().into());
            }
            return Ok(Some(face));
        }
    }
}

fn binary_one_read<R: Read>(reader: &mut R) -> Result<Option<StlFace>> {
    let mut chunk = [0; CHUNKSIZE];
    reader.read_exact(&mut chunk)?;
    let mut buf = [0; FACESIZE];
    buf.copy_from_slice(&chunk[..FACESIZE]);
    Ok(Some(bytemuck::cast(buf)))
}

/// Write STL file in `stl_type` format.
///
/// If `stl_type == StlType::Automatic`, write the binary format.
#[inline]
/// # Errors
/// Returns an error on output I/O failure, or when binary STL is asked to encode more than the format count limit.
pub fn write<I: IntoStlIterator, W: Write>(
    iter: I,
    writer: &mut W,
    stl_type: StlType,
) -> Result<()> {
    match stl_type {
        StlType::Ascii => write_ascii(iter, writer),
        _ => write_binary(iter, writer),
    }
}

/// Writes ASCII STL data.
fn write_ascii<I: IntoStlIterator, W: Write>(iter: I, writer: &mut W) -> Result<()> {
    let mut iter = iter.into_iter();
    // Some readers expect the `solid` keyword to be followed by whitespace
    // before the (optional) solid name; emit the trailing space even when
    // the name is empty.
    writer.write_all(b"solid \n")?;
    iter.try_for_each::<_, Result<()>>(|face| {
        writer.write_fmt(format_args!(
            "  facet normal {:e} {:e} {:e}\n",
            face.normal[0], face.normal[1], face.normal[2]
        ))?;
        writer.write_all(b"    outer loop\n")?;
        face.vertices.iter().try_for_each(|pt| {
            writer.write_fmt(format_args!(
                "      vertex {:e} {:e} {:e}\n",
                pt[0], pt[1], pt[2]
            ))
        })?;
        writer.write_all(b"    endloop\n  endfacet\n")?;
        Ok(())
    })?;
    writer.write_all(b"endsolid\n")?;
    Ok(())
}

/// Writes binary STL data.
#[inline]
fn write_binary<I: IntoStlIterator, W: Write>(iter: I, writer: &mut W) -> Result<()> {
    let mut iter = iter.into_iter();
    let len = u32::try_from(iter.len())
        .map_err(|_| invalid_data("binary STL supports at most u32::MAX faces"))?;
    writer.write_all(&[0u8; 80])?;
    writer.write_all(&len.to_le_bytes())?;
    iter.try_for_each(|face| {
        writer.write_all(bytemuck::cast_slice(&[face]))?;
        writer.write_all(&[0u8, 0u8])?;
        Ok(())
    })
}

/// By implementing [`IntoStlIterator`] for a type you define how it will be
/// converted to an iterator.
///
/// This is common for types which describe a collection of some kind.
pub trait IntoStlIterator {
    /// Which kind of iterator are we turning this into?
    type IntoIter: ExactSizeIterator<Item = StlFace>;
    /// Creates an iterator from a value.
    fn into_iter(self) -> Self::IntoIter;
}

/// Generate an STL faces from from a [`PolygonMesh`].
#[derive(Debug)]
pub struct PolygonMeshStlFaceIterator<'a> {
    positions: &'a [Point3],
    faces: faces::TriangleIterator<'a, Vertex>,
    len: usize,
}

impl Iterator for PolygonMeshStlFaceIterator<'_> {
    type Item = StlFace;
    fn next(&mut self) -> Option<StlFace> {
        self.faces.next().map(|face| {
            let p = array![i => self.positions[face[i].pos]; 3];
            let n = (p[1] - p[0]).cross(p[2] - p[0]).normalize();
            let normal = [stl_f32(n.x), stl_f32(n.y), stl_f32(n.z)];
            let vertices = array![i => [stl_f32(p[i].x), stl_f32(p[i].y), stl_f32(p[i].z)]; 3];
            StlFace { normal, vertices }
        })
    }
    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) { (self.len, Some(self.len)) }
}

impl ExactSizeIterator for PolygonMeshStlFaceIterator<'_> {}

impl<'a> IntoStlIterator for &'a PolygonMesh {
    type IntoIter = PolygonMeshStlFaceIterator<'a>;
    fn into_iter(self) -> Self::IntoIter {
        let iter = self.faces().triangle_iter();
        Self::IntoIter {
            positions: self.positions(),
            len: iter.len(),
            faces: iter,
        }
    }
}

impl<I> IntoStlIterator for I
where
    I: IntoIterator<Item = StlFace>,
    I::IntoIter: ExactSizeIterator,
{
    type IntoIter = I::IntoIter;
    fn into_iter(self) -> I::IntoIter { self.into_iter() }
}

fn signup_vector(vector: [f32; 3], map: &mut HashMap<[i64; 3], usize>) -> usize {
    let vector = array![i => quantize_coordinate(vector[i]); 3];
    let len = map.len();
    *map.entry(vector).or_insert_with(|| len)
}

fn decode_vector<T: From<[f64; 3]>>((code, _): ([i64; 3], usize)) -> T {
    array![i => dequantize_coordinate(code[i]); 3].into()
}

impl FromIterator<StlFace> for PolygonMesh {
    fn from_iter<I: IntoIterator<Item = StlFace>>(iter: I) -> Self {
        let mut positions = HashMap::<[i64; 3], usize>::default();
        let mut normals = HashMap::<[i64; 3], usize>::default();
        let closure = |face: StlFace| {
            let n = signup_vector(face.normal, &mut normals);
            let p = array![i => signup_vector(face.vertices[i], &mut positions); 3];
            array![i => (p[i], None, Some(n)).into(); 3]
        };
        let faces: Vec<[Vertex; 3]> = iter.into_iter().map(closure).collect();
        let faces = Faces::from_tri_and_quad_faces(faces, Vec::new());
        let mut positions: Vec<([i64; 3], usize)> = positions.into_iter().collect();
        positions.sort_by_key(|a| a.1);
        let positions: Vec<Point3> = positions.into_iter().map(decode_vector).collect();
        let mut normals: Vec<([i64; 3], usize)> = normals.into_iter().collect();
        normals.sort_by_key(|a| a.1);
        let normals: Vec<Vector3> = normals.into_iter().map(decode_vector).collect();
        Self::new_unchecked(
            StandardAttributes {
                positions,
                uv_coords: Vec::new(),
                normals,
            },
            faces,
        )
    }
}

/// Read STL file and parse to [`PolygonMesh`].
#[inline]
/// # Errors
/// Returns an error if the STL stream is malformed, truncated, or cannot be read.
pub fn read<R: Read>(reader: R, stl_type: StlType) -> Result<PolygonMesh> {
    StlReader::new(reader, stl_type)?.collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stl_float_narrowing_preserves_nonfinite_values() {
        assert!(stl_f32(f64::NAN).is_nan());
        assert_eq!(stl_f32(f64::INFINITY), f32::INFINITY);
        assert_eq!(stl_f32(f64::NEG_INFINITY), f32::NEG_INFINITY);
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn binary_writer_rejects_face_counts_above_format_limit_before_writing() {
        #[derive(Clone, Copy)]
        struct TooManyFaces;

        impl Iterator for TooManyFaces {
            type Item = StlFace;

            fn next(&mut self) -> Option<Self::Item> { None }
            fn size_hint(&self) -> (usize, Option<usize>) {
                let len = usize::try_from(u64::from(u32::MAX) + 1).unwrap();
                (len, Some(len))
            }
        }

        impl ExactSizeIterator for TooManyFaces {}

        let mut output = Vec::new();
        let error = write(TooManyFaces, &mut output, StlType::Binary).unwrap_err();
        assert!(error.to_string().contains("at most"));
        assert!(output.is_empty());
    }
}
