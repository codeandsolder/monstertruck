use crate::{
    Faces, Point3, PolygonMesh, StandardAttributes, StandardVertex, Vector2, Vector3, errors,
};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
type Vertex = StandardVertex;
type Result<T> = std::result::Result<T, errors::Error>;

/// Writes obj data to output stream
/// # Examples
/// ```
/// # fn main() -> anyhow::Result<()> {
/// use monstertruck_mesh::*;
/// let positions = vec![
///     Point3::new(0.0, 0.0, 0.0),
///     Point3::new(1.0, 0.0, 0.0),
///     Point3::new(0.0, 1.0, 0.0),
///     Point3::new(0.0, 0.0, 1.0),
///     Point3::new(1.0, 1.0, 0.0),
///     Point3::new(1.0, 0.0, 1.0),
///     Point3::new(0.0, 1.0, 1.0),
///     Point3::new(1.0, 1.0, 1.0),
/// ];
/// let normals = vec![
///     Vector3::new(1.0, 0.0, 0.0),
///     Vector3::new(0.0, 1.0, 0.0),
///     Vector3::new(0.0, 0.0, 1.0),
///     Vector3::new(-1.0, 0.0, 0.0),
///     Vector3::new(0.0, -1.0, 0.0),
///     Vector3::new(0.0, 0.0, -1.0),
/// ];
/// let faces = Faces::from_iter(&[
///     [(0, None, Some(5)), (1, None, Some(5)), (2, None, Some(5))],
///     [(4, None,Some(5)), (2, None, Some(5)), (1, None, Some(5))],
///     [(1, None, Some(4)), (0, None, Some(4)), (3, None, Some(4))],
///     [(1, None, Some(4)), (3, None, Some(4)), (5, None, Some(4))],
///     [(1, None, Some(0)), (5, None, Some(0)), (4, None, Some(0))],
///     [(4, None, Some(0)), (5, None, Some(0)), (7, None, Some(0))],
///     [(2, None, Some(1)), (4, None, Some(1)), (7, None, Some(1))],
///     [(2, None, Some(1)), (7, None, Some(1)), (6, None, Some(1))],
///     [(0, None, Some(3)), (2, None, Some(3)), (6, None, Some(3))],
///     [(0, None, Some(3)), (6, None, Some(3)), (3, None, Some(3))],
///     [(3, None, Some(2)), (6, None, Some(2)), (7, None, Some(2))],
///     [(3, None, Some(2)), (7, None, Some(2)), (5, None, Some(2))],
/// ]);
/// let mesh = PolygonMesh::new(
///     StandardAttributes {
///         positions,
///         normals,
///         ..Default::default()
///     },
///     faces,
/// ).expect("mesh is valid");
/// obj::write(&mesh, std::fs::File::create("meshdata.obj")?);
/// # Ok(())
/// # }
/// ```
/// # Errors
/// Returns an error if writing the OBJ stream fails.
pub fn write<W: Write>(mesh: &PolygonMesh, writer: W) -> Result<()> {
    sub_write(mesh, &mut BufWriter::new(writer))
}

/// Writes obj data to output stream
/// # Errors
/// Returns an error if writing any OBJ group or mesh data fails.
pub fn write_vec<W: Write>(mesh: &[PolygonMesh], writer: W) -> Result<()> {
    let mut writer = BufWriter::new(writer);
    for (i, mesh) in mesh.iter().enumerate() {
        writer.write_fmt(format_args!("g {i}\n"))?;
        sub_write(mesh, &mut writer)?;
    }
    Ok(())
}

fn write2vec<V: std::ops::Index<usize, Output = f64>, W: Write>(
    writer: &mut BufWriter<W>,
    vecs: &[V],
    prefix: &str,
) -> Result<()> {
    for vec in vecs {
        writer.write_fmt(format_args!("{prefix} {:.10e} {:.10e}\n", vec[0], vec[1]))?;
    }
    Ok(())
}

fn write3vec<V: std::ops::Index<usize, Output = f64>, W: Write>(
    writer: &mut BufWriter<W>,
    vecs: &[V],
    prefix: &str,
) -> Result<()> {
    for vec in vecs {
        writer.write_fmt(format_args!(
            "{prefix} {:.10e} {:.10e} {:.10e}\n",
            vec[0], vec[1], vec[2]
        ))?;
    }
    Ok(())
}

impl Vertex {
    fn write<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        match (self.uv, self.nor) {
            (None, None) => writer.write_fmt(format_args!("{}", self.pos + 1)),
            (Some(uv), None) => writer.write_fmt(format_args!("{}/{}", self.pos + 1, uv + 1)),
            (None, Some(nor)) => writer.write_fmt(format_args!("{}//{}", self.pos + 1, nor + 1)),
            (Some(uv), Some(nor)) => {
                writer.write_fmt(format_args!("{}/{}/{}", self.pos + 1, uv + 1, nor + 1))
            }
        }
    }
}

impl Faces {
    fn write<W: Write>(&self, writer: &mut W) -> Result<()> {
        for face in self.face_iter() {
            writer.write_all(b"f")?;
            for v in face {
                writer.write_all(b" ")?;
                v.write(writer)?;
            }
            writer.write_all(b"\n")?;
        }
        Ok(())
    }
}

fn next_float(args: &mut std::str::SplitWhitespace<'_>) -> Result<f64> {
    let token = args
        .next()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "missing field"))?;
    Ok(token.parse::<f64>()?)
}

fn invalid_data(message: impl Into<String>) -> errors::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, message.into()).into()
}

fn resolve_obj_index(token: Option<&str>, current_len: usize) -> Result<Option<usize>> {
    let Some(token) = token else {
        return Ok(None);
    };
    if token.is_empty() {
        return Ok(None);
    }

    let index = token.parse::<i64>()?;
    if index == 0 {
        return Err(invalid_data("OBJ indices are 1-based; zero is invalid"));
    }
    if index > 0 {
        let zero_based = usize::try_from(index - 1)
            .map_err(|_| invalid_data(format!("OBJ index {index} does not fit usize")))?;
        return Ok(Some(zero_based));
    }

    let distance = usize::try_from(index.unsigned_abs())
        .map_err(|_| invalid_data(format!("OBJ relative index {index} does not fit usize")))?;
    current_len.checked_sub(distance).map(Some).ok_or_else(|| {
        invalid_data(format!(
            "OBJ relative index {index} exceeds the {current_len} values defined so far"
        ))
    })
}

fn sub_write<W: Write>(mesh: &PolygonMesh, writer: &mut BufWriter<W>) -> Result<()> {
    write3vec(writer, mesh.positions(), "v")?;
    write2vec(writer, mesh.uv_coords(), "vt")?;
    write3vec(writer, mesh.normals(), "vn")?;
    mesh.faces.write(writer)
}

/// Reads mesh data from wavefront obj file.
/// # Errors
/// Returns an error for I/O failures, malformed numeric or index fields, invalid relative indices, or out-of-range mesh indices.
pub fn read<R: Read>(reader: R) -> Result<PolygonMesh> {
    let mut positions = Vec::new();
    let mut uv_coords = Vec::new();
    let mut normals = Vec::new();
    let mut faces = Faces::default();
    let reader = BufReader::new(reader);
    for line in reader.lines() {
        let line = line?;
        let mut args = line.split_whitespace();
        if let Some(first_str) = args.next() {
            if first_str == "v" {
                let x = next_float(&mut args)?;
                let y = next_float(&mut args)?;
                let z = next_float(&mut args)?;
                positions.push(Point3::new(x, y, z));
            } else if first_str == "vt" {
                let u = next_float(&mut args)?;
                let v = next_float(&mut args)?;
                uv_coords.push(Vector2::new(u, v));
            } else if first_str == "vn" {
                let x = next_float(&mut args)?;
                let y = next_float(&mut args)?;
                let z = next_float(&mut args)?;
                normals.push(Vector3::new(x, y, z));
            } else if first_str == "f" {
                let mut face = Vec::new();
                for vert_str in args {
                    if vert_str.starts_with('#') {
                        break;
                    }
                    let mut iter = vert_str.split('/');
                    let pos =
                        resolve_obj_index(iter.next(), positions.len())?.ok_or_else(|| {
                            invalid_data("OBJ face vertex is missing its position index")
                        })?;
                    let uv = resolve_obj_index(iter.next(), uv_coords.len())?;
                    let nor = resolve_obj_index(iter.next(), normals.len())?;
                    if iter.next().is_some() {
                        return Err(invalid_data(format!(
                            "OBJ face vertex {vert_str:?} has more than three index components"
                        )));
                    }
                    face.push(Vertex { pos, uv, nor });
                }
                faces.push(face);
            }
        }
    }
    PolygonMesh::new(
        StandardAttributes {
            positions,
            uv_coords,
            normals,
        },
        faces,
    )
}
