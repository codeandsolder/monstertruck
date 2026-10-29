use crate::{
    Arc, AttrVertex, BufferHandler, BufferUsages, CreateBuffers, Device, DeviceHandler, Instance,
    InstanceCreator, Point3, PolygonInstance, PolygonMesh, PolygonShaders, PolygonState,
    PolylineCurve, RenderError, RenderId, StructuredMesh, ToInstance, WireFrameInstance,
    WireFrameState, WireShaders, gpu_index, narrow_f64_to_f32,
};
use itertools::{Itertools, iproduct};

impl CreateBuffers for PolygonMesh {
    #[inline]
    fn buffers(
        &self,
        vertex_usage: BufferUsages,
        index_usage: BufferUsages,
        device: &Device,
    ) -> std::result::Result<(BufferHandler, BufferHandler), RenderError> {
        let expanded = self.expands(|attr| AttrVertex {
            position: [
                narrow_f64_to_f32(attr.position.x),
                narrow_f64_to_f32(attr.position.y),
                narrow_f64_to_f32(attr.position.z),
            ],
            uv_coord: attr.uv_coord.map_or([0.0, 0.0], |v| {
                [narrow_f64_to_f32(v.x), narrow_f64_to_f32(v.y)]
            }),
            normal: attr.normal.map_or([0.0, 0.0, 0.0], |v| {
                [
                    narrow_f64_to_f32(v.x),
                    narrow_f64_to_f32(v.y),
                    narrow_f64_to_f32(v.z),
                ]
            }),
        })?;
        let indices = expanded
            .faces()
            .triangle_iter()
            .flatten()
            .map(|index| gpu_index(index, "polygon triangle"))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        Ok((
            BufferHandler::from_slice(expanded.attributes(), device, vertex_usage),
            BufferHandler::from_slice(&indices, device, index_usage),
        ))
    }
}

impl Instance for PolygonInstance {
    type Shaders = PolygonShaders;
    fn standard_shaders(creator: &InstanceCreator) -> PolygonShaders {
        creator.polygon_shaders.clone()
    }
}

impl ToInstance<PolygonInstance> for PolygonMesh {
    type State = PolygonState;
    #[inline]
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &PolygonShaders,
        state: &PolygonState,
    ) -> std::result::Result<PolygonInstance, RenderError> {
        let (vb, ib) = self.buffers(BufferUsages::VERTEX, BufferUsages::INDEX, handler.device())?;
        Ok(PolygonInstance {
            polygon: (Arc::new(vb), Arc::new(ib)),
            state: state.clone(),
            shaders: shaders.clone(),
            id: RenderId::generate(),
        })
    }
}

impl ToInstance<WireFrameInstance> for PolygonMesh {
    type State = WireFrameState;
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &WireShaders,
        state: &WireFrameState,
    ) -> std::result::Result<WireFrameInstance, RenderError> {
        let device = handler.device();
        let positions: Vec<[f32; 3]> = self
            .positions()
            .iter()
            .map(|p| {
                [
                    narrow_f64_to_f32(p.x),
                    narrow_f64_to_f32(p.y),
                    narrow_f64_to_f32(p.z),
                ]
            })
            .collect();
        let mut strips = Vec::<u32>::new();
        for face in self.faces().face_iter() {
            for (a, b) in face.iter().circular_tuple_windows() {
                if a.pos >= positions.len() {
                    return Err(monstertruck_mesh::errors::Error::OutOfRange(*a).into());
                }
                if b.pos >= positions.len() {
                    return Err(monstertruck_mesh::errors::Error::OutOfRange(*b).into());
                }
                strips.extend([
                    gpu_index(a.pos, "polygon wireframe")?,
                    gpu_index(b.pos, "polygon wireframe")?,
                ]);
            }
        }
        let vb = BufferHandler::from_slice(&positions, device, BufferUsages::VERTEX);
        let ib = BufferHandler::from_slice(&strips, device, BufferUsages::INDEX);
        Ok(WireFrameInstance {
            vertices: Arc::new(vb),
            strips: Arc::new(ib),
            state: state.clone(),
            shaders: shaders.clone(),
            id: RenderId::generate(),
        })
    }
}

impl CreateBuffers for StructuredMesh {
    #[inline]
    fn buffers(
        &self,
        vertex_usage: BufferUsages,
        index_usage: BufferUsages,
        device: &Device,
    ) -> std::result::Result<(BufferHandler, BufferHandler), RenderError> {
        let (m, n) = structured_mesh_dimensions(self)?;
        let mut vertices = Vec::new();
        iproduct!(0..m, 0..n).for_each(|(i, j)| {
            vertices.push(AttrVertex {
                position: [
                    narrow_f64_to_f32(self.positions()[i][j].x),
                    narrow_f64_to_f32(self.positions()[i][j].y),
                    narrow_f64_to_f32(self.positions()[i][j].z),
                ],
                uv_coord: self
                    .uv_division()
                    .map_or([0.0, 0.0], |(u_division, v_division)| {
                        [
                            narrow_f64_to_f32(u_division[i]),
                            narrow_f64_to_f32(v_division[j]),
                        ]
                    }),
                normal: self.normals().map_or([0.0, 0.0, 0.0], |normals| {
                    [
                        narrow_f64_to_f32(normals[i][j].x),
                        narrow_f64_to_f32(normals[i][j].y),
                        narrow_f64_to_f32(normals[i][j].z),
                    ]
                }),
            });
        });
        let mut indices = Vec::<u32>::new();
        for (i, j) in iproduct!(1..m, 1..n) {
            for index in [
                (i - 1) * n + j - 1,
                i * n + j - 1,
                (i - 1) * n + j,
                (i - 1) * n + j,
                i * n + j - 1,
                i * n + j,
            ] {
                indices.push(gpu_index(index, "structured triangle")?);
            }
        }
        Ok((
            BufferHandler::from_slice(&vertices, device, vertex_usage),
            BufferHandler::from_slice(&indices, device, index_usage),
        ))
    }
}

impl ToInstance<PolygonInstance> for StructuredMesh {
    type State = PolygonState;
    #[inline]
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &PolygonShaders,
        state: &PolygonState,
    ) -> std::result::Result<PolygonInstance, RenderError> {
        let (vb, ib) = self.buffers(BufferUsages::VERTEX, BufferUsages::INDEX, handler.device())?;
        Ok(PolygonInstance {
            polygon: (Arc::new(vb), Arc::new(ib)),
            state: state.clone(),
            shaders: shaders.clone(),
            id: RenderId::generate(),
        })
    }
}

fn structured_mesh_dimensions(
    mesh: &StructuredMesh,
) -> std::result::Result<(usize, usize), RenderError> {
    let positions = mesh.positions();
    let rows = positions.len();
    let columns = positions.first().map_or(0, Vec::len);
    if positions.iter().any(|row| row.len() != columns) {
        return Err(monstertruck_mesh::errors::Error::IrregularArray.into());
    }

    if let Some((u_division, v_division)) = mesh.uv_division()
        && (u_division.len() != rows || v_division.len() != columns)
    {
        return Err(monstertruck_mesh::errors::Error::DifferentLengthArrays.into());
    }

    if let Some(normals) = mesh.normals()
        && (normals.len() != rows || normals.iter().any(|row| row.len() != columns))
    {
        return Err(monstertruck_mesh::errors::Error::DifferentLengthArrays.into());
    }

    Ok((rows, columns))
}

fn structured_wire_indices(rows: usize, columns: usize) -> Vec<usize> {
    let mut strips = Vec::new();
    if rows == 0 || columns == 0 {
        return strips;
    }

    // First column: vertical edges between adjacent rows.
    for row in 1..rows {
        strips.extend([(row - 1) * columns, row * columns]);
    }
    // First row: horizontal edges between adjacent columns.
    for column in 1..columns {
        strips.extend([column - 1, column]);
    }
    // Remaining grid edges. Each cell contributes its right and bottom edge;
    // the first row/column above provide the missing outer boundary.
    for row in 1..rows {
        for column in 1..columns {
            strips.extend([
                (row - 1) * columns + column,
                row * columns + column,
                row * columns + column - 1,
                row * columns + column,
            ]);
        }
    }
    strips
}

impl ToInstance<WireFrameInstance> for StructuredMesh {
    type State = WireFrameState;
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &WireShaders,
        state: &WireFrameState,
    ) -> std::result::Result<WireFrameInstance, RenderError> {
        let device = handler.device();
        let (rows, columns) = structured_mesh_dimensions(self)?;
        let positions: Vec<[f32; 3]> = self
            .positions()
            .iter()
            .flatten()
            .map(|p| {
                [
                    narrow_f64_to_f32(p.x),
                    narrow_f64_to_f32(p.y),
                    narrow_f64_to_f32(p.z),
                ]
            })
            .collect();
        let strips = structured_wire_indices(rows, columns)
            .into_iter()
            .map(|index| gpu_index(index, "structured wireframe"))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let vb = BufferHandler::from_slice(&positions, device, BufferUsages::VERTEX);
        let ib = BufferHandler::from_slice(&strips, device, BufferUsages::INDEX);
        Ok(WireFrameInstance {
            vertices: Arc::new(vb),
            strips: Arc::new(ib),
            state: state.clone(),
            shaders: shaders.clone(),
            id: RenderId::generate(),
        })
    }
}

impl ToInstance<WireFrameInstance> for PolylineCurve<Point3> {
    type State = WireFrameState;
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &WireShaders,
        state: &WireFrameState,
    ) -> std::result::Result<WireFrameInstance, RenderError> {
        let device = handler.device();
        let positions: Vec<[f32; 3]> = self
            .iter()
            .map(|p| {
                [
                    narrow_f64_to_f32(p.x),
                    narrow_f64_to_f32(p.y),
                    narrow_f64_to_f32(p.z),
                ]
            })
            .collect();
        let mut strips = Vec::new();
        for index in 1..positions.len() {
            strips.extend([
                gpu_index(index - 1, "polyline")?,
                gpu_index(index, "polyline")?,
            ]);
        }
        let vb = BufferHandler::from_slice(&positions, device, BufferUsages::VERTEX);
        let ib = BufferHandler::from_slice(&strips, device, BufferUsages::INDEX);
        Ok(WireFrameInstance {
            vertices: Arc::new(vb),
            strips: Arc::new(ib),
            state: state.clone(),
            shaders: shaders.clone(),
            id: RenderId::generate(),
        })
    }
}

impl ToInstance<WireFrameInstance> for Vec<PolylineCurve<Point3>> {
    type State = WireFrameState;
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &WireShaders,
        state: &WireFrameState,
    ) -> std::result::Result<WireFrameInstance, RenderError> {
        let device = handler.device();
        let positions: Vec<[f32; 3]> = self
            .iter()
            .flat_map(|poly| poly.iter())
            .map(|p| {
                [
                    narrow_f64_to_f32(p.x),
                    narrow_f64_to_f32(p.y),
                    narrow_f64_to_f32(p.z),
                ]
            })
            .collect();
        let mut strips = Vec::new();
        let mut offset = 0usize;
        for polyline in self {
            for index in 1..polyline.len() {
                strips.extend([
                    gpu_index(offset + index - 1, "polyline collection")?,
                    gpu_index(offset + index, "polyline collection")?,
                ]);
            }
            offset += polyline.len();
        }
        let vb = BufferHandler::from_slice(&positions, device, BufferUsages::VERTEX);
        let ib = BufferHandler::from_slice(&strips, device, BufferUsages::INDEX);
        Ok(WireFrameInstance {
            vertices: Arc::new(vb),
            strips: Arc::new(ib),
            state: state.clone(),
            shaders: shaders.clone(),
            id: RenderId::generate(),
        })
    }
}

#[cfg(test)]
mod structured_wire_index_tests {
    use super::{structured_mesh_dimensions, structured_wire_indices};
    use crate::{Point3, StructuredMesh, Vector3};

    #[test]
    fn structured_wire_indices_cover_grid_without_out_of_range_indices() {
        let indices = structured_wire_indices(3, 3);
        assert_eq!(
            indices,
            vec![
                0, 3, 3, 6, // first-column verticals
                0, 1, 1, 2, // first-row horizontals
                1, 4, 3, 4, 2, 5, 4, 5, // middle row
                4, 7, 6, 7, 5, 8, 7, 8, // bottom row
            ]
        );
        assert!(indices.iter().all(|&index| index < 9));
    }

    #[test]
    fn structured_mesh_dimensions_reject_irregular_unchecked_positions() {
        let mesh = StructuredMesh::from_positions_unchecked(vec![
            vec![Point3::new(0.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ]);
        assert!(structured_mesh_dimensions(&mesh).is_err());
    }

    #[test]
    fn structured_mesh_dimensions_reject_unchecked_uv_shape_mismatch() {
        let mesh = StructuredMesh::from_positions_and_uvs_unchecked(
            vec![vec![Point3::new(0.0, 0.0, 0.0)]],
            (vec![0.0, 1.0], vec![0.0]),
        );
        assert!(structured_mesh_dimensions(&mesh).is_err());
    }

    #[test]
    fn structured_mesh_dimensions_reject_unchecked_normal_shape_mismatch() {
        let mesh = StructuredMesh::from_positions_and_normals_unchecked(
            vec![vec![Point3::new(0.0, 0.0, 0.0)]],
            vec![vec![
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(0.0, 0.0, 1.0),
            ]],
        );
        assert!(structured_mesh_dimensions(&mesh).is_err());
    }
}
