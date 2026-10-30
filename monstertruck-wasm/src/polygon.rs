use crate::{IntoWasm, Shell, Solid};
use derive_more::{Deref, DerefMut, From, Into};
use monstertruck_meshing::prelude::{Vector2, Vector3, Zero, obj, stl};
use wasm_bindgen::prelude::wasm_bindgen;

/// Wasm wrapper by Polygonmesh
#[wasm_bindgen]
#[derive(Clone, Debug, Into, From, Deref, DerefMut)]
pub struct PolygonMesh(monstertruck_meshing::prelude::PolygonMesh);

impl IntoWasm for monstertruck_meshing::prelude::PolygonMesh {
    type WasmWrapper = PolygonMesh;
}

/// STL type.
#[wasm_bindgen]
#[repr(C)]
#[derive(Copy, Clone, Debug)]
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
    Automatic,
    /// ASCII format.
    Ascii,
    /// Binary format.
    Binary,
}

impl From<StlType> for stl::StlType {
    fn from(stl_type: StlType) -> Self {
        match stl_type {
            StlType::Automatic => Self::Automatic,
            StlType::Ascii => Self::Ascii,
            StlType::Binary => Self::Binary,
        }
    }
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "the JavaScript/WebGL interleaved vertex buffer is explicitly f32, so f64 geometry is intentionally narrowed at this serialization boundary"
)]
const fn render_scalar(value: f64) -> f32 { value as f32 }

/// Buffer for rendering polygon
#[wasm_bindgen]
#[derive(Debug, Clone, Default)]
pub struct PolygonBuffer {
    vertices: Vec<f32>,
    indices: Vec<u32>,
}

#[wasm_bindgen]
impl PolygonMesh {
    /// input from obj format
    #[must_use]
    pub fn from_obj(data: &[u8]) -> Option<Self> {
        obj::read::<&[u8]>(data)
            .map_err(|e| {
                gloo::console::error!(format!("{e}"));
            })
            .ok()
            .map(IntoWasm::into_wasm)
    }
    /// input from STL format
    #[must_use]
    pub fn from_stl(data: &[u8], stl_type: StlType) -> Option<Self> {
        stl::read::<&[u8]>(data, stl_type.into())
            .map_err(|e| {
                gloo::console::error!(format!("{e}"));
            })
            .ok()
            .map(IntoWasm::into_wasm)
    }
    /// output obj format
    #[must_use]
    pub fn to_obj(&self) -> Option<Vec<u8>> {
        let mut res = Vec::new();
        obj::write(&self.0, &mut res)
            .map_err(|e| {
                gloo::console::error!(format!("{e}"));
            })
            .ok()?;
        Some(res)
    }
    /// output stl format
    #[must_use]
    pub fn to_stl(&self, stl_type: StlType) -> Option<Vec<u8>> {
        let mut res = Vec::new();
        stl::write(&self.0, &mut res, stl_type.into())
            .map_err(|e| {
                gloo::console::error!(format!("{e}"));
            })
            .ok()?;
        Some(res)
    }
    /// Returns polygon buffer
    #[must_use]
    pub fn to_buffer(&self) -> Option<PolygonBuffer> {
        let exp = self
            .0
            .expands(|attr| {
                let position = attr.position;
                let uv_coord = attr.uv_coord.unwrap_or_else(Vector2::zero);
                let normal = attr.normal.unwrap_or_else(Vector3::zero);
                [
                    render_scalar(position[0]),
                    render_scalar(position[1]),
                    render_scalar(position[2]),
                    render_scalar(uv_coord[0]),
                    render_scalar(uv_coord[1]),
                    render_scalar(normal[0]),
                    render_scalar(normal[1]),
                    render_scalar(normal[2]),
                ]
            })
            .ok()?;
        Some(PolygonBuffer {
            vertices: exp.attributes().iter().flatten().copied().collect(),
            indices: exp
                .faces()
                .triangle_iter()
                .flatten()
                .map(u32::try_from)
                .collect::<Result<Vec<_>, _>>()
                .ok()?,
        })
    }
    /// meshing shell
    #[must_use]
    pub fn from_shell(shell: &Shell, tol: f64) -> Self { shell.to_polygon(tol) }
    /// meshing solid
    #[must_use]
    pub fn from_solid(solid: &Solid, tol: f64) -> Self { solid.to_polygon(tol) }
    /// Returns the bonding box
    #[must_use]
    pub fn bounding_box(&self) -> Vec<f64> {
        let bdd = self.0.bounding_box();
        let min = bdd.min();
        let max = bdd.max();
        vec![min[0], min[1], min[2], max[0], max[1], max[2]]
    }
    /// merge two polygons: `self` and `other`.
    pub fn merge(&mut self, other: Self) { self.0.merge(other.0); }
}

#[wasm_bindgen]
impl PolygonBuffer {
    /// vertex buffer. One attribute contains `position: [f32; 3]`, `uv_coord: [f32; 2]` and `normal: [f32; 3]`.
    #[must_use]
    pub fn vertex_buffer(&self) -> Vec<f32> { self.vertices.clone() }
    /// the length (bytes) of vertex buffer. (Num of attributes) * 8 components * 4 bytes.
    #[must_use]
    #[expect(
        clippy::missing_const_for_fn,
        reason = "wasm-bindgen rejects exported const functions"
    )]
    pub fn vertex_buffer_size(&self) -> usize { self.vertices.len() * 4 }
    /// index buffer. `u32`.
    #[must_use]
    pub fn index_buffer(&self) -> Vec<u32> { self.indices.clone() }
    /// the length (bytes) of index buffer. (Num of triangles) * 3 vertices * 4 bytes.
    #[must_use]
    #[expect(
        clippy::missing_const_for_fn,
        reason = "wasm-bindgen rejects exported const functions"
    )]
    pub fn index_buffer_size(&self) -> usize { self.indices.len() * 4 }
}
