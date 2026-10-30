// wgpu's `Surface` auto-trait resolution exceeds the default recursion limit of
// 128, which nightly reports as `recursion_depth_exceeding_limit`
// (rust-lang/rust#159228) and will eventually be a hard error. The lint attaches
// to the whole crate and cannot be silenced per function; raising the limit is
// what its own help text suggests.
#![recursion_limit = "256"]
//! Shape and polygon mesh visualization built on `monstertruck-gpu`.
//!
//! # Examples
//!
//! ```
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use monstertruck_gpu::*;
//! use monstertruck_render::*;
//! # let mesh = PolygonMesh::default();
//!
//! let handler = pollster::block_on(DeviceHandler::default_device())?;
//! let mut scene = Scene::new(handler, &Default::default());
//!
//! // Create render instances from a mesh
//! let creator = scene.instance_creator();
//! let instance: PolygonInstance = creator.create_instance(
//!     &mesh,
//!     &PolygonState {
//!         material: Material {
//!             albedo: Vector4::new(0.8, 0.2, 0.2, 1.0),
//!             roughness: 0.3,
//!             ..Default::default()
//!         },
//!         ..Default::default()
//!     },
//! )?;
//! scene.add_object(&instance);
//! # Ok(())
//! # }
//! ```

#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

use bytemuck::{Pod, Zeroable};
use image::DynamicImage;
use monstertruck_gpu::{
    BufferHandler, DeviceHandler, PreBindGroupLayoutEntry, RenderId, Rendered, Scene,
    SceneDescriptor, bind_group_util, impl_render_id, matrix4_to_raw_f32,
    wgpu::{
        BindGroup, BindGroupLayout, BindingResource, BindingType, BlendState, BufferAddress,
        BufferBindingType, BufferUsages, ColorTargetState, ColorWrites, CompareFunction,
        DepthBiasState, DepthStencilState, Device, Extent3d, Face, FilterMode, FragmentState,
        FrontFace, MipmapFilterMode, MultisampleState, Origin3d, PipelineCompilationOptions,
        PipelineLayout, PolygonMode, PrimitiveState, PrimitiveTopology, RenderPipeline,
        RenderPipelineDescriptor, Sampler, SamplerBindingType, SamplerDescriptor, ShaderModule,
        ShaderModuleDescriptor, ShaderSource, ShaderStages, StencilState, TexelCopyBufferLayout,
        TexelCopyTextureInfo, Texture, TextureAspect, TextureDescriptor, TextureDimension,
        TextureFormat, TextureSampleType, TextureUsages, TextureView, TextureViewDescriptor,
        TextureViewDimension, VertexAttribute, VertexBufferLayout, VertexFormat, VertexState,
        VertexStepMode,
    },
};
use std::sync::Arc;

/// Re-exports `monstertruck_mesh`.
pub mod polymesh {
    pub use monstertruck_mesh::*;
}
pub use polymesh::*;

/// Errors while converting CPU-side render data into GPU instances.
#[derive(Debug, thiserror::Error)]
pub enum RenderError {
    /// The source polygon mesh is internally inconsistent.
    #[error(transparent)]
    Mesh(#[from] monstertruck_mesh::errors::Error),
    /// A CPU-side index cannot be represented by wgpu's `u32` index format.
    #[error("{what} index {value} exceeds the u32 GPU index range")]
    IndexTooLarge {
        /// Description of the index being converted.
        what: &'static str,
        /// CPU-side index value.
        value: usize,
    },
    /// Image dimensions cannot be represented by the selected GPU device.
    #[error(
        "unsupported texture dimensions {width}×{height}; both axes must be non-zero and at most {max_dimension}"
    )]
    InvalidTextureDimensions {
        /// Texture width in pixels.
        width: u32,
        /// Texture height in pixels.
        height: u32,
        /// Device limit for a 2D texture axis.
        max_dimension: u32,
    },
    /// A texture row cannot be represented by wgpu's `u32` byte-stride field.
    #[error("texture row is too large: width {width} × {pixel_size} bytes per pixel")]
    TextureRowTooLarge {
        /// Texture width in pixels.
        width: u32,
        /// Size of one pixel in bytes.
        pixel_size: usize,
    },
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "wgpu vertex/material/shader interfaces intentionally use f32 while mesh/modeling data uses f64"
)]
const fn narrow_f64_to_f32(value: f64) -> f32 { value as f32 }

fn gpu_index(value: usize, what: &'static str) -> std::result::Result<u32, RenderError> {
    u32::try_from(value).map_err(|_| RenderError::IndexTooLarge { what, value })
}

/// Material information.
///
/// Each instance is rendered based on the microfacet theory.
#[derive(Debug, Clone, Copy)]
pub struct Material {
    /// albedo, base color, [0, 1]-normalized rgba. Default is `Vector4::new(1.0, 1.0, 1.0, 1.0)`.
    /// Transparent by alpha is not yet supported in the current standard shader.
    pub albedo: Vector4,
    /// roughness of the surface: [0, 1]. Default is 0.5.
    pub roughness: f64,
    /// ratio of specular: [0, 1]. Default is 0.25.
    pub reflectance: f64,
    /// ratio of ambient: [0, 1]. Default is 0.02.
    pub ambient_ratio: f64,
    /// ratio of blending background color: [0, 1]. Default is 0.0.
    pub background_ratio: f64,
    /// alpha blend flag
    pub alpha_blend: bool,
}

/// Configures of instances.
#[derive(Clone, Debug)]
pub struct PolygonState {
    /// instance matrix
    pub matrix: Matrix4,
    /// material of instance
    pub material: Material,
    /// texture of instance
    pub texture: Option<Arc<Texture>>,
    /// If this parameter is true, the backface culling will be activated.
    pub backface_culling: bool,
}

/// Configures of `WireFrameInstance`.
#[derive(Clone, Debug)]
pub struct WireFrameState {
    /// instance matrix
    pub matrix: Matrix4,
    /// color of instance
    pub color: Vector4,
}

/// shaders for rendering polygons
#[derive(Debug, Clone)]
pub struct PolygonShaders {
    vertex_module: Arc<ShaderModule>,
    vertex_entry: &'static str,
    fragment_module: Arc<ShaderModule>,
    fragment_entry: &'static str,
    tex_fragment_module: Arc<ShaderModule>,
    tex_fragment_entry: &'static str,
}

/// shaders for rendering wireframes
#[derive(Debug, Clone)]
pub struct WireShaders {
    vertex_module: Arc<ShaderModule>,
    vertex_entry: &'static str,
    fragment_module: Arc<ShaderModule>,
    fragment_entry: &'static str,
}

/// Instance of polygon
///
/// One can duplicate polygons with different postures and materials
/// that have the same mesh data.
/// To save memory, mesh data on the GPU can be used again.
///
/// The duplicated polygon by `Clone::clone` has the same mesh data and descriptor
/// with original, however, its render id is different from the one of original.
#[derive(Debug)]
pub struct PolygonInstance {
    polygon: (Arc<BufferHandler>, Arc<BufferHandler>),
    state: PolygonState,
    shaders: PolygonShaders,
    id: RenderId,
}

/// Wire frame rendering
#[derive(Debug)]
pub struct WireFrameInstance {
    vertices: Arc<BufferHandler>,
    strips: Arc<BufferHandler>,
    state: WireFrameState,
    shaders: WireShaders,
    id: RenderId,
}

/// Constroctor for instances
#[derive(Debug, Clone)]
pub struct InstanceCreator {
    handler: DeviceHandler,
    polygon_shaders: PolygonShaders,
    wire_shaders: WireShaders,
}

/// for creating `InstanceCreator`
pub trait CreatorCreator {
    /// create `InstanceCreator`
    fn instance_creator(&self) -> InstanceCreator;
}

/// The trait for Buffer Objects.
pub trait CreateBuffers {
    /// Creates buffer handlers of attributes and indices.
    ///
    /// # Errors
    /// Returns [`RenderError`] when source mesh data is invalid or an index cannot fit the GPU index format.
    fn buffers(
        &self,
        vertex_usage: BufferUsages,
        index_usage: BufferUsages,
        device: &Device,
    ) -> std::result::Result<(BufferHandler, BufferHandler), RenderError>;
}

/// The trait for generating `Instance` from `Self`.
pub trait ToInstance<I: Instance> {
    /// Configuration descriptor for instance.
    type State;
    /// Creates `Instance` from `self`.
    ///
    /// # Errors
    /// Returns [`RenderError`] when the source cannot be represented safely as GPU render data.
    fn to_instance(
        &self,
        handler: &DeviceHandler,
        shaders: &I::Shaders,
        desc: &Self::State,
    ) -> std::result::Result<I, RenderError>;
}

/// Instance for rendering
pub trait Instance {
    #[doc(hidden)]
    type Shaders;
    /// Get standard shaders from instance creator.
    #[doc(hidden)]
    fn standard_shaders(creator: &InstanceCreator) -> Self::Shaders;
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Zeroable, Pod)]
struct AttrVertex {
    pub position: [f32; 3],
    pub uv_coord: [f32; 2],
    pub normal: [f32; 3],
}

/// utility for creating `Texture`
pub mod image2texture;
mod instance_creator;
mod instance_descriptor;
mod polygon_instance;
mod polyrend;
mod wireframe_instance;
