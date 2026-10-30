//! Graphics utility crate built on wgpu.
//!
//! # Examples
//!
//! ```
//! use std::f64::consts::TAU;
//! use monstertruck_core::cgmath64::*;
//! use monstertruck_gpu::*;
//!
//! let handler = pollster::block_on(DeviceHandler::default_device()).expect("GPU device is available");
//! let scene = Scene::new(handler, &SceneDescriptor {
//!     studio: StudioConfig {
//!         camera: Camera {
//!             matrix: Matrix4::look_at_rh(
//!                 Point3::new(1.0, 1.0, 1.0),
//!                 Point3::origin(),
//!                 Vector3::unit_y(),
//!             ).invert().unwrap(),
//!             method: ProjectionMethod::perspective(Rad(TAU / 8.0)),
//!             near_clip: 0.1,
//!             far_clip: 100.0,
//!         },
//!         ..Default::default()
//!     },
//!     ..Default::default()
//! });
//! ```

// wgpu's `Surface` auto-trait resolution exceeds the default recursion limit of
// 128, so `cargo clippy --all-targets -- -D warnings` fails here with
// "overflow evaluating the requirement `wgpu::Surface<'_>: Send`" -- on nightly
// as the `recursion_depth_exceeding_limit` lint (rust-lang/rust#159228), which
// is on its way to becoming a hard error. The lint attaches to the whole crate
// and cannot be silenced per function; raising the limit is what its own help
// text suggests.
#![recursion_limit = "256"]
#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

use bytemuck::{Pod, Zeroable};
use derive_more::{Deref, DerefMut};
use monstertruck_core::{
    bounding_box::BoundingBox,
    cgmath64::{
        Angle, EuclideanSpace, InnerSpace, Matrix3, Matrix4, Point2, Point3, Rad, SquareMatrix,
        Transform, Vector3, ortho, perspective,
    },
};
use std::sync::Arc;
pub use wgpu;
use wgpu::util::{BufferInitDescriptor, DeviceExt};
use wgpu::{
    Adapter, Backends, BindGroup, BindGroupDescriptor, BindGroupEntry, BindGroupLayout,
    BindGroupLayoutDescriptor, BindGroupLayoutEntry, BindingResource, BindingType, Buffer,
    BufferAsyncError, BufferBinding, BufferBindingType, BufferDescriptor, BufferUsages, Color,
    CommandEncoder, CommandEncoderDescriptor, CompositeAlphaMode, CreateSurfaceError,
    CurrentSurfaceTexture, Device, DeviceDescriptor, ExperimentalFeatures, Extent3d, IndexFormat,
    Instance, InstanceDescriptor, Limits, LoadOp, MapMode, MemoryHints, Operations, Origin3d,
    PipelineLayout, PipelineLayoutDescriptor, PollError, PollType, PowerPreference, PresentMode,
    Queue, RenderPassColorAttachment, RenderPassDepthStencilAttachment, RenderPassDescriptor,
    RenderPipeline, RequestAdapterError, RequestAdapterOptions, RequestDeviceError, ShaderStages,
    StoreOp, Surface, SurfaceConfiguration, TexelCopyBufferInfo, TexelCopyBufferLayout,
    TexelCopyTextureInfo, Texture, TextureAspect, TextureDescriptor, TextureDimension,
    TextureFormat, TextureUsages, TextureView, TextureViewDescriptor,
};

#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant as TimeInstant;
#[cfg(target_arch = "wasm32")]
use web_time::Instant as TimeInstant;

/// maximum number of light
pub const LIGHT_MAX: usize = 255;

#[inline]
#[expect(
    clippy::cast_possible_truncation,
    reason = "wgpu shader data is f32; this helper is the intentional IEEE-754 f64-to-f32 narrowing boundary"
)]
pub(crate) const fn narrow_f64_to_f32(value: f64) -> f32 { value as f32 }

/// Narrows an f64 transform matrix to the raw f32 column layout consumed by wgpu.
///
/// The conversion intentionally uses Rust's floating-point narrowing semantics:
/// finite out-of-range values become infinities and NaNs remain NaNs.
#[inline]
#[must_use]
pub fn matrix4_to_raw_f32(matrix: Matrix4) -> [[f32; 4]; 4] {
    [
        matrix.x.map(narrow_f64_to_f32).into(),
        matrix.y.map(narrow_f64_to_f32).into(),
        matrix.z.map(narrow_f64_to_f32).into(),
        matrix.w.map(narrow_f64_to_f32).into(),
    ]
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Zeroable, Pod)]
struct CameraInfo {
    camera_matrix: [[f32; 4]; 4],
    camera_projection: [[f32; 4]; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Zeroable, Pod)]
struct LightInfo {
    position: [f32; 4],
    color: [f32; 4],
    kind: [u32; 4],
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Zeroable, Pod)]
struct SceneInfo {
    background_color: [f32; 4],
    resolution: [u32; 2],
    time: f32,
    num_of_lights: u32,
}

/// Runtime failures from GPU initialization and readback.
#[derive(Debug, thiserror::Error)]
pub enum GpuError {
    /// A window surface could not be created.
    #[error("failed to create GPU surface: {0}")]
    CreateSurface(#[from] CreateSurfaceError),
    /// No compatible adapter could be requested.
    #[error("failed to request GPU adapter: {0}")]
    RequestAdapter(#[from] RequestAdapterError),
    /// The selected adapter could not create a device.
    #[error("failed to request GPU device: {0}")]
    RequestDevice(#[from] RequestDeviceError),
    /// Waiting for submitted GPU work failed.
    #[error("GPU polling failed: {0}")]
    Poll(#[from] PollError),
    /// Mapping a GPU buffer for CPU readback failed.
    #[error("GPU buffer mapping failed: {0}")]
    BufferMap(#[from] BufferAsyncError),
    /// The asynchronous map callback disappeared before delivering a result.
    #[error("GPU buffer mapping callback was dropped")]
    BufferMapChannelClosed,
    /// Window initialization completed without the surface that was requested.
    #[error("window GPU initialization did not produce a surface")]
    MissingSurface,
    /// The public camera transform is singular and cannot be projected.
    #[error("camera transform is singular")]
    SingularCameraMatrix,
    /// A value required by the GPU API does not fit its 32-bit count field.
    #[error("{what} exceeds the GPU API's u32 count range")]
    CountTooLarge {
        /// Name of the count that overflowed.
        what: &'static str,
    },
    /// A requested GPU buffer copy would overrun the destination.
    #[error(
        "destination GPU buffer is too small: source={source_size} bytes, destination={destination_size} bytes"
    )]
    BufferTooSmall {
        /// Source buffer size in bytes.
        source_size: u64,
        /// Destination buffer size in bytes.
        destination_size: u64,
    },
    /// A CPU-visible render readback cannot fit in this platform's address space.
    #[error("render readback requires {bytes} bytes, which cannot fit in memory on this platform")]
    ReadbackTooLarge {
        /// Number of unpadded image bytes requested by the readback.
        bytes: u64,
    },
    /// The configured render texture does not have a simple copyable texel footprint.
    #[error("render texture format {format:?} is not supported for CPU readback")]
    UnsupportedReadbackFormat {
        /// Texture format that cannot be represented by the readback path.
        format: TextureFormat,
    },
}

/// safe handler of GPU buffer
/// [`Buffer`](https://docs.rs/wgpu/0.10.1/wgpu/struct.Buffer.html)
#[derive(Debug)]
pub struct BufferHandler {
    buffer: Buffer,
    size: u64,
    stride: u64,
}

/// Utility for [`BindGroupLayoutEntry`]
///
/// The member variables of this struct are the ones of [`BindGroupLayoutEntry`]
/// with only `binding` removed. We can create `BindGroupLayout` by
/// giving its iterator to the function monstertruck_gpu::[`create_bind_group_layout`].
/// # Examples
/// ```
/// use std::sync::{Arc, Mutex};
/// use monstertruck_gpu::*;
/// use wgpu::*;
/// let handler = pollster::block_on(DeviceHandler::default_device()).expect("GPU device is available");
/// let entries = [
///     PreBindGroupLayoutEntry {
///         visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
///         ty: BindingType::Buffer {
///             ty: BufferBindingType::Uniform,
///             has_dynamic_offset: false,
///             min_binding_size: None,
///         },
///         count: None,
///     },
///     PreBindGroupLayoutEntry {
///         visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
///         ty: BindingType::Buffer {
///             ty: BufferBindingType::Uniform,
///             has_dynamic_offset: false,
///             min_binding_size: None,
///         },
///         count: None,
///     },
/// ];
/// let layout: BindGroupLayout = bind_group_util::create_bind_group_layout(handler.device(), &entries);
/// ```
///
/// [`BindGroupLayoutEntry`]: https://docs.rs/wgpu/0.10.1/wgpu/struct.BindGroupLayoutEntry.html
#[doc(hidden)]
#[derive(Debug)]
pub struct PreBindGroupLayoutEntry {
    pub visibility: ShaderStages,
    pub ty: BindingType,
    pub count: Option<std::num::NonZeroU32>,
}

/// A collection of GPU buffers used by `wgpu` for rendering
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct RenderObject {
    vertex_buffer: Arc<BufferHandler>,
    index_buffer: Option<Arc<BufferHandler>>,
    pipeline: Arc<RenderPipeline>,
    bind_group_layout: Arc<BindGroupLayout>,
    bind_group: Arc<BindGroup>,
    visible: bool,
}

/// setting of camera
#[derive(Clone, Copy, Debug, PartialEq, derive_more::From)]
pub enum ProjectionMethod {
    /// perspective projection
    Perspective {
        /// field of view
        fov: Rad<f64>,
    },
    /// parallel projection
    Parallel {
        /// the size of screen
        screen_size: f64,
    },
}

/// Camera
///
/// A [`Scene`](./struct.Scene.html) holds only one `Camera`.
/// # Examples
/// ## Perspective camera
/// ```
/// use std::f64::consts::PI;
/// use monstertruck_core::{cgmath64::*, tolerance::Tolerance};
/// use monstertruck_gpu::*;
/// let matrix = Matrix4::look_at_rh(
///     Point3::new(1.0, 1.0, 1.0),
///     Point3::origin(),
///     Vector3::new(0.0, 1.0, 0.0),
/// );
/// let camera = Camera {
///     // depends on the difference of the style with cgmath,
///     // the matrix must be inverted
///     matrix: matrix.invert().unwrap(),
///     method: ProjectionMethod::perspective(Rad(PI / 4.0)),
///     near_clip: 0.1,
///     far_clip: 1.0,
/// };
/// assert!(camera.eye_direction().near(&-Vector3::new(1.0, 1.0, 1.0).normalize()));
/// ```
/// ## Parallel camera
/// ```
/// use monstertruck_core::{cgmath64::*, tolerance::Tolerance};
/// use monstertruck_gpu::*;
/// let matrix = Matrix4::look_at_rh(
///     Point3::new(1.0, 1.0, 1.0),
///     Point3::origin(),
///     Vector3::new(0.0, 1.0, 0.0),
/// );
/// let camera = Camera {
///     // depends on the difference of the style with cgmath,
///     // the matrix must be inverted
///     matrix: matrix.invert().unwrap(),
///     method: ProjectionMethod::parallel(1.0),
///     near_clip: 0.1,
///     far_clip: 1.0,
/// };
/// assert!(camera.head_direction().near(&Vector3::new(-0.5, 1.0, -0.5).normalize()));
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Camera {
    /// camera matrix
    ///
    /// This matrix must be in the Euclidean momentum group, the semi-direct product of O(3) and R^3.
    pub matrix: Matrix4,
    /// projection method: perspective or parallel
    pub method: ProjectionMethod,
    /// distance from near clipping plane
    pub near_clip: f64,
    /// distance from far clipping plane
    pub far_clip: f64,
}

/// Rays corresponding to a point on the screen, defined by the camera.
#[derive(Clone, Copy, Debug)]
pub struct Ray {
    origin: Point3,
    direction: Vector3,
}

/// the kinds of light sources: point or uniform
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LightType {
    /// point light source
    Point,
    /// uniform light source
    Uniform,
}

/// Light
///
/// There is no limit to the number of lights that can be added to a [`Scene`](./struct.Scene.html).
/// The information about the lights is sent to the shader as a storage buffer
/// (cf: [`Scene::lights_buffer()`](./struct.Scene.html#method.lights_buffer)).
#[derive(Clone, Debug, PartialEq)]
pub struct Light {
    /// position of light
    pub position: Point3,
    /// [0, 1] range RGB color of light
    pub color: Vector3,
    /// type of light source: point or uniform
    pub kind: LightType,
}

/// Chain that holds [`Device`], [`Queue`] and [`SurfaceConfiguration`].
///
/// This struct is used for creating [`Scene`].
/// [`Device`] and [`Queue`] must be wrapped `Arc`,
/// and [`SurfaceConfiguration`] `Arc<Mutex>`.
///
/// [`Device`]: https://docs.rs/wgpu/0.10.1/wgpu/struct.Device.html
/// [`Queue`]: https://docs.rs/wgpu/0.10.1/wgpu/struct.Queue.html
/// [`SurfaceConfiguration`]: https://docs.rs/wgpu/0.10.1/wgpu/struct.SurfaceConfiguration.html
/// [`Scene`]: ./struct.Scene.html
#[derive(Debug, Clone)]
pub struct DeviceHandler {
    adapter: Adapter,
    device: Device,
    queue: Queue,
}

#[derive(Clone, Debug)]
struct WindowHandler {
    window: Arc<winit::window::Window>,
    surface: Arc<Surface<'static>>,
}

/// The unique Id for `Rendered` struct.
///
/// This structure is not used explicitly by users for modeling by `monstertruck-modeling` and `monstertruck-render`.
/// If you want to define a new drawing element in which "Rendered" will be implemented, you need to add
/// this structure to the member variables of that drawing element.
///
/// This structure is assigned a unique value each time it is generated by `RenderId::generate()`.
/// This property allows us to map a `Rendred` entity to data in GPU memory held by `Scene`.
/// ```
/// use monstertruck_gpu::RenderId;
/// assert_ne!(RenderId::generate(), RenderId::generate());
/// ```
#[derive(Clone, Copy, Hash, PartialEq, Eq, Debug)]
pub struct RenderId(usize);

/// Configuration for studio to shoot the scene.
#[derive(Debug, Clone)]
pub struct StudioConfig {
    /// camera of the scene. Default is `Camera::default()`.
    pub camera: Camera,
    /// All lights in the scene. Default is `vec![Light::default()]`.
    pub lights: Vec<Light>,
    /// background color. Default is `Color::BLACK`.
    pub background: Color,
}

/// Configuration for buffer preparation
#[derive(Clone, Debug, Copy)]
pub struct BackendBufferConfig {
    /// depth test flag. Default is `true`.
    pub depth_test: bool,
    /// sample count for anti-aliasing by MSAA. 1, 2, 4, 8, or 16. Default is `1`.
    pub sample_count: u32,
}

/// Configuration for rendering texture
#[derive(Clone, Debug, Copy)]
pub struct RenderTextureConfig {
    /// canvas size `(width, height)`. Default is `(1024, 768)`.
    pub canvas_size: (u32, u32),
    /// texture format. Default is `TextureFormat::Rgba8Unorm`.
    pub format: TextureFormat,
}

/// Configures of [`Scene`](./struct.Scene.html).
#[derive(Debug, Clone, Default)]
pub struct SceneDescriptor {
    /// Configures for studio i.e. camera, lights, and background.
    pub studio: StudioConfig,
    /// Configures buffer preparation, depth and MSAA.
    pub backend_buffer: BackendBufferConfig,
    /// Configuration for rendering texture
    pub render_texture: RenderTextureConfig,
}

/// Configures of [`WindowScene`](./struct.WindowScene.html).
/// Compared to the structure `SceneDescriptor`, this excludes `render_texture`, which can be obtained from window.
#[derive(Debug, Clone, Default)]
pub struct WindowSceneDescriptor {
    /// Configures for studio i.e. camera, lights, and background.
    pub studio: StudioConfig,
    /// Configures buffer preparation, depth and MSAA.
    pub backend_buffer: BackendBufferConfig,
}

/// Wraps `wgpu` and provides an intuitive graphics API.
///
/// `Scene` is the most important in `monstertruck-gpu`.
/// This structure holds information about rendering and
/// serves as a bridge to the actual rendering of `Rendered` objects.
#[derive(Clone, Debug)]
pub struct Scene {
    device_handler: DeviceHandler,
    objects: SliceHashMap<RenderId, RenderObject>,
    bind_group_layout: BindGroupLayout,
    forward_depth: Option<Texture>,
    sampling_buffer: Option<Texture>,
    descriptor: SceneDescriptor,
    clock: TimeInstant,
}

/// Utility for wrapper
#[derive(Clone, Debug, Deref, DerefMut)]
pub struct WindowScene {
    #[deref]
    #[deref_mut]
    scene: Scene,
    window_handler: WindowHandler,
}

/// Rendered objects in the scene.
pub trait Rendered {
    /// Returns the render id.
    ///
    /// [`RenderId`](./struct.RenderId.html) is a key that maps `self` to a drawing element.
    /// Each object must have a `RenderId` to ensure that there are no duplicates.
    fn render_id(&self) -> RenderId;
    /// Creates the pair (vertex buffer, index buffer).
    fn vertex_buffer(
        &self,
        device_handler: &DeviceHandler,
    ) -> (Arc<BufferHandler>, Option<Arc<BufferHandler>>);
    /// Creates the bind group layout.
    fn bind_group_layout(&self, device_handler: &DeviceHandler) -> Arc<BindGroupLayout>;
    /// Creates the bind group in `set = 1`.
    fn bind_group(
        &self,
        device_handler: &DeviceHandler,
        layout: &BindGroupLayout,
    ) -> Arc<BindGroup>;
    /// Creates the render pipeline.
    fn pipeline(
        &self,
        device_handler: &DeviceHandler,
        layout: &PipelineLayout,
        scene_descriptor: &SceneDescriptor,
    ) -> Arc<RenderPipeline>;
    #[doc(hidden)]
    fn render_object(&self, scene: &Scene) -> RenderObject {
        let (vertex_buffer, index_buffer) = self.vertex_buffer(scene.device_handler());
        let bind_group_layout = self.bind_group_layout(scene.device_handler());
        let bind_group = self.bind_group(scene.device_handler(), &bind_group_layout);
        let pipeline_layout = scene
            .device()
            .create_pipeline_layout(&PipelineLayoutDescriptor {
                bind_group_layouts: &[Some(&scene.bind_group_layout), Some(&bind_group_layout)],
                ..Default::default()
            });
        let pipeline = self.pipeline(scene.device_handler(), &pipeline_layout, &scene.descriptor);
        RenderObject {
            vertex_buffer,
            index_buffer,
            bind_group_layout,
            bind_group,
            pipeline,
            visible: true,
        }
    }
}

mod buffer_handler;
mod camera;
mod light;
#[doc(hidden)]
pub mod rendered_macros;
mod scene;
mod slice_hashmap;
use slice_hashmap::SliceHashMap;

#[doc(hidden)]
pub mod bind_group_util {
    #[inline]
    fn binding_index(index: usize) -> u32 {
        assert!(
            u32::try_from(index).is_ok(),
            "bind-group binding index exceeds wgpu's u32 range"
        );
        #[expect(
            clippy::cast_possible_truncation,
            reason = "the preceding assertion proves the usize binding index fits u32"
        )]
        let index = index as u32;
        index
    }

    use crate::{
        BindGroup, BindGroupDescriptor, BindGroupEntry, BindGroupLayout, BindGroupLayoutDescriptor,
        BindGroupLayoutEntry, BindingResource, Device, PreBindGroupLayoutEntry,
    };
    #[doc(hidden)]
    pub fn create_bind_group<'a, T: IntoIterator<Item = BindingResource<'a>>>(
        device: &Device,
        layout: &BindGroupLayout,
        resources: T,
    ) -> BindGroup {
        let entries: &Vec<BindGroupEntry<'_>> = &resources
            .into_iter()
            .enumerate()
            .map(move |(i, resource)| BindGroupEntry {
                binding: binding_index(i),
                resource,
            })
            .collect();
        device.create_bind_group(&BindGroupDescriptor {
            layout,
            entries,
            label: None,
        })
    }

    #[doc(hidden)]
    pub fn create_bind_group_layout<'a, T: IntoIterator<Item = &'a PreBindGroupLayoutEntry>>(
        device: &Device,
        entries: T,
    ) -> BindGroupLayout {
        let vec: Vec<_> = entries
            .into_iter()
            .enumerate()
            .map(|(i, e)| BindGroupLayoutEntry {
                binding: binding_index(i),
                visibility: e.visibility,
                ty: e.ty,
                count: e.count,
            })
            .collect();
        device.create_bind_group_layout(&BindGroupLayoutDescriptor {
            label: None,
            entries: &vec,
        })
    }
}
