use crate::{
    Adapter, Arc, BackendBufferConfig, Backends, BindGroup, BindGroupLayout, BindingType,
    BufferBindingType, BufferDescriptor, BufferHandler, BufferUsages, Camera, Color,
    CommandEncoderDescriptor, CompositeAlphaMode, CurrentSurfaceTexture, Device, DeviceDescriptor,
    DeviceHandler, ExperimentalFeatures, Extent3d, GpuError, IndexFormat, Instance,
    InstanceDescriptor, LIGHT_MAX, Light, LightInfo, Limits, LoadOp, MapMode, MemoryHints,
    Operations, Origin3d, PipelineLayoutDescriptor, PollType, PowerPreference,
    PreBindGroupLayoutEntry, PresentMode, Queue, RenderId, RenderPassColorAttachment,
    RenderPassDepthStencilAttachment, RenderPassDescriptor, RenderTextureConfig, Rendered,
    RequestAdapterOptions, Scene, SceneDescriptor, SceneInfo, ShaderStages, StoreOp, StudioConfig,
    Surface, SurfaceConfiguration, TexelCopyBufferInfo, TexelCopyBufferLayout,
    TexelCopyTextureInfo, Texture, TextureAspect, TextureDescriptor, TextureDimension,
    TextureFormat, TextureUsages, TextureView, TextureViewDescriptor, TimeInstant, WindowHandler,
    WindowScene, WindowSceneDescriptor, Zeroable, bind_group_util, narrow_f64_to_f32,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use winit::window::Window;

static MAXID: AtomicUsize = AtomicUsize::new(0);

impl RenderId {
    /// Generate the unique `RenderId`.
    #[inline]
    pub fn generate() -> Self { Self(MAXID.fetch_add(1, Ordering::Relaxed)) }
}

async fn init_default_device(
    window: Option<Arc<Window>>,
) -> Result<(DeviceHandler, Option<WindowHandler>), GpuError> {
    #[cfg(not(feature = "webgl"))]
    let backends = Backends::PRIMARY;
    #[cfg(feature = "webgl")]
    let backends = Backends::all();
    let instance = Instance::new(InstanceDescriptor {
        backends,
        ..InstanceDescriptor::new_without_display_handle()
    });

    let surface = match window.as_ref() {
        Some(window) => Some(instance.create_surface(Arc::clone(window))?),
        None => None,
    };

    let adapter = instance
        .request_adapter(&RequestAdapterOptions {
            #[cfg(not(feature = "webgl"))]
            power_preference: PowerPreference::HighPerformance,
            #[cfg(feature = "webgl")]
            power_preference: PowerPreference::LowPower,
            compatible_surface: surface.as_ref(),
            force_fallback_adapter: false,
        })
        .await?;

    let (device, queue) = adapter
        .request_device(&DeviceDescriptor {
            required_features: wgpu::Features::default(),
            memory_hints: MemoryHints::Performance,
            #[cfg(not(feature = "webgl"))]
            required_limits: Limits::downlevel_defaults().using_resolution(adapter.limits()),
            #[cfg(feature = "webgl")]
            required_limits: Limits::downlevel_webgl2_defaults(),
            label: None,
            trace: wgpu::Trace::default(),
            experimental_features: ExperimentalFeatures::disabled(),
        })
        .await?;
    let device_handler = DeviceHandler::new(adapter, device, queue);
    let window_handler = match (window, surface) {
        (Some(window), Some(surface)) => Some(WindowHandler {
            window,
            surface: Arc::new(surface),
        }),
        (None, None) => None,
        _ => return Err(GpuError::MissingSurface),
    };
    Ok((device_handler, window_handler))
}

impl DeviceHandler {
    /// constructor
    #[inline]
    #[must_use]
    pub const fn new(adapter: Adapter, device: Device, queue: Queue) -> Self {
        Self {
            adapter,
            device,
            queue,
        }
    }
    /// Returns the reference of the adapter.
    #[inline]
    #[must_use]
    pub const fn adapter(&self) -> &Adapter { &self.adapter }
    /// Returns the reference of the device.
    #[inline]
    #[must_use]
    pub const fn device(&self) -> &Device { &self.device }
    /// Returns the reference of the queue.
    #[inline]
    #[must_use]
    pub const fn queue(&self) -> &Queue { &self.queue }

    /// Creates default device handler.
    ///
    /// # Errors
    /// Returns a [`GpuError`] if no compatible adapter/device can be initialized.
    pub async fn default_device() -> Result<Self, GpuError> {
        Ok(init_default_device(None).await?.0)
    }
}

impl Default for StudioConfig {
    #[inline]
    fn default() -> Self {
        Self {
            background: Color::BLACK,
            camera: Camera::default(),
            lights: vec![Light::default()],
        }
    }
}

impl Default for BackendBufferConfig {
    #[inline]
    fn default() -> Self {
        Self {
            depth_test: true,
            sample_count: 1,
        }
    }
}

impl Default for RenderTextureConfig {
    #[inline]
    fn default() -> Self {
        Self {
            canvas_size: (1024, 768),
            format: TextureFormat::Rgba8Unorm,
        }
    }
}

impl RenderTextureConfig {
    /// Returns compatible `SurfaceConfiguration`.
    #[inline]
    #[must_use]
    pub const fn compatible_surface_config(self) -> SurfaceConfiguration {
        SurfaceConfiguration {
            usage: TextureUsages::RENDER_ATTACHMENT,
            format: self.format,
            width: self.canvas_size.0,
            height: self.canvas_size.1,
            alpha_mode: CompositeAlphaMode::Auto,
            present_mode: PresentMode::Fifo,
            view_formats: Vec::new(),
            desired_maximum_frame_latency: 2,
        }
    }
}

impl SceneDescriptor {
    /// Creates a `UNIFORM` buffer of camera.
    ///
    /// The bind group provides [`Scene`] holds this uniform buffer.
    ///
    /// # Shader Example
    /// ```glsl
    /// layout(set = 0, binding = 0) uniform Camera {
    ///     mat4 camera_matrix;     // the camera matrix
    ///     mat4 camera_projection; // the projection into the normalized view volume
    /// };
    /// ```
    /// Returns `None` when the camera transform is singular.
    #[inline]
    #[must_use]
    pub fn camera_buffer(&self, device: &Device) -> Option<BufferHandler> {
        let (width, height) = self.render_texture.canvas_size;
        let as_rat = f64::from(width) / f64::from(height);
        self.studio.camera.buffer(as_rat, device)
    }

    /// Creates a `STORAGE` buffer of all lights.
    ///
    /// The bind group provides [`Scene`] holds this uniform buffer.
    ///
    /// # Shader Example
    /// ```glsl
    /// struct Light {
    ///     vec4 position;      // the position of light, position.w == 1.0
    ///     vec4 color;         // the color of light, color.w == 1.0
    ///     uvec4 light_type;   // Point => uvec4(0, 0, 0, 0), Uniform => uvec4(1, 0, 0, 0)
    /// };
    ///
    /// layout(set = 0, binding = 1) buffer Lights {
    ///     Light lights[];
    /// };
    /// ```
    #[inline]
    pub fn lights_buffer(&self, device: &Device) -> BufferHandler {
        let mut light_vec: Vec<_> = self.studio.lights.iter().map(Light::light_info).collect();
        light_vec.resize(LIGHT_MAX, LightInfo::zeroed());
        BufferHandler::from_slice(&light_vec, device, BufferUsages::UNIFORM)
    }

    #[inline]
    fn sampling_buffer(
        device: &Device,
        render_texture: RenderTextureConfig,
        sample_count: u32,
    ) -> Texture {
        device.create_texture(&TextureDescriptor {
            size: Extent3d {
                width: render_texture.canvas_size.0,
                height: render_texture.canvas_size.1,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count,
            dimension: TextureDimension::D2,
            format: render_texture.format,
            usage: TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
            label: None,
        })
    }

    #[inline]
    fn depth_texture(device: &Device, size: (u32, u32), sample_count: u32) -> Texture {
        device.create_texture(&TextureDescriptor {
            size: Extent3d {
                width: size.0,
                height: size.1,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count,
            dimension: TextureDimension::D2,
            format: TextureFormat::Depth32Float,
            usage: TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
            label: None,
        })
    }

    fn backend_buffers(&self, device: &Device) -> (Option<Texture>, Option<Texture>) {
        let forward_depth = if self.backend_buffer.depth_test {
            Some(Self::depth_texture(
                device,
                self.render_texture.canvas_size,
                self.backend_buffer.sample_count,
            ))
        } else {
            None
        };
        let sampling_buffer = if self.backend_buffer.sample_count > 1 {
            Some(Self::sampling_buffer(
                device,
                self.render_texture,
                self.backend_buffer.sample_count,
            ))
        } else {
            None
        };
        (forward_depth, sampling_buffer)
    }
}

/// Mutable reference of `SceneDescriptor` in `Scene`.
///
/// When this struct is dropped, the backend buffers of scene will be updated.
#[derive(Debug)]
pub struct SceneDescriptorMut<'a>(&'a mut Scene);

impl std::ops::Deref for SceneDescriptorMut<'_> {
    type Target = SceneDescriptor;
    #[inline]
    fn deref(&self) -> &SceneDescriptor { &self.0.descriptor }
}

impl std::ops::DerefMut for SceneDescriptorMut<'_> {
    #[inline]
    fn deref_mut(&mut self) -> &mut SceneDescriptor { &mut self.0.descriptor }
}

impl Drop for SceneDescriptorMut<'_> {
    fn drop(&mut self) {
        let (forward_depth, sampling_buffer) = self.backend_buffers(self.0.device());
        self.0.forward_depth = forward_depth;
        self.0.sampling_buffer = sampling_buffer;
    }
}

impl Scene {
    #[inline]
    fn camera_bgl_entry() -> PreBindGroupLayoutEntry {
        PreBindGroupLayoutEntry {
            visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
            ty: BindingType::Buffer {
                ty: BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        }
    }

    #[inline]
    fn lights_bgl_entry() -> PreBindGroupLayoutEntry {
        PreBindGroupLayoutEntry {
            visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
            ty: BindingType::Buffer {
                ty: BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        }
    }

    #[inline]
    fn scene_bgl_entry() -> PreBindGroupLayoutEntry {
        PreBindGroupLayoutEntry {
            visibility: ShaderStages::VERTEX | ShaderStages::FRAGMENT,
            ty: BindingType::Buffer {
                ty: BufferBindingType::Uniform,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        }
    }

    #[inline]
    fn init_scene_bind_group_layout(device: &Device) -> BindGroupLayout {
        bind_group_util::create_bind_group_layout(
            device,
            &[
                Self::camera_bgl_entry(),
                Self::lights_bgl_entry(),
                Self::scene_bgl_entry(),
            ],
        )
    }

    /// constructor
    // About `scene_desc`, entity is better than reference for the performance.
    // This is reference because only for as wgpu is.
    #[inline]
    #[must_use]
    pub fn new(device_handler: DeviceHandler, scene_desc: &SceneDescriptor) -> Self {
        let device = device_handler.device();
        let (forward_depth, sampling_buffer) = scene_desc.backend_buffers(device);
        let bind_group_layout = Self::init_scene_bind_group_layout(device);
        Self {
            objects: crate::SliceHashMap::default(),
            bind_group_layout,
            forward_depth,
            sampling_buffer,
            clock: TimeInstant::now(),
            descriptor: scene_desc.clone(),
            device_handler,
        }
    }

    /// Construct scene from default GPU device.
    /// # Arguments
    /// - `size`: (width, height)
    ///
    /// # Errors
    /// Returns a [`GpuError`] if GPU device initialization fails.
    pub async fn from_default_device(scene_desc: &SceneDescriptor) -> Result<Self, GpuError> {
        Ok(Self::new(
            DeviceHandler::default_device().await?,
            scene_desc,
        ))
    }

    /// Creates compatible texture for render attachment.
    ///
    /// # Remarks
    /// The usage of texture is `TextureUsages::RENDER_ATTACHMENT | TetureUsages::COPY_SRC`.
    #[inline]
    #[must_use]
    pub fn compatible_texture(&self) -> Texture {
        let config = self.descriptor.render_texture;
        self.device().create_texture(&TextureDescriptor {
            label: None,
            size: Extent3d {
                width: config.canvas_size.0,
                height: config.canvas_size.1,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: TextureDimension::D2,
            format: config.format,
            view_formats: &[],
            usage: TextureUsages::RENDER_ATTACHMENT | TextureUsages::COPY_SRC,
        })
    }

    /// Returns the reference of its own `DeviceHandler`.
    #[inline]
    #[must_use]
    pub const fn device_handler(&self) -> &DeviceHandler { &self.device_handler }

    /// Returns the reference of the adapter.
    #[inline]
    #[must_use]
    pub const fn dapter(&self) -> &Adapter { &self.device_handler.adapter }

    /// Returns the reference of the device.
    #[inline]
    #[must_use]
    pub const fn device(&self) -> &Device { &self.device_handler.device }

    /// Returns the reference of the queue.
    #[inline]
    #[must_use]
    pub const fn queue(&self) -> &Queue { &self.device_handler.queue }

    /// Returns the elapsed time since the scene was created.
    #[inline]
    #[must_use]
    pub fn elapsed(&self) -> std::time::Duration { self.clock.elapsed() }

    /// Returns the reference of the descriptor.
    #[inline]
    #[must_use]
    pub const fn descriptor(&self) -> &SceneDescriptor { &self.descriptor }

    /// Returns the mutable reference of the descriptor.
    ///
    /// # Remarks
    ///
    /// When the return value is dropped, the depth buffer and sampling buffer are automatically updated.
    /// Use `studio_config_mut` if you only want to update the colors of the camera, lights, and background.
    #[inline]
    pub const fn descriptor_mut(&mut self) -> SceneDescriptorMut<'_> { SceneDescriptorMut(self) }

    /// Returns the reference of the studio configuration.
    #[inline]
    #[must_use]
    pub const fn studio_config(&self) -> &StudioConfig { &self.descriptor.studio }

    /// Returns the mutable reference of the studio configuration.
    #[inline]
    pub const fn studio_config_mut(&mut self) -> &mut StudioConfig { &mut self.descriptor.studio }

    /// Returns the bind group layout in the scene.
    #[inline]
    #[must_use]
    pub const fn bind_group_layout(&self) -> &BindGroupLayout { &self.bind_group_layout }

    /// Creates a `UNIFORM` buffer of the camera.
    ///
    /// The bind group provides [`Scene`] holds this uniform buffer.
    ///
    /// # Shader Example
    /// ```glsl
    /// layout(set = 0, binding = 0) uniform Camera {
    ///     mat4 camera_matrix;     // the camera matrix
    ///     mat4 camera_projection; // the projection into the normalized view volume
    /// };
    /// ```
    /// Returns `None` when the camera transform is singular.
    #[inline]
    #[must_use]
    pub fn camera_buffer(&self) -> Option<BufferHandler> {
        self.descriptor.camera_buffer(self.device())
    }

    /// Creates a `STORAGE` buffer of all lights.
    ///
    /// The bind group provides [`Scene`] holds this uniform buffer.
    ///
    /// # Shader Example
    /// ```glsl
    /// struct Light {
    ///     vec4 position;      // the position of light, position.w == 1.0
    ///     vec4 color;         // the color of light, color.w == 1.0
    ///     uvec4 light_type;   // Point => uvec4(0, 0, 0, 0), Uniform => uvec4(1, 0, 0, 0)
    /// };
    ///
    /// layout(set = 0, binding = 1) buffer Lights {
    ///     Light lights[]; // the number of lights must be gotten from another place
    /// };
    /// ```
    #[inline]
    #[must_use]
    pub fn lights_buffer(&self) -> BufferHandler { self.descriptor.lights_buffer(self.device()) }

    /// Creates a `UNIFORM` buffer of the scene status.
    ///
    /// The bind group provides [`Scene`] holds this uniform buffer.
    ///
    /// # Shader Example
    /// ```glsl
    /// layout(set = 0, binding = 2) uniform Scene {
    ///     vec4 bk_color;  // color of back ground
    ///     float time;     // elapsed time since the scene was created.
    ///     uint nlights;   // the number of lights
    /// };
    /// ```
    ///
    /// # Errors
    /// Returns [`GpuError::CountTooLarge`] if the light count cannot fit the shader's `u32` field.
    #[inline]
    pub fn scene_status_buffer(&self) -> Result<BufferHandler, GpuError> {
        let bk = self.descriptor.studio.background;
        let size = self.descriptor.render_texture.canvas_size;
        let num_of_lights = u32::try_from(self.descriptor.studio.lights.len()).map_err(|_| {
            GpuError::CountTooLarge {
                what: "light count",
            }
        })?;
        let scene_info = SceneInfo {
            background_color: [
                narrow_f64_to_f32(bk.r),
                narrow_f64_to_f32(bk.g),
                narrow_f64_to_f32(bk.b),
                narrow_f64_to_f32(bk.a),
            ],
            resolution: size.into(),
            time: self.elapsed().as_secs_f32(),
            num_of_lights,
        };
        Ok(BufferHandler::from_slice(
            &[scene_info],
            self.device(),
            BufferUsages::UNIFORM,
        ))
    }

    /// Creates bind group.
    /// # Shader Examples
    /// Suppose binded as `set = 0`.
    /// ```glsl
    /// layout(set = 0, binding = 0) uniform Camera {
    ///     mat4 camera_matrix;     // the camera matrix
    ///     mat4 camera_projection; // the projection into the normalized view volume
    /// };
    ///
    /// struct Light {
    ///     vec4 position;      // the position of light, position.w == 1.0
    ///     vec4 color;         // the color of light, color.w == 1.0
    ///     uvec4 light_type;   // Point => uvec4(0, 0, 0, 0), Uniform => uvec4(1, 0, 0, 0)
    /// };
    ///
    /// layout(set = 0, binding = 1) buffer Lights {
    ///     Light lights[];
    /// };
    ///
    /// layout(set = 0, binding = 2) uniform Scene {
    ///     float time;     // elapsed time since the scene was created.
    ///     uint nlights;   // the number of lights
    /// };
    /// ```
    /// # Errors
    /// Returns [`GpuError::SingularCameraMatrix`] when the camera transform cannot be projected.
    #[inline]
    pub fn scene_bind_group(&self) -> Result<BindGroup, GpuError> {
        let camera_buffer = self.camera_buffer().ok_or(GpuError::SingularCameraMatrix)?;
        let scene_status_buffer = self.scene_status_buffer()?;
        Ok(bind_group_util::create_bind_group(
            self.device(),
            &self.bind_group_layout,
            vec![
                camera_buffer.binding_resource(),
                self.lights_buffer().binding_resource(),
                scene_status_buffer.binding_resource(),
            ],
        ))
    }

    /// Adds a render object to the scene.
    ///
    /// If there already exists a render object with the same Id,
    /// replaces the render object and returns false.
    #[inline]
    pub fn add_object<R: Rendered>(&mut self, object: &R) -> bool {
        let render_object = object.render_object(self);
        self.objects
            .insert(object.render_id(), render_object)
            .is_none()
    }
    /// Sets the visibility of a render object.
    ///
    /// If there does not exist the render object in the scene, does nothing and returns `false`.
    #[inline]
    pub fn set_visibility<R: Rendered>(&mut self, object: &R, visible: bool) -> bool {
        self.objects
            .get_mut(&object.render_id())
            .map(|obj| obj.visible = visible)
            .is_some()
    }
    /// Adds render objects to the scene.
    ///
    /// If there already exists a render object with the same Id,
    /// replaces the render object and returns false.
    #[inline]
    pub fn add_objects<'a, R, I>(&mut self, objects: I) -> bool
    where
        R: 'a + Rendered,
        I: IntoIterator<Item = &'a R>, {
        let closure = move |flag, object| flag && self.add_object(object);
        objects.into_iter().fold(true, closure)
    }
    /// Removes a render object from the scene.
    ///
    /// If there does not exist the render object in the scene, does nothing and returns `false`.
    #[inline]
    pub fn remove_object<R: Rendered>(&mut self, object: &R) -> bool {
        self.objects.remove(&object.render_id()).is_some()
    }
    /// Removes render objects from the scene.
    ///
    /// If there exists a render object which does not exist in the scene, returns `false`.
    #[inline]
    pub fn remove_objects<'a, R, I>(&mut self, objects: I) -> bool
    where
        R: 'a + Rendered,
        I: IntoIterator<Item = &'a R>, {
        let closure = move |flag, object| flag && self.remove_object(object);
        objects.into_iter().fold(true, closure)
    }

    /// Removes all render objects from the scene.
    #[inline]
    pub fn clear_objects(&mut self) { self.objects.clear() }

    /// Returns the number of the render objects in the scene.
    #[inline]
    #[must_use]
    pub const fn number_of_objects(&self) -> usize { self.objects.len() }

    /// Synchronizes the information of vertices of `object` in the CPU memory
    /// and that in the GPU memory.
    ///
    /// If there does not exist the render object in the scene, does nothing and returns false.
    #[inline]
    pub fn update_vertex_buffer<R: Rendered>(&mut self, object: &R) -> bool {
        let (handler, objects) = (&self.device_handler, &mut self.objects);
        match objects.get_mut(&object.render_id()) {
            None => false,
            Some(render_object) => {
                let (vb, ib) = object.vertex_buffer(handler);
                render_object.vertex_buffer = vb;
                render_object.index_buffer = ib;
                true
            }
        }
    }

    /// Synchronizes the information of vertices of `objects` in the CPU memory
    /// and that in the GPU memory.
    ///
    /// If there exists a render object which does not exist in the scene, returns false.
    #[inline]
    pub fn update_vertex_buffers<'a, R, I>(&mut self, objects: I) -> bool
    where
        R: 'a + Rendered,
        I: IntoIterator<Item = &'a R>, {
        let closure = move |flag, object: &R| flag && self.update_vertex_buffer(object);
        objects.into_iter().fold(true, closure)
    }

    /// Synchronizes the information of bind group of `object` in the CPU memory
    /// and that in the GPU memory.
    ///
    /// If there does not exist the render object in the scene, does nothing and returns false.
    #[inline]
    pub fn update_bind_group<R: Rendered>(&mut self, object: &R) -> bool {
        let (handler, objects) = (&self.device_handler, &mut self.objects);
        match objects.get_mut(&object.render_id()) {
            Some(render_object) => {
                let bind_group = object.bind_group(handler, &render_object.bind_group_layout);
                render_object.bind_group = bind_group;
                true
            }
            _ => false,
        }
    }
    /// Synchronizes the information of bind group of `object` in the CPU memory
    /// and that in the GPU memory.
    ///
    /// If there exists a render object which does not exist in the scene, returns false.
    #[inline]
    pub fn update_bind_groups<'a, R, I>(&mut self, objects: I) -> bool
    where
        R: 'a + Rendered,
        I: IntoIterator<Item = &'a R>, {
        let closure = move |flag, object: &R| flag && self.update_bind_group(object);
        objects.into_iter().fold(true, closure)
    }
    /// Synchronizes the information of pipeline of `object` in the CPU memory
    /// and that in the GPU memory.
    ///
    /// If there does not exist the render object in the scene, does nothing and returns false.
    #[inline]
    pub fn update_pipeline<R: Rendered>(&mut self, object: &R) -> bool {
        let (handler, objects) = (&self.device_handler, &mut self.objects);
        match objects.get_mut(&object.render_id()) {
            Some(render_object) => {
                let device = handler.device();
                let pipeline_layout = device.create_pipeline_layout(&PipelineLayoutDescriptor {
                    bind_group_layouts: &[
                        Some(&self.bind_group_layout),
                        Some(&render_object.bind_group_layout),
                    ],
                    ..Default::default()
                });
                render_object.pipeline =
                    object.pipeline(handler, &pipeline_layout, &self.descriptor);
                true
            }
            _ => false,
        }
    }
    /// Synchronizes the information of pipeline of `object` in the CPU memory
    /// and that in the GPU memory.
    ///
    /// If there exists a render object which does not exist in the scene, returns false.
    #[inline]
    pub fn update_pipelines<'a, R, I>(&mut self, objects: I) -> bool
    where
        R: 'a + Rendered,
        I: IntoIterator<Item = &'a R>, {
        let closure = move |flag, object: &R| flag && self.update_pipeline(object);
        objects.into_iter().fold(true, closure)
    }
    #[inline]
    const fn depth_stencil_attachment_descriptor(
        depth_view: &TextureView,
    ) -> RenderPassDepthStencilAttachment<'_> {
        RenderPassDepthStencilAttachment {
            view: depth_view,
            depth_ops: Some(Operations {
                load: LoadOp::Clear(1.0),
                store: StoreOp::Store,
            }),
            stencil_ops: None,
        }
    }

    /// Renders the scene to `view`.
    ///
    /// # Errors
    /// Returns [`GpuError::SingularCameraMatrix`] when the scene camera is singular.
    pub fn render(&self, view: &TextureView) -> Result<(), GpuError> {
        let bind_group = self.scene_bind_group()?;
        let depth_view = self
            .forward_depth
            .as_ref()
            .map(|tex| tex.create_view(&TextureViewDescriptor::default()));
        let sampled_view = self
            .sampling_buffer
            .as_ref()
            .map(|tex| tex.create_view(&TextureViewDescriptor::default()));
        let mut encoder = self
            .device()
            .create_command_encoder(&CommandEncoderDescriptor { label: None });
        {
            let (attachment, resolve_target) = sampled_view
                .as_ref()
                .map_or((view, None), |sampled_view| (sampled_view, Some(view)));
            let mut rpass = encoder.begin_render_pass(&RenderPassDescriptor {
                color_attachments: &[Some(RenderPassColorAttachment {
                    view: attachment,
                    resolve_target,
                    ops: Operations {
                        load: LoadOp::Clear(self.descriptor.studio.background),
                        store: StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: depth_view
                    .as_ref()
                    .map(Self::depth_stencil_attachment_descriptor),
                ..Default::default()
            });
            rpass.set_bind_group(0, &bind_group, &[]);
            for (_, object) in &self.objects {
                if !object.visible {
                    continue;
                }
                rpass.set_pipeline(&object.pipeline);
                rpass.set_bind_group(1, Some(object.bind_group.as_ref()), &[]);
                rpass.set_vertex_buffer(0, object.vertex_buffer.buffer.slice(..));
                if let Some(ref index_buffer) = object.index_buffer {
                    rpass.set_index_buffer(index_buffer.buffer.slice(..), IndexFormat::Uint32);
                    let index_count = index_buffer.size / size_of::<u32>() as u64;
                    let index_count =
                        u32::try_from(index_count).map_err(|_| GpuError::CountTooLarge {
                            what: "index count",
                        })?;
                    rpass.draw_indexed(0..index_count, 0, 0..1);
                } else {
                    let vertex_count = object.vertex_buffer.size / object.vertex_buffer.stride;
                    let vertex_count =
                        u32::try_from(vertex_count).map_err(|_| GpuError::CountTooLarge {
                            what: "vertex count",
                        })?;
                    rpass.draw(0..vertex_count, 0..1);
                }
            }
        }
        self.queue().submit(vec![encoder.finish()]);
        Ok(())
    }

    /// Render image to buffer.
    ///
    /// # Errors
    /// Returns a [`GpuError`] for singular camera transforms, invalid readback sizing, device
    /// polling failures, buffer-map failures, or a closed map-completion channel.
    pub async fn render_to_buffer(&self) -> Result<Vec<u8>, GpuError> {
        let texture = self.compatible_texture();
        let view = texture.create_view(&TextureViewDescriptor::default());
        self.render(&view)?;
        let (device, queue) = (self.device(), self.queue());

        let (width, height) = self.descriptor.render_texture.canvas_size;
        let format = self.descriptor.render_texture.format;
        if format.block_dimensions() != (1, 1) {
            return Err(GpuError::UnsupportedReadbackFormat { format });
        }
        let bytes_per_texel = format
            .block_copy_size(Some(TextureAspect::All))
            .ok_or(GpuError::UnsupportedReadbackFormat { format })?;
        let unpadded_bytes_per_row =
            width
                .checked_mul(bytes_per_texel)
                .ok_or(GpuError::CountTooLarge {
                    what: "render readback row pitch",
                })?;
        let alignment = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let padding = (alignment - unpadded_bytes_per_row % alignment) % alignment;
        let padded_bytes_per_row =
            unpadded_bytes_per_row
                .checked_add(padding)
                .ok_or(GpuError::CountTooLarge {
                    what: "aligned render readback row pitch",
                })?;
        let buffer_size = u64::from(padded_bytes_per_row) * u64::from(height);
        let output_size_u64 = u64::from(unpadded_bytes_per_row) * u64::from(height);
        let output_size =
            usize::try_from(output_size_u64).map_err(|_| GpuError::ReadbackTooLarge {
                bytes: output_size_u64,
            })?;
        let buffer = device.create_buffer(&BufferDescriptor {
            label: None,
            mapped_at_creation: false,
            usage: BufferUsages::COPY_DST | BufferUsages::MAP_READ,
            size: buffer_size,
        });
        let mut encoder = device.create_command_encoder(&CommandEncoderDescriptor { label: None });
        encoder.copy_texture_to_buffer(
            TexelCopyTextureInfo {
                texture: &texture,
                mip_level: 0,
                origin: Origin3d::ZERO,
                aspect: TextureAspect::All,
            },
            TexelCopyBufferInfo {
                buffer: &buffer,
                layout: TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(padded_bytes_per_row),
                    rows_per_image: Some(height),
                },
            },
            Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );
        queue.submit(Some(encoder.finish()));
        let buffer_slice = buffer.slice(..);
        let (sender, receiver) = futures_intrusive::channel::shared::oneshot_channel();
        buffer_slice.map_async(MapMode::Read, move |result| {
            // If the receiver was dropped there is nobody left to observe the map result.
            let _ = sender.send(result);
        });
        let poll_type = PollType::Wait {
            submission_index: None,
            timeout: None,
        };
        device.poll(poll_type)?;
        match receiver.receive().await {
            Some(Ok(())) => {
                let mapped = buffer_slice.get_mapped_range();
                if padded_bytes_per_row == unpadded_bytes_per_row {
                    return Ok(mapped.iter().copied().collect());
                }
                let row_bytes = usize::try_from(unpadded_bytes_per_row).map_err(|_| {
                    GpuError::ReadbackTooLarge {
                        bytes: output_size_u64,
                    }
                })?;
                let padded_row_bytes = usize::try_from(padded_bytes_per_row).map_err(|_| {
                    GpuError::ReadbackTooLarge {
                        bytes: output_size_u64,
                    }
                })?;
                let rows = usize::try_from(height).map_err(|_| GpuError::ReadbackTooLarge {
                    bytes: output_size_u64,
                })?;
                let mut output = Vec::with_capacity(output_size);
                for row in mapped.chunks_exact(padded_row_bytes).take(rows) {
                    output.extend_from_slice(&row[..row_bytes]);
                }
                Ok(output)
            }
            Some(Err(error)) => Err(error.into()),
            None => Err(GpuError::BufferMapChannelClosed),
        }
    }
}

impl WindowScene {
    /// Initialize scene compatible with `window`.
    ///
    /// # Errors
    /// Returns a [`GpuError`] if adapter/device/surface initialization fails.
    pub async fn from_window(
        window: Arc<Window>,
        scene_desc: &WindowSceneDescriptor,
    ) -> Result<Self, GpuError> {
        let size = window.inner_size();
        let (device_handler, window_handler) = init_default_device(Some(window)).await?;
        let Some(window_handler) = window_handler else {
            return Err(GpuError::MissingSurface);
        };
        let (device, surface) = (&device_handler.device, &window_handler.surface);
        let render_texture = RenderTextureConfig {
            canvas_size: size.into(),
            format: TextureFormat::Bgra8Unorm,
        };
        let config = render_texture.compatible_surface_config();
        surface.configure(device, &config);

        Ok(Self {
            scene: Scene::new(
                device_handler,
                &SceneDescriptor {
                    studio: scene_desc.studio.clone(),
                    backend_buffer: scene_desc.backend_buffer,
                    render_texture,
                },
            ),
            window_handler,
        })
    }
    /// Get the reference of initializing window.
    #[inline]
    #[must_use]
    pub const fn window(&self) -> &Arc<Window> { &self.window_handler.window }
    /// Get the reference of surface.
    #[inline]
    #[must_use]
    pub const fn surface(&self) -> &Arc<Surface<'_>> { &self.window_handler.surface }
    /// Adjusts the size of the backend buffers (depth or sampling buffer) to the size of the window.
    pub fn size_alignment(&mut self) {
        let size = self.window().inner_size();
        let canvas_size = self.scene.descriptor.render_texture.canvas_size;
        if canvas_size != (size.width, size.height) {
            let mut desc = self.scene.descriptor_mut();
            desc.render_texture.canvas_size = size.into();
            let config = desc.render_texture.compatible_surface_config();
            drop(desc);
            self.surface().configure(self.device(), &config);
        }
    }
    /// Render scene to initializing window.
    ///
    /// # Errors
    /// Returns a [`GpuError`] when rendering fails, including a singular camera transform.
    pub fn render_frame(&mut self) -> Result<(), GpuError> {
        self.size_alignment();
        let surface = self.surface();
        // wgpu 29 replaced the `Result<SurfaceTexture, SurfaceError>` return of
        // `get_current_texture` with the `CurrentSurfaceTexture` enum. A
        // `Suboptimal` texture is still usable; any other non-`Success` status
        // means we reconfigure the surface and retry once.
        let usable = |attempt: CurrentSurfaceTexture| match attempt {
            CurrentSurfaceTexture::Success(texture)
            | CurrentSurfaceTexture::Suboptimal(texture) => Some(texture),
            _ => None,
        };
        let surface_texture = if let Some(texture) = usable(surface.get_current_texture()) {
            texture
        } else {
            let config = self
                .scene
                .descriptor
                .render_texture
                .compatible_surface_config();
            surface.configure(self.device(), &config);
            let Some(texture) = usable(surface.get_current_texture()) else {
                return Ok(());
            };
            texture
        };
        let view = surface_texture
            .texture
            .create_view(&TextureViewDescriptor::default());
        self.render(&view)?;
        surface_texture.present();
        Ok(())
    }
}
