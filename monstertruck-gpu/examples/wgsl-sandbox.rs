// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each example and test target is its own crate.
#![recursion_limit = "256"]
//! A sample of creating a render object by implementing "Rendered" in a new structure.
//!
//! One can use xyr WGSL shader in the following way:
//!
//! - Enter the shader path as an argument when executing the program.
//! - Drag and drop the shader into the window.
//!
//! The rule of shaders:
//!
//! - One can draw a image by implementing the function:
//!
//! ```wgsl
//! vec4<f32> main_image(coord: vec2<f32>, env: Environment);
//! ```
//!
//! - The parameter `coord` is the fragment coordinate. The origin is the lower left.
//! - The parameter `env` has the environment information. The declaration of struct is the following:
//!
//! ```wgsl
//! struct Environment {
//!     resolution: vec2<f32>;  // the resolution of the image
//!     mouse: vec4<f32>;       // the mouse information behaving the same as `iMouse` in Shadertoy.
//!     time: f32;              // the number of seconds since the application started.
//! };
//! ```
//!
//! Also, see the sample `newton-cuberoot.wgsl`, default shader, in `examples`.
#![allow(
    deprecated,
    reason = "example intentionally demonstrates the legacy event-loop integration while winit compatibility is maintained"
)]

use monstertruck_gpu::{
    BufferHandler, DeviceHandler, RenderId, Rendered, SceneDescriptor, WindowScene,
    WindowSceneDescriptor, bind_group_util, impl_render_id, wgpu,
};
use std::sync::Arc;
use wgpu::{
    BindGroup, BindGroupLayout, BindGroupLayoutDescriptor, BindGroupLayoutEntry, BindingType,
    BlendState, BufferBindingType, BufferUsages, ColorTargetState, ColorWrites, DepthBiasState,
    DepthStencilState, Device, Face, FragmentState, FrontFace, MultisampleState,
    PipelineCompilationOptions, PipelineLayout, PolygonMode, PrimitiveState, PrimitiveTopology,
    RenderPipeline, RenderPipelineDescriptor, ShaderModule, ShaderModuleDescriptor, ShaderSource,
    ShaderStages, StencilState, TextureFormat, VertexAttribute, VertexBufferLayout, VertexFormat,
    VertexState, VertexStepMode,
};
use winit::event::{ElementState, Event, StartCause, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};

const DEFAULT_SHADER: &str = include_str!("newton-cuberoot.wgsl");

/// minimum example for implementing `Rendered`.
mod plane {
    use super::{
        Arc, BindGroup, BindGroupLayout, BindGroupLayoutDescriptor, BindGroupLayoutEntry,
        BindingType, BlendState, BufferBindingType, BufferHandler, BufferUsages, ColorTargetState,
        ColorWrites, DepthBiasState, DepthStencilState, Device, DeviceHandler, Face, FragmentState,
        FrontFace, MultisampleState, PipelineCompilationOptions, PipelineLayout, PolygonMode,
        PrimitiveState, PrimitiveTopology, RenderId, RenderPipeline, RenderPipelineDescriptor,
        Rendered, SceneDescriptor, ShaderModule, ShaderModuleDescriptor, ShaderSource,
        ShaderStages, StencilState, TextureFormat, VertexAttribute, VertexBufferLayout,
        VertexFormat, VertexState, VertexStepMode, bind_group_util, impl_render_id, wgpu,
    };

    /// Canvas to draw by fragment shader.
    pub struct Plane {
        module: ShaderModule,
        pub mouse: [f32; 4],
        id: RenderId,
    }

    const BASE_PREFIX: &str = "struct SceneInfo {
    background_color: vec4<f32>,
    resolution: vec2<u32>,
    time: f32,
    nlights: u32,
}

struct Mouse {
    mouse: vec4<f32>,
}

@group(0)
@binding(2)
var<uniform> info__: SceneInfo;

@group(1)
@binding(0)
var<uniform> mouse__: Mouse;

struct Environment {
    resolution: vec2<f32>,
    mouse: vec4<f32>,
    time: f32,
}

";

    const BASE_SHADER: &str = "@vertex
fn vs_main(@location(0) idx: u32) -> @builtin(position) vec4<f32> {
    var vertex: array<vec2<f32>, 4>;
    vertex[0] = vec2<f32>(-1.0, -1.0);
    vertex[1] = vec2<f32>(1.0, -1.0);
    vertex[2] = vec2<f32>(-1.0, 1.0);
    vertex[3] = vec2<f32>(1.0, 1.0);
    return vec4<f32>(vertex[idx], 0.0, 1.0);
}

@fragment
fn fs_main(@builtin(position) position: vec4<f32>) -> @location(0) vec4<f32> {
    var env: Environment;
    env.resolution = vec2<f32>(info__.resolution);
    env.mouse = mouse__.mouse;
    env.time = info__.time;
    let coord = vec2<f32>(
        position.x,
        env.resolution.y - position.y,
    );
    return main_image(coord, env);
}
";

    // Implementation of Rendered for Plane.
    impl Rendered for Plane {
        // `Rendered::render_id()` can be implemented by macro.
        impl_render_id!(id);

        // Vertices: [0, 1, 2, 2, 1, 3] as [u32; 6].
        // There is not the index buffer.
        fn vertex_buffer(
            &self,
            handler: &DeviceHandler,
        ) -> (Arc<BufferHandler>, Option<Arc<BufferHandler>>) {
            let vertex_buffer = BufferHandler::from_slice(
                &[0, 1, 2, 2, 1, 3],
                handler.device(),
                BufferUsages::VERTEX,
            );
            (Arc::new(vertex_buffer), None)
        }

        // bind group is only one uniform buffer: resolution
        fn bind_group_layout(&self, handler: &DeviceHandler) -> Arc<BindGroupLayout> {
            Arc::new(
                handler
                    .device()
                    .create_bind_group_layout(&BindGroupLayoutDescriptor {
                        label: None,
                        entries: &[BindGroupLayoutEntry {
                            binding: 0,
                            visibility: ShaderStages::FRAGMENT,
                            ty: BindingType::Buffer {
                                ty: BufferBindingType::Uniform,
                                has_dynamic_offset: false,
                                min_binding_size: None,
                            },
                            count: None,
                        }],
                    }),
            )
        }
        // bind group is only one uniform buffer: resolution
        fn bind_group(&self, handler: &DeviceHandler, layout: &BindGroupLayout) -> Arc<BindGroup> {
            Arc::new(bind_group_util::create_bind_group(
                handler.device(),
                layout,
                vec![
                    BufferHandler::from_slice(&self.mouse, handler.device(), BufferUsages::UNIFORM)
                        .binding_resource(),
                ],
            ))
        }

        // Describe pipeline
        fn pipeline(
            &self,
            handler: &DeviceHandler,
            layout: &PipelineLayout,
            scene_desc: &SceneDescriptor,
        ) -> Arc<RenderPipeline> {
            let SceneDescriptor {
                backend_buffer,
                render_texture,
                ..
            } = scene_desc;
            Arc::new(
                handler
                    .device()
                    .create_render_pipeline(&RenderPipelineDescriptor {
                        layout: Some(layout),
                        vertex: VertexState {
                            module: &self.module,
                            entry_point: Some("vs_main"),
                            buffers: &[VertexBufferLayout {
                                array_stride: 4,
                                step_mode: VertexStepMode::Vertex,
                                attributes: &[VertexAttribute {
                                    format: VertexFormat::Uint32,
                                    offset: 0,
                                    shader_location: 0,
                                }],
                            }],
                            compilation_options: PipelineCompilationOptions::default(),
                        },
                        fragment: Some(FragmentState {
                            module: &self.module,
                            entry_point: Some("fs_main"),
                            targets: &[Some(ColorTargetState {
                                format: render_texture.format,
                                blend: Some(BlendState::REPLACE),
                                write_mask: ColorWrites::ALL,
                            })],
                            compilation_options: PipelineCompilationOptions::default(),
                        }),
                        primitive: PrimitiveState {
                            topology: PrimitiveTopology::TriangleList,
                            front_face: FrontFace::Ccw,
                            cull_mode: Some(Face::Back),
                            polygon_mode: PolygonMode::Fill,
                            ..PrimitiveState::default()
                        },
                        depth_stencil: Some(DepthStencilState {
                            format: TextureFormat::Depth32Float,
                            depth_write_enabled: Some(true),
                            depth_compare: Some(wgpu::CompareFunction::Less),
                            stencil: StencilState::default(),
                            bias: DepthBiasState::default(),
                        }),
                        multisample: MultisampleState {
                            count: backend_buffer.sample_count,
                            mask: !0,
                            alpha_to_coverage_enabled: false,
                        },
                        label: None,
                        multiview_mask: None,
                        cache: None,
                    }),
            )
        }
    }

    impl Plane {
        /// constructor
        /// # Arguments
        /// - device: Device, provided by wgpu.
        /// - shader: the inputted fragment shader
        pub fn new(device: &Device, shader: &str) -> Option<Self> {
            Some(Self {
                module: create_module(device, shader)?,
                mouse: [0.0; 4],
                id: RenderId::generate(),
            })
        }

        pub fn set_shader(&mut self, device: &Device, shader: &str) {
            if let Some(module) = create_module(device, shader) {
                self.module = module;
            }
        }
    }

    fn create_module(device: &Device, shader: &str) -> Option<ShaderModule> {
        use ::naga::{
            front::wgsl,
            valid::{Capabilities, ValidationFlags, Validator},
        };
        let mut source = BASE_PREFIX.to_string();
        source += shader;
        source += BASE_SHADER;

        Validator::new(ValidationFlags::all(), Capabilities::empty())
            .validate(
                &wgsl::parse_str(&source)
                    .map_err(|error| println!("WGSL Parse Error: {error}"))
                    .ok()?,
            )
            .map_err(|error| println!("WGSL Validation Error: {error}"))
            .ok()?;

        Some(device.create_shader_module(ShaderModuleDescriptor {
            source: ShaderSource::Wgsl(source.into()),
            label: None,
        }))
    }
}
use plane::Plane;

#[allow(
    clippy::future_not_send,
    reason = "winit event loops and windows are deliberately thread-affine"
)]
async fn run(
    event_loop: EventLoop<()>,
    window: winit::window::Window,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut scene =
        WindowScene::from_window(Arc::new(window), &WindowSceneDescriptor::default()).await?;
    let args: Vec<_> = std::env::args().collect();
    let source = if args.len() > 1 {
        match std::fs::read_to_string(&args[1]) {
            Ok(code) => code,
            Err(error) => {
                println!("{error:?}");
                DEFAULT_SHADER.to_string()
            }
        }
    } else {
        DEFAULT_SHADER.to_string()
    };
    let Some(mut plane) = Plane::new(scene.device(), &source).or_else(|| {
        eprintln!("input shader is invalid; falling back to the built-in shader");
        Plane::new(scene.device(), DEFAULT_SHADER)
    }) else {
        return Err(std::io::Error::other("built-in WGSL shader is invalid").into());
    };
    // Adds a plane to the scene!
    scene.add_object(&plane);

    let mut dragging = false;
    let mut clicked = false;
    let mut cursor = [0.0; 2];
    let routine = move |ev: Event<()>, target: &ActiveEventLoop| {
        let control_flow = match ev {
            Event::NewEvents(StartCause::Poll) => {
                scene.window().request_redraw();
                ControlFlow::Poll
            }
            Event::WindowEvent { event, .. } => match event {
                WindowEvent::RedrawRequested => {
                    scene.update_bind_group(&plane);
                    if clicked {
                        plane.mouse[3] = -plane.mouse[3];
                        clicked = false;
                    }
                    if let Err(error) = scene.render_frame() {
                        eprintln!("render failed: {error}");
                    }
                    ControlFlow::Poll
                }
                WindowEvent::CloseRequested => {
                    target.exit();
                    ControlFlow::Poll
                }
                WindowEvent::DroppedFile(path) => {
                    match std::fs::read_to_string(path) {
                        Ok(code) => {
                            plane.set_shader(scene.device(), &code);
                            scene.update_pipeline(&plane);
                        }
                        Err(error) => println!("{error:?}"),
                    }
                    ControlFlow::Poll
                }
                WindowEvent::MouseInput { state, .. } => {
                    dragging = state == ElementState::Pressed;
                    clicked = dragging;
                    if dragging {
                        plane.mouse[0] = cursor[0];
                        plane.mouse[1] = cursor[1];
                        plane.mouse[2] = cursor[0];
                        plane.mouse[3] = cursor[1];
                    } else {
                        plane.mouse[2] = -plane.mouse[2];
                    }
                    ControlFlow::Poll
                }
                WindowEvent::CursorMoved { position, .. } => {
                    let height = f64::from(scene.descriptor().render_texture.canvas_size.1);
                    cursor = [
                        shader_coordinate(position.x),
                        shader_coordinate(height - position.y),
                    ];
                    if dragging {
                        plane.mouse[0] = cursor[0];
                        plane.mouse[1] = cursor[1];
                    }
                    ControlFlow::Poll
                }
                _ => ControlFlow::Poll,
            },
            _ => ControlFlow::Poll,
        };
        target.set_control_flow(control_flow);
    };
    #[cfg(not(target_arch = "wasm32"))]
    event_loop.run(routine)?;
    #[cfg(target_arch = "wasm32")]
    {
        use winit::platform::web::EventLoopExtWebSys;
        event_loop.spawn(routine);
    }
    Ok(())
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "winit f64 cursor coordinates are intentionally narrowed to f32 shader uniforms"
)]
const fn shader_coordinate(value: f64) -> f32 { value as f32 }

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let event_loop = winit::event_loop::EventLoop::new()?;
    let mut wa = winit::window::Window::default_attributes();
    wa.title = "WGSL Sandbox".to_string();
    let window = event_loop.create_window(wa)?;
    #[cfg(not(target_arch = "wasm32"))]
    pollster::block_on(run(event_loop, window))?;
    #[cfg(target_arch = "wasm32")]
    {
        std::panic::set_hook(Box::new(console_error_panic_hook::hook));
        console_log::init()?;
        use winit::platform::web::WindowExtWebSys;
        // On wasm, append the canvas to the document body
        let appended = web_sys::window()
            .and_then(|win| win.document())
            .and_then(|doc| doc.body())
            .and_then(|body| {
                let canvas = window.canvas()?;
                canvas.set_width(500);
                canvas.set_height(500);
                body.append_child(&web_sys::Element::from(canvas)).ok()
            });
        if appended.is_none() {
            return Err(std::io::Error::other("could not append canvas to document body").into());
        }
        wasm_bindgen_futures::spawn_local(async move {
            if let Err(error) = run(event_loop, window).await {
                eprintln!("WGSL sandbox failed: {error}");
            }
        });
    }
    Ok(())
}
