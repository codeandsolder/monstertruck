// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each test target is its own crate.
#![recursion_limit = "256"]
mod common;

use monstertruck_gpu::*;
use monstertruck_render::*;
use wgpu::*;

fn render_wireframe_without_depth_format(
    backend: Backends,
    format: TextureFormat,
    bytes_per_texel: usize,
) -> common::TestResult {
    let instance = wgpu::Instance::new(InstanceDescriptor {
        backends: backend,
        ..InstanceDescriptor::new_without_display_handle()
    });
    let handler = common::init_device(&instance)?;
    let mut scene = Scene::new(
        handler,
        &SceneDescriptor {
            backend_buffer: BackendBufferConfig {
                depth_test: false,
                ..Default::default()
            },
            render_texture: RenderTextureConfig {
                canvas_size: (32, 32),
                format,
            },
            ..Default::default()
        },
    );
    let lines = vec![(Point3::new(-0.5, 0.0, 0.0), Point3::new(0.5, 0.0, 0.0))];
    let wireframe: WireFrameInstance = scene
        .instance_creator()
        .create_instance(&lines, &WireFrameState::default())?;

    let pixels = common::render_one(&mut scene, &wireframe)?;
    assert_eq!(pixels.len(), 32 * 32 * bytes_per_texel);
    Ok(())
}

fn render_wireframe_without_depth(backend: Backends, _: &str) -> common::TestResult {
    render_wireframe_without_depth_format(backend, TextureFormat::Rgba8Unorm, 4)
}

fn render_wireframe_rgba16(backend: Backends, _: &str) -> common::TestResult {
    render_wireframe_without_depth_format(backend, TextureFormat::Rgba16Float, 8)
}

#[test]
fn wireframe_pipeline_matches_depth_disabled_scene() -> common::TestResult {
    common::os_alt_exec_test(render_wireframe_without_depth)
}

#[test]
fn render_readback_uses_configured_texture_format_size() -> common::TestResult {
    common::os_alt_exec_test(render_wireframe_rgba16)
}
