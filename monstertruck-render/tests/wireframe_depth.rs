// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each test target is its own crate.
#![recursion_limit = "256"]
mod common;

use monstertruck_gpu::*;
use monstertruck_render::*;
use wgpu::*;

fn render_wireframe_without_depth(backend: Backends, _: &str) -> common::TestResult {
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
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let lines = vec![(Point3::new(-0.5, 0.0, 0.0), Point3::new(0.5, 0.0, 0.0))];
    let wireframe: WireFrameInstance = scene
        .instance_creator()
        .create_instance(&lines, &WireFrameState::default())?;

    let pixels = common::render_one(&mut scene, &wireframe)?;
    assert_eq!(pixels.len(), 32 * 32 * 4);
    Ok(())
}

#[test]
fn wireframe_pipeline_matches_depth_disabled_scene() -> common::TestResult {
    common::os_alt_exec_test(render_wireframe_without_depth)
}
