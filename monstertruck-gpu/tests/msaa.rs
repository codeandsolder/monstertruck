// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each test target is its own crate.
#![recursion_limit = "256"]
mod common;
use common::Plane;
use monstertruck_gpu::{BackendBufferConfig, RenderTextureConfig, Scene, SceneDescriptor};
use wgpu::{Backends, TextureFormat};

pub const PICTURE_WIDTH: u32 = 512;
pub const PICTURE_HEIGHT: u32 = 512;

fn exec_msaa_test(backend: Backends, out_dir: &str) -> common::TestResult {
    std::fs::create_dir_all(out_dir)?;
    let handler = common::init_device(backend)?;
    let mut scene = Scene::new(
        handler,
        &SceneDescriptor {
            backend_buffer: BackendBufferConfig {
                sample_count: 1,
                ..BackendBufferConfig::default()
            },
            render_texture: RenderTextureConfig {
                canvas_size: (PICTURE_WIDTH, PICTURE_HEIGHT),
                format: TextureFormat::Rgba8Unorm,
            },
            ..SceneDescriptor::default()
        },
    );
    let plane = new_plane!("shaders/trapezoid.wgsl", "vs_main", "fs_main");
    let buffer0 = common::render_one(&mut scene, &plane)?;
    common::save_rgba8(
        format!("{out_dir}sample_count_one.png"),
        &buffer0,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    scene.descriptor_mut().backend_buffer.sample_count = 4;
    let buffer1 = common::render_one(&mut scene, &plane)?;
    common::save_rgba8(
        format!("{out_dir}sample_count_two.png"),
        &buffer1,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    assert!(!common::same_buffer(&buffer0, &buffer1));
    Ok(())
}

#[test]
fn msaa_test() -> common::TestResult { common::os_alt_exec_test(exec_msaa_test) }
