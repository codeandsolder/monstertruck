// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each test target is its own crate.
#![recursion_limit = "256"]
mod common;
use common::Plane;
use monstertruck_gpu::{RenderId, RenderTextureConfig, Scene, SceneDescriptor};
use wgpu::{Backends, TextureFormat};

const PICTURE_WIDTH: u32 = 256;
const PICTURE_HEIGHT: u32 = 256;

fn exec_math_util_test(backend: Backends, out_dir: &str) -> common::TestResult {
    std::fs::create_dir_all(out_dir)?;
    let desc = SceneDescriptor {
        render_texture: RenderTextureConfig {
            canvas_size: (PICTURE_WIDTH, PICTURE_HEIGHT),
            format: TextureFormat::Rgba8Unorm,
        },
        ..SceneDescriptor::default()
    };
    let mut scene = Scene::new(common::init_device(backend)?, &desc);
    let plane = new_plane!("shaders/unicolor.wgsl", "vs_main", "fs_main");
    let buffer0 = common::render_one(&mut scene, &plane)?;
    let shader = include_str!("../wgsl-utils/math.wgsl").to_string()
        + include_str!("shaders/math-util.wgsl");
    let plane = Plane {
        shader: &shader,
        vs_entpt: "vs_main",
        fs_entpt: "fs_main",
        id: RenderId::generate(),
    };
    let buffer1 = common::render_one(&mut scene, &plane)?;
    let plane = Plane {
        shader: &shader,
        vs_entpt: "vs_main",
        fs_entpt: "fs_main_anti",
        id: RenderId::generate(),
    };
    let buffer2 = common::render_one(&mut scene, &plane)?;
    common::save_rgba8(
        format!("{out_dir}unicolor.png"),
        &buffer0,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    common::save_rgba8(
        format!("{out_dir}math-util.png"),
        &buffer1,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    common::save_rgba8(
        format!("{out_dir}anti-math-util.png"),
        &buffer2,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    assert!(common::same_buffer(&buffer0, &buffer1));
    assert!(!common::same_buffer(&buffer0, &buffer2));
    Ok(())
}

#[test]
fn math_util_test() -> common::TestResult { common::os_alt_exec_test(exec_math_util_test) }
