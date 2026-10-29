// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each test target is its own crate.
#![recursion_limit = "256"]
mod common;
use common::Plane;
use monstertruck_core::cgmath64::{Matrix4, Point3, Rad, Vector3, Vector4};
use monstertruck_gpu::{
    Camera, Light, LightType, ProjectionMethod, RenderTextureConfig, Scene, SceneDescriptor,
    StudioConfig,
};
use wgpu::{Backends, Color, TextureFormat};

const PICTURE_WIDTH: u32 = 256;
const PICTURE_HEIGHT: u32 = 256;

const CAMERA_MATRIX: Matrix4 = Matrix4::from_cols(
    Vector4::new(1.0, 2.1, 3.2, 4.3),
    Vector4::new(5.4, 6.5, 7.6, 8.7),
    Vector4::new(9.8, 10.9, 11.0, 12.0),
    Vector4::new(13.0, 14.0, 15.0, 16.23),
);
const CAMERA_FOV: Rad<f64> = Rad(std::f64::consts::PI / 4.0);
const CAMERA_NEARCLIP: f64 = 0.1;
const CAMERA_FARCLIP: f64 = 10.0;

const POINT_LIGHT: Light = Light {
    position: Point3::new(0.1, 0.2, 0.3),
    color: Vector3::new(0.4, 0.5, 0.6),
    kind: LightType::Point,
};
const UNIFORM_LIGHT: Light = Light {
    position: Point3::new(1.1, 1.2, 1.3),
    color: Vector3::new(1.4, 1.5, 1.6),
    kind: LightType::Uniform,
};

fn exec_bind_group_test(backend: Backends, out_dir: &str) -> common::TestResult {
    std::fs::create_dir_all(out_dir)?;
    let camera = Camera {
        matrix: CAMERA_MATRIX,
        method: ProjectionMethod::perspective(CAMERA_FOV),
        near_clip: CAMERA_NEARCLIP,
        far_clip: CAMERA_FARCLIP,
    };
    println!(
        "camera projection:\n{:?}",
        camera.projection(f64::from(PICTURE_WIDTH) / f64::from(PICTURE_HEIGHT))
    );
    let desc = SceneDescriptor {
        studio: StudioConfig {
            camera,
            lights: vec![POINT_LIGHT, UNIFORM_LIGHT],
            background: Color {
                r: 0.1,
                g: 0.2,
                b: 0.3,
                a: 0.4,
            },
        },
        render_texture: RenderTextureConfig {
            canvas_size: (PICTURE_WIDTH, PICTURE_HEIGHT),
            format: TextureFormat::Rgba8Unorm,
        },
        ..SceneDescriptor::default()
    };
    let handler = common::init_device(backend)?;
    let mut scene = Scene::new(handler, &desc);
    let plane = new_plane!("shaders/unicolor.wgsl", "vs_main", "fs_main");
    let buffer0 = common::render_one(&mut scene, &plane)?;
    let plane = new_plane!("shaders/bindgroup.wgsl", "vs_main", "fs_main");
    let buffer1 = common::render_one(&mut scene, &plane)?;
    let plane = new_plane!("shaders/bindgroup.wgsl", "vs_main", "fs_main_anti");
    let buffer2 = common::render_one(&mut scene, &plane)?;
    common::save_rgba8(
        format!("{out_dir}unicolor.png"),
        &buffer0,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    common::save_rgba8(
        format!("{out_dir}bindgroup.png"),
        &buffer1,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    common::save_rgba8(
        format!("{out_dir}anti-bindgroup.png"),
        &buffer2,
        PICTURE_WIDTH,
        PICTURE_HEIGHT,
    )?;
    assert!(common::same_buffer(&buffer0, &buffer1));
    assert!(!common::same_buffer(&buffer0, &buffer2));
    Ok(())
}

#[test]
fn bind_group_test() -> common::TestResult { common::os_alt_exec_test(exec_bind_group_test) }
