// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each example and test target is its own crate.
#![recursion_limit = "256"]
//! Material Samples
//!
//! - The more right the model, the higher the reflectance.
//! - The upper the model, the higher the roughness.
//!
//! The most right and lowest model is black because it does not diffuse light
//! and no roughness in microfacet.

mod app;
use app::{App, AppResult, async_trait};
use monstertruck_gpu::{
    BackendBufferConfig, Camera, Light, LightType, ProjectionMethod, StudioConfig, WindowScene,
    WindowSceneDescriptor, wgpu,
};
use monstertruck_meshing::prelude::*;
use monstertruck_modeling::{Solid, builder};
use monstertruck_render::{CreatorCreator, Material, PolygonInstance, PolygonState};
use std::f64::consts::PI;
use std::sync::Arc;

const N: u32 = 5;
const BACKGROUND: [f64; 4] = [45.0 / 255.0, 36.0 / 255.0, 42.0 / 255.0, 1.0];
const BOXCOLOR: [f64; 4] = [208.0 / 255.0, 176.0 / 255.0, 107.0 / 255.0, 1.0];

struct MyApp {
    scene: WindowScene,
    instances: Vec<PolygonInstance>,
    matrices: Vec<Matrix4>,
}

#[async_trait(?Send)]
impl App for MyApp {
    async fn init(window: Arc<winit::window::Window>) -> AppResult<Self> {
        let side_length = f64::from(N + 1) * 1.5;
        let camera_dist = side_length / 2.0 / (PI / 8.0).tan();
        let half_side = side_length / 2.0;
        let half_camera_dist = camera_dist / 2.0;
        let sample_count = 4;
        let scene_desc = WindowSceneDescriptor {
            studio: StudioConfig {
                camera: Camera {
                    matrix: Matrix4::from_translation(camera_dist * Vector3::unit_z()),
                    method: ProjectionMethod::perspective(Rad(PI / 4.0)),
                    near_clip: 0.1,
                    far_clip: 100.0,
                },
                lights: vec![
                    Light {
                        position: Point3::new(-half_side, -half_side, half_camera_dist),
                        color: Vector3::new(0.5, 0.5, 0.5),
                        kind: LightType::Point,
                    },
                    Light {
                        position: Point3::new(-half_side, half_side, half_camera_dist),
                        color: Vector3::new(0.5, 0.5, 0.5),
                        kind: LightType::Point,
                    },
                    Light {
                        position: Point3::new(half_side, -half_side, half_camera_dist),
                        color: Vector3::new(0.5, 0.5, 0.5),
                        kind: LightType::Point,
                    },
                    Light {
                        position: Point3::new(half_side, half_side, half_camera_dist),
                        color: Vector3::new(0.5, 0.5, 0.5),
                        kind: LightType::Point,
                    },
                ],
                background: wgpu::Color {
                    r: BACKGROUND[0],
                    g: BACKGROUND[1],
                    b: BACKGROUND[2],
                    a: BACKGROUND[3],
                },
            },
            backend_buffer: BackendBufferConfig {
                sample_count,
                ..Default::default()
            },
        };
        let mut scene = WindowScene::from_window(window, &scene_desc).await?;
        let vertex = builder::vertex(Point3::new(-0.5, -0.5, -0.5));
        let edge = builder::extrude(&vertex, Vector3::unit_x());
        let face = builder::extrude(&edge, Vector3::unit_y());
        let cube: Solid = builder::extrude(&face, Vector3::unit_z());
        let mesh = cube.triangulation(0.01).to_polygon();
        let instance: PolygonInstance = scene
            .instance_creator()
            .create_instance(&mesh, &PolygonState::default())?;
        let mut matrices = Vec::new();
        let instances: Vec<_> = (0..N)
            .flat_map(move |row| (0..N).map(move |column| (row, column)))
            .map(|(row, column)| {
                let mut instance = instance.clone_instance();
                let (reflectance, roughness) = (
                    f64::from(row) / f64::from(N - 1),
                    f64::from(column) / f64::from(N - 1),
                );
                let matrix = Matrix4::from_translation(Vector3::new(
                    1.5_f64.mul_add(f64::from(row + 1), -(side_length / 2.0)),
                    1.5_f64.mul_add(f64::from(column + 1), -(side_length / 2.0)),
                    0.0,
                ));
                matrices.push(matrix);
                *instance.instance_state_mut() = PolygonState {
                    matrix,
                    material: Material {
                        albedo: Vector4::from(BOXCOLOR),
                        reflectance,
                        roughness,
                        ambient_ratio: 0.02,
                        background_ratio: 0.0,
                        alpha_blend: false,
                    },
                    ..Default::default()
                };
                instance
            })
            .collect();
        for shape in &instances {
            scene.add_object(shape);
        }
        Ok(Self {
            scene,
            instances,
            matrices,
        })
    }
    fn render(&mut self) {
        let time = self.scene.elapsed().as_secs_f64();
        for (index, shape) in self.instances.iter_mut().enumerate() {
            let pair_sign = if (index / 2) % 2 == 0 { 1.0 } else { -1.0 };
            let axis = if index % 2 == 0 {
                pair_sign * Vector3::unit_y()
            } else {
                -pair_sign * Vector3::unit_x()
            };
            shape.instance_state_mut().matrix =
                self.matrices[index] * Matrix4::from_axis_angle(axis, Rad(time * PI / 2.0));
            self.scene.update_bind_group(shape);
        }
        if let Err(error) = self.scene.render_frame() {
            eprintln!("render failed: {error}");
        }
    }
}

fn main() { MyApp::run() }
