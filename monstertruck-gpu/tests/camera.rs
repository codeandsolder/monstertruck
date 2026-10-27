// See the note in src/lib.rs: wgpu overflows the default auto-trait recursion
// limit, and each test target is its own crate.
#![recursion_limit = "256"]
#![cfg(not(target_arch = "wasm32"))]

use monstertruck_core::{bounding_box::*, cgmath64::*, prop_assert_near, tolerance::*};
use monstertruck_gpu::*;
use proptest::prelude::*;
use std::f64::consts::PI;

proptest! {
    #[test]
    fn parallel_view_fitting(
        rot_axis in (-1.0f64..=1.0, 0.0..2.0 * PI),
        rot_angle in 0.0..2.0 * PI,
        aspect in 0.5f64..=10.0f64,
        near_clip in 0.01f64..=100.0,
        points in prop::array::uniform32(prop::array::uniform3(-100.0f64..=100.0f64)),
    ) {
        let z = rot_axis.0;
        let r = (1.0 - z * z).sqrt();
        let rot_axis = Vector3::new(r * f64::cos(rot_axis.1), r * f64::sin(rot_axis.1), z);
        let rot = Matrix3::from_axis_angle(rot_axis, Rad(rot_angle));
        let points: [Point3; 32] = points.map(|p| p.into());

        let same_plane = points.windows(4).all(|p| {
            Matrix3::from_cols(p[1] - p[0], p[2] - p[0], p[3] - p[0])
                .determinant()
                .so_small()
        });
        prop_assume!(!same_plane);

        let camera = Camera::parallel_view_fitting(rot, aspect, near_clip, &points).expect("rotation is invertible");
        let nvolume = points
            .iter()
            .map(|p| camera.projection(aspect).transform_point(*p))
            .collect::<BoundingBox<Point3>>();

        let (min, max) = (nvolume.min(), nvolume.max());
        prop_assert_near!(min.x.min(min.y), -1.0);
        prop_assert_near!(max.x.max(max.y), 1.0);
        prop_assert_near!(min.z, 0.0);
        prop_assert_near!(max.z, 1.0);
    }

    #[test]
    fn perspective_view_fitting(
        rot_axis in (-1.0f64..=1.0, 0.0..2.0 * PI),
        rot_angle in 0.0..2.0 * PI,
        aspect in 0.5f64..=10.0f64,
        fov in 0.1f64..=3.05,
        points in prop::array::uniform32(prop::array::uniform3(-100.0f64..=100.0f64)),
    ) {
        let z = rot_axis.0;
        let r = (1.0 - z * z).sqrt();
        let rot_axis = Vector3::new(r * f64::cos(rot_axis.1), r * f64::sin(rot_axis.1), z);
        let rot = Matrix3::from_axis_angle(rot_axis, Rad(rot_angle));
        let points: [Point3; 32] = points.map(|p| p.into());

        let same_plane = points.windows(4).all(|p| {
            Matrix3::from_cols(p[1] - p[0], p[2] - p[0], p[3] - p[0])
                .determinant()
                .so_small()
        });
        prop_assume!(!same_plane);

        let camera = Camera::perspective_view_fitting(rot, aspect, Rad(fov), &points).expect("rotation is invertible");
        prop_assume!(camera.near_clip > TOLERANCE);

        let nvolume = points
            .iter()
            .map(|p| camera.projection(aspect).transform_point(*p))
            .collect::<BoundingBox<Point3>>();

        let (min, max) = (nvolume.min(), nvolume.max());
        prop_assert_near!(min.x.min(min.y), -1.0);
        prop_assert_near!(max.x.max(max.y), 1.0);
        prop_assert_near!(min.z, 0.0);
        prop_assert_near!(max.z, 1.0);
    }
}

#[test]
fn singular_camera_is_undefined_instead_of_aborting() {
    let camera = Camera {
        matrix: Matrix4::from_value(0.0),
        ..Camera::default()
    };
    assert!(camera.try_projection(1.0).is_none());
    assert!(camera.projection(1.0).x.x.is_nan());
    let ray = camera.ray(Point2::origin());
    assert!(ray.direction().x.is_nan());
}

#[test]
fn singular_view_basis_is_rejected() {
    let direction = Matrix3::from_value(0.0);
    let points = [Point3::origin(), Point3::new(1.0, 1.0, 1.0)];
    assert!(Camera::parallel_view_fitting(direction, 1.0, 0.1, &points).is_none());
    assert!(Camera::perspective_view_fitting(direction, 1.0, Rad(1.0), &points).is_none());
}
