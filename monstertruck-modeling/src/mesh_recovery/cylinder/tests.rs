use super::*;
use crate::{ParametricSurface, ParametricSurface3D, Surface};
use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
use monstertruck_topology::shell::ShellCondition;

fn regular_cylinder(divisions: usize, radius: f64, height: f64, axis: Vector3) -> PolygonMesh {
    segmented_cylinder(divisions, radius, height, axis, 1, false)
}

fn segmented_cylinder(
    divisions: usize,
    radius: f64,
    height: f64,
    axis: Vector3,
    levels: usize,
    staggered: bool,
) -> PolygonMesh {
    let origin = Point3::new(2.0, -3.0, -4.0);
    let u = Vector3::unit_x();
    let v = axis.cross(u);
    let positions = (0..=levels)
        .flat_map(|level| {
            (0..divisions).map(move |i| {
                let offset = if staggered && level % 2 == 1 {
                    0.5
                } else {
                    0.0
                };
                let angle = std::f64::consts::TAU * (i as f64 + offset) / divisions as f64;
                origin
                    + axis * (height * level as f64 / levels as f64)
                    + u * (radius * angle.cos())
                    + v * (radius * angle.sin())
            })
        })
        .chain([origin, origin + axis * height])
        .collect::<Vec<_>>();
    let mut triangles = Vec::new();
    let upper = levels * divisions;
    let lower_center = (levels + 1) * divisions;
    let upper_center = lower_center + 1;
    for i in 0..divisions {
        let j = (i + 1) % divisions;
        triangles.extend([[lower_center, j, i], [upper_center, upper + i, upper + j]]);
        for level in 0..levels {
            let bottom = level * divisions;
            let top = bottom + divisions;
            triangles.extend([
                [bottom + i, bottom + j, top + j],
                [bottom + i, top + j, top + i],
            ]);
        }
    }
    PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        Faces::from_iter(triangles),
    )
    .unwrap()
}

#[test]
fn revolved_line_with_shared_circular_caps_makes_closed_native_cylinder() {
    let solid = build_cylinder_solid(
        Point3::new(2.0, -3.0, -4.0),
        Vector3::unit_z(),
        Vector3::unit_x(),
        2.5,
        7.0,
    )
    .unwrap();
    assert_eq!(solid.boundaries().len(), 1);
    let shell = &solid.boundaries()[0];
    assert_eq!(shell.shell_condition(), ShellCondition::Closed);
    assert_eq!(shell.len(), 4);
    assert_eq!(
        shell
            .face_iter()
            .filter(|f| matches!(f.surface(), Surface::Plane(_)))
            .count(),
        2
    );
    assert_eq!(
        shell
            .face_iter()
            .filter(|f| matches!(f.surface(), Surface::Revolution(_)))
            .count(),
        2
    );
    for face in shell.face_iter() {
        let surface = face.oriented_surface();
        let normal = surface.normal(0.5, 0.5);
        match surface {
            Surface::Revolution(_) => {
                let point = surface.evaluate(0.5, 0.5);
                let radial = point - Point3::new(2.0, -3.0, -4.0);
                let radial = radial - Vector3::unit_z() * radial.z;
                assert!(
                    normal.dot(radial) > 0.0,
                    "lateral cylinder normal points inward"
                );
            }
            Surface::Plane(_) => {
                let point = face.boundaries()[0][0].front().point();
                let expected = if point.z < 0.0 {
                    -Vector3::unit_z()
                } else {
                    Vector3::unit_z()
                };
                assert!(
                    normal.dot(expected) > 0.0,
                    "cylindrical cap normal points inward"
                );
            }
            _ => unreachable!(),
        }
    }
}

#[test]
fn recovers_true_analytic_cylinder_from_128_triangles_with_cap_center_vertices() {
    let input = regular_cylinder(32, 3.0, 7.0, Vector3::unit_z());
    let result = recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("closed faceted cylinder should fit");
    assert!((result.radius - 3.0).abs() < 1.0e-9);
    assert!((result.height - 7.0).abs() < 1.0e-9);
    assert!(result.max_radial_residual < 1.0e-9);
    assert!(result.max_surface_deviation > 0.0);
    assert!(result.max_surface_deviation < 0.05);
    assert_eq!(result.cap_source_face_indices[0].len(), 32);
    assert_eq!(result.cap_source_face_indices[1].len(), 32);
    assert_eq!(result.lateral_source_face_indices.len(), 64);
    let shell = &result.solid.boundaries()[0];
    assert_eq!(shell.shell_condition(), ShellCondition::Closed);
    assert_eq!(shell.len(), 4);
    assert_eq!(
        shell
            .face_iter()
            .filter(|f| matches!(f.surface(), Surface::Plane(_)))
            .count(),
        2
    );
    assert_eq!(
        shell
            .face_iter()
            .filter(|f| matches!(f.surface(), Surface::Revolution(_)))
            .count(),
        2
    );
}

#[test]
fn recovers_axis_along_y_with_nonzero_world_origin() {
    let input = regular_cylinder(24, 1.25, 4.0, Vector3::unit_y());
    let result = recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("arbitrarily aligned cap normal should be supported");
    assert!((result.radius - 1.25).abs() < 1.0e-9);
    assert!((result.height - 4.0).abs() < 1.0e-9);
    assert!((result.axis - Vector3::unit_y()).magnitude() < 1.0e-9);
}

#[test]
fn does_not_reinterpret_an_ellipse_as_a_cylinder() {
    let mut input = regular_cylinder(32, 2.0, 4.0, Vector3::unit_z());
    for vertex in input.positions_mut() {
        vertex.x = 2.0 + (vertex.x - 2.0) * 1.2;
    }
    assert!(
        recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn chord_error_budget_controls_polygonal_to_smooth_conversion() {
    let input = regular_cylinder(8, 3.0, 7.0, Vector3::unit_z());
    let options = CylinderRecoveryOptions {
        min_ring_vertices: 8,
        max_surface_deviation: 0.05,
        ..Default::default()
    };
    assert!(
        recover_cylindrical_solid(&input, options)
            .unwrap()
            .is_none()
    );
    let options = CylinderRecoveryOptions {
        max_surface_deviation: 0.3,
        ..options
    };
    let result = recover_cylindrical_solid(&input, options).unwrap().unwrap();
    assert!(result.max_surface_deviation > 0.2);
    assert!(result.max_surface_deviation < 0.3);
}

#[test]
fn does_not_hide_noncoaxial_ring_offset() {
    let mut input = regular_cylinder(32, 2.0, 4.0, Vector3::unit_z());
    for vertex in input.positions_mut().iter_mut().skip(32).take(32) {
        vertex.x += 0.01;
    }
    input.positions_mut()[65].x += 0.01;
    assert!(
        recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn open_cylindrical_sleeve_is_not_claimed_as_a_closed_solid() {
    let mut input = regular_cylinder(24, 1.25, 4.0, Vector3::unit_z());
    let triangles = (0..24)
        .flat_map(|i| {
            let j = (i + 1) % 24;
            [[i, j, 24 + j], [i, 24 + j, 24 + i]]
        })
        .collect::<Vec<_>>();
    input = PolygonMesh::new(
        StandardAttributes {
            positions: input.positions().to_vec(),
            ..Default::default()
        },
        Faces::from_iter(triangles),
    )
    .unwrap();
    assert!(
        recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn rejects_invalid_cylinder_fit_parameters() {
    let input = regular_cylinder(24, 1.0, 2.0, Vector3::unit_z());
    let options = CylinderRecoveryOptions {
        max_surface_deviation: f64::NAN,
        ..Default::default()
    };
    assert!(matches!(
        recover_cylindrical_solid(&input, options),
        Err(MeshRecoveryError::InvalidCylinderOptions)
    ));
}

#[test]
fn recovers_three_axial_sections_even_with_staggered_intermediate_rings() {
    let input = segmented_cylinder(32, 2.0, 5.0, Vector3::unit_z(), 3, true);
    let result = recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("intermediate axial levels must not prevent analytic recovery");
    assert_eq!(result.lateral_source_face_indices.len(), 192);
    assert_eq!(result.cap_source_face_indices[0].len(), 32);
    assert_eq!(result.cap_source_face_indices[1].len(), 32);
    assert!((result.radius - 2.0).abs() < 1.0e-9);
    assert!((result.height - 5.0).abs() < 1.0e-9);
    assert!(result.max_surface_deviation < 0.05);
    assert_eq!(result.solid.boundaries()[0].len(), 4);
}

#[test]
fn rejects_radially_bulged_middle_ring() {
    let mut input = segmented_cylinder(32, 2.0, 5.0, Vector3::unit_z(), 3, false);
    // Both cap rings remain perfect circles, but the interior is not cylindrical.
    for point in input.positions_mut().iter_mut().skip(32).take(32) {
        point.x = 2.0 + (point.x - 2.0) * 1.03;
        point.y = -3.0 + (point.y + 3.0) * 1.03;
    }
    assert!(
        recover_cylindrical_solid(&input, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}
