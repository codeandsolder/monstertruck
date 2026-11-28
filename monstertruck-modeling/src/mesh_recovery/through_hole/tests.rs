use super::*;
use crate::mesh_recovery::signed_source_volume;
use crate::{ParametricSurface, ParametricSurface3D, Surface};
use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
use monstertruck_topology::shell::ShellCondition;

/// A square prism of material with a cylindrical hole. Its circular source
/// boundary is an N-gon; each cap is a triangulated annulus and both the
/// exterior walls and the faceted hole wall use the same exact edge indices.
fn perforated_prism(samples: usize, radius: f64, half_width: f64, height: f64) -> PolygonMesh {
    let origin = Point3::new(3.0, -4.0, 2.0);
    let rings = (0..4)
        .flat_map(|ring| {
            (0..samples).map(move |i| {
                let angle =
                    std::f64::consts::TAU * crate::usize_to_f64(i) / crate::usize_to_f64(samples);
                let (sine, cosine) = angle.sin_cos();
                let scale = if ring < 2 {
                    radius
                } else {
                    half_width / cosine.abs().max(sine.abs())
                };
                Point3::new(
                    origin.x + scale * cosine,
                    origin.y + scale * sine,
                    origin.z + if ring % 2 == 0 { 0.0 } else { height },
                )
            })
        })
        .collect::<Vec<_>>();
    let mut triangles = Vec::new();
    for i in 0..samples {
        let j = (i + 1) % samples;
        let ib = i;
        let it = samples + i;
        let ob = samples * 2 + i;
        let ot = samples * 3 + i;
        let ibj = j;
        let itj = samples + j;
        let obj = samples * 2 + j;
        let otj = samples * 3 + j;
        triangles.extend([
            [ib, ibj, obj],
            [ib, obj, ob],
            [it, ot, otj],
            [it, otj, itj],
            [ib, it, itj],
            [ib, itj, ibj],
            [ob, obj, otj],
            [ob, otj, ot],
        ]);
    }
    PolygonMesh::new(
        StandardAttributes {
            positions: rings,
            ..Default::default()
        },
        Faces::from_iter(triangles),
    )
    .unwrap()
}

#[test]
fn thirty_two_sided_source_hole_becomes_shared_circular_analytic_trim() {
    let mesh = perforated_prism(32, 1.0, 3.0, 4.0);
    let recovered = recover_cylindrical_through_hole(&mesh, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("single cylindrical through-hole must be reconstructible");
    assert!((recovered.radius - 1.0).abs() < 1.0e-9);
    assert!((recovered.height - 4.0).abs() < 1.0e-9);
    assert!(recovered.source_signed_volume > 0.0);
    assert!(recovered.max_radial_residual < 1.0e-8);
    assert!(recovered.max_surface_deviation < 0.05);
    assert_eq!(recovered.cylindrical_source_face_indices.len(), 64);
    assert_eq!(recovered.cap_source_face_indices[0].len(), 64);
    assert_eq!(recovered.cap_source_face_indices[1].len(), 64);
    let shell = &recovered.solid.boundaries()[0];
    assert_eq!(shell.shell_condition(), ShellCondition::Closed);
    assert!(shell.is_connected());
    assert!(recovered.solid.is_geometric_consistent());
    let planar = shell
        .face_iter()
        .filter(|f| matches!(f.surface(), Surface::Plane(_)))
        .count();
    let cylindrical = shell
        .face_iter()
        .filter(|f| matches!(f.surface(), Surface::Revolution(_)))
        .count();
    assert_eq!((planar, cylindrical), (6, 2));
    assert_eq!(shell.len(), 8);
    // Topologically closed alone is not enough: the hole wall normals must
    // point into the void, even though the exterior normals point outwards.
    for face in shell.face_iter() {
        if matches!(face.surface(), Surface::Revolution(_)) {
            let surface = face.oriented_surface();
            let point = surface.evaluate(0.4, 0.5);
            let delta = point - recovered.cap_centers[0];
            let radial = delta - recovered.axis * delta.dot(&recovered.axis);
            assert!(surface.normal(0.4, 0.5).dot(&radial) < 0.0);
        }
    }
    // Each annular end cap retains exactly one hole bounded by two native
    // circular arcs, identical in topological ID to the cylinder side edges.
    assert_eq!(
        shell
            .face_iter()
            .filter(|f| f.boundaries().len() == 2)
            .count(),
        2
    );
}

#[test]
fn non_circular_hole_is_not_promoted_to_cylinder() {
    let mut mesh = perforated_prism(32, 1.0, 3.0, 4.0);
    for point in mesh.positions_mut().iter_mut().take(64) {
        point.x = 3.0 + (point.x - 3.0) * 1.1;
    }
    assert!(
        recover_cylindrical_through_hole(&mesh, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn low_polygon_count_requires_an_explicit_faceting_budget() {
    let mesh = perforated_prism(8, 1.0, 3.0, 4.0);
    let tight = CylinderRecoveryOptions {
        min_ring_vertices: 8,
        max_surface_deviation: 0.02,
        ..Default::default()
    };
    assert!(
        recover_cylindrical_through_hole(&mesh, tight)
            .unwrap()
            .is_none()
    );
    let loose = CylinderRecoveryOptions {
        max_surface_deviation: 0.1,
        ..tight
    };
    let recovered = recover_cylindrical_through_hole(&mesh, loose)
        .unwrap()
        .expect("coarse cylindrical hole fits within explicit budget");
    assert!(recovered.max_surface_deviation > 0.07);
    assert!(recovered.max_surface_deviation < 0.1);
}

#[test]
fn shifted_upper_hole_is_not_coaxial() {
    let mut mesh = perforated_prism(32, 1.0, 3.0, 4.0);
    for point in mesh.positions_mut().iter_mut().skip(32).take(32) {
        point.x += 0.04;
    }
    assert!(
        recover_cylindrical_through_hole(&mesh, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn inverted_mesh_cannot_infer_an_outward_material_through_hole() {
    let mut mesh = perforated_prism(32, 1.0, 3.0, 4.0);
    for polygon in mesh.face_iter_mut() {
        polygon.reverse();
    }
    assert!(
        recover_cylindrical_through_hole(&mesh, CylinderRecoveryOptions::default())
            .unwrap()
            .is_none()
    );
}

#[test]
fn inclined_source_axis_preserves_analytic_through_hole_and_winding() {
    let mut mesh = perforated_prism(24, 1.3, 4.0, 2.5);
    // Proper rigid 90 degree rotation about X, followed by translation.
    for point in mesh.positions_mut() {
        let (_, y, z) = (point.x, point.y, point.z);
        point.x += 11.0;
        point.y = z + 8.0;
        point.z = 7.0 - y;
    }
    let recovered = recover_cylindrical_through_hole(&mesh, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("arbitrary source orientations must recover");
    assert!((recovered.radius - 1.3).abs() < 1.0e-8);
    assert!((recovered.height - 2.5).abs() < 1.0e-8);
    assert!(recovered.axis.dot(&Vector3::y()) > 0.9999);
    assert!(recovered.solid.is_geometric_consistent());
}

#[test]
fn three_level_wall_with_intermediate_ring_still_recovers_as_one_carrier() {
    let samples = 24;
    let original = perforated_prism(samples, 1.0, 3.0, 4.0);
    let mut positions = original.positions().to_vec();
    let middle = positions.len();
    for i in 0..samples {
        let mut point = positions[i];
        point.z += 2.0;
        positions.push(point);
    }
    let mut triangles = Vec::new();
    for (index, polygon) in original.face_iter().enumerate() {
        if index % 8 != 4 && index % 8 != 5 {
            triangles.push(polygon.iter().map(|v| v.pos).collect::<Vec<_>>());
        }
    }
    for i in 0..samples {
        let j = (i + 1) % samples;
        let bottom = i;
        let top = samples + i;
        let mb = middle + i;
        let mj = middle + j;
        triangles.extend([
            vec![bottom, mb, mj],
            vec![bottom, mj, j],
            vec![mb, top, samples + j],
            vec![mb, samples + j, mj],
        ]);
    }
    let source = PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        Faces::from_iter(triangles),
    )
    .unwrap();
    let recovered = recover_cylindrical_through_hole(&source, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("middle wall ring must not become a spurious B-rep seam");
    assert_eq!(recovered.cylindrical_source_face_indices.len(), samples * 4);
    assert_eq!(recovered.solid.boundaries()[0].len(), 8);
    assert!(recovered.solid.is_geometric_consistent());
}

#[test]
fn signed_mesh_volume_detects_globally_inverted_material_and_is_translation_stable() {
    let source = perforated_prism(32, 1.0, 3.0, 4.0);
    let (positions, polygons) = positions_and_polygons(&source, KERNEL_EPSILON).unwrap();
    let signed = signed_source_volume(&positions, &polygons).unwrap();
    let expected_hole_area = 16.0 * (std::f64::consts::TAU / 32.0).sin();
    let expected_material_volume = (36.0 - expected_hole_area) * 4.0;
    assert!((signed - expected_material_volume).abs() < 1.0e-8);
    let offset = positions
        .iter()
        .map(|p| *p + Vector3::new(1.0e9, -1.0e9, 1.0e9))
        .collect::<Vec<_>>();
    let translated = signed_source_volume(&offset, &polygons).unwrap();
    assert!((translated - signed).abs() < 1.0e-4);
    let mut reversed = polygons;
    for polygon in &mut reversed {
        polygon.vertices.reverse();
    }
    let inverted_volume = signed_source_volume(&positions, &reversed).unwrap();
    assert!((inverted_volume + signed).abs() < 1.0e-8);
}
