//! End-to-end recovery: tessellated STL -> native analytic planar solid -> STEP.
//!
//! This fixture intentionally introduces an extra midpoint in just one side
//! of a cube's shared seam so the input is not a conforming triangulation.
//! The exported STEP must contain six analytic faces, not thirteen triangles.

use monstertruck_io::step::load::{Table, step_geometry};
use monstertruck_io::step::save::{CompleteStepDisplay, StepModel};
use monstertruck_mesh::stl::{self, StlFace, StlType};
use monstertruck_modeling::mesh_recovery::{PlanarRecoveryOptions, recover_planar_mesh};
use monstertruck_modeling::{Solid, ToCompressedTrimmedParameterCurves};
use monstertruck_topology::shell::ShellCondition;

#[test]
fn nonconforming_stl_cube_recovers_analytic_step_solid() {
    let positions: [[f32; 3]; 9] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
        [0.5, 0.0, 1.0],
    ];
    let faces: [[usize; 3]; 13] = [
        [0, 2, 1],
        [0, 3, 2],
        [4, 8, 6],
        [8, 5, 6],
        [4, 6, 7],
        [0, 1, 5],
        [0, 5, 4],
        [1, 2, 6],
        [1, 6, 5],
        [2, 3, 7],
        [2, 7, 6],
        [3, 0, 4],
        [3, 4, 7],
    ];
    let facets = faces.map(|indices| StlFace {
        normal: [0.0, 0.0, 0.0],
        vertices: indices.map(|id| positions[id]),
    });
    let mut stl_bytes = Vec::new();
    stl::write(facets, &mut stl_bytes, StlType::Binary).unwrap();
    let polygon_mesh = stl::read(stl_bytes.as_slice(), StlType::Automatic).unwrap();

    let recovery = recover_planar_mesh(&polygon_mesh, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovery.shell.len(), 6);
    assert_eq!(recovery.shell.shell_condition(), ShellCondition::Closed);
    assert!(
        recovery
            .patches
            .iter()
            .all(|patch| patch.boundary_loops.len() == 1 && patch.boundary_loops[0].len() == 4)
    );
    let solid = Solid::new(vec![recovery.shell]).unwrap();
    let compressed = solid.compress_with_parameter_curves(1.0e-6);
    assert_eq!(compressed.boundaries[0].vertices.len(), 8);
    assert_eq!(compressed.boundaries[0].edges.len(), 12);
    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();
    assert!(step.contains("MANIFOLD_SOLID_BREP"));
    assert_eq!(step.matches("ADVANCED_FACE").count(), 6);
    let table = Table::from_step(&step).unwrap();
    let shell = table.shell.values().next().expect("STEP shell exists");
    let roundtrip = table.to_compressed_trimmed_shell(shell).unwrap();
    assert_eq!(roundtrip.vertices.len(), 8);
    assert_eq!(roundtrip.edges.len(), 12);
    assert_eq!(roundtrip.faces.len(), 6);
    assert!(roundtrip.faces.iter().all(|face| matches!(
        face.surface,
        step_geometry::Surface::Elementary(step_geometry::ElementarySurface::Plane(_))
    )));
}

#[test]
fn sampled_cylindrical_stl_roundtrips_as_analytic_step_cylinders() {
    use monstertruck_modeling::mesh_recovery::{
        CylinderRecoveryOptions, recover_cylindrical_solid,
    };
    use monstertruck_modeling::{InnerSpace, ParametricSurface, Vector3};
    let sides = 32;
    let rings = (0..2)
        .flat_map(|level| {
            (0..sides).map(move |i| {
                let angle = std::f32::consts::TAU * i as f32 / sides as f32;
                [2.0 * angle.cos(), 2.0 * angle.sin(), 3.0 * level as f32]
            })
        })
        .collect::<Vec<_>>();
    let mut triangles = Vec::new();
    for i in 0..sides {
        let j = (i + 1) % sides;
        let bottom = [0.0; 3];
        let top = [0.0, 0.0, 3.0];
        triangles.extend([
            [bottom, rings[j], rings[i]],
            [top, rings[sides + i], rings[sides + j]],
            [rings[i], rings[j], rings[sides + j]],
            [rings[i], rings[sides + j], rings[sides + i]],
        ]);
    }
    let facets = triangles
        .into_iter()
        .map(|vertices| StlFace {
            normal: [0.0; 3],
            vertices,
        })
        .collect::<Vec<_>>();
    let mut encoded = Vec::new();
    stl::write(facets, &mut encoded, StlType::Binary).unwrap();
    let source = stl::read(encoded.as_slice(), StlType::Automatic).unwrap();
    let recovered = recover_cylindrical_solid(&source, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("STL cylinder should fit with bounded chord deviation");
    assert!((recovered.radius - 2.0).abs() < 1.0e-5);
    assert!((recovered.height - 3.0).abs() < 1.0e-5);
    assert!(recovered.max_surface_deviation < 0.05);
    assert_eq!(recovered.lateral_source_face_indices.len(), 64);
    let compressed = recovered.solid.compress_with_parameter_curves(1.0e-6);
    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();
    assert!(step.contains("MANIFOLD_SOLID_BREP"));
    assert!(step.contains("CYLINDRICAL_SURFACE("));
    assert!(!step.contains("B_SPLINE_SURFACE"));
    let table = Table::from_step(&step).unwrap();
    let shell = table.shell.values().next().unwrap();
    let reloaded = table.to_compressed_trimmed_shell(shell).unwrap();
    assert_eq!(reloaded.faces.len(), 4);
    let mut planes = 0;
    let mut cylinders = 0;
    for face in &reloaded.faces {
        match &face.surface {
            step_geometry::Surface::Elementary(step_geometry::ElementarySurface::Plane(_)) => {
                planes += 1;
            }
            step_geometry::Surface::Elementary(
                step_geometry::ElementarySurface::CylindricalSurface(cylinder),
            ) => {
                cylinders += 1;
                for axial in [0.0, 0.5, 1.0] {
                    for angle in [0.0, 1.0, 2.0] {
                        let p = cylinder.evaluate(axial, angle);
                        let radial = Vector3::new(p.x, p.y, 0.0);
                        assert!((radial.magnitude() - 2.0).abs() < 1.0e-5);
                    }
                }
            }
            other => panic!("unexpected STEP surface: {other:?}"),
        }
    }
    assert_eq!((planes, cylinders), (2, 2));
}
