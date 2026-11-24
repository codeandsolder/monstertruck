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

#[test]
fn faceted_through_hole_stl_roundtrips_as_mixed_planar_cylindrical_step() {
    use monstertruck_modeling::mesh_recovery::{
        CylinderRecoveryOptions, recover_cylindrical_through_hole,
    };
    // A checked-in binary STL fixture is used here rather than generating
    // the mesh from the recovery code under test or linking synthetic CAD.
    // 256 triangles: two triangulated annular caps, four outer planar sides,
    // and a 32-segment faceted cylindrical hole.
    let source = stl::read(
        &include_bytes!("../fixtures/mesh-recovery/through-hole-32.stl")[..],
        StlType::Automatic,
    )
    .unwrap();
    let recovered = recover_cylindrical_through_hole(&source, CylinderRecoveryOptions::default())
        .unwrap()
        .expect("circular STL hole must recover a true analytic through-hole");
    assert!((recovered.radius - 1.0).abs() < 1.0e-5);
    assert!((recovered.height - 4.0).abs() < 1.0e-5);
    assert!(recovered.max_surface_deviation < 0.05);
    assert_eq!(recovered.cylindrical_source_face_indices.len(), 64);
    assert!(recovered.solid.is_geometric_consistent());
    assert_eq!(
        recovered.solid.boundaries()[0].shell_condition(),
        ShellCondition::Closed
    );

    let compressed = recovered.solid.compress_with_parameter_curves(1.0e-6);
    assert_eq!(compressed.boundaries[0].faces.len(), 8);
    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();
    assert!(step.contains("MANIFOLD_SOLID_BREP"));
    assert!(step.contains("CYLINDRICAL_SURFACE("));
    assert!(!step.contains("B_SPLINE_SURFACE"));
    let parsed = Table::from_step(&step).unwrap();
    let native_shell = parsed.shell.values().next().unwrap();
    let reloaded = parsed.to_compressed_trimmed_shell(native_shell).unwrap();
    assert_eq!(reloaded.faces.len(), 8);
    assert_eq!(
        reloaded
            .faces
            .iter()
            .filter(|f| f.boundaries.len() == 2)
            .count(),
        2
    );
    let (mut planes, mut cylinders) = (0, 0);
    for face in &reloaded.faces {
        match &face.surface {
            step_geometry::Surface::Elementary(step_geometry::ElementarySurface::Plane(_)) => {
                planes += 1;
            }
            step_geometry::Surface::Elementary(
                step_geometry::ElementarySurface::CylindricalSurface(_),
            ) => {
                cylinders += 1;
                assert_eq!(face.boundaries.len(), 1);
            }
            other => panic!("unexpected recovered STEP face surface: {other:?}"),
        }
    }
    assert_eq!((planes, cylinders), (6, 2));
}

#[test]
fn shared_periodic_trim_graph_roundtrips_as_analytic_step_solid()
-> Result<(), Box<dyn std::error::Error>> {
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_modeling::mesh_recovery::{
        EvidenceOptions, NativeTrimOptions, RegionBoundaryIndex, SharedTrimOptions, SurfaceCarrier,
        SurfaceRegionBudget, analyze_mesh_evidence,
    };
    use monstertruck_modeling::{Point3, Shell, Vector3};

    let sides = 16_usize;
    let mut positions = Vec::new();
    for z in [0.0, 3.0] {
        for i in 0..sides {
            let (sin, cos) = (std::f64::consts::TAU * i as f64 / sides as f64).sin_cos();
            positions.push(Point3::new(2.0 * cos, 2.0 * sin, z));
        }
    }
    let mut faces = (0..sides)
        .map(|i| {
            let j = (i + 1) % sides;
            vec![i, j, j + sides, i + sides]
        })
        .collect::<Vec<_>>();
    faces.push((0..sides).rev().collect());
    faces.push((sides..2 * sides).collect());
    let mesh = PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        Faces::from_iter(faces),
    )?;
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
    let boundary = RegionBoundaryIndex::new(&evidence)?;
    let mut supports = evidence.source_planar_support[..sides].to_vec();
    supports.sort_unstable();
    supports.dedup();
    let proposal = boundary
        .measure_surface_region(
            &supports,
            SurfaceCarrier::LinearRevolution {
                origin: Point3::new(0.0, 0.0, 0.0),
                axis: Vector3::unit_z(),
                radius_at_origin: 2.0,
                slope: 0.0,
            },
            SurfaceRegionBudget {
                vertex: 1.0e-8,
                chord: 0.15,
            },
        )
        .ok_or("cylinder proposal")?;
    let partition = boundary.select_surface_partition(&[proposal], &evidence, 100)?;
    let graph = boundary.assemble_shared_trims(
        &partition,
        &evidence.positions,
        SharedTrimOptions {
            absolute_tolerance: 1.0e-8,
            relative_tolerance: 1.0e-4,
            max_interval_edges: 128,
        },
    )?;
    let native = boundary.construct_native_trims(
        &partition,
        &graph,
        &evidence.positions,
        NativeTrimOptions {
            max_vertex_displacement: 1.0e-8,
            intersection_tolerance: 1.0e-10,
            angular_resolution: 1.0e-8,
        },
    )?;
    assert!(native.regions.iter().all(|r| !r.requires_chart_cut));
    let mut shell = Shell::new();
    for region in native.regions {
        shell.extend(region.faces?);
    }
    assert_eq!(shell.len(), 4);
    assert_eq!(shell.shell_condition(), ShellCondition::Closed);
    let solid = Solid::new(vec![shell])?;
    assert!(solid.is_geometric_consistent());
    let compressed = solid.compress_with_parameter_curves(1.0e-6);
    let step =
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string();
    assert!(step.contains("MANIFOLD_SOLID_BREP"));
    assert!(step.contains("CYLINDRICAL_SURFACE("));
    assert!(!step.contains("B_SPLINE_SURFACE"));
    let table = Table::from_step(&step)?;
    let shell = table.shell.values().next().ok_or("STEP shell")?;
    let restored = table.to_compressed_trimmed_shell(shell)?;
    assert_eq!(restored.faces.len(), 4);
    assert_eq!(
        restored
            .faces
            .iter()
            .filter(|face| matches!(
                face.surface,
                step_geometry::Surface::Elementary(
                    step_geometry::ElementarySurface::CylindricalSurface(_)
                )
            ))
            .count(),
        2
    );
    Ok(())
}
