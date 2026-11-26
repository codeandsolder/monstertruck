use super::{MeshRecoveryError, PlanarRecoveryOptions, recover_planar_mesh};
use crate::{Point3, Solid};
use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
use monstertruck_topology::shell::ShellCondition;

fn mesh(positions: &[[f64; 3]], faces: &[[usize; 3]]) -> PolygonMesh {
    let positions = positions
        .iter()
        .map(|&[x, y, z]| Point3::new(x, y, z))
        .collect();
    PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        Faces::from_iter(faces),
    )
    .expect("valid polygon indices")
}

fn square_vertices() -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ]
}

fn assert_separate_triangle_regions(positions: &[[f64; 3]; 6]) {
    let source = mesh(positions, &[[0, 1, 2], [3, 4, 5]]);
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 2);
    assert!(
        recovered
            .patches
            .iter()
            .all(|p| p.boundary_loops[0].len() == 3)
    );
}

#[test]
fn two_coplanar_triangles_become_one_planar_brep_face() {
    let source = mesh(&square_vertices(), &[[0, 1, 2], [0, 2, 3]]);
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 1);
    assert_eq!(recovered.shell.len(), 1);
    assert_eq!(recovered.patches[0].source_face_indices, [0, 1]);
    assert_eq!(recovered.patches[0].boundary_loops.len(), 1);
    assert_eq!(recovered.patches[0].boundary_loops[0].len(), 4);
    assert_eq!(recovered.patches[0].max_planar_residual, 0.0);
    assert!(matches!(
        recovered.shell[0].oriented_surface(),
        crate::Surface::Plane(_)
    ));
    assert_eq!(recovered.shell.shell_condition(), ShellCondition::Oriented);
}

#[test]
fn triangulated_cube_becomes_six_shared_edge_faces_and_closed_solid() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        &[
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.shell.len(), 6);
    assert_eq!(recovered.canonical_positions.len(), 8);
    assert!(
        recovered
            .patches
            .iter()
            .all(|patch| patch.source_face_indices.len() == 2)
    );
    assert!(
        recovered
            .patches
            .iter()
            .all(|patch| patch.boundary_loops[0].len() == 4)
    );
    assert_eq!(recovered.shell.shell_condition(), ShellCondition::Closed);
    assert!(Solid::new(vec![recovered.shell]).is_ok());
}

#[test]
fn ring_retains_oppositely_wound_hole() {
    let mut triangles = Vec::new();
    for i in 0..4 {
        let next = (i + 1) % 4;
        triangles.push([i, next, next + 4]);
        triangles.push([i, next + 4, i + 4]);
    }
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [4.0, 4.0, 0.0],
            [0.0, 4.0, 0.0],
            [1.0, 1.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 3.0, 0.0],
            [1.0, 3.0, 0.0],
        ],
        &triangles,
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.shell.len(), 1);
    assert_eq!(recovered.patches[0].boundary_loops.len(), 2);
    assert_eq!(recovered.shell[0].boundaries().len(), 2);
    assert!(
        recovered.patches[0]
            .boundary_loops
            .iter()
            .all(|loop_| loop_.len() == 4)
    );
}

#[test]
fn coincident_but_distinct_vertex_indices_are_welded() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        &[[0, 1, 2], [4, 5, 3]],
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 1);
    assert_eq!(recovered.patches[0].boundary_loops[0].len(), 4);
}

#[test]
fn coplanar_disconnected_islands_are_separate_faces() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ],
        &[[0, 1, 2], [3, 4, 5]],
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 2);
    assert!(
        recovered
            .patches
            .iter()
            .all(|patch| patch.boundary_loops.len() == 1)
    );
}

#[test]
fn a_dihedral_edge_is_not_merged_into_one_plane() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.1],
        ],
        &[[0, 1, 2], [1, 3, 2]],
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 2);
    assert_eq!(recovered.shell.shell_condition(), ShellCondition::Oriented);
}

#[test]
fn rejects_non_manifold_edges() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 1, 2], [1, 0, 3], [0, 1, 4]],
    );
    assert!(matches!(
        recover_planar_mesh(&source, PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::NonManifoldEdge(_)),
    ));
}

#[test]
fn rejects_inconsistent_winding() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        &[[0, 1, 2], [0, 1, 3]],
    );
    assert!(matches!(
        recover_planar_mesh(&source, PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::InconsistentWinding(_)),
    ));
}

#[test]
fn rejects_invalid_tolerances() {
    let source = mesh(&square_vertices(), &[[0, 1, 2]]);
    let options = PlanarRecoveryOptions {
        distance_tolerance: f64::NAN,
        ..Default::default()
    };
    assert!(matches!(
        recover_planar_mesh(&source, options),
        Err(MeshRecoveryError::InvalidTolerance)
    ));
}

#[test]
fn rejects_degenerate_face_after_exact_position_welding() {
    let source = mesh(
        &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        &[[0, 1, 2]],
    );
    assert!(matches!(
        recover_planar_mesh(&source, PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::DegeneratePolygon(0)),
    ));
}

#[test]
fn empty_polygon_mesh_produces_empty_shell() {
    let source = mesh(&[], &[]);
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert!(recovered.patches.is_empty());
    assert!(recovered.shell.is_empty());
}

#[test]
fn binary_stl_import_recovers_the_original_square_face() {
    use monstertruck_mesh::stl::{self, StlFace, StlType};
    let facets = [
        StlFace {
            normal: [0.0, 0.0, 1.0],
            vertices: [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        },
        StlFace {
            normal: [0.0, 0.0, 1.0],
            vertices: [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        },
    ];
    let mut binary = Vec::new();
    stl::write(facets, &mut binary, StlType::Binary).unwrap();
    let source = stl::read(binary.as_slice(), StlType::Automatic).unwrap();
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.shell.len(), 1);
    assert_eq!(recovered.patches[0].boundary_loops[0].len(), 4);
}

#[test]
fn rejects_warped_quadrangle_instead_of_pretending_it_is_a_plane() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.1),
    ];
    let source = PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        Faces::from_iter([[0, 1, 2, 3]]),
    )
    .unwrap();
    assert!(matches!(
        recover_planar_mesh(&source, PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::NonPlanarPolygon(0)),
    ));
}

#[test]
fn rejects_bad_indices_even_from_an_unchecked_mesh_constructor() {
    let source = PolygonMesh::new_unchecked(
        StandardAttributes {
            positions: square_vertices()
                .iter()
                .map(|&[x, y, z]| Point3::new(x, y, z))
                .collect(),
            ..Default::default()
        },
        Faces::from_iter([[0, 1, 50]]),
    );
    assert!(matches!(
        recover_planar_mesh(&source, PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::InvalidPositionIndex {
            polygon: 0,
            position: 50
        }),
    ));
}

#[test]
fn rejects_loose_surface_fit_that_would_break_brep_incidence() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.1],
        ],
        &[[0, 1, 2], [1, 3, 2]],
    );
    let options = PlanarRecoveryOptions {
        distance_tolerance: 0.5,
        angular_tolerance_radians: 0.5,
    };
    assert!(matches!(
        recover_planar_mesh(&source, options),
        Err(MeshRecoveryError::ExcessivePlanarResidual { .. }),
    ));
}

#[test]
fn dense_planar_grid_collapses_without_losing_boundary_segments() {
    const N: usize = 64;
    let mut positions = Vec::new();
    for y in 0..=N {
        for x in 0..=N {
            positions.push([x as f64, y as f64, 0.0]);
        }
    }
    let mut triangles = Vec::new();
    for y in 0..N {
        for x in 0..N {
            let a = y * (N + 1) + x;
            let b = a + 1;
            let c = a + N + 1;
            let d = c + 1;
            triangles.extend([[a, b, d], [a, d, c]]);
        }
    }
    let source = mesh(&positions, &triangles);
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 1);
    assert_eq!(recovered.patches[0].source_face_indices.len(), 2 * N * N);
    assert_eq!(recovered.patches[0].boundary_loops[0].len(), 4);
}

#[test]
fn geometric_t_junction_in_cube_is_split_and_then_globally_simplified() {
    // The upper face has a midpoint on a seam; the adjoining side uses
    // the original unsplit long edge. The patch topology must agree.
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
            [0.5, 0.0, 1.0],
        ],
        &[
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
        ],
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 6);
    assert_eq!(recovered.shell.shell_condition(), ShellCondition::Closed);
    assert!(
        recovered
            .patches
            .iter()
            .all(|patch| { patch.boundary_loops.len() == 1 && patch.boundary_loops[0].len() == 4 })
    );
    // No orphaned half-edge or T-junction survives as a B-rep edge.
    for face in recovered.shell.face_iter() {
        assert_eq!(face.boundaries()[0].len(), 4);
    }
    assert!(Solid::new(vec![recovered.shell]).is_ok());
}

#[test]
fn near_but_not_incident_vertex_does_not_weld_to_a_seam() {
    // This point is 1e-8 above the plane, not connected to its carrier edge.
    assert_separate_triangle_regions(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.0, 1e-8],
        [0.5, 0.0, 1.0],
        [0.5, 1.0, 1e-8],
    ]);
}

#[test]
fn a_nearly_straight_real_boundary_feature_is_not_erased() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [0.5, 1.0e-8, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        &[[0, 1, 3], [1, 2, 3]],
    );
    let recovered = recover_planar_mesh(&source, PlanarRecoveryOptions::default()).unwrap();
    assert_eq!(recovered.patches.len(), 1);
    assert_eq!(recovered.patches[0].boundary_loops[0].len(), 4);
    assert!(recovered.patches[0].boundary_loops[0].contains(&1));
}

#[test]
fn t_junction_repair_rejects_a_third_incident_polygon() {
    let source = mesh(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
        ],
        &[[0, 1, 2], [1, 0, 3], [0, 4, 5]],
    );
    assert!(matches!(
        recover_planar_mesh(&source, PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::NonManifoldEdge(_)),
    ));
}

#[test]
fn point_only_contact_on_an_edge_does_not_create_a_false_seam() {
    // A vertex on the edge with only transverse spokes is a point contact,
    // not a T-junction with overlapping collinear source edges.
    assert_separate_triangle_regions(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.0, 1.0],
        [0.5, 1.0, 1.0],
    ]);
}

#[test]
fn large_coordinate_offsets_do_not_turn_a_nearby_parallel_edge_into_a_seam() {
    // At 1e12, 8*epsilon*coordinate is 0.0018: larger than this real gap.
    let base = 1.0e12;
    assert_separate_triangle_regions(&[
        [base, base, 0.0],
        [base + 1.0, base, 0.0],
        [base, base + 1.0, 0.0],
        [base + 0.5, base - 0.001, 0.0],
        [base + 0.5, base - 0.5, 0.0],
        [base + 0.75, base - 0.001, 0.0],
    ]);
}
