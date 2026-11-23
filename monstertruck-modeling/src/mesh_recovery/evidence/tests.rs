use super::*;
use crate::Vector3;
use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

fn make_mesh(positions: &[[f64; 3]], triangles: &[[usize; 3]]) -> PolygonMesh {
    let coordinates = positions
        .iter()
        .map(|&[x, y, z]| Point3::new(x, y, z))
        .collect();
    PolygonMesh::new(
        StandardAttributes {
            positions: coordinates,
            ..Default::default()
        },
        Faces::from_iter(triangles.iter().copied()),
    )
    .unwrap()
}

fn unit_box_mesh() -> PolygonMesh {
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let quads = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ];
    let triangles = quads
        .into_iter()
        .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
        .collect::<Vec<_>>();
    make_mesh(&corners, &triangles)
}

#[test]
fn cube_diagonal_triangulation_disappears_before_parallel_grouping() {
    let input = unit_box_mesh();
    let evidence = analyze_mesh_evidence(&input, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.positions.len(), 8);
    assert_eq!(evidence.suppressed_interior_edge_count, 6);
    assert_eq!(evidence.edges.len(), 12);
    assert_eq!(evidence.chains.len(), 12);
    assert_eq!(evidence.parallel_bundles.len(), 3);
    for bundle in &evidence.parallel_bundles {
        assert_eq!(bundle.placements.len(), 4);
        assert!(
            bundle
                .placements
                .iter()
                .all(|member| member.angular_residual_radians <= 1.0e-12)
        );
        let separations = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| bundle.transverse_separation(i, j).unwrap()))
            .collect::<Vec<_>>();
        assert_eq!(
            separations
                .iter()
                .filter(|&&d| (d - 1.0).abs() < 1.0e-12)
                .count(),
            4
        );
        assert_eq!(
            separations
                .iter()
                .filter(|&&d| (d - 2.0_f64.sqrt()).abs() < 1.0e-12)
                .count(),
            2
        );
    }
}

#[test]
fn exact_collinear_source_subdivisions_form_one_chain_with_provenance() {
    let input = make_mesh(
        &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        &[[0, 1, 4], [1, 2, 3], [1, 3, 4]],
    );
    let evidence = analyze_mesh_evidence(&input, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.suppressed_interior_edge_count, 2);
    assert_eq!(evidence.edges.len(), 5);
    assert_eq!(evidence.chains.len(), 4);
    assert_eq!(evidence.parallel_bundles.len(), 2);
    let subdivided = evidence
        .chains
        .iter()
        .find(|chain| chain.edge_indices.len() == 2)
        .expect("two collinear source segments must be one exact chain");
    assert_eq!(subdivided.vertices, [0, 1, 2]);
    assert_eq!(subdivided.kind, FeatureEdgeKind::Boundary);
}

#[test]
fn separate_nearby_parallel_lines_are_hypothesized_not_welded() {
    let small = 0.000_000_001;
    let input = make_mesh(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, small],
            [2.0, 0.0, small],
            [0.0, 1.0, small],
        ],
        &[[0, 1, 2], [3, 4, 5]],
    );
    let evidence = analyze_mesh_evidence(&input, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.positions.len(), 6);
    assert_eq!(evidence.edges.len(), 6);
    assert_eq!(evidence.chains.len(), 6);
    let grouped = evidence
        .parallel_bundles
        .iter()
        .find(|bundle| bundle.direction.dot(Vector3::unit_x()).abs() > 0.99999)
        .expect("same horizontal direction");
    assert_eq!(grouped.placements.len(), 2);
    assert!((grouped.transverse_separation(0, 1).unwrap() - small).abs() < 1.0e-15);
}

#[test]
fn faceted_through_hole_exposes_parallel_generatrices_without_cylinder_assumption() {
    let source = monstertruck_mesh::stl::read(
        &include_bytes!(
            "../../../../monstertruck-io/tests/fixtures/mesh-recovery/through-hole-32.stl"
        )[..],
        monstertruck_mesh::stl::StlType::Automatic,
    )
    .unwrap();
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    let axial = evidence
        .parallel_bundles
        .iter()
        .find(|bundle| bundle.direction.dot(Vector3::unit_z()).abs() > 0.99999)
        .expect("faceted wall has parallel axial sample edges");
    assert!(axial.placements.len() >= 32);
    assert!(axial.placements.iter().all(|member| {
        let span = member.longitudinal_interval[1] - member.longitudinal_interval[0];
        (span - 4.0).abs() < 1.0e-5
    }));
    assert!(
        evidence
            .edges
            .iter()
            .any(|edge| edge.kind == FeatureEdgeKind::Dihedral)
    );
    assert!(evidence.suppressed_interior_edge_count > 0);
}

#[test]
fn direction_cluster_does_not_grow_by_transitive_angular_drift() {
    let angles: [f64; 3] = [0.0, 0.75_f64.to_radians(), 1.5_f64.to_radians()];
    let mut coordinates = Vec::new();
    let mut triangles = Vec::new();
    for (i, angle) in angles.into_iter().enumerate() {
        let x = i as f64 * 20.0;
        let base = 10.0 - i as f64;
        let first = coordinates.len();
        coordinates.extend([
            [x, 0.0, 0.0],
            [x + base * angle.cos(), base * angle.sin(), 0.0],
            [x, 2.0, 0.0],
        ]);
        triangles.push([first, first + 1, first + 2]);
    }
    let input = make_mesh(&coordinates, &triangles);
    let evidence = analyze_mesh_evidence(
        &input,
        EvidenceOptions {
            max_parallel_angle_radians: 1.0_f64.to_radians(),
            ..Default::default()
        },
    )
    .unwrap();
    let seed = evidence
        .parallel_bundles
        .iter()
        .find(|bundle| {
            bundle.placements.iter().any(|placement| {
                let chain = &evidence.chains[placement.chain_index];
                chain.vertices.contains(&0) && chain.vertices.contains(&1)
            })
        })
        .unwrap();
    assert!(seed.placements.iter().any(|p| {
        let chain = &evidence.chains[p.chain_index];
        chain.vertices.contains(&3) && chain.vertices.contains(&4)
    }));
    assert!(!seed.placements.iter().any(|p| {
        let chain = &evidence.chains[p.chain_index];
        chain.vertices.contains(&6) && chain.vertices.contains(&7)
    }));
}

#[test]
fn invalid_evidence_options_are_rejected() {
    let input = unit_box_mesh();
    for bad in [
        EvidenceOptions {
            min_dihedral_angle_radians: f64::NAN,
            ..Default::default()
        },
        EvidenceOptions {
            min_dihedral_angle_radians: -1.0,
            ..Default::default()
        },
        EvidenceOptions {
            max_parallel_angle_radians: f64::INFINITY,
            ..Default::default()
        },
        EvidenceOptions {
            max_parallel_angle_radians: std::f64::consts::PI,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            analyze_mesh_evidence(&input, bad),
            Err(MeshRecoveryError::InvalidEvidenceOptions),
        ));
    }
}

fn chain_at(evidence: &MeshEvidence, a: usize, b: usize) -> usize {
    evidence
        .chains
        .iter()
        .position(|chain| {
            chain.vertices.len() == 2 && chain.vertices.contains(&a) && chain.vertices.contains(&b)
        })
        .expect("cube edge is an exact single-segment recovered chain")
}

#[test]
fn line_pair_relations_distinguish_coplanar_intersecting_and_skew_supports() {
    let cube = analyze_mesh_evidence(&unit_box_mesh(), EvidenceOptions::default()).unwrap();
    let x_front = chain_at(&cube, 0, 1);
    let x_back = chain_at(&cube, 3, 2);
    let y_front = chain_at(&cube, 0, 3);
    let z_back = chain_at(&cube, 2, 6);

    let parallel = cube.chain_pair(x_front, x_back).unwrap();
    assert!(parallel.angle_radians < 1.0e-12);
    assert!(parallel.exactly_coplanar);
    assert!(!parallel.exactly_coincident_carriers);
    assert!((parallel.supporting_line_distance - 1.0).abs() < 1.0e-12);

    let intersecting = cube.chain_pair(x_front, y_front).unwrap();
    assert!((intersecting.angle_radians - std::f64::consts::FRAC_PI_2).abs() < 1.0e-12);
    assert!(intersecting.exactly_coplanar);
    assert!(intersecting.supporting_line_distance < 1.0e-12);

    let skew = cube.chain_pair(x_front, z_back).unwrap();
    assert!((skew.angle_radians - std::f64::consts::FRAC_PI_2).abs() < 1.0e-12);
    assert!(!skew.exactly_coplanar);
    assert!((skew.supporting_line_distance - 1.0).abs() < 1.0e-12);
    assert!(cube.chain_pair(usize::MAX, 0).is_none());
}

fn assert_chain_cover_is_exact(evidence: &MeshEvidence) {
    let mut used = vec![0; evidence.edges.len()];
    for chain in &evidence.chains {
        assert_eq!(chain.vertices.len(), chain.edge_indices.len() + 1);
        for (i, &id) in chain.edge_indices.iter().enumerate() {
            used[id] += 1;
            let observed = evidence.edges[id].vertices;
            let mut expected = [chain.vertices[i], chain.vertices[i + 1]];
            expected.sort_unstable();
            assert_eq!(observed, expected);
        }
        if chain.vertices.len() > 2 {
            let first = evidence.positions[chain.vertices[0]];
            let last = evidence.positions[*chain.vertices.last().unwrap()];
            assert!(
                chain.vertices[1..chain.vertices.len() - 1]
                    .iter()
                    .all(|&id| seams::on_segment(first, last, evidence.positions[id]).is_some())
            );
        }
    }
    assert!(used.iter().all(|&count| count == 1));
}

#[test]
fn dense_planar_grid_preserves_linear_boundary_evidence_without_diagonals() {
    let n = 64;
    let stride = n + 1;
    let positions = (0..=n)
        .flat_map(|y| (0..=n).map(move |x| [crate::usize_to_f64(x), crate::usize_to_f64(y), 0.0]))
        .collect::<Vec<_>>();
    let mut triangles = Vec::with_capacity(2 * n * n);
    for y in 0..n {
        for x in 0..n {
            let a = y * stride + x;
            let b = a + 1;
            let c = a + stride;
            let d = c + 1;
            triangles.extend([[a, b, d], [a, d, c]]);
        }
    }
    let source = make_mesh(&positions, &triangles);
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.suppressed_interior_edge_count, 3 * n * n - 2 * n);
    assert_eq!(evidence.edges.len(), 4 * n);
    assert_eq!(evidence.chains.len(), 4);
    assert_eq!(evidence.parallel_bundles.len(), 2);
    assert!(
        evidence
            .chains
            .iter()
            .all(|chain| chain.edge_indices.len() == n)
    );
    assert_chain_cover_is_exact(&evidence);
}

#[test]
fn distant_world_offset_preserves_pairwise_line_relationships() {
    let mut source = unit_box_mesh();
    for point in source.positions_mut() {
        point.x += 1.0e9;
        point.y -= 1.0e9;
        point.z += 1.0e9;
    }
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.parallel_bundles.len(), 3);
    let x_front = chain_at(&evidence, 0, 1);
    let x_back = chain_at(&evidence, 3, 2);
    let relation = evidence.chain_pair(x_front, x_back).unwrap();
    assert!(relation.exactly_coplanar);
    assert!((relation.supporting_line_distance - 1.0).abs() < 1.0e-12);
    assert_chain_cover_is_exact(&evidence);
}

#[test]
fn ambiguous_multiple_body_incidence_is_preserved_without_brep_certification() {
    let source = make_mesh(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 1, 2], [1, 0, 3], [0, 1, 4]],
    );
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.non_manifold_edge_count, 1);
    assert_eq!(evidence.inconsistent_winding_edge_count, 0);
    let junction = evidence
        .edges
        .iter()
        .find(|edge| edge.kind == FeatureEdgeKind::NonManifoldJunction)
        .unwrap();
    assert_eq!(junction.vertices, [0, 1]);
    assert_eq!(junction.source_face_indices, [0, 1, 2]);
    assert_eq!(junction.source_forward_winding, [true, false, true]);
    assert_chain_cover_is_exact(&evidence);
    assert!(matches!(
        super::super::recover_planar_mesh(&source, super::super::PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::NonManifoldEdge(_))
    ));
}

#[test]
fn conflicting_mesh_windings_are_descriptive_not_implicitly_corrected() {
    let source = make_mesh(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        &[[0, 1, 2], [0, 1, 3]],
    );
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.non_manifold_edge_count, 0);
    assert_eq!(evidence.inconsistent_winding_edge_count, 1);
    let conflict = evidence
        .edges
        .iter()
        .find(|edge| edge.kind == FeatureEdgeKind::WindingConflict)
        .unwrap();
    assert_eq!(conflict.vertices, [0, 1]);
    assert_eq!(conflict.source_forward_winding, [true, true]);
    assert!(matches!(
        super::super::recover_planar_mesh(&source, super::super::PlanarRecoveryOptions::default()),
        Err(MeshRecoveryError::InconsistentWinding(_))
    ));
}

fn observed_test_edge(a: usize, b: usize) -> FeatureEdge {
    FeatureEdge {
        vertices: [a.min(b), a.max(b)],
        source_face_indices: vec![],
        source_forward_winding: vec![],
        kind: FeatureEdgeKind::Dihedral,
        dihedral_angle_radians: Some(std::f64::consts::FRAC_PI_2),
    }
}

#[test]
fn collinear_tracks_continue_through_high_degree_mesh_vertices() {
    let positions = [
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    // Two geometrically independent lines cross at a topological branch.
    let edges = [
        observed_test_edge(0, 1),
        observed_test_edge(1, 2),
        observed_test_edge(3, 1),
        observed_test_edge(1, 4),
    ];
    let chains = trace_chains(&positions, &edges);
    assert_eq!(chains.len(), 2);
    assert!(chains.iter().all(|c| c.edge_indices.len() == 2));
    let mut covered = chains
        .iter()
        .flat_map(|chain| chain.edge_indices.iter().copied())
        .collect::<Vec<_>>();
    covered.sort_unstable();
    assert_eq!(covered, [0, 1, 2, 3]);
    assert!(
        chains
            .iter()
            .any(|c| c.vertices.contains(&0) && c.vertices.contains(&2))
    );
    assert!(
        chains
            .iter()
            .any(|c| c.vertices.contains(&3) && c.vertices.contains(&4))
    );
}

#[test]
fn overlapping_collinear_ambiguity_prevents_unjustified_continuation() {
    let positions = [
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
    ];
    let edges = [
        observed_test_edge(0, 1),
        observed_test_edge(1, 2),
        observed_test_edge(1, 3),
    ];
    let chains = trace_chains(&positions, &edges);
    assert_eq!(chains.len(), 3);
    assert!(chains.iter().all(|c| c.edge_indices.len() == 1));
}
