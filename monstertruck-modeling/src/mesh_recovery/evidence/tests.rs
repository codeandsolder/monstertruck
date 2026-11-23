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

#[test]
fn connected_planar_supports_recover_cube_planes_with_measured_area() {
    let source = unit_box_mesh();
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.source_planar_support.len(), 12);
    assert_eq!(evidence.planar_supports.len(), 6);
    for support in &evidence.planar_supports {
        assert_eq!(support.source_face_indices.len(), 2);
        assert!((support.area - 1.0).abs() < 1.0e-12);
        assert!(support.max_vertex_plane_residual < 1.0e-12);
        assert!(support.max_normal_angle_radians < 1.0e-12);
    }
    for (index, _) in evidence.edges.iter().enumerate() {
        let [a, b] = evidence.edge_planar_support_pair(index).unwrap();
        assert_ne!(a, b);
    }
}

#[test]
fn planar_supports_do_not_depend_on_dihedral_feature_selection() {
    let source = unit_box_mesh();
    let full = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    let fewer_features = analyze_mesh_evidence(
        &source,
        EvidenceOptions {
            min_dihedral_angle_radians: 1.75,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        full.planar_supports.len(),
        fewer_features.planar_supports.len()
    );
    assert_eq!(
        full.source_planar_support,
        fewer_features.source_planar_support
    );
    assert_eq!(fewer_features.edges.len(), 0);
}

#[test]
fn partial_or_ambiguous_incidence_does_not_fuse_planar_supports() {
    let source = make_mesh(
        &[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
            [1.0, 0.0, 1.0],
        ],
        &[[0, 1, 2], [1, 0, 3], [0, 1, 4]],
    );
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.non_manifold_edge_count, 1);
    assert_eq!(evidence.planar_supports.len(), 3);
    assert_eq!(evidence.source_planar_support.len(), 3);
    assert!(evidence.edges.iter().enumerate().any(|(id, edge)| {
        edge.kind == FeatureEdgeKind::NonManifoldJunction
            && evidence.edge_planar_support_pair(id).is_none()
    }));
}

#[test]
fn rejects_invalid_independent_planar_support_bounds() {
    let mesh = unit_box_mesh();
    for options in [
        PlanarRecoveryOptions {
            distance_tolerance: f64::NAN,
            ..Default::default()
        },
        PlanarRecoveryOptions {
            angular_tolerance_radians: f64::INFINITY,
            ..Default::default()
        },
        PlanarRecoveryOptions {
            distance_tolerance: -1.0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            analyze_mesh_evidence(
                &mesh,
                EvidenceOptions {
                    planar_support: options,
                    ..Default::default()
                }
            ),
            Err(MeshRecoveryError::InvalidEvidenceOptions)
        ));
    }
}

#[test]
fn concave_polygon_support_uses_oriented_area_not_absolute_fan_triangles() {
    let coordinates = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 2.0, 0.0],
        [1.0, 0.5, 0.0],
        [0.0, 2.0, 0.0],
    ]
    .map(|[x, y, z]| Point3::new(x, y, z))
    .to_vec();
    let source = PolygonMesh::new(
        StandardAttributes {
            positions: coordinates,
            ..Default::default()
        },
        Faces::from_iter([vec![0, 1, 2, 3, 4]]),
    )
    .unwrap();
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.source_planar_support, [0]);
    assert_eq!(evidence.planar_supports.len(), 1);
    assert!((evidence.planar_supports[0].area - 2.5).abs() < 1.0e-12);
    assert_eq!(evidence.edges.len(), 5);
}

#[test]
fn cube_support_interfaces_form_exact_original_shared_edges() {
    let mesh = unit_box_mesh();
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.planar_supports.len(), 6);
    assert_eq!(evidence.support_interfaces.len(), 12);
    assert_eq!(evidence.support_interface_indices.len(), 6);
    let mut seen = std::collections::HashSet::new();
    for interface in &evidence.support_interfaces {
        assert!(interface.supports[0] < interface.supports[1]);
        assert_eq!(interface.topology, SupportBoundaryTopology::OpenPath);
        assert_eq!(interface.source_edges.len(), 1);
        assert_eq!(interface.source_polygon_pairs.len(), 1);
        assert!(seen.insert(interface.source_edges[0]));
        assert!((interface.boundary_length - 1.0).abs() < 1.0e-12);
        assert!((interface.mean_turn_radians - std::f64::consts::FRAC_PI_2).abs() < 1.0e-12);
        assert!(interface.max_seed_line_offset < 1.0e-12);
        assert!(interface.max_seed_direction_spread_radians < 1.0e-12);
        for (region, &polygon) in interface
            .supports
            .iter()
            .zip(&interface.source_polygon_pairs[0])
        {
            assert_eq!(evidence.source_planar_support[polygon], *region);
        }
    }
    assert!(
        evidence
            .support_interface_indices
            .iter()
            .all(|ids| ids.len() == 4)
    );
}

#[test]
fn ambiguous_source_edge_incidence_does_not_produce_fake_interfaces() {
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
    assert!(evidence.support_interfaces.is_empty());
    assert!(evidence.support_interface_indices.iter().all(Vec::is_empty));
}

#[test]
fn cylinder_wall_reveals_repeated_coherent_normal_rotation_without_primitive_fitting() {
    let source = monstertruck_mesh::stl::read(
        &include_bytes!(
            "../../../../monstertruck-io/tests/fixtures/mesh-recovery/through-hole-32.stl"
        )[..],
        monstertruck_mesh::stl::StlType::Automatic,
    )
    .unwrap();
    let evidence = analyze_mesh_evidence(&source, EvidenceOptions::default()).unwrap();
    let mut coherent_axial = 0;
    for (center, interfaces) in evidence.support_interface_indices.iter().enumerate() {
        for (i, &first) in interfaces.iter().enumerate() {
            for &second in &interfaces[i + 1..] {
                let Some(observation) = evidence.support_normal_rotation(center, first, second)
                else {
                    continue;
                };
                let z = observation.reference_direction.dot(Vector3::unit_z()).abs();
                if z < 0.99999 || !observation.coherent_rotation {
                    continue;
                }
                if observation.interface_axis_residual_radians > 1.0e-4
                    || observation.max_normal_axis_dot > 1.0e-4
                    || observation.turn_magnitude_mismatch_radians > 1.0e-4
                {
                    continue;
                }
                let expected_turn = std::f64::consts::TAU / 32.0;
                if (observation.turn_radians[0] - expected_turn).abs() < 1.0e-3
                    && (observation.turn_radians[1] - expected_turn).abs() < 1.0e-3
                {
                    coherent_axial += 1;
                }
            }
        }
    }
    assert!(
        coherent_axial >= 24,
        "found {coherent_axial} coherent wall strips"
    );
}

#[test]
fn normal_rotation_requires_distinct_center_incident_interfaces() {
    let evidence = analyze_mesh_evidence(&unit_box_mesh(), EvidenceOptions::default()).unwrap();
    assert!(evidence.support_normal_rotation(0, 0, 0).is_none());
    assert!(evidence.support_normal_rotation(usize::MAX, 0, 1).is_none());
    assert!(evidence.support_normal_rotation(0, usize::MAX, 1).is_none());
    for (index, interface) in evidence.support_interfaces.iter().enumerate() {
        if !interface.supports.contains(&0) {
            assert!(evidence.support_normal_rotation(0, index, 0).is_none());
        }
    }
}

fn faceted_conical_wall(
    bottom_radius: f64,
    top_radius: f64,
    sectors: usize,
    xy_scale: [f64; 2],
    shift: [f64; 3],
) -> PolygonMesh {
    let mut vertices = Vec::with_capacity(2 * sectors);
    for z in [0.0, 2.0] {
        let radius = if z == 0.0 { bottom_radius } else { top_radius };
        for index in 0..sectors {
            let phi = std::f64::consts::TAU * index as f64 / sectors as f64;
            vertices.push([
                radius.mul_add(xy_scale[0] * phi.cos(), shift[0]),
                radius.mul_add(xy_scale[1] * phi.sin(), shift[1]),
                z + shift[2],
            ]);
        }
    }
    let mut triangles = Vec::with_capacity(2 * sectors);
    for index in 0..sectors {
        let next = (index + 1) % sectors;
        triangles.push([index, next, sectors + next]);
        triangles.push([index, sectors + next, sectors + index]);
    }
    make_mesh(&vertices, &triangles)
}

#[test]
fn global_curvature_tracks_recover_full_cylindrical_normal_ring() {
    let mesh = faceted_conical_wall(1.0, 1.0, 32, [1.0, 1.0], [0.0; 3]);
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default()).unwrap();
    let tracks = trace_curvature_tracks(&evidence, CurvatureTrackOptions::default()).unwrap();
    let ring = tracks
        .iter()
        .find(|track| {
            track.support_indices.len() >= 30
                && track.seed_axis.dot(Vector3::unit_z()).abs() > 0.9999
                && track.normal_axis_projection_spread < 1.0e-6
        })
        .unwrap_or_else(|| panic!("no full cylindrical normal ring in {} tracks", tracks.len()));
    assert!(ring.closed);
    assert_eq!(ring.distinct_support_count, 32);
    assert_eq!(ring.support_indices.first(), ring.support_indices.last());
    assert_eq!(ring.interface_indices.len() + 1, ring.support_indices.len());
    assert!((ring.accumulated_normal_turn_radians - std::f64::consts::TAU).abs() < 1.0e-4);
    assert!(ring.mean_normal_axis_projection.abs() < 1.0e-6);
    assert!(ring.max_seed_axis_error_radians < 1.0e-4);
}

#[test]
fn normal_space_axis_tracks_conical_surfaces_without_cylinder_assumption() {
    let mesh = faceted_conical_wall(1.0, 0.5, 24, [1.0, 1.0], [0.0; 3]);
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default()).unwrap();
    let tracks = trace_curvature_tracks(&evidence, CurvatureTrackOptions::default()).unwrap();
    let ring = tracks
        .iter()
        .find(|track| {
            track.support_indices.len() >= 24
                && track.seed_axis.dot(Vector3::unit_z()).abs() > 0.9999
                && track.normal_axis_projection_spread < 1.0e-5
        })
        .unwrap_or_else(|| panic!("no conical ring in {} tracks", tracks.len()));
    assert!(ring.closed);
    assert_eq!(ring.distinct_support_count, 24);
    assert!(ring.mean_normal_axis_projection.abs() > 0.05);
    assert!(ring.max_seed_axis_error_radians < 1.0e-4);
}

#[test]
fn normal_space_tracks_validate_options_and_do_not_fuse_ambiguous_source_contacts() {
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
    assert!(
        trace_curvature_tracks(&evidence, CurvatureTrackOptions::default())
            .unwrap()
            .is_empty()
    );
    for bad_options in [
        CurvatureTrackOptions {
            min_normal_turn_radians: f64::NAN,
            ..Default::default()
        },
        CurvatureTrackOptions {
            min_normal_axis_strength: -1.0,
            ..Default::default()
        },
        CurvatureTrackOptions {
            max_seed_axis_error_radians: std::f64::consts::PI,
            ..Default::default()
        },
        CurvatureTrackOptions {
            max_normal_axis_projection_spread: f64::INFINITY,
            ..Default::default()
        },
        CurvatureTrackOptions {
            min_distinct_supports: 2,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            trace_curvature_tracks(&evidence, bad_options),
            Err(MeshRecoveryError::InvalidEvidenceOptions)
        ));
    }
}

#[test]
fn original_source_polygon_positions_are_available_without_mesh_reparsing() {
    let evidence = analyze_mesh_evidence(&unit_box_mesh(), EvidenceOptions::default()).unwrap();
    assert_eq!(evidence.source_polygon_vertex_offsets.len(), 13);
    assert_eq!(evidence.source_polygon_vertex_indices.len(), 36);
    assert_eq!(evidence.source_polygon_vertices(0).unwrap().len(), 3);
    assert_eq!(evidence.source_polygon_vertices(11).unwrap().len(), 3);
    assert!(evidence.source_polygon_vertices(12).is_none());
    assert!(evidence.source_polygon_vertices(usize::MAX).is_none());
}

#[test]
fn spatial_evidence_fits_shifted_cylinder_using_both_source_planes_and_vertices() {
    let shift = [1.0e7, -2.0e7, 3.0e7];
    let mesh = faceted_conical_wall(1.0, 1.0, 32, [1.0, 1.0], shift);
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default()).unwrap();
    let tracks = trace_curvature_tracks(&evidence, CurvatureTrackOptions::default()).unwrap();
    let ring = tracks
        .iter()
        .find(|track| track.closed && track.distinct_support_count == 32)
        .unwrap();
    let spatial = measure_curvature_spatial(&evidence, ring).unwrap();
    assert_eq!(spatial.unique_source_vertex_count, 64);
    assert_eq!(spatial.axial_stations.len(), 2);
    assert!(spatial.max_meridian_slope_change.is_none());
    assert!(spatial.transverse_normal_rank_ratio > 0.9);
    assert!((spatial.axis_origin.x - shift[0]).abs() < 1.0e-7);
    assert!((spatial.axis_origin.y - shift[1]).abs() < 1.0e-7);
    assert!((spatial.mean_sample_radius - 1.0).abs() < 1.0e-7);
    assert!(spatial.max_support_plane_offset_residual < 1.0e-7);
    assert!(spatial.max_sample_radial_residual < 1.0e-7);
    assert!(spatial.relative_max_sample_radial_residual < 1.0e-7);
    assert!(spatial.fitted_sample_radius_slope.unwrap().abs() < 1.0e-7);
    assert!(spatial.normal_implied_radius_slope.unwrap().abs() < 1.0e-7);
    let apothem = (std::f64::consts::PI / 32.0).cos();
    assert!((spatial.signed_support_plane_offset.abs() - apothem).abs() < 1.0e-7);
}

#[test]
fn spatial_evidence_fits_cone_slope_independently_from_observed_normals() {
    let mesh = faceted_conical_wall(1.0, 0.5, 24, [1.0, 1.0], [12.0, -7.0, 3.0]);
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default()).unwrap();
    let tracks = trace_curvature_tracks(&evidence, CurvatureTrackOptions::default()).unwrap();
    let ring = tracks
        .iter()
        .find(|track| track.closed && track.distinct_support_count == 24)
        .unwrap();
    let spatial = measure_curvature_spatial(&evidence, ring).unwrap();
    assert_eq!(spatial.unique_source_vertex_count, 48);
    assert_eq!(spatial.axial_stations.len(), 2);
    assert!(spatial.max_meridian_slope_change.is_none());
    assert!(spatial.transverse_normal_rank_ratio > 0.9);
    assert!(spatial.max_support_plane_offset_residual < 1.0e-9);
    assert!(spatial.max_sample_radial_residual < 1.0e-9);
    assert!((spatial.fitted_sample_radius_slope.unwrap().abs() - 0.25).abs() < 1.0e-9);
    // The facet plane describes the apothem, not the original circumradius.
    // Their expected slopes differ by cos(half a polygon segment angle);
    // this is a measurable faceting effect, not a bad conical fit.
    let polygon_apothem_factor = (std::f64::consts::PI / 24.0).cos();
    assert!(
        (spatial.normal_implied_radius_slope.unwrap().abs() - 0.25 * polygon_apothem_factor).abs()
            < 1.0e-9
    );
    assert!(
        (spatial.radius_slope_disagreement.unwrap() - 0.25 * (1.0 - polygon_apothem_factor)).abs()
            < 1.0e-9
    );
}

#[test]
fn spatial_evidence_refutes_axis_coherent_elliptical_extrusion_as_a_cylinder() {
    // The normals still turn monotonically through 360 degrees around Z,
    // but no cylinder axis and fixed radius can explain the source vertices.
    let mesh = faceted_conical_wall(1.0, 1.0, 48, [1.6, 1.0], [0.0; 3]);
    let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default()).unwrap();
    let tracks = trace_curvature_tracks(&evidence, CurvatureTrackOptions::default()).unwrap();
    let ring = tracks
        .iter()
        .find(|track| track.closed && track.distinct_support_count == 48)
        .unwrap();
    assert!(ring.max_seed_axis_error_radians < 1.0e-8);
    let spatial = measure_curvature_spatial(&evidence, ring).unwrap();
    assert_eq!(spatial.unique_source_vertex_count, 96);
    assert_eq!(spatial.axial_stations.len(), 2);
    assert!(spatial.max_meridian_slope_change.is_none());
    assert!(spatial.transverse_normal_rank_ratio > 0.1);
    assert!(spatial.relative_max_sample_radial_residual > 0.1);
    assert!(spatial.max_support_plane_offset_residual > 0.1);
}
