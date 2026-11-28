use super::*;
use crate::mesh_recovery::{NativeRecoveryOptions, recover_native_mesh};

fn collapsed_cube(reversed: bool) -> PolygonMesh {
    let points = [
        [-2., -2., 0.],
        [2., -2., 0.],
        [2., 2., 0.],
        [-2., 2., 0.],
        [-2., -2., 1.],
        [2., -2., 1.],
        [2., 2., 1.],
        [-2., 2., 1.],
        [-1., 0., 1.],
        [1., 0., 1.],
        [-1., 0., 0.9],
        [0., -1., 0.9],
        [1., 0., 0.9],
        [0., 1., 0.9],
    ];
    let mut faces = vec![
        vec![0, 3, 2, 1],
        vec![0, 1, 5, 4],
        vec![1, 2, 6, 5],
        vec![2, 3, 7, 6],
        vec![3, 0, 4, 7],
        vec![4, 5, 9],
        vec![4, 9, 8],
        vec![5, 6, 9],
        vec![6, 7, 8],
        vec![6, 8, 9],
        vec![7, 4, 8],
        vec![10, 11, 12, 13],
    ];
    let floor = [10, 11, 12, 13];
    let top = [8, 8, 9, 9];
    for ((a, b), (c, d)) in planar::edges(&floor).zip(planar::edges(&top)) {
        faces.push(vec![b, a, c]);
        if c != d {
            faces.push(vec![b, c, d]);
        }
    }
    if reversed {
        for face in &mut faces {
            face.reverse();
        }
    }
    PolygonMesh::new_unchecked(
        StandardAttributes {
            positions: points
                .into_iter()
                .map(|[x, y, z]| Point3::new(x, y, z))
                .collect(),
            ..Default::default()
        },
        Faces::from_iter(faces),
    )
}

#[test]
fn collapsed_aperture_preserves_floor_and_reopens_contact() -> Result<(), Box<dyn std::error::Error>>
{
    for reversed in [false, true] {
        let mesh = collapsed_cube(reversed);
        let source_parts = partition_recovery_mesh(&mesh)?;
        assert_eq!(source_parts.len(), 2);
        let mut candidates = collapsed_aperture_candidates(&mesh, 0.2)?;
        assert_eq!(candidates.len(), 1);
        let candidate = candidates.pop().ok_or("candidate")?;
        assert!((candidate.evidence.depth - 0.1).abs() < KERNEL_EPSILON);
        assert!(candidate.evidence.deviation_upper_bound < 0.101);
        assert_eq!(candidate.evidence.floor_ring.len(), 4);
        assert_eq!(candidate.evidence.contact_edges.len(), 1);
        let removed = candidate
            .evidence
            .cap_polygons
            .iter()
            .chain(&candidate.evidence.wall_polygons)
            .copied()
            .collect::<BTreeSet<_>>();
        let original = mesh.face_iter().collect::<Vec<_>>();
        let mut observed = BTreeSet::new();
        for (face, origin) in candidate.mesh.face_iter().zip(&candidate.polygon_origins) {
            if let AperturePolygonOrigin::Source(id) = origin {
                assert!(observed.insert(*id));
                assert_eq!(
                    face.iter()
                        .map(|v| candidate.mesh.positions()[v.pos])
                        .collect::<Vec<_>>(),
                    original[*id]
                        .iter()
                        .map(|v| mesh.positions()[v.pos])
                        .collect::<Vec<_>>()
                );
            }
        }
        assert_eq!(
            observed,
            (0..original.len())
                .filter(|id| !removed.contains(id))
                .collect()
        );
        let parts = partition_recovery_mesh(&candidate.mesh)?;
        assert_eq!(parts.len(), 1);
        assert!(parts[0].independently_closed);
        let recovery =
            recover_native_mesh(&candidate.mesh, &NativeRecoveryOptions::default(), |_| {})?;
        assert_eq!(
            recovery
                .shell
                .ok_or("closed aperture")?
                .connected_components()
                .len(),
            1
        );
    }
    Ok(())
}

#[test]
fn collapsed_aperture_requires_budget_and_unobstructed_prism()
-> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = collapsed_cube(false);
    assert!(collapsed_aperture_candidates(&mesh, 0.09)?.is_empty());
    assert!(collapsed_aperture_candidates(&mesh, f64::NAN).is_err());
    mesh.merge(PolygonMesh::new_unchecked(
        StandardAttributes {
            positions: vec![
                Point3::new(0.1, 0.1, 0.92),
                Point3::new(0.2, 0.1, 0.92),
                Point3::new(0.1, 0.2, 0.92),
                Point3::new(0.1, 0.1, 0.96),
            ],
            ..Default::default()
        },
        Faces::from_iter([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]),
    ));
    assert!(collapsed_aperture_candidates(&mesh, 0.2)?.is_empty());
    Ok(())
}

#[test]
fn planar_candidate_preserves_multiple_holes_and_concave_boundary() -> Result<(), &'static str> {
    let points = [
        [0., 0.],
        [8., 0.],
        [8., 8.],
        [5., 8.],
        [5., 4.],
        [3., 4.],
        [3., 8.],
        [0., 8.],
        [1., 1.],
        [1., 2.],
        [2., 2.],
        [2., 1.],
        [6., 5.],
        [6., 6.],
        [7., 6.],
        [7., 5.],
    ];
    let loops = vec![(0..8).collect(), (8..12).collect(), (12..16).collect()];
    let triangles = planar::triangulate(&points, &loops).ok_or("triangulation")?;
    let area = triangles
        .iter()
        .map(|p| planar::area_twice(&points, p))
        .sum::<f64>()
        * 0.5;
    assert!((area - 54.0).abs() < 1e-12);
    Ok(())
}
