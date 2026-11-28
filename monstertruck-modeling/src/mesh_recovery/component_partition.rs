//! Partition recovery work only at certified independent source components.
//! Source points and polygon cycles are preserved exactly. Components with
//! a shared exact vertex are never separated by this optimization; all
//! unresolved components remain together for the original topology repair.

use super::{
    MeshRecoveryError, component_groups::Groups, evidence::certify_independent_source_group,
    input_polygon_from_ids, welded_source_polygons,
};
use crate::{KERNEL_EPSILON, Point3};
use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
use rustc_hash::FxHashMap;
use std::{borrow::Cow, collections::BTreeMap};

/// A geometry-only recovery partition with original polygon provenance.
#[derive(Debug)]
pub struct RecoveryMeshComponent<'a> {
    /// Original mesh if no split is needed; otherwise compact geometry with
    /// unchanged positions and directed polygon cycles.
    pub mesh: Cow<'a, PolygonMesh>,
    /// Original polygon index for every local polygon, in local face order.
    pub original_polygons: Vec<usize>,
    /// Every source polygon belongs to an independently closed, consistently
    /// oriented component with a certified nonzero signed volume. This does
    /// not classify material, exclude intersections, or guarantee recovery.
    pub independently_closed: bool,
}

fn vertex_groups(vertex_count: usize, polygons: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut groups = Groups::new(polygons.len());
    let mut owner = vec![None; vertex_count];
    for (id, polygon) in polygons.iter().enumerate() {
        for &vertex in polygon {
            if let Some(first) = owner[vertex] {
                groups.join(first, id);
            } else {
                owner[vertex] = Some(id);
            }
        }
    }
    let mut members = BTreeMap::<usize, Vec<usize>>::new();
    for id in 0..polygons.len() {
        members.entry(groups.root(id)).or_default().push(id);
    }
    members.into_values().collect()
}

fn certify_group(positions: &[Point3], polygons: &[Vec<usize>], ids: &[usize]) -> bool {
    let measured = ids
        .iter()
        .map(|&id| input_polygon_from_ids(positions, polygons[id].clone(), id, KERNEL_EPSILON))
        .collect::<Result<Vec<_>, _>>();
    measured.is_ok_and(|measured| certify_independent_source_group(positions, &measured))
}

fn compact_mesh(positions: &[Point3], polygons: &[Vec<usize>], ids: &[usize]) -> PolygonMesh {
    let mut members = ids
        .iter()
        .flat_map(|&id| polygons[id].iter().copied())
        .collect::<Vec<_>>();
    members.sort_unstable();
    members.dedup();
    let local_positions = members.iter().map(|&id| positions[id]).collect();
    let remap = members
        .into_iter()
        .enumerate()
        .map(|(local, global)| (global, local))
        .collect::<FxHashMap<_, _>>();
    let faces = Faces::from_iter(ids.iter().map(|&id| {
        polygons[id]
            .iter()
            .map(|vertex| remap[vertex])
            .collect::<Vec<_>>()
    }));
    // Each index was allocated together with its position. The original
    // polygon order is already triangles, quads, then remaining n-gons;
    // a sorted subsequence therefore retains the original local mapping.
    PolygonMesh::new_unchecked(
        StandardAttributes {
            positions: local_positions,
            ..Default::default()
        },
        faces,
    )
}

/// Separate independently closed source components before global carrier
/// discovery. Exact-coordinate vertex incidence (including signed zero) is
/// shared with ordinary mesh recovery. No tolerance welding or CAD labels
/// are involved. Open, nonmanifold, degenerate and zero-volume components
/// stay in a common unresolved pool so cross-component source repair remains
/// possible. Disconnected native shells are not automatically material solids.
///
/// # Errors
/// Returns the original source-ingestion error for nonfinite coordinates or
/// invalid source position indices.
pub fn partition_recovery_mesh(
    mesh: &PolygonMesh,
) -> Result<Vec<RecoveryMeshComponent<'_>>, MeshRecoveryError> {
    let (positions, polygons) = welded_source_polygons(mesh)?;
    let mut partitions = Vec::new();
    let mut unresolved = Vec::new();
    for ids in vertex_groups(positions.len(), &polygons) {
        if certify_group(&positions, &polygons, &ids) {
            partitions.push((ids, true));
        } else {
            unresolved.extend(ids);
        }
    }
    if !unresolved.is_empty() || partitions.is_empty() {
        unresolved.sort_unstable();
        partitions.push((unresolved, false));
    }
    if partitions.len() == 1 {
        let (original_polygons, independently_closed) = partitions.remove(0);
        return Ok(vec![RecoveryMeshComponent {
            mesh: Cow::Borrowed(mesh),
            original_polygons,
            independently_closed,
        }]);
    }
    Ok(partitions
        .into_iter()
        .map(
            |(original_polygons, independently_closed)| RecoveryMeshComponent {
                mesh: Cow::Owned(compact_mesh(&positions, &polygons, &original_polygons)),
                original_polygons,
                independently_closed,
            },
        )
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Shell,
        mesh_recovery::{AnalyticRecoveryOptions, recover_analytic_mesh},
    };
    use monstertruck_topology::shell::ShellCondition;

    fn tetrahedron(offset: f64, reversed: bool) -> PolygonMesh {
        let points = [
            Point3::new(offset, 0.0, 0.0),
            Point3::new(offset + 1.0, 0.0, 0.0),
            Point3::new(offset, 1.0, 0.0),
            Point3::new(offset, 0.0, 1.0),
        ];
        let mut faces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        if reversed {
            faces.iter_mut().for_each(|face| face.reverse());
        }
        PolygonMesh::new_unchecked(
            StandardAttributes {
                positions: points.to_vec(),
                ..Default::default()
            },
            Faces::from_iter(faces),
        )
    }

    #[test]
    fn disjoint_closed_components_preserve_source_cycles_and_inward_winding()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut mesh = tetrahedron(0.0, false);
        mesh.merge(tetrahedron(3.0, true));
        let parts = partition_recovery_mesh(&mesh)?;
        assert_eq!(parts.len(), 2);
        let mut shell = Shell::new();
        let original = mesh.face_iter().collect::<Vec<_>>();
        for part in parts {
            assert!(part.independently_closed);
            assert_eq!(part.original_polygons.len(), 4);
            for (local, source) in part.mesh.face_iter().zip(&part.original_polygons) {
                let a = local
                    .iter()
                    .map(|v| part.mesh.positions()[v.pos])
                    .collect::<Vec<_>>();
                let b = original[*source]
                    .iter()
                    .map(|v| mesh.positions()[v.pos])
                    .collect::<Vec<_>>();
                assert_eq!(a, b);
            }
            let recovery = recover_analytic_mesh(&part.mesh, &AnalyticRecoveryOptions::default())?;
            shell.append(&mut recovery.complete_shell().ok_or("native component")?);
        }
        assert_eq!(shell.shell_condition(), ShellCondition::Closed);
        assert_eq!(shell.connected_components().len(), 2);
        Ok(())
    }

    #[test]
    fn touching_vertex_and_unresolved_source_groups_are_not_split()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut touching = tetrahedron(0.0, false);
        touching.merge(tetrahedron(1.0, false));
        let parts = partition_recovery_mesh(&touching)?;
        assert_eq!(parts.len(), 1);
        assert!(matches!(parts[0].mesh, Cow::Borrowed(_)));
        let mut mixed = tetrahedron(0.0, false);
        let open = |x| {
            PolygonMesh::new_unchecked(
                StandardAttributes {
                    positions: vec![
                        Point3::new(x, 0., 0.),
                        Point3::new(x + 1., 0., 0.),
                        Point3::new(x, 1., 0.),
                    ],
                    ..Default::default()
                },
                Faces::from_iter([[0, 1, 2]]),
            )
        };
        mixed.merge(open(4.));
        mixed.merge(open(8.));
        let parts = partition_recovery_mesh(&mixed)?;
        assert_eq!(parts.len(), 2);
        assert!(parts[0].independently_closed);
        assert!(!parts[1].independently_closed);
        assert_eq!(parts[1].original_polygons, vec![4, 5]);
        Ok(())
    }
}
