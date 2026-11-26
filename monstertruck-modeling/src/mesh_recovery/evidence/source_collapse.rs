//! Exact zero-area source facet cancellation with oriented-chain provenance.
//! Noncollinear degeneracies are not treated as zero merely because an
//! approximate area calculation vanishes.

use super::super::{
    EdgeKey, InputPolygon, MeshRecoveryError, edge_key, seams, source_edge_incidence,
};
use crate::Point3;
use std::collections::BTreeMap;

/// An original polygon whose complete oriented boundary collapses on an
/// exactly collinear set of source vertices. Its filled area is zero.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CollapsedSourceFacet {
    /// Original input polygon ID, before any source normalization.
    pub polygon: usize,
    /// Original directed boundary walk in the initial exact-welded table.
    pub vertices: Vec<usize>,
    /// All distinct vertices ordered along one exact common line.
    /// Subdivision at these samples makes the directed edge chain zero.
    pub ordered_vertices: Vec<usize>,
}

impl CollapsedSourceFacet {
    pub(super) fn certify(
        positions: &[Point3],
        polygon: usize,
        vertices: Vec<usize>,
    ) -> Option<Self> {
        let &first = vertices.first()?;
        let mut ordered_vertices = vertices.clone();
        ordered_vertices.sort_unstable();
        ordered_vertices.dedup();
        if let Some(second) = ordered_vertices.iter().copied().find(|&id| id != first) {
            if !ordered_vertices
                .iter()
                .all(|&id| seams::on_line(positions[first], positions[second], positions[id]))
            {
                return None;
            }
            let delta = positions[second] - positions[first];
            let axis = (0..3).max_by(|&a, &b| delta[a].abs().total_cmp(&delta[b].abs()))?;
            ordered_vertices
                .sort_unstable_by(|&a, &b| positions[a][axis].total_cmp(&positions[b][axis]));
        }
        // Split the closed walk at every intermediate sample. A difference
        // array records its signed coverage on each successive interval:
        // increasing edges contribute +1 and decreasing edges contribute -1.
        let index = ordered_vertices
            .iter()
            .enumerate()
            .map(|(i, &id)| (id, i))
            .collect::<BTreeMap<_, _>>();
        let mut changes = vec![0_i64; ordered_vertices.len()];
        for (&a, &b) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
            changes[*index.get(&a)?] += 1;
            changes[*index.get(&b)?] -= 1;
        }
        let mut coverage = 0;
        if changes.into_iter().any(|change| {
            coverage += change;
            coverage != 0
        }) {
            return None;
        }
        Some(Self {
            polygon,
            vertices,
            ordered_vertices,
        })
    }
}

// A regular original source edge paired with a collapsed facet must still
// have an opposite partner over every subsegment after exact seam repair.
pub(super) struct RequiredRepair {
    active_polygon: usize,
    source_edge: EdgeKey,
    collapsed_polygon: usize,
}

pub(super) fn required_repairs(source: &[Vec<usize>], active_ids: &[usize]) -> Vec<RequiredRepair> {
    let mut active_by_source = vec![None; source.len()];
    for (active, &original) in active_ids.iter().enumerate() {
        active_by_source[original] = Some(active);
    }
    let incidence = source_edge_incidence(source.iter().map(Vec::as_slice));
    let mut repairs = Vec::new();
    for (edge, uses) in incidence {
        if let [a, b] = uses.as_slice()
            && a.1 != b.1
        {
            let pair = match (active_by_source[a.0], active_by_source[b.0]) {
                (Some(active), None) => Some((active, b.0)),
                (None, Some(active)) => Some((active, a.0)),
                _ => None,
            };
            if let Some((active_polygon, collapsed_polygon)) = pair {
                repairs.push(RequiredRepair {
                    active_polygon,
                    source_edge: edge,
                    collapsed_polygon,
                });
            }
        }
    }
    repairs.sort_unstable_by_key(|r| (r.collapsed_polygon, r.active_polygon, r.source_edge));
    repairs
}

pub(super) fn verify_repairs(
    positions: &[Point3],
    polygons: &[InputPolygon],
    repairs: &[RequiredRepair],
) -> Result<(), MeshRecoveryError> {
    if repairs.is_empty() {
        return Ok(());
    }
    let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
    for repair in repairs {
        let [start, end] = repair.source_edge.map(|id| positions[id]);
        let inside = |id: usize| {
            let p = positions[id];
            p == start || p == end || seams::on_segment(start, end, p).is_some()
        };
        let vertices = &polygons[repair.active_polygon].vertices;
        let mut observed = false;
        for (&a, &b) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
            if !inside(a) || !inside(b) {
                continue;
            }
            observed = true;
            let uses = incidence.get(&edge_key(a, b).0).map(Vec::as_slice);
            if !matches!(uses, Some([a, b]) if a.1 != b.1) {
                return Err(MeshRecoveryError::DegeneratePolygon(
                    repair.collapsed_polygon,
                ));
            }
        }
        if !observed {
            return Err(MeshRecoveryError::DegeneratePolygon(
                repair.collapsed_polygon,
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::Point3;
    use crate::mesh_recovery::{MeshRecoveryError, analyze_mesh_evidence, recover_analytic_mesh};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_topology::shell::ShellCondition;

    #[test]
    fn removes_exact_collinear_triangle_and_repairs_its_subdivided_neighbor()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(2.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                    Point3::new(1.0, 0.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([
                [0, 2, 1],
                [0, 1, 4], // zero area, joins the long edge to two short uses
                [0, 4, 3],
                [4, 1, 3],
                [0, 3, 2],
                [1, 2, 3],
            ]),
        )?;
        let recovered = recover_analytic_mesh(&mesh, Default::default())?;
        let provenance = &recovered.evidence.provenance;
        assert_eq!(provenance.original_polygon_count, 6);
        assert_eq!(provenance.original_polygon_ids, [0, 2, 3, 4, 5]);
        let [collapsed] = provenance.collapsed_facets.as_slice() else {
            return Err("one exact collapsed facet".into());
        };
        assert_eq!(collapsed.polygon, 1);
        assert_eq!(collapsed.vertices, [0, 1, 4]);
        assert_eq!(collapsed.ordered_vertices, [0, 4, 1]);
        assert_eq!(
            recovered.evidence.source_polygon_vertices(0),
            Some([0, 2, 1, 4].as_slice())
        );
        assert!(
            (recovered
                .evidence
                .closed_source_volume()
                .ok_or("closed source")?
                - 1.0 / 3.0)
                .abs()
                < 1.0e-12
        );
        let shell = recovered.complete_shell().ok_or("closed tetrahedron")?;
        assert_eq!(shell.shell_condition(), ShellCondition::Closed);
        assert_eq!(shell.len(), 4);
        Ok(())
    }

    #[test]
    fn rejects_noncollinear_zero_area_and_unrepaired_collapsed_seam()
    -> Result<(), Box<dyn std::error::Error>> {
        let attributes = StandardAttributes {
            positions: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(2.0, 1.0, 0.0),
            ],
            ..Default::default()
        };
        let missing =
            PolygonMesh::new(attributes.clone(), Faces::from_iter([[0, 3, 1], [0, 1, 2]]))?;
        assert!(matches!(
            analyze_mesh_evidence(&missing, Default::default()),
            Err(MeshRecoveryError::DegeneratePolygon(1))
        ));
        let bowtie = PolygonMesh::new(attributes, Faces::from_iter([[0, 4, 1, 3]]))?;
        assert!(matches!(
            analyze_mesh_evidence(&bowtie, Default::default()),
            Err(MeshRecoveryError::DegeneratePolygon(0))
        ));
        Ok(())
    }
}
