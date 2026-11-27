//! Globally shared source trim paths for an arbitrary polygon partition.
//! A path terminates at every region junction or ambiguous source vertex,
//! so both adjacent regions inherit identical source edge subdivisions.

use super::{RegionBoundaryEdge, RegionBoundaryError, RegionBoundaryIndex};

/// One maximal consistently oriented source contact path for a partition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionContactPath {
    /// Region on the left, and the right region if the mesh is not open.
    /// For two-sided contacts the smaller region ID is always first.
    pub regions: (usize, Option<usize>),
    /// Original source edges, oriented by the first region's polygon winding.
    /// The second region uses this exact path in reverse.
    pub edges: Vec<RegionBoundaryEdge>,
    /// True for an entire cyclic contact with no forced junction anchor.
    /// A path that returns to a junction is still anchored there.
    pub unanchored_cycle: bool,
}

/// All regular shared source trim paths in one complete region partition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionContactGraph {
    /// Paths with every regular inter-region or open-mesh edge used once.
    pub paths: Vec<RegionContactPath>,
    /// Nonmanifold or winding-conflicting source edges, never silently
    /// assigned a two-region contact. Their vertices terminate regular paths.
    pub ambiguous_edges: Vec<[usize; 2]>,
}

#[derive(Clone, Copy)]
struct ContactEdge {
    regions: (usize, Option<usize>),
    edge: RegionBoundaryEdge,
}

#[derive(Clone, Copy, Default)]
struct Junction {
    incoming: Option<usize>,
    outgoing: Option<usize>,
    blocked: bool,
}

impl RegionBoundaryIndex<'_> {
    /// Extract shared trim paths for any complete source-polygon partition.
    ///
    /// Region IDs need not be consecutive. Every original edge is considered,
    /// including tessellation diagonals that separate assigned regions.
    /// A path continues only at a vertex with one incoming and one outgoing
    /// edge, matching region identities, and no ambiguous source incidence.
    /// This is source topology, not a smooth analytic-curve classification.
    ///
    /// # Errors
    /// Returns an error if the partition does not assign every original
    /// polygon or if the stored source evidence is inconsistent.
    pub fn partition_contacts(
        &self,
        region_by_polygon: &[usize],
    ) -> Result<RegionContactGraph, RegionBoundaryError> {
        if region_by_polygon.len() != self.evidence.source_planar_support.len() {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let mut keys = self.incidence.keys().copied().collect::<Vec<_>>();
        keys.sort_unstable();
        let mut contacts = Vec::new();
        let mut ambiguous_edges = Vec::new();
        let mut junctions = vec![Junction::default(); self.evidence.positions.len()];
        for key in keys {
            let uses = self
                .incidence
                .get(&key)
                .ok_or(RegionBoundaryError::InvalidEvidence)?;
            let (first, second) = match uses.as_slice() {
                [a] => (*a, None),
                [a, b] if a.1 != b.1 => {
                    let left = region_by_polygon[a.0];
                    let right = region_by_polygon[b.0];
                    if left == right {
                        continue;
                    }
                    if left < right {
                        (*a, Some(*b))
                    } else {
                        (*b, Some(*a))
                    }
                }
                _ => {
                    ambiguous_edges.push(key);
                    for vertex in key {
                        junctions
                            .get_mut(vertex)
                            .ok_or(RegionBoundaryError::InvalidEvidence)?
                            .blocked = true;
                    }
                    continue;
                }
            };
            let vertices = if first.1 { key } else { [key[1], key[0]] };
            let other_polygon = second.map(|other| other.0);
            let id = contacts.len();
            let incoming = junctions
                .get_mut(vertices[1])
                .ok_or(RegionBoundaryError::InvalidEvidence)?;
            incoming.blocked |= incoming.incoming.replace(id).is_some();
            let outgoing = junctions
                .get_mut(vertices[0])
                .ok_or(RegionBoundaryError::InvalidEvidence)?;
            outgoing.blocked |= outgoing.outgoing.replace(id).is_some();
            contacts.push(ContactEdge {
                regions: (
                    region_by_polygon[first.0],
                    other_polygon.map(|polygon| region_by_polygon[polygon]),
                ),
                edge: RegionBoundaryEdge {
                    vertices,
                    interior_polygon: first.0,
                    exterior_polygon: other_polygon,
                },
            });
        }
        Ok(RegionContactGraph {
            paths: trace_contact_paths(&contacts, &junctions),
            ambiguous_edges,
        })
    }
}

fn trace_contact_paths(contacts: &[ContactEdge], junctions: &[Junction]) -> Vec<RegionContactPath> {
    let mut next = vec![None; contacts.len()];
    let mut previous = vec![None; contacts.len()];
    for junction in junctions {
        if junction.blocked {
            continue;
        }
        if let (Some(incoming), Some(outgoing)) = (junction.incoming, junction.outgoing)
            && contacts[incoming].regions == contacts[outgoing].regions
        {
            next[incoming] = Some(outgoing);
            previous[outgoing] = Some(incoming);
        }
    }
    let mut visited = vec![false; contacts.len()];
    let mut paths = Vec::new();
    // First walk anchored paths, then cycles without junctions. The original
    // canonical source edge order makes both deterministic.
    let anchored = (0..contacts.len()).filter(|&id| previous[id].is_none());
    for start in anchored.chain(0..contacts.len()) {
        if visited[start] {
            continue;
        }
        let mut edges = Vec::new();
        let mut current = start;
        let unanchored_cycle = loop {
            visited[current] = true;
            edges.push(contacts[current].edge);
            match next[current] {
                Some(successor) if !visited[successor] => current = successor,
                Some(successor) => break successor == start,
                None => break false,
            }
        };
        paths.push(RegionContactPath {
            regions: contacts[start].regions,
            edges,
            unanchored_cycle,
        });
    }
    paths
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Point3,
        mesh_recovery::{EvidenceOptions, analyze_mesh_evidence},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    #[test]
    fn shared_region_paths_stop_at_partition_junctions() -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [1, 0, 3]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let graph = index.partition_contacts(&[7, 2])?;
        assert!(graph.ambiguous_edges.is_empty());
        assert_eq!(graph.paths.len(), 3);
        assert_eq!(graph.paths.iter().map(|p| p.edges.len()).sum::<usize>(), 5);
        let shared = graph
            .paths
            .iter()
            .find(|p| p.regions == (2, Some(7)))
            .ok_or("missing contact")?;
        assert_eq!(shared.edges.len(), 1);
        assert_eq!(shared.edges[0].interior_polygon, 1);
        assert_eq!(shared.edges[0].exterior_polygon, Some(0));
        assert!(graph.paths.iter().all(|p| !p.unanchored_cycle));
        let joined = index.partition_contacts(&[7, 7])?;
        assert_eq!(joined.paths.len(), 1);
        assert_eq!(joined.paths[0].edges.len(), 4);
        assert!(joined.paths[0].unanchored_cycle);
        Ok(())
    }

    #[test]
    fn ambiguous_incidence_is_not_a_shared_trim() -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                    Point3::new(0.0, -1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [1, 0, 3], [0, 1, 4]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let graph = index.partition_contacts(&[0, 1, 2])?;
        assert_eq!(graph.ambiguous_edges.len(), 1);
        assert_eq!(graph.paths.len(), 3);
        assert!(
            graph
                .paths
                .iter()
                .all(|p| p.edges.len() == 2 && !p.unanchored_cycle && p.regions.1.is_none())
        );
        assert_eq!(
            index.partition_contacts(&[0, 1]),
            Err(RegionBoundaryError::InvalidPartition)
        );
        Ok(())
    }
}
