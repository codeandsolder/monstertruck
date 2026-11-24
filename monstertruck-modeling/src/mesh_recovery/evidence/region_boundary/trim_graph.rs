//! Assemble shared trim geometry into the original oriented region loops.
//! No vertex projection, periodic chart cut, or independent face-edge fit
//! is performed here. Every regular source contact retains one owner.

use super::{
    RegionBoundaryEdge, RegionBoundaryError, RegionBoundaryIndex, SharedTrimGeometry,
    SharedTrimOptions, SurfaceRegionPartition,
};
use rustc_hash::FxHashMap as HashMap;

/// One global trim with its source provenance, shared by both incident faces.
#[derive(Clone, Debug)]
pub struct SharedAnalyticTrim {
    /// Index of the original contact path.
    pub contact_path: usize,
    /// Region IDs in the canonical contact orientation.
    pub regions: (usize, Option<usize>),
    /// Original edges in the trim's directed order.
    pub source_edges: Vec<RegionBoundaryEdge>,
    /// Independently checked geometry, or an explicitly unresolved path.
    pub geometry: SharedTrimGeometry,
}

/// An oriented reference into the global trim table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrientedTrim {
    /// Index of the unique shared trim.
    pub trim: usize,
    /// Whether the region traverses the trim in reverse order.
    pub reversed: bool,
}

/// Winding-preserving loops of shared trim references.
#[derive(Debug)]
pub struct RegionTrimLoops {
    /// Cyclic sequences of globally shared edges.
    pub loops: Vec<Vec<OrientedTrim>>,
    /// Original boundary-sector pinches, retained for later chart planning.
    pub pinched_vertices: Vec<usize>,
}

/// The shared trim boundary of one selected or fallback region.
#[derive(Debug)]
pub struct RegionTrimBoundary {
    /// Active region ID in the surface partition.
    pub region: usize,
    /// Failed source boundaries are represented rather than hidden.
    pub boundary: Result<RegionTrimLoops, RegionBoundaryError>,
}

/// Topological trim graph before geometric vertex reconciliation and charts.
#[derive(Debug)]
pub struct SharedTrimGraph {
    /// Every regular inter-region or open boundary source edge belongs to
    /// exactly one entry, even if analytic fitting was unsuccessful.
    pub trims: Vec<SharedAnalyticTrim>,
    /// Every active region, including boundary-free and failed regions.
    pub regions: Vec<RegionTrimBoundary>,
    /// Original nonmanifold or winding-conflicting source edges.
    pub ambiguous_edges: Vec<[usize; 2]>,
    /// Contact paths whose mixed-curve search was restricted.
    pub limited_paths: usize,
}

#[derive(Clone, Copy)]
struct EdgeLocation {
    trim: usize,
    offset: usize,
}

struct TrimTable {
    trims: Vec<SharedAnalyticTrim>,
    locations: HashMap<[usize; 2], EdgeLocation>,
    limited_paths: usize,
}

impl RegionBoundaryIndex<'_> {
    /// Assemble shared trims and region loops from a complete partition.
    ///
    /// A region with ambiguous incidence or a noncyclic source boundary
    /// retains its individual error. For every valid region, each trim use
    /// consumes its full original directed edge sequence exactly once.
    ///
    /// # Errors
    /// Inconsistent partition membership or contact data returns an error.
    /// Invalid geometric fit options leave the original contacts unresolved.
    pub fn assemble_shared_trims(
        &self,
        partition: &SurfaceRegionPartition,
        options: SharedTrimOptions,
    ) -> Result<SharedTrimGraph, RegionBoundaryError> {
        if partition.region_by_support.len() != self.evidence.planar_supports.len()
            || self.evidence.source_planar_support.len() != partition.region_by_polygon.len()
            || self
                .evidence
                .source_planar_support
                .iter()
                .zip(&partition.region_by_polygon)
                .any(|(&support, &region)| {
                    partition.region_by_support.get(support) != Some(&region)
                })
            || self.partition_contacts(&partition.region_by_polygon)? != partition.contacts
        {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let mut groups = vec![Vec::new(); partition.carriers.len()];
        for (support, &region) in partition.region_by_support.iter().enumerate() {
            if partition
                .carriers
                .get(region)
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(RegionBoundaryError::InvalidPartition);
            }
            groups[region].push(support);
        }
        if groups
            .iter()
            .zip(&partition.carriers)
            .any(|(group, carrier)| group.is_empty() == carrier.is_some())
        {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let TrimTable {
            trims,
            locations,
            limited_paths,
        } = self.shared_trim_table(partition, options)?;
        let mut regions = Vec::new();
        for (region, supports) in groups
            .iter()
            .enumerate()
            .filter(|(_, group)| !group.is_empty())
        {
            let boundary = self.boundary(supports).and_then(|boundary| {
                let loops = boundary
                    .loops
                    .iter()
                    .map(|cycle| shared_cycle(cycle, region, &trims, &locations))
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(RegionTrimLoops {
                    loops,
                    pinched_vertices: boundary.pinched_vertices,
                })
            });
            regions.push(RegionTrimBoundary { region, boundary });
        }
        validate_uses(&trims, &regions)?;
        Ok(SharedTrimGraph {
            trims,
            regions,
            ambiguous_edges: partition.contacts.ambiguous_edges.clone(),
            limited_paths,
        })
    }

    fn shared_trim_table(
        &self,
        partition: &SurfaceRegionPartition,
        options: SharedTrimOptions,
    ) -> Result<TrimTable, RegionBoundaryError> {
        let mut trims = Vec::new();
        let mut locations = HashMap::default();
        let mut limited_paths = 0;
        for (contact_path, path) in partition.contacts.paths.iter().enumerate() {
            let sides = std::iter::once(path.regions.0)
                .chain(path.regions.1)
                .map(|id| {
                    partition
                        .carriers
                        .get(id)
                        .copied()
                        .flatten()
                        .ok_or(RegionBoundaryError::InvalidPartition)
                })
                .collect::<Result<Vec<_>, _>>()?;
            let segmentation =
                path.segment_analytic_trims(&self.evidence.positions, &sides, options);
            let pieces = if let Some(segmentation) = segmentation {
                limited_paths += usize::from(segmentation.limited);
                segmentation.pieces
            } else {
                vec![super::SharedTrimPiece {
                    first_edge: 0,
                    edge_count: path.edges.len(),
                    geometry: SharedTrimGeometry::Unresolved,
                }]
            };
            let mut covered = vec![false; path.edges.len()];
            for piece in pieces {
                if piece.edge_count == 0
                    || piece.edge_count > path.edges.len()
                    || piece.first_edge >= path.edges.len()
                {
                    return Err(RegionBoundaryError::InvalidEvidence);
                }
                let trim = trims.len();
                let mut source_edges = Vec::with_capacity(piece.edge_count);
                for offset in 0..piece.edge_count {
                    let id = (piece.first_edge + offset) % path.edges.len();
                    if std::mem::replace(&mut covered[id], true) {
                        return Err(RegionBoundaryError::InvalidEvidence);
                    }
                    let edge = path.edges[id];
                    let key = super::edge_key(edge.vertices[0], edge.vertices[1]).0;
                    if locations
                        .insert(key, EdgeLocation { trim, offset })
                        .is_some()
                    {
                        return Err(RegionBoundaryError::InvalidEvidence);
                    }
                    source_edges.push(edge);
                }
                trims.push(SharedAnalyticTrim {
                    contact_path,
                    regions: path.regions,
                    source_edges,
                    geometry: piece.geometry,
                });
            }
            if covered.contains(&false) {
                return Err(RegionBoundaryError::InvalidEvidence);
            }
        }
        Ok(TrimTable {
            trims,
            locations,
            limited_paths,
        })
    }
}

fn locate(
    edge: &RegionBoundaryEdge,
    region: usize,
    trims: &[SharedAnalyticTrim],
    locations: &HashMap<[usize; 2], EdgeLocation>,
) -> Result<(EdgeLocation, bool), RegionBoundaryError> {
    let key = super::edge_key(edge.vertices[0], edge.vertices[1]).0;
    let &location = locations
        .get(&key)
        .ok_or(RegionBoundaryError::InvalidEvidence)?;
    let trim = &trims[location.trim];
    let source = trim.source_edges[location.offset];
    if trim.regions.0 == region && source == *edge {
        Ok((location, false))
    } else if trim.regions.1 == Some(region)
        && source.vertices == [edge.vertices[1], edge.vertices[0]]
        && source.exterior_polygon == Some(edge.interior_polygon)
        && edge.exterior_polygon == Some(source.interior_polygon)
    {
        Ok((location, true))
    } else {
        Err(RegionBoundaryError::InvalidEvidence)
    }
}

fn shared_cycle(
    cycle: &[RegionBoundaryEdge],
    region: usize,
    trims: &[SharedAnalyticTrim],
    locations: &HashMap<[usize; 2], EdgeLocation>,
) -> Result<Vec<OrientedTrim>, RegionBoundaryError> {
    let ordered = cycle
        .iter()
        .map(|edge| locate(edge, region, trims, locations))
        .collect::<Result<Vec<_>, _>>()?;
    let start = ordered
        .iter()
        .position(|(location, reversed)| {
            location.offset
                == if *reversed {
                    trims[location.trim].source_edges.len() - 1
                } else {
                    0
                }
        })
        .ok_or(RegionBoundaryError::NonCyclicBoundary)?;
    let mut result = Vec::new();
    let mut consumed = 0;
    while consumed < ordered.len() {
        let (first, reversed) = ordered[(start + consumed) % ordered.len()];
        let length = trims[first.trim].source_edges.len();
        if length > ordered.len() - consumed {
            return Err(RegionBoundaryError::NonCyclicBoundary);
        }
        for i in 0..length {
            let (location, direction) = ordered[(start + consumed + i) % ordered.len()];
            if location.trim != first.trim
                || direction != reversed
                || location.offset != if reversed { length - 1 - i } else { i }
            {
                return Err(RegionBoundaryError::NonCyclicBoundary);
            }
        }
        result.push(OrientedTrim {
            trim: first.trim,
            reversed,
        });
        consumed += length;
    }
    Ok(result)
}

fn validate_uses(
    trims: &[SharedAnalyticTrim],
    regions: &[RegionTrimBoundary],
) -> Result<(), RegionBoundaryError> {
    let mut uses = vec![[0_usize; 2]; trims.len()];
    let valid = regions
        .iter()
        .filter(|region| region.boundary.is_ok())
        .map(|region| region.region)
        .collect::<rustc_hash::FxHashSet<_>>();
    for region in regions {
        if let Ok(boundary) = &region.boundary {
            for edge in boundary.loops.iter().flatten() {
                uses[edge.trim][usize::from(edge.reversed)] += 1;
            }
        }
    }
    for (trim, count) in trims.iter().zip(uses) {
        let expected = [
            usize::from(valid.contains(&trim.regions.0)),
            usize::from(trim.regions.1.is_some_and(|region| valid.contains(&region))),
        ];
        if count != expected {
            return Err(RegionBoundaryError::InvalidEvidence);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Point3,
        mesh_recovery::{EvidenceOptions, analyze_mesh_evidence},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    fn options() -> SharedTrimOptions {
        SharedTrimOptions {
            absolute_tolerance: 1.0e-8,
            relative_tolerance: 1.0e-4,
            max_interval_edges: 128,
        }
    }

    #[test]
    fn adjacent_faces_reference_the_same_trim_in_opposite_directions()
    -> Result<(), Box<dyn std::error::Error>> {
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
        let partition = index.select_surface_partition(&[], &evidence, 100)?;
        let graph = index.assemble_shared_trims(&partition, options())?;
        assert_eq!(graph.trims.len(), 5);
        assert_eq!(graph.regions.len(), 2);
        let shared = graph
            .trims
            .iter()
            .position(|edge| edge.regions.1.is_some())
            .ok_or("missing shared edge")?;
        let first = graph.regions[0]
            .boundary
            .as_ref()
            .map_err(|error| *error)?
            .loops[0]
            .iter()
            .find(|edge| edge.trim == shared)
            .ok_or("missing first use")?;
        let second = graph.regions[1]
            .boundary
            .as_ref()
            .map_err(|error| *error)?
            .loops[0]
            .iter()
            .find(|edge| edge.trim == shared)
            .ok_or("missing second use")?;
        assert_ne!(first.reversed, second.reversed);
        Ok(())
    }

    #[test]
    fn closed_planar_cycle_retains_the_complete_boundary() -> Result<(), Box<dyn std::error::Error>>
    {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(2.0, 0.0, 0.0),
                    Point3::new(2.0, 1.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 4], [1, 3, 4], [1, 2, 3]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let partition = index.select_surface_partition(&[], &evidence, 100)?;
        let graph = index.assemble_shared_trims(&partition, options())?;
        assert_eq!(graph.trims.len(), 4);
        assert_eq!(
            graph
                .trims
                .iter()
                .map(|trim| trim.source_edges.len())
                .sum::<usize>(),
            5
        );
        assert_eq!(
            graph.regions[0]
                .boundary
                .as_ref()
                .map_err(|error| *error)?
                .loops[0]
                .len(),
            4
        );
        Ok(())
    }
}
