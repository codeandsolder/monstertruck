//! Topological candidate expansion across original source vertex fans.
//!
//! Connect selected sectors through explicitly permitted regions. This is
//! only a proposal: the expanded patch still needs independent junction,
//! approximation and native-shell validation.

use super::{RegionBoundaryError, RegionBoundaryIndex, SurfaceRegionPartition};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

struct RegionIncidence {
    vertex_polygons: Vec<Vec<usize>>,
    region_vertices: Vec<Vec<usize>>,
    region_sizes: Vec<usize>,
}

impl RegionIncidence {
    fn new(
        index: &RegionBoundaryIndex<'_>,
        partition: &SurfaceRegionPartition,
    ) -> Result<Self, RegionBoundaryError> {
        let evidence = index.evidence;
        if partition.region_by_polygon.len() != evidence.source_planar_support.len() {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let mut vertex_polygons = vec![Vec::new(); evidence.positions.len()];
        let mut region_vertices = vec![Vec::new(); partition.carriers.len()];
        let mut region_sizes = vec![0; partition.carriers.len()];
        for (polygon, &region) in partition.region_by_polygon.iter().enumerate() {
            if partition.carriers.get(region).is_none_or(Option::is_none) {
                return Err(RegionBoundaryError::InvalidPartition);
            }
            let vertices = evidence
                .source_polygon_vertices(polygon)
                .ok_or(RegionBoundaryError::InvalidEvidence)?;
            region_sizes[region] += 1;
            for &vertex in vertices {
                vertex_polygons
                    .get_mut(vertex)
                    .ok_or(RegionBoundaryError::InvalidEvidence)?
                    .push(polygon);
                region_vertices[region].push(vertex);
            }
        }
        for vertices in &mut region_vertices {
            vertices.sort_unstable();
            vertices.dedup();
        }
        Ok(Self {
            vertex_polygons,
            region_vertices,
            region_sizes,
        })
    }
}

struct FanSelection<'a, 'e> {
    index: &'a RegionBoundaryIndex<'e>,
    partition: &'a SurfaceRegionPartition,
    incidence: RegionIncidence,
    allowed: BTreeSet<usize>,
    selected: BTreeSet<usize>,
}

impl FanSelection<'_, '_> {
    fn region(&self, polygon: usize) -> Option<usize> {
        self.partition.region_by_polygon.get(polygon).copied()
    }

    // Follow the directed source fan, not region adjacency: one region
    // may occupy multiple disjoint sectors at the same vertex.
    fn next(&self, polygon: usize, vertex: usize) -> Result<Option<usize>, RegionBoundaryError> {
        self.index
            .opposite_polygon(polygon, self.index.outgoing_edge(polygon, vertex)?)
    }

    fn gap(&self, start: usize, vertex: usize) -> Option<BTreeSet<usize>> {
        let mut polygon = start;
        let mut added = BTreeSet::new();
        for _ in 0..self.incidence.vertex_polygons[vertex].len() {
            let region = self.region(polygon)?;
            if self.selected.contains(&region) {
                return (!added.is_empty()).then_some(added);
            }
            if !self.allowed.contains(&region) {
                return None;
            }
            added.insert(region);
            polygon = self.next(polygon, vertex).ok()??;
        }
        None
    }

    fn connectors(&self, vertex: usize) -> BTreeSet<usize> {
        let mut gaps = Vec::new();
        for &polygon in &self.incidence.vertex_polygons[vertex] {
            let Some(region) = self.region(polygon) else {
                return BTreeSet::new();
            };
            if !self.selected.contains(&region) {
                continue;
            }
            match self.next(polygon, vertex) {
                Ok(Some(next))
                    if self
                        .region(next)
                        .is_some_and(|id| self.selected.contains(&id)) => {}
                Ok(Some(next)) => gaps.push(self.gap(next, vertex)),
                Ok(None) => gaps.push(None),
                // No sector order can be established at this source edge.
                Err(_) => return BTreeSet::new(),
            }
        }
        if gaps.len() < 2 {
            return BTreeSet::new();
        }
        // A closed fan with only admissible gaps needs all but one gap
        // to connect its selected sectors. Leave the most expensive gap
        // outside; ties are deterministic proposals, never evidence.
        if gaps.iter().all(Option::is_some) {
            let largest = gaps
                .iter()
                .enumerate()
                .max_by_key(|(_, gap)| {
                    gap.iter()
                        .flatten()
                        .map(|&id| self.incidence.region_sizes[id])
                        .sum::<usize>()
                })
                .map(|(id, _)| id);
            if let Some(id) = largest {
                gaps.swap_remove(id);
            }
        }
        gaps.into_iter().flatten().flatten().collect()
    }
}

impl RegionBoundaryIndex<'_> {
    /// Propose the monotone closure of explicitly selected source regions
    /// through permitted connector regions at shared source vertices.
    ///
    /// Only bounded exterior sectors between selected sectors can add
    /// regions. Open, ambiguous and forbidden fans cannot authorize an
    /// addition. Each addition retains the full original region ownership;
    /// all touched vertices are reconsidered until the selection stabilizes.
    /// The returned sorted region IDs may still form multiple connected
    /// components or unresolved pinches. This is not a geometric certificate.
    ///
    /// # Errors
    /// Returns an error for inconsistent source references or unknown
    /// selected/connector regions.
    pub fn complete_patch_fans(
        &self,
        partition: &SurfaceRegionPartition,
        required: &[usize],
        connectors: &[usize],
    ) -> Result<Vec<usize>, RegionBoundaryError> {
        let incidence = RegionIncidence::new(self, partition)?;
        let selected = required.iter().copied().collect::<BTreeSet<_>>();
        let allowed = connectors.iter().copied().collect::<BTreeSet<_>>();
        if selected
            .iter()
            .chain(&allowed)
            .any(|&id| incidence.region_sizes.get(id).is_none_or(|&size| size == 0))
        {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let mut pending = selected
            .iter()
            .flat_map(|&id| incidence.region_vertices[id].iter().copied())
            .collect::<BTreeSet<_>>();
        let mut queue = pending.iter().copied().collect::<VecDeque<_>>();
        let mut selection = FanSelection {
            index: self,
            partition,
            incidence,
            allowed,
            selected,
        };
        while let Some(vertex) = queue.pop_front() {
            pending.remove(&vertex);
            for region in selection.connectors(vertex) {
                if selection.selected.insert(region) {
                    queue.extend(
                        selection.incidence.region_vertices[region]
                            .iter()
                            .copied()
                            .filter(|&id| pending.insert(id)),
                    );
                }
            }
        }
        Ok(selection.selected.into_iter().collect())
    }
}

impl SurfaceRegionPartition {
    /// Edge-connected components of an explicitly selected region set.
    /// Point-only adjacency is not a connection between source sectors.
    ///
    /// # Errors
    /// A selected ID does not refer to an active region.
    pub fn connected_region_groups(&self, regions:&[usize]) -> Result<Vec<Vec<usize>>,RegionBoundaryError> {
        let mut remaining = regions.iter().copied().collect::<BTreeSet<_>>();
        if remaining.iter().any(|&id|self.carriers.get(id).is_none_or(Option::is_none)) {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let mut adjacency=BTreeMap::<usize,Vec<usize>>::new();
        for path in &self.contacts.paths {
            if let (a,Some(b))=path.regions
                && remaining.contains(&a) && remaining.contains(&b) {
                adjacency.entry(a).or_default().push(b);
                adjacency.entry(b).or_default().push(a);
            }
        }
        let mut groups=Vec::new();
        while let Some(start)=remaining.pop_first() {
            let mut group=Vec::new();
            let mut pending=vec![start];
            while let Some(region)=pending.pop() {
                group.push(region);
                if let Some(neighbors)=adjacency.get(&region) {
                    pending.extend(neighbors.iter().copied().filter(|id|remaining.remove(id)));
                }
            }
            group.sort_unstable();
            groups.push(group);
        }
        Ok(groups)
    }
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
    fn fan_closure_uses_only_bounded_permitted_sectors() -> Result<(), Box<dyn std::error::Error>> {
        let mut positions = vec![Point3::new(0.0, 0.0, 0.0)];
        for i in 0..5 {
            let angle = f64::from(i) * std::f64::consts::TAU / 5.0;
            positions.push(Point3::new(
                angle.cos(),
                angle.sin(),
                0.4 + 0.2 * f64::from(i),
            ));
        }
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let partition = index.select_surface_partition(&[], &evidence, 100)?;
        let ids = &partition.region_by_polygon;
        assert_eq!(ids.iter().copied().collect::<BTreeSet<_>>().len(), 5);
        let required = [ids[0], ids[2]];
        let expected = [ids[0], ids[1], ids[2]]
            .into_iter()
            .collect::<BTreeSet<_>>();
        for allowed in [vec![ids[1]], vec![ids[1], ids[3], ids[4]]] {
            let selected = index.complete_patch_fans(&partition, &required, &allowed)?;
            assert_eq!(selected.into_iter().collect::<BTreeSet<_>>(), expected);
        }
        let selected = index.complete_patch_fans(&partition, &required, &[ids[3]])?;
        assert_eq!(
            selected.into_iter().collect::<BTreeSet<_>>(),
            required.into_iter().collect()
        );
        assert!(
            index
                .complete_patch_fans(&partition, &[usize::MAX], &[])
                .is_err()
        );
        Ok(())
    }
}
