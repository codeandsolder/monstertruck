//! Refine planar evidence only where source-backed curved hypotheses need
//! original facet ownership. Preserve every unaffected connected planar
//! domain, its original seed plane, and exact source incidence.

use super::{
    MeshEvidence, PlanarSupport, build_support_interfaces, measure_planar_support, pinched_planar,
    regular_polygon_neighbors,
};
use crate::KERNEL_EPSILON;
use crate::mesh_recovery::{
    InputPolygon, MeshRecoveryError, input_polygon_from_ids, source_edge_incidence,
};
use rustc_hash::FxHashSet;

// Selected facets become singletons. The remainder is partitioned into
// regular-edge-connected children, never joined across a vertex pinch.
fn connected_remainder(
    support: &PlanarSupport,
    selected: &[bool],
    neighbors: &[Vec<usize>],
) -> Vec<Vec<usize>> {
    let mut remaining = support
        .source_face_indices
        .iter()
        .copied()
        .filter(|&id| !selected[id])
        .collect::<FxHashSet<_>>();
    let mut groups = Vec::new();
    for &seed in &support.source_face_indices {
        if selected[seed] {
            groups.push(vec![seed]);
            continue;
        }
        if !remaining.remove(&seed) {
            continue;
        }
        let mut group = vec![seed];
        let mut cursor = 0;
        while let Some(&id) = group.get(cursor) {
            for &other in &neighbors[id] {
                if remaining.remove(&other) {
                    group.push(other);
                }
            }
            cursor += 1;
        }
        group.sort_unstable();
        groups.push(group);
    }
    groups
}

fn assign_membership(
    supports: &[PlanarSupport],
    polygon_count: usize,
) -> Result<Vec<usize>, MeshRecoveryError> {
    let mut membership = vec![usize::MAX; polygon_count];
    for (id, support) in supports.iter().enumerate() {
        for &polygon in &support.source_face_indices {
            let slot = membership
                .get_mut(polygon)
                .ok_or(MeshRecoveryError::DegeneratePolygon(polygon))?;
            if *slot != usize::MAX {
                return Err(MeshRecoveryError::DegeneratePolygon(polygon));
            }
            *slot = id;
        }
    }
    if membership.contains(&usize::MAX) {
        return Err(MeshRecoveryError::InvalidEvidenceOptions);
    }
    Ok(membership)
}

impl MeshEvidence {
    /// Isolate the specified source facets as individual supports. Facets are
    /// indexed in the active source polygon array, not the pre-normalization
    /// input. All other polygons retain their existing seed plane and remain
    /// grouped wherever regular edge connectivity and simple face domains
    /// permit. Preserve source positions, winding, feature observations and
    /// complete original polygon provenance.
    ///
    /// # Errors
    /// Invalid source polygon references or inconsistent source incidence.
    pub fn refine_source_facets(&self, facets: &[usize]) -> Result<Self, MeshRecoveryError> {
        let count = self.source_planar_support.len();
        let mut selected = vec![false; count];
        for &id in facets {
            *selected
                .get_mut(id)
                .ok_or(MeshRecoveryError::InvalidEvidenceOptions)? = true;
        }
        let polygons = (0..count)
            .map(|id| {
                let vertices = self
                    .source_polygon_vertices(id)
                    .ok_or(MeshRecoveryError::DegeneratePolygon(id))?;
                input_polygon_from_ids(&self.positions, vertices.to_vec(), id, KERNEL_EPSILON)
            })
            .collect::<Result<Vec<InputPolygon>, _>>()?;
        let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
        let neighbors = regular_polygon_neighbors(&incidence, count);
        let mut supports = Vec::with_capacity(self.planar_supports.len());
        for support in &self.planar_supports {
            if support.source_face_indices.len() == 1
                || !support.source_face_indices.iter().any(|&id| selected[id])
            {
                supports.push(support.clone());
                continue;
            }
            let mut children = Vec::new();
            for group in connected_remainder(support, &selected, &neighbors) {
                if group.len() == 1 {
                    children.push(group);
                } else {
                    children.extend(pinched_planar::split(
                        &polygons, &incidence, &neighbors, group,
                    ));
                }
            }
            children.sort_unstable_by_key(|group| group[0]);
            for group in children {
                let seed = &polygons[group[0]];
                // Only an explicitly isolated facet gets its own fitted
                // plane. Retain the already certified parent carrier for
                // all remaining children, including a singleton remainder.
                let (origin, normal) = if group.len() == 1 && selected[group[0]] {
                    (seed.origin, seed.normal)
                } else {
                    (support.origin, support.normal)
                };
                supports.push(measure_planar_support(
                    &polygons,
                    &self.positions,
                    group,
                    origin,
                    normal,
                ));
            }
        }
        let membership = assign_membership(&supports, count)?;
        let (support_interfaces, support_interface_indices) = build_support_interfaces(
            &self.positions,
            &polygons,
            &incidence,
            &membership,
            supports.len(),
        );
        let mut refined = self.clone();
        refined.planar_supports = supports;
        refined.source_planar_support = membership;
        refined.support_interfaces = support_interfaces;
        refined.support_interface_indices = support_interface_indices;
        Ok(refined)
    }
}
