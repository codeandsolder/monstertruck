//! Budget-aware planar normalization of complete hybrid recovery results.

use super::{HybridPatchEvidence, HybridPolygonFaces, HybridRecovery};
use crate::{
    Edge, EdgeId, Face, KERNEL_EPSILON, Shell,
    mesh_recovery::{PlanarFanOptions, simplify_planar_fan_shell},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Range,
};

/// Summary of a committed, budget-preserving hybrid simplification.
#[derive(Debug, Default, Clone, Copy)]
pub struct HybridPlanarSimplification {
    /// Faces before the transaction.
    pub faces_before: usize,
    /// Faces after the transaction.
    pub faces_after: usize,
    /// Number of certified planar face groups.
    pub merged_groups: usize,
    /// Number of candidate groups for which no proof was established.
    pub rejected_groups: usize,
    /// Number of whole patch proposals that failed validation.
    pub rejected_patches: usize,
    /// Maximum additional deviation of any simplified patch.
    pub added_deviation_upper_bound: f64,
}

struct PatchReplacement {
    regions: BTreeMap<usize, Vec<Face>>,
    polygons: Vec<HybridPolygonFaces>,
    added_deviation: f64,
    merged_groups: usize,
    rejected_groups: usize,
}

struct PatchInput<'a> {
    recovery: &'a HybridRecovery,
    region_indices: &'a BTreeMap<usize, usize>,
    protected: &'a [EdgeId],
}

impl PatchInput<'_> {
    fn propose(&self, patch: &HybridPatchEvidence) -> Option<PatchReplacement> {
        let limit = (patch.budget.vertex - patch.vertex_displacement)
            .min(patch.budget.chord - patch.deviation_upper_bound)
            .min(KERNEL_EPSILON.min(patch.budget.chord) - patch.seam_deviation_upper_bound);
        if !limit.is_finite() || limit < 0.0 {
            return None;
        }
        let mut shell = Shell::new();
        let mut owners = Vec::new();
        let mut ranges = BTreeMap::<usize, Range<usize>>::new();
        for id in &patch.regions {
            let index = *self.region_indices.get(id)?;
            let faces = self.recovery.native.regions[index].faces.as_ref().ok()?;
            let first = shell.len();
            shell.extend(faces.iter().cloned());
            owners.extend(std::iter::repeat_n(index, faces.len()));
            if ranges.insert(*id, first..shell.len()).is_some() {
                return None;
            }
        }
        let result = simplify_planar_fan_shell(
            &shell,
            PlanarFanOptions {
                owner_by_face: Some(&owners),
                protected_edges: self.protected,
                ..PlanarFanOptions::new(limit)
            },
        )?;
        let added_deviation = result.deviation_upper_bound;
        if patch.vertex_displacement + added_deviation > patch.budget.vertex
            || patch.deviation_upper_bound + added_deviation > patch.budget.chord
            || patch.seam_deviation_upper_bound + added_deviation
                > KERNEL_EPSILON.min(patch.budget.chord)
        {
            return None;
        }
        // Assign every simplified face to its original region. A union
        // cannot cross owners. Keep region-local indices for provenance.
        let mut new_owners = vec![None; result.shell.len()];
        for (&owner, &new) in owners.iter().zip(&result.original_to_simplified) {
            let slot = new_owners.get_mut(new)?;
            if slot.is_some_and(|existing| existing != owner) {
                return None;
            }
            *slot = Some(owner);
        }
        let mut regions = BTreeMap::<usize, Vec<Face>>::new();
        let mut local_indices = Vec::with_capacity(new_owners.len());
        for (face, owner) in result.shell.iter().zip(new_owners) {
            let native = regions.entry(owner?).or_default();
            local_indices.push(native.len());
            native.push(face.clone());
        }
        let polygons = patch
            .polygons
            .iter()
            .map(|polygon| {
                let range = ranges.get(&polygon.region)?;
                let faces = polygon
                    .faces
                    .iter()
                    .map(|&index| {
                        let old = range.start.checked_add(index)?;
                        if old >= range.end {
                            return None;
                        }
                        local_indices
                            .get(*result.original_to_simplified.get(old)?)
                            .copied()
                    })
                    .collect::<Option<BTreeSet<_>>>()?;
                if faces.is_empty() {
                    return None;
                }
                Some(HybridPolygonFaces {
                    polygon: polygon.polygon,
                    region: polygon.region,
                    faces: faces.into_iter().collect(),
                })
            })
            .collect::<Option<Vec<_>>>()?;
        Some(PatchReplacement {
            regions,
            polygons,
            added_deviation,
            merged_groups: result.merged_groups,
            rejected_groups: result.rejected_groups,
        })
    }
}

impl HybridRecovery {
    /// Simplify only certified approximate fans of a complete hybrid shell.
    /// Preserve every original native trim, region owner and source polygon
    /// provenance. Charge the additional geometric deviation to the original
    /// vertex, chord and seam bounds; no source tolerance is relaxed.
    /// Unproved patch proposals are left unchanged. None means the global
    /// transaction was invalid, in which case the recovery is unchanged.
    pub fn simplify_planar_fans(&mut self) -> Option<HybridPlanarSimplification> {
        let original = self.complete_shell()?;
        let region_indices = self
            .native
            .regions
            .iter()
            .enumerate()
            .map(|(index, region)| (region.region, index))
            .collect::<BTreeMap<_, _>>();
        if region_indices.len() != self.native.regions.len() {
            return None;
        }
        let protected = self
            .native
            .trims
            .iter()
            .filter_map(|trim| trim.as_ref().ok())
            .flatten()
            .map(Edge::id)
            .collect::<Vec<_>>();
        let input = PatchInput {
            recovery: self,
            region_indices: &region_indices,
            protected: &protected,
        };
        let mut report = HybridPlanarSimplification {
            faces_before: original.len(),
            faces_after: original.len(),
            ..Default::default()
        };
        let mut replacements = BTreeMap::<usize, Vec<Face>>::new();
        let mut patches = BTreeMap::<usize, PatchReplacement>::new();
        for (index, patch) in self.patches.iter().enumerate() {
            let Some(mut proposal) = input.propose(patch) else {
                report.rejected_patches += 1;
                continue;
            };
            report.rejected_groups += proposal.rejected_groups;
            if proposal.merged_groups == 0 {
                continue;
            }
            report.merged_groups += proposal.merged_groups;
            report.added_deviation_upper_bound = report
                .added_deviation_upper_bound
                .max(proposal.added_deviation);
            for (region, faces) in std::mem::take(&mut proposal.regions) {
                if replacements.insert(region, faces).is_some() {
                    return None;
                }
            }
            patches.insert(index, proposal);
        }
        if replacements.is_empty() {
            return Some(report);
        }
        let mut shell = Shell::new();
        for (index, region) in self.native.regions.iter().enumerate() {
            let faces = match replacements.get(&index) {
                Some(faces) => faces,
                None => region.faces.as_ref().ok()?,
            };
            shell.extend(faces.iter().cloned());
        }
        if shell.shell_condition() != original.shell_condition()
            || shell.connected_components().len() != original.connected_components().len()
            || !shell.singular_vertices().is_empty()
            || !shell.is_geometric_consistent()
        {
            return None;
        }
        report.faces_after = shell.len();
        for (index, faces) in replacements {
            self.native.regions[index].faces = Ok(faces);
        }
        for (index, proposal) in patches {
            let patch = &mut self.patches[index];
            patch.polygons = proposal.polygons;
            patch.vertex_displacement += proposal.added_deviation;
            patch.deviation_upper_bound += proposal.added_deviation;
            patch.seam_deviation_upper_bound += proposal.added_deviation;
        }
        Some(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Point3, builder,
        mesh_recovery::{
            BoundaryPatchOptions, NativeRegionFace, NativeTrimGraph, PrescribedSourceEdge,
            SurfaceRegionBudget, approximate_boundary_patch,
        },
    };
    use monstertruck_topology::shell::ShellCondition;

    #[test]
    fn hybrid_planar_simplification_keeps_region_and_trim_provenance()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(0., 0., 0.),
            Point3::new(2., 0., 0.),
            Point3::new(2., 2., 0.),
            Point3::new(0., 2., 0.),
        ];
        let vertices = points.map(builder::vertex);
        let diagonal = builder::line(&vertices[0], &vertices[2]);
        let source: [&[usize]; 2] = [&[0, 1, 2], &[0, 2, 3]];
        let budget = SurfaceRegionBudget {
            vertex: 1e-5,
            chord: 0.1,
        };
        let patch = approximate_boundary_patch(
            &points,
            &source,
            &[PrescribedSourceEdge {
                vertices: [0, 2],
                edge: diagonal.clone(),
            }],
            BoundaryPatchOptions::new(budget),
        )?;
        let regions = patch
            .source_face_ranges
            .iter()
            .enumerate()
            .map(|(region, range)| NativeRegionFace {
                region,
                faces: Ok(patch.shell[range.clone()].to_vec()),
                requires_chart_cut: false,
                chart_obstacles: Vec::new(),
                chart_seams: 0,
            })
            .collect();
        let polygons = (0..2)
            .map(|polygon| HybridPolygonFaces {
                polygon,
                region: polygon,
                faces: vec![0, 1, 2],
            })
            .collect();
        let mut result = HybridRecovery {
            native: NativeTrimGraph {
                vertices: vertices.into_iter().map(Some).collect(),
                trims: vec![Ok(vec![diagonal.clone()])],
                regions,
                max_vertex_displacement: 0.,
            },
            patches: vec![HybridPatchEvidence {
                budget,
                regions: vec![0, 1],
                polygons,
                coalesced_cells: 0,
                vertex_displacement: patch.vertex_displacement,
                deviation_upper_bound: patch.deviation_upper_bound,
                seam_deviation_upper_bound: patch.seam_deviation_upper_bound,
            }],
            lifted_trims: BTreeMap::new(),
            ambiguous_edges: false,
        };
        let report = result.simplify_planar_fans().ok_or("transaction")?;
        assert_eq!((report.faces_before, report.faces_after), (6, 2));
        assert_eq!(report.merged_groups, 2);
        let shell = result.complete_shell().ok_or("complete shell")?;
        assert_eq!(shell.shell_condition(), ShellCondition::Oriented);
        assert_eq!(
            shell
                .edge_iter()
                .filter(|e| e.id() == diagonal.id())
                .count(),
            2
        );
        assert!(
            result.patches[0]
                .polygons
                .iter()
                .all(|polygon| polygon.faces == [0])
        );
        assert!(result.patches[0].vertex_displacement <= budget.vertex);
        assert!(result.patches[0].deviation_upper_bound <= budget.chord);
        Ok(())
    }
}
