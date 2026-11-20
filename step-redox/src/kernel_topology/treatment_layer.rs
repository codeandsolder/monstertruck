//! Dependency and conflict scheduling of independent inverse treatment groups.
use super::{
    KernelPatchEdgeRef, KernelTopologyId, KernelTopologyPatch, PlannedTreatmentTransaction,
    combine_treatment_transactions, merge_topology_patches, plan_treatment_component,
    treatment_components,
};
use crate::{
    edge_treatments::EdgeTreatmentEvidence, inverse_edge_treatment::TreatmentTopology,
    reduced_brep::ReducedSolid,
};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};

pub struct PlannedTreatmentLayer {
    pub transaction: PlannedTreatmentTransaction,
    pub evidence: Vec<EdgeTreatmentEvidence>,
}

/// Candidate footprint in the ORIGINAL topology of this layer. Different
/// splices on one face are compatible, provided their removed runs and
/// surviving references don't conflict.
#[derive(Default)]
struct Footprint {
    removed_faces: BTreeSet<KernelTopologyId>,
    edited_faces: BTreeSet<KernelTopologyId>,
    replaced_faces: BTreeSet<KernelTopologyId>,
    removed_edges: BTreeSet<KernelTopologyId>,
    retained_edges: BTreeSet<KernelTopologyId>,
    rewritten_edges: BTreeSet<KernelTopologyId>,
    splice_edges: BTreeSet<(KernelTopologyId, KernelTopologyId)>,
}

impl Footprint {
    fn new(patch: &KernelTopologyPatch) -> Self {
        let mut result = Self {
            removed_faces: patch.remove_faces.iter().copied().collect(),
            removed_edges: patch.require_unused_edges.iter().copied().collect(),
            rewritten_edges: patch.rewrite_edges.iter().map(|edge| edge.edge).collect(),
            ..Self::default()
        };
        let mut add_edge = |reference: KernelPatchEdgeRef| {
            if let KernelPatchEdgeRef::Existing(id) = reference {
                result.retained_edges.insert(id);
            }
        };
        for splice in &patch.face_splices {
            result.edited_faces.insert(splice.face);
            result
                .splice_edges
                .extend(splice.remove_edges.iter().map(|&edge| (splice.face, edge)));
            for edge in &splice.replacement {
                add_edge(edge.edge);
            }
        }
        for replacement in &patch.face_boundary_replacements {
            result.replaced_faces.insert(replacement.face);
            result.edited_faces.insert(replacement.face);
            for edge in replacement.boundaries.iter().flatten() {
                add_edge(edge.edge);
            }
        }
        result.retained_edges.extend(&result.rewritten_edges);
        result
    }

    fn conflicts(&self, other: &Self) -> bool {
        !self.removed_faces.is_disjoint(&other.edited_faces)
            || !other.removed_faces.is_disjoint(&self.edited_faces)
            || !self.removed_faces.is_disjoint(&other.removed_faces)
            || !self.replaced_faces.is_disjoint(&other.edited_faces)
            || !other.replaced_faces.is_disjoint(&self.edited_faces)
            || !self.removed_edges.is_disjoint(&other.retained_edges)
            || !other.removed_edges.is_disjoint(&self.retained_edges)
            || !self.splice_edges.is_disjoint(&other.splice_edges)
    }
}

struct Candidate {
    plan: PlannedTreatmentTransaction,
    evidence: Vec<EdgeTreatmentEvidence>,
    footprint: Footprint,
}

impl Candidate {
    fn new(plan: PlannedTreatmentTransaction, evidence: Vec<EdgeTreatmentEvidence>) -> Self {
        let footprint = Footprint::new(&plan.patch);
        Self {
            plan,
            evidence,
            footprint,
        }
    }
}

/// If two independently plausible plans contradict one another, retain both
/// on the solid rather than assigning a winner by STEP entity order. Only
/// the compatible groups will be redetected after the next native edit.
fn compatible_candidates(candidates: Vec<Candidate>) -> (Vec<Candidate>, Vec<String>) {
    let mut blocked = BTreeSet::new();
    let mut conflicts = Vec::new();
    for (i, left) in candidates.iter().enumerate() {
        for (j, right) in candidates.iter().enumerate().skip(i + 1) {
            let reason = if left.footprint.conflicts(&right.footprint) {
                Some("overlapping removal or boundary edits".to_owned())
            } else if !left
                .footprint
                .rewritten_edges
                .is_disjoint(&right.footprint.rewritten_edges)
            {
                merge_topology_patches(vec![left.plan.patch.clone(), right.plan.patch.clone()])
                    .err()
                    .map(|error| format!("{error:#}"))
            } else {
                None
            };
            if let Some(reason) = reason {
                blocked.extend([i, j]);
                let ids = |candidate: &Candidate| {
                    candidate
                        .evidence
                        .iter()
                        .map(|item| item.treatment_face_id)
                        .collect::<Vec<_>>()
                };
                conflicts.push(format!("{:?} and {:?}: {reason}", ids(left), ids(right)));
            }
        }
    }
    (
        candidates
            .into_iter()
            .enumerate()
            .filter_map(|(i, candidate)| (!blocked.contains(&i)).then_some(candidate))
            .collect(),
        conflicts,
    )
}

fn carrier_users(evidence: &[EdgeTreatmentEvidence]) -> BTreeMap<u64, BTreeSet<u64>> {
    let mut users = BTreeMap::<u64, BTreeSet<u64>>::new();
    for item in evidence {
        for carrier in &item.carriers {
            users
                .entry(carrier.face_id)
                .or_default()
                .insert(item.treatment_face_id);
        }
    }
    users
}

fn ready(group: &[&EdgeTreatmentEvidence], users: &BTreeMap<u64, BTreeSet<u64>>) -> bool {
    let ids = group
        .iter()
        .map(|item| item.treatment_face_id)
        .collect::<BTreeSet<_>>();
    ids.iter().all(|id| {
        users
            .get(id)
            .is_none_or(|dependents| dependents.is_subset(&ids))
    })
}

/// Partition the compatible plans by original vertex/edge incidence. An
/// unrelated operation on another boundary run of the same carrier face
/// need not share its rollback fate. Adjacent endpoint rewrites must remain
/// together because individually they may leave an open boundary.
fn neighborhoods(
    solid: &ReducedSolid,
    section: usize,
    provenance: &impl TreatmentTopology,
    candidates: Vec<Candidate>,
) -> Result<Vec<Vec<Candidate>>> {
    let mut touched = Vec::with_capacity(candidates.len());
    for candidate in &candidates {
        let fp = &candidate.footprint;
        let mut edges = fp.removed_edges.clone();
        edges.extend(&fp.retained_edges);
        edges.extend(fp.splice_edges.iter().map(|(_, edge)| *edge));
        touched.push(edges);
    }
    let all_edges = touched.iter().flatten().copied().collect::<BTreeSet<_>>();
    let mut endpoints = BTreeMap::new();
    for edge in solid
        .faces()
        .flat_map(crate::reduced_brep::ReducedFace::loop_edges)
    {
        let Some(id) = provenance.edge(section, edge.edge_id) else {
            continue;
        };
        if all_edges.contains(&id) && !endpoints.contains_key(&id) {
            let a = provenance
                .vertex(section, edge.start_vertex)
                .context("missing neighborhood endpoint")?;
            let b = provenance
                .vertex(section, edge.end_vertex)
                .context("missing neighborhood endpoint")?;
            endpoints.insert(id, [a, b]);
        }
    }
    let mut owners = BTreeMap::<KernelTopologyId, Vec<usize>>::new();
    for (i, edges) in touched.iter().enumerate() {
        let mut resources = edges.clone();
        for edge in edges {
            resources.extend(endpoints.get(edge).context("missing neighborhood edge")?);
        }
        for resource in resources {
            owners.entry(resource).or_default().push(i);
        }
    }
    let mut adjacency = vec![Vec::new(); candidates.len()];
    for indices in owners.values() {
        if let Some((&first, rest)) = indices.split_first() {
            for &other in rest {
                adjacency[first].push(other);
                adjacency[other].push(first);
            }
        }
    }
    let mut seen = BTreeSet::new();
    let mut candidates = candidates.into_iter().map(Some).collect::<Vec<_>>();
    let mut result = Vec::new();
    for start in 0..candidates.len() {
        if !seen.insert(start) {
            continue;
        }
        let mut stack = vec![start];
        let mut component = Vec::new();
        while let Some(i) = stack.pop() {
            component.push(
                candidates[i]
                    .take()
                    .context("duplicated neighborhood visit")?,
            );
            for &other in &adjacency[i] {
                if seen.insert(other) {
                    stack.push(other);
                }
            }
        }
        result.push(component);
    }
    Ok(result)
}

/// Plan one outer inverse layer. An atomic component is ready only when
/// no other candidate uses one of its faces as a carrier. Independently
/// proven, compatible plans can be committed together. Conflicting plans
/// remain on the solid, and are only reconsidered against fresh geometry.
pub fn plan_treatment_layer(
    solid: &ReducedSolid,
    section: usize,
    evidence: &[EdgeTreatmentEvidence],
    provenance: &impl TreatmentTopology,
) -> Result<Vec<PlannedTreatmentLayer>> {
    let users = carrier_users(evidence);
    let mut candidates = Vec::new();
    let mut refusals = Vec::new();
    for (closed, group) in treatment_components(solid, evidence) {
        if !ready(&group, &users) {
            continue;
        }
        match plan_treatment_component(solid, section, closed, &group, provenance) {
            Ok(plan) => candidates.push(Candidate::new(plan, group.into_iter().cloned().collect())),
            Err(error) => refusals.push(format!(
                "{:?}: {error:#}",
                group
                    .iter()
                    .map(|item| item.treatment_face_id)
                    .collect::<Vec<_>>()
            )),
        }
    }
    let (candidates, conflicts) = compatible_candidates(candidates);
    refusals.extend(conflicts);
    if candidates.is_empty() {
        if refusals.is_empty() {
            return Ok(Vec::new());
        }
        bail!(
            "no independently plannable outer edge treatment: {}",
            refusals.join("; ")
        );
    }
    neighborhoods(solid, section, provenance, candidates)?
        .into_iter()
        .map(|component| {
            let mut plans = Vec::with_capacity(component.len());
            let mut evidence = Vec::new();
            for candidate in component {
                plans.push(candidate.plan);
                evidence.extend(candidate.evidence);
            }
            Ok(PlannedTreatmentLayer {
                transaction: combine_treatment_transactions(plans)?,
                evidence,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel_topology::{KernelPatchCurve, KernelPatchVertexRef, KernelRetainedEdge};

    fn endpoint_edit(face: u64, front: bool, chosen: [f64; 3]) -> Candidate {
        let other = if front { [2.0, 0.0, 0.0] } else { [0.0; 3] };
        let (vertices, points) = if front {
            (
                (
                    KernelPatchVertexRef::New(0),
                    KernelPatchVertexRef::Existing(KernelTopologyId::from_raw(12)),
                ),
                [chosen, other],
            )
        } else {
            (
                (
                    KernelPatchVertexRef::Existing(KernelTopologyId::from_raw(11)),
                    KernelPatchVertexRef::New(0),
                ),
                [other, chosen],
            )
        };
        Candidate::new(
            PlannedTreatmentTransaction {
                coupled_selectors: Vec::new(),
                patch: KernelTopologyPatch {
                    remove_faces: vec![KernelTopologyId::from_raw(face)],
                    new_vertices: vec![chosen],
                    new_edges: Vec::new(),
                    rewrite_edges: vec![KernelRetainedEdge {
                        edge: KernelTopologyId::from_raw(100),
                        vertices,
                        curve: KernelPatchCurve::Line {
                            start_mm: points[0],
                            end_mm: points[1],
                        },
                    }],
                    face_splices: Vec::new(),
                    face_boundary_replacements: Vec::new(),
                    require_unused_edges: Vec::new(),
                },
            },
            Vec::new(),
        )
    }

    #[test]
    fn native_treatment_layers_conflicting_junctions_do_not_block_independent_edits() {
        let first = endpoint_edit(1, true, [1.0, 0.0, 0.0]);
        let second = endpoint_edit(2, false, [1.0, 0.0, 0.0]);
        let mut independent = endpoint_edit(3, true, [1.0, 0.0, 0.0]);
        independent.plan.patch.rewrite_edges[0].edge = KernelTopologyId::from_raw(200);
        independent.footprint = Footprint::new(&independent.plan.patch);
        let (accepted, conflicts) = compatible_candidates(vec![first, second, independent]);
        assert_eq!(conflicts.len(), 1);
        assert!(conflicts[0].contains("collapses"));
        assert_eq!(accepted.len(), 1);
        assert_eq!(
            accepted[0].plan.patch.remove_faces,
            vec![KernelTopologyId::from_raw(3)]
        );
    }

    #[test]
    fn native_treatment_layers_compatible_opposite_endpoints_merge() -> Result<()> {
        let first = endpoint_edit(1, true, [0.5, 0.0, 0.0]);
        let second = endpoint_edit(2, false, [1.5, 0.0, 0.0]);
        let (accepted, conflicts) = compatible_candidates(vec![first, second]);
        assert!(conflicts.is_empty());
        let combined =
            combine_treatment_transactions(accepted.into_iter().map(|c| c.plan).collect())?;
        assert_eq!(combined.patch.rewrite_edges.len(), 1);
        assert_eq!(
            combined.patch.rewrite_edges[0].curve,
            KernelPatchCurve::Line {
                start_mm: [0.5, 0.0, 0.0],
                end_mm: [1.5, 0.0, 0.0],
            }
        );
        Ok(())
    }
}
