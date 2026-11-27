//! Native-consistency feedback for competing source partitions.
//! Compare complete region neighborhoods, not individual faces: changing
//! ownership at one junction changes the shared trim constraints on all
//! incident regions. Neither alternative changes original source ownership.

use super::{
    AnalyticRecoveryError, AnalyticRecoveryOptions, MeshEvidence, NativeTrimGraph,
    RegionBoundaryIndex, SelectedSurfaceProposal, SurfaceCarrier, SurfaceRegionPartition,
    assemble_native,
};
use std::collections::BTreeSet;
mod feedback;
mod ownership;
use ownership::OwnershipTrial;

#[derive(Clone, Debug)]
pub(super) struct PartitionAlternative {
    pub partition: SurfaceRegionPartition,
    pub selected: Vec<SelectedSurfaceProposal>,
}

struct Groups {
    parents: Vec<usize>,
}

impl Groups {
    fn new(count: usize) -> Self {
        Self {
            parents: (0..count).collect(),
        }
    }

    fn root(&mut self, mut id: usize) -> usize {
        while self.parents[id] != id {
            self.parents[id] = self.parents[self.parents[id]];
            id = self.parents[id];
        }
        id
    }

    fn join(&mut self, a: usize, b: usize) {
        let a = self.root(a);
        let b = self.root(b);
        self.parents[a.max(b)] = a.min(b);
    }
}

// First close each group under both region partitions. Then include every
// region meeting a changed group at an original vertex. Changes sharing an
// incident region are inseparable, even if their own polygons do not touch.
fn affected_groups(
    evidence: &MeshEvidence,
    first: &PartitionAlternative,
    second: &PartitionAlternative,
) -> Vec<usize> {
    let count = evidence.planar_supports.len();
    let mut groups = Groups::new(count);
    for partition in [&first.partition, &second.partition] {
        for (support, &region) in partition.region_by_support.iter().enumerate() {
            groups.join(support, region);
        }
    }
    let common = first
        .partition
        .selected_proposals
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .intersection(
            &second
                .partition
                .selected_proposals
                .iter()
                .copied()
                .collect(),
        )
        .copied()
        .collect::<BTreeSet<_>>();
    let mut changed = vec![false; count];
    for alternative in [first, second] {
        for (&id, proposal) in alternative
            .partition
            .selected_proposals
            .iter()
            .zip(&alternative.selected)
        {
            if !common.contains(&id) {
                changed[groups.root(proposal.supports[0])] = true;
            }
        }
    }
    let atomic = (0..count).map(|id| groups.root(id)).collect::<Vec<_>>();
    let mut incident = vec![Vec::new(); evidence.positions.len()];
    for (polygon, &support) in evidence.source_planar_support.iter().enumerate() {
        if let Some(vertices) = evidence.source_polygon_vertices(polygon) {
            for &vertex in vertices {
                incident[vertex].push(atomic[support]);
            }
        }
    }
    for at_vertex in incident {
        if let Some(&seed) = at_vertex.iter().find(|&&id| changed[id]) {
            for id in at_vertex {
                groups.join(seed, id);
            }
        }
    }
    (0..count).map(|id| groups.root(id)).collect()
}

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
struct Coverage {
    native_polygons: usize,
    analytic_polygons: usize,
}

fn group_coverage(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    native: &NativeTrimGraph,
    groups: &[usize],
) -> Vec<Coverage> {
    let count = evidence.planar_supports.len();
    let mut valid = vec![false; count];
    for region in &native.regions {
        valid[region.region] = region.faces.is_ok() && !region.requires_chart_cut;
    }
    let mut result = vec![Coverage::default(); count];
    for (support, &region) in partition.region_by_support.iter().enumerate() {
        if valid[region] {
            let weight = evidence.planar_supports[support].source_face_indices.len();
            let coverage = &mut result[groups[support]];
            coverage.native_polygons += weight;
            if !matches!(
                partition.carriers[region],
                Some(SurfaceCarrier::Plane { .. })
            ) {
                coverage.analytic_polygons += weight;
            }
        }
    }
    result
}

fn total_coverage(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    native: &NativeTrimGraph,
) -> Coverage {
    let groups = vec![0; evidence.planar_supports.len()];
    group_coverage(evidence, partition, native, &groups)
        .first()
        .copied()
        .unwrap_or_default()
}

// Each chosen neighborhood contains complete regions from either partition.
// Region IDs remain original support IDs, so they cannot collide across
// independently chosen neighborhoods.
fn combine(
    evidence: &MeshEvidence,
    first: &PartitionAlternative,
    second: &PartitionAlternative,
    take_second: &[bool],
) -> Result<PartitionAlternative, AnalyticRecoveryError> {
    let count = evidence.planar_supports.len();
    let mut region_by_support = Vec::with_capacity(count);
    let mut carriers = vec![None; count];
    for (support, &second_choice) in take_second.iter().enumerate() {
        let partition = if second_choice {
            &second.partition
        } else {
            &first.partition
        };
        let region = partition.region_by_support[support];
        region_by_support.push(region);
        carriers[region] = partition.carriers[region];
    }
    let region_by_polygon = evidence
        .source_planar_support
        .iter()
        .map(|&id| region_by_support[id])
        .collect::<Vec<_>>();
    let contacts = RegionBoundaryIndex::new(evidence)?.partition_contacts(&region_by_polygon)?;
    let mut chosen = Vec::new();
    for (index, alternative) in [first, second].into_iter().enumerate() {
        for (&id, proposal) in alternative
            .partition
            .selected_proposals
            .iter()
            .zip(&alternative.selected)
        {
            if take_second[proposal.supports[0]] == (index == 1) {
                chosen.push((id, proposal.clone()));
            }
        }
    }
    chosen.sort_unstable_by_key(|(id, _)| *id);
    let (selected_proposals, selected) = chosen.into_iter().unzip();
    let consolidated_supports = count - carriers.iter().flatten().count();
    Ok(PartitionAlternative {
        partition: SurfaceRegionPartition {
            selected_proposals,
            region_by_support,
            region_by_polygon,
            carriers,
            contacts,
            consolidated_supports,
            exact_components: second.partition.exact_components,
            limited_components: second.partition.limited_components,
        },
        selected,
    })
}

pub(super) struct NativeSelection {
    pub alternative: PartitionAlternative,
    pub assembly: super::NativeAssembly,
    pub compared_neighborhoods: usize,
    pub expanded_neighborhoods: usize,
    pub ownership: Option<OwnershipTrial>,
    pub ownership_rejection: Option<super::JointFitAttempt>,
}

enum RawChoice {
    Core,
    Expanded,
    Combined(Box<(PartitionAlternative, super::NativeAssembly)>),
}

struct RawDecision {
    choice: RawChoice,
    compared_neighborhoods: usize,
    expanded_neighborhoods: usize,
}

/// Compare independent ownership alternatives. Reconcile displaced source
/// junctions on a competing expanded partition before accepting it, without
/// allowing new native failures or invalidating previously successful supports.
pub(super) fn select_native(
    evidence: &MeshEvidence,
    core: PartitionAlternative,
    expanded: PartitionAlternative,
    options: AnalyticRecoveryOptions,
) -> Result<NativeSelection, AnalyticRecoveryError> {
    let core_assembly = assemble_native(evidence, core.partition.clone(), options)?;
    let expanded_assembly = assemble_native(evidence, expanded.partition.clone(), options)?;
    let decision = select_raw(
        evidence,
        &core,
        &expanded,
        &core_assembly,
        &expanded_assembly,
        options,
    )?;
    let baseline = match &decision.choice {
        RawChoice::Core => &core_assembly,
        RawChoice::Expanded => &expanded_assembly,
        RawChoice::Combined(combined) => &combined.1,
    };
    let outcome = ownership::attempt(
        evidence,
        &core,
        &expanded,
        &expanded_assembly,
        baseline,
        options,
    )?;
    let mut ownership_rejection = None;
    if let Some(trial) = outcome {
        match trial {
            Ok((alternative, assembly, fit)) => {
                return Ok(NativeSelection {
                    alternative,
                    assembly,
                    compared_neighborhoods: decision.compared_neighborhoods,
                    expanded_neighborhoods: decision.compared_neighborhoods,
                    ownership: Some(fit),
                    ownership_rejection: None,
                });
            }
            Err(rejection) => ownership_rejection = Some(rejection),
        }
    }
    let (alternative, assembly) = match decision.choice {
        RawChoice::Core => (core, core_assembly),
        RawChoice::Expanded => (expanded, expanded_assembly),
        RawChoice::Combined(combined) => *combined,
    };
    Ok(NativeSelection {
        alternative,
        assembly,
        compared_neighborhoods: decision.compared_neighborhoods,
        expanded_neighborhoods: decision.expanded_neighborhoods,
        ownership: None,
        ownership_rejection,
    })
}

fn select_raw(
    evidence: &MeshEvidence,
    core: &PartitionAlternative,
    expanded: &PartitionAlternative,
    core_assembly: &super::NativeAssembly,
    expanded_assembly: &super::NativeAssembly,
    options: AnalyticRecoveryOptions,
) -> Result<RawDecision, AnalyticRecoveryError> {
    let groups = affected_groups(evidence, core, expanded);
    let core_scores = group_coverage(evidence, &core.partition, &core_assembly.native, &groups);
    let expanded_scores = group_coverage(
        evidence,
        &expanded.partition,
        &expanded_assembly.native,
        &groups,
    );
    let take_expanded = groups
        .iter()
        .map(|&group| expanded_scores[group] > core_scores[group])
        .collect::<Vec<_>>();
    let mut differing = BTreeSet::new();
    for (support, &group) in groups.iter().enumerate() {
        if core.partition.region_by_support[support]
            != expanded.partition.region_by_support[support]
        {
            differing.insert(group);
        }
    }
    // Different carriers on the same exact membership also count.
    for (alternative, other) in [(core, expanded), (expanded, core)] {
        for (&id, proposal) in alternative
            .partition
            .selected_proposals
            .iter()
            .zip(&alternative.selected)
        {
            if other
                .partition
                .selected_proposals
                .binary_search(&id)
                .is_err()
            {
                differing.insert(groups[proposal.supports[0]]);
            }
        }
    }
    let compared_neighborhoods = differing.len();
    let expanded_neighborhoods = differing
        .iter()
        .filter(|&&group| take_expanded[group])
        .count();
    let (choice, expanded_neighborhoods) = if expanded_neighborhoods == 0 {
        (RawChoice::Core, 0)
    } else if expanded_neighborhoods == compared_neighborhoods {
        (RawChoice::Expanded, expanded_neighborhoods)
    } else {
        let combined = combine(evidence, core, expanded, &take_expanded)?;
        let assembly = assemble_native(evidence, combined.partition.clone(), options)?;
        // Reconciliation and native chart seams may affect other regions.
        // Validate the actual combined assembly against both candidates.
        let score = total_coverage(evidence, &combined.partition, &assembly.native);
        let core_score = total_coverage(evidence, &core.partition, &core_assembly.native);
        let expanded_score =
            total_coverage(evidence, &expanded.partition, &expanded_assembly.native);
        if score >= core_score.max(expanded_score) {
            (
                RawChoice::Combined(Box::new((combined, assembly))),
                expanded_neighborhoods,
            )
        } else if expanded_score > core_score {
            (RawChoice::Expanded, compared_neighborhoods)
        } else {
            (RawChoice::Core, 0)
        }
    };
    Ok(RawDecision {
        choice,
        compared_neighborhoods,
        expanded_neighborhoods,
    })
}
