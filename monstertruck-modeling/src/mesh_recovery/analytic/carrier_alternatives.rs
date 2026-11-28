//! Preserve competing geometry until actual shared native contacts can
//! discriminate between source-supported carriers on the same ownership.

use super::{
    AnalyticRecoveryError, AnalyticRecoveryOptions, NativeAssembly, SelectedSurfaceProposal,
    assemble_native,
    native_selection::{PartitionAlternative, accept_repair},
};
use crate::mesh_recovery::{MeshEvidence, SurfaceRegionProposal};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug)]
struct IndexedCarrier {
    id: usize,
    proposal: SelectedSurfaceProposal,
}

#[derive(Debug)]
pub(super) struct CarrierAlternatives {
    by_support: BTreeMap<Vec<usize>, Vec<IndexedCarrier>>,
}

#[derive(Default)]
pub(super) struct CarrierTrials {
    pub trials: usize,
    pub accepted: usize,
}

impl CarrierAlternatives {
    pub(super) fn new(proposals: &[SurfaceRegionProposal<'_>]) -> Self {
        let mut membership = BTreeMap::<&[usize], Vec<usize>>::new();
        for (id, p) in proposals.iter().enumerate() {
            membership.entry(p.supports()).or_default().push(id);
        }
        let by_support = membership
            .into_iter()
            .filter(|(_, ids)| ids.len() > 1)
            .map(|(supports, ids)| {
                let mut candidates = ids
                    .into_iter()
                    .map(|id| IndexedCarrier {
                        id,
                        proposal: SelectedSurfaceProposal::from_measured(&proposals[id]),
                    })
                    .collect::<Vec<_>>();
                candidates.sort_unstable_by(|a, b| {
                    b.proposal
                        .independently_tested_profile
                        .cmp(&a.proposal.independently_tested_profile)
                        .then_with(|| {
                            a.proposal
                                .relative_vertex_error
                                .total_cmp(&b.proposal.relative_vertex_error)
                        })
                        .then_with(|| {
                            a.proposal
                                .relative_chord_bound
                                .total_cmp(&b.proposal.relative_chord_bound)
                        })
                        .then(a.id.cmp(&b.id))
                });
                (supports.to_vec(), candidates)
            })
            .collect();
        Self { by_support }
    }

    pub(super) fn repair(
        &self,
        evidence: &MeshEvidence,
        alternative: &mut PartitionAlternative,
        assembly: &mut NativeAssembly,
        options: &AnalyticRecoveryOptions,
    ) -> Result<CarrierTrials, AnalyticRecoveryError> {
        let mut result = CarrierTrials::default();
        if self.by_support.is_empty() || options.max_carrier_trials == 0 {
            return Ok(result);
        }
        let mut rejected = BTreeSet::new();
        while result.trials < options.max_carrier_trials {
            let implicated = implicated_regions(assembly);
            let mut candidates = alternative
                .partition
                .selected_proposals
                .iter()
                .zip(&alternative.selected)
                .enumerate()
                .filter_map(|(slot, (&current, proposal))| {
                    let region =
                        assembly.reconciled.partition.region_by_support[proposal.supports[0]];
                    let weight = *implicated.get(&region)?;
                    let variants = self.by_support.get(&proposal.supports)?;
                    Some(
                        variants
                            .iter()
                            .filter(move |variant| {
                                variant.id != current
                                    && (!proposal.independently_tested_profile
                                        || variant.proposal.independently_tested_profile)
                            })
                            .map(move |variant| (std::cmp::Reverse(weight), slot, variant)),
                    )
                })
                .flatten()
                .filter(|(_, _, variant)| !rejected.contains(&variant.id))
                .collect::<Vec<_>>();
            candidates.sort_by_key(|(weight, slot, _)| (*weight, *slot));
            let Some((_, slot, variant)) = candidates.first().copied() else {
                break;
            };
            result.trials += 1;
            let proposed = replace_carrier(alternative, slot, variant);
            let candidate = match assemble_native(evidence, proposed.partition.clone(), options) {
                Ok(candidate) => candidate,
                Err(
                    AnalyticRecoveryError::Boundary(_)
                    | AnalyticRecoveryError::Native(_)
                    | AnalyticRecoveryError::SourceChordBudget(_),
                ) => {
                    rejected.insert(variant.id);
                    continue;
                }
                Err(error) => return Err(error),
            };
            if accept_repair(evidence, assembly, &candidate).is_ok() {
                *alternative = proposed;
                *assembly = candidate;
                result.accepted += 1;
                rejected.clear();
            } else {
                rejected.insert(variant.id);
            }
        }
        Ok(result)
    }
}

fn implicated_regions(assembly: &NativeAssembly) -> BTreeMap<usize, usize> {
    let invalid = assembly
        .native
        .regions
        .iter()
        .filter(|region| region.faces.is_err() || region.requires_chart_cut)
        .map(|region| region.region)
        .collect::<BTreeSet<_>>();
    let mut weights = invalid
        .iter()
        .map(|&id| (id, 1_usize))
        .collect::<BTreeMap<_, _>>();
    for (trim, native) in assembly.trims.trims.iter().zip(&assembly.native.trims) {
        let incident = std::iter::once(trim.regions.0).chain(trim.regions.1);
        if native.is_err() || incident.clone().any(|region| invalid.contains(&region)) {
            for id in incident {
                *weights.entry(id).or_default() += trim.source_edges.len();
            }
        }
    }
    weights
}

fn replace_carrier(
    original: &PartitionAlternative,
    slot: usize,
    variant: &IndexedCarrier,
) -> PartitionAlternative {
    let mut candidate = original.clone();
    let region = candidate.partition.region_by_support[variant.proposal.supports[0]];
    candidate.partition.carriers[region] = Some(variant.proposal.carrier);
    candidate.partition.selected_proposals[slot] = variant.id;
    candidate.selected[slot] = variant.proposal.clone();
    // Selection IDs must stay sorted and paired with the corresponding
    // evidence records, including after a same-ownership replacement.
    let mut pairs = candidate
        .partition
        .selected_proposals
        .into_iter()
        .zip(candidate.selected)
        .collect::<Vec<_>>();
    pairs.sort_unstable_by_key(|(id, _)| *id);
    (candidate.partition.selected_proposals, candidate.selected) = pairs.into_iter().unzip();
    candidate
}

#[cfg(test)]
mod tests {
    use super::super::{
        native_complete_shell, recover_analytic_mesh, tests::sampled_cylinder_mesh,
    };
    use super::*;
    use crate::{
        Point3, Vector3,
        mesh_recovery::{
            RegionBoundaryIndex, SurfaceCarrier, SurfaceRegionBudget,
            carrier_identity::{GeometryBuckets, GeometryCandidate},
        },
    };

    fn cylinder(radius: f64, axial_offset: f64, sign: f64) -> SurfaceCarrier {
        SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, axial_offset),
            axis: Vector3::new(0.0, 0.0, sign),
            radius_at_origin: radius,
            slope: 0.0,
        }
    }

    struct Hypothesis {
        supports: Vec<usize>,
        carrier: SurfaceCarrier,
        error: f64,
    }

    impl GeometryCandidate for Hypothesis {
        fn supports(&self) -> &[usize] { &self.supports }
        fn carrier(&self) -> SurfaceCarrier { self.carrier }
        fn error(&self) -> f64 { self.error }
    }

    #[test]
    fn only_geometrically_equivalent_candidates_share_a_bucket()
    -> Result<(), Box<dyn std::error::Error>> {
        let recovered = recover_analytic_mesh(
            &sampled_cylinder_mesh()?,
            &AnalyticRecoveryOptions::default(),
        )?;
        let supports = &recovered.discovery.selected[0].supports;
        let mut buckets = GeometryBuckets::default();
        for (carrier, error) in [
            (cylinder(1.0, 0.0, 1.0), 0.2),
            (cylinder(1.0, 2.0, -1.0), 0.1),
            (cylinder(1.001, 0.0, 1.0), 0.3),
        ] {
            buckets.insert(
                Hypothesis {
                    supports: supports.clone(),
                    carrier,
                    error,
                },
                &recovered.evidence,
                1.0e-8,
            );
        }
        let candidates = buckets.into_candidates();
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].error, 0.1);
        assert_eq!(candidates[1].error, 0.3);
        Ok(())
    }

    #[test]
    fn native_contacts_discriminate_same_ownership_carriers()
    -> Result<(), Box<dyn std::error::Error>> {
        let options = AnalyticRecoveryOptions::default();
        let recovered = recover_analytic_mesh(&sampled_cylinder_mesh()?, &options)?;
        let evidence = &recovered.evidence;
        let supports = &recovered.discovery.selected[0].supports;
        let boundary = RegionBoundaryIndex::new(evidence)?;
        let budget = SurfaceRegionBudget {
            vertex: 0.002,
            chord: 0.02,
        };
        let proposals = [
            boundary
                .measure_surface_region(supports, cylinder(1.001, 0.0, 1.0), budget)
                .ok_or("offset cylinder must be source admissible")?,
            boundary
                .measure_surface_region(supports, cylinder(1.0, 0.0, 1.0), budget)
                .ok_or("correct cylinder must be source admissible")?,
        ];
        let partition = boundary.select_surface_partition(&proposals[..1], evidence, 100)?;
        let ownership = partition.region_by_support.clone();
        let mut alternative = PartitionAlternative {
            selected: vec![SelectedSurfaceProposal::from_measured(&proposals[0])],
            partition,
        };
        let mut assembly = assemble_native(evidence, alternative.partition.clone(), &options)?;
        assert!(native_complete_shell(&assembly.trims, &assembly.native).is_none());
        let alternatives = CarrierAlternatives::new(&proposals);
        let mut disabled = options.clone();
        disabled.max_carrier_trials = 0;
        assert_eq!(
            alternatives
                .repair(evidence, &mut alternative, &mut assembly, &disabled)?
                .trials,
            0
        );
        let trials = alternatives.repair(evidence, &mut alternative, &mut assembly, &options)?;
        assert_eq!(trials.trials, 1);
        assert_eq!(trials.accepted, 1);
        assert_eq!(alternative.partition.selected_proposals, vec![1]);
        assert_eq!(alternative.partition.region_by_support, ownership);
        assert!(native_complete_shell(&assembly.trims, &assembly.native).is_some());
        Ok(())
    }
}
