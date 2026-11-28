//! Treat an expanded ownership hypothesis and its displaced source junctions
//! as one transaction, rather than ranking the unfinished native assembly.
use super::{PartitionAlternative, feedback, total_coverage};
use crate::{
    KERNEL_EPSILON,
    mesh_recovery::{
        MeshEvidence, NativeTrimOptions, ReconciledPartition, RegionBoundaryIndex, SurfaceCarrier,
        SurfaceRegionPartition,
        analytic::{
            AnalyticRecoveryError, AnalyticRecoveryOptions, JointFitAttempt, JointFitReport,
            NativeAssembly,
            coplanar_contacts::{failed_edges, valid_regions},
            joint_recovery::native_displacement,
            native_complete_shell,
        },
        evidence::{
            FrozenSourceFit, JointFitMethod, JointFitOptions, JointFitRejection,
            fit_ownership_junctions,
        },
    },
};
use std::collections::BTreeSet;

pub(in crate::mesh_recovery::analytic) struct OwnershipTrial {
    pub report: JointFitReport,
    pub frozen: FrozenSourceFit,
}

pub(super) struct RejectedOwnership {
    pub attempt: JointFitAttempt,
    pub pending: Option<Box<NativeAssembly>>,
}

type TrialResult =
    Result<(PartitionAlternative, NativeAssembly, OwnershipTrial), RejectedOwnership>;

fn improves_contacts(baseline: &NativeAssembly, candidate: &NativeAssembly) -> bool {
    let before = failed_edges(baseline);
    let after = failed_edges(candidate);
    before.len() > after.len()
        && after.is_subset(&before)
        && baseline.trims.ambiguous_edges == candidate.trims.ambiguous_edges
}

fn quality(
    evidence: &MeshEvidence,
    assembly: &NativeAssembly,
) -> (std::cmp::Reverse<usize>, super::Coverage) {
    (
        std::cmp::Reverse(failed_edges(assembly).len()),
        total_coverage(evidence, &assembly.reconciled.partition, &assembly.native),
    )
}

pub(super) const fn rejection(
    regions: Vec<usize>,
    rejection: JointFitRejection,
) -> JointFitAttempt {
    JointFitAttempt {
        regions,
        method: JointFitMethod::OwnershipJunctions,
        rejection,
    }
}

fn changed_regions(
    core: &PartitionAlternative,
    expanded: &PartitionAlternative,
) -> BTreeSet<usize> {
    expanded
        .partition
        .selected_proposals
        .iter()
        .zip(&expanded.selected)
        .filter(|&(id, p)| {
            !core.partition.selected_proposals.contains(id)
                && !matches!(p.carrier, SurfaceCarrier::Plane { .. })
        })
        .map(|(_, p)| expanded.partition.region_by_support[p.supports[0]])
        .collect()
}

const fn speculative_error(
    error: AnalyticRecoveryError,
) -> Result<JointFitRejection, AnalyticRecoveryError> {
    match error {
        AnalyticRecoveryError::Boundary(error) => Ok(JointFitRejection::CandidateBoundary(error)),
        AnalyticRecoveryError::Native(error) => Ok(JointFitRejection::CandidateNative(error)),
        AnalyticRecoveryError::SourceChordBudget(support) => {
            Ok(JointFitRejection::SourceChordBudget(support))
        }
        other => Err(other),
    }
}

pub(super) fn attempt(
    evidence: &MeshEvidence,
    core: &PartitionAlternative,
    expanded: &PartitionAlternative,
    raw: &NativeAssembly,
    baseline: &NativeAssembly,
    options: &AnalyticRecoveryOptions,
) -> Result<Option<TrialResult>, AnalyticRecoveryError> {
    if options.max_joint_fit_trials == 0
        || options.max_reconstruction_adjustment == 0.0
        || failed_edges(raw).is_empty()
    {
        return Ok(None);
    }
    let mut regions = changed_regions(core, expanded);
    if regions.is_empty() {
        return Ok(None);
    }
    let mut partition = raw.reconciled.partition.clone();
    let mut best = None::<(PartitionAlternative, NativeAssembly, OwnershipTrial)>;
    let mut last_rejection = None;
    let mut pending = None::<Box<NativeAssembly>>;
    let mut additions = 0;
    for _ in 0..options.max_joint_fit_trials {
        let (assembly, mut fit) = match fit_partition(evidence, raw, &partition, &regions, options)
        {
            Ok(Ok(value)) => value,
            Ok(Err(reason)) => {
                last_rejection = Some(reason);
                break;
            }
            Err(error) => {
                last_rejection = Some(rejection(
                    regions.iter().copied().collect(),
                    speculative_error(error)?,
                ));
                break;
            }
        };
        fit.report.ownership_additions = additions;
        let next = match feedback::expand(evidence, &partition, &assembly, &regions, options) {
            Ok(next) => next,
            Err(error) => {
                last_rejection = Some(rejection(
                    regions.iter().copied().collect(),
                    speculative_error(error)?,
                ));
                None
            }
        };
        match accept(evidence, baseline, &assembly) {
            Ok(repaired_edges) => {
                fit.report.repaired_edges = repaired_edges;
                if best.as_ref().is_none_or(|(_, previous, _)| {
                    let candidate = quality(evidence, &assembly);
                    let old = quality(evidence, previous);
                    (candidate.1, candidate.0) > (old.1, old.0)
                }) {
                    let mut alternative = expanded.clone();
                    alternative.partition = partition.clone();
                    best = Some((alternative, assembly, fit));
                }
            }
            Err(reason) => {
                last_rejection = Some(rejection(fit.report.regions, reason));
                // A failed strict chart is not a reason to discard a
                // source-bounded candidate with strictly better contacts.
                // Keep it separately for an explicitly budgeted hybrid trial.
                if improves_contacts(baseline, &assembly)
                    && pending
                        .as_ref()
                        .is_none_or(|old| quality(evidence, &assembly) > quality(evidence, old))
                {
                    pending = Some(Box::new(assembly));
                }
            }
        }
        let Some((next_partition, next_regions, added)) = next else {
            break;
        };
        partition = next_partition;
        regions = next_regions;
        additions += added;
    }
    Ok(best
        .map(Ok)
        .or_else(|| last_rejection.map(|attempt| Err(RejectedOwnership { attempt, pending }))))
}

fn fit_partition(
    evidence: &MeshEvidence,
    raw: &NativeAssembly,
    partition: &SurfaceRegionPartition,
    regions: &BTreeSet<usize>,
    options: &AnalyticRecoveryOptions,
) -> Result<Result<(NativeAssembly, OwnershipTrial), JointFitAttempt>, AnalyticRecoveryError> {
    let policy = JointFitOptions {
        maximum_adjustment: options.max_reconstruction_adjustment,
        declared_uncertainty: options.source_position_uncertainty,
        angular_resolution: options.axis_relative_resolution,
    };
    let (candidate, frozen) = match fit_ownership_junctions(
        evidence,
        partition,
        &raw.reconciled.positions,
        regions,
        policy,
    ) {
        Ok(candidate) => candidate,
        Err(reason) => return Ok(Err(rejection(regions.iter().copied().collect(), reason))),
    };
    let trim = options.trim_options(evidence);
    let boundary = RegionBoundaryIndex::new(evidence)?;
    if let Some(support) = options.source_chord_violation(&boundary, &candidate.partition) {
        return Ok(Err(rejection(
            regions.iter().copied().collect(),
            JointFitRejection::SourceChordBudget(support),
        )));
    }
    let mut trims =
        boundary.assemble_shared_trims(&candidate.partition, &candidate.positions, trim)?;
    let native_options = NativeTrimOptions {
        max_vertex_displacement: trim.absolute_tolerance,
        intersection_tolerance: KERNEL_EPSILON * 0.1,
        angular_resolution: options.axis_relative_resolution,
    };
    let mut native = boundary.construct_native_trims(
        &candidate.partition,
        &trims,
        &candidate.positions,
        native_options,
    )?;
    boundary.repair_planar_crossings(
        &candidate.partition,
        &candidate.positions,
        &mut trims,
        &mut native,
        native_options,
    )?;
    let max_displacement =
        match native_displacement(evidence, raw, &candidate, &native, trim.absolute_tolerance) {
            Ok(displacement) => displacement,
            Err(reason) => return Ok(Err(rejection(regions.iter().copied().collect(), reason))),
        };
    let report = JointFitReport {
        method: JointFitMethod::OwnershipJunctions,
        regions: regions
            .iter()
            .copied()
            .chain(candidate.carrier_changes.iter().map(|&(id, _)| id))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect(),
        refitted_regions: candidate.refitted_regions,
        repaired_edges: 0,
        ownership_additions: 0,
        max_displacement,
        max_allowance: candidate.allowances.iter().copied().fold(0.0, f64::max),
        max_additional_carrier_residual: candidate.max_additional_carrier_residual,
    };
    let assembly = NativeAssembly {
        reconciled: ReconciledPartition {
            partition: candidate.partition,
            positions: candidate.positions,
            planes: raw
                .reconciled
                .planes
                .iter()
                .filter_map(|plane| {
                    let region = partition.region_by_support[plane.supports[0]];
                    if !matches!(
                        partition.carriers[region],
                        Some(SurfaceCarrier::Plane { .. })
                    ) || plane
                        .supports
                        .iter()
                        .any(|&support| partition.region_by_support[support] != region)
                    {
                        return None;
                    }
                    let mut record = plane.clone();
                    record.region = region;
                    Some(record)
                })
                .collect(),
            max_displacement: raw
                .reconciled
                .max_displacement
                .max(candidate.maximum_displacement),
            rejected: raw.reconciled.rejected,
        },
        trims,
        native,
    };
    Ok(Ok((assembly, OwnershipTrial { report, frozen })))
}

pub(super) fn accept(
    evidence: &MeshEvidence,
    baseline: &NativeAssembly,
    candidate: &NativeAssembly,
) -> Result<usize, JointFitRejection> {
    let before = failed_edges(baseline);
    let after = failed_edges(candidate);
    if let Some(&edge) = after.difference(&before).next() {
        return Err(JointFitRejection::NewFailedEdge(edge));
    }
    if baseline.trims.ambiguous_edges != candidate.trims.ambiguous_edges {
        return Err(JointFitRejection::AmbiguousIncidence);
    }
    let old = &baseline.reconciled.partition;
    let new = &candidate.reconciled.partition;
    let old_valid = valid_regions(baseline);
    let new_valid = valid_regions(candidate);
    for (&a, &b) in old.region_by_support.iter().zip(&new.region_by_support) {
        if old_valid.contains(&a) && !new_valid.contains(&b) {
            let error = candidate
                .native
                .regions
                .iter()
                .find(|r| r.region == b)
                .and_then(|r| r.faces.as_ref().err().copied());
            return Err(JointFitRejection::InvalidatedRegion { region: b, error });
        }
    }
    if after.is_empty()
        && candidate.trims.ambiguous_edges.is_empty()
        && new_valid.len() == candidate.native.regions.len()
        && native_complete_shell(&candidate.trims, &candidate.native).is_none()
    {
        return Err(JointFitRejection::InconsistentShell);
    }
    if total_coverage(evidence, new, &candidate.native)
        <= total_coverage(evidence, old, &baseline.native)
        && after.len() == before.len()
    {
        return Err(JointFitRejection::NoImprovement);
    }
    Ok(before.len() - after.len())
}
