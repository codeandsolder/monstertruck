//! Bounded transactions on failed analytic contact neighborhoods. Source
//! ownership is unchanged, and successful analytic regions must stay valid.
use super::{
    AnalyticRecoveryError, AnalyticRecoveryOptions, NativeAssembly, NativeTrimOptions,
    ReconciledPartition, RegionBoundaryIndex,
    coplanar_contacts::{failed_edges, valid_regions},
};
use crate::{
    KERNEL_EPSILON,
    mesh_recovery::{
        MeshEvidence,
        evidence::{
            FrozenSourceFit, JointFitCandidate, JointFitContext, JointFitMethod, JointFitOptions,
            JointFitRejection, fit_neighborhood,
        },
    },
};
use std::collections::{BTreeMap, BTreeSet};

/// Accepted source-aware joint analytic fit. These are measured source
/// adjustments, not an inferred tolerance of the original CAD model.
#[derive(Clone, Debug)]
pub struct JointFitReport {
    /// Strategy that produced the accepted reconstruction.
    pub method: JointFitMethod,
    /// Region IDs participating in the common constraint system.
    pub regions: Vec<usize>,
    /// Number of fitted carrier parameter sets that changed.
    pub refitted_regions: usize,
    /// Original failed source edges repaired by this transaction.
    pub repaired_edges: usize,
    /// Additional original supports absorbed at witnessed failed contacts.
    pub ownership_additions: usize,
    /// Maximum actual displacement of an affected source vertex, including native assembly.
    pub max_displacement: f64,
    /// Largest local allowance supplied by sampled chord evidence or caller uncertainty.
    pub max_allowance: f64,
    /// Bound on added normalized carrier residual over all original source facets.
    pub max_additional_carrier_residual: f64,
}

/// One rejected trial, with the original region IDs and exact rejection stage.
#[derive(Clone, Debug)]
pub struct JointFitAttempt {
    /// Adjustable carrier region IDs.
    pub regions: Vec<usize>,
    /// Strategy attempted.
    pub method: JointFitMethod,
    /// Reason this trial was not committed.
    pub rejection: JointFitRejection,
}

fn failed_neighborhoods(assembly: &NativeAssembly) -> Vec<Vec<usize>> {
    let mut pairs = BTreeMap::<[usize; 2], usize>::new();
    for (trim, native) in assembly.trims.trims.iter().zip(&assembly.native.trims) {
        let Some(right) = trim.regions.1 else {
            continue;
        };
        if native.is_err() {
            let mut pair = [trim.regions.0, right];
            pair.sort_unstable();
            *pairs.entry(pair).or_default() += trim.source_edges.len();
        }
    }
    let mut weighted = pairs.into_iter().collect::<Vec<_>>();
    weighted.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut adjacency = BTreeMap::<usize, Vec<usize>>::new();
    for &(pair, _) in &weighted {
        adjacency.entry(pair[0]).or_default().push(pair[1]);
        adjacency.entry(pair[1]).or_default().push(pair[0]);
    }
    let mut result = Vec::new();
    let mut visited = BTreeSet::new();
    for &(pair, _) in &weighted {
        if !visited.insert(pair[0]) {
            continue;
        }
        let mut group = vec![pair[0]];
        let mut cursor = 0;
        while let Some(&region) = group.get(cursor) {
            for &neighbor in &adjacency[&region] {
                if visited.insert(neighbor) {
                    group.push(neighbor);
                }
            }
            cursor += 1;
        }
        if group.len() <= 6 {
            group.sort_unstable();
            result.push(group);
        }
    }
    // Large connected failure graphs and rejected whole-neighborhood fits
    // can still have independently repairable local contacts.
    for (pair, _) in weighted {
        let pair = pair.to_vec();
        if !result.contains(&pair) {
            result.push(pair);
        }
    }
    result
}

fn dependent_neighborhoods(
    assembly: &NativeAssembly,
    context: &JointFitContext,
) -> Vec<Vec<usize>> {
    let partition = &assembly.reconciled.partition;
    let failed = failed_edges(assembly);
    let mut stars = BTreeMap::<usize, BTreeSet<usize>>::new();
    for [a, b] in failed {
        for id in [a, b] {
            stars
                .entry(id)
                .or_insert_with(|| context.dependent_at(id, partition).collect());
        }
    }
    // Merge stars through their shared dependent planes, not through
    // independent curved carriers that may span a whole component.
    let mut adjacency = BTreeMap::<usize, BTreeSet<usize>>::new();
    let mut frequency = BTreeMap::<usize, usize>::new();
    for star in stars.values() {
        let Some(&first) = star.first() else {
            continue;
        };
        for &id in star {
            adjacency.entry(first).or_default().insert(id);
            adjacency.entry(id).or_default().insert(first);
            *frequency.entry(id).or_default() += 1;
        }
    }
    let mut groups = Vec::new();
    let mut seen = BTreeSet::new();
    for &seed in adjacency.keys() {
        if !seen.insert(seed) {
            continue;
        }
        let mut group = vec![seed];
        let mut cursor = 0;
        while let Some(&id) = group.get(cursor) {
            for &neighbor in &adjacency[&id] {
                if seen.insert(neighbor) {
                    group.push(neighbor);
                }
            }
            cursor += 1;
        }
        group.sort_unstable();
        let score = group.iter().map(|id| frequency[id]).sum::<usize>();
        groups.push((score, group));
    }
    groups.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut result = groups
        .into_iter()
        .map(|(_, group)| group)
        .collect::<Vec<_>>();
    // A wide dependent component may contain an unrepairable junction.
    // Retry vertex-local stars instead of repeatedly selecting the same
    // large independent carrier.
    let mut local = stars
        .into_values()
        .filter(|star| !star.is_empty())
        .map(|star| star.into_iter().collect::<Vec<_>>())
        .collect::<Vec<_>>();
    local.sort_unstable_by(|a, b| b.len().cmp(&a.len()).then(a.cmp(b)));
    local.dedup();
    for star in local {
        if !result.contains(&star) {
            result.push(star);
        }
    }
    result
}

pub(super) fn native_displacement(
    evidence: &MeshEvidence,
    assembly: &NativeAssembly,
    candidate: &JointFitCandidate,
    native: &crate::mesh_recovery::NativeTrimGraph,
    numerical_allowance: f64,
) -> Result<f64, JointFitRejection> {
    let mut maximum = candidate.maximum_displacement;
    for (id, (vertex, old)) in native
        .vertices
        .iter()
        .zip(&assembly.native.vertices)
        .enumerate()
    {
        let Some(vertex) = vertex else {
            continue;
        };
        let source = evidence.positions[id];
        let displacement = (vertex.point() - source).norm();
        let allowed = if candidate.allowances[id] > 0.0 {
            maximum = maximum.max(displacement);
            candidate.allowances[id]
        } else {
            // No fresh approximation allowance for unrelated vertices.
            let old = old.as_ref().map_or(
                assembly.reconciled.positions[id],
                monstertruck_topology::Vertex::point,
            );
            (old - source).norm() + numerical_allowance
        };
        if !displacement.is_finite() || displacement > allowed {
            return Err(JointFitRejection::NativeDisplacement {
                vertex: id,
                allowance: allowed,
            });
        }
    }
    Ok(maximum)
}

// A full native rebuild is only useful if at least one previously
// failed source edge can now have an accepted trim. Check affected paths
// only; a newly unresolved edge would necessarily reject the transaction.
fn can_improve_trims(
    assembly: &NativeAssembly,
    candidate: &JointFitCandidate,
    regions: &[usize],
    options: super::SharedTrimOptions,
) -> bool {
    let failed = failed_edges(assembly);
    let mut possible = false;
    for path in &candidate.partition.contacts.paths {
        if !regions.contains(&path.regions.0)
            && !path.regions.1.is_some_and(|id| regions.contains(&id))
            && !path
                .edges
                .iter()
                .any(|e| e.vertices.iter().any(|&id| candidate.allowances[id] > 0.0))
        {
            continue;
        }
        let carriers = std::iter::once(path.regions.0)
            .chain(path.regions.1)
            .filter_map(|id| candidate.partition.carriers[id])
            .collect::<Vec<_>>();
        let Some(segmentation) =
            path.segment_analytic_trims(&candidate.positions, &carriers, options)
        else {
            continue;
        };
        for piece in segmentation.pieces {
            for offset in 0..piece.edge_count {
                let mut edge = path.edges[(piece.first_edge + offset) % path.edges.len()].vertices;
                edge.sort_unstable();
                if matches!(
                    piece.geometry,
                    crate::mesh_recovery::SharedTrimGeometry::Unresolved
                ) {
                    if !failed.contains(&edge) {
                        return false;
                    }
                } else if failed.contains(&edge) {
                    possible = true;
                }
            }
        }
    }
    possible
}

fn try_fit(
    evidence: &MeshEvidence,
    assembly: &mut NativeAssembly,
    regions: &[usize],
    options: AnalyticRecoveryOptions,
    context: &mut JointFitContext,
    method: JointFitMethod,
) -> Result<Result<JointFitReport, JointFitRejection>, AnalyticRecoveryError> {
    let candidate = match fit_neighborhood(
        evidence,
        &assembly.reconciled.partition,
        &assembly.reconciled.positions,
        regions,
        context,
        method,
        Some(&assembly.native),
    ) {
        Ok(candidate) => candidate,
        Err(reason) => return Ok(Err(reason)),
    };
    if candidate.refitted_regions == 0 && candidate.positions == assembly.reconciled.positions {
        return Ok(Err(JointFitRejection::Unchanged));
    }
    let boundary = RegionBoundaryIndex::new(evidence)?;
    let trim = options.trim_options(evidence);
    if !can_improve_trims(assembly, &candidate, regions, trim) {
        return Ok(Err(JointFitRejection::TrimPreflight));
    }
    let trims = boundary.assemble_shared_trims(&candidate.partition, &candidate.positions, trim)?;
    let native = boundary.construct_native_trims(
        &candidate.partition,
        &trims,
        &candidate.positions,
        NativeTrimOptions {
            max_vertex_displacement: trim.absolute_tolerance,
            intersection_tolerance: KERNEL_EPSILON * 0.1,
            angular_resolution: options.axis_relative_resolution,
        },
    )?;
    let max_displacement = match native_displacement(
        evidence,
        assembly,
        &candidate,
        &native,
        trim.absolute_tolerance,
    ) {
        Ok(displacement) => displacement,
        Err(reason) => return Ok(Err(reason)),
    };
    let refitted_regions = candidate.refitted_regions;
    let max_additional_carrier_residual = candidate.max_additional_carrier_residual;
    let max_allowance = candidate.allowances.iter().copied().fold(0.0_f64, f64::max);
    let mut next = NativeAssembly {
        reconciled: ReconciledPartition {
            partition: candidate.partition,
            positions: candidate.positions,
            planes: Vec::new(),
            max_displacement: candidate
                .maximum_displacement
                .max(assembly.reconciled.max_displacement),
            rejected: assembly.reconciled.rejected,
        },
        trims,
        native,
    };
    let before = failed_edges(assembly);
    let after = failed_edges(&next);
    if after.len() >= before.len() {
        return Ok(Err(JointFitRejection::NoImprovement));
    }
    if let Some(&edge) = after.difference(&before).next() {
        return Ok(Err(JointFitRejection::NewFailedEdge(edge)));
    }
    if next.trims.ambiguous_edges != assembly.trims.ambiguous_edges {
        return Ok(Err(JointFitRejection::AmbiguousIncidence));
    }
    if let Some(&region) = valid_regions(assembly)
        .difference(&valid_regions(&next))
        .next()
    {
        let error = next
            .native
            .regions
            .iter()
            .find(|r| r.region == region)
            .and_then(|r| r.faces.as_ref().err().copied());
        return Ok(Err(JointFitRejection::InvalidatedRegion { region, error }));
    }
    // Existing plane records describe prior accepted transactions; keep
    // them as provenance, not a claim that the carriers were never refit.
    next.reconciled.planes = std::mem::take(&mut assembly.reconciled.planes);
    context.commit(&candidate.carrier_changes);
    *assembly = next;
    Ok(Ok(JointFitReport {
        method,
        regions: regions.to_vec(),
        refitted_regions,
        repaired_edges: before.len() - after.len(),
        ownership_additions: 0,
        max_displacement,
        max_allowance,
        max_additional_carrier_residual,
    }))
}

pub(super) fn repair(
    evidence: &MeshEvidence,
    assembly: &mut NativeAssembly,
    options: AnalyticRecoveryOptions,
    frozen: Option<&FrozenSourceFit>,
) -> Result<(Vec<JointFitReport>, Vec<JointFitAttempt>), AnalyticRecoveryError> {
    let mut reports = Vec::new();
    let mut rejected = Vec::new();
    if options.max_reconstruction_adjustment == 0.0 || options.max_joint_fit_trials == 0 {
        return Ok((reports, rejected));
    }
    if failed_neighborhoods(assembly).is_empty() {
        return Ok((reports, rejected));
    }
    let policy = JointFitOptions {
        maximum_adjustment: options.max_reconstruction_adjustment,
        declared_uncertainty: options.source_position_uncertainty,
        angular_resolution: options.axis_relative_resolution,
    };
    let partition = &assembly.reconciled.partition;
    let context = frozen.map_or_else(
        || JointFitContext::new(evidence, partition, policy),
        |fit| fit.context(evidence, partition, policy),
    );
    let Some(mut context) = context else {
        return Ok((reports, rejected));
    };
    let mut trials = 0;
    while trials < options.max_joint_fit_trials {
        let neighborhoods = failed_neighborhoods(assembly);
        let fillets = neighborhoods
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .filter(|&id| context.anchored_fillet_possible(&assembly.reconciled.partition, id));
        let mut plans = fillets
            .map(|id| (vec![id], &[JointFitMethod::AnchoredFillet][..]))
            .collect::<Vec<_>>();
        let general = &[
            JointFitMethod::FixedCarriers,
            JointFitMethod::DependentPlanes,
            JointFitMethod::JointRefit,
        ][..];
        plans.extend(
            neighborhoods
                .iter()
                .filter(|r| r.len() > 2)
                .map(|r| (r.clone(), general)),
        );
        plans.extend(
            dependent_neighborhoods(assembly, &context)
                .into_iter()
                .map(|r| (r, &[JointFitMethod::DependentPlanes][..])),
        );
        plans.extend(
            neighborhoods
                .into_iter()
                .filter(|r| r.len() <= 2)
                .map(|r| (r, general)),
        );
        let mut improved = false;
        for (regions, methods) in plans {
            if trials >= options.max_joint_fit_trials {
                break;
            }
            trials += 1;
            for &method in methods {
                match try_fit(evidence, assembly, &regions, options, &mut context, method)? {
                    Ok(report) => {
                        reports.push(report);
                        improved = true;
                        break;
                    }
                    Err(rejection) => rejected.push(JointFitAttempt {
                        regions: regions.clone(),
                        method,
                        rejection,
                    }),
                }
            }
            if improved {
                break;
            }
        }
        if !improved {
            break;
        }
    }
    Ok((reports, rejected))
}
