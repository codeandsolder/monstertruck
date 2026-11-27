//! Failed curved/planar contact edges identify ownership alternatives. A
//! complete support may be absorbed only by a source-witnessed curved carrier;
//! no source label or endpoint misfit alone establishes a claim.
use crate::mesh_recovery::{
    MeshEvidence, RegionBoundaryIndex, SurfaceCarrier, SurfaceRegionPartition,
    analytic::{AnalyticRecoveryError, AnalyticRecoveryOptions, NativeAssembly},
    evidence::{OwnershipGrowthOptions, support_claim},
};
use std::collections::{BTreeMap, BTreeSet};

fn failed_planar_supports(evidence: &MeshEvidence, assembly: &NativeAssembly) -> BTreeSet<usize> {
    let partition = &assembly.reconciled.partition;
    let mut result = BTreeSet::new();
    for (trim, native) in assembly.trims.trims.iter().zip(&assembly.native.trims) {
        if native.is_ok() {
            continue;
        }
        let Some(right) = trim.regions.1 else {
            continue;
        };
        let Some(left_carrier) = partition.carriers[trim.regions.0] else {
            continue;
        };
        let Some(right_carrier) = partition.carriers[right] else {
            continue;
        };
        let left_planar = matches!(left_carrier, SurfaceCarrier::Plane { .. });
        if left_planar == matches!(right_carrier, SurfaceCarrier::Plane { .. }) {
            continue;
        }
        for edge in &trim.source_edges {
            let polygon = if left_planar {
                Some(edge.interior_polygon)
            } else {
                edge.exterior_polygon
            };
            if let Some(polygon) = polygon {
                result.insert(evidence.source_planar_support[polygon]);
            }
        }
    }
    result
}

fn best_claim(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    support: usize,
    policy: OwnershipGrowthOptions,
) -> Option<usize> {
    let mut claims = BTreeMap::<usize, f64>::new();
    for &interface in &evidence.support_interface_indices[support] {
        let [a, b] = evidence.support_interfaces[interface].supports;
        let neighbor = if a == support { b } else { a };
        let region = partition.region_by_support[neighbor];
        let carrier = partition.carriers[region]?;
        if matches!(carrier, SurfaceCarrier::Plane { .. }) {
            continue;
        }
        if let Some(score) = support_claim(evidence, support, neighbor, carrier, policy) {
            claims
                .entry(region)
                .and_modify(|previous| *previous = previous.min(score))
                .or_insert(score);
        }
    }
    let mut claims = claims.into_iter().collect::<Vec<_>>();
    claims.sort_unstable_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let &(region, score) = claims.first()?;
    if claims
        .get(1)
        .is_some_and(|&(_, second)| second - score <= crate::KERNEL_EPSILON * 0.01)
    {
        return None;
    }
    Some(region)
}

pub(super) fn expand(
    evidence: &MeshEvidence,
    original: &SurfaceRegionPartition,
    assembly: &NativeAssembly,
    regions: &BTreeSet<usize>,
    options: AnalyticRecoveryOptions,
) -> Result<Option<(SurfaceRegionPartition, BTreeSet<usize>, usize)>, AnalyticRecoveryError> {
    let policy = OwnershipGrowthOptions {
        maximum_adjustment: options.max_reconstruction_adjustment,
        declared_uncertainty: options.source_position_uncertainty,
        angular_tolerance: options.profile_normal_tolerance_radians,
        linear_chord_fraction: options.max_linear_relative_chord_error,
        circular_chord_fraction: options.max_circular_relative_chord_error,
    };
    let mut labels = original.region_by_support.clone();
    let mut targets = regions.clone();
    let mut added = 0;
    for support in failed_planar_supports(evidence, assembly) {
        let Some(target) = best_claim(evidence, original, support, policy) else {
            continue;
        };
        if labels[support] != target {
            labels[support] = target;
            targets.insert(target);
            added += 1;
        }
    }
    if added == 0 {
        return Ok(None);
    }
    let boundary = RegionBoundaryIndex::new(evidence)?;
    let mut partition =
        boundary.repartition_connected(&labels, |label| original.carriers[label])?;
    partition.exact_components = original.exact_components;
    partition.limited_components = original.limited_components;
    let regions = labels
        .iter()
        .zip(&partition.region_by_support)
        .filter_map(|(old, &new)| targets.contains(old).then_some(new))
        .collect();
    Ok(Some((partition, regions, added)))
}
