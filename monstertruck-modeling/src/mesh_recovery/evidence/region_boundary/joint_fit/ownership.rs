//! A competing partition can recognize source facets correctly but initially
//! fail native incidence. Reconcile only its displaced junctions, eliminating
//! incident triangular interpolation planes instead of treating them as
//! independent constraints. Source curvature evidence is frozen before fitting.
use super::{
    JointFitCandidate, JointFitContext, JointFitOptions, JointFitRejection, MeshEvidence,
    Neighborhood, Point3, SurfaceCarrier, SurfaceRegionPartition, fixed_sample_position,
    plane_from_junctions, preserves_polygon_orientation, satisfies,
};
use crate::KERNEL_EPSILON;
use std::collections::BTreeSet;
mod growth;
pub(in crate::mesh_recovery) use growth::{OwnershipGrowthOptions, support_claim};

// A witnessed nonlinear departure from endpoint interpolation. An endpoint
// misfit by itself cannot create its own displacement budget. Sample several
// interior locations so an inflection at the midpoint does not erase evidence.
fn chord_uncertainty(carrier: SurfaceCarrier, p: Point3, q: Point3) -> f64 {
    if matches!(carrier, SurfaceCarrier::Plane { .. }) {
        return 0.0;
    }
    let (Some((_, a)), Some((_, b))) = (carrier.signed_constraint(p), carrier.signed_constraint(q))
    else {
        return 0.0;
    };
    [0.25, 0.5, 0.75]
        .into_iter()
        .filter_map(|t| {
            let (_, residual) = carrier.signed_constraint(p + (q - p) * t)?;
            Some((residual - (b - a).mul_add(t, a)).abs())
        })
        .fold(0.0, f64::max)
}

/// The same source allowances survive subsequent joint fits. Region changes
/// are indexed by an original support, so a later exact ownership correction
/// cannot replenish the consumed part of a carrier budget.
pub(in crate::mesh_recovery) struct FrozenSourceFit {
    allowances: Vec<f64>,
    changes: Vec<(usize, f64)>,
}

impl FrozenSourceFit {
    pub(in crate::mesh_recovery) fn context(
        &self,
        evidence: &MeshEvidence,
        partition: &SurfaceRegionPartition,
        options: JointFitOptions,
    ) -> Option<JointFitContext> {
        let mut context =
            JointFitContext::with_chord_witness(evidence, partition, options, |_, _, _| 0.0)?;
        if context.samples.len() != self.allowances.len() {
            return None;
        }
        for (sample, &allowance) in context.samples.iter_mut().zip(&self.allowances) {
            sample.allowance = allowance;
        }
        for &(support, change) in &self.changes {
            context.carrier_change[*partition.region_by_support.get(support)?] += change;
        }
        Some(context)
    }
}

struct Junctions {
    positions: Vec<Point3>,
    affected: Vec<usize>,
    dependent: BTreeSet<usize>,
    allowances: Vec<f64>,
    max_displacement: f64,
}

fn place_junctions(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    prior: &[Point3],
    regions: &BTreeSet<usize>,
    context: &JointFitContext,
) -> Result<Junctions, JointFitRejection> {
    let mut result = Junctions {
        positions: prior.to_vec(),
        affected: Vec::new(),
        dependent: BTreeSet::new(),
        allowances: vec![0.0; prior.len()],
        max_displacement: 0.0,
    };
    for sample in &context.samples {
        if !sample.boundary || !sample.regions.iter().any(|id| regions.contains(id)) {
            continue;
        }
        let point = prior[sample.id];
        let needs_placement = sample.regions.iter().any(|&id| {
            regions.contains(&id)
                && partition.carriers[id].is_some_and(|carrier| {
                    !matches!(carrier, SurfaceCarrier::Plane { .. })
                        && !satisfies(point, &[carrier], KERNEL_EPSILON * 0.01)
                })
        });
        if !needs_placement {
            continue;
        }
        let dependent = context
            .dependent_at(sample.id, partition)
            .collect::<BTreeSet<_>>();
        let fixed = sample
            .regions
            .iter()
            .filter(|id| !dependent.contains(id))
            .map(|&id| partition.carriers[id])
            .collect::<Option<Vec<_>>>()
            .ok_or(JointFitRejection::Neighborhood)?;
        let source = evidence.positions[sample.id];
        let next = fixed_sample_position(
            source,
            point,
            None,
            &fixed,
            sample.allowance,
            context.options.angular_resolution,
        )
        .ok_or(JointFitRejection::Sample {
            vertex: sample.id,
            allowance: sample.allowance,
        })?;
        let displacement = (next - source).norm();
        if !displacement.is_finite() || displacement > sample.allowance {
            return Err(JointFitRejection::Sample {
                vertex: sample.id,
                allowance: sample.allowance,
            });
        }
        result.positions[sample.id] = next;
        result.affected.push(sample.id);
        result.dependent.extend(dependent);
        result.allowances[sample.id] = sample.allowance;
        result.max_displacement = result.max_displacement.max(displacement);
    }
    Ok(result)
}

pub(in crate::mesh_recovery) fn fit_ownership_junctions(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    prior: &[Point3],
    regions: &BTreeSet<usize>,
    options: JointFitOptions,
) -> Result<(JointFitCandidate, FrozenSourceFit), JointFitRejection> {
    let context =
        JointFitContext::with_chord_witness(evidence, partition, options, chord_uncertainty)
            .ok_or(JointFitRejection::Neighborhood)?;
    let placed = place_junctions(evidence, partition, prior, regions, &context)?;
    let mut candidate = partition.clone();
    let mut changes = Vec::new();
    let dependent = placed.dependent.iter().copied().collect::<Vec<_>>();
    if !dependent.is_empty() {
        let neighborhood = Neighborhood::new(evidence, partition, &dependent, &context)
            .ok_or(JointFitRejection::Neighborhood)?;
        for &region in &dependent {
            let (carrier, bound) = plane_from_junctions(&neighborhood, region, &placed.positions)?;
            candidate.carriers[region] = Some(carrier);
            if bound > 0.0 {
                changes.push((region, bound));
            }
        }
    }
    for &id in &placed.affected {
        let sample = &context.samples[id];
        let local = sample
            .regions
            .iter()
            .map(|&region| candidate.carriers[region])
            .collect::<Option<Vec<_>>>()
            .ok_or(JointFitRejection::Neighborhood)?;
        if !satisfies(placed.positions[id], &local, KERNEL_EPSILON * 0.01) {
            return Err(JointFitRejection::Sample {
                vertex: id,
                allowance: sample.allowance,
            });
        }
    }
    let supports = placed
        .affected
        .iter()
        .flat_map(|&id| context.vertex_supports[id].iter().copied())
        .collect();
    if !preserves_polygon_orientation(evidence, &placed.positions, &supports) {
        return Err(JointFitRejection::Orientation);
    }
    let frozen = FrozenSourceFit {
        allowances: context.samples.iter().map(|s| s.allowance).collect(),
        changes: changes.clone(),
    };
    Ok((
        JointFitCandidate {
            partition: candidate,
            positions: placed.positions,
            allowances: placed.allowances,
            maximum_displacement: placed.max_displacement,
            refitted_regions: changes.len(),
            max_additional_carrier_residual: changes
                .iter()
                .map(|&(_, bound)| bound)
                .fold(0.0, f64::max),
            carrier_changes: changes,
        },
        frozen,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;

    #[test]
    fn endpoint_error_is_not_a_source_curvature_allowance() {
        let carrier = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let p = Point3::new(1.001, 0.0, 0.0);
        let q = Point3::new(1.001, 0.0, 1.0);
        assert_eq!(
            chord_uncertainty(carrier, p, q).to_bits(),
            0.0_f64.to_bits()
        );
        let q = Point3::new(0.0, 1.001, 1.0);
        assert!(chord_uncertainty(carrier, p, q) > 0.29);
    }

    #[test]
    fn noisy_expansion_places_only_its_shared_junction_and_refits_dependent_plane()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{EvidenceOptions, RegionBoundaryIndex, analyze_mesh_evidence};
        use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
        let perturbed = Point3::new(1.00002, 0.0, 1.0);
        let outside = Point3::new(1.2, -0.5, 0.5);
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 1.0, 1.0),
                    perturbed,
                    outside,
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 2, 3], [0, 3, 4]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let side = evidence.source_planar_support[0].min(evidence.source_planar_support[1]);
        let planar = evidence.source_planar_support[2];
        let mut labels = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        labels[evidence.source_planar_support[0]] = side;
        labels[evidence.source_planar_support[1]] = side;
        let carrier = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let partition = index.repartition_connected(&labels, |id| {
            Some(if id == side {
                carrier
            } else {
                SurfaceCarrier::from_planar_support(&evidence.planar_supports[id])
            })
        })?;
        let policy = JointFitOptions {
            maximum_adjustment: 0.01,
            declared_uncertainty: 0.0,
            angular_resolution: 1.0e-8,
        };
        let (candidate, frozen) = fit_ownership_junctions(
            &evidence,
            &partition,
            &evidence.positions,
            &BTreeSet::from([side]),
            policy,
        )?;
        let id = evidence
            .positions
            .iter()
            .position(|&p| p == perturbed)
            .ok_or("source junction")?;
        let boundary = candidate.positions[id];
        assert!((boundary - Point3::new(1.0, 0.0, 1.0)).norm() <= KERNEL_EPSILON * 0.01);
        let plane = candidate.partition.carriers[planar].ok_or("dependent plane")?;
        assert!(satisfies(
            boundary,
            &[carrier, plane],
            KERNEL_EPSILON * 0.01
        ));
        assert!(satisfies(outside, &[plane], KERNEL_EPSILON * 0.01));
        assert_eq!(candidate.refitted_regions, 1);
        assert!(candidate.maximum_displacement <= candidate.allowances[id]);
        assert!(
            candidate
                .positions
                .iter()
                .zip(&evidence.positions)
                .enumerate()
                .all(|(i, (a, b))| i == id || a == b)
        );
        let next = frozen
            .context(&evidence, &candidate.partition, policy)
            .ok_or("frozen context")?;
        assert_eq!(
            next.samples[id].allowance.to_bits(),
            candidate.allowances[id].to_bits()
        );
        assert!(next.carrier_change[planar] > 0.0);
        let denied = fit_ownership_junctions(
            &evidence,
            &partition,
            &evidence.positions,
            &BTreeSet::from([side]),
            JointFitOptions {
                maximum_adjustment: 1.0e-6,
                ..policy
            },
        );
        assert!(matches!(denied, Err(JointFitRejection::Sample { .. })));
        Ok(())
    }
}
