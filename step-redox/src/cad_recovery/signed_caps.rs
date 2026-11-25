//! Constructive subtraction from translated, oppositely oriented source caps.
//!
//! The cap scanner identifies candidate filled and empty intervals, not
//! recovered solids. Replaying the host and its cutters must independently
//! reproduce the complete oriented source boundary before admission.

use super::{
    CadFragment, CadFragmentSource, CheapestProvenFragment, ConstructiveStep, MaterialFamily,
    SourceRecovery, compose_constructive_sequence_model, solid::recover_source_cap_extrusion_model,
};
use crate::{
    brep,
    cad_ir::{BooleanOp, CadModel, NodeId, ProofStatus, Provenance},
    cad_kernel::CadKernel,
    solid_extrusions::{
        CapInterval, RecoveredSolidExtrusion, RecoveredSourceCapExtrusion, SignedCapPairProposal,
        cap_contains, detect_signed_cap_pairs_for_solid_with_index,
        source_boundary_edge_witness_points,
    },
};
use anyhow::Result;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::{HashMap, HashSet};

fn cap(prism: &RecoveredSolidExtrusion) -> RecoveredSourceCapExtrusion {
    RecoveredSourceCapExtrusion {
        solid_id: prism.solid_id,
        cap_face_id: prism.cap_face_ids[0],
        profile_curves: prism.profile_curves.clone(),
        inner_profile_loops: prism.inner_profile_loops.clone(),
        origin_mm: prism.origin_mm,
        x_axis: prism.x_axis,
        y_axis: prism.y_axis,
        z_axis: prism.z_axis,
        height_mm: prism.height_mm,
        max_residual_mm: prism.max_residual_mm,
    }
}

fn model(prism: &RecoveredSolidExtrusion) -> Result<(CadModel, NodeId)> {
    let (mut model, root) = recover_source_cap_extrusion_model(&cap(prism))?;
    model.set_provenance(
        root,
        Provenance {
            source_entity_ids: std::iter::once(prism.solid_id)
                .chain(prism.cap_face_ids)
                .collect(),
            proof: ProofStatus::HeuristicCandidate,
            max_residual_mm: None,
        },
    )?;
    Ok((model, root))
}

fn cutter_cap_witnesses<K: CadKernel>(
    source: &SourceRecovery<'_, K>,
    proposal: &SignedCapPairProposal,
) -> Option<Vec<[f64; 3]>> {
    let mut points = Vec::new();
    for id in proposal.prism.cap_face_ids {
        for ring in brep::face_loops(id, source.entities, source.index)? {
            for edge in ring.edges {
                points.extend([edge.start_mm, edge.end_mm]);
                points.extend(source_boundary_edge_witness_points(&edge));
            }
        }
    }
    (!points.is_empty()).then_some(points)
}

/// Each endpoint of an empty interval has at most one adjacent matching cap.
/// Select mutual unique shortest pairings; a tied endpoint remains ambiguous
/// rather than picking an entity ID. Disjoint intervals may share the same
/// profile and need not share a global axis.
fn adjacent_pairs<'a>(candidates: &[&'a SignedCapPairProposal]) -> Vec<&'a SignedCapPairProposal> {
    let mut nearest = HashMap::<u64, (usize, f64, bool)>::new();
    let tolerance = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
    for (index, candidate) in candidates.iter().enumerate() {
        let height = candidate.prism.height_mm;
        for face in candidate.prism.cap_face_ids {
            nearest
                .entry(face)
                .and_modify(|old| {
                    if height < old.1 - tolerance {
                        *old = (index, height, false);
                    } else if (height - old.1).abs() <= tolerance {
                        old.2 = true;
                    }
                })
                .or_insert((index, height, false));
        }
    }
    candidates
        .iter()
        .enumerate()
        .filter_map(|(index, &candidate)| {
            candidate
                .prism
                .cap_face_ids
                .iter()
                .all(|face| {
                    nearest
                        .get(face)
                        .is_some_and(|&(closest, _, tied)| closest == index && !tied)
                })
                .then_some(candidate)
        })
        .collect()
}

fn compose(
    host: &SignedCapPairProposal,
    cutters: &[&SignedCapPairProposal],
) -> Result<(CadModel, NodeId)> {
    let mut models = Vec::with_capacity(cutters.len() + 1);
    models.push(model(&host.prism)?);
    for cutter in cutters {
        models.push(model(&cutter.prism)?);
    }
    let steps = models
        .iter()
        .enumerate()
        .map(|(index, (model, root))| ConstructiveStep {
            model,
            root: *root,
            operation: (index != 0).then_some(BooleanOp::Difference),
        })
        .collect::<Vec<_>>();
    compose_constructive_sequence_model(&steps)
}

pub(super) fn recover<K: CadKernel>(
    source: &SourceRecovery<'_, K>,
    ceiling: Option<u64>,
) -> Result<Option<CadFragment>> {
    let proposals = detect_signed_cap_pairs_for_solid_with_index(
        source.solid_id,
        source.entities,
        source.index,
    );
    let cutters = proposals
        .iter()
        .filter(|p| p.interval == CapInterval::Void)
        .filter_map(|p| Some((p, cutter_cap_witnesses(source, p)?)))
        .collect::<Vec<_>>();
    if cutters.is_empty() {
        return Ok(None);
    }
    let mut best = CheapestProvenFragment::default();
    for host in proposals
        .iter()
        .filter(|p| p.interval == CapInterval::Material)
    {
        let host_cap = cap(&host.prism);
        // Actual boundary witnesses of a cutter must be contained by the
        // positive prism if this simple subtraction grammar is sufficient.
        // This is only candidate selection, never a material admission proof.
        // For spline profiles the analytic ray predicate is not available;
        // leave the decision to the complete source proof.
        let spline_host = host
            .prism
            .profile_loops()
            .flatten()
            .any(crate::profile_curves::RecoveredProfileCurve::is_spline);
        let eligible = cutters
            .iter()
            .filter(|(_, points)| spline_host || points.iter().all(|&p| cap_contains(&host_cap, p)))
            .map(|(proposal, _)| *proposal)
            .collect::<Vec<_>>();
        let selected = adjacent_pairs(&eligible);
        if selected.is_empty() {
            continue;
        }
        let (mut model, root) = compose(host, &selected)?;
        if !source
            .proves(
                &model,
                root,
                best.ceiling(ceiling, MaterialFamily::SignedCaps),
            )
            .unwrap_or(false)
        {
            continue;
        }
        let cutter_cap_face_ids = selected
            .iter()
            .map(|c| c.prism.cap_face_ids)
            .collect::<Vec<_>>();
        let mut ids = HashSet::from([source.solid_id]);
        ids.extend(host.prism.cap_face_ids);
        ids.extend(cutter_cap_face_ids.iter().flatten().copied());
        let mut source_entity_ids = ids.into_iter().collect::<Vec<_>>();
        source_entity_ids.sort_unstable();
        model.set_provenance(
            root,
            Provenance {
                source_entity_ids,
                proof: ProofStatus::WithinTolerance,
                max_residual_mm: Some(crate::boundary_proof::ANALYTIC_BOUNDARY_PROOF_TOL_MM),
            },
        )?;
        best.consider(
            MaterialFamily::SignedCaps,
            Ok(Some(CadFragment {
                source: CadFragmentSource::SignedCapAssembly {
                    solid_id: source.solid_id,
                    host_cap_face_ids: host.prism.cap_face_ids,
                    cutter_cap_face_ids,
                },
                model,
                root,
            })),
        )?;
    }
    best.finish()
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;
    use crate::{
        cad_ir::{CadNode, Primitive, RigidTransform},
        cad_kernel::monstertruck::MonstertruckKernel,
        step_graph::build_index,
    };
    use anyhow::Context;
    use monstertruck_io::step::load::ParsedStepSource;

    #[test]
    fn source_side_channel_recovers_as_editable_prism_difference() -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let mut design = CadModel::new();
        let host = design.add_node(CadNode::Primitive(Primitive::Box {
            size_mm: [10.0, 8.0, 6.0],
        }));
        let cutter = design.add_node(CadNode::Primitive(Primitive::Box {
            size_mm: [2.0, 3.0, 2.0],
        }));
        let cutter = design.add_node(CadNode::Transform {
            child: cutter,
            transform: RigidTransform::translation_mm([4.0, -1.0, 2.0]),
        });
        let root = design.add_node(CadNode::Boolean {
            op: BooleanOp::Difference,
            children: vec![host, cutter],
        });
        design.add_root(root)?;
        let evaluated = kernel.evaluate(&design, root)?;
        assert!(kernel.summarize(&evaluated).is_single_solid());
        let step = kernel.to_step(&evaluated)?;
        let exchange = ParsedStepSource::parse(step.as_bytes())?.exchange;
        let entities = &exchange.data[0].entities;
        let index = build_index(entities);
        let ids = brep::manifold_solid_ids(entities);
        let [id] = ids.as_slice() else {
            anyhow::bail!("expected one exported solid");
        };
        let source = SourceRecovery::new(*id, entities, &index, &kernel)
            .context("original source boundary not certified")?;
        let fragment =
            recover(&source, None)?.context("source channel not automatically recovered")?;
        assert!(matches!(
            fragment.source,
            CadFragmentSource::SignedCapAssembly { .. }
        ));
        assert!(source.proves(&fragment.model, fragment.root, None)?);
        Ok(())
    }
}
