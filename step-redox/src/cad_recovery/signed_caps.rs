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

mod convex;

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
    additions: &[AdditiveFeature],
    cutters: &[&SignedCapPairProposal],
    rounds: &[super::channel_round_ends::RoundEndProposal],
) -> Result<(CadModel, NodeId)> {
    let host_model = model(&host.prism)?;
    let mut models = Vec::with_capacity(cutters.len() + rounds.len());
    for cutter in cutters {
        let extended = convex::through_openings(&host.prism, &cutter.prism);
        models.push(model(extended.as_ref().unwrap_or(&cutter.prism))?);
    }
    for round in rounds {
        models.push(round.model()?);
    }
    let mut steps = vec![ConstructiveStep {
        model: &host_model.0,
        root: host_model.1,
        operation: None,
    }];
    steps.extend(additions.iter().map(|addition| ConstructiveStep {
        model: &addition.model,
        root: addition.root,
        operation: Some(BooleanOp::Union),
    }));
    steps.extend(models.iter().map(|(model, root)| ConstructiveStep {
        model,
        root: *root,
        operation: Some(BooleanOp::Difference),
    }));
    compose_constructive_sequence_model(&steps)
}

struct AdditiveFeature {
    model: CadModel,
    root: NodeId,
    face_ids: Vec<u64>,
}

fn additive_features(
    solid: &crate::reduced_brep::ReducedSolid,
    host: &SignedCapPairProposal,
    proposals: &[SignedCapPairProposal],
    extensions: &[crate::solid_lofts::RuledEndExtension],
) -> Result<Vec<AdditiveFeature>> {
    let host_cap = cap(&host.prism);
    let mut features = Vec::new();
    for extension in extensions {
        if !cap_contains(&host_cap, extension.shoulder_point_mm)
            || cap_contains(&host_cap, extension.tip_point_mm)
        {
            continue;
        }
        let CadFragmentSource::RuledEndExtension {
            tip_face_id,
            wall_face_ids,
            ..
        } = &extension.fragment.source
        else {
            continue;
        };
        let face_ids = std::iter::once(*tip_face_id)
            .chain(wall_face_ids.iter().copied())
            .collect();
        features.push(AdditiveFeature {
            model: extension.fragment.model.clone(),
            root: extension.fragment.root,
            face_ids,
        });
    }
    for lip in super::tapered_lip::detect(solid, &host.prism, proposals) {
        let (model, root) = lip.model()?;
        features.push(AdditiveFeature {
            model,
            root,
            face_ids: lip.source_face_ids().to_vec(),
        });
    }
    Ok(features)
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
    let mut best = CheapestProvenFragment::default();
    let reduced = crate::reduced_brep::ReducedSolid::from_step(
        source.solid_id,
        source.entities,
        source.index,
    );
    let extensions = reduced.as_ref().map_or_else(
        Vec::new,
        crate::solid_lofts::detect_extended_ruled_end_extensions,
    );
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
        let mut selected: Vec<&SignedCapPairProposal> = Vec::new();
        for candidate in adjacent_pairs(&eligible) {
            if !selected
                .iter()
                .any(|other| convex::equivalent(&candidate.prism, &other.prism))
            {
                selected.push(candidate);
            }
        }
        let additions = match reduced
            .as_ref()
            .map(|solid| additive_features(solid, host, &proposals, &extensions))
            .transpose()
        {
            Ok(features) => features.unwrap_or_default(),
            Err(error) => {
                best.consider(MaterialFamily::SignedCaps, Err(error))?;
                continue;
            }
        };
        if selected.is_empty() && additions.is_empty() {
            continue;
        }
        let rounds = reduced.as_ref().map_or_else(Vec::new, |solid| {
            selected
                .iter()
                .flat_map(|cutter| super::channel_round_ends::detect(solid, &cutter.prism))
                .collect::<Vec<_>>()
        });
        // Compare the plain interval grammar with all independently witnessed
        // additions and terminations. No arbitrary feature-subset search.
        for (additions, rounds) in [
            (&[][..], &[][..]),
            (&[][..], rounds.as_slice()),
            (additions.as_slice(), &[][..]),
            (additions.as_slice(), rounds.as_slice()),
        ]
        .into_iter()
        .enumerate()
        .filter_map(|(variant, pair)| {
            let (additions, rounds) = pair;
            if (variant & 1 != 0 && rounds.is_empty())
                || (variant & 2 != 0 && additions.is_empty())
                || (selected.is_empty() && additions.is_empty())
            {
                return None;
            }
            Some(pair)
        }) {
            let (mut model, root) = match compose(host, additions, &selected, rounds) {
                Ok(candidate) => candidate,
                Err(error) => {
                    best.consider(MaterialFamily::SignedCaps, Err(error))?;
                    continue;
                }
            };
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
            let round_end_face_ids = rounds
                .iter()
                .flat_map(|round| round.round_face_ids().iter().copied())
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect();
            let additive_face_ids = additions
                .iter()
                .flat_map(|feature| feature.face_ids.iter().copied())
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let mut ids = HashSet::from([source.solid_id]);
            ids.extend(additive_face_ids.iter().copied());
            ids.extend(host.prism.cap_face_ids);
            ids.extend(cutter_cap_face_ids.iter().flatten().copied());
            for round in rounds {
                ids.extend(round.source_entity_ids().iter().copied());
            }
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
                        round_end_face_ids,
                        additive_face_ids,
                    },
                    model,
                    root,
                })),
            )?;
        }
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
