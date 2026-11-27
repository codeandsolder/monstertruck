//! Dependency-ordered recovery across native boundary-feature and treatment
//! transactions. STEP projections are read-only analytic recognizer inputs;
//! every edit is applied to the same backend-owned solid.
use super::{
    CadFragment, CadFragmentSource, CadKernel, RecoveredBoundaryFeatureFamily,
    compose_boundary_feature_fragment, prepare_native_boundary_pass, recognize_native_edited_host,
    recognize_single_healed_host, treatment_layers,
};
use crate::{
    SourceEntity,
    brep::face_loops,
    cad_ir::CadModel,
    cad_kernel::{EvaluationSeed, ReplayMode},
    cad_recovery::SourceRecovery,
    kernel_source_load::{KernelTopologyId, KernelTopologyProvenance},
    kernel_topology::{KernelPatchEdgeRef, KernelPatchEdgeUse, SourceTopologyEditor},
    planar_features::apply_boundary_feature_peel_plan,
    reduced_brep::ReducedSolid,
};
use anyhow::{Context, Result, ensure};
use std::{cell::RefCell, collections::BTreeSet};

enum ReplayStep {
    Boundary(Vec<RecoveredBoundaryFeatureFamily>),
    Treatment(Vec<treatment_layers::ReplayLayer>),
}

#[derive(Default)]
pub(super) struct Journal {
    steps: Vec<ReplayStep>,
}

#[derive(Clone, Copy)]
pub(super) enum NextPeel {
    Boundary,
    Treatment,
}

/// The original unpeeled source, not the intermediate host, is the
/// final qualification boundary for a replayed constructive history.
pub(super) struct ReplayTarget<'a, 's, K: CadKernel> {
    pub(super) source: SourceEntity,
    pub(super) proof: Option<&'a SourceRecovery<'s, K>>,
    pub(super) mode: ReplayMode,
}

impl Journal {
    pub(super) fn boundary(&mut self, families: Vec<RecoveredBoundaryFeatureFamily>) {
        if !families.is_empty() {
            self.steps.push(ReplayStep::Boundary(families));
        }
    }

    pub(super) fn treatment(&mut self, layers: Vec<treatment_layers::ReplayLayer>) {
        if !layers.is_empty() {
            self.steps.push(ReplayStep::Treatment(layers));
        }
    }

    fn next_boundary_pass(&self) -> usize {
        self.steps
            .iter()
            .filter_map(|step| match step {
                ReplayStep::Boundary(families) => Some(families.iter()),
                ReplayStep::Treatment(_) => None,
            })
            .flatten()
            .map(|family| family.pass_index)
            .max()
            .map_or(0, |last| last + 1)
    }

    fn replay(&self, mut host: CadFragment, solid_id: u64) -> Result<CadFragment> {
        let passes = self.next_boundary_pass();
        let mut feature_families = 0;
        let mut feature_instances = 0;
        host.model.provenance.clear();
        for step in self.steps.iter().rev() {
            match step {
                ReplayStep::Boundary(families) => {
                    feature_families += families.len();
                    feature_instances += families
                        .iter()
                        .map(|family| family.members.len())
                        .sum::<usize>();
                    host = compose_boundary_feature_fragment(solid_id, &host, families)?;
                }
                ReplayStep::Treatment(layers) => {
                    host = treatment_layers::replay(host, solid_id, layers.clone())?;
                }
            }
        }
        // A mixed composite is still a boundary-feature composite, even
        // when its outermost authored operation is an edge treatment.
        if feature_families != 0 {
            host.source = CadFragmentSource::BoundaryFeatureComposite {
                solid_id,
                passes,
                feature_families,
                feature_instances,
            };
        }
        Ok(host)
    }

    fn evaluate_replay<K: CadKernel>(
        &self,
        kernel: &K,
        host: CadFragment,
        target: &ReplayTarget<'_, '_, K>,
        seed: Option<EvaluationSeed<'_, K>>,
    ) -> Result<CadFragment> {
        let fragment = self.replay(host, target.source.id)?;
        let evaluated =
            kernel.evaluate_replay(&fragment.model, fragment.root, seed, target.mode)?;
        ensure!(
            kernel.summarize(&evaluated).is_single_solid(),
            "native feature replay did not produce one solid"
        );
        if let Some(proof) = target.proof {
            ensure!(
                proof.matches_evaluated_boundary(&evaluated)?,
                "native feature replay does not match the original source boundary"
            );
        }
        Ok(fragment)
    }

    fn recognize_replay<K: CadKernel + SourceTopologyEditor>(
        &self,
        kernel: &K,
        edited: &K::EditedSource,
        target: &ReplayTarget<'_, '_, K>,
    ) -> Result<Option<CadFragment>> {
        let mut native_error = None;
        if let Some(host) = recognize_native_edited_host(kernel, edited)? {
            match self.evaluate_replay(kernel, host, target, None) {
                Ok(fragment) => return Ok(Some(fragment)),
                Err(error) => native_error = Some(error),
            }
        }
        // Serialize the edited host once, only when native recognition is
        // insufficient. Each material candidate must prove this intermediate
        // boundary and pass the full forward replay before ranking.
        let snapshot = kernel.edited_source_to_step(edited)?;
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(snapshot.as_bytes())?.exchange,
        )?;
        let replay_error = RefCell::new(native_error);
        let acceptance = |model: &CadModel, root, seed: Option<EvaluationSeed<'_, K>>| {
            let host = CadFragment {
                source: CadFragmentSource::SameDomainNormalized {
                    solid_id: target.source.id,
                },
                model: model.clone(),
                root,
            };
            match self.evaluate_replay(kernel, host, target, seed) {
                Ok(_) => true,
                Err(error) => {
                    replay_error.borrow_mut().get_or_insert(error);
                    false
                }
            }
        };
        if let Some(host) =
            recognize_single_healed_host(document.as_ref(), kernel, Some(&acceptance), target.mode)?
        {
            // Qualification already evaluated this exact DAG and journal.
            // Recompose only the IR to retain the actual host's source kind;
            // do not repeat kernel evaluation or the original boundary proof.
            return self.replay(host, target.source.id).map(Some);
        }
        if let Some(error) = replay_error.into_inner() {
            return Err(error
                .context("no recognized host could replay the complete native feature journal"));
        }
        Ok(None)
    }

    pub(super) fn finish<K: CadKernel + SourceTopologyEditor>(
        mut self,
        kernel: &K,
        edited: &mut K::EditedSource,
        target: &ReplayTarget<'_, '_, K>,
        provenance: &KernelTopologyProvenance,
        mut next: NextPeel,
        mut stop: Option<String>,
    ) -> Result<Option<CadFragment>> {
        let mut replay_error = None;
        loop {
            // Each successful inverse transaction strictly decreases native
            // face count. Try complete forward replay before peeling further.
            match self.recognize_replay(kernel, edited, target) {
                Ok(Some(fragment)) => return Ok(Some(fragment)),
                Ok(None) => (),
                Err(error) => replay_error = Some(error),
            }
            match next {
                NextPeel::Boundary => {
                    let families = peel_projected_features(
                        kernel,
                        edited,
                        provenance,
                        self.next_boundary_pass(),
                    )?;
                    if families.is_empty() {
                        break;
                    }
                    self.boundary(families);
                    next = NextPeel::Treatment;
                    stop = None;
                }
                NextPeel::Treatment => {
                    let (layers, reason) = treatment_layers::peel_native(
                        kernel,
                        edited,
                        target.source.data_section,
                        provenance,
                    )?;
                    if layers.is_empty() {
                        stop = reason.or(stop);
                        break;
                    }
                    self.treatment(layers);
                    kernel.normalize_edited_source(edited)?;
                    next = NextPeel::Boundary;
                    stop = reason;
                }
            }
        }
        if let Some(reason) = stop {
            anyhow::bail!("native inverse fixed point: {reason}");
        }
        if let Some(error) = replay_error {
            return Err(error
                .context("no recognized host could replay the complete native feature journal"));
        }
        Ok(None)
    }
}

/// Plan feature reductions against a *new* semantic projection after native
/// treatment inversion. The exporter supplies stable-ID correspondence for
/// new edges as well as the original retained faces.
fn peel_projected_features<K: CadKernel + SourceTopologyEditor>(
    kernel: &K,
    edited: &mut K::EditedSource,
    original: &KernelTopologyProvenance,
    first_pass: usize,
) -> Result<Vec<RecoveredBoundaryFeatureFamily>> {
    let projection = kernel.edited_source_step_projection(edited)?;
    let source = projection.source;
    let mut entities = projection
        .exchange
        .data
        .into_iter()
        .next()
        .context("missing native projection DATA section")?
        .entities;
    let index = crate::step_graph::build_index(&entities);
    let old_ids = treatment_layers::original_face_ids(original)?;
    let mut families = Vec::new();
    let mut native = kernel.edited_source_reduced_brep(edited)?;
    let mut pass = first_pass;
    loop {
        let Some(mut prepared) =
            prepare_native_boundary_pass(source, &entities, &index, &projection.topology, pass)?
        else {
            return Ok(families);
        };
        // STEP surface sense may be folded differently by export. Never use
        // it to reorient retained native coedges. Only the identity/shell
        // plan comes from the projection; the replacement wires are native.
        retain_native_wires(
            &mut prepared,
            &native,
            &entities,
            &index,
            &projection.topology,
        )?;
        let mut expected = native.faces().map(|face| face.id).collect::<BTreeSet<_>>();
        for id in &prepared.patch.remove_faces {
            ensure!(
                expected.remove(&id.raw()),
                "projected feature removes a missing native face"
            );
        }
        ensure!(
            expected.len() < native.faces().count(),
            "projected feature does not decrease face count"
        );
        kernel.apply_edited_source_patches(edited, &[prepared.patch])?;
        apply_boundary_feature_peel_plan(&mut entities, &index, &prepared.plan)?;
        native = kernel.edited_source_reduced_brep(edited)?;
        ensure!(
            native.faces().map(|face| face.id).collect::<BTreeSet<_>>() == expected,
            "projected feature transaction changed unexpected native face identities"
        );
        for mut family in prepared.families {
            for member in std::iter::once(&mut family.source_member).chain(&mut family.members) {
                member.face_ids = member
                    .face_ids
                    .iter()
                    .filter_map(|id| {
                        let native = projection.topology.face_for_source_id(0, *id)?;
                        old_ids.get(&native.raw()).copied()
                    })
                    .collect();
                // Snapshot bound IDs are private to its exchange. They must
                // never leak into original document provenance.
                member.interface_bound_ids.clear();
            }
            family.original_source_ids = true;
            families.push(family);
        }
        pass += 1;
    }
}

fn retain_native_wires(
    prepared: &mut super::NativeBoundaryPass,
    native: &ReducedSolid,
    entities: &[step_p21::ast::EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
    topology: &KernelTopologyProvenance,
) -> Result<()> {
    for edit in &prepared.plan.host_bound_removals {
        let id = topology
            .face_for_source_id(0, edit.host_face_id)
            .context("projected feature host has no native identity")?;
        let face = native
            .faces()
            .find(|face| face.id == id.raw())
            .context("projected feature host is missing from current topology")?;
        let snapshot_loops = face_loops(edit.host_face_id, entities, index)
            .context("projected feature host has unreadable bounds")?;
        ensure!(
            snapshot_loops.len() == face.loops.len(),
            "projected native wire count mismatch"
        );
        let mut boundaries = Vec::new();
        for (snapshot, wire) in snapshot_loops.iter().zip(&face.loops) {
            ensure!(
                snapshot.edges.len() == wire.edges.len(),
                "projected native coedge count mismatch"
            );
            for (a, b) in snapshot.edges.iter().zip(&wire.edges) {
                let mapped = topology.edge_for_source_id(0, a.edge_id);
                ensure!(
                    mapped.is_some_and(|id| id.raw() == b.edge_id)
                        && (a.parameter_forward == a.curve_same_sense)
                            == (b.parameter_forward == b.curve_same_sense),
                    "projected native coedge incidence mismatch"
                );
            }
            if !edit.bound_ids.contains(&snapshot.bound_id) {
                boundaries.push(
                    wire.edges
                        .iter()
                        .map(|edge| KernelPatchEdgeUse {
                            edge: KernelPatchEdgeRef::Existing(KernelTopologyId::from_raw(
                                edge.edge_id,
                            )),
                            orientation: edge.parameter_forward == edge.curve_same_sense,
                        })
                        .collect(),
                );
            }
        }
        ensure!(
            !boundaries.is_empty(),
            "projected feature removes all native host wires"
        );
        let replacement = prepared
            .patch
            .face_boundary_replacements
            .iter_mut()
            .find(|replacement| replacement.face == id)
            .context("projected feature has no native host replacement")?;
        replacement.boundaries = boundaries;
    }
    Ok(())
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;
    use crate::{
        cad_ir::{
            BooleanOp, CadModel, CadNode, EdgeSelector, EdgeTreatmentEdge, EdgeTreatmentProfile,
            Profile2d, RigidTransform,
        },
        cad_kernel::monstertruck::MonstertruckKernel,
        cad_recovery::SourceRecovery,
        kernel_source_load::SourceStepLoadInspector,
    };

    fn peel_single_source<K: SourceTopologyEditor>(
        kernel: &K,
        document: &K::Document,
    ) -> Result<(
        SourceEntity,
        KernelTopologyProvenance,
        treatment_layers::PeelOutcome<K::EditedSource>,
    )> {
        let load = kernel.inspect_source_document(document)?;
        let [solid] = load.solids.as_slice() else {
            anyhow::bail!("expected one source body");
        };
        let source = solid.source;
        let provenance = solid
            .topology_provenance
            .clone()
            .context("missing original identity journal")?;
        let section = &document.as_ref().data[source.data_section];
        let index = crate::step_graph::build_index(&section.entities);
        let evidence =
            super::super::replayable_source_treatments(&section.entities, &index, source.id);
        let outcome = treatment_layers::peel(kernel, document, source, &provenance, &evidence)?;
        Ok((source, provenance, outcome))
    }

    #[test]
    fn drafted_round_host_peels_nonconvex_planar_recesses() -> Result<()> {
        let bytes =
            include_bytes!("../../validation/fixtures/native_drafted_round_network_0b5a3544.step");
        let kernel = MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange,
        )?;
        let (
            _,
            provenance,
            treatment_layers::PeelOutcome {
                mut edited,
                layers,
                stop,
            },
        ) = peel_single_source(&kernel, &document)?;
        assert_eq!(layers.len(), 2);
        assert!(stop.is_none(), "{stop:?}");
        kernel.normalize_edited_source(&mut edited)?;
        // Independently prove the materialized tool witnesses against
        // evaluated CAD, without treating the unsupported source solid as
        // a successful recovery.
        let projection = kernel.edited_source_step_projection(&edited)?;
        let mut exchange = projection.exchange;
        let entities = &mut exchange.data[0].entities;
        let index = crate::step_graph::build_index(entities);
        let report = crate::planar_features::diagnose_boundary_features(entities);
        let tools = crate::planar_features::materialize_boundary_feature_tools_from_report(
            entities, &index, &report, 1,
        )?;
        let index = crate::step_graph::build_index(entities);
        let ids = tools.iter().map(|tool| tool.solid_id).collect::<Vec<_>>();
        let tolerances =
            crate::solid_revolutions::source_tolerance_by_representation_item(entities, &index);
        let candidates = super::super::recover_constructive_fragments_with_tolerances(
            entities,
            &index,
            &tolerances,
            &ids,
        )?;
        let mut proven = BTreeSet::new();
        for candidate in &candidates {
            let id = super::super::source_solid_id(&candidate.source);
            if let Some(proof) = SourceRecovery::new(id, entities, &index, &kernel)
                && proof.proves(&candidate.model, candidate.root, None)?
            {
                proven.insert(id);
            }
        }
        let families = peel_projected_features(&kernel, &mut edited, &provenance, 0)?;
        let faces = kernel.inspect_edited_source(&edited)?.faces;
        eprintln!(
            "strict tool proofs={}, peeled families={}, remaining faces={faces}",
            proven.len(),
            families.len()
        );
        assert!(
            proven.len() >= 2,
            "expected strictly proven materialized tools"
        );
        assert!(
            families.len() >= 5,
            "families={}, faces={faces}",
            families.len()
        );
        assert!(
            faces < 136,
            "the previously unrecognized planar tools must be removed"
        );
        // The remaining topology must still be exportable, even though
        // full-source reconstruction is not established.
        kernel.edited_source_step_projection(&edited)?;
        Ok(())
    }

    #[test]
    fn native_treatment_inversion_exposes_mapped_boundary_feature() -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let mut model = CadModel::new();
        let base = model.add_node(CadNode::Extrude {
            profile: Profile2d::polygon(vec![[0., 0.], [10., 0.], [10., 8.], [0., 8.]])?,
            vector_mm: [0., 0., 2.],
        });
        let boss = model.add_node(CadNode::Extrude {
            profile: Profile2d::polygon(vec![[3., 3.], [5., 3.], [5., 5.], [3., 5.]])?,
            vector_mm: [0., 0., 2.],
        });
        let raised = model.add_node(CadNode::Transform {
            child: boss,
            transform: RigidTransform::translation_mm([0., 0., 2.]),
        });
        let joined = model.add_node(CadNode::Boolean {
            op: BooleanOp::Union,
            children: vec![base, raised],
        });
        let root = model.add_node(CadNode::EdgeTreatment {
            child: joined,
            profile: EdgeTreatmentProfile::Round,
            edges: vec![EdgeTreatmentEdge {
                selector: EdgeSelector::Segment {
                    start_mm: [3., 3., 4.],
                    end_mm: [5., 3., 4.],
                },
                radius_mm: 0.25,
            }],
        });
        model.add_root(root)?;
        let evaluated = kernel.evaluate(&model, root)?;
        let text = kernel.to_step(&evaluated)?;
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(text.as_bytes())?.exchange,
        )?;
        let (
            source,
            provenance,
            treatment_layers::PeelOutcome {
                mut edited,
                layers,
                stop,
            },
        ) = peel_single_source(&kernel, &document)?;
        assert!(!layers.is_empty(), "treatment inversion stop: {stop:?}");
        let projection = kernel.edited_source_step_projection(&edited)?;
        assert!(
            projection.topology.edges.iter().any(|entry| {
                !provenance
                    .edges
                    .iter()
                    .any(|original| original.topology_id == entry.topology_id)
            }),
            "inverse-created edges must receive their current native identities in the projection"
        );
        let expected_original = treatment_layers::original_face_ids(&provenance)?;
        let families = peel_projected_features(&kernel, &mut edited, &provenance, 0)?;
        assert_eq!(families.len(), 1, "expected a single restored boss family");
        assert!(families[0].members.iter().all(|member| {
            member.interface_bound_ids.is_empty()
                && !member.face_ids.is_empty()
                && member
                    .face_ids
                    .iter()
                    .all(|id| expected_original.values().any(|original| original == id))
        }));
        assert_eq!(kernel.inspect_edited_source(&edited)?.faces, 6);
        let host = super::super::recognize_native_edited_host(&kernel, &edited)?
            .context("peeled box was not recognized")?;
        let mut journal = Journal::default();
        journal.treatment(layers);
        journal.boundary(families);
        let candidate = journal.replay(host, source.id)?;
        assert!(matches!(
            candidate.model.node(candidate.root)?,
            CadNode::EdgeTreatment { .. }
        ));
        let section = &document.as_ref().data[source.data_section];
        let index = crate::step_graph::build_index(&section.entities);
        let proof = SourceRecovery::new(source.id, &section.entities, &index, &kernel)
            .context("treated boss has no source proof")?;
        let rebuilt = kernel.evaluate(&candidate.model, candidate.root)?;
        assert!(proof.matches_evaluated_boundary(&rebuilt)?);
        Ok(())
    }
}
