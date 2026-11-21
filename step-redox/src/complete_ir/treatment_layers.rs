//! Fixed-point inverse treatment layers on one backend-owned native solid.
use super::{
    CadFragment, ReplayTreatment, TreatmentReplayGroup, apply_treatment_replay_groups,
    make_treatment_replay_groups, replayable_edge_treatment, select_replayable_treatments,
};
use crate::{
    SourceEntity,
    edge_treatments::{EdgeTreatmentEvidence, detect_native_edge_treatments},
    inverse_edge_treatment::NativeTreatmentTopology,
    kernel_source_load::KernelTopologyProvenance,
    kernel_topology::{
        PlannedTreatmentLayer, PlannedTreatmentTransaction, SourceTopologyEditor,
        plan_treatment_layer, source_reduced_solid,
    },
};
use anyhow::{Context, Result, ensure};
use std::collections::BTreeMap;

pub(super) struct ReplayLayer {
    groups: Vec<TreatmentReplayGroup>,
}

impl ReplayLayer {
    fn new(
        evidence: &[EdgeTreatmentEvidence],
        plan: PlannedTreatmentTransaction,
        source_ids: impl Fn(u64) -> Result<u64>,
    ) -> Result<(crate::kernel_topology::KernelTopologyPatch, Self)> {
        let replay = evidence
            .iter()
            .map(|item| {
                let (profile, edge) = replayable_edge_treatment(item)
                    .context("planned treatment has no forward operation")?;
                Ok((item, profile, edge))
            })
            .collect::<Result<Vec<ReplayTreatment<'_>>>>()?;
        let mut groups = make_treatment_replay_groups(replay, plan.coupled_groups)?;
        ensure!(
            groups
                .iter()
                .map(|group| group.face_ids.len())
                .sum::<usize>()
                == plan.patch.remove_faces.len(),
            "inverse patch and forward replay face provenance counts disagree"
        );
        for group in &mut groups {
            group.face_ids = group
                .face_ids
                .iter()
                .copied()
                .map(&source_ids)
                .collect::<Result<_>>()?;
            group.face_ids.sort_unstable();
        }
        Ok((plan.patch, Self { groups }))
    }
}

fn original_face_ids(provenance: &KernelTopologyProvenance) -> Result<BTreeMap<u64, u64>> {
    let mut ids = BTreeMap::new();
    for face in &provenance.faces {
        let Some(source) = face.shell_member.or(face.source) else {
            continue;
        };
        let id = face.topology_id.raw();
        if let Some(previous) = ids.insert(id, source.id) {
            ensure!(
                previous == source.id,
                "native face #{id} has ambiguous original provenance"
            );
        }
    }
    Ok(ids)
}

/// A layer commits all independent valid neighborhoods. Failed neighborhoods
/// are retained transactionally and reconsidered only against a fresh scan.
fn apply_layer<K: SourceTopologyEditor>(
    kernel: &K,
    edited: &mut K::EditedSource,
    planned: Vec<PlannedTreatmentLayer>,
    source_ids: impl Fn(u64) -> Result<u64>,
) -> Result<Option<ReplayLayer>> {
    if planned.is_empty() {
        return Ok(None);
    }
    let mut prepared = Vec::new();
    let mut refusals = Vec::new();
    for item in planned {
        match ReplayLayer::new(&item.evidence, item.transaction, &source_ids) {
            Ok(value) => prepared.push(value),
            Err(error) => refusals.push(format!("{error:#}")),
        }
    }
    let patches = prepared
        .iter()
        .map(|(patch, _)| patch.clone())
        .collect::<Vec<_>>();
    if !patches.is_empty() && kernel.apply_edited_source_patches(edited, &patches).is_ok() {
        return Ok(Some(ReplayLayer {
            groups: prepared
                .into_iter()
                .flat_map(|(_, layer)| layer.groups)
                .collect(),
        }));
    }
    let mut groups = Vec::new();
    for (patch, layer) in prepared {
        match kernel.apply_edited_source_patches(edited, &[patch]) {
            Ok(()) => groups.extend(layer.groups),
            Err(error) => refusals.push(format!(
                "{:?}: {error:#}",
                layer
                    .groups
                    .iter()
                    .flat_map(|group| &group.face_ids)
                    .collect::<Vec<_>>()
            )),
        }
    }
    if groups.is_empty() {
        anyhow::bail!(
            "no native treatment neighborhood committed: {}",
            refusals.join("; ")
        );
    }
    Ok(Some(ReplayLayer { groups }))
}

pub(super) struct PeelOutcome<T> {
    pub edited: T,
    pub layers: Vec<ReplayLayer>,
    pub stop: Option<String>,
}

/// Apply first-layer STEP-backed evidence, then redetect from the exact
/// edited native geometry after each committed removal. No intermediate
/// serialization, reimport, or stale face/carrier geometry.
pub(super) fn peel<K: SourceTopologyEditor>(
    kernel: &K,
    document: &K::Document,
    source: SourceEntity,
    provenance: &KernelTopologyProvenance,
    initial: &[EdgeTreatmentEvidence],
) -> Result<PeelOutcome<K::EditedSource>> {
    let mut edited = kernel.edit_source_topology(document, source, &[])?;
    let mut layers = Vec::new();
    let mut stop = None;
    if !initial.is_empty() {
        let reduced = source_reduced_solid(document.as_ref(), source)?;
        let first = plan_treatment_layer(&reduced, source.data_section, initial, provenance)
            .and_then(|planned| apply_layer(kernel, &mut edited, planned, Ok));
        match first {
            Ok(Some(layer)) => layers.push(layer),
            Ok(None) => {}
            Err(error) => stop = Some(format!("{error:#}")),
        }
    }
    if !layers.is_empty() {
        let original = original_face_ids(provenance)?;
        loop {
            let faces = kernel.inspect_edited_source(&edited)?.faces;
            let current = kernel.edited_source_reduced_brep(&edited)?;
            let evidence = select_replayable_treatments(detect_native_edge_treatments(&current));
            let lookup = NativeTreatmentTopology::new(&current);
            let next = plan_treatment_layer(&current, source.data_section, &evidence, &lookup)
                .and_then(|planned| {
                    apply_layer(kernel, &mut edited, planned, |id| {
                        original.get(&id).copied().with_context(|| {
                            format!("native treatment face #{id} has no original STEP provenance")
                        })
                    })
                });
            match next {
                Ok(Some(layer)) => {
                    ensure!(
                        kernel.inspect_edited_source(&edited)?.faces < faces,
                        "native inverse treatment layer did not decrease face count"
                    );
                    layers.push(layer);
                }
                Ok(None) => break,
                Err(error) => {
                    stop = Some(format!("{error:#}"));
                    break;
                }
            }
        }
    }
    Ok(PeelOutcome {
        edited,
        layers,
        stop,
    })
}

pub(super) fn replay(
    host: CadFragment,
    solid_id: u64,
    layers: Vec<ReplayLayer>,
) -> Result<CadFragment> {
    apply_treatment_replay_groups(
        host,
        solid_id,
        layers.into_iter().rev().flat_map(|layer| layer.groups),
    )
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;
    use crate::{
        cad_kernel::{CadKernel, monstertruck::MonstertruckKernel},
        kernel_source_load::SourceStepLoadInspector,
    };

    #[test]
    fn native_treatment_layers_reverse_replay_order() -> Result<()> {
        use crate::cad_ir::{
            CadModel, CadNode, EdgeSelector, EdgeTreatmentEdge, EdgeTreatmentProfile, Primitive,
        };
        use crate::cad_recovery::CadFragmentSource;
        let mut model = CadModel::new();
        let root = model.add_node(CadNode::Primitive(Primitive::Box { size_mm: [5.0; 3] }));
        model.add_root(root)?;
        let host = CadFragment {
            source: CadFragmentSource::SolidExtrusion {
                solid_id: 42,
                cap_face_ids: [1, 2],
                side_face_ids: Vec::new(),
            },
            model,
            root,
        };
        let layer = |id, radius_mm| ReplayLayer {
            groups: vec![TreatmentReplayGroup {
                face_ids: vec![id],
                profile: EdgeTreatmentProfile::Round,
                edges: vec![EdgeTreatmentEdge {
                    selector: EdgeSelector::Segment {
                        start_mm: [0.0; 3],
                        end_mm: [0.0, 0.0, 5.0],
                    },
                    radius_mm,
                }],
                max_residual_mm: 0.0,
            }],
        };
        // Outer feature removed first; its carrier exists only after the
        // second inverse layer has been replayed forward.
        let recovered = replay(host, 42, vec![layer(100, 0.1), layer(50, 0.5)])?;
        assert_eq!(
            recovered.source,
            CadFragmentSource::EdgeTreatmentComposite {
                solid_id: 42,
                treatment_face_ids: vec![50, 100]
            }
        );
        let CadNode::EdgeTreatment { child, edges, .. } = recovered.model.node(recovered.root)?
        else {
            anyhow::bail!("missing outer replay");
        };
        assert!((edges[0].radius_mm - 0.1).abs() < f64::EPSILON);
        assert_eq!(
            recovered.model.provenance[&recovered.root].source_entity_ids,
            vec![42, 100]
        );
        assert_eq!(
            recovered.model.provenance[child].source_entity_ids,
            vec![42, 50]
        );
        let CadNode::EdgeTreatment { child: inner, .. } = recovered.model.node(*child)? else {
            anyhow::bail!("missing inner replay");
        };
        assert_eq!(*inner, root);
        Ok(())
    }

    fn assert_te_coupled_rounds(layers: &[ReplayLayer], stop: Option<&str>) -> Result<()> {
        // Native face-sense composition must permit the four radius-2
        // vertical rounds of the reversed planar carrier shell to peel.
        let removed = layers
            .iter()
            .flat_map(|layer| &layer.groups)
            .flat_map(|group| group.face_ids.iter().copied())
            .collect::<std::collections::BTreeSet<_>>();
        for id in [8019, 8047, 8075, 8102] {
            assert!(
                removed.contains(&id),
                "reversed native carrier blocked treatment #{id}"
            );
        }
        assert!(layers.first().is_some_and(|layer| {
            layer
                .groups
                .iter()
                .map(|group| group.face_ids.len())
                .sum::<usize>()
                >= 64
        }));
        // The last four spherical corner patches belong to the same
        // eight-edge network as the surrounding cylindrical fillets.
        for id in [
            8430, 8459, 8487, 8513, 8527, 8554, 8581, 9149, 8472, 9160, 8567, 9171,
        ] {
            assert!(
                removed.contains(&id),
                "spherical round network retained #{id}: {:?}",
                stop
            );
        }
        let spherical_faces = std::collections::BTreeSet::from([
            8430, 8459, 8487, 8513, 8527, 8554, 8581, 9149, 8472, 9160, 8567, 9171,
        ]);
        let network = layers
            .iter()
            .flat_map(|layer| &layer.groups)
            .find(|group| {
                group
                    .face_ids
                    .iter()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>()
                    == spherical_faces
            })
            .context("spherical junctions must share one replay operation")?;
        assert_eq!(network.edges.len(), 8);
        assert_eq!(network.profile, crate::cad_ir::EdgeTreatmentProfile::Round);
        assert!(
            network
                .edges
                .iter()
                .all(|edge| matches!(edge.selector, crate::cad_ir::EdgeSelector::Segment { .. }))
        );
        Ok(())
    }

    #[test]
    fn native_treatment_layers_real_te_outer_layer() -> Result<()> {
        let bytes = crate::inverse_edge_treatment::te_treatment_fixture()?;
        let kernel = MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(&bytes)?.exchange,
        )?;
        let load = kernel.inspect_source_document(&document)?;
        let solid = load
            .solids
            .iter()
            .find(|s| s.source.id == 9354)
            .context("missing TE")?;
        let provenance = solid
            .topology_provenance
            .as_ref()
            .context("missing provenance")?;
        let source = solid.source;
        let initial = super::super::replayable_source_treatments(document.as_ref(), source)?;
        let result = peel(&kernel, &document, source, provenance, &initial)?;
        let summary = kernel.inspect_edited_source(&result.edited)?;
        let original = original_face_ids(provenance)?;
        let reduced = kernel.edited_source_reduced_brep(&result.edited)?;
        assert!(reduced.faces().all(|face| original.contains_key(&face.id)));
        let counts = result
            .layers
            .iter()
            .map(|layer| {
                layer
                    .groups
                    .iter()
                    .map(|group| group.face_ids.len())
                    .sum::<usize>()
            })
            .collect::<Vec<_>>();
        eprintln!(
            "TE layers {counts:?}, remaining {}, stop {:?}",
            summary.faces, result.stop
        );
        assert_eq!(summary.boundary_shells, 1);
        assert_eq!(summary.faces, 309 - counts.iter().sum::<usize>());
        assert_te_coupled_rounds(&result.layers, result.stop.as_deref())?;
        assert!(summary.faces <= 225);
        assert!(
            counts.len() >= 2,
            "TE must peel at least two native layers: {:?}",
            result.stop
        );
        // The same real fixture must now admit a constructive peeled host,
        // not merely operate on its native face topology.
        let mut edited = result.edited;
        kernel.normalize_edited_source(&mut edited)?;
        let snapshot = kernel.edited_source_to_step(&edited)?;
        let exchange =
            monstertruck_io::step::load::ParsedStepSource::parse(snapshot.as_bytes())?.exchange;
        let host = super::super::recognize_single_healed_host(&exchange, &kernel)?
            .context("normalized TE host did not pass analytic boundary proof")?;
        assert!(matches!(
            host.source,
            crate::cad_recovery::CadFragmentSource::SectionalExtrusionUnion { .. }
        ));
        let rebuilt = kernel.evaluate(&host.model, host.root)?;
        assert!(kernel.summarize(&rebuilt).is_single_solid());
        let recovered = replay(host, source.id, result.layers)?;
        let reconstructed = kernel.evaluate(&recovered.model, recovered.root)?;
        let reconstructed_summary = kernel.summarize(&reconstructed);
        assert!(
            reconstructed_summary.is_single_solid(),
            "{reconstructed_summary:?}"
        );
        assert_eq!(reconstructed_summary.faces, 309);
        // Every reconstructed support, bounded edge locus, and oriented
        // boundary chain must agree with the original, including the
        // cylindrical, toroidal, and spherical treatment surfaces.
        assert!(
            crate::cad_recovery::rebuilt_kernel_boundary_matches_source(
                &kernel,
                &reconstructed,
                source.id,
                &document.as_ref().data[source.data_section].entities,
            )?,
            "replayed TE differs from the original analytic boundary"
        );
        Ok(())
    }

    #[test]
    fn real_te_is_automatically_admitted_with_source_boundary_proof() -> Result<()> {
        let bytes = crate::inverse_edge_treatment::te_treatment_fixture()?;
        let kernel = MonstertruckKernel::project_standard();
        let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(&bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.unrecovered_solid_count, 0,
            "{:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.constructively_recovered_solid_count, 1);
        let [selected] = report.edge_treatment_composites.as_slice() else {
            anyhow::bail!(
                "expected recovered treatment composite, failures: {:?}",
                report.kernel_recovery_failures
            );
        };
        assert_eq!(selected.source_solid.id, 9354);
        assert_eq!(
            report.source_boundary_verified_solids,
            [selected.source_solid]
        );
        let fragment = &selected.fragment;
        fragment.model.validate()?;
        let reconstructed = kernel.evaluate(&fragment.model, fragment.root)?;
        let summary = kernel.summarize(&reconstructed);
        assert!(summary.is_single_solid(), "{summary:?}");
        let original = monstertruck_io::step::load::ParsedStepSource::parse(&bytes)?.exchange;
        assert!(crate::cad_recovery::rebuilt_kernel_boundary_matches_source(
            &kernel,
            &reconstructed,
            selected.source_solid.id,
            &original.data[selected.source_solid.data_section].entities,
        )?);
        // An editable result must respond to a real parameter change, not
        // just reproduce a fixed source snapshot. The terminal round has no
        // downstream dependent selectors and can be reduced independently.
        let mut edited = fragment.model.clone();
        let Some(crate::cad_ir::CadNode::EdgeTreatment {
            profile: crate::cad_ir::EdgeTreatmentProfile::Round,
            edges,
            ..
        }) = edited.nodes.get_mut(fragment.root.0)
        else {
            anyhow::bail!("expected an editable terminal round");
        };
        ensure!(!edges.is_empty(), "terminal round has no selectors");
        for edge in edges {
            edge.radius_mm *= 0.75;
        }
        edited.provenance.remove(&fragment.root);
        edited.validate()?;
        let shape = kernel.evaluate(&edited, fragment.root)?;
        let summary = kernel.summarize(&shape);
        assert!(summary.is_single_solid(), "edited round: {summary:?}");
        let step = kernel.to_step(&shape)?;
        let exchange =
            monstertruck_io::step::load::ParsedStepSource::parse(step.as_bytes())?.exchange;
        let [section] = exchange.data.as_slice() else {
            anyhow::bail!("edited round export must contain one DATA section");
        };
        let solids = section
            .entities
            .iter()
            .filter_map(|entity| {
                crate::step_graph::simple_record(entity)
                    .is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
                    .then_some(crate::step_graph::entity_id(entity))
            })
            .collect::<Vec<_>>();
        let [id] = solids.as_slice() else {
            anyhow::bail!("edited round export must contain one solid");
        };
        let comparison = crate::boundary_proof::compare_analytic_boundaries(
            selected.source_solid.id,
            &original.data[selected.source_solid.data_section].entities,
            *id,
            &section.entities,
        );
        assert!(
            comparison.first.supported && comparison.second.supported,
            "edited boundary must remain certifiable: {:?}",
            comparison.second.reason
        );
        assert!(
            !comparison.equivalent,
            "radius edit must change the certified boundary"
        );
        Ok(())
    }
}
