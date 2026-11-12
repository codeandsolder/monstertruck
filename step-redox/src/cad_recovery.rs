use crate::cad_ir::{
    BooleanOp, CadModel, CadNode, FusedPeriodicChain, IndexedCount, NodeId, ProofStatus, Provenance,
};
use crate::cad_kernel::CadKernel;
use crate::math3::dot;
use crate::periodic_chains::PeriodicChainPattern;
use crate::step_graph::{build_index, entity_id, simple_record};
use anyhow::{Result, bail};
use monstertruck_io::step::load::ParsedStepSource;
use serde::Serialize;

mod solid;
pub use solid::{
    recover_open_rectangular_sweep_fragment, recover_open_rectangular_sweep_fragments,
    recover_radial_slot_revolution_fragment, recover_radial_slot_revolution_fragments,
    recover_revolved_round_tail_fragment, recover_revolved_round_tail_fragments,
    recover_round_sweep_fragment, recover_round_sweep_fragments, recover_solid_extrusion_fragment,
    recover_solid_extrusion_fragments, recover_solid_revolution_fragment,
    recover_solid_revolution_fragments,
};

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CadFragment {
    pub source: CadFragmentSource,
    pub model: CadModel,
    pub root: NodeId,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum CadFragmentSource {
    SolidExtrusion {
        solid_id: u64,
        cap_face_ids: [u64; 2],
        side_face_ids: Vec<u64>,
    },
    SolidLoft {
        solid_id: u64,
        cap_face_ids: [u64; 2],
        side_face_ids: Vec<u64>,
    },
    SolidRevolution {
        solid_id: u64,
        face_ids: Vec<u64>,
    },
    RadialSlotRevolution {
        solid_id: u64,
        base_face_ids: Vec<u64>,
        slot_face_ids: [u64; 3],
    },
    SolidSweep {
        solid_id: u64,
        face_ids: Vec<u64>,
        closed: bool,
    },
    RevolvedRoundTail {
        solid_id: u64,
        head_face_ids: Vec<u64>,
        tail_face_ids: Vec<u64>,
    },
    PeriodicChain {
        solid_id: u64,
        sites: usize,
    },
    ExtrusionCoreUnion {
        solid_id: u64,
        cap_face_ids: Vec<[u64; 2]>,
    },
    /// Exact ruled end caps added to an extrusion-core union.
    CapExtrusionRuledEnds {
        solid_id: u64,
        cap_face_id: u64,
        tip_face_ids: Vec<u64>,
    },
    ExtrusionCoreRuledEnds {
        solid_id: u64,
        cap_face_ids: Vec<[u64; 2]>,
        tip_face_ids: Vec<u64>,
    },
    /// A source-derived partial solid: a terminal cap and one or more
    /// witnessed planar wall rings; never a recovered complete source body.
    RuledEndExtension {
        solid_id: u64,
        tip_face_id: u64,
        wall_face_ids: Vec<u64>,
    },
    BoundaryFeatureComposite {
        solid_id: u64,
        passes: usize,
        feature_families: usize,
        feature_instances: usize,
    },
    EdgeTreatmentComposite {
        solid_id: u64,
        treatment_face_ids: Vec<u64>,
    },
    SameDomainNormalized {
        solid_id: u64,
    },
}

/// Compose two or more already-recovered CAD fragments into a union model.
///
/// This helper only composes executable CAD models; it deliberately does not create a
/// `CadFragmentSource` claim for the union. Callers must independently prove that the input
/// fragments form the intended source solid before admitting the composite into Complete IR.
///
/// # Errors
/// Returns an error if fewer than two fragments are provided, an input model is invalid, or a
/// fragment root cannot be remapped into the composed model.
pub fn compose_union_model(fragments: &[CadFragment]) -> Result<(CadModel, NodeId)> {
    let models = fragments
        .iter()
        .map(|fragment| (&fragment.model, fragment.root))
        .collect::<Vec<_>>();
    compose_union_models(&models)
}

/// One independently editable CAD feature in a source-ordered constructive
/// sequence. The first feature is positive material (`operation = None`);
/// each following feature combines with the accumulated body in this exact
/// order. Difference is intentionally not commutative.
pub struct ConstructiveStep<'a> {
    pub model: &'a CadModel,
    pub root: NodeId,
    pub operation: Option<BooleanOp>,
}

/// Compose source-derived features into one executable CAD tree, retaining
/// each feature's provenance and preserving the exact Boolean operation order.
///
/// The component models, not imported STEP solids, are the owners of all
/// sketch/profile/loft geometry. Evaluating this model starts at the leaf CAD
/// nodes and never reloads an intermediate B-rep.
///
/// # Errors
/// Refuses empty sequences, an operation on the first feature, a missing
/// subsequent operation, invalid feature roots, cycles or invalid provenance.
pub fn compose_constructive_sequence_model(
    steps: &[ConstructiveStep<'_>],
) -> Result<(CadModel, NodeId)> {
    let (model, stages) = compose_constructive_sequence_model_with_stages(steps)?;
    let root = stages
        .last()
        .copied()
        .ok_or_else(|| anyhow::anyhow!("missing constructive root"))?;
    Ok((model, root))
}

/// Compose an ordered CAD feature sequence and return every cumulative stage
/// root, so an incremental evaluator can inspect each feature result without
/// reconstructing or re-evaluating earlier Boolean subtrees.
///
/// Only the final cumulative stage is registered as the model's root. All
/// intermediate stage IDs remain valid and can be evaluated independently.
///
/// # Errors
/// Refuses empty sequences, invalid feature roots or provenance, and missing
/// or misplaced Boolean operations.
pub fn compose_constructive_sequence_model_with_stages(
    steps: &[ConstructiveStep<'_>],
) -> Result<(CadModel, Vec<NodeId>)> {
    let Some((first, rest)) = steps.split_first() else {
        bail!("constructive sequence must contain at least one feature");
    };
    if first.operation.is_some() {
        bail!("first constructive feature must create positive material");
    }
    let mut model = CadModel::new();
    let mut stages = Vec::with_capacity(steps.len());
    for step in std::iter::once(first).chain(rest) {
        step.model.validate()?;
        step.model.node(step.root)?;
        let offset = model.append_model(step.model)?;
        let owned_root = NodeId(step.root.0.checked_add(offset).ok_or_else(|| {
            anyhow::anyhow!("constructive node id overflow while composing features")
        })?);
        let root = match (stages.last().copied(), step.operation) {
            (None, None) => owned_root,
            (Some(body), Some(operation)) => model.add_node(CadNode::Boolean {
                op: operation,
                children: vec![body, owned_root],
            }),
            (None, Some(_)) => bail!("first constructive feature cannot be a Boolean"),
            (Some(_), None) => bail!("subsequent constructive feature needs a Boolean operation"),
        };
        stages.push(root);
    }
    let root = stages
        .last()
        .copied()
        .ok_or_else(|| anyhow::anyhow!("missing constructive root"))?;
    model.add_root(root)?;
    model.validate()?;
    Ok((model, stages))
}

fn compose_union_models(parts: &[(&CadModel, NodeId)]) -> Result<(CadModel, NodeId)> {
    if parts.len() < 2 {
        bail!("a CAD union composite needs at least two models");
    }

    let mut model = CadModel::new();
    let mut children = Vec::with_capacity(parts.len());
    for &(part, part_root) in parts {
        part.validate()?;
        let offset = model.append_model(part)?;
        let root = NodeId(
            part_root
                .0
                .checked_add(offset)
                .ok_or_else(|| anyhow::anyhow!("CAD root id overflow while composing union"))?,
        );
        model.node(root)?;
        children.push(root);
    }

    let root = add_bounded_union_tree(&mut model, children);
    model.add_root(root)?;
    model.validate()?;
    Ok((model, root))
}

const MAX_FLAT_UNION_CHILDREN: usize = 16;
const AXIS_GROUP_DOT_TOLERANCE: f64 = 1.0e-10;

fn add_bounded_union_tree(model: &mut CadModel, children: Vec<NodeId>) -> NodeId {
    if children.len() == 1 {
        return children[0];
    }
    if children.len() <= MAX_FLAT_UNION_CHILDREN {
        return model.add_node(CadNode::Boolean {
            op: BooleanOp::Union,
            children,
        });
    }

    let mut level = children;
    while level.len() > 1 {
        let mut next = Vec::with_capacity(level.len().div_ceil(2));
        for pair in level.chunks(2) {
            if let [first, second] = pair {
                next.push(model.add_node(CadNode::Boolean {
                    op: BooleanOp::Union,
                    children: vec![*first, *second],
                }));
            } else {
                next.push(pair[0]);
            }
        }
        level = next;
    }
    level[0]
}

fn compose_extrusion_core_union_model(
    fragments: &[CadFragment],
    axes: &[[f64; 3]],
) -> Result<(CadModel, NodeId)> {
    if fragments.len() != axes.len() || fragments.len() < 2 {
        bail!("extrusion-core union composition needs matching fragments and axes");
    }
    if fragments.len() <= MAX_FLAT_UNION_CHILDREN {
        return compose_union_model(fragments);
    }

    let mut model = CadModel::new();
    let mut groups = Vec::<([f64; 3], Vec<NodeId>)>::new();
    for (fragment, &axis) in fragments.iter().zip(axes) {
        fragment.model.validate()?;
        let offset = model.append_model(&fragment.model)?;
        let root = NodeId(
            fragment
                .root
                .0
                .checked_add(offset)
                .ok_or_else(|| anyhow::anyhow!("CAD root id overflow while composing union"))?,
        );
        model.node(root)?;

        if let Some((_, roots)) = groups
            .iter_mut()
            .find(|(group_axis, _)| dot(*group_axis, axis).abs() >= 1.0 - AXIS_GROUP_DOT_TOLERANCE)
        {
            roots.push(root);
        } else {
            groups.push((axis, vec![root]));
        }
    }

    let group_roots = groups
        .into_iter()
        .map(|(_, roots)| add_bounded_union_tree(&mut model, roots))
        .collect::<Vec<_>>();
    let root = add_bounded_union_tree(&mut model, group_roots);
    model.add_root(root)?;
    model.validate()?;
    Ok((model, root))
}
fn normalize_extrusion_cap_pairs(cap_face_ids: &[[u64; 2]]) -> Result<Vec<[u64; 2]>> {
    let mut pairs = cap_face_ids
        .iter()
        .copied()
        .map(|mut pair| {
            pair.sort_unstable();
            pair
        })
        .collect::<Vec<_>>();
    if pairs.iter().any(|pair| pair[0] == 0 || pair[0] == pair[1]) {
        bail!("extrusion-core union cap pairs must contain two distinct nonzero face ids");
    }
    pairs.sort_unstable();
    if pairs.windows(2).any(|pair| pair[0] == pair[1]) {
        bail!("extrusion-core union cap pairs must be unique");
    }
    Ok(pairs)
}

pub(crate) fn rebuilt_kernel_boundary_matches_source<K: CadKernel>(
    kernel: &K,
    evaluated: &K::Evaluated,
    source_solid_id: u64,
    source_entities: &[ruststep::ast::EntityInstance],
) -> Result<bool> {
    let rebuilt_step = kernel.to_step(evaluated)?;
    let rebuilt = ParsedStepSource::parse(rebuilt_step.as_bytes())?.exchange;
    let [data] = rebuilt.data.as_slice() else {
        return Ok(false);
    };
    let rebuilt_entities = &data.entities;
    let rebuilt_solid_ids = rebuilt_entities
        .iter()
        .filter_map(|entity| {
            simple_record(entity)
                .is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
                .then_some(entity_id(entity))
        })
        .collect::<Vec<_>>();
    let [rebuilt_solid_id] = rebuilt_solid_ids.as_slice() else {
        return Ok(false);
    };
    Ok(crate::boundary_proof::compare_analytic_boundaries(
        source_solid_id,
        source_entities,
        *rebuilt_solid_id,
        rebuilt_entities,
    )
    .equivalent)
}

/// Recover and admit a union of sharp extrusion cores only when an independent
/// analytic boundary proof establishes that the evaluated union is the source solid.
///
/// Cap-pair recovery is source-derived and fail-closed. The supplied CAD kernel
/// merely proposes executable geometry: after evaluation and STEP serialization,
/// the result must contain exactly one manifold solid and its complete oriented
/// plane/cylinder boundary signature must equal the source signature. Unsupported
/// carriers, unproven edge loci, disconnected kernel results, or any signature
/// mismatch return `Ok(None)` rather than manufacturing semantic IR.
///
/// `data_section` is the zero-based DATA-section index in `source_step`.
///
/// # Errors
/// Returns an error for malformed API input, a missing DATA section, STEP parsing failures,
/// invalid CAD IR, or a backend evaluation/serialization failure.
pub fn recover_proven_extrusion_core_union_fragment<K: CadKernel>(
    source_step: &[u8],
    data_section: usize,
    solid_id: u64,
    cap_face_ids: &[[u64; 2]],
    kernel: &K,
) -> Result<Option<CadFragment>> {
    if solid_id == 0 || cap_face_ids.is_empty() {
        bail!("extrusion-core construction admission needs a solid and at least one cap pair");
    }

    let canonical_pairs = normalize_extrusion_cap_pairs(cap_face_ids)?;

    let source = ParsedStepSource::parse(source_step)?.exchange;
    let Some(section) = source.data.get(data_section) else {
        bail!("extrusion-core union DATA section {data_section} does not exist");
    };
    let entities = &section.entities;
    let index = build_index(entities);
    let Some(source_face_ids) = crate::brep::solid_face_ids(solid_id, entities, &index) else {
        return Ok(None);
    };

    let mut fragments = Vec::with_capacity(canonical_pairs.len());
    let mut axes = Vec::with_capacity(canonical_pairs.len());
    let mut max_residual_mm = crate::boundary_proof::ANALYTIC_BOUNDARY_PROOF_TOL_MM;
    for &pair in &canonical_pairs {
        let Some(core) =
            crate::solid_extrusions::recover_any_extrusion_core_from_cap_faces_with_index(
                solid_id, pair, entities, &index,
            )
        else {
            return Ok(None);
        };
        max_residual_mm = max_residual_mm.max(core.max_residual_mm);
        axes.push(core.z_axis);
        fragments.push(recover_solid_extrusion_fragment(&core)?);
    }

    // Source-exact planar end treatments can extend beyond every constant
    // section without being separate translated-cap extrusion candidates.
    // Each loft is built from a genuine terminal source cap, its four
    // incident planar walls and their exact shoulder vertices.
    let mut tips = Vec::new();
    if let Some(reduced) = crate::reduced_brep::ReducedSolid::from_step(solid_id, entities, &index)
    {
        for extension in crate::solid_lofts::detect_ruled_end_extensions(&reduced) {
            tips.push(extension.tip_face_id);
            axes.push(extension.axis);
            fragments.push(extension.fragment);
        }
    }
    tips.sort_unstable();
    let admitted_single = (fragments.len() == 1).then(|| fragments[0].source.clone());
    let (mut model, root) = if fragments.len() == 1 {
        let one = fragments
            .pop()
            .ok_or_else(|| anyhow::anyhow!("missing core"))?;
        (one.model, one.root)
    } else {
        compose_extrusion_core_union_model(&fragments, &axes)?
    };
    let evaluated = kernel.evaluate(&model, root)?;
    let summary = kernel.summarize(&evaluated);
    if !summary.geometrically_consistent || summary.shells != 1 {
        return Ok(None);
    }
    if !rebuilt_kernel_boundary_matches_source(kernel, &evaluated, solid_id, entities)? {
        return Ok(None);
    }

    let mut source_entity_ids = Vec::with_capacity(source_face_ids.len() + 1);
    source_entity_ids.push(solid_id);
    source_entity_ids.extend(source_face_ids);
    source_entity_ids.sort_unstable();
    source_entity_ids.dedup();
    model.set_provenance(
        root,
        Provenance {
            source_entity_ids,
            proof: ProofStatus::WithinTolerance,
            max_residual_mm: Some(max_residual_mm),
        },
    )?;
    model.validate()?;

    Ok(Some(CadFragment {
        source: admitted_single.unwrap_or(if tips.is_empty() {
            CadFragmentSource::ExtrusionCoreUnion {
                solid_id,
                cap_face_ids: canonical_pairs,
            }
        } else {
            CadFragmentSource::ExtrusionCoreRuledEnds {
                solid_id,
                cap_face_ids: canonical_pairs,
                tip_face_ids: tips,
            }
        }),
        model,
        root,
    }))
}

/// Propose a partial extrusion from each actual source cap, together with all
/// source-derived terminal lofts, and accept only an independently exact union.
///
/// This handles source solids whose rear cap is intact but whose front cap is
/// partly occluded by ruled protrusions. The source cap, generator station and
/// ruled carrier intersections are all derived from the STEP B-rep.
fn recover_auto_proven_partial_cap_ruled_ends_fragment<K: CadKernel>(
    solid_id: u64,
    entities: &[ruststep::ast::EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
    kernel: &K,
) -> Result<Option<CadFragment>> {
    let Some(reduced) = crate::reduced_brep::ReducedSolid::from_step(solid_id, entities, index)
    else {
        return Ok(None);
    };
    let extensions = crate::solid_lofts::detect_extended_ruled_end_extensions(&reduced);
    if extensions.is_empty() {
        return Ok(None);
    }
    let cap_candidates =
        crate::solid_extrusions::partial_source_cap_extrusion_candidates_with_index(
            solid_id, entities, index,
        );
    for cap in cap_candidates {
        // The source-witnessed shoulder must coincide with the extrusion's
        // end station; the terminal cap must extend *beyond* that station.
        // This is a topological attachment test, not arbitrary subset search.
        const STATION_TOL_MM: f64 = 1.0e-7;
        let aligned = extensions
            .iter()
            .filter(|ext| {
                if dot(ext.axis, cap.z_axis).abs() < 1.0 - AXIS_GROUP_DOT_TOLERANCE {
                    return false;
                }
                let shoulder = dot(
                    crate::math3::sub(ext.shoulder_point_mm, cap.origin_mm),
                    cap.z_axis,
                );
                let tip = dot(
                    crate::math3::sub(ext.tip_point_mm, cap.origin_mm),
                    cap.z_axis,
                );
                (shoulder - cap.height_mm).abs() <= STATION_TOL_MM
                    && tip > cap.height_mm + STATION_TOL_MM
            })
            .collect::<Vec<_>>();
        if aligned.is_empty() {
            continue;
        }
        let tip_face_ids = aligned
            .iter()
            .map(|ext| ext.tip_face_id)
            .collect::<Vec<_>>();
        let cap_face_id = cap.cap_face_id;
        let (cap_model, cap_root) = solid::recover_source_cap_extrusion_model(&cap)?;
        let mut parts = Vec::with_capacity(1 + aligned.len());
        parts.push((&cap_model, cap_root));
        parts.extend(
            aligned
                .iter()
                .map(|ext| (&ext.fragment.model, ext.fragment.root)),
        );
        let (mut model, root) = compose_union_models(&parts)?;
        let Ok(shape) = kernel.evaluate(&model, root) else {
            continue;
        };
        if !kernel.summarize(&shape).is_single_solid() {
            continue;
        }
        if !rebuilt_kernel_boundary_matches_source(kernel, &shape, solid_id, entities)
            .unwrap_or(false)
        {
            continue;
        }
        let mut source_entity_ids = vec![solid_id, cap_face_id];
        source_entity_ids.extend(&tip_face_ids);
        for extension in aligned {
            if let CadFragmentSource::RuledEndExtension { wall_face_ids, .. } =
                &extension.fragment.source
            {
                source_entity_ids.extend(wall_face_ids);
            }
        }
        source_entity_ids.sort_unstable();
        source_entity_ids.dedup();
        model.set_provenance(
            root,
            Provenance {
                source_entity_ids,
                proof: ProofStatus::WithinTolerance,
                max_residual_mm: Some(crate::boundary_proof::ANALYTIC_BOUNDARY_PROOF_TOL_MM),
            },
        )?;
        model.validate()?;
        return Ok(Some(CadFragment {
            source: CadFragmentSource::CapExtrusionRuledEnds {
                solid_id,
                cap_face_id,
                tip_face_ids,
            },
            model,
            root,
        }));
    }
    Ok(None)
}

/// Detect all source-derived extrusion cores, discard those that contradict
/// known source boundary geometry, and admit only the complete surviving union.
///
/// A candidate core cannot have an original source boundary *strictly inside*
/// its volume. Source line-edge witnesses prove such contradictions without
/// assuming that the core's cap pair explains the entire source solid.
/// This geometric filter replaces omission enumeration; no subsets are
/// guessed. Final admission independently compares the complete source and
/// evaluated union's oriented analytic boundaries.
///
/// # Errors
/// Returns an error for malformed STEP, invalid CAD IR, or a backend error
/// evaluating the remaining source-derived construction.
pub fn recover_auto_proven_extrusion_core_union_fragment<K: CadKernel>(
    source_step: &[u8],
    data_section: usize,
    solid_id: u64,
    kernel: &K,
) -> Result<Option<CadFragment>> {
    let source = ParsedStepSource::parse(source_step)?.exchange;
    let Some(section) = source.data.get(data_section) else {
        bail!("automatic extrusion-core union DATA section {data_section} does not exist");
    };
    let entities = &section.entities;
    let index = build_index(entities);
    let proposed =
        crate::solid_extrusions::detect_extrusion_core_candidates_with_index(entities, &index)
            .into_iter()
            .filter(|candidate| candidate.solid_id == solid_id)
            .collect();
    let mut filtered = crate::solid_extrusions::discard_boundary_penetrating_extrusion_cores(
        solid_id, proposed, entities, &index,
    );
    // Prefer an independently source-proven full core (possibly with exact
    // ruled ends) to an overlapping union of redundant partial interpretations.
    // Volume is a discovery-order hint only, never admission evidence.
    filtered.sort_by(|a, b| {
        crate::solid_extrusions::analytic_extrusion_volume(b)
            .unwrap_or(0.0)
            .total_cmp(&crate::solid_extrusions::analytic_extrusion_volume(a).unwrap_or(0.0))
            .then_with(|| a.cap_face_ids.cmp(&b.cap_face_ids))
    });
    let mut seen = std::collections::HashSet::new();
    let mut pairs = Vec::new();
    for candidate in filtered {
        let mut pair = candidate.cap_face_ids;
        pair.sort_unstable();
        if seen.insert(pair) {
            pairs.push(pair);
        }
    }
    for &pair in &pairs {
        // Failure of one candidate's tentative kernel evaluation does not
        // imply the next independent candidate is invalid. Only the final
        // comparison to the original source can admit semantic geometry.
        if let Ok(Some(fragment)) = recover_proven_extrusion_core_union_fragment(
            source_step,
            data_section,
            solid_id,
            &[pair],
            kernel,
        ) {
            return Ok(Some(fragment));
        }
    }
    if pairs.len() >= 2 {
        pairs.sort_unstable();
        if let Ok(Some(fragment)) = recover_proven_extrusion_core_union_fragment(
            source_step,
            data_section,
            solid_id,
            &pairs,
            kernel,
        ) {
            return Ok(Some(fragment));
        }
    }
    recover_auto_proven_partial_cap_ruled_ends_fragment(solid_id, entities, &index, kernel)
}

/// Recover a first-class fused periodic-chain IR fragment from a proven source partition.
///
/// This node is descriptive, not yet executable by Monstertruck/KCL. It references the exact
/// source solid while recording the proven pitch/count grammar; detailed face ownership stays
/// in the detector/debug evidence rather than bloating high-level IR. Backends that do not
/// implement fused periodic chains must reject the node explicitly instead of lowering it as
/// a pattern of independent solids.
///
/// # Errors
/// Returns an error if the chain is not fully proven or its partition is inconsistent.
pub fn recover_periodic_chain_fragment(chain: &PeriodicChainPattern) -> Result<CadFragment> {
    validate_periodic_chain(chain)?;

    let adjacent_site_edges_per_boundary = chain
        .adjacent_site_edge_counts
        .first()
        .copied()
        .unwrap_or(0);

    let mut model = CadModel::new();
    let root = model.add_node(CadNode::PeriodicChain(FusedPeriodicChain {
        source_solid_id: chain.solid_id,
        axis: chain.axis,
        pitch_mm: chain.pitch_mm,
        sites: chain.sites,
        first_site_center_mm: chain.site_centers_mm[0],
        interior_site_face_count: chain.interior_site_face_count,
        site_face_count_overrides: count_overrides(
            &chain.site_face_counts,
            chain.interior_site_face_count,
        ),
        interior_gap_face_count: chain.interior_gap_face_count,
        gap_face_count_overrides: count_overrides(
            &chain.gap_face_counts,
            chain.interior_gap_face_count,
        ),
        stretch_face_count: chain.stretch_face_ids.len(),
        fixed_negative_face_count: chain.fixed_negative_face_ids.len(),
        fixed_positive_face_count: chain.fixed_positive_face_ids.len(),
        adjacent_site_edges_per_boundary,
        repeat_coverage_ratio: chain.repeat_coverage_ratio,
    }));

    model.set_provenance(
        root,
        Provenance {
            source_entity_ids: vec![chain.solid_id],
            // This is a structural partition/grammar proof, not a reconstructed
            // geometry comparison, so there is no geometric residual to invent.
            proof: ProofStatus::StructurallyProven,
            max_residual_mm: None,
        },
    )?;
    model.add_root(root)?;
    model.validate()?;

    Ok(CadFragment {
        source: CadFragmentSource::PeriodicChain {
            solid_id: chain.solid_id,
            sites: chain.sites,
        },
        model,
        root,
    })
}

/// Recover compact CAD-IR fragments for every periodic chain whose detector proof is complete.
///
/// Diagnostic candidates that are not `read_only_proven` remain available through the
/// periodic-chain report, but they are deliberately excluded from canonical CAD IR.
///
/// # Errors
/// Returns the first inconsistency found in a chain already marked as proven.
pub fn recover_periodic_chain_fragments(
    chains: &[PeriodicChainPattern],
) -> Result<Vec<CadFragment>> {
    chains
        .iter()
        .filter(|chain| chain.read_only_proven)
        .map(recover_periodic_chain_fragment)
        .collect()
}

fn count_overrides(counts: &[usize], baseline: usize) -> Vec<IndexedCount> {
    counts
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(index, count)| (count != baseline).then_some(IndexedCount { index, count }))
        .collect()
}

fn validate_periodic_chain(chain: &PeriodicChainPattern) -> Result<()> {
    if !chain.read_only_proven
        || !chain.complete_partition
        || chain.faces_without_geometry != 0
        || chain.nonmanifold_edges != 0
        || chain.nonlocal_cross_site_edges != 0
        || !chain.fixed_middle_face_ids.is_empty()
    {
        bail!("periodic-chain source proof is not complete and editable");
    }
    if chain.solid_id == 0 || chain.sites < 2 {
        bail!("periodic-chain source identity/site count is invalid");
    }
    if !chain.pitch_mm.is_finite() || chain.pitch_mm <= 0.0 {
        bail!("periodic-chain source pitch is invalid");
    }
    if chain.axis.iter().any(|value| !value.is_finite()) {
        bail!("periodic-chain source axis is non-finite");
    }
    let axis_norm = chain.axis[2]
        .mul_add(
            chain.axis[2],
            chain.axis[1].mul_add(chain.axis[1], chain.axis[0] * chain.axis[0]),
        )
        .sqrt();
    if (axis_norm - 1.0).abs() > 1.0e-9 {
        bail!("periodic-chain source axis is not unit length");
    }
    if chain.site_centers_mm.iter().any(|value| !value.is_finite()) {
        bail!("periodic-chain source site centers are non-finite");
    }
    let pitch_tolerance = 1.0e-5_f64.max(chain.pitch_mm.abs() * 1.0e-9);
    if chain
        .site_centers_mm
        .windows(2)
        .any(|pair| ((pair[1] - pair[0]).abs() - chain.pitch_mm).abs() > pitch_tolerance)
    {
        bail!("periodic-chain source site centers disagree with pitch");
    }
    if chain.site_face_ids.len() != chain.sites
        || chain.site_face_counts.len() != chain.sites
        || chain.site_centers_mm.len() != chain.sites
        || chain.gap_face_ids.len() + 1 != chain.sites
        || chain.gap_face_counts.len() + 1 != chain.sites
        || chain.adjacent_site_edge_counts.len() + 1 != chain.sites
    {
        bail!("periodic-chain proof arrays disagree with site count");
    }
    if chain.interior_site_face_count == 0 || chain.site_face_counts.contains(&0) {
        bail!("periodic-chain proof contains an empty site");
    }
    if chain
        .site_face_ids
        .iter()
        .zip(&chain.site_face_counts)
        .any(|(faces, &count)| faces.len() != count)
        || chain
            .gap_face_ids
            .iter()
            .zip(&chain.gap_face_counts)
            .any(|(faces, &count)| faces.len() != count)
    {
        bail!("periodic-chain face counts disagree with face partitions");
    }
    if chain
        .adjacent_site_edge_counts
        .first()
        .is_some_and(|first| {
            chain
                .adjacent_site_edge_counts
                .iter()
                .any(|count| count != first)
        })
    {
        bail!("periodic-chain adjacent-site seam width is inconsistent");
    }

    let mut owned = std::collections::HashSet::new();
    for face in chain
        .site_face_ids
        .iter()
        .flatten()
        .chain(chain.gap_face_ids.iter().flatten())
        .chain(chain.stretch_face_ids.iter())
        .chain(chain.fixed_negative_face_ids.iter())
        .chain(chain.fixed_positive_face_ids.iter())
    {
        if *face == 0 || !owned.insert(*face) {
            bail!("periodic-chain source face partition is invalid");
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "cad_recovery/mixed_fixture_models.rs"]
pub(crate) mod mixed_fixture_models;

#[cfg(test)]
mod tests {
    #![expect(
        clippy::redundant_closure,
        reason = "legacy test fixture code is allowed narrow lint exceptions; production paths remain strict"
    )]
    use super::*;

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn mixed_source_cap_model(face: u64, height: f64) -> Result<(CadModel, NodeId)> {
        let bytes = include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");
        let parsed = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &parsed.data[0].entities;
        let index = build_index(entities);
        let candidate = crate::solid_extrusions::source_cap_extrusion_proposal_at_station(
            171,
            face,
            [1.0, 0.0, 0.0],
            height,
            entities,
            &index,
        )
        .ok_or_else(|| anyhow::anyhow!("face {face} is not a recoverable source sketch"))?;
        solid::recover_source_cap_extrusion_model(&candidate)
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn mixed_four_section_lobe_model() -> Result<(CadModel, NodeId)> {
        use crate::reduced_brep::ReducedSolid;
        let bytes = include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");
        let parsed = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &parsed.data[0].entities;
        let index = build_index(entities);
        let reduced = ReducedSolid::from_step(171, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("source mixed B-rep is missing"))?;
        let extensions = crate::solid_lofts::detect_extended_ruled_end_extensions(&reduced);
        let lobe = extensions
            .iter()
            .find(|entry| {
                entry
                    .fragment
                    .model
                    .nodes
                    .iter()
                    .any(|node| matches!(node,CadNode::Loft{sections} if sections.len()==4))
            })
            .ok_or_else(|| anyhow::anyhow!("source-derived four-section lobe absent"))?;
        Ok((lobe.fragment.model.clone(), lobe.fragment.root))
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn mixed_source_cap_and_four_section_lobe_replay_from_constructive_ir() -> Result<()> {
        use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
        let (base, base_root) = mixed_source_cap_model(620, 1.27)?;
        let (lobe, lobe_root) = mixed_four_section_lobe_model()?;
        let steps = [
            ConstructiveStep {
                model: &base,
                root: base_root,
                operation: None,
            },
            ConstructiveStep {
                model: &lobe,
                root: lobe_root,
                operation: Some(BooleanOp::Union),
            },
        ];
        let (model, root) = compose_constructive_sequence_model(&steps)?;
        assert_eq!(model.roots, [root]);
        assert!(
            model
                .nodes
                .iter()
                .any(|node| matches!(node,CadNode::Loft{sections} if sections.len()==4))
        );
        assert!(
            model
                .nodes
                .iter()
                .any(|node| matches!(node, CadNode::Extrude { .. }))
        );
        let kernel = MonstertruckKernel::standard_mm(1e-5);
        let body = kernel.evaluate(&model, root)?;
        assert!(kernel.summarize(&body).is_single_solid());
        assert_eq!(kernel.summarize(&body).faces, 31);
        if let Ok(folder) = std::env::var("MIXED_CONSTRUCTIVE_IR_OUTPUT_DIR") {
            std::fs::write(
                std::path::Path::new(&folder).join("mixed-cap-lobe-cad-ir.step"),
                kernel.to_step(&body)?,
            )?;
        }
        Ok(())
    }

    /// Full source-measured CAD sketches exercised as one editable six-part
    /// constructive tree. The cap and lobe are extracted from the original
    /// STEP; the remaining analytic fixture models still need automatic
    /// feature detection before this can be admitted as Complete IR.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn mixed_complete_constructive_ir_replays_without_intermediate_step_import() -> Result<()> {
        use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
        let (base, base_root) = mixed_source_cap_model(620, 1.27)?;
        let (lobe, lobe_root) = mixed_four_section_lobe_model()?;
        let (mut lip, lip_inner) = mixed_fixture_models::curved_lip()?;
        let (void, void_root) = mixed_source_cap_model(391, 0.3048)?;
        let (mut rounds, upper_inner, lower_root) = mixed_fixture_models::rounding_caps()?;
        // These local intersections have a different numerical conditioning
        // from the outer source-scale union/difference sequence. The scoped
        // tolerance controls the solver, NOT admissible geometric error.
        let lip_root = lip.add_node(CadNode::SolverTolerance {
            child: lip_inner,
            tolerance_mm: 0.01,
        });
        let upper_root = rounds.add_node(CadNode::SolverTolerance {
            child: upper_inner,
            tolerance_mm: 0.01,
        });
        let steps = [
            ConstructiveStep {
                model: &base,
                root: base_root,
                operation: None,
            },
            ConstructiveStep {
                model: &lobe,
                root: lobe_root,
                operation: Some(BooleanOp::Union),
            },
            ConstructiveStep {
                model: &lip,
                root: lip_root,
                operation: Some(BooleanOp::Union),
            },
            ConstructiveStep {
                model: &void,
                root: void_root,
                operation: Some(BooleanOp::Difference),
            },
            ConstructiveStep {
                model: &rounds,
                root: upper_root,
                operation: Some(BooleanOp::Difference),
            },
            ConstructiveStep {
                model: &rounds,
                root: lower_root,
                operation: Some(BooleanOp::Difference),
            },
        ];
        let (model, stage_roots) = compose_constructive_sequence_model_with_stages(&steps)?;
        let root = *stage_roots
            .last()
            .ok_or_else(|| anyhow::anyhow!("missing CAD stage"))?;
        assert_eq!(model.roots, [root]);
        let kernel = MonstertruckKernel::standard_mm(1e-5);
        assert!(
            model
                .nodes
                .iter()
                .any(|node| matches!(node,CadNode::Loft{sections} if sections.len()==4))
        );
        assert!(model.nodes.iter().any(|node| matches!(
            node,
            CadNode::Boolean {
                op: BooleanOp::Intersection,
                ..
            }
        )));
        let output_folder = std::env::var("MIXED_CONSTRUCTIVE_IR_OUTPUT_DIR").ok();
        if let Some(folder) = &output_folder {
            for (name, component, component_root) in [
                ("lip", &lip, lip_root),
                ("void", &void, void_root),
                ("upper-round", &rounds, upper_root),
                ("lower-round", &rounds, lower_root),
            ] {
                let evaluated = kernel
                    .evaluate(component, component_root)
                    .map_err(|err| anyhow::anyhow!("standalone editable {name} failed: {err:#}"))?;
                anyhow::ensure!(
                    kernel.summarize(&evaluated).is_single_solid(),
                    "standalone editable {name} is not one solid"
                );
                std::fs::write(
                    std::path::Path::new(folder)
                        .join(format!("mixed-editable-component-{name}.step")),
                    kernel.to_step(&evaluated)?,
                )?;
            }
        }
        if let Some(folder) = &output_folder {
            std::fs::write(
                std::path::Path::new(folder).join("mixed-complete-editable-cad-ir.json"),
                serde_json::to_vec_pretty(&model)?,
            )?;
        }
        let stage_names = ["base", "loft", "lip", "void", "upper-round", "lower-round"];
        let expected_faces = [17, 31, 37, 47, 47, 47];
        let mut evaluated = None;
        let mut session = kernel.session(&model)?;
        for (index, &stage_root) in stage_roots.iter().enumerate() {
            let stage = index + 1;
            let start = std::time::Instant::now();
            eprintln!(
                "editable CAD constructive stage {stage} ({}) starting",
                stage_names[index]
            );
            let result = session.evaluate(stage_root).map_err(|err| {
                anyhow::anyhow!(
                    "constructive stage {} ({}) failed: {err:#}",
                    stage,
                    stage_names[stage - 1]
                )
            })?;
            eprintln!(
                "editable CAD constructive stage {stage} ({}) evaluated in {:.3} s",
                stage_names[index],
                start.elapsed().as_secs_f64(),
            );
            let summary = kernel.summarize(&result);
            anyhow::ensure!(
                summary.is_single_solid(),
                "constructive stage {} ({}) is not one solid: {summary:?}",
                stage,
                stage_names[stage - 1]
            );
            anyhow::ensure!(
                summary.faces == expected_faces[stage - 1],
                "constructive stage {} ({}) faces={} expected={}",
                stage,
                stage_names[stage - 1],
                summary.faces,
                expected_faces[stage - 1]
            );
            if let Some(folder) = &output_folder {
                std::fs::write(
                    std::path::Path::new(folder).join(format!(
                        "mixed-editable-stage-{}.step",
                        stage_names[stage - 1]
                    )),
                    kernel.to_step(&result)?,
                )?;
            }
            evaluated = Some(result);
        }
        let evaluated = evaluated.ok_or_else(|| anyhow::anyhow!("no constructive result"))?;
        let summary = kernel.summarize(&evaluated);
        assert!(
            summary.is_single_solid(),
            "complete editable model must produce one body: {summary:?}"
        );
        assert_eq!(summary.faces, 47);
        if let Ok(folder) = std::env::var("MIXED_CONSTRUCTIVE_IR_OUTPUT_DIR") {
            let path = std::path::Path::new(&folder);
            std::fs::write(
                path.join("mixed-complete-editable-cad-ir.step"),
                kernel.to_step(&evaluated)?,
            )?;
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_source_cap_proposals_for_mixed_treatment() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");
        let source = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &source.data[0].entities;
        let index = build_index(entities);
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for &(face, length, filename) in &[(620, 1.27, "base"), (391, 0.3048, "void")] {
            let proposal = crate::solid_extrusions::source_cap_extrusion_proposal_at_station(
                171,
                face,
                [1.0, 0.0, 0.0],
                length,
                entities,
                &index,
            )
            .ok_or_else(|| anyhow::anyhow!("source cap {face} is not a valid profile"))?;
            let (model, root) = solid::recover_source_cap_extrusion_model(&proposal)?;
            let result = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&result).is_single_solid());
            if let Ok(dir) = std::env::var("MIXED_SOURCE_PRISM_OUTPUT_DIR") {
                std::fs::write(
                    std::path::Path::new(&dir).join(format!("native-{filename}-prism.step")),
                    kernel.to_step(&result)?,
                )?;
            }
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn original_bored_cap_and_two_bounded_lobes_recover_automatically() -> Result<()> {
        let source =
            include_bytes!("../validation/fixtures/native_mixed_carrier_chamfer_06d018ae.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let fragment = recover_auto_proven_extrusion_core_union_fragment(source, 0, 19, &kernel)?
            .ok_or_else(|| {
            anyhow::anyhow!("real bored-cap source must admit exact editable CAD")
        })?;
        assert_eq!(
            fragment.source,
            CadFragmentSource::CapExtrusionRuledEnds {
                solid_id: 19,
                cap_face_id: 400,
                tip_face_ids: vec![301, 416],
            },
        );
        assert_eq!(
            fragment
                .model
                .nodes
                .iter()
                .filter(|node| matches!(node, CadNode::Loft { .. }))
                .count(),
            2,
        );
        assert_eq!(
            fragment
                .model
                .nodes
                .iter()
                .filter(|node| matches!(node, CadNode::Extrude { .. }))
                .count(),
            1,
        );
        fragment.model.validate()?;
        let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&evaluated).is_single_solid());
        let exchange = ParsedStepSource::parse(source)?.exchange;
        assert!(rebuilt_kernel_boundary_matches_source(
            &kernel,
            &evaluated,
            19,
            &exchange.data[0].entities
        )?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn source_single_spine_plus_ruled_ends_has_exact_geometry() -> Result<()> {
        let source = include_bytes!("../validation/fixtures/native_beveled_l_hook.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let recovered =
            recover_proven_extrusion_core_union_fragment(source, 0, 84, &[[50, 245]], &kernel)?
                .ok_or_else(|| {
                    anyhow::anyhow!("complete spine and both end treatments must admit")
                })?;
        assert_eq!(
            recovered.source,
            CadFragmentSource::ExtrusionCoreRuledEnds {
                solid_id: 84,
                cap_face_ids: vec![[50, 245]],
                tip_face_ids: vec![330, 355],
            }
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn source_ruled_end_extensions_recover_beveled_elbow() -> Result<()> {
        let source = include_bytes!("../validation/fixtures/native_beveled_l_hook.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let recovered = recover_auto_proven_extrusion_core_union_fragment(source, 0, 84, &kernel)?
            .ok_or_else(|| anyhow::anyhow!("beveled source must admit its ruled end extensions"))?;
        assert_eq!(
            recovered.source,
            CadFragmentSource::ExtrusionCoreRuledEnds {
                solid_id: 84,
                cap_face_ids: vec![[50, 245]],
                tip_face_ids: vec![330, 355],
            }
        );
        assert_eq!(
            recovered
                .model
                .nodes
                .iter()
                .filter(|node| matches!(node, CadNode::Loft { .. }))
                .count(),
            2,
            "both source-derived tapered ends must remain independently editable",
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_open_cap_hook_automatic_material_core_union() -> Result<()> {
        let source = include_bytes!("../validation/fixtures/native_open_cap_hook_union.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let recovered = recover_auto_proven_extrusion_core_union_fragment(source, 0, 30, &kernel)?
            .ok_or_else(|| anyhow::anyhow!("source-derived material cores must reconstruct"))?;
        assert_eq!(
            recovered.source,
            CadFragmentSource::ExtrusionCoreUnion {
                solid_id: 30,
                cap_face_ids: vec![[65, 168], [116, 199]],
            }
        );
        Ok(())
    }

    /// Both source shapes were originally reconstructed using overlapping
    /// extrusion-core unions. Complete analytic boundary proof and independent
    /// OCCT comparisons establish that each is a *single* exact extrusion.
    /// Keep the simpler editable representation rather than pinning history.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_curved_source_automatic_complete_extrusion_regressions() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let samples: &[(&[u8], u64, [u64; 2])] = &[
            (
                include_bytes!("../validation/fixtures/native_curved_extrusion_43d97bd3.step"),
                4,
                [89, 92],
            ),
            (
                include_bytes!("../validation/fixtures/native_curved_extrusion_b917e6ba.step"),
                78,
                [5, 54],
            ),
        ];
        for &(source, solid_id, expected_pair) in samples {
            let recovered =
                recover_auto_proven_extrusion_core_union_fragment(source, 0, solid_id, &kernel)?
                    .ok_or_else(|| {
                        anyhow::anyhow!("curved source solid #{solid_id} was not admitted")
                    })?;
            match &recovered.source {
                CadFragmentSource::SolidExtrusion {
                    solid_id: actual_id,
                    cap_face_ids,
                    ..
                } if *actual_id == solid_id => {
                    // The two source-cap IDs are a set for provenance; the
                    // extrusion's signed axis may choose either as its start.
                    let mut actual_pair = *cap_face_ids;
                    actual_pair.sort_unstable();
                    anyhow::ensure!(
                        actual_pair == expected_pair,
                        "unexpected source caps for solid #{solid_id}: {actual_pair:?}",
                    );
                }
                other => bail!(
                    "expected simpler exact source extrusion for #{solid_id}, found {other:?}",
                ),
            }
            recovered.model.validate()?;
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_open_cap_hook_explicit_core_union() -> Result<()> {
        let source = include_bytes!("../validation/fixtures/native_open_cap_hook_union.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let recovered = recover_proven_extrusion_core_union_fragment(
            source,
            0,
            30,
            &[[65, 168], [116, 199]],
            &kernel,
        )?;
        anyhow::ensure!(
            recovered.is_some(),
            "source-derived two-core union must recover"
        );
        Ok(())
    }

    fn primitive_fragment(solid_id: u64) -> Result<CadFragment> {
        let mut model = CadModel::new();
        let root = model.add_node(CadNode::Primitive(crate::cad_ir::Primitive::Box {
            size_mm: [1.0, 1.0, 1.0],
        }));
        model.add_root(root)?;
        Ok(CadFragment {
            source: CadFragmentSource::SolidExtrusion {
                solid_id,
                cap_face_ids: [solid_id + 1, solid_id + 2],
                side_face_ids: Vec::new(),
            },
            model,
            root,
        })
    }

    #[test]
    fn constructive_sequence_preserves_order_and_feature_provenance() -> Result<()> {
        let fragments = (0..6)
            .map(|index| -> Result<CadModel> {
                let mut model = CadModel::new();
                let body = model.add_node(CadNode::Primitive(crate::cad_ir::Primitive::Box {
                    size_mm: [index as f64 + 1.0, 1.0, 1.0],
                }));
                model.add_root(body)?;
                model.set_provenance(
                    body,
                    crate::cad_ir::Provenance {
                        source_entity_ids: vec![index + 1],
                        proof: crate::cad_ir::ProofStatus::WithinTolerance,
                        max_residual_mm: Some(0.00014),
                    },
                )?;
                Ok(model)
            })
            .collect::<Result<Vec<_>>>()?;
        let operations = [
            None,
            Some(BooleanOp::Union),
            Some(BooleanOp::Union),
            Some(BooleanOp::Difference),
            Some(BooleanOp::Difference),
            Some(BooleanOp::Difference),
        ];
        let steps = fragments
            .iter()
            .zip(operations)
            .map(|(model, operation)| ConstructiveStep {
                model,
                root: NodeId(0),
                operation,
            })
            .collect::<Vec<_>>();
        let (model, root) = compose_constructive_sequence_model(&steps)?;
        assert_eq!(model.roots.as_slice(), &[root]);
        assert_eq!(model.provenance.len(), 6);
        let mut current = root;
        for (index, expected) in [
            (5, BooleanOp::Difference),
            (4, BooleanOp::Difference),
            (3, BooleanOp::Difference),
            (2, BooleanOp::Union),
            (1, BooleanOp::Union),
        ] {
            let CadNode::Boolean { op, children } = model.node(current)? else {
                bail!("missing intermediate Boolean operation");
            };
            assert_eq!(*op, expected);
            assert_eq!(children.len(), 2);
            let leaf = children[1];
            assert_eq!(model.provenance[&leaf].source_entity_ids, vec![index + 1]);
            current = children[0];
        }
        assert_eq!(model.provenance[&current].source_entity_ids, vec![1]);
        assert!(model.node(root).is_ok());
        assert!(compose_constructive_sequence_model(&[]).is_err());
        let mut invalid = steps;
        invalid[0].operation = Some(BooleanOp::Union);
        assert!(compose_constructive_sequence_model(&invalid).is_err());
        invalid[0].operation = None;
        invalid[2].operation = None;
        assert!(compose_constructive_sequence_model(&invalid).is_err());
        Ok(())
    }

    #[test]
    fn union_model_composes_fragment_roots_without_claiming_semantics() -> Result<()> {
        let fragments = [primitive_fragment(10)?, primitive_fragment(20)?];
        let (model, root) = compose_union_model(&fragments)?;
        let CadNode::Boolean { op, children } = model.node(root)? else {
            bail!("expected union root");
        };
        assert_eq!(*op, BooleanOp::Union);
        assert_eq!(children.as_slice(), &[NodeId(0), NodeId(1)]);
        assert_eq!(model.roots.as_slice(), &[root]);
        Ok(())
    }

    #[test]
    fn union_model_keeps_bounded_fragment_roots_flat() -> Result<()> {
        let fragments = (0..16)
            .map(primitive_fragment)
            .collect::<Result<Vec<_>>>()?;
        let (model, root) = compose_union_model(&fragments)?;
        assert_eq!(root, NodeId(16));
        let CadNode::Boolean { op, children } = model.node(root)? else {
            bail!("expected flat union node");
        };
        assert_eq!(*op, BooleanOp::Union);
        assert_eq!(
            children.as_slice(),
            &(0..16).map(NodeId).collect::<Vec<_>>()
        );
        Ok(())
    }

    #[test]
    fn union_model_balances_large_fragment_sets() -> Result<()> {
        let fragments = (0..17)
            .map(|index| primitive_fragment(index))
            .collect::<Result<Vec<_>>>()?;
        let (model, root) = compose_union_model(&fragments)?;
        assert_eq!(root, NodeId(32));
        for node_index in 17..=32 {
            let CadNode::Boolean { op, children } = model.node(NodeId(node_index))? else {
                bail!("expected binary union node");
            };
            assert_eq!(*op, BooleanOp::Union);
            assert_eq!(children.len(), 2);
        }
        Ok(())
    }

    #[test]
    fn extrusion_core_union_groups_large_sets_by_axis() -> Result<()> {
        let fragments = (0..17)
            .map(|index| primitive_fragment(index))
            .collect::<Result<Vec<_>>>()?;
        let axes = vec![[0.0, 0.0, 1.0]; 9]
            .into_iter()
            .chain(vec![[1.0, 0.0, 0.0]; 8])
            .collect::<Vec<_>>();
        let (model, root) = compose_extrusion_core_union_model(&fragments, &axes)?;
        assert_eq!(root, NodeId(19));
        let CadNode::Boolean { op, children } = model.node(root)? else {
            bail!("expected grouped union root");
        };
        assert_eq!(*op, BooleanOp::Union);
        assert_eq!(children.as_slice(), &[NodeId(17), NodeId(18)]);
        for (node, expected_len) in [(NodeId(17), 9), (NodeId(18), 8)] {
            let CadNode::Boolean { op, children } = model.node(node)? else {
                bail!("expected axis-group union node");
            };
            assert_eq!(*op, BooleanOp::Union);
            assert_eq!(children.len(), expected_len);
        }
        Ok(())
    }

    fn periodic_chain_fixture() -> PeriodicChainPattern {
        PeriodicChainPattern {
            solid_id: 42,
            axis: [1.0, 0.0, 0.0],
            pitch_mm: 2.54,
            sites: 4,
            site_centers_mm: vec![0.0, 2.54, 5.08, 7.62],
            lattice_face_family_votes: 8,
            complete_partition: true,
            read_only_proven: true,
            faces_without_geometry: 0,
            nonmanifold_edges: 0,
            cross_site_edges: 6,
            nonlocal_cross_site_edges: 0,
            adjacent_site_edge_counts: vec![2, 2, 2],
            site_face_counts: vec![2, 2, 2, 3],
            interior_site_face_count: 2,
            gap_face_counts: vec![1, 1, 1],
            interior_gap_face_count: 1,
            site_face_ids: vec![vec![10, 11], vec![20, 21], vec![30, 31], vec![40, 41, 42]],
            gap_face_ids: vec![vec![50], vec![51], vec![52]],
            stretch_face_ids: vec![60],
            fixed_negative_face_ids: vec![70],
            fixed_middle_face_ids: Vec::new(),
            fixed_positive_face_ids: vec![71],
            repeat_coverage_ratio: 0.9,
            edge_category_counts: std::collections::BTreeMap::new(),
            site_adjacency_signatures: vec![std::collections::BTreeMap::new(); 4],
            gap_adjacency_signatures: vec![std::collections::BTreeMap::new(); 3],
        }
    }

    #[test]
    fn proven_fused_periodic_chain_becomes_first_class_ir() -> Result<()> {
        let chain = periodic_chain_fixture();
        let fragment = recover_periodic_chain_fragment(&chain)?;
        assert_eq!(
            fragment.source,
            CadFragmentSource::PeriodicChain {
                solid_id: 42,
                sites: 4,
            }
        );
        let CadNode::PeriodicChain(node) = fragment.model.node(fragment.root)? else {
            bail!("expected fused periodic-chain root");
        };
        assert_eq!(node.sites, 4);
        assert_eq!(node.pitch_mm, 2.54);
        assert_eq!(node.adjacent_site_edges_per_boundary, 2);
        assert_eq!(node.interior_site_face_count, 2);
        assert_eq!(node.interior_gap_face_count, 1);
        assert_eq!(
            node.site_face_count_overrides,
            vec![IndexedCount { index: 3, count: 3 }]
        );
        assert_eq!(node.gap_face_count_overrides, Vec::<IndexedCount>::new());
        assert_eq!(node.stretch_face_count, 1);
        assert_eq!(
            fragment.model.provenance[&fragment.root].proof,
            ProofStatus::StructurallyProven
        );
        assert_eq!(
            fragment.model.provenance[&fragment.root].max_residual_mm,
            None
        );
        assert!(fragment.model.complexity_score(fragment.root)? < 200);
        Ok(())
    }

    #[test]
    fn periodic_chain_recovery_fails_closed_on_incomplete_or_ambiguous_proof() -> Result<()> {
        let mut chain = periodic_chain_fixture();
        chain.read_only_proven = false;
        assert!(recover_periodic_chain_fragment(&chain).is_err());

        let mut chain = periodic_chain_fixture();
        chain.site_face_ids[1][0] = chain.site_face_ids[0][0];
        assert!(recover_periodic_chain_fragment(&chain).is_err());

        let mut chain = periodic_chain_fixture();
        chain.adjacent_site_edge_counts[1] = 3;
        assert!(recover_periodic_chain_fragment(&chain).is_err());

        let mut diagnostic = periodic_chain_fixture();
        diagnostic.read_only_proven = false;
        let fragments = recover_periodic_chain_fragments(&[periodic_chain_fixture(), diagnostic])?;
        assert_eq!(fragments.len(), 1);
        Ok(())
    }
}
