use crate::SourceEntity;
use crate::cad_ir::{
    BooleanOp, CadModel, CadNode, EdgeSelector, EdgeTreatmentEdge, EdgeTreatmentProfile, NodeId,
    ProofStatus, Provenance, RigidTransform,
};
use crate::cad_kernel::CadKernel;
use crate::cad_recovery::{
    CadFragment, CadFragmentSource, CandidateAcceptance, SourceRecovery,
    recover_open_rectangular_sweep_fragments, recover_periodic_chain_fragments,
    recover_proven_material_fragment, recover_radial_slot_revolution_fragments,
    recover_revolved_round_tail_fragments, recover_round_sweep_fragments,
    recover_solid_extrusion_fragment, recover_solid_extrusion_fragments,
    recover_solid_revolution_fragments,
};
#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
use crate::cad_recovery::{
    rebuilt_kernel_boundary_matches_source, recover_auto_proven_extrusion_core_union_fragment,
};
use crate::kernel_source_load::{
    KernelSourceLoadEvidence, KernelTopologyProvenance, SourceStepLoadInspector,
};
#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
use crate::kernel_topology::plan_source_coplanar_normalization_patch;
use crate::kernel_topology::{
    SourceTopologyEditor, plan_source_boundary_feature_peel_from_entities,
};
use crate::patterns::InstancePattern;
use crate::planar_features::{
    BoundaryFeatureInstanceEvidence, BoundaryFeaturePeelPlan, FeaturePolarity,
};
use crate::solid_revolutions::SolidSurfaceSignature;
use anyhow::{Context, Result, bail};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use step_p21::ast::{EntityInstance, Exchange};

mod native_peels;
mod treatment_layers;

pub const COMPLETE_IR_SCHEMA: &str = "step-redox-complete-ir-v13";
pub const KERNEL_COMPLETE_IR_SCHEMA: &str = "step-redox-complete-ir-v21";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReducedBrepMutationMode {
    ParserOwned,
    KernelOwned,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SurfaceSignatureSummary {
    pub source: SourceEntity,
    pub face_count: usize,
    pub support_counts: BTreeMap<String, usize>,
    pub unique_edge_count: usize,
    pub edge_use_count: usize,
    pub closed_two_manifold: bool,
}

impl SurfaceSignatureSummary {
    fn from_source(source: SourceEntity, value: &SolidSurfaceSignature) -> Self {
        Self {
            source,
            face_count: value.face_count,
            support_counts: value.support_counts.clone(),
            unique_edge_count: value.unique_edge_count,
            edge_use_count: value.edge_use_count,
            closed_two_manifold: value.closed_two_manifold,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct SourceSolidSummary {
    pub source: SourceEntity,
    pub faces: usize,
    pub edges: usize,
    pub closure_entities: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct BodyScanSummary {
    pub solids: Vec<SourceSolidSummary>,
}

/// One recovered CAD fragment qualified by its document-global source solid.
///
/// Raw STEP entity IDs inside `fragment.source` and node provenance are local to
/// `source_solid.data_section`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SourceCadFragment {
    pub source_solid: SourceEntity,
    pub fragment: CadFragment,
}

/// Instance-pattern diagnostics qualified by the DATA section owning their raw STEP IDs.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SourceInstancePattern {
    pub data_section: usize,
    pub pattern: InstancePattern,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KernelRecoveryFailure {
    pub source: SourceEntity,
    pub error: String,
}

/// A candidate verified against a source proof view in which only uniquely
/// contradictory, fully proven spline-strip face-sense flags were reversed.
///
/// This is deliberately distinct from verification of the authored source.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceOrientationReconciliation {
    pub source: SourceEntity,
    pub repairs: Vec<crate::boundary_proof::StripFaceSenseRepair>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CompleteIr {
    pub schema: &'static str,
    pub surface_signature_summaries: Vec<SurfaceSignatureSummary>,
    pub instance_pattern_diagnostics: Vec<SourceInstancePattern>,
    pub solid_extrusions: Vec<SourceCadFragment>,
    pub solid_lofts: Vec<SourceCadFragment>,
    pub solid_revolutions: Vec<SourceCadFragment>,
    pub radial_slot_revolutions: Vec<SourceCadFragment>,
    pub solid_sweeps: Vec<SourceCadFragment>,
    pub periodic_chains: Vec<SourceCadFragment>,
    pub boundary_feature_composites: Vec<SourceCadFragment>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub edge_treatment_composites: Vec<SourceCadFragment>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub same_domain_normalized_solids: Vec<SourceCadFragment>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub extrusion_core_unions: Vec<SourceCadFragment>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub kernel_recovery_failures: Vec<KernelRecoveryFailure>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kernel_source_load: Option<KernelSourceLoadEvidence>,
    pub unrecovered_solids: Vec<SourceSolidSummary>,
    pub body_scan: BodyScanSummary,
    pub body_count: usize,
    /// Selected models independently compared with their original oriented
    /// analytic source boundaries, not merely structurally reconstructed.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub source_boundary_verified_solids: Vec<SourceEntity>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub source_orientation_reconciled_solids: Vec<SourceOrientationReconciliation>,
    pub constructively_recovered_solid_count: usize,
    pub semantically_recovered_solid_count: usize,
    pub unrecovered_solid_count: usize,
}

impl CompleteIr {
    /// All selected whole-solid fragments, regardless of semantic family.
    /// The source ID belongs to its original DATA section.
    pub fn selected_fragments(&self) -> impl Iterator<Item = &SourceCadFragment> {
        [
            self.solid_extrusions.as_slice(),
            self.solid_lofts.as_slice(),
            self.solid_revolutions.as_slice(),
            self.radial_slot_revolutions.as_slice(),
            self.solid_sweeps.as_slice(),
            self.periodic_chains.as_slice(),
            self.boundary_feature_composites.as_slice(),
            self.edge_treatment_composites.as_slice(),
            self.same_domain_normalized_solids.as_slice(),
            self.extrusion_core_unions.as_slice(),
        ]
        .into_iter()
        .flatten()
    }
}

#[derive(Debug, Clone)]
struct RecoveredBoundaryFeatureFamily {
    pass_index: usize,
    polarity: FeaturePolarity,
    tool: CadFragment,
    source_member: BoundaryFeatureInstanceEvidence,
    members: Vec<BoundaryFeatureInstanceEvidence>,
    /// False when feature IDs came from a kernel-generated intermediate STEP.
    original_source_ids: bool,
}

#[derive(Debug)]
struct RecoverySection<'a> {
    data_section: usize,
    entities: &'a [EntityInstance],
    // A session lends its cached index; standalone kernel recovery owns one.
    index: Cow<'a, HashMap<u64, usize>>,
}

impl RecoverySection<'_> {
    fn source_entity(&self, id: u64) -> Option<SourceEntity> {
        crate::step_document::StepSection::new(self.data_section, self.entities, &self.index)
            .entity(id)
            .map(crate::step_document::StepEntity::source)
    }

    fn local_target_ids(&self, targets: &BTreeSet<SourceEntity>) -> BTreeSet<u64> {
        targets
            .iter()
            .filter(|source| source.data_section == self.data_section)
            .map(|source| source.id)
            .collect()
    }

    fn qualify_fragment(&self, fragment: CadFragment) -> Result<SourceCadFragment> {
        let solid_id = source_solid_id(&fragment.source);
        let source = self.source_entity(solid_id).ok_or_else(|| {
            anyhow::anyhow!(
                "recovered CAD fragment refers to missing solid #{solid_id} in DATA section {}",
                self.data_section
            )
        })?;
        Ok(SourceCadFragment {
            source_solid: source,
            fragment,
        })
    }

    fn qualify_fragments<I>(&self, fragments: I) -> Result<Vec<SourceCadFragment>>
    where I: IntoIterator<Item = CadFragment> {
        fragments
            .into_iter()
            .map(|fragment| self.qualify_fragment(fragment))
            .collect()
    }
}

#[derive(Debug)]
struct LocalRecoveryTargets {
    solid_ids: BTreeSet<u64>,
    shell_to_solid: BTreeMap<u64, u64>,
}

impl LocalRecoveryTargets {
    fn from_solids(
        solid_ids: &BTreeSet<u64>,
        entities: &[EntityInstance],
        index: &HashMap<u64, usize>,
    ) -> Result<Self> {
        let mut shell_to_solid = BTreeMap::new();
        for &solid_id in solid_ids {
            let Some(shell_id) = crate::brep::manifold_shell(solid_id, entities, index) else {
                continue;
            };
            if let Some(previous) = shell_to_solid.insert(shell_id, solid_id)
                && previous != solid_id
            {
                bail!(
                    "local recovery shell #{shell_id} belongs to both solid #{previous} and #{solid_id}"
                );
            }
        }
        Ok(Self {
            solid_ids: solid_ids.clone(),
            shell_to_solid,
        })
    }

    const fn is_empty(&self) -> bool { self.solid_ids.is_empty() }

    fn shell_ids(&self) -> HashSet<u64> { self.shell_to_solid.keys().copied().collect() }

    fn solid_for_shell(&self, shell_id: u64) -> Option<u64> {
        self.shell_to_solid.get(&shell_id).copied()
    }

    fn retain_unrecovered(&mut self, recovered: &BTreeMap<u64, CadFragment>) {
        self.solid_ids
            .retain(|solid_id| !recovered.contains_key(solid_id));
        self.shell_to_solid
            .retain(|_, solid_id| self.solid_ids.contains(solid_id));
    }
}

#[derive(Debug, Clone)]
enum LocalReductionStep {
    BoundaryFeature(RecoveredBoundaryFeatureFamily),
}

/// Parser-only semantic peel journal. The original Part 21 graph remains immutable;
/// a sequence of proven deletion plans projects a temporary entity view for the
/// existing recognizers, and is never mistaken for a separately owned B-rep.
/// The kernel path does not construct this journal.
#[derive(Debug, Default)]
struct ParserPeelJournal {
    plans: Vec<BoundaryFeaturePeelPlan>,
    replay: BTreeMap<u64, Vec<LocalReductionStep>>,
}

impl ParserPeelJournal {
    fn view<'a>(
        &self,
        source: &'a [EntityInstance],
        index: &HashMap<u64, usize>,
    ) -> Result<Cow<'a, [EntityInstance]>> {
        if self.plans.is_empty() {
            return Ok(Cow::Borrowed(source));
        }
        Ok(Cow::Owned(
            crate::planar_features::materialize_boundary_feature_peel_plans(
                source,
                index,
                &self.plans,
            )?,
        ))
    }

    fn commit(
        &mut self,
        steps: Vec<(u64, LocalReductionStep)>,
        plan: BoundaryFeaturePeelPlan,
        current: &[EntityInstance],
        index: &HashMap<u64, usize>,
    ) -> Result<bool> {
        if steps.is_empty() {
            return Ok(false);
        }
        if steps.len() != plan.stats.families || plan.stats.families == 0 {
            bail!(
                "parser peel journal has {} replay steps but {} proven families",
                steps.len(),
                plan.stats.families
            );
        }
        // Prove that the next projection is actually materializable before
        // committing either topology removals or semantic replay. A failed
        // transaction leaves the journal (and original input) unchanged.
        let _ = crate::planar_features::materialize_boundary_feature_peel_plans(
            current,
            index,
            std::slice::from_ref(&plan),
        )?;
        self.plans.push(plan);
        for (solid_id, step) in steps {
            self.replay.entry(solid_id).or_default().push(step);
        }
        Ok(true)
    }

    fn boundary_features(&self, solid_id: u64) -> Vec<RecoveredBoundaryFeatureFamily> {
        self.replay
            .get(&solid_id)
            .into_iter()
            .flatten()
            .map(|step| match step {
                LocalReductionStep::BoundaryFeature(feature) => feature.clone(),
            })
            .collect()
    }

    fn next_boundary_pass_index(&self) -> usize {
        self.replay
            .values()
            .flatten()
            .map(|step| match step {
                LocalReductionStep::BoundaryFeature(feature) => feature.pass_index,
            })
            .max()
            .map_or(0, |last| last + 1)
    }
}

#[derive(Debug)]
struct PreparedBoundaryReduction {
    steps: Vec<(u64, LocalReductionStep)>,
    family_indices: Vec<usize>,
}

/// Stronger evidence takes precedence over structural ranking. Original
/// proof is never reported for a metadata-reconciled source.
#[derive(Debug, Clone, PartialEq, Eq)]
enum BoundaryVerification {
    Structural,
    Reconciled(Vec<crate::boundary_proof::StripFaceSenseRepair>),
    Original,
}

impl BoundaryVerification {
    const fn rank(&self) -> u8 {
        match self {
            Self::Structural => 0,
            Self::Reconciled(_) => 1,
            Self::Original => 2,
        }
    }

    const fn is_certified(&self) -> bool { !matches!(self, Self::Structural) }
}

#[derive(Debug)]
struct SelectedWholeSolid {
    fragment: CadFragment,
    complexity: u64,
    /// This candidate passed the independent, oriented source-boundary check,
    /// not only a source-topology pattern recognizer.
    boundary_verification: BoundaryVerification,
}

#[derive(Debug)]
struct WholeSolidRecoveryState {
    source_solids: BTreeSet<SourceEntity>,
    selected: BTreeMap<SourceEntity, SelectedWholeSolid>,
    revision: u64,
}

impl WholeSolidRecoveryState {
    fn new(body_scan: &[SourceSolidSummary]) -> Self {
        Self {
            source_solids: body_scan.iter().map(|solid| solid.source).collect(),
            selected: BTreeMap::new(),
            revision: 0,
        }
    }

    fn add_candidates<I>(&mut self, candidates: I) -> Result<usize>
    where I: IntoIterator<Item = SourceCadFragment> {
        self.add_candidates_with_proof(candidates, &BoundaryVerification::Structural)
    }

    fn add_candidates_with_proof<I>(
        &mut self,
        candidates: I,
        boundary_verification: &BoundaryVerification,
    ) -> Result<usize>
    where
        I: IntoIterator<Item = SourceCadFragment>,
    {
        let mut changes = 0;
        for candidate in candidates {
            let SourceCadFragment {
                source_solid: source,
                fragment,
            } = candidate;
            let local_solid_id = source_solid_id(&fragment.source);
            if local_solid_id != source.id {
                bail!(
                    "recovered CAD fragment for DATA section {} solid #{} carries local solid #{local_solid_id}",
                    source.data_section,
                    source.id,
                );
            }
            if !self.source_solids.contains(&source) {
                bail!(
                    "recovered CAD fragment refers to missing source solid {:?} #{}",
                    source.data_section,
                    source.id
                );
            }
            let complexity = fragment.model.complexity_score(fragment.root)?;
            let replace = self.selected.get(&source).is_none_or(|existing| {
                if boundary_verification.rank() == existing.boundary_verification.rank() {
                    whole_solid_candidate_is_preferred(
                        &existing.fragment.source,
                        existing.complexity,
                        &fragment.source,
                        complexity,
                    )
                } else {
                    boundary_verification.rank() > existing.boundary_verification.rank()
                }
            });
            if replace {
                self.selected.insert(
                    source,
                    SelectedWholeSolid {
                        fragment,
                        complexity,
                        boundary_verification: boundary_verification.clone(),
                    },
                );
                self.revision += 1;
                changes += 1;
            }
        }
        Ok(changes)
    }

    fn constructively_unresolved_sources(&self) -> BTreeSet<SourceEntity> {
        self.source_solids
            .iter()
            .copied()
            .filter(|source| {
                self.selected
                    .get(source)
                    .is_none_or(|selected| !is_constructive_source(&selected.fragment.source))
            })
            .collect()
    }

    #[cfg(test)]
    fn semantically_unrecovered_sources(&self) -> BTreeSet<SourceEntity> {
        self.source_solids
            .iter()
            .copied()
            .filter(|source| !self.selected.contains_key(source))
            .collect()
    }

    fn selected_fragments(&self) -> BTreeMap<SourceEntity, &CadFragment> {
        self.selected
            .iter()
            .map(|(&source, selected)| (source, &selected.fragment))
            .collect()
    }
}

#[derive(Debug)]
struct CompleteIrBuildState {
    surface_signature_summaries: Vec<SurfaceSignatureSummary>,
    instance_pattern_diagnostics: Vec<SourceInstancePattern>,
    body_scan: Vec<SourceSolidSummary>,
    recovery: WholeSolidRecoveryState,
    direct_candidates: BTreeMap<SourceEntity, Vec<CadFragment>>,
}

impl CompleteIrBuildState {
    fn finish(
        self,
        schema: &'static str,
        kernel_recovery_failures: Vec<KernelRecoveryFailure>,
        kernel_source_load: Option<KernelSourceLoadEvidence>,
    ) -> CompleteIr {
        let selected = self.recovery.selected_fragments();
        let recovered_ids = selected.keys().copied().collect::<BTreeSet<_>>();
        let unrecovered_solids = self
            .body_scan
            .iter()
            .filter(|solid| !recovered_ids.contains(&solid.source))
            .cloned()
            .collect::<Vec<_>>();
        let constructively_recovered_solid_count = selected
            .values()
            .filter(|fragment| is_constructive_source(&fragment.source))
            .count();

        CompleteIr {
            schema,
            surface_signature_summaries: self.surface_signature_summaries,
            instance_pattern_diagnostics: self.instance_pattern_diagnostics,
            solid_extrusions: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::SolidExtrusion { .. })
            }),
            solid_lofts: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::SolidLoft { .. })
            }),
            solid_revolutions: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::SolidRevolution { .. })
            }),
            radial_slot_revolutions: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::RadialSlotRevolution { .. })
            }),
            solid_sweeps: selected_category(&selected, |source| {
                matches!(
                    source,
                    CadFragmentSource::SolidSweep { .. }
                        | CadFragmentSource::RevolvedRoundTail { .. }
                )
            }),
            periodic_chains: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::PeriodicChain { .. })
            }),
            boundary_feature_composites: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::BoundaryFeatureComposite { .. })
            }),
            edge_treatment_composites: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::EdgeTreatmentComposite { .. })
            }),
            same_domain_normalized_solids: selected_category(&selected, |source| {
                matches!(source, CadFragmentSource::SameDomainNormalized { .. })
            }),
            extrusion_core_unions: selected_category(&selected, |source| {
                matches!(
                    source,
                    CadFragmentSource::ExtrusionCoreUnion { .. }
                        | CadFragmentSource::SectionalExtrusionUnion { .. }
                        | CadFragmentSource::SectionalExtrusionTransitions { .. }
                        | CadFragmentSource::SectionalLoftUnion { .. }
                        | CadFragmentSource::FormedSheetExtrusionUnion { .. }
                        | CadFragmentSource::SourceCapExtrusionUnion { .. }
                        | CadFragmentSource::ExtrusionCoreRuledEnds { .. }
                        | CadFragmentSource::CapExtrusionRuledEnds { .. }
                )
            }),
            kernel_recovery_failures,
            kernel_source_load,
            unrecovered_solid_count: unrecovered_solids.len(),
            body_count: self.body_scan.len(),
            source_boundary_verified_solids: self
                .recovery
                .selected
                .iter()
                .filter_map(|(&source, selected)| {
                    matches!(
                        selected.boundary_verification,
                        BoundaryVerification::Original
                    )
                    .then_some(source)
                })
                .collect(),
            source_orientation_reconciled_solids: self
                .recovery
                .selected
                .iter()
                .filter_map(|(&source, selected)| {
                    let BoundaryVerification::Reconciled(repairs) = &selected.boundary_verification
                    else {
                        return None;
                    };
                    Some(SourceOrientationReconciliation {
                        source,
                        repairs: repairs.clone(),
                    })
                })
                .collect(),
            constructively_recovered_solid_count,
            semantically_recovered_solid_count: selected.len(),
            body_scan: BodyScanSummary {
                solids: self.body_scan,
            },
            unrecovered_solids,
        }
    }
}

type WholeSolidCandidateProvider = for<'a> fn(
    &[RecoverySection<'a>],
    &BTreeSet<SourceEntity>,
    ReducedBrepMutationMode,
) -> Result<Vec<SourceCadFragment>>;

const WHOLE_SOLID_CANDIDATE_PROVIDERS: [WholeSolidCandidateProvider; 2] = [
    recover_boundary_feature_candidates,
    recover_reduced_brep_candidates,
];
const MAX_LOCAL_REDUCTION_ROUNDS: usize = 32;

/// Build the compact complete semantic IR for a normalized STEP exchange.
///
/// This report deliberately separates compact canonical IR from source diagnostics.
/// Every recovered manifold solid is represented by the cheapest proven whole-solid
/// semantic fragment. Bodies without a proven constructive representation are listed
/// explicitly in `unrecovered_solids` and never enter `cad_ir`. Instance-pattern
/// evidence remains diagnostic until its repeated child can be linked to proven CAD IR.
///
/// # Errors
/// Returns an error if STEP parsing fails or a detector's proven result cannot be converted into valid CAD IR.
pub fn recover_complete_ir_bytes(input: &[u8]) -> Result<CompleteIr> {
    let session = crate::AnalysisSession::from_bytes(input)?;
    session.complete_ir()
}

/// Build Complete IR and opportunistically recover otherwise-unrecovered solids
/// through a supplied CAD kernel plus the independent analytic boundary proof.
///
/// The parser-only [`recover_complete_ir_bytes`] path reports schema v10 and never
/// invokes a CAD backend. This explicit kernel-backed path reports its versioned schema,
/// includes parser-neutral source-load evidence,
/// records admitted kernel-owned topology-normalized composites and extrusion-core unions separately, and leaves any unsupported
/// or unproven body in `unrecovered_solids`. Backend failures are diagnostic per
/// solid rather than aborting the rest of the report.
///
/// # Errors
/// Returns an error if STEP parsing/Complete IR construction fails or if the selected kernel cannot inspect the original STEP source.
/// Independently certify each source-derived constructive candidate against
/// the original closed, oriented analytic STEP boundary. Parser recognition is
/// a proposal, not proof of occupied material; rejected candidates do not
/// prevent another valid representation of the same source solid from winning.
fn recover_verified_direct_constructive_fragment<K: CadKernel>(
    proof: &SourceRecovery<'_, K>,
    candidates: Vec<CadFragment>,
    ceiling: Option<u64>,
) -> Result<Option<CadFragment>> {
    let mut candidates = candidates
        .into_iter()
        .map(|fragment| Ok((fragment.model.complexity_score(fragment.root)?, fragment)))
        .collect::<Result<Vec<_>>>()?;
    candidates.sort_by_key(|(complexity, _)| *complexity);
    for (cost, fragment) in candidates {
        if ceiling.is_some_and(|limit| cost >= limit) {
            break;
        }
        if proof
            .proves(&fragment.model, fragment.root, None)
            .unwrap_or(false)
        {
            return Ok(Some(fragment));
        }
    }
    Ok(None)
}

/// Recover the full source document using native kernel ownership, then
/// independently certify executable constructive candidates against STEP geometry
/// when the kernel supports exporting evaluated shapes.
///
/// # Errors
/// Returns an error if source STEP parsing, native inspection, candidate
/// construction, or kernel-backed recovery cannot be initialized. Individual
/// candidate failures remain attached to the returned recovery report.
pub fn recover_complete_ir_with_kernel_bytes<
    K: CadKernel + SourceStepLoadInspector + SourceTopologyEditor,
>(
    input: &[u8],
    kernel: &K,
) -> Result<CompleteIr> {
    let exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
    let document = kernel.resolve_source_exchange(exchange)?;
    let exchange = document.as_ref();
    // Reuse the same section-local evidence for the initial parser proposal
    // pass and every independent proof provider, including multi-DATA inputs.
    let indexes: Vec<_> = exchange
        .data
        .iter()
        .map(|section| crate::step_graph::build_index(&section.entities))
        .collect();
    let source_tolerances: Vec<_> = exchange
        .data
        .iter()
        .zip(&indexes)
        .map(|(section, index)| {
            crate::solid_revolutions::source_tolerance_by_representation_item(
                &section.entities,
                index,
            )
        })
        .collect();
    let sections: Vec<_> = indexes
        .iter()
        .zip(&source_tolerances)
        .map(|(index, source_tolerances)| SectionEvidence {
            index,
            source_tolerances,
            surface_signatures: None,
        })
        .collect();
    let mut build = recover_complete_ir_build_state(
        exchange,
        ReducedBrepMutationMode::KernelOwned,
        Some(&sections),
    )?;
    let source_load = kernel.inspect_source_document(&document)?;
    let mut failures = Vec::new();
    let export = kernel.supports_step_export();
    // All candidates for one source share immutable boundary evidence.
    // A proof-only orientation reconciliation retains its own provenance.
    // Drop it before proceeding to the next solid: do not retain the full
    // assembly's spline signatures merely to reuse the evidence later.
    for source in build.recovery.source_solids.clone() {
        let entities = &exchange.data[source.data_section].entities;
        let index = &indexes[source.data_section];
        let proof = if export {
            SourceRecovery::new(source.id, entities, index, kernel).or_else(|| {
                SourceRecovery::reconcile_strip_senses(source.id, entities, index, kernel)
            })
        } else {
            None
        };
        if !build.recovery.selected.contains_key(&source) {
            let candidate = recover_kernel_boundary_feature_fragment(
                source,
                &source_load,
                &document,
                kernel,
                index,
                proof.as_ref(),
            );
            failures.extend(build.recovery.admit_kernel_candidate(source, candidate));
        }
        if !build.recovery.selected.contains_key(&source) {
            let candidate = recover_kernel_topology_normalized_fragment(
                &document,
                source,
                &source_load,
                kernel,
                index,
                proof.as_ref(),
            );
            failures.extend(build.recovery.admit_kernel_candidate(source, candidate));
        }
        let candidates = build.direct_candidates.remove(&source).unwrap_or_default();
        let Some(proof) = proof else {
            continue;
        };
        // A direct source feature is the canonical representation on an
        // equal score. Prove it before searching alternative decompositions;
        // strictly cheaper certified material candidates can still replace it.
        let candidate = recover_verified_direct_constructive_fragment(
            &proof,
            candidates,
            build.recovery.certified_ceiling(source),
        )
        .map(|candidate| {
            candidate.map(|fragment| KernelRecoveryCandidate::certified(fragment, &proof))
        });
        failures.extend(build.recovery.admit_kernel_candidate(source, candidate));
        let candidate =
            recover_proven_material_fragment(&proof, build.recovery.certified_ceiling(source)).map(
                |candidate| {
                    candidate.map(|fragment| KernelRecoveryCandidate::certified(fragment, &proof))
                },
            );
        failures.extend(build.recovery.admit_kernel_candidate(source, candidate));
    }
    Ok(build.finish(KERNEL_COMPLETE_IR_SCHEMA, failures, Some(source_load)))
}

/// Peel fully proven boundary-feature leaves on kernel-owned topology. The parser
/// constructs tool/selection evidence, never a mutable recovery shape.
fn recover_kernel_boundary_feature_fragment<
    K: CadKernel + SourceStepLoadInspector + SourceTopologyEditor,
>(
    source: SourceEntity,
    source_load: &KernelSourceLoadEvidence,
    source_document: &K::Document,
    kernel: &K,
    index: &HashMap<u64, usize>,
    proof: Option<&SourceRecovery<'_, K>>,
) -> Result<Option<KernelRecoveryCandidate>> {
    let Some(load) = source_load.solids.iter().find(|item| item.source == source) else {
        return Ok(None);
    };
    if load.failure.is_some() {
        return Ok(None);
    }
    let Some(provenance) = load.topology_provenance.as_ref() else {
        return Ok(None);
    };
    let source_exchange = source_document.as_ref();
    let section = source_exchange
        .data
        .get(source.data_section)
        .ok_or_else(|| anyhow::anyhow!("kernel peel lost its source DATA section"))?;
    let mut projected_entities = section.entities.clone();
    // Only shell membership and host FACE_BOUND references are ever changed
    // in this semantic projection; the original entity index stays valid.
    let original_outer_faces =
        crate::brep::solid_face_ids(source.id, &projected_entities, index)
            .ok_or_else(|| anyhow::anyhow!("source solid has no readable outer-shell faces"))?;
    let mut edited = None::<K::EditedSource>;
    let mut replay = Vec::<RecoveredBoundaryFeatureFamily>::new();

    for pass in 0..MAX_LOCAL_REDUCTION_ROUNDS {
        let Some(native_pass) =
            prepare_native_boundary_pass(source, &projected_entities, index, provenance, pass)?
        else {
            break;
        };
        let NativeBoundaryPass {
            patch,
            plan,
            families,
        } = native_pass;
        // The native geometry transaction commits first. Then update only
        // shell membership and FACE_BOUND references in the disposable source
        // projection. Its preflight is atomic, and we never clone/reimport a
        // second geometric source at any pass.
        if let Some(existing) = &mut edited {
            kernel.apply_edited_source_patches(existing, &[patch])?;
        } else {
            edited = Some(kernel.edit_source_topology(source_document, source, &[patch])?);
        }
        crate::planar_features::apply_boundary_feature_peel_plan(
            &mut projected_entities,
            index,
            &plan,
        )?;
        let current = edited
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("kernel peel dropped its edited source"))?;
        validate_native_peel_projection(
            source,
            &original_outer_faces,
            &projected_entities,
            index,
            provenance,
            load.boundaries.len(),
            &kernel.inspect_edited_source(current)?,
        )?;

        // Advance semantic replay only after the authoritative native
        // transaction has passed the source-identity consistency check.
        replay.extend(families.into_iter().map(|mut feature| {
            feature.original_source_ids = pass == 0;
            feature
        }));
        if pass + 1 == MAX_LOCAL_REDUCTION_ROUNDS {
            bail!("kernel boundary-feature peeling did not reach a fixed point");
        }
    }
    let Some(edited) = edited else {
        return Ok(None);
    };
    finish_native_boundary_feature_replay(kernel, source, edited, provenance, replay, proof)
}

struct NativeBoundaryPass {
    patch: crate::kernel_topology::KernelTopologyPatch,
    plan: BoundaryFeaturePeelPlan,
    families: Vec<RecoveredBoundaryFeatureFamily>,
}

fn prepare_native_boundary_pass(
    source: SourceEntity,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    provenance: &KernelTopologyProvenance,
    pass_index: usize,
) -> Result<Option<NativeBoundaryPass>> {
    let targets = LocalRecoveryTargets::from_solids(&BTreeSet::from([source.id]), entities, index)?;
    if targets.is_empty() {
        return Ok(None);
    }
    let shells = targets.shell_ids();
    let report = crate::planar_features::diagnose_boundary_features_for_shells_with_index(
        entities, index, &shells,
    );
    let selected = (0..report.families.len())
        .filter(|&i| {
            targets
                .solid_for_shell(report.families[i].shell_id)
                .is_some()
        })
        .collect::<Vec<_>>();
    if selected.is_empty() {
        return Ok(None);
    }
    let prepared = prepare_boundary_feature_reduction(
        entities, index, &targets, &report, &selected, pass_index,
    )?;
    if prepared.family_indices.is_empty() {
        return Ok(None);
    }
    let (patch, plan) = plan_source_boundary_feature_peel_from_entities(
        source,
        entities,
        index,
        &report,
        &prepared.family_indices,
        provenance,
    )?;
    let families = prepared
        .steps
        .into_iter()
        .map(|(_, LocalReductionStep::BoundaryFeature(feature))| feature)
        .collect();
    Ok(Some(NativeBoundaryPass {
        patch,
        plan,
        families,
    }))
}

/// Ensure projected shell membership exactly matches the retained native
/// `StableIds` after every peel; a semantic-view drift must never be replayed.
fn validate_native_peel_projection(
    source: SourceEntity,
    original_outer_faces: &[u64],
    projected_entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    provenance: &KernelTopologyProvenance,
    boundary_shells: usize,
    native: &crate::kernel_topology::KernelEditedSourceEvidence,
) -> Result<()> {
    let projected_outer_faces =
        crate::brep::solid_face_ids(source.id, projected_entities, index)
            .ok_or_else(|| anyhow::anyhow!("projected source shell is not readable"))?;
    let native_face_ids = |faces: &[u64]| -> Result<BTreeSet<_>> {
        faces
            .iter()
            .map(|&face| {
                provenance
                    .face_for_source_id(source.data_section, face)
                    .ok_or_else(|| {
                        anyhow::anyhow!("source face #{face} has no original native identity")
                    })
            })
            .collect()
    };
    let original_outer_ids = native_face_ids(original_outer_faces)?;
    let projected_outer_ids = native_face_ids(&projected_outer_faces)?;
    if !projected_outer_ids.is_subset(&original_outer_ids) {
        bail!("kernel peel projection introduced source faces not in the original outer shell");
    }
    // Outer-shell faces removed by feature peeling disappear from the native
    // shape. Face identities from every other boundary shell (including voids)
    // must survive unchanged; they are not in the editable source projection.
    let removed_ids = original_outer_ids
        .difference(&projected_outer_ids)
        .copied()
        .collect::<BTreeSet<_>>();
    let expected_ids = provenance
        .faces
        .iter()
        .map(|face| face.topology_id)
        .filter(|id| !removed_ids.contains(id))
        .collect::<BTreeSet<_>>();
    let retained_ids = native
        .retained_source_topology
        .faces
        .iter()
        .map(|face| face.topology_id)
        .collect::<BTreeSet<_>>();
    if native.source != source
        || native.boundary_shells != boundary_shells
        || native.faces != expected_ids.len()
        || expected_ids != retained_ids
    {
        bail!(
            "kernel peel diverged from source face membership: projected_outer={} native_faces={} expected_retained={} actual_retained={}",
            projected_outer_faces.len(),
            native.faces,
            expected_ids.len(),
            retained_ids.len()
        );
    }
    Ok(())
}

fn finish_native_boundary_feature_replay<
    K: CadKernel + SourceStepLoadInspector + SourceTopologyEditor,
>(
    kernel: &K,
    source: SourceEntity,
    mut edited: K::EditedSource,
    provenance: &KernelTopologyProvenance,
    replay: Vec<RecoveredBoundaryFeatureFamily>,
    proof: Option<&SourceRecovery<'_, K>>,
) -> Result<Option<KernelRecoveryCandidate>> {
    let mut journal = native_peels::Journal::default();
    journal.boundary(replay);
    let Some((fragment, evaluated)) = journal.finish(
        kernel,
        &mut edited,
        source,
        provenance,
        native_peels::NextPeel::Treatment,
        None,
    )?
    else {
        return Ok(None);
    };
    // Certify the entire evaluated result against the original, unpeeled
    // source. Unsupported boundaries keep their existing structural status;
    // successful proof also supplies the safe candidate-complexity ceiling.
    Ok(Some(match proof {
        Some(proof)
            if proof
                .matches_evaluated_boundary(&evaluated)
                .unwrap_or(false) =>
        {
            KernelRecoveryCandidate::certified(fragment, proof)
        }
        _ => KernelRecoveryCandidate::structural(fragment),
    }))
}

fn recover_kernel_topology_normalized_fragment<
    K: CadKernel + SourceStepLoadInspector + SourceTopologyEditor,
>(
    document: &K::Document,
    source: SourceEntity,
    source_load: &KernelSourceLoadEvidence,
    kernel: &K,
    index: &HashMap<u64, usize>,
    proof: Option<&SourceRecovery<'_, K>>,
) -> Result<Option<KernelRecoveryCandidate>> {
    let exchange = document.as_ref();
    let Some(load) = source_load
        .solids
        .iter()
        .find(|solid| solid.source == source)
    else {
        return Ok(None);
    };
    if load.failure.is_some() {
        return Ok(None);
    }
    let Some(provenance) = load.topology_provenance.clone() else {
        return Ok(None);
    };

    let evidence = replayable_source_treatments(
        &exchange.data[source.data_section].entities,
        index,
        source.id,
    );
    let treatment_layers::PeelOutcome {
        mut edited,
        layers,
        stop,
    } = treatment_layers::peel(kernel, document, source, &provenance, &evidence)
        .context("peel dependency-ordered native edge treatments")?;
    if kernel.inspect_edited_source(&edited)?.faces == 0 {
        return Ok(None);
    }
    let same_domain_merges = kernel
        .normalize_edited_source(&mut edited)
        .context("normalize native source faces after inverse treatments")?;
    let unchanged = layers.is_empty() && same_domain_merges == 0;
    // Even an unedited source can have a native constructive representation
    // that the parser-only path did not recognize (e.g. an exact U-directed
    // rational spline extrusion). Do not needlessly serialize that source or
    // mislabel direct recovery as a topology-normalization operation.
    if unchanged {
        let Some(mut host) = recognize_native_edited_host(kernel, &edited)? else {
            if let Some(reason) = stop {
                bail!("inverse treatment fixed point: {reason}");
            }
            return Ok(None);
        };
        host.model.provenance.clear();
        return kernel_replay_is_single_solid(kernel, &host)
            .map(|accepted| accepted.then(|| KernelRecoveryCandidate::structural(host)));
    }
    let mut journal = native_peels::Journal::default();
    journal.treatment(layers);
    let Some((mut fragment, evaluated)) = journal.finish(
        kernel,
        &mut edited,
        source,
        &provenance,
        native_peels::NextPeel::Boundary,
        stop,
    )?
    else {
        return Ok(None);
    };
    if !matches!(
        fragment.source,
        CadFragmentSource::EdgeTreatmentComposite { .. }
            | CadFragmentSource::BoundaryFeatureComposite { .. }
    ) {
        fragment.source = CadFragmentSource::SameDomainNormalized {
            solid_id: source.id,
        };
        return Ok(Some(KernelRecoveryCandidate::structural(fragment)));
    }
    if !kernel.supports_step_export() {
        return Ok(Some(KernelRecoveryCandidate::structural(fragment)));
    }
    // A closed replay is only a candidate. Export-capable kernels must
    // certify the reconstructed occupied boundary against the original,
    // not merely the intermediate peeled host.
    let Some(proof) = proof else {
        bail!("native edge-treatment replay lacks a supported source boundary proof");
    };
    if !proof.matches_evaluated_boundary(&evaluated)? {
        bail!("native edge-treatment replay did not prove source boundary equivalence");
    }
    Ok(Some(KernelRecoveryCandidate::certified(fragment, proof)))
}

fn replayable_source_treatments(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    solid_id: u64,
) -> Vec<crate::edge_treatments::EdgeTreatmentEvidence> {
    let target = HashSet::from([solid_id]);
    let mut evidence = select_replayable_treatments(
        crate::edge_treatments::detect_edge_treatment_evidence_for_solids_with_index(
            entities, index, &target,
        ),
    );
    evidence.sort_by_key(|item| item.treatment_face_id);
    evidence
}

/// Prefer native exact-trim semantic recognition while retaining the general
/// STEP-record recognizer for shapes not yet consumed by the native adapter.
fn recognize_edited_host_with_fallback<K: SourceTopologyEditor>(
    kernel: &K,
    edited: &K::EditedSource,
    fallback: impl FnOnce() -> Result<Option<CadFragment>>,
) -> Result<Option<CadFragment>> {
    recognize_native_edited_host(kernel, edited)?.map_or_else(fallback, |host| Ok(Some(host)))
}

/// Recognize exact analytic hosts from the edited native B-rep, with no STEP
/// serialization. Unknown surface/curve classes remain unclassified in the
/// native view; the existing general STEP recognizer handles those separately.
fn recognize_native_edited_host<K: SourceTopologyEditor>(
    kernel: &K,
    edited: &K::EditedSource,
) -> Result<Option<CadFragment>> {
    let reduced = kernel.edited_source_reduced_brep(edited)?;
    if let Some(extrusion) = crate::solid_extrusions::detect_solid_extrusion_from_reduced(&reduced)
    {
        return recover_solid_extrusion_fragment(&extrusion).map(Some);
    }
    if let Some(loft) = crate::solid_lofts::detect_sectioned_planar_loft(&reduced) {
        return Ok(Some(loft));
    }
    crate::solid_revolutions::detect_native_solid_revolution_from_reduced(&reduced)
        .as_ref()
        .map(crate::cad_recovery::recover_solid_revolution_fragment)
        .transpose()
}

fn recognize_single_healed_host<K: CadKernel>(
    healed: &Exchange,
    kernel: &K,
) -> Result<Option<CadFragment>> {
    recognize_single_healed_host_matching(healed, kernel, None)
}

fn recognize_single_healed_host_matching<K: CadKernel>(
    healed: &Exchange,
    kernel: &K,
    acceptance: Option<&CandidateAcceptance<'_>>,
) -> Result<Option<CadFragment>> {
    let [section] = healed.data.as_slice() else {
        return Ok(None);
    };
    let index = crate::step_graph::build_index(&section.entities);
    let solid_ids = section
        .entities
        .iter()
        .filter_map(|entity| {
            let record = crate::step_graph::simple_record(entity)?;
            (record.name == "MANIFOLD_SOLID_BREP").then(|| crate::step_graph::entity_id(entity))
        })
        .collect::<Vec<_>>();
    let [solid_id] = solid_ids.as_slice() else {
        return Ok(None);
    };
    let mut candidates =
        recover_constructive_fragments_from_entities_with_index(&section.entities, &index)?;
    // Monstertruck's STEP re-export is an internal handoff, not a requirement that the raw-entity
    // recognizers understand every equivalent exporter decomposition. ReducedSolid is allowed here only
    // as an immutable analytic view of the topology the kernel has already normalized.
    if let Some(reduced) =
        crate::reduced_brep::ReducedSolid::from_step(*solid_id, &section.entities, &index)
    {
        if let Some(extrusion) =
            crate::solid_extrusions::detect_solid_extrusion_from_reduced(&reduced)
        {
            candidates.push(recover_solid_extrusion_fragment(&extrusion)?);
        }
        if let Some(revolution) = crate::solid_revolutions::detect_solid_revolution_from_reduced(
            &reduced,
            &section.entities,
            &index,
        ) {
            candidates.push(crate::cad_recovery::recover_solid_revolution_fragment(
                &revolution,
            )?);
        }
    }
    let refs = candidates
        .iter()
        .filter(|fragment| acceptance.is_none_or(|accept| accept(&fragment.model, fragment.root)))
        .collect::<Vec<_>>();
    let selected = select_cheapest_whole_solid_fragments(&refs)?;
    if let Some(host) = selected
        .get(solid_id)
        .filter(|host| is_constructive_source(&host.source))
    {
        return Ok(Some((*host).clone()));
    }
    // Apply all independently certified material recognizers to this
    // immutable final native snapshot, not just the sectional family.
    let Some(proof) = SourceRecovery::new(*solid_id, &section.entities, &index, kernel) else {
        return Ok(None);
    };
    let proof = if let Some(accept) = acceptance {
        proof.requiring(accept)
    } else {
        proof
    };
    recover_proven_material_fragment(&proof, None)
}

fn kernel_replay_is_single_solid<K: CadKernel>(kernel: &K, fragment: &CadFragment) -> Result<bool> {
    let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
    Ok(kernel.summarize(&evaluated).is_single_solid())
}

/// Whole-solid evidence belongs to the individual candidate, not to the
/// provider: native topology recovery may prove a treatment composite with
/// an independent source-boundary comparison, or return a structural host.
struct KernelRecoveryCandidate {
    fragment: CadFragment,
    boundary_verification: BoundaryVerification,
}

impl KernelRecoveryCandidate {
    const fn structural(fragment: CadFragment) -> Self {
        Self {
            fragment,
            boundary_verification: BoundaryVerification::Structural,
        }
    }

    fn certified<K: CadKernel>(fragment: CadFragment, proof: &SourceRecovery<'_, K>) -> Self {
        let repairs = proof.strip_face_sense_repairs();
        let boundary_verification = if repairs.is_empty() {
            BoundaryVerification::Original
        } else {
            BoundaryVerification::Reconciled(repairs.to_vec())
        };
        Self {
            fragment,
            boundary_verification,
        }
    }
}

impl WholeSolidRecoveryState {
    fn certified_ceiling(&self, source: SourceEntity) -> Option<u64> {
        self.selected
            .get(&source)
            .filter(|selected| {
                selected.boundary_verification.is_certified()
                    && is_constructive_source(&selected.fragment.source)
            })
            .map(|selected| selected.complexity)
    }

    fn admit_kernel_candidate(
        &mut self,
        source: SourceEntity,
        candidate: Result<Option<KernelRecoveryCandidate>>,
    ) -> Option<KernelRecoveryFailure> {
        let result = match candidate {
            Ok(Some(KernelRecoveryCandidate {
                fragment,
                boundary_verification,
            })) if source_solid_id(&fragment.source) == source.id => self
                .add_candidates_with_proof(
                    [SourceCadFragment {
                        source_solid: source,
                        fragment,
                    }],
                    &boundary_verification,
                )
                .map(|_| ()),
            Ok(Some(KernelRecoveryCandidate { fragment, .. })) => Err(anyhow::anyhow!(
                "kernel recovery returned local solid #{} while recovering DATA section {:?} solid #{}",
                source_solid_id(&fragment.source),
                source.data_section,
                source.id
            )),
            Ok(None) => Ok(()),
            Err(error) => Err(error),
        };
        result.err().map(|error| KernelRecoveryFailure {
            source,
            error: format!("{error:#}"),
        })
    }
}

#[cfg(test)]
fn apply_kernel_recovery<F>(
    recovery: &mut WholeSolidRecoveryState,
    mut recover: F,
    all_sources: bool,
) -> Vec<KernelRecoveryFailure>
where
    F: FnMut(SourceEntity, Option<u64>) -> Result<Option<KernelRecoveryCandidate>>,
{
    let targets = if all_sources {
        recovery.source_solids.clone()
    } else {
        recovery.semantically_unrecovered_sources()
    };
    let mut failures = Vec::new();
    for source in targets {
        let candidate = recover(source, recovery.certified_ceiling(source));
        failures.extend(recovery.admit_kernel_candidate(source, candidate));
    }
    failures
}

/// Borrowed, per-DATA read-only evidence supplied by an `AnalysisSession`.
/// Never qualify signatures using a document-global raw entity ID.
pub(crate) struct SectionEvidence<'a> {
    pub index: &'a HashMap<u64, usize>,
    pub source_tolerances: &'a HashMap<u64, f64>,
    pub surface_signatures: Option<&'a [SolidSurfaceSignature]>,
}

/// Borrow section-local indexes and source tolerances from the owning
/// `AnalysisSession`. Kernel-owned recovery still builds evidence independently.
pub(crate) fn recover_complete_ir_exchange_with_context(
    exchange: &Exchange,
    sections: &[SectionEvidence<'_>],
) -> Result<CompleteIr> {
    anyhow::ensure!(
        exchange.data.len() == sections.len(),
        "one cached evidence pair is required per STEP DATA section"
    );
    recover_complete_ir_exchange_inner(exchange, Some(sections))
}

fn recover_complete_ir_exchange_inner(
    exchange: &Exchange,
    sections: Option<&[SectionEvidence<'_>]>,
) -> Result<CompleteIr> {
    Ok(
        recover_complete_ir_build_state(exchange, ReducedBrepMutationMode::ParserOwned, sections)?
            .finish(COMPLETE_IR_SCHEMA, Vec::new(), None),
    )
}

fn recover_complete_ir_build_state<'a>(
    exchange: &'a Exchange,
    mutation_mode: ReducedBrepMutationMode,
    sections: Option<&[SectionEvidence<'a>]>,
) -> Result<CompleteIrBuildState> {
    let mut signatures = Vec::new();
    let mut patterns = Vec::new();
    let mut body_scan = Vec::new();
    let mut direct_candidates = Vec::new();
    let mut recovery_sections = Vec::with_capacity(exchange.data.len());

    for (data_section, section) in exchange.data.iter().enumerate() {
        let cached = sections.and_then(|items| items.get(data_section));
        let index = cached.map_or_else(
            || Cow::Owned(crate::step_graph::build_index(&section.entities)),
            |context| Cow::Borrowed(context.index),
        );
        let recovery_section = RecoverySection {
            data_section,
            entities: &section.entities,
            index,
        };

        let section_signatures = cached
            .and_then(|context| context.surface_signatures)
            .map_or_else(
                || {
                    Cow::Owned(
                        crate::solid_revolutions::detect_solid_surface_signatures_with_index(
                            recovery_section.entities,
                            &recovery_section.index,
                        ),
                    )
                },
                Cow::Borrowed,
            );
        for signature in section_signatures.iter() {
            let source = recovery_section.source_entity(signature.solid_id).ok_or_else(|| {
                anyhow::anyhow!(
                    "surface signature refers to missing solid #{} in DATA section {data_section}",
                    signature.solid_id
                )
            })?;
            signatures.push(SurfaceSignatureSummary::from_source(source, signature));
        }

        let mut section_candidates = if let Some(context) = cached {
            recover_constructive_fragments_with_tolerances(
                recovery_section.entities,
                &recovery_section.index,
                context.source_tolerances,
                &crate::brep::manifold_solid_ids(recovery_section.entities),
            )?
        } else {
            recover_constructive_fragments_from_entities_with_index(
                recovery_section.entities,
                &recovery_section.index,
            )?
        };
        let periodic_chains = crate::periodic_chains::detect_periodic_chains_with_index(
            recovery_section.entities,
            &recovery_section.index,
        );
        section_candidates.extend(recover_periodic_chain_fragments(&periodic_chains)?);
        direct_candidates.extend(recovery_section.qualify_fragments(section_candidates)?);

        patterns.extend(
            crate::patterns::detect_instance_patterns_with_index(
                recovery_section.entities,
                &recovery_section.index,
                DEFAULT_STEP_RECOVERY_TOLERANCES.pattern_geometry_mm,
                4,
            )
            .into_iter()
            .map(|pattern| SourceInstancePattern {
                data_section,
                pattern,
            }),
        );
        body_scan.extend(scan_source_solids(&recovery_section));
        recovery_sections.push(recovery_section);
    }

    body_scan.sort_by_key(|solid| solid.source);
    signatures.sort_by_key(|signature| signature.source);

    let surface_signature_summaries = signatures;
    let instance_pattern_diagnostics = patterns;
    let (recovery, direct_candidates) = recover_and_retain_direct_candidates(
        &recovery_sections,
        &body_scan,
        direct_candidates,
        mutation_mode,
    )?;

    Ok(CompleteIrBuildState {
        surface_signature_summaries,
        instance_pattern_diagnostics,
        body_scan,
        recovery,
        direct_candidates,
    })
}

type DirectCandidatesBySource = BTreeMap<SourceEntity, Vec<CadFragment>>;

fn recover_and_retain_direct_candidates(
    recovery_sections: &[RecoverySection<'_>],
    body_scan: &[SourceSolidSummary],
    direct_candidates: Vec<SourceCadFragment>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<(WholeSolidRecoveryState, DirectCandidatesBySource)> {
    let result = if mutation_mode == ReducedBrepMutationMode::KernelOwned {
        let recovery = recover_whole_solids(
            recovery_sections,
            body_scan,
            direct_candidates.iter().cloned(),
            mutation_mode,
        )?;
        let mut by_source = BTreeMap::<SourceEntity, Vec<CadFragment>>::new();
        for candidate in direct_candidates {
            if is_constructive_source(&candidate.fragment.source) {
                by_source
                    .entry(candidate.source_solid)
                    .or_default()
                    .push(candidate.fragment);
            }
        }
        (recovery, by_source)
    } else {
        (
            recover_whole_solids(
                recovery_sections,
                body_scan,
                direct_candidates,
                mutation_mode,
            )?,
            BTreeMap::new(),
        )
    };

    Ok(result)
}

fn recover_whole_solids(
    sections: &[RecoverySection<'_>],
    body_scan: &[SourceSolidSummary],
    direct_candidates: impl IntoIterator<Item = SourceCadFragment>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<WholeSolidRecoveryState> {
    let mut recovery = WholeSolidRecoveryState::new(body_scan);
    recovery.add_candidates(direct_candidates)?;

    for provider in WHOLE_SOLID_CANDIDATE_PROVIDERS {
        let targets = recovery.constructively_unresolved_sources();
        if targets.is_empty() {
            break;
        }
        recovery.add_candidates(provider(sections, &targets, mutation_mode)?)?;
    }
    Ok(recovery)
}

fn recover_boundary_feature_candidates(
    sections: &[RecoverySection<'_>],
    targets: &BTreeSet<SourceEntity>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<Vec<SourceCadFragment>> {
    if mutation_mode == ReducedBrepMutationMode::KernelOwned {
        return Ok(Vec::new());
    }
    let mut candidates = Vec::new();
    for section in sections {
        let local_targets = section.local_target_ids(targets);
        if local_targets.is_empty() {
            continue;
        }
        candidates.extend(section.qualify_fragments(
            recover_boundary_feature_composites_for_section(
                section,
                &local_targets,
                mutation_mode,
            )?,
        )?);
    }
    Ok(candidates)
}

fn recover_reduced_brep_candidates(
    sections: &[RecoverySection<'_>],
    targets: &BTreeSet<SourceEntity>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<Vec<SourceCadFragment>> {
    let mut candidates = Vec::new();
    for section in sections {
        let local_targets = section.local_target_ids(targets);
        if local_targets.is_empty() {
            continue;
        }
        candidates.extend(section.qualify_fragments(recover_reduced_brep_fragments(
            section.entities,
            &section.index,
            Some(&local_targets),
            mutation_mode,
        )?)?);
    }
    Ok(candidates)
}

fn recover_boundary_feature_composites_for_section(
    section: &RecoverySection<'_>,
    target_solid_ids: &BTreeSet<u64>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<Vec<CadFragment>> {
    let entities = section.entities;
    let initial_index = &section.index;
    let mut representative_targets = BTreeSet::new();
    let mut recovery_groups = HashMap::new();
    // Canonical ShapeKey equality requires equal face counts. Build that
    // inexpensive invariant first: unique counts cannot have a matching
    // source solid to instance, so they do not need full geometry/topology
    // fingerprinting. Ambiguous counts still receive the exact identity
    // proof used previously; this also covers the common singleton case.
    let mut by_face_count = BTreeMap::<usize, Vec<u64>>::new();
    for &solid_id in target_solid_ids {
        if let Some(faces) =
            crate::shape_identity::manifold_solid_face_ids(solid_id, entities, initial_index)
        {
            by_face_count.entry(faces.len()).or_default().push(solid_id);
        } else {
            // Full solid_identity would also reject an unresolved face list.
            representative_targets.insert(solid_id);
        }
    }
    let mut keyed_targets = BTreeMap::<crate::shape_identity::ShapeKey, Vec<_>>::new();
    for face_group in by_face_count.into_values() {
        if let [solid_id] = face_group.as_slice() {
            representative_targets.insert(*solid_id);
            continue;
        }
        for solid_id in face_group {
            if let Some(identity) =
                crate::shape_identity::solid_identity(solid_id, entities, initial_index)
            {
                keyed_targets
                    .entry(identity.key.clone())
                    .or_default()
                    .push(identity);
            } else {
                representative_targets.insert(solid_id);
            }
        }
    }
    for (_, mut group) in keyed_targets {
        group.sort_by_key(|identity| identity.root);
        let Some(representative) = group.first().map(|identity| identity.root) else {
            continue;
        };
        representative_targets.insert(representative);
        recovery_groups.insert(representative, group);
    }

    let mut active_targets =
        LocalRecoveryTargets::from_solids(&representative_targets, entities, initial_index)?;
    if active_targets.is_empty() {
        return Ok(Vec::new());
    }

    let mut journal = ParserPeelJournal::default();
    let recovered = recover_local_reductions_to_fixed_point(
        &mut journal,
        entities,
        initial_index,
        &mut active_targets,
        mutation_mode,
    )?;
    let out = recovered.into_values().collect::<Vec<_>>();

    let mut expanded = Vec::new();
    for fragment in out {
        let representative = source_solid_id(&fragment.source);
        if let Some(group) = recovery_groups.get(&representative)
            && let Some(source_identity) = group.first()
        {
            for target_identity in group.iter().skip(1) {
                if let Some(instance) = instantiate_boundary_feature_fragment(
                    &fragment,
                    source_identity,
                    target_identity,
                    entities,
                    initial_index,
                )? {
                    expanded.push(instance);
                }
            }
        }
        expanded.push(fragment);
    }
    Ok(expanded)
}

fn recover_local_reductions_to_fixed_point(
    journal: &mut ParserPeelJournal,
    source: &[EntityInstance],
    index: &HashMap<u64, usize>,
    active_targets: &mut LocalRecoveryTargets,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<BTreeMap<u64, CadFragment>> {
    let mut recovered = BTreeMap::new();
    for _ in 0..MAX_LOCAL_REDUCTION_ROUNDS {
        let view = journal.view(source, index)?;
        if !recover_boundary_feature_pass(journal, view.as_ref(), index, active_targets)? {
            return Ok(recovered);
        }
        let healed = journal.view(source, index)?;
        for fragment in recover_replayed_boundary_hosts(
            journal,
            healed.as_ref(),
            index,
            &active_targets.solid_ids,
            mutation_mode,
        )? {
            recovered.insert(source_solid_id(&fragment.source), fragment);
        }
        active_targets.retain_unrecovered(&recovered);
        if active_targets.is_empty() {
            return Ok(recovered);
        }
    }
    bail!("local recovery did not reach a fixed point after {MAX_LOCAL_REDUCTION_ROUNDS} rounds")
}

fn recover_replayed_boundary_hosts(
    journal: &ParserPeelJournal,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    target_solid_ids: &BTreeSet<u64>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<Vec<CadFragment>> {
    if target_solid_ids.is_empty() {
        return Ok(Vec::new());
    }
    let mut host_candidates =
        recover_constructive_fragments_from_entities_with_index(entities, index)?;
    host_candidates.extend(recover_reduced_brep_fragments(
        entities,
        index,
        Some(target_solid_ids),
        mutation_mode,
    )?);
    let host_refs = host_candidates.iter().collect::<Vec<_>>();
    let selected_hosts = select_cheapest_whole_solid_fragments(&host_refs)?;
    let mut out = Vec::new();
    for &solid_id in target_solid_ids {
        let families = journal.boundary_features(solid_id);
        if families.is_empty() {
            continue;
        }
        let Some(host) = selected_hosts.get(&solid_id) else {
            continue;
        };
        if !is_constructive_source(&host.source) {
            continue;
        }
        out.push(compose_boundary_feature_fragment(
            solid_id, host, &families,
        )?);
    }
    Ok(out)
}

fn instantiate_boundary_feature_fragment(
    fragment: &CadFragment,
    source_identity: &crate::shape_identity::SolidIdentity,
    target_identity: &crate::shape_identity::SolidIdentity,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Result<Option<CadFragment>> {
    let CadFragmentSource::BoundaryFeatureComposite {
        passes,
        feature_families,
        feature_instances,
        ..
    } = fragment.source
    else {
        return Ok(None);
    };
    let Some(quarter) = crate::shape_identity::relative_z90_quarter(
        source_identity,
        target_identity,
        entities,
        index,
    ) else {
        return Ok(None);
    };
    let translation =
        relative_quarter_translation(source_identity.center, target_identity.center, quarter);
    let mut model = fragment.model.clone();
    model.roots.clear();
    let root = model.add_node(CadNode::Transform {
        transform: quarter_turn_transform(quarter, translation),
        child: fragment.root,
    });
    model.set_provenance(
        root,
        Provenance {
            source_entity_ids: vec![target_identity.root],
            proof: ProofStatus::StructurallyProven,
            max_residual_mm: None,
        },
    )?;
    model.add_root(root)?;
    model.validate()?;
    Ok(Some(CadFragment {
        source: CadFragmentSource::BoundaryFeatureComposite {
            solid_id: target_identity.root,
            passes,
            feature_families,
            feature_instances,
        },
        model,
        root,
    }))
}

fn relative_quarter_translation(
    source_center: [f64; 3],
    target_center: [f64; 3],
    quarter: u8,
) -> [f64; 3] {
    let (x, y) = match quarter % 4 {
        0 => (source_center[0], source_center[1]),
        1 => (-source_center[1], source_center[0]),
        2 => (-source_center[0], -source_center[1]),
        _ => (source_center[1], -source_center[0]),
    };
    [
        target_center[0] - x,
        target_center[1] - y,
        target_center[2] - source_center[2],
    ]
}

fn recover_boundary_feature_pass(
    journal: &mut ParserPeelJournal,
    current: &[EntityInstance],
    index: &HashMap<u64, usize>,
    targets: &LocalRecoveryTargets,
) -> Result<bool> {
    let pass_index = journal.next_boundary_pass_index();
    let target_shell_ids = targets.shell_ids();
    let report = crate::planar_features::diagnose_boundary_features_for_shells_with_index(
        current,
        index,
        &target_shell_ids,
    );
    let target_family_indices = (0..report.families.len())
        .filter(|&candidate| {
            targets
                .shell_to_solid
                .contains_key(&report.families[candidate].shell_id)
        })
        .collect::<Vec<_>>();
    if target_family_indices.is_empty() {
        return Ok(false);
    }

    let prepared = prepare_boundary_feature_reduction(
        current,
        index,
        targets,
        &report,
        &target_family_indices,
        pass_index,
    )?;
    if prepared.family_indices.is_empty() {
        return Ok(false);
    }
    let plan = crate::planar_features::plan_boundary_feature_family_peel(
        current,
        index,
        &report,
        &prepared.family_indices,
    )?;
    journal.commit(prepared.steps, plan, current, index)
}

fn prepare_boundary_feature_reduction(
    source_entities: &[EntityInstance],
    source_index: &HashMap<u64, usize>,
    targets: &LocalRecoveryTargets,
    report: &crate::planar_features::BoundaryFeatureDiagnostics,
    target_family_indices: &[usize],
    pass_index: usize,
) -> Result<PreparedBoundaryReduction> {
    // This clone is the disposable proof-only tool materialization. The source
    // is immutable, and neither consumer receives a mutable topology owner.
    let mut tool_entities = source_entities.to_vec();
    let materialized = crate::planar_features::materialize_boundary_feature_tools_from_report(
        &mut tool_entities,
        source_index,
        report,
        1,
    )?;
    let tool_ids = materialized
        .iter()
        .map(|tool| tool.solid_id)
        .collect::<Vec<_>>();
    let tool_index = crate::step_graph::build_index(&tool_entities);
    let tolerances = crate::solid_revolutions::source_tolerance_by_representation_item(
        &tool_entities,
        &tool_index,
    );
    let tool_candidates = recover_constructive_fragments_with_tolerances(
        &tool_entities,
        &tool_index,
        &tolerances,
        &tool_ids,
    )?;
    let tool_refs = tool_candidates.iter().collect::<Vec<_>>();
    let selected_tools = select_cheapest_whole_solid_fragments(&tool_refs)?;
    let mut peelable_family_indices = Vec::new();
    let mut pending_steps = Vec::new();

    for &family_index in target_family_indices {
        let family = &report.families[family_index];
        let Some(source_solid_id) = targets.solid_for_shell(family.shell_id) else {
            continue;
        };
        let Some(source_member) = family.members.first().cloned() else {
            continue;
        };
        let Some(materialization) = materialized
            .iter()
            .find(|tool| tool.family_index == family_index)
        else {
            continue;
        };
        let Some(tool) = selected_tools.get(&materialization.solid_id) else {
            continue;
        };
        if !is_constructive_source(&tool.source) {
            continue;
        }
        let tool = (*tool).clone();
        if !validate_boundary_feature_tool(&tool, family, source_entities, source_index) {
            continue;
        }
        pending_steps.push((
            source_solid_id,
            LocalReductionStep::BoundaryFeature(RecoveredBoundaryFeatureFamily {
                pass_index,
                polarity: family.polarity,
                tool,
                source_member,
                members: family.members.clone(),
                original_source_ids: true,
            }),
        ));
        peelable_family_indices.push(family_index);
    }
    let selected = select_closed_feature_families(&report.families, &peelable_family_indices);
    let (family_indices, steps) = peelable_family_indices
        .into_iter()
        .zip(pending_steps)
        .filter(|(index, _)| selected.contains(index))
        .unzip();
    Ok(PreparedBoundaryReduction {
        steps,
        family_indices,
    })
}

/// Closed constructive tools can already account for nested feature hypotheses.
/// Only defer a parent if its dependent hypothesis extends outside its patch.
fn select_closed_feature_families(
    families: &[crate::planar_features::BoundaryFeatureFamilyEvidence],
    proven: &[usize],
) -> BTreeSet<usize> {
    let faces = families
        .iter()
        .map(|family| {
            family
                .members
                .iter()
                .flat_map(|member| member.face_ids.iter().copied())
                .collect::<BTreeSet<_>>()
        })
        .collect::<Vec<_>>();
    let hosts = families
        .iter()
        .map(|family| {
            family
                .host_face_ids
                .iter()
                .copied()
                .collect::<BTreeSet<_>>()
        })
        .collect::<Vec<_>>();
    let mut eligible = proven
        .iter()
        .copied()
        .filter(|&i| {
            families.iter().enumerate().all(|(j, _)| {
                i == j || faces[i].is_disjoint(&hosts[j]) || faces[j].is_subset(&faces[i])
            })
        })
        .collect::<Vec<_>>();
    eligible.sort_by_key(|&i| (std::cmp::Reverse(faces[i].len()), i));
    let mut selected = BTreeSet::new();
    let mut removed = BTreeSet::new();
    let mut retained = BTreeSet::new();
    for i in eligible {
        if faces[i].is_disjoint(&removed)
            && faces[i].is_disjoint(&retained)
            && hosts[i].is_disjoint(&removed)
        {
            removed.extend(&faces[i]);
            retained.extend(&hosts[i]);
            selected.insert(i);
        }
    }
    selected
}

fn validate_boundary_feature_tool(
    tool: &CadFragment,
    family: &crate::planar_features::BoundaryFeatureFamilyEvidence,
    source_entities: &[EntityInstance],
    source_index: &HashMap<u64, usize>,
) -> bool {
    if !matches!(tool.source, CadFragmentSource::SolidExtrusion { .. }) {
        return true;
    }
    let Some(host_outwards) =
        feature_host_outward_directions(&family.host_face_ids, source_entities, source_index)
    else {
        return false;
    };
    validate_boundary_tool_for_boolean(tool, &host_outwards).is_ok()
}

fn recover_constructive_fragments_from_entities_with_index(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Result<Vec<CadFragment>> {
    let source_tolerances =
        crate::solid_revolutions::source_tolerance_by_representation_item(entities, index);
    recover_constructive_fragments_with_tolerances(
        entities,
        index,
        &source_tolerances,
        &crate::brep::manifold_solid_ids(entities),
    )
}

fn recover_constructive_fragments_with_tolerances(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    source_tolerances: &HashMap<u64, f64>,
    solid_ids: &[u64],
) -> Result<Vec<CadFragment>> {
    let extrusions =
        crate::solid_extrusions::detect_solid_extrusions_with_index(entities, index, solid_ids);
    let revolutions = crate::solid_revolutions::detect_solid_revolutions_with_tolerances(
        entities,
        index,
        source_tolerances,
        solid_ids,
    );
    let radial_slots = crate::solid_revolutions::detect_radial_slot_revolutions_with_tolerances(
        entities,
        index,
        source_tolerances,
        solid_ids,
    );
    let closed_sweeps = crate::solid_sweeps::detect_round_sweeps(entities, index, solid_ids);
    let revolved_round_tails =
        crate::solid_sweeps::detect_revolved_round_tails(entities, index, solid_ids);
    let open_sweeps =
        crate::solid_sweeps::sheet::detect_open_rectangular_sweeps(entities, index, solid_ids);

    let mut out = recover_solid_extrusion_fragments(&extrusions)?;
    out.extend(recover_solid_revolution_fragments(&revolutions)?);
    out.extend(crate::solid_lofts::detect_sectioned_planar_lofts(
        entities, index, solid_ids,
    ));
    out.extend(recover_radial_slot_revolution_fragments(&radial_slots)?);
    out.extend(recover_round_sweep_fragments(&closed_sweeps)?);
    out.extend(recover_revolved_round_tail_fragments(
        &revolved_round_tails,
    )?);
    out.extend(recover_open_rectangular_sweep_fragments(&open_sweeps)?);
    Ok(out)
}

fn recover_reduced_brep_fragments(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    target_solid_ids: Option<&BTreeSet<u64>>,
    mutation_mode: ReducedBrepMutationMode,
) -> Result<Vec<CadFragment>> {
    if mutation_mode == ReducedBrepMutationMode::KernelOwned {
        return Ok(Vec::new());
    }
    let candidate_solid_ids = entities
        .iter()
        .filter_map(|entity| {
            let record = crate::step_graph::simple_record(entity)?;
            (record.name == "MANIFOLD_SOLID_BREP").then(|| crate::step_graph::entity_id(entity))
        })
        .filter(|solid_id| target_solid_ids.is_none_or(|targets| targets.contains(solid_id)))
        .collect::<HashSet<_>>();
    if candidate_solid_ids.is_empty() {
        return Ok(Vec::new());
    }

    let evidence = select_replayable_treatments(
        crate::edge_treatments::detect_edge_treatment_evidence_for_solids_with_index(
            entities,
            index,
            &candidate_solid_ids,
        ),
    );
    let mut by_solid = BTreeMap::<u64, Vec<_>>::new();
    for item in &evidence {
        by_solid.entry(item.solid_id).or_default().push(item);
    }

    let mut out = Vec::new();
    for solid_id in candidate_solid_ids {
        let Some(mut reduced) =
            crate::reduced_brep::ReducedSolid::from_step(solid_id, entities, index)
        else {
            continue;
        };
        let mut replay = Vec::new();
        if let Some(items) = by_solid.get_mut(&solid_id) {
            items.sort_by_key(|item| item.treatment_face_id);
            for item in items.iter().copied() {
                let Some((profile, edge)) = replayable_edge_treatment(item) else {
                    continue;
                };
                if reduced.remove_edge_treatment(item).is_some() {
                    replay.push((item, profile, edge));
                }
            }
        }
        let coalesced = reduced.coalesce_coplanar_faces();
        if replay.is_empty() && !coalesced {
            continue;
        }

        if let Some(extrusion) =
            crate::solid_extrusions::detect_solid_extrusion_from_reduced(&reduced)
        {
            let fragment = recover_solid_extrusion_fragment(&extrusion)?;
            out.push(replay_edge_treatments(fragment, solid_id, replay.clone())?);
        }
        if let Some(revolution) = crate::solid_revolutions::detect_solid_revolution_from_reduced(
            &reduced, entities, index,
        ) {
            let fragment = crate::cad_recovery::recover_solid_revolution_fragment(&revolution)?;
            out.push(replay_edge_treatments(fragment, solid_id, replay)?);
        }
    }
    Ok(out)
}

type ReplayTreatment<'a> = (
    &'a crate::edge_treatments::EdgeTreatmentEvidence,
    EdgeTreatmentProfile,
    EdgeTreatmentEdge,
);

fn replay_edge_treatments(
    fragment: CadFragment,
    solid_id: u64,
    replay: Vec<ReplayTreatment<'_>>,
) -> Result<CadFragment> {
    replay_grouped_edge_treatments(fragment, solid_id, replay, Vec::new())
}

/// Owned forward operation. Geometry is captured in the inverse layer's
/// current coordinate space, while provenance is independently mapped back
/// to original STEP member identities.
#[derive(Clone)]
struct TreatmentReplayGroup {
    face_ids: Vec<u64>,
    profile: EdgeTreatmentProfile,
    edges: Vec<EdgeTreatmentEdge>,
    max_residual_mm: f64,
}

fn replay_grouped_edge_treatments(
    fragment: CadFragment,
    solid_id: u64,
    replay: Vec<ReplayTreatment<'_>>,
    coupled_groups: Vec<crate::kernel_topology::CoupledTreatmentReplay>,
) -> Result<CadFragment> {
    apply_treatment_replay_groups(
        fragment,
        solid_id,
        make_treatment_replay_groups(replay, coupled_groups)?,
    )
}

/// Every connected perimeter is one operation, not a sequence of
/// independently replayed rounds that would destroy their common junctions.
fn make_treatment_replay_groups(
    replay: Vec<ReplayTreatment<'_>>,
    coupled_groups: Vec<crate::kernel_topology::CoupledTreatmentReplay>,
) -> Result<Vec<TreatmentReplayGroup>> {
    let mut remaining = BTreeMap::new();
    for item in replay {
        anyhow::ensure!(
            remaining.insert(item.0.treatment_face_id, item).is_none(),
            "multiple forward treatment interpretations for the same face"
        );
    }
    let mut groups = Vec::<TreatmentReplayGroup>::new();
    for planned in coupled_groups {
        let mut group = planned
            .selectors
            .into_iter()
            .map(|(id, selector)| {
                let mut item = remaining
                    .remove(&id)
                    .with_context(|| format!("coupled replay has no source treatment #{id}"))?;
                item.2.selector = selector;
                Ok(item)
            })
            .collect::<Result<Vec<_>>>()?;
        anyhow::ensure!(!group.is_empty(), "empty coupled replay group");
        group.sort_by_key(|(item, _, _)| item.treatment_face_id);
        let mut prepared = prepare_treatment_replay_group(group)?;
        if !planned.junction_faces.is_empty() {
            prepared.max_residual_mm = prepared
                .max_residual_mm
                .max(crate::numeric::NUMERICAL_GEOMETRY_TOL_MM);
        }
        prepared.face_ids.extend(planned.junction_faces);
        prepared.face_ids.sort_unstable();
        anyhow::ensure!(
            prepared.face_ids.windows(2).all(|pair| pair[0] != pair[1]),
            "coupled treatment contains a duplicated junction"
        );
        groups.push(prepared);
    }
    for item in remaining.into_values() {
        groups.push(prepare_treatment_replay_group(vec![item])?);
    }
    groups.sort_by_key(|group| group.face_ids[0]);
    Ok(groups)
}

fn prepare_treatment_replay_group(group: Vec<ReplayTreatment<'_>>) -> Result<TreatmentReplayGroup> {
    let profile = group[0].1;
    anyhow::ensure!(
        group.iter().all(|(_, p, _)| *p == profile),
        "one coupled treatment perimeter has mixed forward profiles"
    );
    let mut face_ids = Vec::with_capacity(group.len());
    let mut max_residual_mm = 0.0_f64;
    let mut edges = Vec::<EdgeTreatmentEdge>::new();
    for (item, _, edge) in group {
        face_ids.push(item.treatment_face_id);
        max_residual_mm = max_residual_mm.max(item.max_boundary_residual_mm);
        if let Some(existing) = edges
            .iter()
            .find(|existing| existing.selector == edge.selector)
        {
            anyhow::ensure!(
                crate::numeric::near_imported_feature_mm(
                    (existing.radius_mm - edge.radius_mm).abs()
                ),
                "split treatment has conflicting forward dimensions"
            );
        } else {
            edges.push(edge);
        }
    }
    Ok(TreatmentReplayGroup {
        face_ids,
        profile,
        edges,
        max_residual_mm,
    })
}

fn apply_treatment_replay_groups(
    mut fragment: CadFragment,
    solid_id: u64,
    groups: impl IntoIterator<Item = TreatmentReplayGroup>,
) -> Result<CadFragment> {
    let mut treatment_face_ids = BTreeSet::new();
    fragment.model.roots.clear();
    let mut root = fragment.root;
    for group in groups {
        for &id in &group.face_ids {
            anyhow::ensure!(
                treatment_face_ids.insert(id),
                "source treatment replayed twice"
            );
        }
        root = fragment.model.add_node(CadNode::EdgeTreatment {
            child: root,
            profile: group.profile,
            edges: group.edges,
        });
        let mut source_entity_ids = vec![solid_id];
        source_entity_ids.extend(group.face_ids);
        fragment.model.set_provenance(
            root,
            Provenance {
                source_entity_ids,
                proof: ProofStatus::WithinTolerance,
                max_residual_mm: Some(group.max_residual_mm),
            },
        )?;
    }
    fragment.model.add_root(root)?;
    fragment.root = root;
    if !treatment_face_ids.is_empty() {
        fragment.source = CadFragmentSource::EdgeTreatmentComposite {
            solid_id,
            treatment_face_ids: treatment_face_ids.into_iter().collect(),
        };
    }
    fragment.model.validate()?;
    Ok(fragment)
}

/// Local analytic evidence is not proof that a face should be removed.
/// A symmetric ring of mutually interchangeable planar candidates is a host
/// profile, not a simultaneous batch of edge treatments. Keep its authored
/// profile intact, rather than choosing alternating chamfers by source ID.
fn select_replayable_treatments(
    mut evidence: Vec<crate::edge_treatments::EdgeTreatmentEvidence>,
) -> Vec<crate::edge_treatments::EdgeTreatmentEvidence> {
    evidence.retain(|item| replayable_edge_treatment(item).is_some());
    let profile_faces = crate::edge_treatments::reciprocal_planar_profile_faces(&evidence);
    evidence.retain(|item| !profile_faces.contains(&(item.solid_id, item.treatment_face_id)));
    let mut counts = HashMap::new();
    for item in &evidence {
        *counts
            .entry((item.solid_id, item.treatment_face_id))
            .or_insert(0_usize) += 1;
    }
    evidence.retain(|item| counts[&(item.solid_id, item.treatment_face_id)] == 1);
    evidence
}

fn replayable_edge_treatment(
    item: &crate::edge_treatments::EdgeTreatmentEvidence,
) -> Option<(EdgeTreatmentProfile, EdgeTreatmentEdge)> {
    if !item.carrier_pair_unique {
        return None;
    }
    let selector = match &item.virtual_edge {
        crate::edge_treatments::EdgeTreatmentVirtualEdge::Segment { start_mm, end_mm } => {
            EdgeSelector::Segment {
                start_mm: *start_mm,
                end_mm: *end_mm,
            }
        }
        crate::edge_treatments::EdgeTreatmentVirtualEdge::Circle {
            center_mm,
            normal,
            radius_mm,
        } => EdgeSelector::Circle {
            center_mm: *center_mm,
            normal: *normal,
            radius_mm: *radius_mm,
        },
    };
    let (profile, radius_mm) = match item.kind {
        crate::edge_treatments::EdgeTreatmentKind::PlanarChamfer => {
            let [first, second] = item.inferred_setbacks_mm?;
            if !first.is_finite()
                || !second.is_finite()
                || first <= 0.0
                || second <= 0.0
                || !crate::numeric::near_imported_feature_mm((first - second).abs())
            {
                return None;
            }
            // Numerical point tolerance is not a symmetry certificate:
            // a tiny setback difference over a short bevel can alter its
            // plane normal far beyond the direction quantum. Coalesce only
            // arithmetic roundoff, not measured geometric asymmetry.
            if (first - second).abs() <= f64::EPSILON * first.max(second) {
                (EdgeTreatmentProfile::Chamfer, first.midpoint(second))
            } else {
                // Preserve the source carrier position, rather than fitting a
                // nominal common size that moves later dependent selectors.
                (
                    EdgeTreatmentProfile::CarrierChamfer {
                        first_carrier_normal: item.first_setback_carrier_normal?,
                        second_setback_mm: second,
                    },
                    first,
                )
            }
        }
        crate::edge_treatments::EdgeTreatmentKind::CylindricalFillet
        | crate::edge_treatments::EdgeTreatmentKind::ToroidalFillet => {
            let radius = item.inferred_radius_mm?;
            if !radius.is_finite() || radius <= 0.0 {
                return None;
            }
            (EdgeTreatmentProfile::Round, radius)
        }
        crate::edge_treatments::EdgeTreatmentKind::ConicalChamfer => return None,
    };
    Some((
        profile,
        EdgeTreatmentEdge {
            selector,
            radius_mm,
        },
    ))
}

fn feature_host_outward_directions(
    host_face_ids: &[u64],
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<Vec<[f64; 3]>> {
    let mut out = Vec::with_capacity(host_face_ids.len());
    for &face_id in host_face_ids {
        let surface_id = crate::brep::face_surface(face_id, entities, index)?;
        let crate::brep::SurfaceSupport::Plane(plane) =
            crate::brep::surface_support(surface_id, entities, index)
        else {
            return None;
        };
        let sign = match crate::brep::face_sense(face_id, entities, index)?.as_str() {
            "T" => 1.0,
            "F" => -1.0,
            _ => return None,
        };
        out.push([
            sign * plane.normal[0],
            sign * plane.normal[1],
            sign * plane.normal[2],
        ]);
    }
    (!out.is_empty()).then_some(out)
}

/// The healed host and the exact closed tool already describe the occupied
/// material. Preserve both source cap stations; do not enlarge the tool just
/// to force volumetric overlap.
fn validate_boundary_tool_for_boolean(
    tool: &CadFragment,
    host_outwards: &[[f64; 3]],
) -> Result<()> {
    let CadNode::Transform { transform, child } = tool.model.node(tool.root)? else {
        bail!("boundary-feature extrusion tool root is not a rigid transform");
    };
    let z_axis = [
        transform.matrix[0][2],
        transform.matrix[1][2],
        transform.matrix[2][2],
    ];
    let norm = crate::math3::norm(z_axis);
    if !norm.is_finite()
        || (norm - 1.0).abs() > DEFAULT_STEP_RECOVERY_TOLERANCES.source_identity_axis_unit
    {
        bail!("boundary-feature extrusion tool axis is not unit length");
    }
    for &outward in host_outwards {
        let length = crate::math3::norm(outward);
        if !length.is_finite() || length <= 0.0 {
            bail!("boundary-feature host normal is invalid");
        }
        let alignment = crate::math3::dot(outward, z_axis) / length;
        if (alignment.abs() - 1.0).abs()
            > DEFAULT_STEP_RECOVERY_TOLERANCES.source_axis_direction_cosine
        {
            bail!("boundary-feature interface is not an extrusion end cap");
        }
    }
    let CadNode::Extrude { vector_mm, .. } = tool.model.node(*child)? else {
        bail!("boundary-feature tool transform does not wrap an extrusion");
    };
    if vector_mm[0].abs() > DEFAULT_STEP_RECOVERY_TOLERANCES.source_identity_transform
        || vector_mm[1].abs() > DEFAULT_STEP_RECOVERY_TOLERANCES.source_identity_transform
        || !vector_mm[2].is_finite()
        || vector_mm[2] <= 0.0
    {
        bail!("boundary-feature extrusion vector is not canonical +Z");
    }
    tool.model.validate()
}

fn compose_boundary_feature_fragment(
    source_solid_id: u64,
    host: &CadFragment,
    families: &[RecoveredBoundaryFeatureFamily],
) -> Result<CadFragment> {
    let mut model = CadModel::new();
    let host_offset = model.append_model(&host.model)?;
    let mut current = shifted_fragment_root(host.root, host_offset)?;

    let mut batch = None::<FeatureReplayBatch>;
    for family in families.iter().rev() {
        let tool_offset = model.append_model(&family.tool.model)?;
        let tool_root = shifted_fragment_root(family.tool.root, tool_offset)?;
        let exemplar_sources = if family.original_source_ids {
            feature_member_source_ids(source_solid_id, &family.source_member)
        } else {
            vec![source_solid_id]
        };
        for node_index in tool_offset..tool_offset + family.tool.model.nodes.len() {
            if let Some(provenance) = model.provenance.get_mut(&NodeId(node_index)) {
                provenance.source_entity_ids.clone_from(&exemplar_sources);
            }
        }

        let mut instance_roots = Vec::with_capacity(family.members.len());
        for member in &family.members {
            let quarter =
                (family.source_member.normalized_quarter + 4 - member.normalized_quarter) % 4;
            let translation = crate::planar_features::rigid_translation(
                family.source_member.center_mm,
                member.center_mm,
                quarter,
            );
            let root = if quarter == 0
                && translation.iter().all(|value| {
                    value.abs() <= DEFAULT_STEP_RECOVERY_TOLERANCES.source_identity_transform
                }) {
                tool_root
            } else {
                let transform = model.add_node(CadNode::Transform {
                    transform: quarter_turn_transform(quarter, translation),
                    child: tool_root,
                });
                model.set_provenance(
                    transform,
                    Provenance {
                        source_entity_ids: if family.original_source_ids {
                            feature_member_source_ids(source_solid_id, member)
                        } else {
                            vec![source_solid_id]
                        },
                        proof: ProofStatus::StructurallyProven,
                        max_residual_mm: None,
                    },
                )?;
                transform
            };
            instance_roots.push(root);
        }

        let op = match family.polarity {
            FeaturePolarity::Additive => BooleanOp::Union,
            FeaturePolarity::Subtractive => BooleanOp::Difference,
        };
        if batch
            .as_ref()
            .is_some_and(|(pass, previous, _, _)| *pass != family.pass_index || *previous != op)
        {
            current = finish_feature_batch(&mut model, current, batch.take())?;
        }
        let (_, _, children, sources) = batch.get_or_insert_with(|| {
            (
                family.pass_index,
                op,
                Vec::new(),
                BTreeSet::from([source_solid_id]),
            )
        });
        children.extend(instance_roots);
        if family.original_source_ids {
            for member in &family.members {
                sources.extend(feature_member_source_ids(source_solid_id, member));
            }
        }
    }
    current = finish_feature_batch(&mut model, current, batch)?;

    model.add_root(current)?;
    model.validate()?;
    let passes = families
        .iter()
        .map(|family| family.pass_index)
        .max()
        .map_or(0, |last| last + 1);
    Ok(CadFragment {
        source: CadFragmentSource::BoundaryFeatureComposite {
            solid_id: source_solid_id,
            passes,
            feature_families: families.len(),
            feature_instances: families.iter().map(|family| family.members.len()).sum(),
        },
        model,
        root: current,
    })
}

type FeatureReplayBatch = (usize, BooleanOp, Vec<NodeId>, BTreeSet<u64>);

fn finish_feature_batch(
    model: &mut CadModel,
    host: NodeId,
    batch: Option<FeatureReplayBatch>,
) -> Result<NodeId> {
    let Some((_, op, tools, sources)) = batch else {
        return Ok(host);
    };
    let children = std::iter::once(host).chain(tools).collect();
    let root = model.add_node(CadNode::Boolean { op, children });
    model.set_provenance(
        root,
        Provenance {
            source_entity_ids: sources.into_iter().collect(),
            proof: ProofStatus::StructurallyProven,
            max_residual_mm: None,
        },
    )?;
    Ok(root)
}

fn shifted_fragment_root(root: NodeId, offset: usize) -> Result<NodeId> {
    Ok(NodeId(root.0.checked_add(offset).ok_or_else(|| {
        anyhow::anyhow!("CAD root id overflow while composing fragments")
    })?))
}

fn feature_member_source_ids(
    source_solid_id: u64,
    member: &BoundaryFeatureInstanceEvidence,
) -> Vec<u64> {
    let mut out = Vec::with_capacity(member.face_ids.len() + member.interface_bound_ids.len() + 1);
    out.push(source_solid_id);
    out.extend(member.face_ids.iter().copied());
    out.extend(member.interface_bound_ids.iter().copied());
    out.sort_unstable();
    out.dedup();
    out
}

fn quarter_turn_transform(quarter: u8, translation: [f64; 3]) -> RigidTransform {
    let (xx, xy, yx, yy) = match quarter % 4 {
        0 => (1.0, 0.0, 0.0, 1.0),
        1 => (0.0, -1.0, 1.0, 0.0),
        2 => (-1.0, 0.0, 0.0, -1.0),
        3 => (0.0, 1.0, -1.0, 0.0),
        _ => unreachable!(),
    };
    RigidTransform {
        matrix: [
            [xx, xy, 0.0, translation[0]],
            [yx, yy, 0.0, translation[1]],
            [0.0, 0.0, 1.0, translation[2]],
            [0.0, 0.0, 0.0, 1.0],
        ],
    }
}

fn scan_source_solids(section: &RecoverySection<'_>) -> Vec<SourceSolidSummary> {
    let mut out = Vec::new();
    for entity in section.entities {
        let Some(record) = crate::step_graph::simple_record(entity) else {
            continue;
        };
        if record.name != "MANIFOLD_SOLID_BREP" {
            continue;
        }
        let solid_id = crate::step_graph::entity_id(entity);
        let Some(source) = section.source_entity(solid_id) else {
            continue;
        };
        let closure = crate::shape_identity::semantic_solid_closure(
            solid_id,
            section.entities,
            &section.index,
        )
        .unwrap_or_else(|| {
            crate::step_entities::closure_from(solid_id, section.entities, &section.index)
        });
        let mut faces = 0usize;
        let mut edges = 0usize;
        for id in &closure {
            let Some(&entity_index) = section.index.get(id) else {
                continue;
            };
            let Some(child) = crate::step_graph::simple_record(&section.entities[entity_index])
            else {
                continue;
            };
            match child.name.as_str() {
                "ADVANCED_FACE" => faces += 1,
                "EDGE_CURVE" => edges += 1,
                _ => {}
            }
        }
        out.push(SourceSolidSummary {
            source,
            faces,
            edges,
            closure_entities: closure.len(),
        });
    }
    out
}

const fn source_solid_id(source: &CadFragmentSource) -> u64 {
    match source {
        CadFragmentSource::SolidExtrusion { solid_id, .. }
        | CadFragmentSource::SolidLoft { solid_id, .. }
        | CadFragmentSource::SolidRevolution { solid_id, .. }
        | CadFragmentSource::RadialSlotRevolution { solid_id, .. }
        | CadFragmentSource::SolidSweep { solid_id, .. }
        | CadFragmentSource::RevolvedRoundTail { solid_id, .. }
        | CadFragmentSource::PeriodicChain { solid_id, .. }
        | CadFragmentSource::ExtrusionCoreUnion { solid_id, .. }
        | CadFragmentSource::SectionalExtrusionUnion { solid_id, .. }
        | CadFragmentSource::SectionalExtrusionTransitions { solid_id, .. }
        | CadFragmentSource::SectionalLoftUnion { solid_id, .. }
        | CadFragmentSource::FormedSheetExtrusionUnion { solid_id, .. }
        | CadFragmentSource::SourceCapExtrusionUnion { solid_id, .. }
        | CadFragmentSource::ExtrusionCoreRuledEnds { solid_id, .. }
        | CadFragmentSource::CapExtrusionRuledEnds { solid_id, .. }
        | CadFragmentSource::RuledEndExtension { solid_id, .. }
        | CadFragmentSource::BoundaryFeatureComposite { solid_id, .. }
        | CadFragmentSource::EdgeTreatmentComposite { solid_id, .. }
        | CadFragmentSource::SameDomainNormalized { solid_id } => *solid_id,
    }
}

const fn is_constructive_source(source: &CadFragmentSource) -> bool {
    !matches!(
        source,
        CadFragmentSource::PeriodicChain { .. } | CadFragmentSource::RuledEndExtension { .. }
    )
}

const fn whole_solid_candidate_is_preferred(
    existing_source: &CadFragmentSource,
    existing_complexity: u64,
    candidate_source: &CadFragmentSource,
    candidate_complexity: u64,
) -> bool {
    match (
        is_constructive_source(existing_source),
        is_constructive_source(candidate_source),
    ) {
        (true, false) => false,
        (false, true) => true,
        _ => candidate_complexity < existing_complexity,
    }
}

fn select_cheapest_whole_solid_fragments<'a>(
    candidates: &[&'a CadFragment],
) -> Result<BTreeMap<u64, &'a CadFragment>> {
    let mut selected = BTreeMap::<u64, (&CadFragment, u64)>::new();
    for &fragment in candidates {
        let solid_id = source_solid_id(&fragment.source);
        let complexity = fragment.model.complexity_score(fragment.root)?;
        let replace =
            selected
                .get(&solid_id)
                .is_none_or(|(existing_fragment, existing_complexity)| {
                    whole_solid_candidate_is_preferred(
                        &existing_fragment.source,
                        *existing_complexity,
                        &fragment.source,
                        complexity,
                    )
                });
        if replace {
            selected.insert(solid_id, (fragment, complexity));
        }
    }
    Ok(selected
        .into_iter()
        .map(|(solid_id, (fragment, _))| (solid_id, fragment))
        .collect())
}

fn selected_category<F>(
    selected: &BTreeMap<SourceEntity, &CadFragment>,
    predicate: F,
) -> Vec<SourceCadFragment>
where
    F: Fn(&CadFragmentSource) -> bool,
{
    selected
        .iter()
        .filter(|(_, fragment)| predicate(&fragment.source))
        .map(|(&source, fragment)| SourceCadFragment {
            source_solid: source,
            fragment: (*fragment).clone(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AnalysisSession;
    use crate::cad_ir::{CadModel, CadNode, Primitive};

    /// Compare the canonical-identity cost with the proof-preserving
    /// boundary provider on real STEP models. A diagnostic, not a timing gate.
    #[test]
    #[ignore = "manual source-identity prefilter benchmark"]
    fn boundary_identity_prefilter_benchmark() -> Result<()> {
        use std::hint::black_box;
        use std::time::Instant;

        for (label, bytes, runs) in [
            (
                "radial_slot",
                include_bytes!("../validation/fixtures/production_radial_slot.step").as_slice(),
                12_u32,
            ),
            (
                "split_ring",
                include_bytes!("../validation/fixtures/native_split_ring_hollow_cylinder.step")
                    .as_slice(),
                4_u32,
            ),
        ] {
            let exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
            let [data] = exchange.data.as_slice() else {
                anyhow::bail!("fixture must contain exactly one DATA section");
            };
            let entities = &data.entities;
            let index = crate::step_graph::build_index(entities);
            let section = RecoverySection {
                data_section: 0,
                entities,
                index: Cow::Borrowed(&index),
            };
            let targets = scan_source_solids(&section)
                .into_iter()
                .map(|summary| summary.source.id)
                .collect::<BTreeSet<_>>();
            let session = AnalysisSession::from_bytes(bytes)?;
            let measure = |task: &mut dyn FnMut() -> Result<()>| -> Result<f64> {
                let start = Instant::now();
                for _ in 0..runs {
                    task()?;
                }
                Ok(start.elapsed().as_secs_f64() * 1_000.0 / f64::from(runs))
            };
            let identity_ms = measure(&mut || {
                for &solid_id in &targets {
                    black_box(crate::shape_identity::solid_identity(
                        solid_id, entities, &index,
                    ));
                }
                Ok(())
            })?;
            let boundary_ms = measure(&mut || {
                black_box(recover_boundary_feature_composites_for_section(
                    &section,
                    &targets,
                    ReducedBrepMutationMode::ParserOwned,
                )?);
                Ok(())
            })?;
            let complete_ir_ms = measure(&mut || {
                black_box(session.complete_ir()?);
                Ok(())
            })?;
            let parsing_ms = measure(&mut || {
                black_box(AnalysisSession::from_bytes(bytes)?);
                Ok(())
            })?;
            println!(
                "boundary identity prefilter: {label}, entities={}, targets={}, runs={runs}, canonical_identity_ms={identity_ms:.3}, boundary_provider_ms={boundary_ms:.3}, parsed_complete_ir_ms={complete_ir_ms:.3}, parse_ms={parsing_ms:.3}",
                entities.len(),
                targets.len(),
            );
        }
        Ok(())
    }

    #[test]
    fn singleton_boundary_recovery_preserves_nested_peels_and_data_section_identity() -> Result<()>
    {
        let bytes = include_bytes!("../validation/fixtures/boundary_boss_nested.step");
        let original = recover_complete_ir_bytes(bytes)?;
        let [expected] = original.boundary_feature_composites.as_slice() else {
            anyhow::bail!("nested-boss fixture must recover exactly one boundary composite");
        };
        assert!(matches!(
            expected.fragment.source,
            CadFragmentSource::BoundaryFeatureComposite {
                passes: 2,
                feature_families: 2,
                feature_instances: 2,
                ..
            }
        ));
        let mut exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let [source] = exchange.data.as_slice() else {
            anyhow::bail!("nested-boss fixture must contain one DATA section");
        };
        exchange.data.push(source.clone());
        let doubled =
            recover_complete_ir_build_state(&exchange, ReducedBrepMutationMode::ParserOwned, None)?
                .finish(COMPLETE_IR_SCHEMA, Vec::new(), None);
        assert_eq!(doubled.body_count, 2);
        assert_eq!(doubled.boundary_feature_composites.len(), 2);
        for (data_section, composite) in doubled.boundary_feature_composites.iter().enumerate() {
            assert_eq!(composite.source_solid.data_section, data_section);
            assert_eq!(composite.fragment, expected.fragment);
        }
        assert_eq!(doubled.unrecovered_solid_count, 0);
        Ok(())
    }

    /// Clone an independent source graph into the same DATA section while
    /// preserving its internal references. Used to exercise real multi-solid
    /// grouping, not the independent-index behavior of two DATA sections.
    fn append_disjoint_source_entities(
        into: &mut Vec<EntityInstance>,
        source: &[EntityInstance],
    ) -> Result<HashMap<u64, u64>> {
        let offset = into
            .iter()
            .map(crate::step_graph::entity_id)
            .max()
            .unwrap_or(0)
            .checked_add(1)
            .ok_or_else(|| anyhow::anyhow!("STEP test entity ID overflow"))?;
        let aliases = source
            .iter()
            .map(|entity| {
                let id = crate::step_graph::entity_id(entity);
                id.checked_add(offset)
                    .map(|replacement| (id, replacement))
                    .ok_or_else(|| anyhow::anyhow!("STEP test entity ID overflow"))
            })
            .collect::<Result<HashMap<_, _>>>()?;
        let mut cloned = source.to_vec();
        for entity in &mut cloned {
            crate::step_graph::rewrite_entity_refs(entity, &aliases);
            let id = match entity {
                EntityInstance::Simple { id, .. } | EntityInstance::Complex { id, .. } => id,
            };
            *id = aliases[id];
        }
        into.extend(cloned);
        Ok(aliases)
    }

    #[test]
    fn boundary_identity_prefilter_preserves_multi_solid_equivalence_and_distinct_shapes()
    -> Result<()> {
        let boss_bytes = include_bytes!("../validation/fixtures/boundary_boss_nested.step");
        let conical_bytes = include_bytes!("../validation/fixtures/native_conical_frustum.step");
        let mut boss = monstertruck_io::step::load::ParsedStepSource::parse(boss_bytes)?.exchange;
        let conical = monstertruck_io::step::load::ParsedStepSource::parse(conical_bytes)?.exchange;
        let [boss_data] = boss.data.as_mut_slice() else {
            anyhow::bail!("expected exactly one boss DATA section");
        };
        let original = boss_data.entities.clone();
        let boss_id = original
            .iter()
            .filter(|entity| {
                crate::step_graph::simple_record(entity)
                    .is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
            })
            .map(crate::step_graph::entity_id)
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing boss source solid"))?;
        let boss_aliases = append_disjoint_source_entities(&mut boss_data.entities, &original)?;
        let [conical_data] = conical.data.as_slice() else {
            anyhow::bail!("expected exactly one conical DATA section");
        };
        let conical_aliases =
            append_disjoint_source_entities(&mut boss_data.entities, &conical_data.entities)?;
        let conical_id = conical_data
            .entities
            .iter()
            .filter(|entity| {
                crate::step_graph::simple_record(entity)
                    .is_some_and(|record| record.name == "MANIFOLD_SOLID_BREP")
            })
            .map(crate::step_graph::entity_id)
            .next()
            .ok_or_else(|| anyhow::anyhow!("missing conical source solid"))?;
        let index = crate::step_graph::build_index(&boss_data.entities);
        let boss_faces =
            crate::shape_identity::manifold_solid_face_ids(boss_id, &boss_data.entities, &index)
                .ok_or_else(|| anyhow::anyhow!("missing boss faces"))?;
        let conical_faces = crate::shape_identity::manifold_solid_face_ids(
            conical_aliases[&conical_id],
            &boss_data.entities,
            &index,
        )
        .ok_or_else(|| anyhow::anyhow!("missing conical faces"))?;
        assert_ne!(boss_faces.len(), conical_faces.len());

        let section = RecoverySection {
            data_section: 0,
            entities: &boss_data.entities,
            index: Cow::Borrowed(&index),
        };
        let targets = BTreeSet::from([
            boss_id,
            boss_aliases[&boss_id],
            conical_aliases[&conical_id],
        ]);
        let fragments = recover_boundary_feature_composites_for_section(
            &section,
            &targets,
            ReducedBrepMutationMode::ParserOwned,
        )?;
        let recovered_boss_ids = fragments
            .iter()
            .filter(|fragment| {
                matches!(
                    fragment.source,
                    CadFragmentSource::BoundaryFeatureComposite {
                        passes: 2,
                        feature_families: 2,
                        feature_instances: 2,
                        ..
                    }
                )
            })
            .map(|fragment| source_solid_id(&fragment.source))
            .collect::<BTreeSet<_>>();
        assert_eq!(
            recovered_boss_ids,
            BTreeSet::from([boss_id, boss_aliases[&boss_id]])
        );
        for fragment in &fragments {
            fragment.model.validate()?;
        }
        Ok(())
    }

    #[test]
    fn cached_session_complete_ir_matches_independent_section_indexes() -> Result<()> {
        for bytes in [
            include_bytes!("../validation/fixtures/native_conical_frustum.step").as_slice(),
            include_bytes!("../validation/fixtures/production_radial_slot.step").as_slice(),
        ] {
            let session = AnalysisSession::from_bytes(bytes)?;
            let cached = session.complete_ir()?;
            let exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
            let rebuilt = recover_complete_ir_build_state(
                &exchange,
                ReducedBrepMutationMode::ParserOwned,
                None,
            )?
            .finish(COMPLETE_IR_SCHEMA, Vec::new(), None);
            assert_eq!(
                serde_json::to_value(&cached)?,
                serde_json::to_value(&rebuilt)?
            );
        }
        Ok(())
    }

    #[test]
    fn owned_analysis_session_outlives_its_step_input_and_omits_rewrite_diagnostics() -> Result<()>
    {
        let session = {
            let input =
                include_bytes!("../validation/fixtures/production_radial_slot.step").to_vec();
            AnalysisSession::from_bytes(&input)?
        };
        let report = session.complete_ir()?;
        assert_eq!(report.schema, COMPLETE_IR_SCHEMA);
        assert_eq!(report.body_count, 1);
        let serialized = serde_json::to_value(&report)?;
        assert!(serialized.get("idempotent_compact_rewrite").is_none());
        assert!(serialized.get("compact_second_pass_stats").is_none());
        Ok(())
    }

    #[test]
    fn explicitly_scoped_tool_recognition_matches_full_section_candidates() -> Result<()> {
        // Include both authored and newly materialized bodies, preserving the
        // original representation contexts and every geometry reference.
        let bytes = include_bytes!("../validation/fixtures/native_relief_blind_636f2bd5.step");
        let mut parsed = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let entities = &mut parsed.data[0].entities;
        let original_index = crate::step_graph::build_index(entities);
        let original_ids = crate::brep::manifold_solid_ids(entities);
        let report = crate::planar_features::diagnose_boundary_features(entities);
        let tools = crate::planar_features::materialize_boundary_feature_tools_from_report(
            entities,
            &original_index,
            &report,
            1,
        )?;
        let tool_ids = tools.iter().map(|tool| tool.solid_id).collect::<Vec<_>>();
        assert!(!tool_ids.is_empty());
        assert!(tool_ids.iter().all(|id| !original_ids.contains(id)));
        let index = crate::step_graph::build_index(entities);
        let tolerances =
            crate::solid_revolutions::source_tolerance_by_representation_item(entities, &index);
        let full = recover_constructive_fragments_with_tolerances(
            entities,
            &index,
            &tolerances,
            &crate::brep::manifold_solid_ids(entities),
        )?;
        let scoped = recover_constructive_fragments_with_tolerances(
            entities,
            &index,
            &tolerances,
            &tool_ids,
        )?;
        let expected = full
            .into_iter()
            .filter(|fragment| tool_ids.contains(&source_solid_id(&fragment.source)))
            .collect::<Vec<_>>();
        assert!(!expected.is_empty());
        assert_eq!(scoped, expected);
        assert_eq!(
            recover_constructive_fragments_with_tolerances(entities, &index, &tolerances, &[])?,
            []
        );
        Ok(())
    }

    const fn test_source_in(data_section: usize, solid_id: u64) -> SourceEntity {
        SourceEntity {
            data_section,
            ordinal: 0,
            id: solid_id,
        }
    }

    const fn test_source(solid_id: u64) -> SourceEntity { test_source_in(0, solid_id) }

    fn qualified_fragment_in(data_section: usize, fragment: CadFragment) -> SourceCadFragment {
        SourceCadFragment {
            source_solid: test_source_in(data_section, source_solid_id(&fragment.source)),
            fragment,
        }
    }

    fn qualified_fragment(fragment: CadFragment) -> SourceCadFragment {
        qualified_fragment_in(0, fragment)
    }

    fn report_with_unrecovered(solid_ids: &[u64]) -> CompleteIr {
        let solids = solid_ids
            .iter()
            .map(|&solid_id| SourceSolidSummary {
                source: test_source(solid_id),
                faces: 0,
                edges: 0,
                closure_entities: 0,
            })
            .collect::<Vec<_>>();
        CompleteIr {
            schema: COMPLETE_IR_SCHEMA,
            surface_signature_summaries: Vec::new(),
            instance_pattern_diagnostics: Vec::new(),
            solid_extrusions: Vec::new(),
            solid_lofts: Vec::new(),
            solid_revolutions: Vec::new(),
            radial_slot_revolutions: Vec::new(),
            solid_sweeps: Vec::new(),
            periodic_chains: Vec::new(),
            boundary_feature_composites: Vec::new(),
            edge_treatment_composites: Vec::new(),
            same_domain_normalized_solids: Vec::new(),
            extrusion_core_unions: Vec::new(),
            kernel_recovery_failures: Vec::new(),
            kernel_source_load: None,
            unrecovered_solids: solids.clone(),
            body_scan: BodyScanSummary { solids },
            body_count: solid_ids.len(),
            source_boundary_verified_solids: Vec::new(),
            source_orientation_reconciled_solids: Vec::new(),
            constructively_recovered_solid_count: 0,
            semantically_recovered_solid_count: 0,
            unrecovered_solid_count: solid_ids.len(),
        }
    }

    fn build_state_with_unrecovered(solid_ids: &[u64]) -> CompleteIrBuildState {
        let report = report_with_unrecovered(solid_ids);
        let body_scan = report.body_scan.solids;
        let recovery = WholeSolidRecoveryState::new(&body_scan);
        CompleteIrBuildState {
            surface_signature_summaries: report.surface_signature_summaries,
            instance_pattern_diagnostics: report.instance_pattern_diagnostics,
            body_scan,
            recovery,
            direct_candidates: BTreeMap::new(),
        }
    }

    fn apply_test_kernel_recovery<F>(solid_ids: &[u64], mut recover: F) -> CompleteIr
    where F: FnMut(u64) -> Result<Option<CadFragment>> {
        let mut build = build_state_with_unrecovered(solid_ids);
        let failures = apply_kernel_recovery(
            &mut build.recovery,
            |source, _| {
                recover(source.id)
                    .map(|candidate| candidate.map(KernelRecoveryCandidate::structural))
            },
            false,
        );
        build.finish(KERNEL_COMPLETE_IR_SCHEMA, failures, None)
    }

    fn constructive_fragment(solid_id: u64, boxes: usize) -> Result<CadFragment> {
        if boxes == 0 {
            bail!("test fragment requires at least one box");
        }
        let mut model = CadModel::new();
        let children = (0..boxes)
            .map(|_| {
                model.add_node(CadNode::Primitive(Primitive::Box {
                    size_mm: [1.0, 1.0, 1.0],
                }))
            })
            .collect::<Vec<_>>();
        let root = if children.len() == 1 {
            children[0]
        } else {
            model.add_node(CadNode::Assembly { children })
        };
        model.add_root(root)?;
        Ok(CadFragment {
            source: CadFragmentSource::SolidExtrusion {
                solid_id,
                cap_face_ids: [1, 2],
                side_face_ids: Vec::new(),
            },
            model,
            root,
        })
    }

    struct SummaryOnlyKernel {
        summary: crate::cad_kernel::KernelSummary,
    }

    impl CadKernel for SummaryOnlyKernel {
        type Evaluated = crate::cad_kernel::KernelSummary;

        fn evaluate(&self, _: &CadModel, _: NodeId) -> Result<Self::Evaluated> { Ok(self.summary) }

        fn summarize(&self, evaluated: &Self::Evaluated) -> crate::cad_kernel::KernelSummary {
            *evaluated
        }

        fn supports_step_export(&self) -> bool { false }

        fn to_step(&self, _: &Self::Evaluated) -> Result<String> {
            bail!("kernel replay must not export and reparse a previously evaluated solid")
        }
    }

    #[test]
    fn kernel_replay_admits_single_native_solid_without_step_roundtrip() -> Result<()> {
        use crate::cad_kernel::KernelSummary;
        let fragment = constructive_fragment(12, 1)?;
        for (summary, expected) in [
            (
                KernelSummary {
                    geometrically_consistent: true,
                    shells: 1,
                    faces: 6,
                },
                true,
            ),
            (
                KernelSummary {
                    geometrically_consistent: true,
                    shells: 0,
                    faces: 0,
                },
                false,
            ),
            (
                KernelSummary {
                    geometrically_consistent: true,
                    shells: 1,
                    faces: 0,
                },
                false,
            ),
            (
                KernelSummary {
                    geometrically_consistent: true,
                    shells: 2,
                    faces: 12,
                },
                false,
            ),
            (
                KernelSummary {
                    geometrically_consistent: false,
                    shells: 1,
                    faces: 6,
                },
                false,
            ),
        ] {
            let kernel = SummaryOnlyKernel { summary };
            assert_eq!(kernel_replay_is_single_solid(&kernel, &fragment)?, expected);
        }
        Ok(())
    }

    fn test_treatment_evidence(
        solid_id: u64,
        treatment_face_id: u64,
    ) -> crate::edge_treatments::EdgeTreatmentEvidence {
        use crate::edge_treatments::{
            EdgeTreatmentCarrier, EdgeTreatmentEvidence, EdgeTreatmentEvidenceStrength,
            EdgeTreatmentKind, EdgeTreatmentVirtualEdge,
        };
        EdgeTreatmentEvidence {
            solid_id,
            treatment_face_id,
            kind: EdgeTreatmentKind::PlanarChamfer,
            strength: EdgeTreatmentEvidenceStrength::AnalyticTangency,
            carrier_pair_unique: true,
            treatment_support: "plane".to_owned(),
            carriers: [
                EdgeTreatmentCarrier {
                    face_id: treatment_face_id + 1,
                    support: "plane".to_owned(),
                    shared_edge_ids: vec![treatment_face_id + 10],
                },
                EdgeTreatmentCarrier {
                    face_id: treatment_face_id + 2,
                    support: "plane".to_owned(),
                    shared_edge_ids: vec![treatment_face_id + 11],
                },
            ],
            inferred_radius_mm: None,
            inferred_setbacks_mm: Some([0.5, 0.5]),
            first_setback_carrier_normal: Some([1.0, 0.0, 0.0]),
            inferred_angle_rad: Some(std::f64::consts::FRAC_PI_4),
            carrier_angle_rad: std::f64::consts::FRAC_PI_2,
            virtual_edge: EdgeTreatmentVirtualEdge::Segment {
                start_mm: [0.0, 0.0, 0.0],
                end_mm: [0.0, 0.0, 1.0],
            },
            max_boundary_residual_mm: 1.0e-8,
        }
    }

    fn reciprocal_chamfer_candidate_ring(
        solid_id: u64,
        count: u32,
    ) -> Vec<crate::edge_treatments::EdgeTreatmentEvidence> {
        use crate::edge_treatments::EdgeTreatmentVirtualEdge;
        (0..count)
            .map(|i| {
                let mut item = test_treatment_evidence(solid_id, u64::from(i));
                let previous = (i + count - 1) % count;
                let next = (i + 1) % count;
                item.carriers[0].face_id = u64::from(previous);
                item.carriers[0].shared_edge_ids = vec![100 + u64::from(previous)];
                item.carriers[1].face_id = u64::from(next);
                item.carriers[1].shared_edge_ids = vec![100 + u64::from(i)];
                let angle = std::f64::consts::TAU * f64::from(i) / f64::from(count);
                let (sine, cosine) = angle.sin_cos();
                item.virtual_edge = EdgeTreatmentVirtualEdge::Segment {
                    start_mm: [cosine, sine, 0.0],
                    end_mm: [cosine, sine, 1.0],
                };
                item
            })
            .collect()
    }

    #[test]
    fn reciprocal_planar_profile_ring_retains_all_host_faces_not_arbitrary_chamfers() {
        let mut candidates = reciprocal_chamfer_candidate_ring(42, 8);
        let ordinary = test_treatment_evidence(42, 99);
        candidates.push(ordinary.clone());
        // Reused raw entity numbers in a different solid do not connect rings.
        let other_solid = test_treatment_evidence(43, 2);
        candidates.push(other_solid.clone());
        assert_eq!(
            select_replayable_treatments(candidates),
            [ordinary, other_solid]
        );
    }

    #[test]
    fn reciprocal_planar_profile_ring_does_not_discard_nonuniform_or_open_chamfer_sets() {
        let ring = reciprocal_chamfer_candidate_ring(42, 8);
        let mut unequal = ring.clone();
        unequal[0].inferred_setbacks_mm = Some([0.2, 0.2]);
        assert_eq!(select_replayable_treatments(unequal.clone()), unequal);
        let mut open = ring.clone();
        open.pop();
        assert_eq!(select_replayable_treatments(open.clone()), open);
        let mut unrelated_edge = ring.clone();
        unrelated_edge[0].carriers[0].shared_edge_ids = vec![999];
        assert_eq!(
            select_replayable_treatments(unrelated_edge.clone()),
            unrelated_edge
        );
        let mut reverse = ring;
        if let crate::edge_treatments::EdgeTreatmentVirtualEdge::Segment { start_mm, end_mm } =
            &mut reverse[0].virtual_edge
        {
            std::mem::swap(start_mm, end_mm);
        }
        reverse.reverse();
        assert!(select_replayable_treatments(reverse).is_empty());
    }

    #[test]
    fn supplier_chamfer_micron_setback_drift_preserves_both_dimensions() -> Result<()> {
        let mut evidence = test_treatment_evidence(42, 77);
        evidence.inferred_setbacks_mm = Some([0.5, 0.506_7]);
        let (profile, edge) = replayable_edge_treatment(&evidence)
            .ok_or_else(|| anyhow::anyhow!("6.7 micrometre difference must replay"))?;
        assert_eq!(
            profile,
            EdgeTreatmentProfile::CarrierChamfer {
                first_carrier_normal: [1.0, 0.0, 0.0],
                second_setback_mm: 0.506_7,
            }
        );
        assert_eq!(edge.radius_mm, 0.5);
        evidence.inferred_setbacks_mm = Some([0.5, 0.500_000_01]);
        let (profile, edge) = replayable_edge_treatment(&evidence)
            .context("sub-tolerance setback asymmetry must be preserved")?;
        assert_eq!(
            profile,
            EdgeTreatmentProfile::CarrierChamfer {
                first_carrier_normal: [1.0, 0.0, 0.0],
                second_setback_mm: 0.500_000_01,
            }
        );
        assert_eq!(edge.radius_mm, 0.5);
        evidence.inferred_setbacks_mm = Some([0.5, 0.510_1]);
        assert!(replayable_edge_treatment(&evidence).is_none());
        evidence.inferred_setbacks_mm = Some([0.5, f64::NAN]);
        assert!(replayable_edge_treatment(&evidence).is_none());
        Ok(())
    }

    #[test]
    fn replayed_edge_treatments_are_first_class_composites() -> Result<()> {
        let solid_id = 42;
        let evidence = test_treatment_evidence(solid_id, 77);
        let fragment = replay_edge_treatments(
            constructive_fragment(solid_id, 1)?,
            solid_id,
            vec![(
                &evidence,
                EdgeTreatmentProfile::Chamfer,
                EdgeTreatmentEdge {
                    selector: EdgeSelector::Segment {
                        start_mm: [0.0, 0.0, 0.0],
                        end_mm: [0.0, 0.0, 1.0],
                    },
                    radius_mm: 0.5,
                },
            )],
        )?;
        assert_eq!(
            fragment.source,
            CadFragmentSource::EdgeTreatmentComposite {
                solid_id,
                treatment_face_ids: vec![77],
            }
        );
        assert!(matches!(
            fragment.model.node(fragment.root)?,
            CadNode::EdgeTreatment { .. }
        ));
        assert_eq!(
            fragment.model.provenance[&fragment.root].source_entity_ids,
            vec![solid_id, 77]
        );
        Ok(())
    }

    #[test]
    fn coupled_replay_retains_bounded_arcs_and_an_atomic_treatment_node() -> Result<()> {
        let mut first = test_treatment_evidence(42, 77);
        first.kind = crate::edge_treatments::EdgeTreatmentKind::ToroidalFillet;
        first.inferred_radius_mm = Some(0.2);
        first.virtual_edge = crate::edge_treatments::EdgeTreatmentVirtualEdge::Circle {
            center_mm: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            radius_mm: 2.0,
        };
        let mut second = test_treatment_evidence(42, 78);
        second.kind = crate::edge_treatments::EdgeTreatmentKind::CylindricalFillet;
        second.inferred_radius_mm = Some(0.2);
        let replay = [&first, &second]
            .into_iter()
            .map(|e| {
                let (profile, edge) =
                    replayable_edge_treatment(e).context("missing replayable treatment")?;
                Ok((e, profile, edge))
            })
            .collect::<Result<Vec<_>>>()?;
        let arc = EdgeSelector::CircleArc {
            center_mm: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            radius_mm: 2.0,
            start_mm: [2.0, 0.0, 0.0],
            sweep_rad: std::f64::consts::FRAC_PI_2,
        };
        let segment = EdgeSelector::Segment {
            start_mm: [0.0, 2.0, 0.0],
            end_mm: [-1.0, 2.0, 0.0],
        };
        let fragment = replay_grouped_edge_treatments(
            constructive_fragment(42, 1)?,
            42,
            replay,
            vec![crate::kernel_topology::CoupledTreatmentReplay {
                selectors: vec![(78, segment.clone()), (77, arc.clone())],
                junction_faces: Vec::new(),
            }],
        )?;
        let CadNode::EdgeTreatment { edges, profile, .. } = fragment.model.node(fragment.root)?
        else {
            bail!("missing atomic treatment node");
        };
        assert_eq!(*profile, EdgeTreatmentProfile::Round);
        assert_eq!(edges.len(), 2);
        assert_eq!(edges[0].selector, arc);
        assert_eq!(edges[1].selector, segment);
        assert_eq!(
            fragment.model.provenance[&fragment.root].source_entity_ids,
            vec![42, 77, 78]
        );
        Ok(())
    }

    fn boundary_reduction_step(solid_id: u64) -> Result<(u64, LocalReductionStep)> {
        let member = BoundaryFeatureInstanceEvidence {
            face_ids: vec![10],
            center_mm: [0.0; 3],
            normalized_quarter: 0,
            interface_bound_ids: vec![11],
        };
        Ok((
            solid_id,
            LocalReductionStep::BoundaryFeature(RecoveredBoundaryFeatureFamily {
                pass_index: 0,
                polarity: FeaturePolarity::Subtractive,
                tool: constructive_fragment(99, 1)?,
                source_member: member.clone(),
                members: vec![member],
                original_source_ids: true,
            }),
        ))
    }

    fn descriptive_fragment(solid_id: u64, boxes: usize) -> Result<CadFragment> {
        if boxes == 0 {
            bail!("test fragment requires at least one box");
        }
        let mut model = CadModel::new();
        let children = (0..boxes)
            .map(|_| {
                model.add_node(CadNode::Primitive(Primitive::Box {
                    size_mm: [1.0, 1.0, 1.0],
                }))
            })
            .collect::<Vec<_>>();
        let root = if children.len() == 1 {
            children[0]
        } else {
            model.add_node(CadNode::Assembly { children })
        };
        model.add_root(root)?;
        Ok(CadFragment {
            source: CadFragmentSource::PeriodicChain { solid_id, sites: 2 },
            model,
            root,
        })
    }

    fn extrusion_core_union_fragment(solid_id: u64) -> Result<CadFragment> {
        let mut model = CadModel::new();
        let root = model.add_node(CadNode::Primitive(Primitive::Box {
            size_mm: [1.0, 1.0, 1.0],
        }));
        model.add_root(root)?;
        Ok(CadFragment {
            source: CadFragmentSource::ExtrusionCoreUnion {
                solid_id,
                cap_face_ids: vec![[1, 2], [3, 4]],
            },
            model,
            root,
        })
    }

    #[test]
    fn parser_only_boss_still_recovers_the_constructive_feature_composite() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/boundary_boss_extrusion.step");
        let report = recover_complete_ir_bytes(input)?;
        assert_eq!(report.schema, COMPLETE_IR_SCHEMA);
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.semantically_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        let [composite] = report.boundary_feature_composites.as_slice() else {
            bail!("parser-only recovery lost its proven boss feature");
        };
        assert!(matches!(
            composite.fragment.source,
            CadFragmentSource::BoundaryFeatureComposite {
                passes: 1,
                feature_families: 1,
                feature_instances: 1,
                ..
            }
        ));
        composite.fragment.model.validate()?;
        Ok(())
    }

    #[test]
    fn parser_only_nested_boss_preserves_two_proven_fixed_point_passes() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/boundary_boss_nested.step");
        let report = recover_complete_ir_bytes(input)?;
        assert_eq!(report.schema, COMPLETE_IR_SCHEMA);
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        let [composite] = report.boundary_feature_composites.as_slice() else {
            bail!("nested boss lost its parser-only boundary composite");
        };
        assert!(matches!(
            composite.fragment.source,
            CadFragmentSource::BoundaryFeatureComposite {
                passes: 2,
                feature_families: 2,
                feature_instances: 2,
                ..
            }
        ));
        composite.fragment.model.validate()?;
        Ok(())
    }

    #[test]
    fn parser_peel_journal_is_atomic_and_keeps_original_step_immutable() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/boundary_boss_extrusion.step");
        let exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let entities = &exchange.data[0].entities;
        let original = entities.clone();
        let index = crate::step_graph::build_index(entities);
        let report = crate::planar_features::diagnose_boundary_features(entities);
        let selected = crate::planar_features::leaf_feature_family_indices(&report.families);
        assert_eq!(selected.len(), 1);
        let plan = crate::planar_features::plan_boundary_feature_family_peel(
            entities, &index, &report, &selected,
        )?;
        assert_eq!(plan.stats.families, 1);
        assert_eq!(plan.stats.faces_removed_from_shells, 5);
        assert_eq!(plan.stats.interface_bounds_healed, 1);

        let mut journal = ParserPeelJournal::default();
        assert!(matches!(journal.view(entities, &index)?, Cow::Borrowed(_)));
        assert!(!journal.commit(Vec::new(), plan.clone(), entities, &index)?);

        let mut invalid = plan.clone();
        invalid.host_bound_removals[0].bound_ids.push(u64::MAX);
        assert!(
            journal
                .commit(vec![boundary_reduction_step(7)?], invalid, entities, &index)
                .is_err()
        );
        assert_eq!(journal.plans.len(), 0);
        assert_eq!(journal.replay.len(), 0);

        assert!(journal.commit(
            vec![boundary_reduction_step(7)?],
            plan.clone(),
            entities,
            &index
        )?);
        assert_eq!(journal.boundary_features(7).len(), 1);
        assert_eq!(journal.next_boundary_pass_index(), 1);
        let projected = journal.view(entities, &index)?;
        assert!(matches!(projected, Cow::Owned(_)));
        assert_eq!(entities, &original);
        let shell = plan.shell_face_removals[0].shell_id;
        let host = plan.host_bound_removals[0].host_face_id;
        let refs = crate::brep::ref_list_param;
        let original_shell = refs(shell, 1, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("source shell unreadable"))?;
        let projected_shell = refs(shell, 1, projected.as_ref(), &index)
            .ok_or_else(|| anyhow::anyhow!("projected shell unreadable"))?;
        let original_host = refs(host, 1, entities, &index)
            .ok_or_else(|| anyhow::anyhow!("source host unreadable"))?;
        let projected_host = refs(host, 1, projected.as_ref(), &index)
            .ok_or_else(|| anyhow::anyhow!("projected host unreadable"))?;
        assert_eq!(original_shell.len(), projected_shell.len() + 5);
        assert_eq!(original_host.len(), projected_host.len() + 1);
        assert!(
            journal
                .commit(
                    vec![boundary_reduction_step(7)?],
                    plan,
                    projected.as_ref(),
                    &index
                )
                .is_err()
        );
        assert_eq!(journal.plans.len(), 1);
        assert_eq!(journal.boundary_features(7).len(), 1);
        Ok(())
    }

    #[test]
    fn recovery_state_unifies_selection_and_unresolved_tracking() -> Result<()> {
        let report = report_with_unrecovered(&[7, 8]);
        let mut state = WholeSolidRecoveryState::new(&report.body_scan.solids);
        assert_eq!(
            state.constructively_unresolved_sources(),
            BTreeSet::from([test_source(7), test_source(8)])
        );

        assert_eq!(
            state.add_candidates([qualified_fragment(descriptive_fragment(7, 1)?)])?,
            1
        );
        assert_eq!(
            state.constructively_unresolved_sources(),
            BTreeSet::from([test_source(7), test_source(8)])
        );
        assert_eq!(
            state.semantically_unrecovered_sources(),
            BTreeSet::from([test_source(8)])
        );

        assert_eq!(
            state.add_candidates([qualified_fragment(constructive_fragment(7, 4)?)])?,
            1
        );
        assert_eq!(
            state.constructively_unresolved_sources(),
            BTreeSet::from([test_source(8)])
        );
        let revision = state.revision;

        assert_eq!(
            state.add_candidates([qualified_fragment(constructive_fragment(7, 6)?)])?,
            0
        );
        assert_eq!(state.revision, revision);
        assert_eq!(
            state.add_candidates([qualified_fragment(constructive_fragment(7, 1)?)])?,
            1
        );
        assert!(state.revision > revision);
        assert_eq!(
            state.selected_fragments()[&test_source(7)]
                .model
                .complexity_score(state.selected_fragments()[&test_source(7)].root)?,
            1
        );
        Ok(())
    }

    #[test]
    fn source_boundary_verified_candidate_outranks_cheaper_structural_guess() -> Result<()> {
        let report = report_with_unrecovered(&[7]);
        let mut state = WholeSolidRecoveryState::new(&report.body_scan.solids);
        let verified = constructive_fragment(7, 5)?;
        let expected_complexity = verified.model.complexity_score(verified.root)?;
        state.add_candidates([qualified_fragment(constructive_fragment(7, 1)?)])?;
        let failures = apply_kernel_recovery(
            &mut state,
            |_, ceiling| {
                assert_eq!(
                    ceiling, None,
                    "structural guess cannot prune certified evidence"
                );
                Ok(Some(KernelRecoveryCandidate {
                    fragment: verified.clone(),
                    boundary_verification: BoundaryVerification::Original,
                }))
            },
            true,
        );
        assert!(failures.is_empty(), "{failures:?}");
        assert!(matches!(
            state.selected[&test_source(7)].boundary_verification,
            BoundaryVerification::Original
        ));
        assert_eq!(
            state.selected[&test_source(7)].complexity,
            expected_complexity
        );
        let failures = apply_kernel_recovery(
            &mut state,
            |_, ceiling| {
                assert_eq!(ceiling, Some(expected_complexity));
                Ok(None)
            },
            true,
        );
        assert!(failures.is_empty(), "{failures:?}");
        // A later unverified shortcut cannot displace an exact proven body.
        state.add_candidates([qualified_fragment(constructive_fragment(7, 1)?)])?;
        assert_eq!(
            state.selected[&test_source(7)].complexity,
            expected_complexity
        );
        Ok(())
    }

    #[test]
    fn recovery_state_keeps_same_local_id_from_distinct_data_sections() -> Result<()> {
        let body_scan = vec![
            SourceSolidSummary {
                source: test_source_in(0, 7),
                faces: 0,
                edges: 0,
                closure_entities: 0,
            },
            SourceSolidSummary {
                source: test_source_in(1, 7),
                faces: 0,
                edges: 0,
                closure_entities: 0,
            },
        ];
        let mut state = WholeSolidRecoveryState::new(&body_scan);
        assert_eq!(
            state.add_candidates([
                qualified_fragment_in(0, constructive_fragment(7, 1)?),
                qualified_fragment_in(1, constructive_fragment(7, 1)?),
            ])?,
            2
        );
        assert_eq!(state.selected_fragments().len(), 2);
        assert!(state.semantically_unrecovered_sources().is_empty());
        assert!(
            state
                .selected_fragments()
                .contains_key(&test_source_in(0, 7))
        );
        assert!(
            state
                .selected_fragments()
                .contains_key(&test_source_in(1, 7))
        );
        Ok(())
    }

    #[test]
    fn recovery_state_rejects_candidate_for_unknown_solid() -> Result<()> {
        let report = report_with_unrecovered(&[7]);
        let mut state = WholeSolidRecoveryState::new(&report.body_scan.solids);
        let result = state.add_candidates([qualified_fragment(constructive_fragment(8, 1)?)]);
        let Err(error) = result else {
            bail!("unknown source solid unexpectedly succeeded");
        };
        assert!(error.to_string().contains("missing source solid"));
        Ok(())
    }

    const DUPLICATE_SOLID_ID_MULTI_DATA: &[u8] = br"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('duplicate-id regression'),'1');
FILE_NAME('duplicate.step','2026-11-02T00:00:00',(''),(''),'step-redox','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=MANIFOLD_SOLID_BREP('',#2);
#2=CLOSED_SHELL('',());
ENDSEC;
DATA;
#1=MANIFOLD_SOLID_BREP('',#2);
#2=CLOSED_SHELL('',());
ENDSEC;
END-ISO-10303-21;
";

    #[test]
    fn duplicate_solid_ids_in_distinct_data_sections_remain_distinct_sources() -> Result<()> {
        let report = recover_complete_ir_bytes(DUPLICATE_SOLID_ID_MULTI_DATA)?;
        assert_eq!(report.schema, COMPLETE_IR_SCHEMA);
        assert_eq!(report.body_count, 2);
        assert_eq!(report.unrecovered_solid_count, 2);
        let sources = report
            .body_scan
            .solids
            .iter()
            .map(|solid| solid.source)
            .collect::<Vec<_>>();
        assert_eq!(
            sources,
            vec![
                SourceEntity {
                    data_section: 0,
                    ordinal: 0,
                    id: 1,
                },
                SourceEntity {
                    data_section: 1,
                    ordinal: 0,
                    id: 1,
                },
            ]
        );
        assert_eq!(
            report
                .unrecovered_solids
                .iter()
                .map(|solid| solid.source)
                .collect::<Vec<_>>(),
            sources
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_complete_ir_handles_duplicate_solid_ids_across_data_sections() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(DUPLICATE_SOLID_ID_MULTI_DATA, &kernel)?;
        assert_eq!(report.schema, KERNEL_COMPLETE_IR_SCHEMA);
        assert_eq!(report.body_count, 2);
        assert_eq!(report.unrecovered_solid_count, 2);
        let evidence = report
            .kernel_source_load
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("kernel v14 omitted source-load evidence"))?;
        assert_eq!(evidence.sections.len(), 2);
        assert_eq!(evidence.solids.len(), 2);
        assert_eq!(
            evidence
                .solids
                .iter()
                .map(|solid| (solid.source.data_section, solid.source.id))
                .collect::<Vec<_>>(),
            vec![(0, 1), (1, 1)]
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_reuses_qualified_direct_candidates_for_distinct_complete_data_sections() -> Result<()>
    {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let box_exchange = |size_mm| -> Result<Exchange> {
            let mut model = CadModel::new();
            let root = model.add_node(CadNode::Primitive(Primitive::Box { size_mm }));
            model.add_root(root)?;
            let step = kernel.to_step(&kernel.evaluate(&model, root)?)?;
            Ok(monstertruck_io::step::load::ParsedStepSource::parse(step.as_bytes())?.exchange)
        };
        let mut exchange = box_exchange([2.0, 3.0, 4.0])?;
        exchange.data.extend(box_exchange([5.0, 7.0, 8.0])?.data);
        let build =
            recover_complete_ir_build_state(&exchange, ReducedBrepMutationMode::KernelOwned, None)?;
        let sources = build
            .recovery
            .source_solids
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let [first, second] = sources.as_slice() else {
            bail!("expected two source boxes")
        };
        assert_eq!(first.id, second.id);
        assert_ne!(first.data_section, second.data_section);
        for source in &sources {
            let entities = &exchange.data[source.data_section].entities;
            let index = crate::step_graph::build_index(entities);
            let direct = recover_constructive_fragments_from_entities_with_index(entities, &index)?;
            assert_eq!(build.direct_candidates[source], direct);
            assert!(!direct.is_empty());
        }
        let report = recover_complete_ir_with_kernel_bytes(
            crate::write_exchange(&exchange)?.as_bytes(),
            &kernel,
        )?;
        assert_eq!(report.source_boundary_verified_solids, sources);
        assert_eq!(report.constructively_recovered_solid_count, 2);
        let fragments = report.selected_fragments().collect::<Vec<_>>();
        assert_eq!(fragments.len(), 2);
        for selected in fragments {
            let shape = kernel.evaluate(&selected.fragment.model, selected.fragment.root)?;
            let own = selected.source_solid;
            let other = if own == *first { *second } else { *first };
            assert!(rebuilt_kernel_boundary_matches_source(
                &kernel,
                &shape,
                own.id,
                &exchange.data[own.data_section].entities
            )?);
            assert!(!rebuilt_kernel_boundary_matches_source(
                &kernel,
                &shape,
                other.id,
                &exchange.data[other.data_section].entities
            )?);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn assert_boss_native_edit(
        document: &monstertruck_io::step::load::ResolvedStepDocument,
        source: SourceEntity,
        provenance: &KernelTopologyProvenance,
        edited: &crate::kernel_topology::KernelEditedSourceEvidence,
    ) -> Result<()> {
        let original_face_count = document
            .import_solid(source.data_section, source.id)?
            .solid
            .boundaries
            .iter()
            .map(|shell| shell.faces.len())
            .sum::<usize>();
        assert_eq!(edited.source, source);
        assert_eq!(edited.boundary_shells, 1);
        assert_eq!(edited.faces + 5, original_face_count);
        assert_eq!(
            edited.retained_source_topology.faces.len() + 5,
            provenance.faces.len()
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_v15_square_boss_peel_has_complete_source_and_kernel_proofs() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/boundary_boss_extrusion.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange,
        )?;
        let source_load = kernel.inspect_source_document(&document)?;
        let [loaded] = source_load.solids.as_slice() else {
            bail!("single-boss fixture must import as exactly one source solid");
        };
        let provenance = loaded
            .topology_provenance
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("boss fixture did not retain kernel StableIds"))?;
        let exchange = document.as_ref();
        let section = &exchange.data[loaded.source.data_section];
        let index = crate::step_graph::build_index(&section.entities);
        let shell = crate::brep::manifold_shell(loaded.source.id, &section.entities, &index)
            .ok_or_else(|| anyhow::anyhow!("boss fixture has no solid shell"))?;
        let report = crate::planar_features::diagnose_boundary_features_for_shells_with_index(
            &section.entities,
            &index,
            &HashSet::from([shell]),
        );
        let selected = crate::planar_features::leaf_feature_family_indices(&report.families);
        assert_eq!(
            report
                .families
                .iter()
                .filter(|family| {
                    family.polarity == FeaturePolarity::Additive && family.faces_per_instance == 5
                })
                .count(),
            1,
            "unexpected boss families: {report:?}"
        );
        let targets = LocalRecoveryTargets::from_solids(
            &BTreeSet::from([loaded.source.id]),
            &section.entities,
            &index,
        )?;
        let prepared = prepare_boundary_feature_reduction(
            &section.entities,
            &index,
            &targets,
            &report,
            &selected,
            0,
        )?;
        assert_eq!(
            prepared.family_indices.len(),
            1,
            "boss evidence did not produce exactly one constructive tool: {report:?}"
        );
        let (patch, _) = plan_source_boundary_feature_peel_from_entities(
            loaded.source,
            &section.entities,
            &index,
            &report,
            &prepared.family_indices,
            provenance,
        )?;
        assert_eq!(patch.remove_faces.len(), 5);
        assert_eq!(patch.face_boundary_replacements.len(), 1);
        assert_eq!(patch.require_unused_edges.len(), 4);
        let mut edited = kernel.edit_source_topology(&document, loaded.source, &[])?;
        kernel.apply_edited_source_patches(&mut edited, std::slice::from_ref(&patch))?;
        let committed = kernel.inspect_edited_source(&edited)?;
        // Trying to remove the same faces twice must leave the first valid
        // native transaction intact, without export or re-import.
        assert!(
            kernel
                .apply_edited_source_patches(&mut edited, &[patch])
                .is_err()
        );
        assert_eq!(kernel.inspect_edited_source(&edited)?, committed);
        assert_boss_native_edit(
            &document,
            loaded.source,
            provenance,
            &kernel.inspect_edited_source(&edited)?,
        )?;
        let snapshot = kernel.edited_source_to_step(&edited)?;
        let healed_document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(snapshot.as_bytes())?.exchange,
        )?;
        let healed_report = kernel.inspect_source_document(&healed_document)?;
        assert_eq!(healed_report.solids.len(), 1);
        let host = recognize_single_healed_host(healed_document.as_ref(), &kernel)?
            .ok_or_else(|| anyhow::anyhow!("boss peel did not expose an analytic box host"))?;
        assert!(matches!(
            host.source,
            CadFragmentSource::SolidExtrusion { .. }
        ));
        Ok(())
    }

    /// Relief-bearing manufacturer parts have hole boundaries that may also be
    /// recognized as subordinate feature hypotheses. A closed extrusion tool
    /// already accounts for these boundaries; do not strand the parent.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_shallow_relief_cohort_peels_complete_tools() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let fixtures: &[(&[u8], usize, usize)] = &[
            (
                include_bytes!("../validation/fixtures/native_relief_box_0b5fc9e7.step"),
                230,
                6,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_housing_0036d9d5.step"),
                263,
                30,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_housing_0112d499.step"),
                264,
                40,
            ),
            (
                // The correctly oriented subtractive tools now remove 55
                // additional faces in this first native transaction.
                include_bytes!("../validation/fixtures/native_relief_housing_0914dda4.step"),
                385,
                90,
            ),
        ];
        for &(bytes, original_faces, healed_faces) in fixtures {
            let document = kernel.resolve_source_exchange(
                monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange,
            )?;
            let load = kernel.inspect_source_document(&document)?;
            let [loaded] = load.solids.as_slice() else {
                bail!("expected one native relief solid");
            };
            let provenance = loaded
                .topology_provenance
                .as_ref()
                .context("relief source lacks native identity")?;
            let section = &document.as_ref().data[loaded.source.data_section];
            let index = crate::step_graph::build_index(&section.entities);
            let original = crate::brep::solid_face_ids(loaded.source.id, &section.entities, &index)
                .context("relief source has no faces")?;
            assert_eq!(original.len(), original_faces);
            let targets = LocalRecoveryTargets::from_solids(
                &BTreeSet::from([loaded.source.id]),
                &section.entities,
                &index,
            )?;
            let shells = targets.shell_ids();
            let report = crate::planar_features::diagnose_boundary_features_for_shells_with_index(
                &section.entities,
                &index,
                &shells,
            );
            let prepared = prepare_boundary_feature_reduction(
                &section.entities,
                &index,
                &targets,
                &report,
                &(0..report.families.len()).collect::<Vec<_>>(),
                0,
            )?;
            let (patch, plan) = plan_source_boundary_feature_peel_from_entities(
                loaded.source,
                &section.entities,
                &index,
                &report,
                &prepared.family_indices,
                provenance,
            )?;
            let edited = kernel.edit_source_topology(&document, loaded.source, &[patch])?;
            let native = kernel.inspect_edited_source(&edited)?;
            let mut projected = section.entities.clone();
            crate::planar_features::apply_boundary_feature_peel_plan(
                &mut projected,
                &index,
                &plan,
            )?;
            validate_native_peel_projection(
                loaded.source,
                &original,
                &projected,
                &index,
                provenance,
                loaded.boundaries.len(),
                &native,
            )?;
            assert_eq!(native.faces, healed_faces);
            let reduced = kernel.edited_source_reduced_brep(&edited)?;
            assert!(reduced.validate_closed_topology().is_some());
            if healed_faces == 6 {
                assert!(recognize_native_edited_host(&kernel, &edited)?.is_some());
            }
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_shallow_relief_replays_complete_source_cohort() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let fixtures: &[(&[u8], usize)] = &[
            (
                include_bytes!("../validation/fixtures/native_relief_box_0b5fc9e7.step"),
                230,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_box_16147ee2.step"),
                230,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_box_b974305f.step"),
                230,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_box_587c2a36.step"),
                234,
            ),
        ];
        for &(bytes, faces) in fixtures {
            let document = kernel.resolve_source_exchange(
                monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange,
            )?;
            let load = kernel.inspect_source_document(&document)?;
            let [source] = load.solids.as_slice() else {
                bail!("expected one relief source")
            };
            let entities = &document.as_ref().data[source.source.data_section].entities;
            let index = crate::step_graph::build_index(entities);
            let proof = SourceRecovery::new(source.source.id, entities, &index, &kernel);
            let fragment = recover_kernel_boundary_feature_fragment(
                source.source,
                &load,
                &document,
                &kernel,
                &index,
                proof.as_ref(),
            )?
            .context("native relief did not reconstruct a whole body")?;
            assert!(
                matches!(
                    fragment.boundary_verification,
                    BoundaryVerification::Original
                ),
                "native source proof refused for {:?}",
                source.source
            );
            let fragment = fragment.fragment;
            assert!(matches!(
                fragment.source,
                CadFragmentSource::BoundaryFeatureComposite {
                    passes: 1,
                    feature_families: 9,
                    feature_instances: 13,
                    ..
                }
            ));
            let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.is_single_solid(), "{summary:?}");
            assert_eq!(summary.faces, faces);
            let exchange = monstertruck_io::step::load::ParsedStepSource::parse(
                kernel.to_step(&evaluated)?.as_bytes(),
            )?
            .exchange;
            let reimported = kernel.resolve_source_exchange(exchange)?;
            let reloaded = kernel.inspect_source_document(&reimported)?;
            assert_eq!(reloaded.solids.len(), 1);
            assert!(reloaded.solids[0].failure.is_none());
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn assert_source_verified_relief(bytes: &[u8]) -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        let fragments = report.selected_fragments().collect::<Vec<_>>();
        let [body] = fragments.as_slice() else {
            bail!("expected exactly one selected whole-source reconstruction");
        };
        assert!(kernel_replay_is_single_solid(&kernel, &body.fragment)?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_shallow_relief_is_admitted_by_complete_ir() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_box_0b5fc9e7.step"
        ))
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_round_at_terminal_shoulder_limit_recovers_exactly() -> Result<()> {
        use crate::cad_ir::{
            CadModel, CadNode, EdgeSelector, EdgeTreatmentEdge, EdgeTreatmentProfile, Profile2d,
            RigidTransform,
        };
        use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};

        let mut model = CadModel::new();
        let section = model.add_node(CadNode::Extrude {
            profile: Profile2d::polygon(vec![
                [0., 0.],
                [1., 0.],
                [1., 4.9],
                [0.4, 5.5],
                [0., 5.5],
            ])?,
            vector_mm: [0., 0., 6.65],
        });
        let child = model.add_node(CadNode::Transform {
            child: section,
            transform: RigidTransform {
                matrix: [
                    [0., 0., 1., 0.],
                    [1., 0., 0., 0.],
                    [0., 1., 0., 0.],
                    [0., 0., 0., 1.],
                ],
            },
        });
        let root = model.add_node(CadNode::EdgeTreatment {
            child,
            profile: EdgeTreatmentProfile::Round,
            edges: vec![EdgeTreatmentEdge {
                selector: EdgeSelector::Segment {
                    start_mm: [6.65, 0., 5.5],
                    end_mm: [6.65, 0.4, 5.5],
                },
                radius_mm: 0.6,
            }],
        });
        model.add_root(root)?;
        let kernel = MonstertruckKernel::project_standard();
        let evaluated = kernel.evaluate(&model, root)?;
        let source = kernel.to_step(&evaluated)?;
        assert_source_verified_relief(source.as_bytes())
    }

    /// A real relief-bearing host with a limiting oblique fillet. Check
    /// the composed provider itself, independently of candidate ranking.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_boundary_feature_peel_exposes_treated_host() -> Result<()> {
        let bytes =
            include_bytes!("../validation/fixtures/native_relief_oblique_fillet_ba0460c9.step");
        let recovered = verified_native_boundary_feature_candidate(bytes)?;
        assert!(matches!(
            recovered.fragment.source,
            CadFragmentSource::BoundaryFeatureComposite { .. }
        ));
        assert!(
            recovered
                .fragment
                .model
                .nodes
                .iter()
                .any(|node| matches!(node, CadNode::EdgeTreatment { .. }))
        );
        Ok(())
    }

    // Real high-face-count relief housing, formerly timing out in repeated
    // core recognition. Require complete independent original-boundary proof,
    // not any particular reconstruction family or operation count.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_large_relief_housing_recovery_is_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_housing_271239a1.step"
        ))
    }

    // The host itself requires intersecting planar Booleans. Their rebuilt
    // transverse edges must retain certified finite bounds for the next
    // shallow contact batch.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_on_boolean_host_is_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_boolean_host_22bfe2fc.step"
        ))
    }

    // A native revolution cutter must retain affine planar caps so that
    // the contact interface can be rebuilt before the lettering is replayed.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_with_circular_recesses_is_source_verified() -> Result<()> {
        for bytes in [
            include_bytes!("../validation/fixtures/native_relief_curved_1c75cccd.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_curved_7a4544b3.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_curved_224fa2f0.step").as_slice(),
        ] {
            assert_source_verified_relief(bytes)?;
        }
        Ok(())
    }

    // A translated full disk can acquire a few ulps of axis displacement
    // on its second semicircle. Require complete original-source proof
    // for circular recesses with differing radius and depth.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_with_translated_revolution_caps_is_source_verified() -> Result<()> {
        for bytes in [
            include_bytes!("../validation/fixtures/native_relief_curved_a7efbe91.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_curved_e319f3f9.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_curved_0036d9d5.step").as_slice(),
        ] {
            assert_source_verified_relief(bytes)?;
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_limiting_hollow_analytic_loft_is_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_mixed_cda28227.step"
        ))
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_mixed_sectional_host_is_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_mixed_65712a81.step"
        ))
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_tapered_hosts_are_source_verified() -> Result<()> {
        for bytes in [
            include_bytes!("../validation/fixtures/native_relief_tapered_4eef029b.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_tapered_a9bc1938.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_tapered_6a408847.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_tapered_e43d7cca.step").as_slice(),
        ] {
            assert_source_verified_relief(bytes)?;
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_turning_spline_relief_is_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_turning_09b26184.step"
        ))
    }

    // The host requires parallel extrusion transition algebra. Replayed
    // spline contact intersections must converge beyond the weld radius;
    // competing analytic-only cores cannot cover its polynomial strips.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_polynomial_crossings_and_core_exclusion_are_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_polynomial_622f6bfb.step"
        ))
    }

    // Adjacent canonical cells occur simultaneously in plane direction,
    // plane offset, and the spline/endpoint coordinates of this original.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_mixed_boundary_cells_are_source_verified() -> Result<()> {
        assert_source_verified_relief(include_bytes!(
            "../validation/fixtures/native_relief_mixed_cells_30821d92.step"
        ))
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_partitioned_spline_caps_are_source_verified() -> Result<()> {
        for bytes in [
            include_bytes!("../validation/fixtures/native_relief_partitioned_63d50d37.step")
                .as_slice(),
            include_bytes!("../validation/fixtures/native_relief_partitioned_67cf7259.step")
                .as_slice(),
        ] {
            assert_source_verified_relief(bytes)?;
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_local_contact_rebuilds_stepped_relief_sources() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for bytes in [
            include_bytes!("../validation/fixtures/native_relief_stepped_66eb641f.step").as_slice(),
            include_bytes!("../validation/fixtures/native_relief_stepped_c7c4ea15.step").as_slice(),
        ] {
            let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(report.constructively_recovered_solid_count, 1);
            assert_eq!(report.unrecovered_solid_count, 0);
            assert_eq!(
                report.source_boundary_verified_solids,
                [report.body_scan.solids[0].source]
            );
            let [body] = report.boundary_feature_composites.as_slice() else {
                bail!("stepped relief must recover an executable boundary feature composite");
            };
            let rebuilt = kernel.evaluate(&body.fragment.model, body.fragment.root)?;
            let summary = kernel.summarize(&rebuilt);
            assert!(summary.is_single_solid(), "{summary:?}");
            assert_eq!(summary.faces, 238);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_planar_contact_rebuilds_exact_relief_recesses() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for (bytes, faces) in [
            (
                include_bytes!("../validation/fixtures/native_relief_recess_202872e8.step")
                    .as_slice(),
                246,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_recess_17a62594.step")
                    .as_slice(),
                250,
            ),
            (
                include_bytes!("../validation/fixtures/native_relief_recess_dc36f700.step")
                    .as_slice(),
                254,
            ),
        ] {
            let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(report.constructively_recovered_solid_count, 1);
            assert_eq!(
                report.source_boundary_verified_solids,
                [report.body_scan.solids[0].source],
                "recess original source proof failed ({faces} faces)"
            );
            let [body] = report.boundary_feature_composites.as_slice() else {
                bail!("recess source must recover an executable feature composite");
            };
            let rebuilt = kernel.evaluate(&body.fragment.model, body.fragment.root)?;
            let summary = kernel.summarize(&rebuilt);
            assert!(summary.is_single_solid(), "{summary:?}");
            assert_eq!(summary.faces, faces);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_blind_spline_recesses_reconstruct_whole_body() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_relief_blind_254dc55c.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        let [body] = report.boundary_feature_composites.as_slice() else {
            bail!("blind spline recesses must recover one editable feature composite");
        };
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        assert!(report.source_boundary_verified_solids.is_empty());
        let [reconciliation] = report.source_orientation_reconciled_solids.as_slice() else {
            bail!("expected an independently verified orientation reconciliation");
        };
        assert_eq!(reconciliation.source, report.body_scan.solids[0].source);
        assert_eq!(
            reconciliation
                .repairs
                .iter()
                .map(|r| r.face_id)
                .collect::<Vec<_>>(),
            [530, 729, 845, 1686]
        );
        assert!(
            reconciliation
                .repairs
                .iter()
                .all(|r| !r.original_same_sense && r.reconciled_same_sense)
        );
        assert!(kernel_replay_is_single_solid(&kernel, &body.fragment)?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_spherical_bosses_and_blind_recesses_are_source_verified() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_relief_blind_636f2bd5.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        let [body] = report.boundary_feature_composites.as_slice() else {
            bail!("expected editable subtractive and spherical feature composite");
        };
        let rebuilt = kernel.evaluate(&body.fragment.model, body.fragment.root)?;
        let summary = kernel.summarize(&rebuilt);
        assert!(summary.is_single_solid(), "{summary:?}");
        assert_eq!(summary.faces, 84);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_relief_housing_reconciles_only_certified_strip_face_senses() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_relief_housing_0914dda4.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert!(report.source_boundary_verified_solids.is_empty());
        let [reconciliation] = report.source_orientation_reconciled_solids.as_slice() else {
            bail!("expected a distinct orientation-reconciled recovery certificate");
        };
        assert_eq!(reconciliation.source, report.body_scan.solids[0].source);
        assert_eq!(
            reconciliation
                .repairs
                .iter()
                .map(|r| r.face_id)
                .collect::<Vec<_>>(),
            [799, 1682, 2395, 3004, 3833, 4324, 4496, 4755, 6122, 8592]
        );
        assert!(
            reconciliation
                .repairs
                .iter()
                .all(|r| r.original_same_sense != r.reconciled_same_sense)
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn healed_host_can_use_independently_proven_sectional_recovery() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_sectional_stepped_655821f7.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let host = recognize_single_healed_host(&exchange, &kernel)?
            .context("no constructive healed sectional host")?;
        assert!(matches!(
            host.source,
            CadFragmentSource::SectionalExtrusionUnion { slabs: 5, .. }
                | CadFragmentSource::SectionalExtrusionTransitions { slabs: 5, .. }
        ));
        let evaluated = kernel.evaluate(&host.model, host.root)?;
        assert!(crate::cad_recovery::rebuilt_kernel_boundary_matches_source(
            &kernel,
            &evaluated,
            99,
            &exchange.data[0].entities
        )?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn healed_host_can_use_proven_material_core_union() -> Result<()> {
        let bytes =
            include_bytes!("../validation/fixtures/native_relief_healed_host_0036d9d5.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let exchange = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
        let host = recognize_single_healed_host(&exchange, &kernel)?
            .context("no constructive healed union host")?;
        assert!(is_constructive_source(&host.source), "{:?}", host.source);
        assert_eq!(source_solid_id(&host.source), 16);
        let evaluated = kernel.evaluate(&host.model, host.root)?;
        assert!(crate::cad_recovery::rebuilt_kernel_boundary_matches_source(
            &kernel,
            &evaluated,
            16,
            &exchange.data[0].entities
        )?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn prove_direct_fixture<K: CadKernel>(
        source: SourceEntity,
        exchange: &Exchange,
        kernel: &K,
    ) -> Result<Option<CadFragment>> {
        let entities = &exchange.data[source.data_section].entities;
        let index = crate::step_graph::build_index(entities);
        let Some(proof) = SourceRecovery::new(source.id, entities, &index, kernel) else {
            return Ok(None);
        };
        let mut build =
            recover_complete_ir_build_state(exchange, ReducedBrepMutationMode::KernelOwned, None)?;
        let candidates = build.direct_candidates.remove(&source).unwrap_or_default();
        recover_verified_direct_constructive_fragment(&proof, candidates, None)
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_v15_replays_square_boss_through_complete_ir() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/boundary_boss_extrusion.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.schema, KERNEL_COMPLETE_IR_SCHEMA);
        assert_eq!(
            report.kernel_recovery_failures,
            Vec::<KernelRecoveryFailure>::new()
        );
        assert_eq!(report.unrecovered_solid_count, 0);
        let exact = recover_auto_proven_extrusion_core_union_fragment(input, 0, 16, &kernel)?
            .context("square boss has no exact simple material model")?;
        let chosen = selected_unique_source_fragment(&report)?;
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        assert!(
            chosen.model.complexity_score(chosen.root)?
                <= exact.model.complexity_score(exact.root)?
        );
        let source = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let shape = kernel.evaluate(&chosen.model, chosen.root)?;
        assert!(kernel.summarize(&shape).is_single_solid());
        assert!(rebuilt_kernel_boundary_matches_source(
            &kernel,
            &shape,
            16,
            &source.data[0].entities
        )?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_boss_replay_candidate_has_independent_original_boundary_proof() -> Result<()> {
        verified_native_boundary_feature_candidate(include_bytes!(
            "../validation/fixtures/boundary_boss_extrusion.step"
        ))?;
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn verified_native_boundary_feature_candidate(input: &[u8]) -> Result<KernelRecoveryCandidate> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange,
        )?;
        let load = kernel.inspect_source_document(&document)?;
        let [solid] = load.solids.as_slice() else {
            bail!("expected exactly one source solid");
        };
        let entities = &document.as_ref().data[solid.source.data_section].entities;
        let index = crate::step_graph::build_index(entities);
        let proof = SourceRecovery::new(solid.source.id, entities, &index, &kernel);
        let recovered = recover_kernel_boundary_feature_fragment(
            solid.source,
            &load,
            &document,
            &kernel,
            &index,
            proof.as_ref(),
        )?
        .context("native boundary-feature recovery refused")?;
        assert!(matches!(
            recovered.boundary_verification,
            BoundaryVerification::Original
        ));
        Ok(recovered)
    }

    /// Exact kernel adapter for native-peeling tests: any intermediate or final
    /// STEP serialization makes the test fail instead of silently falling back.
    #[cfg(feature = "cad-kernel-monstertruck")]
    struct NoSourceSnapshotKernel(crate::cad_kernel::monstertruck::MonstertruckKernel);

    #[cfg(feature = "cad-kernel-monstertruck")]
    impl CadKernel for NoSourceSnapshotKernel {
        type Evaluated =
            <crate::cad_kernel::monstertruck::MonstertruckKernel as CadKernel>::Evaluated;

        fn evaluate(&self, model: &CadModel, root: NodeId) -> Result<Self::Evaluated> {
            self.0.evaluate(model, root)
        }

        fn summarize(&self, evaluated: &Self::Evaluated) -> crate::cad_kernel::KernelSummary {
            self.0.summarize(evaluated)
        }

        fn supports_step_export(&self) -> bool { false }

        fn to_step(&self, _: &Self::Evaluated) -> Result<String> {
            bail!("native feature peeling must not require CAD STEP serialization")
        }
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    impl SourceStepLoadInspector for NoSourceSnapshotKernel {
        type Document = <crate::cad_kernel::monstertruck::MonstertruckKernel as SourceStepLoadInspector>::Document;

        fn resolve_source_exchange(&self, exchange: Exchange) -> Result<Self::Document> {
            self.0.resolve_source_exchange(exchange)
        }

        fn inspect_source_document(
            &self,
            document: &Self::Document,
        ) -> Result<KernelSourceLoadEvidence> {
            self.0.inspect_source_document(document)
        }
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    impl SourceTopologyEditor for NoSourceSnapshotKernel {
        type EditedSource = <crate::cad_kernel::monstertruck::MonstertruckKernel as SourceTopologyEditor>::EditedSource;

        fn edit_source_topology(
            &self,
            document: &Self::Document,
            source: SourceEntity,
            patches: &[crate::kernel_topology::KernelTopologyPatch],
        ) -> Result<Self::EditedSource> {
            self.0.edit_source_topology(document, source, patches)
        }

        fn apply_edited_source_patches(
            &self,
            edited: &mut Self::EditedSource,
            patches: &[crate::kernel_topology::KernelTopologyPatch],
        ) -> Result<()> {
            self.0.apply_edited_source_patches(edited, patches)
        }

        fn normalize_edited_source(&self, edited: &mut Self::EditedSource) -> Result<usize> {
            self.0.normalize_edited_source(edited)
        }

        fn edited_source_reduced_brep(
            &self,
            edited: &Self::EditedSource,
        ) -> Result<crate::reduced_brep::ReducedSolid> {
            self.0.edited_source_reduced_brep(edited)
        }

        fn inspect_edited_source(
            &self,
            edited: &Self::EditedSource,
        ) -> Result<crate::kernel_topology::KernelEditedSourceEvidence> {
            self.0.inspect_edited_source(edited)
        }

        fn edited_source_step_projection(
            &self,
            _: &Self::EditedSource,
        ) -> Result<crate::kernel_topology::KernelEditedSourceProjection> {
            bail!("native feature peeling must not export an intermediate STEP projection")
        }

        fn edited_source_to_step(&self, _: &Self::EditedSource) -> Result<String> {
            bail!("native feature peeling must not export intermediate STEP")
        }
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_nested_boss_peels_two_passes_with_original_native_identity() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/boundary_boss_nested.step");
        let kernel = NoSourceSnapshotKernel(
            crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard(),
        );
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert!(
            report.kernel_recovery_failures.is_empty(),
            "{:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.unrecovered_solid_count, 0);
        let [composite] = report.boundary_feature_composites.as_slice() else {
            bail!("kernel-owned nested boss did not reconstruct the boundary composite");
        };
        assert!(matches!(
            composite.fragment.source,
            CadFragmentSource::BoundaryFeatureComposite {
                passes: 2,
                feature_families: 2,
                feature_instances: 2,
                ..
            }
        ));
        composite.fragment.model.validate()?;
        assert!(kernel_replay_is_single_solid(&kernel, &composite.fragment)?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_curved_extrusions_reconstruct_executable_cad_from_real_step_exchange() -> Result<()> {
        use crate::brep::SurfaceSupport;
        use crate::cad_ir::{Curve2d, Profile2d, ProfileLoop};

        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let rectangle_tail = |right: f64| {
            [
                Curve2d::Line {
                    start_mm: [right, 0.0],
                    end_mm: [right, -2.0],
                },
                Curve2d::Line {
                    start_mm: [right, -2.0],
                    end_mm: [0.0, -2.0],
                },
                Curve2d::Line {
                    start_mm: [0.0, -2.0],
                    end_mm: [0.0, 0.0],
                },
            ]
        };
        let mut profiles = Vec::new();
        profiles.push((
            "circular wall",
            vec![Curve2d::CircleArc {
                center_mm: [0.0, 0.0],
                radius_mm: 2.0,
                start_angle_rad: 0.0,
                end_angle_rad: std::f64::consts::TAU,
            }],
        ));
        let mut bezier = vec![Curve2d::Bezier {
            control_points_mm: vec![[0.0, 0.0], [1.0, 1.0], [2.0, 1.0], [3.0, 0.0]],
        }];
        bezier.extend(rectangle_tail(3.0));
        profiles.push(("Bezier wall", bezier));
        let mut rational = vec![Curve2d::BSpline {
            degree: 2,
            control_points_mm: vec![[0.0, 0.0], [1.0, 1.5], [2.0, 0.0]],
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            weights: Some(vec![1.0, 0.7, 1.0]),
        }];
        rational.extend(rectangle_tail(2.0));
        profiles.push(("rational spline wall", rational));

        for (label, curves) in profiles {
            let mut model = CadModel::new();
            let source_root = model.add_node(CadNode::Extrude {
                profile: Profile2d {
                    loops: vec![ProfileLoop { curves }],
                },
                vector_mm: [0.0, 0.0, 4.0],
            });
            model.add_root(source_root)?;
            let original = kernel.evaluate(&model, source_root)?;
            assert!(kernel.summarize(&original).is_single_solid(), "{label}");

            let step = kernel.to_step(&original)?;
            let document = kernel.resolve_source_exchange(
                monstertruck_io::step::load::ParsedStepSource::parse(step.as_bytes())?.exchange,
            )?;
            let inspected = kernel.inspect_source_document(&document)?;
            let [load] = inspected.solids.as_slice() else {
                bail!("{label}: expected exactly one STEP source solid");
            };
            let native = kernel.edit_source_topology(&document, load.source, &[])?;
            let reduced = kernel.edited_source_reduced_brep(&native)?;
            assert!(
                reduced.faces().any(|face| matches!(
                    &face.surface,
                    SurfaceSupport::SplineExtrusion(_) | SurfaceSupport::Cylinder(_)
                )),
                "{label}: native side support was not recognized: {:?}",
                reduced
                    .faces()
                    .map(|face| &face.surface)
                    .collect::<Vec<_>>()
            );
            let host = recognize_native_edited_host(&kernel, &native)?
                .ok_or_else(|| anyhow::anyhow!("{label}: native CAD reconstruction failed"))?;
            host.model.validate()?;
            let rebuilt = kernel.evaluate(&host.model, host.root)?;
            assert!(
                kernel.summarize(&rebuilt).is_single_solid(),
                "{label}: recovered CAD did not evaluate to a valid solid"
            );
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_curved_u_sweep_recovers_native_rational_extrusion_without_step_snapshot() -> Result<()>
    {
        // A real STEP solid whose caps use a three-pole rational quadratic
        // arc but whose exactly swept side describes the *same* arc with five
        // poles and a different knot basis. Both cap faces are encoded as
        // B-spline surfaces, not as PLANE records.
        let input = include_bytes!("../validation/fixtures/native_rational_u_sweep.step");
        let kernel = NoSourceSnapshotKernel(
            crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard(),
        );
        let result = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(result.body_count, 1);
        assert_eq!(
            result.constructively_recovered_solid_count, 1,
            "{:?}",
            result.kernel_recovery_failures
        );
        assert_eq!(result.unrecovered_solid_count, 0);
        let [recovered] = result.solid_extrusions.as_slice() else {
            bail!("native U-sweep must recover exactly one executable extrusion");
        };
        let fragment = &recovered.fragment;
        assert!(matches!(
            fragment.source,
            CadFragmentSource::SolidExtrusion { .. }
        ));
        fragment.model.validate()?;
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(
            kernel.summarize(&rebuilt).is_single_solid(),
            "native U-sweep CAD did not rebuild a single geometric solid"
        );
        Ok(())
    }

    /// Four distinct real JLCPCB closed solids comprised of two round
    /// elbows, three cylindrical legs and planar circular end caps.
    /// No source-specific dimensions, face count or bend count is encoded
    /// in the detector; these fixtures cover exporter face partitioning.
    /// A real JLCPCB 14-face rectangular bar with independent four-sided
    /// chamfers at both ends: four coplanar, line-edged rectangle stations
    /// spanning a long constant section and two short tapered sections.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_staged_rectangle_recovers_complete_native_ruled_loft() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/native_staged_rectangular_loft.step");
        let parser = recover_complete_ir_bytes(input)?;
        assert_eq!(parser.constructively_recovered_solid_count, 1);
        let kernel = NoSourceSnapshotKernel(
            crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard(),
        );
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.constructively_recovered_solid_count, 1,
            "{:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.unrecovered_solid_count, 0);
        let [loft] = report.solid_lofts.as_slice() else {
            bail!("expected one complete parametric loft");
        };
        let fragment = &loft.fragment;
        let [CadNode::Loft { sections }] = fragment.model.nodes.as_slice() else {
            bail!("expected one editable loft node");
        };
        assert_eq!(sections.len(), 4);
        fragment.model.validate()?;
        let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
        assert!(kernel.summarize(&rebuilt).is_single_solid());
        Ok(())
    }

    /// The real JLC connector's four independent rectangular stations are
    /// recognized without source-specific dimensions. A valid native loft
    /// alone is not an automatic source admission: the exported *oriented*
    /// analytic boundary must independently match the original solid.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_staged_connector_loft_is_source_boundary_proven() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/native_staged_rectangular_loft.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let source = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let [section] = source.data.as_slice() else {
            bail!("JLC connector source needs exactly one DATA section");
        };
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.unrecovered_solid_count, 0,
            "{:?}",
            report.kernel_recovery_failures
        );
        let [loft] = report.solid_lofts.as_slice() else {
            bail!("JLC connector did not admit precisely one editable loft");
        };
        let CadFragmentSource::SolidLoft { solid_id, .. } = loft.fragment.source else {
            bail!("JLC connector selected an unexpected CAD source");
        };
        loft.fragment.model.validate()?;
        let rebuilt = kernel.evaluate(&loft.fragment.model, loft.fragment.root)?;
        assert!(kernel.summarize(&rebuilt).is_single_solid());
        assert!(
            rebuilt_kernel_boundary_matches_source(&kernel, &rebuilt, solid_id, &section.entities)?,
            "automatically recovered JLC staged loft must match its complete oriented source boundary"
        );
        assert!(
            prove_direct_fixture(report.body_scan.solids[0].source, &source, &kernel)?.is_some(),
            "the generic source-admission path must independently certify this JLC loft"
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_open_round_bend_paths_recover_as_executable_cad_sweeps() -> Result<()> {
        let kernel = NoSourceSnapshotKernel(
            crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard(),
        );
        for (name, input) in [
            (
                "a",
                include_bytes!("../validation/fixtures/open_round_two_bend_a.step").as_slice(),
            ),
            (
                "b",
                include_bytes!("../validation/fixtures/open_round_two_bend_b.step").as_slice(),
            ),
            (
                "c",
                include_bytes!("../validation/fixtures/open_round_two_bend_c.step").as_slice(),
            ),
            (
                "d",
                include_bytes!("../validation/fixtures/open_round_two_bend_d.step").as_slice(),
            ),
        ] {
            let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
            assert_eq!(report.body_count, 1, "{name}");
            assert_eq!(
                report.constructively_recovered_solid_count, 1,
                "{name}: {:?}",
                report.kernel_recovery_failures
            );
            assert_eq!(report.unrecovered_solid_count, 0, "{name}");
            let [recovered] = report.solid_sweeps.as_slice() else {
                bail!("{name}: expected precisely one complete sweep");
            };
            assert!(
                matches!(
                    recovered.fragment.source,
                    CadFragmentSource::SolidSweep { closed: false, .. }
                ),
                "{name}"
            );
            let fragment = &recovered.fragment;
            let [crate::cad_ir::CadNode::Sweep { path, .. }] = fragment.model.nodes.as_slice()
            else {
                bail!("{name}: expected one editable sweep node");
            };
            assert!(!path.closed, "{name}");
            assert_eq!(path.segments.len(), 5, "{name}");
            fragment.model.validate()?;
            let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&rebuilt).is_single_solid(), "{name}");
        }
        Ok(())
    }

    /// The ninth real source must be discovered by the general recognizer;
    /// a separately source-constructed fixture is insufficient for Complete IR.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_bored_cap_with_ruled_lobes_is_automatically_complete_ir() -> Result<()> {
        let source =
            include_bytes!("../validation/fixtures/native_mixed_carrier_chamfer_06d018ae.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(source, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.constructively_recovered_solid_count, 1,
            "{:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.unrecovered_solid_count, 0);
        let chosen = selected_unique_source_fragment(&report)?;
        assert_eq!(
            chosen.source,
            CadFragmentSource::CapExtrusionRuledEnds {
                solid_id: 19,
                cap_face_id: 400,
                tip_face_ids: vec![301, 416],
            }
        );
        assert_eq!(
            chosen
                .model
                .nodes
                .iter()
                .filter(|node| matches!(node, CadNode::Loft { .. }))
                .count(),
            2,
        );
        chosen.model.validate()?;
        let rebuilt = kernel.evaluate(&chosen.model, chosen.root)?;
        assert!(kernel.summarize(&rebuilt).is_single_solid());
        let parsed = monstertruck_io::step::load::ParsedStepSource::parse(source)?.exchange;
        assert!(rebuilt_kernel_boundary_matches_source(
            &kernel,
            &rebuilt,
            19,
            &parsed.data[0].entities,
        )?);
        Ok(())
    }

    /// Historical corpus case 99527c was unrecovered before the general
    /// source-derived union improvements. Its independently published STEP
    /// manifold is the acceptance oracle, not any implementation's node order.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_four_core_union_99527c_is_automatically_source_recovered() -> Result<()> {
        use crate::cad_kernel::CadKernel;

        let input = include_bytes!("../validation/fixtures/native_four_core_union_99527c.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.constructively_recovered_solid_count, 1,
            "source solid #56 must automatically recover: {:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.unrecovered_solid_count, 0);
        assert_eq!(report.body_scan.solids[0].source.id, 56);

        let selected = selected_unique_source_fragment(&report)?;
        selected.model.validate()?;
        let rebuilt = kernel.evaluate(&selected.model, selected.root)?;
        assert!(kernel.summarize(&rebuilt).is_single_solid());
        let source = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let [source_section] = source.data.as_slice() else {
            bail!("four-core source did not have one DATA section");
        };
        assert!(
            rebuilt_kernel_boundary_matches_source(
                &kernel,
                &rebuilt,
                56,
                &source_section.entities,
            )?,
            "automatically selected editable model is not source-boundary equivalent"
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_open_cap_and_hook_union_reconstructs_complete_editable_body() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/native_open_cap_hook_union.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.unrecovered_solid_count, 0,
            "{:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.extrusion_core_unions.len(), 1);
        let union = &report.extrusion_core_unions[0].fragment;
        assert!(matches!(
            union.source,
            CadFragmentSource::ExtrusionCoreUnion { solid_id: 30, .. }
        ));
        union.model.validate()?;
        let rebuilt = kernel.evaluate(&union.model, union.root)?;
        assert!(kernel.summarize(&rebuilt).is_single_solid());
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn selected_unique_source_fragment(report: &CompleteIr) -> Result<&CadFragment> {
        let selected = report.selected_fragments().collect::<Vec<_>>();
        let [chosen] = selected.as_slice() else {
            bail!(
                "expected exactly one selected source CAD model, got {}",
                selected.len()
            );
        };
        Ok(&chosen.fragment)
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[derive(Clone, Copy)]
    enum TargetedFamily {
        Loft,
        SourceCapUnion,
        CoreRuledEnds,
        CapRuledEnds,
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    impl TargetedFamily {
        fn matches(self, source: &CadFragmentSource, solid_id: u64) -> bool {
            source_solid_id(source) == solid_id
                && matches!(
                    (self, source),
                    (Self::Loft, CadFragmentSource::SolidLoft { .. })
                        | (
                            Self::SourceCapUnion,
                            CadFragmentSource::SourceCapExtrusionUnion { .. }
                        )
                        | (
                            Self::CoreRuledEnds,
                            CadFragmentSource::ExtrusionCoreRuledEnds { .. }
                        )
                        | (
                            Self::CapRuledEnds,
                            CadFragmentSource::CapExtrusionRuledEnds { .. }
                        )
                )
        }
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    fn verify_targeted_analytic_sources(
        samples: &[(&str, &[u8], u64, TargetedFamily)],
    ) -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for &(label, input, solid_id, family) in samples {
            let exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
            let [section] = exchange.data.as_slice() else {
                bail!("{label}: expected sole source DATA section");
            };
            let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
            assert_eq!(report.body_count, 1, "{label}");
            assert_eq!(
                report.constructively_recovered_solid_count, 1,
                "{label}: {:?}",
                report.kernel_recovery_failures
            );
            let fragment = selected_unique_source_fragment(&report)?;
            assert!(
                family.matches(&fragment.source, solid_id),
                "{label}: unexpected constructive family {:?}",
                fragment.source
            );
            fragment.model.validate()?;
            let rebuilt = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&rebuilt).is_single_solid(), "{label}");
            assert!(
                rebuilt_kernel_boundary_matches_source(
                    &kernel,
                    &rebuilt,
                    solid_id,
                    &section.entities
                )?,
                "{label}: source boundary differs"
            );
        }
        Ok(())
    }

    /// Representative previously unrecovered source bodies from the frozen
    /// catalog census. Admit their whole editable solid, not just a candidate
    /// with a valid topology. Both the loft and constructive union paths
    /// must reproduce the original oriented analytic boundary.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_targeted_analytic_sources_are_automatically_recovered() -> Result<()> {
        let samples: &[(&str, &[u8], u64, TargetedFamily)] = &[
            (
                "42c4c451",
                include_bytes!("../validation/fixtures/native_targeted_42c4c451.step"),
                22,
                TargetedFamily::Loft,
            ),
            (
                "5c8c706c",
                include_bytes!("../validation/fixtures/native_targeted_5c8c706c.step"),
                50,
                TargetedFamily::Loft,
            ),
            (
                "66abf5fd",
                include_bytes!("../validation/fixtures/native_targeted_66abf5fd.step"),
                38,
                TargetedFamily::Loft,
            ),
            (
                "eb76f3ce",
                include_bytes!("../validation/fixtures/native_targeted_eb76f3ce.step"),
                20,
                TargetedFamily::Loft,
            ),
            (
                "58dd1537",
                include_bytes!("../validation/fixtures/native_targeted_58dd1537.step"),
                66,
                TargetedFamily::SourceCapUnion,
            ),
            (
                "b13fc76c",
                include_bytes!("../validation/fixtures/native_targeted_b13fc76c.step"),
                54,
                TargetedFamily::SourceCapUnion,
            ),
        ];
        verify_targeted_analytic_sources(samples)
    }

    /// Nine additional historical failures, all now representable as exact
    /// lofts or extrusions with source-witnessed ruled ends. Keep whole-source
    /// equivalence as the acceptance condition rather than matching only caps.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_low_face_analytic_sources_are_automatically_recovered() -> Result<()> {
        verify_targeted_analytic_sources(&[
            (
                "9288385a",
                include_bytes!("../validation/fixtures/native_targeted_9288385a.step"),
                24,
                TargetedFamily::Loft,
            ),
            (
                "119a071d",
                include_bytes!("../validation/fixtures/native_targeted_119a071d.step"),
                40,
                TargetedFamily::CoreRuledEnds,
            ),
            (
                "14116ee2",
                include_bytes!("../validation/fixtures/native_targeted_14116ee2.step"),
                60,
                TargetedFamily::CoreRuledEnds,
            ),
            (
                "2b38a50f",
                include_bytes!("../validation/fixtures/native_targeted_2b38a50f.step"),
                6,
                TargetedFamily::Loft,
            ),
            (
                "3752a555",
                include_bytes!("../validation/fixtures/native_targeted_3752a555.step"),
                139,
                TargetedFamily::CapRuledEnds,
            ),
            (
                "49a2725d",
                include_bytes!("../validation/fixtures/native_targeted_49a2725d.step"),
                8,
                TargetedFamily::Loft,
            ),
            (
                "530c77f2",
                include_bytes!("../validation/fixtures/native_targeted_530c77f2.step"),
                38,
                TargetedFamily::Loft,
            ),
            (
                "7b59e66c",
                include_bytes!("../validation/fixtures/native_targeted_7b59e66c.step"),
                37,
                TargetedFamily::CoreRuledEnds,
            ),
            (
                "c146449c",
                include_bytes!("../validation/fixtures/native_targeted_c146449c.step"),
                53,
                TargetedFamily::Loft,
            ),
            (
                "ed8dbb8d",
                include_bytes!("../validation/fixtures/native_targeted_ed8dbb8d.step"),
                16,
                TargetedFamily::Loft,
            ),
        ])
    }

    /// Three distinct real JLCPCB solids with long source profiles (including
    /// three cylindrical side patches on one) were structurally recognized
    /// before the generic independent candidate-proof path existed. Their
    /// *entire* native extrusion, rather than a cheaper guess or a subset, must
    /// be admitted only after oriented source-boundary equality succeeds.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_distinct_extrusions_are_source_boundary_proven_in_complete_ir() -> Result<()> {
        let samples: &[(&str, &[u8], u64, [u64; 2])] = &[
            (
                "17-edge planar profile",
                include_bytes!("../validation/fixtures/native_planar_extrusion_bc1b8e95.step")
                    .as_slice(),
                44,
                [336, 354],
            ),
            (
                "11-edge planar profile",
                include_bytes!("../validation/fixtures/native_planar_extrusion_23c65eb0.step")
                    .as_slice(),
                70,
                [136, 241],
            ),
            (
                "nine-edge profile with three circular arcs",
                include_bytes!("../validation/fixtures/native_curved_extrusion_406bfe61.step")
                    .as_slice(),
                205,
                [111, 193],
            ),
        ];
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for &(label, source, solid_id, expected_caps) in samples {
            let exchange = monstertruck_io::step::load::ParsedStepSource::parse(source)?.exchange;
            let [section] = exchange.data.as_slice() else {
                bail!("{label}: source did not contain exactly one DATA section");
            };
            let report = recover_complete_ir_with_kernel_bytes(source, &kernel)?;
            assert_eq!(report.body_count, 1, "{label}");
            assert_eq!(
                report.unrecovered_solid_count, 0,
                "{label}: {:?}",
                report.kernel_recovery_failures
            );
            let selected = selected_unique_source_fragment(&report)?;
            let CadFragmentSource::SolidExtrusion {
                solid_id: recovered_id,
                cap_face_ids,
                ..
            } = &selected.source
            else {
                bail!(
                    "{label}: expected an editable source extrusion, got {:?}",
                    selected.source
                );
            };
            assert_eq!(*recovered_id, solid_id, "{label}");
            let mut caps = *cap_face_ids;
            caps.sort_unstable();
            assert_eq!(
                caps, expected_caps,
                "{label}: source cap provenance changed"
            );
            selected.model.validate()?;
            let rebuilt = kernel.evaluate(&selected.model, selected.root)?;
            assert!(kernel.summarize(&rebuilt).is_single_solid(), "{label}");
            assert!(
                rebuilt_kernel_boundary_matches_source(
                    &kernel,
                    &rebuilt,
                    solid_id,
                    &section.entities,
                )?,
                "{label}: the selected CAD failed strict source-boundary equality"
            );
            let source_entity = report.body_scan.solids[0].source;
            assert!(
                prove_direct_fixture(source_entity, &exchange, &kernel)?.is_some(),
                "{label}: the independent direct-candidate proof did not admit the source"
            );
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_beveled_l_ruled_ends_reconstruct_complete_editable_body() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/native_beveled_l_hook.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.unrecovered_solid_count, 0,
            "{:?}",
            report.kernel_recovery_failures
        );
        assert_eq!(report.constructively_recovered_solid_count, 1);
        // Complete IR is allowed to prefer a cheaper verified primitive,
        // sweep or source-derived core union. Do not pin an inferior
        // representation simply because it happened to be reconstructed
        // first; test the *chosen* constructive model against source geometry.
        let recovered = selected_unique_source_fragment(&report)?;
        recovered.model.validate()?;
        let evaluated = kernel.evaluate(&recovered.model, recovered.root)?;
        assert!(kernel.summarize(&evaluated).is_single_solid());

        let output = kernel.to_step(&evaluated)?;
        let source_exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let rebuilt_exchange =
            monstertruck_io::step::load::ParsedStepSource::parse(output.as_bytes())?.exchange;
        let source_entities = &source_exchange.data[0].entities;
        let rebuilt_entities = &rebuilt_exchange.data[0].entities;
        let output_solid_ids = rebuilt_entities
            .iter()
            .filter_map(|entity| {
                crate::step_graph::simple_record(entity)
                    .filter(|record| record.name == "MANIFOLD_SOLID_BREP")
                    .map(|_| crate::step_graph::entity_id(entity))
            })
            .collect::<Vec<_>>();
        let [output_solid_id] = output_solid_ids.as_slice() else {
            bail!("recovered complete IR did not export exactly one STEP solid");
        };
        let comparison = crate::boundary_proof::compare_analytic_boundaries(
            84,
            source_entities,
            *output_solid_id,
            rebuilt_entities,
        );
        assert!(
            comparison.equivalent,
            "Complete IR selected {:?}, but its exported boundary does not equal the source: {:?}",
            recovered.source, comparison
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_corpus_ssi_orientation_regression_122_recovers_in_source_order() -> Result<()> {
        // Monstertruck PR #57 reversed 24 out of 78 geometrically identical
        // SSI chains. The source-order four-extrusion union must be stable.
        let input = include_bytes!("../validation/fixtures/ssi_oriented_union_case122.step");
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(
            report.unrecovered_solid_count, 0,
            "{:?}",
            report.kernel_recovery_failures
        );
        let explicit = crate::cad_recovery::recover_proven_extrusion_core_union_fragment(
            input,
            0,
            145,
            &[[62, 128], [199, 317], [233, 287], [305, 362]],
            &kernel,
        )?
        .context("source-order four-core SSI union stopped reproducing the source")?;
        let selected = selected_unique_source_fragment(&report)?;
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        assert!(
            selected.model.complexity_score(selected.root)?
                < explicit.model.complexity_score(explicit.root)?
        );
        let original = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let shape = kernel.evaluate(&selected.model, selected.root)?;
        assert!(kernel.summarize(&shape).is_single_solid());
        assert!(rebuilt_kernel_boundary_matches_source(
            &kernel,
            &shape,
            145,
            &original.data[0].entities
        )?);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_v14_owns_coplanar_normalization_without_parser_topology_mutation() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/coplanar_partitioned_extrusion.step");
        let parser = recover_complete_ir_bytes(input)?;
        assert_eq!(parser.schema, COMPLETE_IR_SCHEMA);
        assert_eq!(parser.semantically_recovered_solid_count, 1);
        assert_eq!(parser.unrecovered_solid_count, 0);
        assert_eq!(parser.solid_extrusions.len(), 1);
        assert_eq!(
            parser.same_domain_normalized_solids,
            Vec::<SourceCadFragment>::new()
        );

        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange,
        )?;
        let source_load = kernel.inspect_source_document(&document)?;
        let [solid] = source_load.solids.as_slice() else {
            bail!("coplanar fixture should load as exactly one source solid");
        };
        assert!(solid.failure.is_none());
        let provenance = solid
            .topology_provenance
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("coplanar fixture omitted topology provenance"))?;
        let patch =
            plan_source_coplanar_normalization_patch(document.as_ref(), solid.source, provenance)?
                .ok_or_else(|| {
                    anyhow::anyhow!("coplanar fixture did not produce a normalization patch")
                })?;
        assert_eq!(patch.remove_faces.len(), 1);
        assert_eq!(patch.face_boundary_replacements.len(), 1);
        assert_eq!(patch.require_unused_edges.len(), 1);
        assert_eq!(patch.new_vertices, Vec::<[f64; 3]>::new());
        assert_eq!(
            patch.new_edges,
            Vec::<crate::kernel_topology::KernelNewEdge>::new()
        );
        assert_eq!(
            patch.rewrite_edges,
            Vec::<crate::kernel_topology::KernelRetainedEdge>::new()
        );
        assert_eq!(
            patch.face_splices,
            Vec::<crate::kernel_topology::KernelFaceSplice>::new()
        );

        let mut native = kernel.edit_source_topology(&document, solid.source, &[])?;
        let before = kernel.inspect_edited_source(&native)?;
        assert_eq!(kernel.normalize_edited_source(&mut native)?, 1);
        let after = kernel.inspect_edited_source(&native)?;
        assert_eq!(before.faces, after.faces + 1);
        assert_eq!(after.boundary_shells, before.boundary_shells);
        assert_eq!(
            after.retained_source_topology.faces.len() + 1,
            provenance.faces.len()
        );
        assert_eq!(kernel.normalize_edited_source(&mut native)?, 0);
        let direct = kernel.edited_source_reduced_brep(&native)?;
        let exact_host = crate::solid_extrusions::detect_solid_extrusion_from_reduced(&direct)
            .ok_or_else(|| anyhow::anyhow!("native planar extrusion recognition failed"))?;
        assert_eq!(exact_host.height_mm, 1.0);
        let native_fragment = recognize_native_edited_host(&kernel, &native)?.ok_or_else(|| {
            anyhow::anyhow!("native edited host was not constructively recovered")
        })?;
        native_fragment.model.validate()?;
        let legacy_snapshot = kernel.edited_source_to_step(&native)?;
        let legacy_document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(legacy_snapshot.as_bytes())?
                .exchange,
        )?;
        let legacy = recognize_single_healed_host(legacy_document.as_ref(), &kernel)?
            .ok_or_else(|| anyhow::anyhow!("legacy STEP host was not constructively recovered"))?;
        assert!(kernel_replay_is_single_solid(&kernel, &native_fragment)?);
        assert_eq!(
            kernel.summarize(&kernel.evaluate(&native_fragment.model, native_fragment.root)?),
            kernel.summarize(&kernel.evaluate(&legacy.model, legacy.root)?),
            "native and previously exported host must have equal evaluated topology"
        );

        let report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(report.schema, KERNEL_COMPLETE_IR_SCHEMA);
        assert_eq!(report.semantically_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        // The native normalization is proven independently above. Complete
        // IR may prefer a simpler body if the *original source* independently
        // admits that entire constructive boundary.
        let exact =
            recover_auto_proven_extrusion_core_union_fragment(input, 0, solid.source.id, &kernel)?
                .ok_or_else(|| anyhow::anyhow!("normalized source has no exact simple CAD"))?;
        let chosen = selected_unique_source_fragment(&report)?;
        assert_eq!(
            *chosen, exact,
            "source-proof selection changed the admitted CAD"
        );
        assert_eq!(
            report.kernel_recovery_failures,
            Vec::<KernelRecoveryFailure>::new()
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn native_turning_hosts_recover_from_edited_analytic_brep() -> Result<()> {
        use crate::brep::SurfaceSupport;
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for (name, input) in [
            (
                "conical frustum",
                include_bytes!("../validation/fixtures/native_conical_frustum.step").as_slice(),
            ),
            (
                "line-generated revolution",
                include_bytes!(
                    "../validation/fixtures/native_line_surface_of_revolution_frustum.step"
                )
                .as_slice(),
            ),
            (
                "hollow conical frustum",
                include_bytes!("../validation/fixtures/native_hollow_conical_frustum.step")
                    .as_slice(),
            ),
            (
                "ring torus",
                include_bytes!("../validation/fixtures/native_torus.step").as_slice(),
            ),
            (
                "spherical cap",
                include_bytes!("../validation/fixtures/native_spherical_cap.step").as_slice(),
            ),
        ] {
            let document = kernel.resolve_source_exchange(
                monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange,
            )?;
            let source_load = kernel.inspect_source_document(&document)?;
            let [loaded] = source_load.solids.as_slice() else {
                bail!("{name} should load exactly one solid");
            };
            if let Some(failure) = &loaded.failure {
                bail!("{name} native STEP import failed: {failure:?}");
            }
            let edited = kernel.edit_source_topology(&document, loaded.source, &[])?;
            let native = kernel.edited_source_reduced_brep(&edited)?;
            let analytic = native
                .faces()
                .filter(|face| {
                    matches!(
                        face.surface,
                        SurfaceSupport::Cone(_)
                            | SurfaceSupport::Cylinder(_)
                            | SurfaceSupport::Sphere(_)
                            | SurfaceSupport::Torus(_)
                    )
                })
                .count();
            assert!(analytic > 0, "{name} lost native analytic surfaces");
            let host = recognize_native_edited_host(&kernel, &edited)?
                .ok_or_else(|| anyhow::anyhow!("{name} was not reconstructed from native B-rep"))?;
            host.model.validate()?;
            let evaluated = kernel.evaluate(&host.model, host.root)?;
            let summary = kernel.summarize(&evaluated);
            assert!(summary.is_single_solid(), "{name} produced {summary:?}");
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn kernel_complete_ir_v14_attaches_parser_neutral_topology_provenance() -> Result<()> {
        let input = include_bytes!("../validation/fixtures/production_radial_slot.step");
        let parser_only = recover_complete_ir_bytes(input)?;
        assert_eq!(parser_only.schema, COMPLETE_IR_SCHEMA);
        assert!(parser_only.kernel_source_load.is_none());
        assert!(
            serde_json::to_value(&parser_only)?
                .get("kernel_source_load")
                .is_none()
        );

        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let kernel_report = recover_complete_ir_with_kernel_bytes(input, &kernel)?;
        assert_eq!(kernel_report.schema, KERNEL_COMPLETE_IR_SCHEMA);
        let evidence = kernel_report
            .kernel_source_load
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("kernel v14 omitted source-load evidence"))?;
        assert_eq!(evidence.sections.len(), 1);
        let solid = evidence
            .solids
            .iter()
            .find(|solid| solid.source.id == 15_385)
            .ok_or_else(|| anyhow::anyhow!("missing source-load evidence for solid #15385"))?;
        assert_eq!(solid.source.data_section, 0);
        assert!(solid.failure.is_none());
        assert_eq!(
            solid.boundaries.first().map(|boundary| boundary.role),
            Some(crate::kernel_source_load::KernelBoundaryRole::Outer)
        );
        let provenance = solid
            .topology_provenance
            .as_ref()
            .ok_or_else(|| anyhow::anyhow!("kernel v14 omitted topology provenance"))?;
        assert_ne!(provenance.vertices.len(), 0);
        assert_ne!(provenance.edges.len(), 0);
        assert_ne!(provenance.faces.len(), 0);
        assert!(
            provenance
                .vertices
                .iter()
                .all(|entry| { entry.source.data_section == 0 && entry.topology_id.raw() != 0 })
        );
        assert!(
            provenance
                .edges
                .iter()
                .all(|entry| { entry.source.data_section == 0 && entry.topology_id.raw() != 0 })
        );
        assert!(provenance.faces.iter().all(|entry| {
            entry.topology_id.raw() != 0
                && entry.source.is_none_or(|source| source.data_section == 0)
                && entry
                    .shell_member
                    .is_none_or(|source| source.data_section == 0)
        }));
        assert!(
            serde_json::to_value(&kernel_report)?
                .get("kernel_source_load")
                .is_some()
        );
        Ok(())
    }

    #[test]
    fn kernel_augmentation_skips_recovery_when_nothing_is_unrecovered() {
        let mut calls = 0usize;
        let report = apply_test_kernel_recovery(&[], |_| {
            calls += 1;
            Ok(None)
        });
        assert_eq!(calls, 0);
        assert_eq!(report.schema, KERNEL_COMPLETE_IR_SCHEMA);
        assert_eq!(report.kernel_recovery_failures, Vec::new());
    }

    #[test]
    fn kernel_augmentation_rejects_fragment_for_another_solid() -> Result<()> {
        let mut fragment = Some(extrusion_core_union_fragment(8)?);
        let report = apply_test_kernel_recovery(&[7], |_| Ok(fragment.take()));
        assert_eq!(report.unrecovered_solid_count, 1);
        assert_eq!(report.semantically_recovered_solid_count, 0);
        assert_eq!(report.extrusion_core_unions, Vec::new());
        assert_eq!(report.kernel_recovery_failures.len(), 1);
        assert_eq!(report.kernel_recovery_failures[0].source.id, 7);
        Ok(())
    }

    #[test]
    fn kernel_augmentation_records_backend_failure_per_solid() {
        let report = apply_test_kernel_recovery(&[7], |solid_id| {
            Err(anyhow::anyhow!("kernel refused solid #{solid_id}"))
        });
        assert_eq!(report.unrecovered_solid_count, 1);
        assert_eq!(report.semantically_recovered_solid_count, 0);
        assert_eq!(report.kernel_recovery_failures.len(), 1);
        assert_eq!(report.kernel_recovery_failures[0].source.id, 7);
    }

    #[test]
    fn kernel_augmentation_moves_admitted_union_out_of_unrecovered() -> Result<()> {
        let mut fragment = Some(extrusion_core_union_fragment(7)?);
        let report = apply_test_kernel_recovery(&[7], |_| Ok(fragment.take()));
        assert_eq!(report.unrecovered_solid_count, 0);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.semantically_recovered_solid_count, 1);
        assert_eq!(report.extrusion_core_unions.len(), 1);
        assert_eq!(report.kernel_recovery_failures, Vec::new());
        Ok(())
    }

    #[test]
    fn cheapest_whole_solid_candidate_wins() -> Result<()> {
        let expensive = constructive_fragment(7, 4)?;
        let cheap = constructive_fragment(7, 1)?;
        let selected = select_cheapest_whole_solid_fragments(&[&expensive, &cheap])?;
        assert_eq!(selected.len(), 1);
        assert!(std::ptr::eq(selected[&7], std::ptr::from_ref(&cheap)));
        Ok(())
    }

    #[test]
    fn constructive_candidate_beats_cheaper_descriptive_chain() -> Result<()> {
        let constructive = constructive_fragment(7, 3)?;
        let mut model = CadModel::new();
        let root = model.add_node(CadNode::Primitive(Primitive::Box {
            size_mm: [1.0, 1.0, 1.0],
        }));
        model.add_root(root)?;
        let descriptive = CadFragment {
            source: CadFragmentSource::PeriodicChain {
                solid_id: 7,
                sites: 2,
            },
            model,
            root,
        };
        let selected = select_cheapest_whole_solid_fragments(&[&descriptive, &constructive])?;
        assert_eq!(selected[&7], &constructive);
        Ok(())
    }
}
