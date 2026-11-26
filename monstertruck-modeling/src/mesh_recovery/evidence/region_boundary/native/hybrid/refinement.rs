//! Local source-support refinement for geometrically impossible hybrid seams.
//!
//! A retained analytic carrier can be farther from the middle of a source
//! chord than the patch's entire allowed deviation. Reparameterizing its
//! lifted edge cannot fix that. Release only the adjacent observed support
//! (if explicitly permitted), not the whole independently fitted region.

use super::{HybridInput, HybridRecovery, HybridRecoveryError};
use crate::mesh_recovery::{
    MeshEvidence, NativeTrimOptions, RecoveredAnalyticMesh, RegionBoundaryError,
    RegionBoundaryIndex, RegionPatchOptions, SharedTrimGraph, SharedTrimOptions, SurfaceCarrier,
    SurfaceRegionPartition,
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// One necessary reason for enlarging an approximate patch.
#[derive(Clone, Debug)]
pub struct PatchSupportConflict {
    /// Original support whose carrier constraint had to be released.
    pub support: usize,
    /// Original fitted region containing that support.
    pub parent_region: usize,
    /// Original source edge whose chord cannot be represented on that carrier.
    pub edge: [usize; 2],
    /// Numerical lower witness for the maximum distance to the retained
    /// carrier, in source units, after a roundoff allowance.
    pub distance_lower_witness: f64,
}

/// Explicit policy for local support release; all geometric budgets are
/// identical to the ordinary hybrid transaction.
#[derive(Clone, Copy, Debug)]
pub struct HybridRefinementOptions {
    /// Independent approximation policy.
    pub patch: RegionPatchOptions,
    /// Native construction policy of the original strict reconstruction.
    pub strict: NativeTrimOptions,
    /// Original strict shared-curve fitting policy.
    pub trim: SharedTrimOptions,
    /// Maximum additional source polygons authorized for release.
    pub max_added_polygons: usize,
}

/// Explicit original ownership to approximate before witnessed expansion.
#[derive(Clone, Copy, Debug)]
pub enum HybridRefinementSeed<'a> {
    /// Release every support of the listed original fitted regions.
    Regions(&'a [usize]),
    /// Release only the listed original planar source supports, splitting
    /// their fitted regions when necessary.
    Supports(&'a [usize]),
}

/// The original regions are split, never refitted. New region IDs refer to
/// this partition, not the original discovery graph.
#[derive(Debug)]
pub struct RefinedHybridRecovery {
    /// Complete source ownership, including retained carrier fragments.
    pub partition: SurfaceRegionPartition,
    /// Original region of each active refined region.
    pub parent_by_region: Vec<Option<usize>>,
    /// Refined regions whose carrier constraints were relinquished.
    pub selected_regions: Vec<usize>,
    /// Source-edge witnesses for additional support releases.
    pub conflicts: Vec<PatchSupportConflict>,
    /// Number of extra original source polygons replaced.
    pub added_polygons: usize,
    /// Refined analytic graph and complete original source-edge provenance.
    /// Native lifted trims in the result replace the corresponding graph
    /// geometry; their IDs remain aligned with this graph.
    pub trims: SharedTrimGraph,
    /// Globally shared native result in refined region/trim order.
    pub recovery: HybridRecovery,
}

struct SupportSelection {
    selected: BTreeSet<usize>,
    conflicts: Vec<PatchSupportConflict>,
    added_polygons: usize,
}

fn carrier_distance_witness(
    carrier: SurfaceCarrier,
    a: crate::Point3,
    b: crate::Point3,
) -> Option<f64> {
    // Cylindrical coordinates (z, rho) are nonexpansive. Distance from
    // their image to the full meridian line/circle lower-bounds distance
    // to the represented carrier, including any restricted radial sheet.
    let witness = carrier.measure_segment(a, b)?.witnessed_max;
    let normalized = match carrier {
        SurfaceCarrier::LinearRevolution { slope, .. } => witness / 1.0_f64.hypot(slope),
        _ => witness,
    };
    let margin =
        super::junction_roundoff(a, &[carrier]).max(super::junction_roundoff(b, &[carrier])) * 8.0;
    (normalized.is_finite() && margin.is_finite()).then_some((normalized - margin).max(0.0))
}

fn expand_supports(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    seed: HybridRefinementSeed<'_>,
    permitted: &[usize],
    options: HybridRefinementOptions,
) -> Result<SupportSelection, HybridRecoveryError> {
    let count = evidence.planar_supports.len();
    if partition.region_by_support.len() != count
        || !options.patch.budget.chord.is_finite()
        || options.patch.budget.chord <= 0.0
        || !options.patch.intersection_tolerance.is_finite()
        || options.patch.intersection_tolerance <= 0.0
        || permitted.iter().any(|&id| id >= count)
    {
        return Err(HybridRecoveryError::InvalidInput);
    }
    let permitted = permitted.iter().copied().collect::<BTreeSet<_>>();
    let mut selected = match seed {
        HybridRefinementSeed::Regions(regions) => {
            if regions
                .iter()
                .any(|&id| partition.carriers.get(id).is_none_or(Option::is_none))
            {
                return Err(HybridRecoveryError::InvalidInput);
            }
            let required = regions.iter().copied().collect::<BTreeSet<_>>();
            partition
                .region_by_support
                .iter()
                .enumerate()
                .filter_map(|(support, region)| required.contains(region).then_some(support))
                .collect::<BTreeSet<_>>()
        }
        HybridRefinementSeed::Supports(supports) => {
            if supports.iter().any(|&id| {
                partition
                    .region_by_support
                    .get(id)
                    .and_then(|&r| partition.carriers.get(r))
                    .is_none_or(Option::is_none)
            }) {
                return Err(HybridRecoveryError::InvalidInput);
            }
            supports.iter().copied().collect::<BTreeSet<_>>()
        }
    };
    if selected.is_empty() {
        return Err(HybridRecoveryError::InvalidInput);
    }
    let mut pending = selected.iter().copied().collect::<VecDeque<_>>();
    let mut conflicts = Vec::new();
    let mut added_polygons = 0_usize;
    while let Some(support) = pending.pop_front() {
        for &interface_id in &evidence.support_interface_indices[support] {
            let interface = &evidence.support_interfaces[interface_id];
            let [a, b] = interface.supports;
            let outside = if a == support { b } else { a };
            if selected.contains(&outside) {
                continue;
            }
            let parent_region = partition.region_by_support[outside];
            let carrier =
                partition.carriers[parent_region].ok_or(HybridRecoveryError::InvalidInput)?;
            let conflict = interface.source_edges.iter().find_map(|&edge| {
                let [a, b] = edge.map(|id| evidence.positions[id]);
                let distance_lower_witness = carrier_distance_witness(carrier, a, b)?;
                // Keep a further numerical incidence allowance before
                // calling the obstruction unavoidable.
                (distance_lower_witness
                    > options.patch.budget.chord + options.patch.intersection_tolerance)
                    .then_some(PatchSupportConflict {
                        support: outside,
                        parent_region,
                        edge,
                        distance_lower_witness,
                    })
            });
            let Some(conflict) = conflict else {
                continue;
            };
            if !permitted.contains(&outside) {
                return Err(HybridRecoveryError::ForbiddenSupport(conflict));
            }
            added_polygons = added_polygons
                .checked_add(evidence.planar_supports[outside].source_face_indices.len())
                .ok_or(HybridRecoveryError::InvalidInput)?;
            if added_polygons > options.max_added_polygons {
                return Err(HybridRecoveryError::SupportLimit {
                    required: added_polygons,
                    limit: options.max_added_polygons,
                });
            }
            selected.insert(outside);
            conflicts.push(conflict);
            pending.push_back(outside);
        }
    }
    Ok(SupportSelection {
        selected,
        conflicts,
        added_polygons,
    })
}

struct RefinedPartition {
    partition: SurfaceRegionPartition,
    parents: Vec<Option<usize>>,
    selected_regions: Vec<usize>,
}

fn split_partition(
    index: &RegionBoundaryIndex<'_>,
    evidence: &MeshEvidence,
    original: &SurfaceRegionPartition,
    selected: &BTreeSet<usize>,
) -> Result<RefinedPartition, RegionBoundaryError> {
    let count = evidence.planar_supports.len();
    let mut region_by_support = vec![usize::MAX; count];
    let mut carriers = vec![None; count];
    let mut parents = vec![None; count];
    let mut selected_regions = Vec::new();
    for seed in 0..count {
        if region_by_support[seed] != usize::MAX {
            continue;
        }
        let parent = original.region_by_support[seed];
        let carrier = original.carriers[parent].ok_or(RegionBoundaryError::InvalidPartition)?;
        let approximate = selected.contains(&seed);
        carriers[seed] = Some(carrier);
        parents[seed] = Some(parent);
        if approximate {
            selected_regions.push(seed);
        }
        region_by_support[seed] = seed;
        let mut queue = VecDeque::from([seed]);
        while let Some(support) = queue.pop_front() {
            for &id in &evidence.support_interface_indices[support] {
                let [a, b] = evidence.support_interfaces[id].supports;
                let other = if a == support { b } else { a };
                if region_by_support[other] == usize::MAX
                    && original.region_by_support[other] == parent
                    && selected.contains(&other) == approximate
                {
                    region_by_support[other] = seed;
                    queue.push_back(other);
                }
            }
        }
    }
    let region_by_polygon = evidence
        .source_planar_support
        .iter()
        .map(|&support| region_by_support[support])
        .collect::<Vec<_>>();
    let contacts = index.partition_contacts(&region_by_polygon)?;
    let consolidated_supports = count - carriers.iter().flatten().count();
    Ok(RefinedPartition {
        partition: SurfaceRegionPartition {
            // This is a refinement of existing ownership, not a second
            // proposal selection. Parent IDs carry the discovery provenance.
            selected_proposals: Vec::new(),
            region_by_support,
            region_by_polygon,
            carriers,
            contacts,
            consolidated_supports,
            exact_components: 0,
            limited_components: 0,
        },
        parents,
        selected_regions,
    })
}

fn check_original_retention(
    original: &RecoveredAnalyticMesh,
    refined: &RefinedPartition,
    graph: &crate::mesh_recovery::SharedTrimGraph,
    result: &HybridRecovery,
) -> Result<(), HybridRecoveryError> {
    let selected = refined
        .selected_regions
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let old = original
        .native
        .regions
        .iter()
        .map(|region| (region.region, region))
        .collect::<BTreeMap<_, _>>();
    for region in &result.native.regions {
        if selected.contains(&region.region) {
            continue;
        }
        let parent = refined.parents[region.region].ok_or(HybridRecoveryError::InvalidInput)?;
        let before = old[&parent];
        if before.faces.is_ok()
            && (region.faces.is_err() || (!before.requires_chart_cut && region.requires_chart_cut))
        {
            return Err(HybridRecoveryError::RetainedRegion(region.region));
        }
    }
    // Original valid retained-retained contacts cannot disappear or become
    // unresolved merely because a new support seam changes the segmentation.
    let key = |[a, b]: [usize; 2]| [a.min(b), a.max(b)];
    let mut valid = BTreeSet::new();
    for (trim, native) in graph.trims.iter().zip(&result.native.trims) {
        if native.is_ok() {
            valid.extend(trim.source_edges.iter().map(|edge| key(edge.vertices)));
        }
    }
    for (id, (trim, native)) in original
        .trims
        .trims
        .iter()
        .zip(&original.native.trims)
        .enumerate()
    {
        if native.is_err() {
            continue;
        }
        for edge in &trim.source_edges {
            let inside = refined.partition.region_by_polygon[edge.interior_polygon];
            let outside = edge
                .exterior_polygon
                .map(|p| refined.partition.region_by_polygon[p]);
            if selected.contains(&inside) || outside.is_some_and(|r| selected.contains(&r)) {
                continue;
            }
            if !valid.contains(&key(edge.vertices)) {
                return Err(HybridRecoveryError::RetainedTrim {
                    trim: id,
                    error: crate::mesh_recovery::NativeTrimError::InvalidSource,
                });
            }
        }
    }
    Ok(())
}

impl RecoveredAnalyticMesh {
    /// Extend an explicit approximation request across only those source
    /// support interfaces with a witnessed impossible retained-carrier
    /// distance. Split each affected original region into connected
    /// selected and retained fragments without refitting its carrier.
    ///
    /// Additional support IDs must be explicitly permitted and their total
    /// source polygon count must fit the separate authorization budget.
    /// All geometric approximation bounds and original evidence are unchanged.
    /// A small distance witness is not proof that a lift will succeed.
    ///
    /// # Errors
    /// Rejects unauthorized releases, excessive expansion, unresolved
    /// geometry, and regressions of previously valid retained regions/trims.
    pub fn construct_refined_hybrid(
        &self,
        seed: HybridRefinementSeed<'_>,
        permitted_supports: &[usize],
        options: HybridRefinementOptions,
    ) -> Result<RefinedHybridRecovery, HybridRecoveryError> {
        let index = RegionBoundaryIndex::new(&self.evidence)?;
        let original = &self.reconciled.partition;
        let selection =
            expand_supports(&self.evidence, original, seed, permitted_supports, options)?;
        let selected_parents = selection
            .selected
            .iter()
            .map(|&support| original.region_by_support[support])
            .collect::<BTreeSet<_>>();
        let whole_regions =
            original
                .region_by_support
                .iter()
                .enumerate()
                .all(|(support, region)| {
                    !selected_parents.contains(region) || selection.selected.contains(&support)
                });
        if selection.added_polygons == 0 && whole_regions {
            // No support ownership changed. Keep the exact existing
            // recognition/trim graph and native baseline; there is no
            // reason to fit all unchanged contacts again.
            let selected_regions = selected_parents.into_iter().collect::<Vec<_>>();
            let recovery =
                self.construct_hybrid(&selected_regions, options.patch, options.strict)?;
            return Ok(RefinedHybridRecovery {
                partition: original.clone(),
                parent_by_region: original
                    .carriers
                    .iter()
                    .enumerate()
                    .map(|(id, c)| c.as_ref().map(|_| id))
                    .collect(),
                selected_regions,
                conflicts: selection.conflicts,
                added_polygons: 0,
                trims: self.trims.clone(),
                recovery,
            });
        }
        let refined = split_partition(&index, &self.evidence, original, &selection.selected)?;
        let graph = index.assemble_shared_trims(
            &refined.partition,
            &self.reconciled.positions,
            options.trim,
        )?;
        let baseline = index.construct_native_trims(
            &refined.partition,
            &graph,
            &self.reconciled.positions,
            options.strict,
        )?;
        let recovery = HybridInput {
            evidence: &self.evidence,
            partition: &refined.partition,
            graph: &graph,
            baseline: &baseline,
            positions: &self.reconciled.positions,
            selected: refined.selected_regions.iter().copied().collect(),
            patch: options.patch,
            strict: options.strict,
        }
        .assemble()?;
        check_original_retention(self, &refined, &graph, &recovery)?;
        Ok(RefinedHybridRecovery {
            partition: refined.partition,
            parent_by_region: refined.parents,
            selected_regions: refined.selected_regions,
            conflicts: selection.conflicts,
            added_polygons: selection.added_polygons,
            trims: graph,
            recovery,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::super::patch_lift::tests::spherical_fixture;
    use super::*;
    use crate::{Point3, Vector3, mesh_recovery::SurfaceRegionBudget};

    fn options(chord: f64, limit: usize) -> HybridRefinementOptions {
        HybridRefinementOptions {
            patch: RegionPatchOptions {
                budget: SurfaceRegionBudget {
                    vertex: 1e-8,
                    chord,
                },
                intersection_tolerance: 1e-10,
                angular_resolution: 1e-8,
            },
            strict: NativeTrimOptions {
                max_vertex_displacement: 1e-7,
                intersection_tolerance: 1e-10,
                angular_resolution: 1e-8,
            },
            trim: SharedTrimOptions {
                absolute_tolerance: 1e-7,
                relative_tolerance: 1e-5,
                max_interval_edges: 100,
            },
            max_added_polygons: limit,
        }
    }

    #[test]
    fn cone_witness_is_normalized_spatial_distance() {
        let cone = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radius_at_origin: 1.0,
            slope: 1.0,
        };
        let a = Point3::new(1.0, 0.0, 0.0);
        let b = Point3::new(0.0, 1.0, 0.0);
        let expected = (1.0 - std::f64::consts::FRAC_1_SQRT_2) * std::f64::consts::FRAC_1_SQRT_2;
        let witness = carrier_distance_witness(cone, a, b).expect("finite meridian");
        assert!(witness <= expected + 1e-14);
        assert!(expected - witness < 1e-10);
        assert!(
            carrier_distance_witness(cone, a, Point3::new(2.0, 0.0, 1.0)).expect("generator")
                < 1e-12
        );
    }

    #[test]
    fn support_expansion_requires_permission_and_polygon_budget()
    -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, partition, fallback, sphere) = spherical_fixture()?;
        let permitted = partition
            .region_by_support
            .iter()
            .enumerate()
            .filter_map(|(id, &r)| (r == sphere).then_some(id))
            .collect::<Vec<_>>();
        let unchanged = expand_supports(
            &evidence,
            &partition,
            HybridRefinementSeed::Regions(&fallback),
            &[],
            options(0.1, 0),
        )?;
        assert_eq!(unchanged.added_polygons, 0);
        assert!(matches!(
            expand_supports(
                &evidence,
                &partition,
                HybridRefinementSeed::Regions(&fallback),
                &[],
                options(0.02, 4)
            ),
            Err(HybridRecoveryError::ForbiddenSupport(_))
        ));
        assert!(matches!(
            expand_supports(
                &evidence,
                &partition,
                HybridRefinementSeed::Regions(&fallback),
                &permitted,
                options(0.02, 0)
            ),
            Err(HybridRecoveryError::SupportLimit { .. })
        ));
        let expanded = expand_supports(
            &evidence,
            &partition,
            HybridRefinementSeed::Regions(&fallback),
            &permitted,
            options(0.02, 4),
        )?;
        assert_eq!(expanded.added_polygons, 4);
        assert!(
            expanded
                .conflicts
                .iter()
                .all(|c| c.parent_region == sphere && c.distance_lower_witness > 0.02)
        );
        assert_eq!(expanded.selected.len(), evidence.planar_supports.len());
        Ok(())
    }

    #[test]
    fn explicit_support_seed_does_not_release_entire_parent()
    -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, partition, _, sphere) = spherical_fixture()?;
        let support = evidence.source_planar_support[0];
        assert_eq!(partition.region_by_support[support], sphere);
        let selected = expand_supports(
            &evidence,
            &partition,
            HybridRefinementSeed::Supports(&[support]),
            &[],
            options(0.5, 0),
        )?;
        assert_eq!(selected.selected, BTreeSet::from([support]));
        assert_eq!(selected.added_polygons, 0);
        assert!(matches!(
            expand_supports(
                &evidence,
                &partition,
                HybridRefinementSeed::Supports(&[evidence.planar_supports.len()]),
                &[],
                options(0.5, 0)
            ),
            Err(HybridRecoveryError::InvalidInput)
        ));
        Ok(())
    }

    #[test]
    fn partition_refinement_retains_other_supports_on_their_original_carrier()
    -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, partition, _, sphere) = spherical_fixture()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let selected = BTreeSet::from([evidence.source_planar_support[0]]);
        let refined = split_partition(&index, &evidence, &partition, &selected)?;
        assert_eq!(refined.selected_regions.len(), 1);
        let release = refined.selected_regions[0];
        assert_eq!(refined.parents[release], Some(sphere));
        for (id, &support) in evidence.source_planar_support.iter().enumerate() {
            let region = refined.partition.region_by_polygon[id];
            assert_eq!(
                refined.parents[region],
                Some(partition.region_by_polygon[id])
            );
            assert_eq!(region == release, selected.contains(&support));
        }
        assert_eq!(
            refined.partition.contacts,
            index.partition_contacts(&refined.partition.region_by_polygon)?
        );
        assert!(
            refined
                .partition
                .carriers
                .iter()
                .enumerate()
                .filter(|(r, _)| refined.parents[*r] == Some(sphere))
                .all(|(_, c)| matches!(c, Some(SurfaceCarrier::CircularRevolution { .. })))
        );
        Ok(())
    }
}
