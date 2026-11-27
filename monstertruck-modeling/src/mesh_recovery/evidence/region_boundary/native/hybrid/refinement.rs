//! Local source-support refinement for geometrically impossible hybrid seams.
//!
//! A retained analytic carrier can be farther from the middle of a source
//! chord than the patch's entire allowed deviation. Reparameterizing its
//! lifted edge cannot fix that. Release only the adjacent observed support
//! (if explicitly permitted), not the whole independently fitted region.

use super::{HybridInput, HybridRecovery, HybridRecoveryError, HybridSource};
use crate::mesh_recovery::{
    MeshEvidence, NativeTrimOptions, RecoveredAnalyticMesh, RegionBoundaryError,
    RegionBoundaryIndex, RegionPatchOptions, SharedTrimGraph, SharedTrimOptions, SurfaceCarrier,
    SurfaceRegionPartition,
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Original source locus proving an impossible retained carrier constraint.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PatchSupportLocus {
    /// A finite boundary chord cannot fit its permitted continuous deviation.
    Chord([usize; 2]),
    /// A shared source junction cannot fit its permitted displacement.
    Vertex(usize),
}

/// One necessary reason for enlarging an approximate patch.
#[derive(Clone, Debug)]
pub struct PatchSupportConflict {
    /// Original support whose carrier constraint had to be released.
    pub support: usize,
    /// Original fitted region containing that support.
    pub parent_region: usize,
    /// Original source locus with an impossible retained constraint.
    pub locus: PatchSupportLocus,
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
    /// Source-locus witnesses for additional support releases.
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

// Independent vertex and chord budgets induce necessary local release
// dependencies. Index source incidence once per search, rather than scan
// the complete source mesh for every candidate support or vertex.
pub(super) struct ReleaseContext<'a> {
    evidence: &'a MeshEvidence,
    partition: &'a SurfaceRegionPartition,
    permitted: BTreeSet<usize>,
    options: HybridRefinementOptions,
    support_vertices: Vec<Vec<usize>>,
    vertex_conflicts: Vec<Vec<PatchSupportConflict>>,
}

impl<'a> ReleaseContext<'a> {
    pub(super) fn new(
        evidence: &'a MeshEvidence,
        partition: &'a SurfaceRegionPartition,
        permitted: &[usize],
        options: HybridRefinementOptions,
    ) -> Result<Self, HybridRecoveryError> {
        let count = evidence.planar_supports.len();
        if partition.region_by_support.len() != count
            || !options.patch.budget.vertex.is_finite()
            || options.patch.budget.vertex < 0.0
            || !options.patch.budget.chord.is_finite()
            || options.patch.budget.chord <= 0.0
            || !options.patch.intersection_tolerance.is_finite()
            || options.patch.intersection_tolerance <= 0.0
            || permitted.iter().any(|&id| id >= count)
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        let mut support_vertices = Vec::with_capacity(count);
        let mut vertex_conflicts = vec![Vec::new(); evidence.positions.len()];
        for (support, plane) in evidence.planar_supports.iter().enumerate() {
            let parent_region = partition.region_by_support[support];
            let carrier = partition
                .carriers
                .get(parent_region)
                .copied()
                .flatten()
                .ok_or(HybridRecoveryError::InvalidInput)?;
            let mut ids = BTreeSet::new();
            for &polygon in &plane.source_face_indices {
                ids.extend(
                    evidence
                        .source_polygon_vertices(polygon)
                        .ok_or(HybridRecoveryError::InvalidInput)?
                        .iter()
                        .copied(),
                );
            }
            for &vertex in &ids {
                let point = *evidence
                    .positions
                    .get(vertex)
                    .ok_or(HybridRecoveryError::InvalidInput)?;
                if let Some(distance_lower_witness) = vertex_distance_witness(carrier, point)
                    && distance_lower_witness
                        > options.patch.budget.vertex + options.patch.intersection_tolerance
                {
                    vertex_conflicts[vertex].push(PatchSupportConflict {
                        support,
                        parent_region,
                        locus: PatchSupportLocus::Vertex(vertex),
                        distance_lower_witness,
                    });
                }
            }
            support_vertices.push(ids.into_iter().collect());
        }
        Ok(Self {
            evidence,
            partition,
            permitted: permitted.iter().copied().collect(),
            options,
            support_vertices,
            vertex_conflicts,
        })
    }

    fn initial(
        &self,
        seed: HybridRefinementSeed<'_>,
    ) -> Result<BTreeSet<usize>, HybridRecoveryError> {
        let partition = self.partition;
        let selected = match seed {
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
                if supports.iter().any(|&id| id >= self.support_vertices.len()) {
                    return Err(HybridRecoveryError::InvalidInput);
                }
                supports.iter().copied().collect::<BTreeSet<_>>()
            }
        };
        if selected.is_empty() {
            return Err(HybridRecoveryError::InvalidInput);
        }
        Ok(selected)
    }

    fn chord_conflicts(
        &self,
        support: usize,
        selected: &BTreeSet<usize>,
    ) -> Vec<PatchSupportConflict> {
        self.evidence.support_interface_indices[support]
            .iter()
            .filter_map(|&id| {
                let interface = &self.evidence.support_interfaces[id];
                let [a, b] = interface.supports;
                let outside = if a == support { b } else { a };
                if selected.contains(&outside) {
                    return None;
                }
                let parent_region = self.partition.region_by_support[outside];
                let carrier = self.partition.carriers[parent_region]?;
                interface.source_edges.iter().find_map(|&edge| {
                    let [a, b] = edge.map(|id| self.evidence.positions[id]);
                    let distance_lower_witness = carrier_distance_witness(carrier, a, b)?;
                    (distance_lower_witness
                        > self.options.patch.budget.chord
                            + self.options.patch.intersection_tolerance)
                        .then_some(PatchSupportConflict {
                            support: outside,
                            parent_region,
                            locus: PatchSupportLocus::Chord(edge),
                            distance_lower_witness,
                        })
                })
            })
            .collect()
    }

    fn expand(
        &self,
        seed: HybridRefinementSeed<'_>,
        limit: usize,
    ) -> Result<SupportSelection, HybridRecoveryError> {
        let mut selected = self.initial(seed)?;
        let mut pending = selected.iter().copied().collect::<VecDeque<_>>();
        let mut conflicts = Vec::new();
        let mut added_polygons = 0_usize;
        while let Some(support) = pending.pop_front() {
            // Every retained support incident at an exposed source vertex
            // is constrained to the same displacement ball, even if it
            // meets the patch only at a point. Releasing another support
            // can expose new junctions; compute the transitive closure.
            let vertex_conflicts = self.support_vertices[support]
                .iter()
                .flat_map(|&id| self.vertex_conflicts[id].iter().cloned());
            let chord_conflicts = self.chord_conflicts(support, &selected);
            for conflict in vertex_conflicts.chain(chord_conflicts) {
                let outside = conflict.support;
                if selected.contains(&outside) {
                    continue;
                }
                if !self.permitted.contains(&outside) {
                    return Err(HybridRecoveryError::ForbiddenSupport(conflict));
                }
                added_polygons = added_polygons
                    .checked_add(
                        self.evidence.planar_supports[outside]
                            .source_face_indices
                            .len(),
                    )
                    .ok_or(HybridRecoveryError::InvalidInput)?;
                if added_polygons > limit {
                    return Err(HybridRecoveryError::SupportLimit {
                        required: added_polygons,
                        limit,
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

    // Mandatory closure of one seed is contained in the closure of any
    // seed set containing it. The same index is used for cover costs and
    // the actual transaction, so it cannot underestimate mandatory work.
    pub(super) fn closure(&self, support: usize) -> Result<BTreeSet<usize>, HybridRecoveryError> {
        if !self.permitted.contains(&support) {
            return Err(HybridRecoveryError::NoSupportCandidate);
        }
        let cost = self
            .evidence
            .planar_supports
            .get(support)
            .ok_or(HybridRecoveryError::InvalidInput)?
            .source_face_indices
            .len();
        let limit = self.options.max_added_polygons;
        let remaining = limit
            .checked_sub(cost)
            .ok_or(HybridRecoveryError::SupportLimit {
                required: cost,
                limit,
            })?;
        self.expand(HybridRefinementSeed::Supports(&[support]), remaining)
            .map(|selection| selection.selected)
    }
}

fn vertex_distance_witness(carrier: SurfaceCarrier, point: crate::Point3) -> Option<f64> {
    // The normalized plane/meridian residual is a lower bound on distance
    // to the retained carrier, also at point-only source contacts.
    let residual = carrier.signed_constraint(point).map(|(_, r)| r.abs());
    match residual {
        Some(distance) if distance.is_finite() => Some(
            8.0_f64
                .mul_add(-super::junction_roundoff(point, &[carrier]), distance)
                .max(0.0),
        ),
        _ => carrier_distance_witness(carrier, point, point),
    }
}

#[cfg(test)]
fn expand_supports(
    evidence: &MeshEvidence,
    partition: &SurfaceRegionPartition,
    seed: HybridRefinementSeed<'_>,
    permitted: &[usize],
    options: HybridRefinementOptions,
) -> Result<SupportSelection, HybridRecoveryError> {
    ReleaseContext::new(evidence, partition, permitted, options)?
        .expand(seed, options.max_added_polygons)
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
    let labels = original
        .region_by_support
        .iter()
        .enumerate()
        .map(|(support, &parent)| (parent, selected.contains(&support)))
        .collect::<Vec<_>>();
    if labels.len() != evidence.planar_supports.len() {
        return Err(RegionBoundaryError::InvalidPartition);
    }
    let partition = index.repartition_connected(&labels, |(parent, _)| {
        original.carriers.get(parent).copied().flatten()
    })?;
    let mut parents = vec![None; labels.len()];
    let mut selected_regions = Vec::new();
    for (region, carrier) in partition.carriers.iter().enumerate() {
        if carrier.is_some() {
            let (parent, approximate) = labels[region];
            parents[region] = Some(parent);
            if approximate {
                selected_regions.push(region);
            }
        }
    }
    Ok(RefinedPartition {
        partition,
        parents,
        selected_regions,
    })
}

fn check_original_retention(
    original: &HybridSource<'_>,
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
        HybridSource::primary(self).construct_refined_hybrid(seed, permitted_supports, options)
    }
}

impl HybridSource<'_> {
    pub(super) fn construct_refined_hybrid(
        &self,
        seed: HybridRefinementSeed<'_>,
        permitted_supports: &[usize],
        options: HybridRefinementOptions,
    ) -> Result<RefinedHybridRecovery, HybridRecoveryError> {
        let context = ReleaseContext::new(
            self.evidence,
            &self.reconciled.partition,
            permitted_supports,
            options,
        )?;
        self.construct_refined_with_context(seed, &context, options.max_added_polygons)
    }

    pub(super) fn construct_refined_with_context(
        &self,
        seed: HybridRefinementSeed<'_>,
        context: &ReleaseContext<'_>,
        limit: usize,
    ) -> Result<RefinedHybridRecovery, HybridRecoveryError> {
        let options = context.options;
        let index = RegionBoundaryIndex::new(self.evidence)?;
        let original = &self.reconciled.partition;
        let selection = context.expand(seed, limit)?;
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
        if whole_regions {
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
                added_polygons: selection.added_polygons,
                trims: self.trims.clone(),
                recovery,
            });
        }
        let refined = split_partition(&index, self.evidence, original, &selection.selected)?;
        let mut graph = index.assemble_shared_trims(
            &refined.partition,
            &self.reconciled.positions,
            options.trim,
        )?;
        let mut baseline = index.construct_native_trims(
            &refined.partition,
            &graph,
            &self.reconciled.positions,
            options.strict,
        )?;
        index.repair_planar_crossings(
            &refined.partition,
            &self.reconciled.positions,
            &mut graph,
            &mut baseline,
            options.strict,
        )?;
        let recovery = HybridInput {
            evidence: self.evidence,
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
    fn vertex_budget_conflicts_propagate_independently_of_chord_budget()
    -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, mut partition, fallback, sphere) = spherical_fixture()?;
        let Some(SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        }) = partition.carriers[sphere]
        else {
            return Err("expected sphere carrier".into());
        };
        partition.carriers[sphere] = Some(SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius: profile_radius + 1e-4,
        });
        let permitted = partition
            .region_by_support
            .iter()
            .enumerate()
            .filter_map(|(id, &region)| (region == sphere).then_some(id))
            .collect::<Vec<_>>();
        let policy = options(0.5, 4);
        assert!(matches!(
            expand_supports(
                &evidence,
                &partition,
                HybridRefinementSeed::Regions(&fallback),
                &[],
                policy
            ),
            Err(HybridRecoveryError::ForbiddenSupport(
                PatchSupportConflict {
                    locus: PatchSupportLocus::Vertex(_),
                    ..
                }
            ))
        ));
        let context = ReleaseContext::new(&evidence, &partition, &permitted, policy)?;
        let expanded = context.expand(HybridRefinementSeed::Regions(&fallback), 4)?;
        assert_eq!(expanded.added_polygons, 4);
        assert!(
            expanded
                .conflicts
                .iter()
                .all(|c| matches!(c.locus, PatchSupportLocus::Vertex(_))
                    && c.distance_lower_witness > policy.patch.budget.vertex)
        );
        let unauthorized = partition
            .region_by_support
            .iter()
            .position(|&region| region != sphere)
            .ok_or("fallback support")?;
        assert!(matches!(
            context.closure(unauthorized),
            Err(HybridRecoveryError::NoSupportCandidate)
        ));
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
