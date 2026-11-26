//! Global hybrid reconstruction using one source-junction and trim table.
//!
//! Explicitly selected approximate regions relinquish their fitted carrier
//! constraints. All retained analytic contacts keep their strict geometry.
//! Approximate contact paths use one exact retained-carrier edge per
//! original source segment; both incident faces consume the same identity.

use super::super::patch_lift::{junction_roundoff, lift_native_source_edge};
use super::{
    NativeRegionFace, NativeTrimError, NativeTrimGraph, NativeTrimOptions, RegionAssemblyInput,
    SectionBranch, SharedTrimGeometry, SharedTrimGraph, SurfaceCarrier, SurfaceRegionPartition,
    TrimInterval, build_trim, geometry_constraints, junction_candidates, trim_intervals,
};
use crate::{
    Edge, Face, InnerSpace, Point3, Shell, Vertex, builder,
    mesh_recovery::{
        BoundaryPatchError, BoundaryPatchOptions, MeshEvidence, PrescribedSourceEdge,
        RecoveredAnalyticMesh, RegionBoundaryError, RegionPatchOptions, approximate_boundary_patch,
    },
};
use monstertruck_topology::shell::ShellCondition;
use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Range,
};

mod connectors;
pub use connectors::{HybridConnectorOptions, HybridConnectorRecovery, HybridConnectorRelease};
mod refinement;
mod simplify;
pub use refinement::{HybridRefinementOptions, PatchSupportConflict, RefinedHybridRecovery};
pub use simplify::HybridPlanarSimplification;

/// Provenance of one original polygon in its approximate native region.
#[derive(Debug, Clone)]
pub struct HybridPolygonFaces {
    /// Active original source polygon ID.
    pub polygon: usize,
    /// Original selected region ID.
    pub region: usize,
    /// Native face indices within that region's face list.
    /// A face can cover multiple source polygons after cell coalescence
    /// or certified planar simplification.
    pub faces: Vec<usize>,
}

/// Evidence for one connected, bounded approximate region union.
#[derive(Debug)]
pub struct HybridPatchEvidence {
    budget: crate::mesh_recovery::SurfaceRegionBudget,
    /// Regions explicitly replaced by source-backed approximations.
    pub regions: Vec<usize>,
    /// Source polygon-to-native-face provenance.
    pub polygons: Vec<HybridPolygonFaces>,
    /// Number of certified cells covering multiple original polygons.
    pub coalesced_cells: usize,
    /// Bound on distance of an original source vertex to the native
    /// domain, including any certified planar simplification.
    pub vertex_displacement: f64,
    /// Maximum continuous two-sided facet deviation.
    pub deviation_upper_bound: f64,
    /// Maximum bounded radial seam discrepancy within the patch.
    pub seam_deviation_upper_bound: f64,
}

impl HybridPatchEvidence {
    /// The original immutable vertex and chord budgets of this patch.
    #[must_use]
    pub const fn budget(&self) -> crate::mesh_recovery::SurfaceRegionBudget { self.budget }
}

/// One global hybrid reconstruction. The input strict result is unchanged.
#[derive(Debug)]
pub struct HybridRecovery {
    /// Shared native edge and face identities in original graph order.
    pub native: NativeTrimGraph,
    /// Certified approximate region unions.
    pub patches: Vec<HybridPatchEvidence>,
    /// Original trim IDs replaced by bounded single-carrier source paths.
    /// The value is the maximum deviation from any original source chord.
    pub lifted_trims: BTreeMap<usize, f64>,
    ambiguous_edges: bool,
}

impl HybridRecovery {
    /// Extract a complete oriented, nonsingular, geometrically consistent
    /// shell. This is not an occupied-volume or self-intersection proof.
    #[must_use]
    pub fn complete_shell(&self) -> Option<Shell> {
        if self.ambiguous_edges || self.native.trims.iter().any(Result::is_err) {
            return None;
        }
        let mut shell = Shell::new();
        for region in &self.native.regions {
            if region.requires_chart_cut {
                return None;
            }
            shell.extend(region.faces.as_ref().ok()?.iter().cloned());
        }
        (!shell.is_empty()
            && matches!(
                shell.shell_condition(),
                ShellCondition::Closed | ShellCondition::Oriented
            )
            && shell.singular_vertices().is_empty()
            && shell.is_geometric_consistent())
        .then_some(shell)
    }
}

/// Failure of an explicit global hybrid transaction.
#[derive(Debug, thiserror::Error)]
pub enum HybridRecoveryError {
    /// Inconsistent source references or numerical policy.
    #[error("invalid hybrid recovery inputs")]
    InvalidInput,
    /// Source region selection or topology is invalid.
    #[error(transparent)]
    Boundary(#[from] RegionBoundaryError),
    /// A global junction or trim could not be constructed.
    #[error(transparent)]
    Native(#[from] NativeTrimError),
    /// A selected source patch could not be approximated.
    #[error("hybrid patch starting at region {region}: {error}")]
    Patch {
        /// First original region ID of the connected selected union.
        region: usize,
        /// Failed approximation certificate.
        error: BoundaryPatchError,
    },
    /// A source segment could not be lifted within the selected budgets.
    #[error("hybrid trim {trim} source edge {vertices:?} could not be lifted")]
    LiftedEdge {
        /// Original shared trim ID.
        trim: usize,
        /// Original source endpoint IDs.
        vertices: [usize; 2],
    },
    /// An unavoidable seam conflict would require an unauthorized support release.
    #[error("hybrid support release is not permitted: {0:?}")]
    ForbiddenSupport(PatchSupportConflict),
    /// Additional original polygon ownership would exceed the authorized limit.
    #[error("hybrid support release requires {required} additional polygons, limit {limit}")]
    SupportLimit {
        /// Minimum number found so far.
        required: usize,
        /// Authorized maximum.
        limit: usize,
    },
    /// No admissible transaction within the connector search budgets.
    #[error("hybrid connector search exhausted after {trials} trials: {last}")]
    ConnectorSearch {
        /// Number of attempted native transactions.
        trials: usize,
        /// Most recent geometric rejection.
        last: Box<Self>,
    },
    /// A previously valid retained analytic trim would fail.
    #[error("hybrid recovery would invalidate retained trim {trim}: {error}")]
    RetainedTrim {
        /// Source graph trim ID.
        trim: usize,
        /// Construction or continuous-domain rejection.
        error: NativeTrimError,
    },
    /// A previously valid retained analytic region would fail or require
    /// a previously unnecessary chart decomposition.
    #[error("hybrid recovery would invalidate retained region {0}")]
    RetainedRegion(usize),
}

#[derive(Clone, Copy)]
enum TrimPolicy {
    Retained,
    Lifted(Option<SurfaceCarrier>),
}

struct HybridPlan {
    graph: SharedTrimGraph,
    policy: Vec<TrimPolicy>,
    intervals: Vec<Result<Vec<TrimInterval>, NativeTrimError>>,
    touched: BTreeSet<usize>,
}

struct HybridEdges {
    trims: Vec<Result<Vec<Edge>, NativeTrimError>>,
    bounds: BTreeMap<usize, f64>,
}

struct JunctionConstraints {
    carriers: Vec<SurfaceCarrier>,
    sections: Vec<SectionBranch>,
}

struct HybridInput<'a> {
    evidence: &'a MeshEvidence,
    partition: &'a SurfaceRegionPartition,
    graph: &'a SharedTrimGraph,
    baseline: &'a NativeTrimGraph,
    positions: &'a [Point3],
    selected: BTreeSet<usize>,
    patch: RegionPatchOptions,
    strict: NativeTrimOptions,
}

impl RecoveredAnalyticMesh {
    /// Rebuild selected regions as bounded source-backed approximations
    /// in one shared native graph. No retained carrier is replaced.
    ///
    /// All incident retained constraints are reconciled at newly exposed
    /// source junctions; retained-retained trims preserve their original
    /// analytic geometry. The native options must be the ones used for
    /// the strict result. The separate patch budget applies to original
    /// source coordinates, not already displaced native junctions.
    ///
    /// # Errors
    /// Rejects an uncertified selected patch or any regression of a
    /// previously successful retained trim or region. Other existing
    /// strict recovery failures may remain in the partial result.
    pub fn construct_hybrid(
        &self,
        selected: &[usize],
        patch: RegionPatchOptions,
        strict: NativeTrimOptions,
    ) -> Result<HybridRecovery, HybridRecoveryError> {
        HybridInput {
            evidence: &self.evidence,
            partition: &self.reconciled.partition,
            graph: &self.trims,
            baseline: &self.native,
            positions: &self.reconciled.positions,
            selected: selected.iter().copied().collect(),
            patch,
            strict,
        }
        .assemble()
    }
}

impl HybridInput<'_> {
    fn validate(&self) -> Result<(), HybridRecoveryError> {
        if self.selected.is_empty()
            || self
                .selected
                .iter()
                .any(|&id| self.partition.carriers.get(id).is_none_or(Option::is_none))
            || self.positions.len() != self.evidence.positions.len()
            || self.baseline.vertices.len() != self.positions.len()
            || self.graph.trims.len() != self.baseline.trims.len()
            || self.partition.region_by_polygon.len() != self.evidence.source_planar_support.len()
            || self.partition.region_by_support.len() != self.evidence.planar_supports.len()
            || !self.patch.budget.vertex.is_finite()
            || self.patch.budget.vertex < 0.0
            || [
                self.patch.budget.chord,
                self.patch.intersection_tolerance,
                self.patch.angular_resolution,
                self.strict.max_vertex_displacement,
                self.strict.intersection_tolerance,
                self.strict.angular_resolution,
            ]
            .iter()
            .any(|&v| !v.is_finite() || v <= 0.0)
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        if self
            .evidence
            .source_planar_support
            .iter()
            .zip(&self.partition.region_by_polygon)
            .any(|(&support, &region)| {
                self.partition.region_by_support.get(support) != Some(&region)
                    || self
                        .partition
                        .carriers
                        .get(region)
                        .is_none_or(Option::is_none)
            })
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        Ok(())
    }

    fn plan(&self) -> Result<HybridPlan, HybridRecoveryError> {
        let mut graph = self.graph.clone();
        let mut policy = Vec::with_capacity(graph.trims.len());
        let mut intervals = Vec::with_capacity(graph.trims.len());
        let mut touched = BTreeSet::new();
        for trim in &mut graph.trims {
            let sides = std::iter::once(trim.regions.0).chain(trim.regions.1);
            let (approximate, retained): (Vec<_>, Vec<_>) =
                sides.partition(|id| self.selected.contains(id));
            if approximate.is_empty() {
                policy.push(TrimPolicy::Retained);
                intervals.push(trim_intervals(trim, self.positions));
                continue;
            }
            let carrier = retained
                .first()
                .map(|&id| {
                    self.partition
                        .carriers
                        .get(id)
                        .copied()
                        .flatten()
                        .ok_or(HybridRecoveryError::InvalidInput)
                })
                .transpose()?;
            policy.push(TrimPolicy::Lifted(carrier));
            // A former line/circle may now follow a different exact
            // single-carrier path; do not feed stale ring evidence to
            // periodic chart planning.
            trim.geometry = SharedTrimGeometry::Unresolved;
            let mut pieces = Vec::with_capacity(trim.source_edges.len());
            for edge in &trim.source_edges {
                if edge.vertices.iter().any(|&id| id >= self.positions.len()) {
                    return Err(HybridRecoveryError::InvalidInput);
                }
                touched.extend(edge.vertices);
                pieces.push(TrimInterval {
                    endpoints: edge.vertices,
                    sweep: 0.0,
                });
            }
            if pieces.is_empty() {
                return Err(HybridRecoveryError::InvalidInput);
            }
            intervals.push(Ok(pieces));
        }
        Ok(HybridPlan {
            graph,
            policy,
            intervals,
            touched,
        })
    }

    fn assemble(&self) -> Result<HybridRecovery, HybridRecoveryError> {
        self.validate()?;
        let plan = self.plan()?;
        let vertices = self.vertices(&plan)?;
        let HybridEdges {
            mut trims,
            bounds: lifted_trims,
        } = self.edges(&plan, &vertices)?;
        let mut points = self.positions.to_vec();
        for &id in &plan.touched {
            points[id] = vertices[id]
                .as_ref()
                .ok_or(NativeTrimError::UnresolvedVertex(id))?
                .point();
        }
        // Bounds are against original chords; the chart tests use the
        // displaced source endpoints. Add the maximum endpoint movement
        // to translate the bound to those chords.
        let deviations = (0..trims.len())
            .map(|id| {
                lifted_trims
                    .get(&id)
                    .map(|bound| bound + self.patch.budget.vertex)
            })
            .collect::<Vec<_>>();
        let mut regions = RegionAssemblyInput {
            evidence: self.evidence,
            partition: self.partition,
            graph: &plan.graph,
            positions: &points,
            options: self.strict,
            approximated: &self.selected,
            lifted_deviations: &deviations,
        }
        .assemble(&mut trims)?;
        self.check_retained(&plan, &trims, &regions)?;
        let patches = self.assemble_patches(&plan, &trims, &mut regions)?;
        regions.sort_unstable_by_key(|r| r.region);
        let max_vertex_displacement = vertices
            .iter()
            .zip(&self.evidence.positions)
            .filter_map(|(v, p)| v.as_ref().map(|v| (v.point() - p).magnitude()))
            .fold(0.0_f64, f64::max);
        Ok(HybridRecovery {
            native: NativeTrimGraph {
                vertices,
                trims,
                regions,
                max_vertex_displacement,
            },
            patches,
            lifted_trims,
            ambiguous_edges: !self.graph.ambiguous_edges.is_empty(),
        })
    }

    fn check_retained(
        &self,
        plan: &HybridPlan,
        trims: &[Result<Vec<Edge>, NativeTrimError>],
        regions: &[NativeRegionFace],
    ) -> Result<(), HybridRecoveryError> {
        for (id, (old, new)) in self.baseline.trims.iter().zip(trims).enumerate() {
            if matches!(plan.policy[id], TrimPolicy::Retained)
                && old.is_ok()
                && let Err(error) = new
            {
                return Err(HybridRecoveryError::RetainedTrim {
                    trim: id,
                    error: *error,
                });
            }
        }
        let old = self
            .baseline
            .regions
            .iter()
            .map(|r| (r.region, r))
            .collect::<BTreeMap<_, _>>();
        for region in regions {
            let before = old
                .get(&region.region)
                .ok_or(HybridRecoveryError::InvalidInput)?;
            if before.faces.is_ok()
                && (region.faces.is_err()
                    || (!before.requires_chart_cut && region.requires_chart_cut))
            {
                return Err(HybridRecoveryError::RetainedRegion(region.region));
            }
        }
        Ok(())
    }

    fn constraints(
        &self,
        plan: &HybridPlan,
    ) -> Result<BTreeMap<usize, JunctionConstraints>, HybridRecoveryError> {
        let mut retained = BTreeMap::<usize, BTreeSet<usize>>::new();
        for &id in &plan.touched {
            retained.insert(id, BTreeSet::new());
        }
        // Include all incident retained regions, even point-only
        // contacts. Source topology determines the constraint set.
        for (polygon, &region) in self.partition.region_by_polygon.iter().enumerate() {
            if self.selected.contains(&region) {
                continue;
            }
            for id in self
                .evidence
                .source_polygon_vertices(polygon)
                .ok_or(HybridRecoveryError::InvalidInput)?
            {
                if let Some(regions) = retained.get_mut(id) {
                    regions.insert(region);
                }
            }
        }
        let mut constraints = retained
            .into_iter()
            .map(|(id, regions)| {
                (
                    id,
                    JunctionConstraints {
                        carriers: regions
                            .into_iter()
                            .filter_map(|region| self.partition.carriers[region])
                            .collect(),
                        sections: Vec::new(),
                    },
                )
            })
            .collect::<BTreeMap<_, _>>();
        for (id, trim) in plan.graph.trims.iter().enumerate() {
            if !matches!(plan.policy[id], TrimPolicy::Retained) {
                continue;
            }
            let Ok(intervals) = &plan.intervals[id] else {
                continue;
            };
            let geometry = geometry_constraints(trim);
            let mut seen = BTreeSet::new();
            for vertex in intervals.iter().flat_map(|interval| interval.endpoints) {
                if !seen.insert(vertex) {
                    continue;
                }
                if let Some(constraints) = constraints.get_mut(&vertex) {
                    constraints.carriers.extend(geometry.iter().copied());
                    if let Some(section) = trim.geometry.section() {
                        constraints.sections.push(section);
                    }
                }
            }
        }
        Ok(constraints)
    }

    fn vertices(&self, plan: &HybridPlan) -> Result<Vec<Option<Vertex>>, HybridRecoveryError> {
        let constraints = self.constraints(plan)?;
        let options = NativeTrimOptions {
            max_vertex_displacement: self.patch.budget.vertex,
            intersection_tolerance: self
                .patch
                .intersection_tolerance
                .min(self.strict.intersection_tolerance),
            angular_resolution: self.patch.angular_resolution,
        };
        let mut vertices = self.baseline.vertices.clone();
        for (id, JunctionConstraints { carriers, sections }) in constraints {
            let source = self.evidence.positions[id];
            let branch = sections.as_slice();
            let valid = |p: Point3, tolerance: f64| {
                (p - source).magnitude() <= options.max_vertex_displacement
                    && carriers.iter().all(|carrier| {
                        carrier
                            .signed_constraint(p)
                            .is_some_and(|(_, error)| error.abs() <= tolerance)
                    })
                    && branch.iter().all(|section| {
                        section
                            .projected(p)
                            .is_some_and(|q| (p - q).magnitude() <= tolerance)
                    })
            };
            let fine = NativeTrimOptions {
                intersection_tolerance: junction_roundoff(source, &carriers)
                    .min(options.intersection_tolerance),
                ..options
            };
            let mut point = None;
            for policy in [fine, options] {
                if let Some(old) = &vertices[id]
                    && valid(old.point(), policy.intersection_tolerance)
                {
                    point = Some(old.clone());
                    break;
                }
                let (initial, alternatives) =
                    junction_candidates(source, &carriers, branch, policy);
                if let Some(candidate) = initial
                    .into_iter()
                    .chain(alternatives)
                    .find(|&p| valid(p, policy.intersection_tolerance))
                {
                    point = Some(builder::vertex(candidate));
                    break;
                }
            }
            let point = point.ok_or(NativeTrimError::UnresolvedVertex(id))?;
            vertices[id] = Some(point);
        }
        Ok(vertices)
    }

    fn edges(
        &self,
        plan: &HybridPlan,
        vertices: &[Option<Vertex>],
    ) -> Result<HybridEdges, HybridRecoveryError> {
        let mut trims = Vec::with_capacity(plan.graph.trims.len());
        let mut bounds = BTreeMap::new();
        for (id, trim) in plan.graph.trims.iter().enumerate() {
            match plan.policy[id] {
                TrimPolicy::Retained => {
                    trims.push(build_trim(trim, &plan.intervals[id], vertices, self.strict));
                }
                TrimPolicy::Lifted(carrier) => {
                    let mut edges = Vec::with_capacity(trim.source_edges.len());
                    let mut deviation = 0.0_f64;
                    for edge in &trim.source_edges {
                        let [a, b] = edge.vertices;
                        let front = vertices[a]
                            .as_ref()
                            .ok_or(NativeTrimError::UnresolvedVertex(a))?;
                        let back = vertices[b]
                            .as_ref()
                            .ok_or(NativeTrimError::UnresolvedVertex(b))?;
                        let (native, bound) = lift_native_source_edge(
                            carrier,
                            edge.vertices.map(|id| self.evidence.positions[id]),
                            front,
                            back,
                            RegionPatchOptions {
                                intersection_tolerance: self
                                    .patch
                                    .intersection_tolerance
                                    .min(self.strict.intersection_tolerance),
                                ..self.patch
                            },
                        )
                        .ok_or(HybridRecoveryError::LiftedEdge {
                            trim: id,
                            vertices: edge.vertices,
                        })?;
                        deviation = deviation.max(bound);
                        edges.push(native);
                    }
                    bounds.insert(id, deviation);
                    trims.push(Ok(edges));
                }
            }
        }
        Ok(HybridEdges { trims, bounds })
    }

    fn assemble_patches(
        &self,
        plan: &HybridPlan,
        trims: &[Result<Vec<Edge>, NativeTrimError>],
        regions: &mut Vec<NativeRegionFace>,
    ) -> Result<Vec<HybridPatchEvidence>, HybridRecoveryError> {
        let groups = self
            .partition
            .connected_region_groups(&self.selected.iter().copied().collect::<Vec<_>>())?;
        let by_region = groups
            .iter()
            .enumerate()
            .flat_map(|(group, ids)| ids.iter().map(move |&id| (id, group)))
            .collect::<BTreeMap<_, _>>();
        let mut polygons = vec![Vec::new(); groups.len()];
        for (polygon, region) in self.partition.region_by_polygon.iter().enumerate() {
            if let Some(&group) = by_region.get(region) {
                polygons[group].push(polygon);
            }
        }
        let index = crate::mesh_recovery::RegionBoundaryIndex::new(self.evidence)?;
        let mut supports = vec![Vec::new(); groups.len()];
        for (support, region) in self.partition.region_by_support.iter().enumerate() {
            if let Some(&group) = by_region.get(region) {
                supports[group].push(support);
            }
        }
        let boundaries = supports
            .iter()
            .map(|group| {
                let boundary = index.boundary(group)?;
                if !boundary.pinched_vertices.is_empty() {
                    return Err(HybridRecoveryError::Native(
                        NativeTrimError::InvalidBoundary,
                    ));
                }
                Ok(boundary
                    .loops
                    .iter()
                    .flatten()
                    .flat_map(|edge| edge.vertices)
                    .collect::<BTreeSet<_>>())
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.validate_patch_vertices(&by_region, &boundaries)?;
        let mut prescribed = vec![Vec::new(); groups.len()];
        for (id, trim) in plan.graph.trims.iter().enumerate() {
            if !matches!(plan.policy[id], TrimPolicy::Lifted(_)) {
                continue;
            }
            let group = by_region
                .get(&trim.regions.0)
                .or_else(|| trim.regions.1.and_then(|r| by_region.get(&r)))
                .ok_or(HybridRecoveryError::InvalidInput)?;
            let native = trims[id].as_ref().map_err(|&error| error)?;
            if native.len() != trim.source_edges.len() {
                return Err(HybridRecoveryError::InvalidInput);
            }
            prescribed[*group].extend(trim.source_edges.iter().zip(native).map(
                |(source, edge)| PrescribedSourceEdge {
                    vertices: source.vertices,
                    edge: edge.clone(),
                },
            ));
        }
        groups
            .into_iter()
            .zip(polygons)
            .zip(prescribed)
            .map(|((ids, polygons), edges)| self.assemble_patch(ids, &polygons, &edges, regions))
            .collect()
    }

    // An interior source vertex is allocated within its own patch, not
    // by the shared trim graph. A point-only contact to another component
    // must not silently turn one source identity into two native vertices.
    fn validate_patch_vertices(
        &self,
        by_region: &BTreeMap<usize, usize>,
        boundaries: &[BTreeSet<usize>],
    ) -> Result<(), HybridRecoveryError> {
        let mut interior = BTreeMap::new();
        for (polygon, region) in self.partition.region_by_polygon.iter().enumerate() {
            let Some(&group) = by_region.get(region) else {
                continue;
            };
            for &vertex in self
                .evidence
                .source_polygon_vertices(polygon)
                .ok_or(HybridRecoveryError::InvalidInput)?
            {
                if !boundaries[group].contains(&vertex)
                    && interior
                        .insert(vertex, group)
                        .is_some_and(|old| old != group)
                {
                    return Err(NativeTrimError::InvalidBoundary.into());
                }
            }
        }
        for (polygon, region) in self.partition.region_by_polygon.iter().enumerate() {
            for vertex in self
                .evidence
                .source_polygon_vertices(polygon)
                .ok_or(HybridRecoveryError::InvalidInput)?
            {
                if let Some(owner) = interior.get(vertex)
                    && by_region.get(region) != Some(owner)
                {
                    return Err(NativeTrimError::InvalidBoundary.into());
                }
            }
        }
        Ok(())
    }

    fn assemble_patch(
        &self,
        ids: Vec<usize>,
        polygons: &[usize],
        edges: &[PrescribedSourceEdge],
        regions: &mut Vec<NativeRegionFace>,
    ) -> Result<HybridPatchEvidence, HybridRecoveryError> {
        let source = polygons
            .iter()
            .map(|&id| self.evidence.source_polygon_vertices(id))
            .collect::<Option<Vec<_>>>()
            .ok_or(HybridRecoveryError::InvalidInput)?;
        let owners = polygons
            .iter()
            .map(|&id| self.partition.region_by_polygon[id])
            .collect::<Vec<_>>();
        let approximation = approximate_boundary_patch(
            &self.evidence.positions,
            &source,
            edges,
            BoundaryPatchOptions {
                owner_by_polygon: Some(&owners),
                ..BoundaryPatchOptions::new(self.patch.budget)
            },
        )
        .map_err(|error| HybridRecoveryError::Patch {
            region: ids[0],
            error,
        })?;
        let mut faces = BTreeMap::<usize, Vec<Face>>::new();
        let mut ranges = BTreeMap::<(usize, usize, usize), Range<usize>>::new();
        let mut provenance = Vec::with_capacity(polygons.len());
        for ((&polygon, &region), range) in polygons
            .iter()
            .zip(&owners)
            .zip(&approximation.source_face_ranges)
        {
            let native_range = ranges
                .entry((region, range.start, range.end))
                .or_insert_with(|| {
                    let native = faces.entry(region).or_default();
                    let start = native.len();
                    native.extend(approximation.shell[range.clone()].iter().cloned());
                    start..native.len()
                });
            provenance.push(HybridPolygonFaces {
                polygon,
                region,
                faces: native_range.clone().collect(),
            });
        }
        regions.extend(faces.into_iter().map(|(region, faces)| NativeRegionFace {
            region,
            faces: Ok(faces),
            requires_chart_cut: false,
            chart_seams: 0,
        }));
        Ok(HybridPatchEvidence {
            budget: self.patch.budget,
            regions: ids,
            polygons: provenance,
            coalesced_cells: approximation.coalesced_cells,
            vertex_displacement: approximation.vertex_displacement,
            deviation_upper_bound: approximation.deviation_upper_bound,
            seam_deviation_upper_bound: approximation.seam_deviation_upper_bound,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::patch_lift::tests::spherical_fixture;
    use super::*;
    use crate::mesh_recovery::{RegionBoundaryIndex, SharedTrimOptions, SurfaceRegionBudget};

    #[test]
    fn hybrid_graph_shares_all_native_trims_and_preserves_retained_faces()
    -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, partition, fallback, sphere) = spherical_fixture()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let graph = index.assemble_shared_trims(
            &partition,
            &evidence.positions,
            SharedTrimOptions {
                absolute_tolerance: 1.0e-7,
                relative_tolerance: 1.0e-5,
                max_interval_edges: 100,
            },
        )?;
        let strict = NativeTrimOptions {
            max_vertex_displacement: 1.0e-7,
            intersection_tolerance: 1.0e-10,
            angular_resolution: 1.0e-8,
        };
        let baseline =
            index.construct_native_trims(&partition, &graph, &evidence.positions, strict)?;
        let patch = RegionPatchOptions {
            budget: SurfaceRegionBudget {
                vertex: 1.0e-9,
                chord: 0.1,
            },
            intersection_tolerance: 1.0e-10,
            angular_resolution: 1.0e-8,
        };
        for selected in [&fallback[..1], fallback.as_slice()] {
            let input = HybridInput {
                evidence: &evidence,
                partition: &partition,
                graph: &graph,
                baseline: &baseline,
                positions: &evidence.positions,
                selected: selected.iter().copied().collect(),
                patch,
                strict,
            };
            let result = input.assemble()?;
            let shell = result.complete_shell().ok_or("complete hybrid shell")?;
            assert_eq!(shell.shell_condition(), ShellCondition::Closed);
            let retained = result
                .native
                .regions
                .iter()
                .find(|r| r.region == sphere)
                .ok_or("sphere")?;
            assert!(
                retained
                    .faces
                    .as_ref()
                    .map_err(|&error| error)?
                    .iter()
                    .all(|face| matches!(face.surface(), crate::Surface::Sphere(_)))
            );
            assert!(
                result
                    .patches
                    .iter()
                    .all(|p| p.deviation_upper_bound <= patch.budget.chord)
            );
            let expected = partition
                .region_by_polygon
                .iter()
                .enumerate()
                .filter_map(|(id, r)| selected.contains(r).then_some(id))
                .collect::<BTreeSet<_>>();
            let actual = result
                .patches
                .iter()
                .flat_map(|p| p.polygons.iter().map(|p| p.polygon))
                .collect::<BTreeSet<_>>();
            assert_eq!(actual, expected);
            for &id in result.lifted_trims.keys() {
                let native = result.native.trims[id].as_ref().map_err(|&error| error)?;
                assert_eq!(native.len(), graph.trims[id].source_edges.len());
                for edge in native {
                    assert_eq!(shell.edge_iter().filter(|e| e.id() == edge.id()).count(), 2);
                }
            }
            assert!(
                HybridInput {
                    patch: RegionPatchOptions {
                        budget: SurfaceRegionBudget {
                            chord: 0.02,
                            ..patch.budget
                        },
                        ..patch
                    },
                    ..input
                }
                .assemble()
                .is_err()
            );
        }
        Ok(())
    }
}
