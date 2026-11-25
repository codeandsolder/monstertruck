//! Native shared edges and provisional kernel faces. The source contact
//! graph supplies topological identity; geometry is never independently
//! rebuilt by each incident face. Periodic-domain cuts are separate from
//! the geometric consistency of the native edge and surface carriers.

use super::super::super::sample_fit::normal_basis;
use super::{
    ContactToroidalSectionEvidence, MeshEvidence, RegionBoundaryIndex, SharedAnalyticTrim,
    SharedTrimGeometry, SharedTrimGraph, SurfaceCarrier, SurfaceRegionPartition,
    reconcile::{ReconciliationOptions, reconcile_point},
};
use crate::{
    Curve, Edge, Face, InnerSpace, Invertible, Line, Matrix4, ParametricSurface,
    ParametricSurface3D, Plane, Point3, Processor, RevolutionSurface, Sphere, Surface,
    ToSameGeometry, Torus, TrimmedCurve, UnitCircle, Vector3, Vertex, Wire, builder,
};
use std::{
    collections::BTreeMap,
    f64::consts::{PI, TAU},
};

mod band;
mod hybrid;
mod planar_domain;
pub use hybrid::{
    HybridPatchEvidence, HybridPolygonFaces, HybridRecovery, HybridRecoveryError,
    HybridRefinementOptions, PatchSupportConflict, RefinedHybridRecovery,
};

/// Numerical and displacement policies for native trim vertex reconciliation.
#[derive(Clone, Copy, Debug)]
pub struct NativeTrimOptions {
    /// Maximum additional displacement from the supplied source position table.
    pub max_vertex_displacement: f64,
    /// Residual for all incident analytic constraints at a shared junction.
    pub intersection_tolerance: f64,
    /// Rank resolution of the local normalized constraint system.
    pub angular_resolution: f64,
}

/// Failure of one native construction, not a reason to discard all other
/// independently valid source-backed faces or trim geometries.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum NativeTrimError {
    /// Input membership or source references are invalid.
    #[error("inconsistent native trim source references")]
    InvalidSource,
    /// A trim has no accepted analytic curve.
    #[error("unresolved trim")]
    UnresolvedTrim,
    /// A shared vertex cannot satisfy all incident curve/surface constraints.
    #[error("unresolved shared vertex {0}")]
    UnresolvedVertex(usize),
    /// A curve cannot be represented with consistent native endpoints.
    #[error("native edge geometry inconsistent")]
    EdgeGeometry,
    /// The region has ambiguous or invalid source boundary sectors.
    #[error("invalid region boundary")]
    InvalidBoundary,
    /// A fitted carrier has no regular native chart in the observed domain.
    #[error("unsupported or singular native carrier")]
    UnsupportedCarrier,
    /// Kernel topology construction rejected the oriented shared wires.
    #[error("invalid native face topology")]
    FaceTopology,
    /// A continuous projected native planar trim crosses another trim.
    #[error("self-crossing planar native trim domain")]
    PlanarTrimCrossing,
    /// Shared trims cross on a locally injective revolution chart.
    #[error("crossing shared trim on linear revolution")]
    RevolutionTrimCrossing,
    /// Native curve containment or endpoint consistency failed.
    #[error("native face geometry inconsistent")]
    FaceGeometry,
}

/// One native candidate with original region identity and periodic-domain evidence.
#[derive(Debug)]
pub struct NativeRegionFace {
    /// Original selected region ID.
    pub region: usize,
    /// Kernel-native face using the shared trim edges.
    pub faces: Result<Vec<Face>, NativeTrimError>,
    /// A boundary-free domain, noncontractible loop, pole or ambiguous
    /// outer/inner loop ordering requires an explicit chart decomposition.
    pub requires_chart_cut: bool,
    /// Number of applied internal chart seam sequences.
    pub chart_seams: usize,
}

/// Native graph constructed once from the reconciled global trim graph.
#[derive(Debug)]
pub struct NativeTrimGraph {
    /// Native vertices at original source junction IDs; unused or unresolved
    /// source samples have no native vertex.
    pub vertices: Vec<Option<Vertex>>,
    /// Each shared trim owns one ordered native edge sequence, split only
    /// where topology needs distinct endpoints or an arc exceeds a half turn.
    pub trims: Vec<Result<Vec<Edge>, NativeTrimError>>,
    /// Every original active region, with local failures preserved.
    pub regions: Vec<NativeRegionFace>,
    /// Maximum additional displacement introduced at native junctions.
    pub max_vertex_displacement: f64,
}

#[derive(Clone, Copy)]
struct TrimInterval {
    endpoints: [usize; 2],
    sweep: f64,
}

impl RegionBoundaryIndex<'_> {
    /// Construct native edges and provisional faces from the shared source
    /// graph. This is not a shell-level or periodic chart-domain certificate.
    ///
    /// # Errors
    /// Invalid budgets, position count or partition membership return an error.
    /// Local geometric or topological failures are retained in the result.
    pub fn construct_native_trims(
        &self,
        partition: &SurfaceRegionPartition,
        graph: &SharedTrimGraph,
        positions: &[Point3],
        options: NativeTrimOptions,
    ) -> Result<NativeTrimGraph, NativeTrimError> {
        if positions.len() != self.evidence.positions.len()
            || partition.region_by_support.len() != self.evidence.planar_supports.len()
            || !options.max_vertex_displacement.is_finite()
            || options.max_vertex_displacement <= 0.0
            || !options.intersection_tolerance.is_finite()
            || options.intersection_tolerance <= 0.0
            || !options.angular_resolution.is_finite()
            || options.angular_resolution <= 0.0
        {
            return Err(NativeTrimError::InvalidSource);
        }
        let intervals = graph
            .trims
            .iter()
            .map(|trim| trim_intervals(trim, positions))
            .collect::<Vec<_>>();
        let NativeEdges {
            vertices,
            mut trims,
        } = native_vertices_and_trims(partition, graph, positions, &intervals, options)?;
        let max_vertex_displacement = vertices
            .iter()
            .zip(positions)
            .filter_map(|(vertex, &source)| {
                vertex.as_ref().map(|v| (v.point() - source).magnitude())
            })
            .fold(0.0_f64, f64::max);
        let regions = RegionAssemblyInput {
            evidence: self.evidence,
            partition,
            graph,
            positions,
            options,
            approximated: &std::collections::BTreeSet::new(),
            lifted_deviations: &[],
        }
        .assemble(&mut trims)?;
        Ok(NativeTrimGraph {
            vertices,
            trims,
            regions,
            max_vertex_displacement,
        })
    }
}

struct RegionAssemblyInput<'a> {
    evidence: &'a MeshEvidence,
    partition: &'a SurfaceRegionPartition,
    graph: &'a SharedTrimGraph,
    positions: &'a [Point3],
    options: NativeTrimOptions,
    approximated: &'a std::collections::BTreeSet<usize>,
    lifted_deviations: &'a [Option<f64>],
}

impl RegionAssemblyInput<'_> {
    fn reject_revolution_crossings(&self, trims: &mut [Result<Vec<Edge>, NativeTrimError>]) {
        let mut rejected = std::collections::BTreeSet::new();
        for region in &self.graph.regions {
            if self.approximated.contains(&region.region) {
                continue;
            }
            let Some(SurfaceCarrier::LinearRevolution { origin, axis, .. }) = self
                .partition
                .carriers
                .get(region.region)
                .copied()
                .flatten()
            else {
                continue;
            };
            let Ok(boundary) = &region.boundary else {
                continue;
            };
            let edges = boundary
                .loops
                .iter()
                .flatten()
                .flat_map(|reference| {
                    trims
                        .get(reference.trim)
                        .and_then(|t| t.as_ref().ok())
                        .into_iter()
                        .flatten()
                        .map(move |edge| (reference.trim, edge))
                })
                .collect::<Vec<_>>();
            rejected.extend(planar_domain::revolution_crossings(
                &edges,
                origin,
                axis,
                self.options.intersection_tolerance,
            ));
        }
        for id in rejected {
            trims[id] = Err(NativeTrimError::RevolutionTrimCrossing);
        }
    }

    fn assemble(
        &self,
        trims: &mut [Result<Vec<Edge>, NativeTrimError>],
    ) -> Result<Vec<NativeRegionFace>, NativeTrimError> {
        let mut groups = vec![Vec::new(); self.partition.carriers.len()];
        for (id, &region) in self.partition.region_by_support.iter().enumerate() {
            groups
                .get_mut(region)
                .ok_or(NativeTrimError::InvalidSource)?
                .push(id);
        }
        self.reject_revolution_crossings(trims);
        let plans = self
            .graph
            .regions
            .iter()
            .map(|region| {
                if self.approximated.contains(&region.region) {
                    return None;
                }
                let carrier = self
                    .partition
                    .carriers
                    .get(region.region)
                    .copied()
                    .flatten()?;
                let supports = groups.get(region.region)?;
                let boundary = region.boundary.as_ref().ok()?;
                let carrier = chart_carrier(
                    carrier,
                    boundary,
                    ChartTrimEvidence {
                        graph: self.graph,
                        lifted_deviations: self.lifted_deviations,
                    },
                    self.evidence,
                    supports,
                    self.positions,
                    self.options,
                );
                let (surface, outward) = native_surface(
                    carrier,
                    self.evidence,
                    supports,
                    self.positions,
                    self.options,
                )
                .ok()?;
                if !boundary.pinched_vertices.is_empty() {
                    return None;
                }
                band::BandContext {
                    carrier,
                    surface: &surface,
                    outward,
                    evidence: self.evidence,
                    supports,
                    positions: self.positions,
                    graph: self.graph,
                    options: self.options,
                }
                .plan(boundary, trims)
            })
            .collect::<Vec<_>>();
        let assembly = FaceAssembly {
            evidence: self.evidence,
            partition: self.partition,
            chart: ChartTrimEvidence {
                graph: self.graph,
                lifted_deviations: self.lifted_deviations,
            },
            positions: self.positions,
            options: self.options,
            trims,
            groups: &groups,
        };
        let regions = self
            .graph
            .regions
            .iter()
            .zip(&plans)
            .filter(|(region, _)| !self.approximated.contains(&region.region))
            .map(|(region, plan)| assembly.region(region, plan.as_ref()))
            .collect();
        Ok(regions)
    }
}

#[derive(Clone, Copy)]
struct ChartTrimEvidence<'a> {
    graph: &'a SharedTrimGraph,
    lifted_deviations: &'a [Option<f64>],
}

struct FaceAssembly<'a> {
    evidence: &'a MeshEvidence,
    partition: &'a SurfaceRegionPartition,
    chart: ChartTrimEvidence<'a>,
    positions: &'a [Point3],
    options: NativeTrimOptions,
    trims: &'a [Result<Vec<Edge>, NativeTrimError>],
    groups: &'a [Vec<usize>],
}
struct RegionWires {
    wires: Vec<Wire>,
    areas: Vec<f64>,
    unresolved_chart: bool,
}

impl FaceAssembly<'_> {
    fn region_wires(
        &self,
        boundary: &super::RegionTrimLoops,
        carrier: SurfaceCarrier,
        outward: bool,
    ) -> Result<RegionWires, NativeTrimError> {
        let mut unresolved_chart = false;
        let mut wires = Vec::new();
        let mut areas = Vec::new();
        for cycle in &boundary.loops {
            let mut edges = Vec::new();
            let mut source = Vec::new();
            for reference in cycle {
                let trim = self
                    .chart
                    .graph
                    .trims
                    .get(reference.trim)
                    .ok_or(NativeTrimError::InvalidSource)?;
                let native = self
                    .trims
                    .get(reference.trim)
                    .ok_or(NativeTrimError::InvalidSource)?
                    .as_ref()
                    .map_err(|&error| error)?;
                if reference.reversed {
                    edges.extend(native.iter().rev().map(Edge::inverse));
                    source.extend(trim.source_edges.iter().rev().map(|e| e.vertices[1]));
                } else {
                    edges.extend(native.iter().cloned());
                    source.extend(trim.source_edges.iter().map(|e| e.vertices[0]));
                }
            }
            let area = chart_loop_area(
                carrier,
                self.positions,
                &source,
                self.options.angular_resolution,
            );
            unresolved_chart |= area.is_none();
            areas.push(area.unwrap_or(0.0) * if outward { 1.0 } else { -1.0 });
            wires.push(Wire::from(edges));
        }
        if let SurfaceCarrier::Plane { origin, normal } = carrier
            && planar_domain::has_crossing(
                &wires,
                origin,
                normal,
                self.options.intersection_tolerance,
            )
        {
            return Err(NativeTrimError::PlanarTrimCrossing);
        }
        Ok(RegionWires {
            wires,
            areas,
            unresolved_chart,
        })
    }

    fn region(
        &self,
        region: &super::RegionTrimBoundary,
        plan: Option<&band::BandPlan>,
    ) -> NativeRegionFace {
        let mut requires_chart_cut = false;
        let faces = (|| {
            let carrier = self
                .partition
                .carriers
                .get(region.region)
                .copied()
                .flatten()
                .ok_or(NativeTrimError::InvalidSource)?;
            let boundary = region
                .boundary
                .as_ref()
                .map_err(|_| NativeTrimError::InvalidBoundary)?;
            if !boundary.pinched_vertices.is_empty() {
                return Err(NativeTrimError::InvalidBoundary);
            }
            let supports = self
                .groups
                .get(region.region)
                .ok_or(NativeTrimError::InvalidSource)?;
            let carrier = chart_carrier(
                carrier,
                boundary,
                self.chart,
                self.evidence,
                supports,
                self.positions,
                self.options,
            );
            let (surface, outward) = native_surface(
                carrier,
                self.evidence,
                supports,
                self.positions,
                self.options,
            )?;
            let RegionWires {
                mut wires,
                areas,
                unresolved_chart,
            } = self.region_wires(boundary, carrier, outward)?;
            requires_chart_cut |= unresolved_chart;
            if let Some(plan) = plan {
                let faces = plan.faces(&wires, &surface)?;
                requires_chart_cut = false;
                return Ok(faces);
            }
            let outer = areas
                .iter()
                .enumerate()
                .filter(|(_, area)| **area > 0.0)
                .map(|(id, _)| id)
                .collect::<Vec<_>>();
            requires_chart_cut |= outer.len() != 1;
            if let [id] = outer.as_slice() {
                wires.swap(0, *id);
            }
            let face = Face::new(wires, surface).map_err(|_| NativeTrimError::FaceTopology)?;
            if !face.is_geometric_consistent() {
                return Err(NativeTrimError::FaceGeometry);
            }
            Ok(vec![face])
        })();
        let chart_seams = usize::from(faces.is_ok() && plan.is_some() && !requires_chart_cut) * 2;
        NativeRegionFace {
            region: region.region,
            faces,
            requires_chart_cut,
            chart_seams,
        }
    }
}

fn trim_intervals(
    trim: &SharedAnalyticTrim,
    points: &[Point3],
) -> Result<Vec<TrimInterval>, NativeTrimError> {
    let start = trim
        .source_edges
        .first()
        .ok_or(NativeTrimError::InvalidSource)?
        .vertices[0];
    let end = trim
        .source_edges
        .last()
        .ok_or(NativeTrimError::InvalidSource)?
        .vertices[1];
    match trim.geometry {
        SharedTrimGeometry::Unresolved => Err(NativeTrimError::UnresolvedTrim),
        SharedTrimGeometry::Line(_) => {
            if start == end {
                return Err(NativeTrimError::EdgeGeometry);
            }
            Ok(vec![TrimInterval {
                endpoints: [start, end],
                sweep: 0.0,
            }])
        }
        SharedTrimGeometry::Conic(conic) => {
            angular_intervals(trim, points, start, conic.signed_sweep_radians, |a, b| {
                conic.angular_step(a, b)
            })
        }
        SharedTrimGeometry::ToroidalSection(section) => {
            angular_intervals(trim, points, start, section.signed_sweep_radians, |a, b| {
                section.angular_step(a, b)
            })
        }
        SharedTrimGeometry::Circle(circle) => {
            let mut result = Vec::new();
            let mut first = start;
            let mut previous = start;
            let mut sweep = 0.0;
            for edge in &trim.source_edges {
                let [a, b] = edge.vertices;
                if a != previous {
                    return Err(NativeTrimError::InvalidSource);
                }
                let left = *points.get(a).ok_or(NativeTrimError::InvalidSource)? - circle.center;
                let right = *points.get(b).ok_or(NativeTrimError::InvalidSource)? - circle.center;
                let u = left - circle.normal * left.dot(circle.normal);
                let v = right - circle.normal * right.dot(circle.normal);
                let step = circle.normal.dot(u.cross(v)).atan2(u.dot(v));
                if step * circle.signed_sweep_radians <= 0.0 {
                    return Err(NativeTrimError::EdgeGeometry);
                }
                if (sweep + step).abs() > PI && first != a {
                    result.push(TrimInterval {
                        endpoints: [first, a],
                        sweep,
                    });
                    first = a;
                    sweep = 0.0;
                }
                sweep += step;
                previous = b;
            }
            if first == end {
                return Err(NativeTrimError::EdgeGeometry);
            }
            result.push(TrimInterval {
                endpoints: [first, end],
                sweep,
            });
            Ok(result)
        }
    }
}

fn angular_intervals(
    trim: &SharedAnalyticTrim,
    points: &[Point3],
    start: usize,
    signed_sweep: f64,
    angular_step: impl Fn(Point3, Point3) -> Option<f64>,
) -> Result<Vec<TrimInterval>, NativeTrimError> {
    let mut result = Vec::with_capacity(trim.source_edges.len());
    let mut previous = start;
    for edge in &trim.source_edges {
        let [a, b] = edge.vertices;
        if a != previous {
            return Err(NativeTrimError::InvalidSource);
        }
        let sweep = angular_step(
            *points.get(a).ok_or(NativeTrimError::InvalidSource)?,
            *points.get(b).ok_or(NativeTrimError::InvalidSource)?,
        )
        .ok_or(NativeTrimError::EdgeGeometry)?;
        if sweep * signed_sweep <= 0.0 {
            return Err(NativeTrimError::EdgeGeometry);
        }
        result.push(TrimInterval {
            endpoints: [a, b],
            sweep,
        });
        previous = b;
    }
    Ok(result)
}

fn geometry_constraints(trim: &SharedAnalyticTrim) -> Vec<SurfaceCarrier> {
    match trim.geometry {
        SharedTrimGeometry::Line(line) => <[Vector3; 2]>::from(normal_basis(line.direction))
            .into_iter()
            .map(|normal| SurfaceCarrier::Plane {
                origin: line.endpoints[0],
                normal,
            })
            .collect(),
        SharedTrimGeometry::Circle(circle) => vec![
            SurfaceCarrier::Plane {
                origin: circle.center,
                normal: circle.normal,
            },
            SurfaceCarrier::CircularRevolution {
                axial_center: circle.center,
                axis: circle.normal,
                radial_center: 0.0,
                profile_radius: circle.radius,
            },
        ],
        // A derived bisector plane resolves the branch where the two
        // cylinder gradients become parallel at their intersection.
        SharedTrimGeometry::Conic(conic) => conic.auxiliary_plane().into_iter().collect(),
        SharedTrimGeometry::ToroidalSection(_) | SharedTrimGeometry::Unresolved => Vec::new(),
    }
}

struct NativeEdges {
    vertices: Vec<Option<Vertex>>,
    trims: Vec<Result<Vec<Edge>, NativeTrimError>>,
}

fn native_vertices_and_trims(
    partition: &SurfaceRegionPartition,
    graph: &SharedTrimGraph,
    positions: &[Point3],
    intervals: &[Result<Vec<TrimInterval>, NativeTrimError>],
    options: NativeTrimOptions,
) -> Result<NativeEdges, NativeTrimError> {
    let mut constraints = vec![Vec::new(); positions.len()];
    let mut sections = BTreeMap::<usize, Vec<ContactToroidalSectionEvidence>>::new();
    for (trim, intervals) in graph.trims.iter().zip(intervals) {
        let Ok(intervals) = intervals else {
            continue;
        };
        let mut carriers = geometry_constraints(trim);
        for id in std::iter::once(trim.regions.0).chain(trim.regions.1) {
            carriers.push(
                partition
                    .carriers
                    .get(id)
                    .copied()
                    .flatten()
                    .ok_or(NativeTrimError::InvalidSource)?,
            );
        }
        let ids = intervals
            .iter()
            .flat_map(|i| i.endpoints)
            .collect::<std::collections::BTreeSet<_>>();
        for id in ids {
            constraints
                .get_mut(id)
                .ok_or(NativeTrimError::InvalidSource)?
                .extend(carriers.iter().copied());
            if let SharedTrimGeometry::ToroidalSection(section) = trim.geometry {
                sections.entry(id).or_default().push(section);
            }
        }
    }
    let mut alternatives = BTreeMap::new();
    let mut vertices = positions
        .iter()
        .zip(constraints)
        .enumerate()
        .map(|(id, (&source, carriers))| {
            if carriers.is_empty() {
                return None;
            }
            let (initial, candidates) = junction_candidates(
                source,
                &carriers,
                sections.get(&id).map_or(&[], Vec::as_slice),
                options,
            );
            if !candidates.is_empty() {
                alternatives.insert(id, candidates);
            }
            initial.map(builder::vertex)
        })
        .collect::<Vec<_>>();
    let mut trims = graph
        .trims
        .iter()
        .zip(intervals)
        .map(|(trim, intervals)| build_trim(trim, intervals, &vertices, options))
        .collect::<Vec<_>>();
    refine_native_junctions(
        graph,
        intervals,
        &mut vertices,
        &mut trims,
        alternatives,
        options,
    );
    Ok(NativeEdges { vertices, trims })
}

fn build_trim(
    trim: &SharedAnalyticTrim,
    intervals: &Result<Vec<TrimInterval>, NativeTrimError>,
    vertices: &[Option<Vertex>],
    options: NativeTrimOptions,
) -> Result<Vec<Edge>, NativeTrimError> {
    intervals
        .as_ref()
        .map_err(|&error| error)?
        .iter()
        .map(|&interval| native_edge(trim, interval, vertices, options))
        .collect()
}

// Test candidate vertex changes transactionally against every incident
// trim. An improvement to one intersection must not invalidate a
// previously constructed line, circular arc, conic or section.
fn refine_native_junctions(
    graph: &SharedTrimGraph,
    intervals: &[Result<Vec<TrimInterval>, NativeTrimError>],
    vertices: &mut [Option<Vertex>],
    trims: &mut [Result<Vec<Edge>, NativeTrimError>],
    alternatives: BTreeMap<usize, Vec<Point3>>,
    options: NativeTrimOptions,
) {
    if alternatives.is_empty() {
        return;
    }
    let mut incidence = alternatives
        .keys()
        .map(|&id| (id, Vec::new()))
        .collect::<BTreeMap<_, _>>();
    for (trim, intervals) in intervals.iter().enumerate() {
        if let Ok(intervals) = intervals {
            for id in intervals.iter().flat_map(|interval| interval.endpoints) {
                if let Some(incident) = incidence.get_mut(&id)
                    && incident.last() != Some(&trim)
                {
                    incident.push(trim);
                }
            }
        }
    }
    for (id, candidates) in alternatives {
        let original = vertices[id].take();
        let incident = &incidence[&id];
        for point in candidates {
            vertices[id] = Some(builder::vertex(point));
            let changed = incident
                .iter()
                .map(|&trim| {
                    (
                        trim,
                        build_trim(&graph.trims[trim], &intervals[trim], vertices, options),
                    )
                })
                .collect::<Vec<_>>();
            if changed
                .iter()
                .all(|(id, result)| trims[*id].is_err() || result.is_ok())
            {
                for (id, result) in changed {
                    trims[id] = result;
                }
                break;
            }
            vertices[id] = None;
        }
        if vertices[id].is_none() {
            vertices[id] = original;
        }
    }
}

// A carrier residual by itself does not imply proximity to an
// intersection branch near parallel normals. Use the analytic sections
// as additional point candidates, not as permission to move one side
// independently. Keep the existing solution if no globally compatible
// refinement improves it.
fn junction_candidates(
    source: Point3,
    carriers: &[SurfaceCarrier],
    sections: &[ContactToroidalSectionEvidence],
    options: NativeTrimOptions,
) -> (Option<Point3>, Vec<Point3>) {
    let reconciliation = ReconciliationOptions {
        max_displacement: options.max_vertex_displacement,
        parallel_angle_tolerance: options.angular_resolution,
        intersection_tolerance: options.intersection_tolerance,
        max_normal_angle: options.angular_resolution * 2.0,
    };
    let initial = reconcile_point(source, carriers, reconciliation);
    if sections.is_empty() {
        return (initial, Vec::new());
    }
    let branch_error = |point: Point3| {
        sections
            .iter()
            .map(|section| {
                section
                    .projected(point)
                    .map_or(f64::INFINITY, |p| (p - point).magnitude())
            })
            .fold(0.0_f64, f64::max)
    };
    let seed = initial.unwrap_or(source);
    let candidates = sections.iter().filter_map(|s| s.projected(seed));
    // Do not perturb a junction whose branch constraints already agree.
    if initial.is_some_and(|p| branch_error(p) <= options.intersection_tolerance) {
        return (initial, Vec::new());
    }
    let mut candidates = candidates
        .filter(|&point| {
            (point - source).magnitude() <= options.max_vertex_displacement
                && carriers.iter().all(|carrier| {
                    carrier
                        .signed_constraint(point)
                        .is_some_and(|(_, error)| error.abs() <= options.intersection_tolerance)
                })
                && initial.is_none_or(|original| branch_error(point) < branch_error(original))
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|&a, &b| {
        branch_error(a).total_cmp(&branch_error(b)).then(
            (a - source)
                .magnitude2()
                .total_cmp(&(b - source).magnitude2()),
        )
    });
    (initial, candidates)
}

fn native_edge(
    trim: &SharedAnalyticTrim,
    interval: TrimInterval,
    vertices: &[Option<Vertex>],
    options: NativeTrimOptions,
) -> Result<Edge, NativeTrimError> {
    let [a, b] = interval.endpoints;
    let front = vertices
        .get(a)
        .and_then(Option::as_ref)
        .ok_or(NativeTrimError::UnresolvedVertex(a))?;
    let back = vertices
        .get(b)
        .and_then(Option::as_ref)
        .ok_or(NativeTrimError::UnresolvedVertex(b))?;
    let edge = match trim.geometry {
        SharedTrimGeometry::Line(_) => builder::line(front, back),
        SharedTrimGeometry::Circle(circle) => native_circle_edge(
            front,
            back,
            circle.center,
            circle.normal,
            circle.radius,
            interval.sweep,
            options.intersection_tolerance,
        )?,
        SharedTrimGeometry::Conic(conic) => {
            let curve = conic
                .native_curve(
                    front.point(),
                    back.point(),
                    interval.sweep,
                    options.intersection_tolerance,
                )
                .ok_or(NativeTrimError::EdgeGeometry)?;
            Edge::new(front, back, curve).map_err(|_| NativeTrimError::EdgeGeometry)?
        }
        SharedTrimGeometry::ToroidalSection(section) => {
            let curve = section
                .native_curve(
                    front.point(),
                    back.point(),
                    interval.sweep,
                    options.intersection_tolerance,
                )
                .ok_or(NativeTrimError::EdgeGeometry)?;
            Edge::new(front, back, curve).map_err(|_| NativeTrimError::EdgeGeometry)?
        }
        SharedTrimGeometry::Unresolved => return Err(NativeTrimError::UnresolvedTrim),
    };
    if !edge.is_geometric_consistent() {
        return Err(NativeTrimError::EdgeGeometry);
    }
    Ok(edge)
}

fn native_circle_edge(
    front: &Vertex,
    back: &Vertex,
    center: Point3,
    normal: Vector3,
    radius: f64,
    sweep: f64,
    tolerance: f64,
) -> Result<Edge, NativeTrimError> {
    let delta = front.point() - center;
    let radial = delta - normal * delta.dot(normal);
    if radial.magnitude() <= 0.0 {
        return Err(NativeTrimError::EdgeGeometry);
    }
    let x = radial.normalize() * radius;
    let axis = normal * sweep.signum();
    let y = axis.cross(x);
    let delta = back.point() - center;
    let angle = delta.dot(y).atan2(delta.dot(x)).rem_euclid(TAU);
    if (angle - sweep.abs()).abs() * radius > tolerance * 4.0 {
        return Err(NativeTrimError::EdgeGeometry);
    }
    let transform = Matrix4::from_cols(
        x.extend(0.0),
        y.extend(0.0),
        axis.extend(0.0),
        center.to_homogeneous(),
    );
    let curve = Processor::with_transform(
        TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, angle)),
        transform,
    )
    .to_same_geometry();
    Edge::new(front, back, curve).map_err(|_| NativeTrimError::EdgeGeometry)
}

// A sphere has no intrinsic axis. Align with a complete circular rim
// only for the cap planner; a general circular arc does not define a
// suitable chart axis for the rest of a patch.
fn chart_carrier(
    carrier: SurfaceCarrier,
    boundary: &super::RegionTrimLoops,
    chart: ChartTrimEvidence<'_>,
    evidence: &MeshEvidence,
    supports: &[usize],
    points: &[Point3],
    options: NativeTrimOptions,
) -> SurfaceCarrier {
    let SurfaceCarrier::CircularRevolution {
        axial_center,
        radial_center,
        profile_radius,
        ..
    } = carrier
    else {
        return carrier;
    };
    if radial_center.abs() > crate::KERNEL_EPSILON {
        return carrier;
    }
    let Some(axis) =
        circular_rim_axis(axial_center, boundary, chart.graph, options).or_else(|| {
            compact_sphere_axis(
                axial_center,
                boundary,
                chart,
                evidence,
                supports,
                points,
                options,
            )
        })
    else {
        return carrier;
    };
    SurfaceCarrier::CircularRevolution {
        axial_center,
        axis,
        radial_center: 0.0,
        profile_radius,
    }
}

fn circular_rim_axis(
    center: Point3,
    boundary: &super::RegionTrimLoops,
    graph: &SharedTrimGraph,
    options: NativeTrimOptions,
) -> Option<Vector3> {
    let [cycle] = boundary.loops.as_slice() else {
        return None;
    };
    let SharedTrimGeometry::Circle(first) = graph.trims.get(cycle.first()?.trim)?.geometry else {
        return None;
    };
    let tolerance = options.intersection_tolerance;
    if (first.center - center).cross(first.normal).magnitude() > tolerance {
        return None;
    }
    let mut sweep = 0.0;
    for reference in cycle {
        let SharedTrimGeometry::Circle(circle) = graph.trims.get(reference.trim)?.geometry else {
            return None;
        };
        if !super::trim_intersections::same_circle(&first, &circle, tolerance) {
            return None;
        }
        let orientation = if reference.reversed { -1.0 } else { 1.0 };
        sweep = (circle.signed_sweep_radians * orientation)
            .mul_add(first.normal.dot(circle.normal).signum(), sweep);
    }
    ((sweep.abs() - TAU).abs() <= options.angular_resolution).then_some(first.normal)
}

// A source patch contained in an open hemisphere admits a regular
// longitude/latitude chart with both poles in the separating great circle.
// Use the same halfspace for source facets AND the continuous trim arcs;
// a collection of vertices alone is not sufficient to certify the chart.
fn compact_sphere_axis(
    center: Point3,
    boundary: &super::RegionTrimLoops,
    chart: ChartTrimEvidence<'_>,
    evidence: &MeshEvidence,
    supports: &[usize],
    points: &[Point3],
    options: NativeTrimOptions,
) -> Option<Vector3> {
    let mut vertices = std::collections::BTreeSet::new();
    for &support in supports {
        for &polygon in &evidence.planar_supports.get(support)?.source_face_indices {
            vertices.extend(evidence.source_polygon_vertices(polygon)?.iter().copied());
        }
    }
    let mut direction = Vector3::new(0.0, 0.0, 0.0);
    for &id in &vertices {
        direction += *points.get(id)? - center;
    }
    if !direction.magnitude2().is_finite() || direction.magnitude2() == 0.0 {
        return None;
    }
    direction = direction.normalize();
    let margin = options.max_vertex_displacement + options.intersection_tolerance;
    if vertices
        .iter()
        .any(|&id| (points[id] - center).dot(direction) <= margin)
    {
        return None;
    }
    // Every source facet lies in this convex halfspace. For a minor circle
    // arc the distance from its projected chord is at most its sagitta.
    // Add the original endpoint projection errors and junction movement.
    for reference in boundary.loops.iter().flatten() {
        let trim = chart.graph.trims.get(reference.trim)?;
        let deviation = match trim.geometry {
            SharedTrimGeometry::Line(line) => line.max_vertex_offset,
            SharedTrimGeometry::Circle(circle) => {
                let sagitta = circle.radius * (1.0 - (circle.max_step_radians * 0.5).cos());
                sagitta
                    + circle
                        .max_plane_residual
                        .hypot(circle.max_vertex_radial_residual)
            }
            SharedTrimGeometry::Conic(_) | SharedTrimGeometry::ToroidalSection(_) => return None,
            SharedTrimGeometry::Unresolved => {
                *chart.lifted_deviations.get(reference.trim)?.as_ref()?
            }
        };
        if trim
            .source_edges
            .iter()
            .flat_map(|edge| edge.vertices)
            .any(|id| {
                points
                    .get(id)
                    .is_none_or(|&point| (point - center).dot(direction) <= deviation + margin)
            })
        {
            return None;
        }
    }
    Some(normal_basis(direction).0)
}

fn native_surface(
    carrier: SurfaceCarrier,
    evidence: &MeshEvidence,
    supports: &[usize],
    points: &[Point3],
    options: NativeTrimOptions,
) -> Result<(Surface, bool), NativeTrimError> {
    let outward = region_outward(carrier, evidence, supports, points)?;
    let (mut surface, probe) = match carrier {
        SurfaceCarrier::Plane { origin, normal } => {
            let (u, v) = normal_basis(normal);
            (
                Surface::Plane(Plane::new(origin, origin + u, origin + v)),
                [0.0, 0.0],
            )
        }
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } => {
            let (lo, hi) = axial_span(evidence, supports, points, origin, axis)?;
            let (u, _) = normal_basis(axis);
            if lo >= hi {
                return Err(NativeTrimError::UnsupportedCarrier);
            }
            let front = linear_profile_point(
                origin,
                axis,
                u,
                radius_at_origin,
                slope,
                lo,
                options.intersection_tolerance,
            )
            .ok_or(NativeTrimError::UnsupportedCarrier)?;
            let back = linear_profile_point(
                origin,
                axis,
                u,
                radius_at_origin,
                slope,
                hi,
                options.intersection_tolerance,
            )
            .ok_or(NativeTrimError::UnsupportedCarrier)?;
            let surface = Surface::Revolution(Processor::new(RevolutionSurface::by_revolution(
                Curve::Line(Line(front, back)),
                origin,
                axis,
            )));
            (surface, [0.5, 0.0])
        }
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } => {
            let (u, v) = normal_basis(axis);
            let transform = Matrix4::from_cols(
                u.extend(0.0),
                v.extend(0.0),
                axis.extend(0.0),
                axial_center.to_homogeneous(),
            );
            if radial_center.abs() <= crate::KERNEL_EPSILON {
                (
                    Surface::Sphere(Processor::with_transform(
                        Sphere::new(Point3::new(0.0, 0.0, 0.0), profile_radius),
                        transform,
                    )),
                    [PI * 0.5, 0.0],
                )
            } else {
                let torus = Torus::new(Point3::new(0.0, 0.0, 0.0), radial_center, profile_radius)
                    .map_err(|_| NativeTrimError::UnsupportedCarrier)?;
                (
                    Surface::Torus(Processor::with_transform(torus, transform)),
                    [0.0, 0.0],
                )
            }
        }
    };
    let point = surface.evaluate(probe[0], probe[1]);
    let target = carrier
        .signed_constraint(point)
        .ok_or(NativeTrimError::UnsupportedCarrier)?
        .0;
    if (surface.normal(probe[0], probe[1]).dot(target) > 0.0) != outward {
        surface.invert();
    }
    Ok((surface, outward))
}

fn region_outward(
    carrier: SurfaceCarrier,
    evidence: &MeshEvidence,
    supports: &[usize],
    points: &[Point3],
) -> Result<bool, NativeTrimError> {
    let support = supports
        .iter()
        .map(|&id| &evidence.planar_supports[id])
        .max_by(|a, b| a.area.total_cmp(&b.area))
        .ok_or(NativeTrimError::InvalidSource)?;
    let ids = evidence
        .source_polygon_vertices(
            *support
                .source_face_indices
                .first()
                .ok_or(NativeTrimError::InvalidSource)?,
        )
        .ok_or(NativeTrimError::InvalidSource)?;
    let origin = points[ids[0]];
    let centroid = origin
        + ids.iter().fold(Vector3::new(0.0, 0.0, 0.0), |sum, &id| {
            sum + (points[id] - origin)
        }) / crate::usize_to_f64(ids.len());
    let normal = carrier
        .signed_constraint(centroid)
        .ok_or(NativeTrimError::UnsupportedCarrier)?
        .0;
    Ok(normal.dot(support.normal) > 0.0)
}

// Near a cone tip, roundoff can make the fitted radius slightly
// negative. Snap to the exact apex only within the intersection budget.
fn linear_profile_point(
    origin: Point3,
    axis: Vector3,
    u: Vector3,
    radius_at_origin: f64,
    slope: f64,
    z: f64,
    tolerance: f64,
) -> Option<Point3> {
    let radius = slope.mul_add(z, radius_at_origin);
    if slope != 0.0 {
        let apex_z = -radius_at_origin / slope;
        if (z - apex_z).hypot(radius) <= tolerance {
            return Some(origin + axis * apex_z);
        }
    }
    (radius >= 0.0).then_some(origin + axis * z + u * radius)
}

fn axial_span(
    evidence: &MeshEvidence,
    supports: &[usize],
    points: &[Point3],
    origin: Point3,
    axis: Vector3,
) -> Result<(f64, f64), NativeTrimError> {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &id in supports {
        for &polygon in &evidence.planar_supports[id].source_face_indices {
            for &v in evidence
                .source_polygon_vertices(polygon)
                .ok_or(NativeTrimError::InvalidSource)?
            {
                let z = (points[v] - origin).dot(axis);
                min = min.min(z);
                max = max.max(z);
            }
        }
    }
    Ok((min, max))
}

// Map original directed source loops to a locally unwrapped analytic chart.
// A nonzero winding in either periodic coordinate is not a contractible
// polygonal loop. It needs a topological seam, not a fake closing segment.
fn chart_loop_area(
    carrier: SurfaceCarrier,
    points: &[Point3],
    ids: &[usize],
    resolution: f64,
) -> Option<f64> {
    let first = *points.get(*ids.first()?)?;
    let mut previous = chart_point(carrier, first, resolution)?;
    let origin = previous;
    let mut area = 0.0;
    for &id in ids.iter().skip(1).chain(ids.first()) {
        let mut point = chart_point(carrier, *points.get(id)?, resolution)?;
        for axis in 0..2 {
            if periodic_coordinate(carrier, axis) {
                point[axis] =
                    previous[axis] + (point[axis] - previous[axis] + PI).rem_euclid(TAU) - PI;
            }
        }
        let [a, b] = [previous[0] - origin[0], previous[1] - origin[1]];
        let [c, d] = [point[0] - origin[0], point[1] - origin[1]];
        area += a.mul_add(d, -b * c);
        previous = point;
    }
    ((previous[0] - origin[0]).abs() <= resolution && (previous[1] - origin[1]).abs() <= resolution)
        .then_some(area * 0.5)
}

const fn periodic_coordinate(carrier: SurfaceCarrier, axis: usize) -> bool {
    match carrier {
        SurfaceCarrier::Plane { .. } => false,
        SurfaceCarrier::LinearRevolution { .. } => axis == 0,
        SurfaceCarrier::CircularRevolution { radial_center, .. } => {
            axis == 0 || radial_center > crate::KERNEL_EPSILON
        }
    }
}

fn chart_point(carrier: SurfaceCarrier, point: Point3, resolution: f64) -> Option<[f64; 2]> {
    match carrier {
        SurfaceCarrier::Plane { origin, normal } => {
            let (u, v) = normal_basis(normal);
            let delta = point - origin;
            Some([delta.dot(u), delta.dot(v)])
        }
        SurfaceCarrier::LinearRevolution { origin, axis, .. } => {
            let (u, v) = normal_basis(axis);
            let delta = point - origin;
            let (x, y) = (delta.dot(u), delta.dot(v));
            (x.hypot(y) > resolution).then(|| [y.atan2(x), delta.dot(axis)])
        }
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            ..
        } => {
            let (u, v) = normal_basis(axis);
            let delta = point - axial_center;
            let (x, y) = (delta.dot(u), delta.dot(v));
            let rho = x.hypot(y);
            (rho > resolution).then(|| [y.atan2(x), delta.dot(axis).atan2(rho - radial_center)])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::{SharedTrimOptions, SurfaceRegionBudget};
    use super::*;
    use crate::mesh_recovery::{EvidenceOptions, analyze_mesh_evidence};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    fn trim_options() -> SharedTrimOptions {
        SharedTrimOptions {
            absolute_tolerance: 1.0e-8,
            relative_tolerance: 1.0e-4,
            max_interval_edges: 128,
        }
    }
    fn native_options() -> NativeTrimOptions {
        NativeTrimOptions {
            max_vertex_displacement: 1.0e-8,
            intersection_tolerance: 1.0e-10,
            angular_resolution: 1.0e-8,
        }
    }
    #[test]
    fn native_cube_owns_each_shared_edge_once() -> Result<(), Box<dyn std::error::Error>> {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ];
        let faces = Faces::from_iter([
            [3, 2, 1, 0],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ]);
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            faces,
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let partition = index.select_surface_partition(&[], &evidence, 100)?;
        let trims = index.assemble_shared_trims(&partition, &evidence.positions, trim_options())?;
        let native = index.construct_native_trims(
            &partition,
            &trims,
            &evidence.positions,
            native_options(),
        )?;
        assert_eq!(native.trims.len(), 12);
        assert!(
            native
                .regions
                .iter()
                .all(|r| r.faces.is_ok() && !r.requires_chart_cut)
        );
        let shell = native
            .regions
            .into_iter()
            .flat_map(|r| match r.faces {
                Ok(faces) => faces.into_iter().map(Ok).collect::<Vec<_>>(),
                Err(error) => vec![Err(error)],
            })
            .collect::<Result<crate::Shell, _>>()?;
        assert_eq!(shell.len(), 6);
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Closed
        );
        Ok(())
    }

    #[test]
    fn cylindrical_band_has_shared_native_ring_edges_and_two_chart_patches()
    -> Result<(), Box<dyn std::error::Error>> {
        let ring = 8_usize;
        let mut positions = Vec::new();
        for z in [0.0, 2.0] {
            for i in 0..ring {
                let (sin, cos) =
                    (TAU * crate::usize_to_f64(i) / crate::usize_to_f64(ring)).sin_cos();
                positions.push(Point3::new(cos, sin, z));
            }
        }
        let mut faces = Vec::new();
        for i in 0..ring {
            let next = (i + 1) % ring;
            faces.push(vec![i, next, next + ring, i + ring]);
        }
        faces.push((0..ring).rev().collect());
        faces.push((ring..2 * ring).collect());
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let mut supports = (0..ring)
            .map(|id| evidence.source_planar_support[id])
            .collect::<Vec<_>>();
        supports.sort_unstable();
        supports.dedup();
        let carrier = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let proposal = index
            .measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: 1.0e-8,
                    chord: 0.25,
                },
            )
            .ok_or("cylinder")?;
        let partition = index.select_surface_partition(&[proposal], &evidence, 100)?;
        let trims = index.assemble_shared_trims(&partition, &evidence.positions, trim_options())?;
        let native = index.construct_native_trims(
            &partition,
            &trims,
            &evidence.positions,
            native_options(),
        )?;
        assert_eq!(native.regions.len(), 3);
        assert_eq!(native.trims.len(), 2);
        assert!(
            native
                .trims
                .iter()
                .all(|trim| trim.as_ref().is_ok_and(|edges| edges.len() >= 2))
        );
        assert_eq!(
            native
                .regions
                .iter()
                .filter(|r| r.requires_chart_cut)
                .count(),
            0
        );
        assert_eq!(
            native.regions.iter().map(|r| r.chart_seams).sum::<usize>(),
            2
        );
        let shell = native
            .regions
            .into_iter()
            .flat_map(|r| match r.faces {
                Ok(faces) => faces.into_iter().map(Ok).collect::<Vec<_>>(),
                Err(error) => vec![Err(error)],
            })
            .collect::<Result<crate::Shell, _>>()?;
        assert_eq!(shell.len(), 4);
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Closed
        );
        Ok(())
    }
    #[test]
    fn toroidal_band_splits_the_periodic_meridian_into_shared_chart_seams()
    -> Result<(), Box<dyn std::error::Error>> {
        let ring = 8_usize;
        let stations = 5_usize;
        let mut positions = Vec::new();
        for i in 0..stations {
            let theta = PI * 0.5 * crate::usize_to_f64(i) / crate::usize_to_f64(stations - 1);
            let (sin, cos) = theta.sin_cos();
            for j in 0..ring {
                let (z, rho) = (TAU * crate::usize_to_f64(j) / crate::usize_to_f64(ring)).sin_cos();
                positions.push(Point3::new((3.0 + rho) * cos, (3.0 + rho) * sin, z));
            }
        }
        let mut faces = Vec::new();
        for i in 0..stations - 1 {
            for j in 0..ring {
                let next = (j + 1) % ring;
                faces.push(vec![
                    i * ring + j,
                    (i + 1) * ring + j,
                    (i + 1) * ring + next,
                    i * ring + next,
                ]);
            }
        }
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let supports = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: 3.0,
            profile_radius: 1.0,
        };
        let proposal = index
            .measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: 1.0e-8,
                    chord: 0.5,
                },
            )
            .ok_or("torus")?;
        let partition = index.select_surface_partition(&[proposal], &evidence, 100)?;
        let trims = index.assemble_shared_trims(&partition, &evidence.positions, trim_options())?;
        let native = index.construct_native_trims(
            &partition,
            &trims,
            &evidence.positions,
            native_options(),
        )?;
        assert_eq!(native.regions.len(), 1);
        let region = &native.regions[0];
        assert!(!region.requires_chart_cut);
        assert_eq!(region.chart_seams, 2);
        let faces = region.faces.as_ref().map_err(|err| *err)?;
        assert_eq!(faces.len(), 2);
        let shell = faces.iter().cloned().collect::<crate::Shell>();
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Oriented
        );
        Ok(())
    }
    #[test]
    fn junction_refinement_preserves_existing_native_trim() -> Result<(), Box<dyn std::error::Error>>
    {
        let points = [Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)];
        let circle = super::super::curve_evidence::measure_placed_circle(
            &points,
            &[0, 1],
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            1.0,
            1.0,
        )
        .ok_or("circle")?;
        let graph = SharedTrimGraph {
            trims: vec![SharedAnalyticTrim {
                contact_path: 0,
                regions: (0, Some(1)),
                source_edges: Vec::new(),
                geometry: SharedTrimGeometry::Circle(circle),
            }],
            regions: Vec::new(),
            ambiguous_edges: Vec::new(),
            limited_paths: 0,
        };
        let intervals = [Ok(vec![TrimInterval {
            endpoints: [0, 1],
            sweep: PI * 0.5,
        }])];
        let mut vertices = points.map(|p| Some(builder::vertex(p)));
        let mut trims = [build_trim(
            &graph.trims[0],
            &intervals[0],
            &vertices,
            native_options(),
        )];
        let original = vertices[0].clone().ok_or("vertex")?;
        let (sin, cos) = 1.0e-9_f64.sin_cos();
        let candidate = Point3::new(cos, sin, 0.0);
        assert!(
            (candidate - points[0]).magnitude() > native_options().intersection_tolerance * 4.0
        );
        refine_native_junctions(
            &graph,
            &intervals,
            &mut vertices,
            &mut trims,
            BTreeMap::from([(0, vec![candidate])]),
            native_options(),
        );
        assert_eq!(vertices[0].as_ref(), Some(&original));
        assert!(trims[0].is_ok());
        assert_eq!(
            trims[0].as_ref().map_err(|error| *error)?[0].front(),
            &original
        );
        Ok(())
    }

    #[test]
    fn noncircular_torus_plane_section_is_shared_by_native_faces()
    -> Result<(), Box<dyn std::error::Error>> {
        // A contractible toroidal patch, bounded by two parallel circles,
        // a meridian circle and a noncircular plane intersection.
        let steps = 8;
        let width = steps + 1;
        let mut positions = Vec::new();
        for i in 0..=steps {
            let beta = 0.9 * crate::usize_to_f64(i) / crate::usize_to_f64(steps) - 0.3;
            let (sin, cos) = beta.sin_cos();
            let radius = 0.5_f64.mul_add(cos, 2.0);
            let end = (0.7 / radius).acos();
            for j in 0..=steps {
                let theta =
                    (end - 0.5).mul_add(crate::usize_to_f64(j) / crate::usize_to_f64(steps), 0.5);
                let (sine, cosine) = theta.sin_cos();
                positions.push(Point3::new(radius * cosine, radius * sine, 0.5 * sin));
            }
        }
        let mut faces = Vec::new();
        for i in 0..steps {
            for j in 0..steps {
                let a = i * width + j;
                faces.push(vec![a, a + 1, a + width + 1]);
                faces.push(vec![a, a + width + 1, a + width]);
            }
        }
        let patch_facets = faces.len();
        // Reverse the boundary traversal along the plane side.
        faces.push((0..=steps).rev().map(|i| i * width + steps).collect());
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let supports = evidence
            .planar_supports
            .iter()
            .enumerate()
            .filter(|(_, support)| {
                support
                    .source_face_indices
                    .iter()
                    .all(|&i| i < patch_facets)
            })
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: 2.0,
            profile_radius: 0.5,
        };
        let proposal = index
            .measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: 1.0e-8,
                    chord: 0.1,
                },
            )
            .ok_or("torus patch")?;
        let partition = index.select_surface_partition(&[proposal], &evidence, 100)?;
        let trims = index.assemble_shared_trims(&partition, &evidence.positions, trim_options())?;
        let section = trims
            .trims
            .iter()
            .position(|trim| matches!(trim.geometry, SharedTrimGeometry::ToroidalSection(_)))
            .ok_or("shared torus-plane trim")?;
        assert!(trims.trims[section].regions.1.is_some());
        let native = index.construct_native_trims(
            &partition,
            &trims,
            &evidence.positions,
            native_options(),
        )?;
        let edges = native.trims[section].as_ref().map_err(|error| *error)?;
        assert_eq!(edges.len(), steps);
        assert_eq!(native.regions.len(), 2);
        let shell = native
            .regions
            .into_iter()
            .flat_map(|region| match region.faces {
                Ok(faces) => faces.into_iter().map(Ok).collect::<Vec<_>>(),
                Err(error) => vec![Err(error)],
            })
            .collect::<Result<crate::Shell, _>>()?;
        assert_eq!(shell.len(), 2);
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Oriented
        );
        Ok(())
    }

    #[test]
    fn spherical_patch_with_multiple_distinct_rims_avoids_chart_poles()
    -> Result<(), Box<dyn std::error::Error>> {
        // Three circular boundaries meeting at the original chart's poles:
        // the source patch is an octant, not an annulus or a complete cap.
        let steps = 6_usize;
        let width = steps + 1;
        let mut positions = Vec::new();
        for i in 0..steps {
            let elevation = PI * 0.5 * crate::usize_to_f64(i) / crate::usize_to_f64(steps);
            let (z, rho) = elevation.sin_cos();
            for j in 0..=steps {
                let angle = PI * 0.5 * crate::usize_to_f64(j) / crate::usize_to_f64(steps);
                let (sin, cos) = angle.sin_cos();
                positions.push(Point3::new(rho * cos, rho * sin, z));
            }
        }
        let pole = positions.len();
        positions.push(Point3::new(0.0, 0.0, 1.0));
        let center = positions.len();
        positions.push(Point3::new(0.0, 0.0, 0.0));
        let mut faces = Vec::new();
        for i in 0..steps {
            for j in 0..steps {
                let a = i * width + j;
                if i + 1 == steps {
                    faces.push(vec![a, a + 1, pole]);
                } else {
                    faces.push(vec![a, a + 1, a + 1 + width, a + width]);
                }
            }
        }
        let sphere_polygons = faces.len();
        faces.push(std::iter::once(center).chain((0..=steps).rev()).collect());
        faces.push(
            std::iter::once(center)
                .chain((0..steps).map(|i| i * width))
                .chain(std::iter::once(pole))
                .collect(),
        );
        faces.push(
            [center, pole]
                .into_iter()
                .chain((0..steps).rev().map(|i| i * width + steps))
                .collect(),
        );
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let mut supports = evidence.source_planar_support[..sphere_polygons].to_vec();
        supports.sort_unstable();
        supports.dedup();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: 0.0,
            profile_radius: 1.0,
        };
        let proposal = index
            .measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: 1.0e-8,
                    chord: 0.1,
                },
            )
            .ok_or("sphere patch")?;
        let partition = index.select_surface_partition(&[proposal], &evidence, 100)?;
        let trims = index.assemble_shared_trims(&partition, &evidence.positions, trim_options())?;
        let native = index.construct_native_trims(
            &partition,
            &trims,
            &evidence.positions,
            native_options(),
        )?;
        assert_eq!(native.regions.len(), 4);
        assert!(
            native
                .regions
                .iter()
                .all(|r| !r.requires_chart_cut && r.chart_seams == 0)
        );
        let mut shell = crate::Shell::new();
        for region in native.regions {
            shell.extend(region.faces?);
        }
        assert_eq!(shell.len(), 4);
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Closed
        );
        assert!(shell.is_geometric_consistent());
        Ok(())
    }

    #[test]
    fn spherical_cap_reorients_the_arbitrary_axis_and_shares_its_pole()
    -> Result<(), Box<dyn std::error::Error>> {
        let ring = 8_usize;
        let mut positions = Vec::new();
        for elevation in [0.0, PI * 0.25] {
            let (z, rho) = elevation.sin_cos();
            for j in 0..ring {
                let (sin, cos) =
                    (TAU * crate::usize_to_f64(j) / crate::usize_to_f64(ring)).sin_cos();
                positions.push(Point3::new(rho * cos, rho * sin, z));
            }
        }
        let pole = positions.len();
        positions.push(Point3::new(0.0, 0.0, 1.0));
        let mut faces = Vec::new();
        for j in 0..ring {
            let next = (j + 1) % ring;
            faces.push(vec![j, next, next + ring, j + ring]);
            faces.push(vec![j + ring, next + ring, pole]);
        }
        faces.push((0..ring).rev().collect());
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let mut supports = evidence.source_planar_support[..2 * ring].to_vec();
        supports.sort_unstable();
        supports.dedup();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_y(),
            radial_center: 0.0,
            profile_radius: 1.0,
        };
        let proposal = index
            .measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: 1.0e-8,
                    chord: 0.5,
                },
            )
            .ok_or("sphere")?;
        let partition = index.select_surface_partition(&[proposal], &evidence, 100)?;
        let trims = index.assemble_shared_trims(&partition, &evidence.positions, trim_options())?;
        let native = index.construct_native_trims(
            &partition,
            &trims,
            &evidence.positions,
            native_options(),
        )?;
        assert_eq!(native.regions.len(), 2);
        assert!(native.regions.iter().all(|r| !r.requires_chart_cut));
        assert_eq!(
            native.regions.iter().map(|r| r.chart_seams).sum::<usize>(),
            2
        );
        let mut shell = crate::Shell::new();
        for region in native.regions {
            shell.extend(region.faces?);
        }
        assert_eq!(shell.len(), 3);
        assert_eq!(
            shell.shell_condition(),
            monstertruck_topology::shell::ShellCondition::Closed
        );
        Ok(())
    }
}
