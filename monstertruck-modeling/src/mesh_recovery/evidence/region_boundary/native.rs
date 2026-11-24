//! Native shared edges and provisional kernel faces. The source contact
//! graph supplies topological identity; geometry is never independently
//! rebuilt by each incident face. Periodic-domain cuts are separate from
//! the geometric consistency of the native edge and surface carriers.

use super::super::super::sample_fit::normal_basis;
use super::{
    MeshEvidence, RegionBoundaryIndex, SharedAnalyticTrim, SharedTrimGeometry, SharedTrimGraph,
    SurfaceCarrier, SurfaceRegionPartition,
    reconcile::{ReconciliationOptions, reconcile_point},
};
use crate::{
    Curve, Edge, Face, InnerSpace, Invertible, Line, Matrix4, ParametricSurface,
    ParametricSurface3D, Plane, Point3, Processor, RevolutionSurface, Sphere, Surface,
    ToSameGeometry, Torus, TrimmedCurve, UnitCircle, Vector3, Vertex, Wire, builder,
};
use std::f64::consts::{PI, TAU};

mod band;

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
        let (vertices, max_vertex_displacement) =
            native_vertices(partition, graph, positions, &intervals, options)?;
        let mut trims = graph
            .trims
            .iter()
            .zip(&intervals)
            .map(|(trim, intervals)| {
                intervals
                    .as_ref()
                    .map_err(|&error| error)?
                    .iter()
                    .map(|&interval| native_edge(trim, interval, &vertices, options))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Vec<_>>();
        let mut groups = vec![Vec::new(); partition.carriers.len()];
        for (id, &region) in partition.region_by_support.iter().enumerate() {
            groups
                .get_mut(region)
                .ok_or(NativeTrimError::InvalidSource)?
                .push(id);
        }
        let plans = graph
            .regions
            .iter()
            .map(|region| {
                let carrier = partition.carriers.get(region.region).copied().flatten()?;
                let supports = groups.get(region.region)?;
                let boundary = region.boundary.as_ref().ok()?;
                let carrier = chart_carrier(
                    carrier,
                    boundary,
                    graph,
                    self.evidence,
                    supports,
                    positions,
                    options,
                );
                let (surface, outward) =
                    native_surface(carrier, self.evidence, supports, positions).ok()?;
                if !boundary.pinched_vertices.is_empty() {
                    return None;
                }
                band::BandContext {
                    carrier,
                    surface: &surface,
                    outward,
                    evidence: self.evidence,
                    supports,
                    positions,
                    graph,
                    options,
                }
                .plan(boundary, &mut trims)
            })
            .collect::<Vec<_>>();
        let assembly = FaceAssembly {
            evidence: self.evidence,
            partition,
            graph,
            positions,
            options,
            trims: &trims,
            groups: &groups,
        };
        let regions = graph
            .regions
            .iter()
            .zip(&plans)
            .map(|(region, plan)| assembly.region(region, plan.as_ref()))
            .collect();
        Ok(NativeTrimGraph {
            vertices,
            trims,
            regions,
            max_vertex_displacement,
        })
    }
}

struct FaceAssembly<'a> {
    evidence: &'a MeshEvidence,
    partition: &'a SurfaceRegionPartition,
    graph: &'a SharedTrimGraph,
    positions: &'a [Point3],
    options: NativeTrimOptions,
    trims: &'a [Result<Vec<Edge>, NativeTrimError>],
    groups: &'a [Vec<usize>],
}
impl FaceAssembly<'_> {
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
                self.graph,
                self.evidence,
                supports,
                self.positions,
                self.options,
            );
            let (surface, outward) =
                native_surface(carrier, self.evidence, supports, self.positions)?;
            let mut wires = Vec::new();
            let mut areas = Vec::new();
            for cycle in &boundary.loops {
                let mut edges = Vec::new();
                let mut source = Vec::new();
                for reference in cycle {
                    let trim = self
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
                requires_chart_cut |= area.is_none();
                areas.push(area.unwrap_or(0.0) * if outward { 1.0 } else { -1.0 });
                wires.push(Wire::from(edges));
            }
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
        SharedTrimGeometry::Unresolved => Vec::new(),
    }
}

fn native_vertices(
    partition: &SurfaceRegionPartition,
    graph: &SharedTrimGraph,
    positions: &[Point3],
    intervals: &[Result<Vec<TrimInterval>, NativeTrimError>],
    options: NativeTrimOptions,
) -> Result<(Vec<Option<Vertex>>, f64), NativeTrimError> {
    let mut constraints = vec![Vec::new(); positions.len()];
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
        }
    }
    let mut maximum = 0.0_f64;
    let vertices = positions
        .iter()
        .zip(constraints)
        .map(|(&source, carriers)| {
            if carriers.is_empty() {
                return None;
            }
            let point = reconcile_point(
                source,
                &carriers,
                ReconciliationOptions {
                    max_displacement: options.max_vertex_displacement,
                    parallel_angle_tolerance: options.angular_resolution,
                    intersection_tolerance: options.intersection_tolerance,
                    max_normal_angle: options.angular_resolution * 2.0,
                },
            )?;
            maximum = maximum.max((point - source).magnitude());
            Some(builder::vertex(point))
        })
        .collect();
    Ok((vertices, maximum))
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
    graph: &SharedTrimGraph,
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
    let Some(axis) = circular_rim_axis(axial_center, boundary, graph, options).or_else(|| {
        compact_sphere_axis(
            axial_center,
            boundary,
            graph,
            evidence,
            supports,
            points,
            options,
        )
    }) else {
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
    graph: &SharedTrimGraph,
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
        let trim = graph.trims.get(reference.trim)?;
        let deviation = match trim.geometry {
            SharedTrimGeometry::Line(line) => line.max_vertex_offset,
            SharedTrimGeometry::Circle(circle) => {
                let sagitta = circle.radius * (1.0 - (circle.max_step_radians * 0.5).cos());
                sagitta
                    + circle
                        .max_plane_residual
                        .hypot(circle.max_vertex_radial_residual)
            }
            SharedTrimGeometry::Unresolved => return None,
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
) -> Result<(Surface, bool), NativeTrimError> {
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
    let outward = normal.dot(support.normal) > 0.0;
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
            let profile = |z| origin + axis * z + u * slope.mul_add(z, radius_at_origin);
            if slope.mul_add(lo, radius_at_origin) < 0.0
                || slope.mul_add(hi, radius_at_origin) < 0.0
                || lo >= hi
            {
                return Err(NativeTrimError::UnsupportedCarrier);
            }
            let surface = Surface::Revolution(Processor::new(RevolutionSurface::by_revolution(
                Curve::Line(Line(profile(lo), profile(hi))),
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
