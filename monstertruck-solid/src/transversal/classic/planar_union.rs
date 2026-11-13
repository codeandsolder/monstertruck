//! Planar Boolean overlay retaining original analytic edges.
//!
//! This operates on noded exact-edge arrangements. Repeated boundary walks
//! and overlapping selected cells are unified by nonzero material occupancy,
//! never by parity or by discarding positional face indices.

use super::divide_face::create_parameter_boundary_on_surface;
use super::faces_classification::BooleanRegion;
use super::loops_store::{
    certified_circular_cylinder, common_circular_extrusion_cylinder, common_ruled_extrusion_strip,
};
use super::tolerance::BooleanTolerance;
use super::weld::{node_selected_regions, point_lex_cmp};
use super::{AltCurve, AltRegion, ShapeOpsCurve, ShapeOpsSurface};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::PolylineCurve;
use monstertruck_topology::{Edge, EdgeId, VertexId, Wire};
use rustc_hash::FxHashMap as HashMap;
use std::collections::BTreeSet;

type OverlayLoop<C, S> = (Wire<Point3, AltCurve<C, S>>, Vec<Point2>);

struct ProjectedEdge<C> {
    edge: Edge<Point3, C>,
    uv: Vec<Point2>,
}

struct PlanarArrangement<C, S> {
    surface: S,
    cylinder: Option<CylinderChart>,
    polygons: Vec<Vec<Vec<Point2>>>,
    edges: Vec<ProjectedEdge<AltCurve<C, S>>>,
    span: f64,
}

/// Global metric chart for a certified rational circular extrusion.
///
/// Two independently parameterized rational arcs may inhabit the same exact
/// cylinder but have incompatible local NURBS parameter domains. Projecting
/// one arc through the other's local inverse map can invent winding and
/// intersections. The cylinder's axial coordinate + unwrapped arc length
/// supplies a single chart for the entire certified common carrier.
#[derive(Clone, Copy)]
struct CylinderChart {
    center: Point3,
    axis: Vector3,
    radius: f64,
    basis_u: Vector3,
    basis_v: Vector3,
    handedness: f64,
    seam: f64,
}

impl CylinderChart {
    fn new<S: ShapeOpsSurface>(surface: &S, first_point: Point3, tol: f64) -> Option<Self> {
        let rational = surface.try_into_homogeneous_bspline_surface()?;
        let (center, axis, radius, _) = certified_circular_cylinder(&rational)?;
        let reference = if axis.x.abs() < 0.8 {
            Vector3::unit_x()
        } else {
            Vector3::unit_y()
        };
        let basis_u = axis.cross(reference).normalize();
        let basis_v = axis.cross(basis_u).normalize();
        let axial = (first_point - center).dot(axis);
        let radial = first_point - center - axis * axial;
        if (radial.magnitude() - radius).abs() > tol {
            return None;
        }
        let uv = surface
            .search_parameter(first_point, None, 100)
            .or_else(|| surface.search_nearest_parameter(first_point, None, 100))?;
        let normal = surface.normal(uv.0, uv.1);
        let radial_alignment = normal.dot(radial.normalize()) / normal.magnitude();
        if !radial_alignment.is_finite()
            || radial_alignment.abs()
                < DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_normal_alignment_cosine
        {
            return None;
        }
        // With this right-handed axis basis, (axial, increasing angle) has
        // an inward normal. Choose chart orientation to agree with the
        // reference surface's raw normal, preserving outer/hole winding.
        let handedness = if radial_alignment < 0.0 { 1.0 } else { -1.0 };
        Some(Self {
            center,
            axis,
            radius,
            basis_u,
            basis_v,
            handedness,
            seam: 0.0,
        })
    }

    fn raw_phase(&self, point: Point3) -> f64 {
        let delta = point - self.center;
        (self.handedness * delta.dot(self.basis_v).atan2(delta.dot(self.basis_u)))
            .rem_euclid(std::f64::consts::TAU)
    }

    fn set_seam(&mut self, points: impl IntoIterator<Item = Point3>) -> Option<()> {
        let mut angles = points
            .into_iter()
            .map(|point| self.raw_phase(point))
            .collect::<Vec<_>>();
        if angles.is_empty() {
            return None;
        }
        angles.sort_by(f64::total_cmp);
        let mut largest = (f64::NEG_INFINITY, 0.0);
        for (i, &start) in angles.iter().enumerate() {
            let end = angles
                .get(i + 1)
                .copied()
                .unwrap_or(angles[0] + std::f64::consts::TAU);
            let gap = end - start;
            if gap > largest.0 {
                largest = (gap, start + gap * 0.5);
            }
        }
        self.seam = largest.1.rem_euclid(std::f64::consts::TAU);
        Some(())
    }

    fn project(&self, point: Point3, previous: Option<Point2>, tol: f64) -> Option<Point2> {
        let delta = point - self.center;
        let axial = delta.dot(self.axis);
        let radial = delta - self.axis * axial;
        if !axial.is_finite() || (radial.magnitude() - self.radius).abs() > tol {
            return None;
        }
        let phase = self.raw_phase(point);
        let mut angle = (phase - self.seam).rem_euclid(std::f64::consts::TAU) + self.seam;
        if let Some(previous) = previous {
            angle += std::f64::consts::TAU
                * ((previous.y / self.radius - angle) / std::f64::consts::TAU).round();
        }
        Some(Point2::new(axial, self.radius * angle))
    }

    fn evaluate(&self, uv: Point2) -> Point3 {
        let angle = self.handedness * uv.y / self.radius;
        self.center
            + self.axis * uv.x
            + (self.basis_u * angle.cos() + self.basis_v * angle.sin()) * self.radius
    }
}

fn project_cylinder_wire<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    chart: &CylinderChart,
    wire: &Wire<Point3, AltCurve<C, S>>,
    samples: &mut HashMap<EdgeId<AltCurve<C, S>>, PolylineCurve<Point3>>,
    tol: f64,
) -> Option<PolylineCurve<Point2>> {
    let mut projected = vec![chart.project(wire.front_vertex()?.point(), None, tol)?];
    for edge in wire {
        let poly = samples.entry(edge.id()).or_insert_with(|| {
            let curve = edge.curve();
            PolylineCurve(
                curve
                    .parameter_division(curve.range_tuple(), BooleanTolerance::new(tol).mesh_seed)
                    .1,
            )
        });
        let walk = if edge.orientation() {
            poly.iter().skip(1).copied().collect::<Vec<_>>()
        } else {
            poly.iter().rev().skip(1).copied().collect::<Vec<_>>()
        };
        for world in walk {
            projected.push(chart.project(world, projected.last().copied(), tol)?);
        }
    }
    Some(PolylineCurve(projected))
}

fn signed_area(poly: &[Point2]) -> f64 {
    poly.iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
        .map(|(start, end)| start.x.mul_add(end.y, -end.x * start.y))
        .sum::<f64>()
        * 0.5
}

fn winding(poly: &[Point2], point: Point2) -> isize {
    poly.iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
        .fold(0_isize, |count, (start, end)| {
            let side = (end.x - start.x)
                .mul_add(point.y - start.y, -(end.y - start.y) * (point.x - start.x));
            if start.y <= point.y && end.y > point.y && side > 0.0 {
                count + 1
            } else if end.y <= point.y && start.y > point.y && side < 0.0 {
                count - 1
            } else {
                count
            }
        })
}

fn covered(polygons: &[Vec<Vec<Point2>>], point: Point2) -> bool {
    polygons.iter().any(|region| {
        region
            .iter()
            .map(|poly| winding(poly, point))
            .sum::<isize>()
            != 0
    })
}

/// Return true when existing exact edge identities are insufficient to form
/// an already-noded arrangement. Such crossings must be split, not guessed.
fn needs_exact_noding<C>(edges: &[ProjectedEdge<C>], span: f64, tol: f64) -> bool {
    let epsilon = span * DEFAULT_NUMERICAL_TOLERANCES.boolean_planar_relative_span;
    let cross = |start: Point2, end: Point2, point: Point2| {
        let direction = end - start;
        direction
            .x
            .mul_add(point.y - start.y, -(direction.y * (point.x - start.x)))
    };
    let on_interior = |point: Point2, start: Point2, end: Point2| {
        let direction = end - start;
        let length2 = direction.magnitude2();
        if length2 <= epsilon * epsilon {
            return false;
        }
        let position = (point - start).dot(direction) / length2;
        position > DEFAULT_NUMERICAL_TOLERANCES.boolean_planar_end_fraction
            && position < 1.0 - DEFAULT_NUMERICAL_TOLERANCES.boolean_planar_end_fraction
            && (point - (start + direction * position)).magnitude2() < epsilon * epsilon
    };
    for (index, lhs) in edges.iter().enumerate() {
        for rhs in edges.iter().skip(index + 1) {
            let lhs_uv = &lhs.uv;
            let rhs_uv = &rhs.uv;
            for (uv_endpoint, world_endpoint, other_edge, other_uv) in [
                (lhs_uv[0], lhs.edge.front().point(), &rhs.edge, rhs_uv),
                (
                    lhs_uv[lhs_uv.len() - 1],
                    lhs.edge.back().point(),
                    &rhs.edge,
                    rhs_uv,
                ),
                (rhs_uv[0], rhs.edge.front().point(), &lhs.edge, lhs_uv),
                (
                    rhs_uv[rhs_uv.len() - 1],
                    rhs.edge.back().point(),
                    &lhs.edge,
                    lhs_uv,
                ),
            ] {
                // UV projection is approximate even when analytic edges
                // share a world endpoint. Ignore chart drift only when
                // spatial agreement meets this operation's requested
                // numerical length, not the global kernel epsilon.
                if world_endpoint.distance(other_edge.front().point()) <= tol
                    || world_endpoint.distance(other_edge.back().point()) <= tol
                {
                    continue;
                }
                if other_uv
                    .windows(2)
                    .any(|w| on_interior(uv_endpoint, w[0], w[1]))
                {
                    return true;
                }
            }
            for lhs_segment in lhs_uv.windows(2) {
                let (lhs_start, lhs_end) = (lhs_segment[0], lhs_segment[1]);
                for rhs_segment in rhs_uv.windows(2) {
                    let (rhs_start, rhs_end) = (rhs_segment[0], rhs_segment[1]);
                    let (cross0, cross1) = (
                        cross(lhs_start, lhs_end, rhs_start),
                        cross(lhs_start, lhs_end, rhs_end),
                    );
                    let (cross2, cross3) = (
                        cross(rhs_start, rhs_end, lhs_start),
                        cross(rhs_start, rhs_end, lhs_end),
                    );
                    let epsilon0 = (lhs_end - lhs_start).magnitude() * epsilon;
                    let epsilon1 = (rhs_end - rhs_start).magnitude() * epsilon;
                    if (cross0 > epsilon0 && cross1 < -epsilon0
                        || cross0 < -epsilon0 && cross1 > epsilon0)
                        && (cross2 > epsilon1 && cross3 < -epsilon1
                            || cross2 < -epsilon1 && cross3 > epsilon1)
                    {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// Refine a crossing of projected curve segments on the actual analytic
/// curves, using the projection only as a starting point for Newton's method.
/// The returned point is shared by both original curves, never a polyline
/// intersection supplied as replacement geometry.
fn exact_curve_crossing<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    lhs: &ProjectedEdge<AltCurve<C, S>>,
    rhs: &ProjectedEdge<AltCurve<C, S>>,
    approximate: Point3,
    tol: f64,
) -> Option<Point3> {
    let lhs_curve = lhs.edge.oriented_curve();
    let rhs_curve = rhs.edge.oriented_curve();
    let (lhs_min, lhs_max) = lhs_curve.range_tuple();
    let (rhs_min, rhs_max) = rhs_curve.range_tuple();
    let mut lhs_parameter = lhs_curve.search_nearest_parameter(approximate, None, 100)?;
    let mut rhs_parameter = rhs_curve.search_nearest_parameter(approximate, None, 100)?;
    for _ in 0..32 {
        let (lhs_point, rhs_point) = (
            lhs_curve.evaluate(lhs_parameter),
            rhs_curve.evaluate(rhs_parameter),
        );
        let gap = lhs_point - rhs_point;
        if gap.magnitude() <= tol {
            return Some(lhs_point.midpoint(rhs_point));
        }
        let (lhs_derivative, rhs_derivative) = (
            lhs_curve.derivative(lhs_parameter),
            rhs_curve.derivative(rhs_parameter),
        );
        let diagonal0 = lhs_derivative.dot(lhs_derivative);
        let diagonal1 = rhs_derivative.dot(rhs_derivative);
        let mixed = -lhs_derivative.dot(rhs_derivative);
        let determinant = diagonal0.mul_add(diagonal1, -(mixed * mixed));
        if !determinant.is_finite() || determinant <= f64::EPSILON * diagonal0 * diagonal1 {
            return None;
        }
        let gradient0 = -lhs_derivative.dot(gap);
        let gradient1 = rhs_derivative.dot(gap);
        let step0 = gradient0.mul_add(diagonal1, -(gradient1 * mixed)) / determinant;
        let step1 = diagonal0.mul_add(gradient1, -(mixed * gradient0)) / determinant;
        lhs_parameter = (lhs_parameter + step0).clamp(lhs_min, lhs_max);
        rhs_parameter = (rhs_parameter + step1).clamp(rhs_min, rhs_max);
    }
    None
}

/// Gather real curve intersections and existing endpoints lying in another
/// edge's interior. Every candidate is checked against the exact curve, so
/// near crossings of sampled polylines do not become spurious B-rep vertices.
fn exact_node_points<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    arrangement: &PlanarArrangement<C, S>,
    tol: f64,
) -> Option<Vec<Point3>> {
    let mut nodes = Vec::<Point3>::new();
    for (index, lhs) in arrangement.edges.iter().enumerate() {
        for rhs in arrangement.edges.iter().skip(index + 1) {
            for (point, curve) in [
                (lhs.edge.front().point(), rhs.edge.oriented_curve()),
                (lhs.edge.back().point(), rhs.edge.oriented_curve()),
                (rhs.edge.front().point(), lhs.edge.oriented_curve()),
                (rhs.edge.back().point(), lhs.edge.oriented_curve()),
            ] {
                let on_curve = curve
                    .search_parameter(point, None, 100)
                    .or_else(|| curve.search_nearest_parameter(point, None, 100))
                    .is_some_and(|parameter| curve.evaluate(parameter).distance(point) <= tol);
                if on_curve && !nodes.iter().any(|existing| existing.distance(point) <= tol) {
                    nodes.push(point);
                }
            }
            for first in lhs.uv.windows(2) {
                for second in rhs.uv.windows(2) {
                    let lhs_direction = first[1] - first[0];
                    let rhs_direction = second[1] - second[0];
                    let determinant = lhs_direction
                        .x
                        .mul_add(rhs_direction.y, -(lhs_direction.y * rhs_direction.x));
                    if determinant.abs() <= f64::EPSILON {
                        continue;
                    }
                    let offset = second[0] - first[0];
                    let first_fraction = offset
                        .x
                        .mul_add(rhs_direction.y, -(offset.y * rhs_direction.x))
                        / determinant;
                    let second_fraction = offset
                        .x
                        .mul_add(lhs_direction.y, -(offset.y * lhs_direction.x))
                        / determinant;
                    if !(0.0..=1.0).contains(&first_fraction)
                        || !(0.0..=1.0).contains(&second_fraction)
                    {
                        continue;
                    }
                    let uv = first[0] + lhs_direction * first_fraction;
                    let approximate = arrangement.cylinder.as_ref().map_or_else(
                        || arrangement.surface.evaluate(uv.x, uv.y),
                        |chart| chart.evaluate(uv),
                    );
                    let point = exact_curve_crossing(lhs, rhs, approximate, tol)?;
                    if !nodes.iter().any(|existing| existing.distance(point) <= tol) {
                        nodes.push(point);
                    }
                }
            }
        }
    }
    nodes.sort_by(|lhs, rhs| {
        lhs.x
            .total_cmp(&rhs.x)
            .then_with(|| lhs.y.total_cmp(&rhs.y))
            .then_with(|| lhs.z.total_cmp(&rhs.z))
    });
    (!nodes.is_empty()).then_some(nodes)
}

/// A closed topological wire must not acquire an artificial gap merely
/// because its first UV inverse projection had no continuation hint while
/// the last one did. Preserve actual periodic winding: a full parameter
/// period is a nonlocal jump, larger than the adjacent sampled UV steps.
fn reconcile_local_wire_uv_closure<S: ParametricSurface3D>(
    surface: &S,
    cylinder: Option<&CylinderChart>,
    world_vertex: Point3,
    uv: &mut [Point2],
) {
    if uv.len() < 3 {
        return;
    }
    let first = uv[0];
    let last_index = uv.len() - 1;
    let last = uv[last_index];
    let gap = (first - last).magnitude2();
    let first_step = (uv[1] - first).magnitude2();
    let last_step = (uv[last_index - 1] - last).magnitude2();
    // Reconcile only a local inverse-map discrepancy. A wire which winds
    // around a periodic surface can legitimately start/end at parameter
    // coordinates a whole period apart despite sharing its 3D vertex.
    if !gap.is_finite() || gap >= first_step || gap >= last_step {
        return;
    }
    let evaluate = |uv: Point2| {
        cylinder.map_or_else(|| surface.evaluate(uv.x, uv.y), |chart| chart.evaluate(uv))
    };
    let first_residual = evaluate(first).distance2(world_vertex);
    let last_residual = evaluate(last).distance2(world_vertex);
    let representative = if first_residual <= last_residual {
        first
    } else {
        last
    };
    uv[0] = representative;
    uv[last_index] = representative;
}

/// Share a surface-chart parameter for every occurrence of a welded 3D
/// vertex, including occurrences in independently imported source faces.
///
/// Distinct periodic-chart lifts are retained: an occurrence is synchronized
/// only with candidates closer than both local incident edge steps. Among
/// equivalent inverse-map solutions, choose the smaller world-space residual.
fn reconcile_shared_vertex_uv<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    source: &[AltRegion<C, S>],
    surface: &S,
    cylinder: Option<&CylinderChart>,
    samples: &HashMap<EdgeId<AltCurve<C, S>>, PolylineCurve<Point3>>,
    projected: &mut [Vec<Vec<Point2>>],
) -> Option<()> {
    struct Occurrence {
        region: usize,
        wire: usize,
        index: usize,
        uv: Point2,
        local_step2: f64,
        world: Point3,
    }

    let mut occurrences = HashMap::<VertexId<Point3>, Vec<Occurrence>>::default();
    for (region_index, (region, loops)) in source.iter().zip(projected.iter()).enumerate() {
        for (wire_index, (wire, poly)) in region.absolute_boundaries().iter().zip(loops).enumerate()
        {
            let mut current = 0;
            for edge in wire {
                let segment_count = samples.get(&edge.id())?.len().checked_sub(1)?;
                if segment_count == 0 {
                    return None;
                }
                let next = current + segment_count;
                for (vertex, index) in [(edge.front(), current), (edge.back(), next)] {
                    let uv = *poly.get(index)?;
                    let previous_step = (index > 0).then(|| (uv - poly[index - 1]).magnitude2());
                    let next_step = poly.get(index + 1).map(|next| (*next - uv).magnitude2());
                    let local_step2 = previous_step
                        .into_iter()
                        .chain(next_step)
                        .filter(|step| *step > 0.0 && step.is_finite())
                        .fold(f64::INFINITY, f64::min);
                    occurrences
                        .entry(vertex.id())
                        .or_default()
                        .push(Occurrence {
                            region: region_index,
                            wire: wire_index,
                            index,
                            uv,
                            local_step2,
                            world: vertex.point(),
                        });
                }
                current = next;
            }
            if current + 1 != poly.len() {
                return None;
            }
        }
    }
    for group in occurrences.values() {
        for point in group {
            let representative = group
                .iter()
                .filter(|candidate| {
                    let gap2 = (candidate.uv - point.uv).magnitude2();
                    gap2 == 0.0
                        || candidate.local_step2.is_finite()
                            && point.local_step2.is_finite()
                            && gap2 < candidate.local_step2
                            && gap2 < point.local_step2
                })
                .min_by(|lhs, rhs| {
                    let residual = |occ: &Occurrence| {
                        cylinder
                            .map_or_else(
                                || surface.evaluate(occ.uv.x, occ.uv.y),
                                |chart| chart.evaluate(occ.uv),
                            )
                            .distance2(occ.world)
                    };
                    residual(lhs)
                        .total_cmp(&residual(rhs))
                        .then_with(|| lhs.uv.x.total_cmp(&rhs.uv.x))
                        .then_with(|| lhs.uv.y.total_cmp(&rhs.uv.y))
                })?;
            projected[point.region][point.wire][point.index] = representative.uv;
        }
    }
    Some(())
}

fn project_source<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    source: &[AltRegion<C, S>],
    tol: f64,
    certified_curved: bool,
) -> Option<PlanarArrangement<C, S>> {
    let surface = source.first()?.surface();
    if !matches!(
        surface.try_into_analytic_surface_kind(),
        Some(AnalyticSurfaceKind::Plane(_))
    ) {
        if !certified_curved {
            return None;
        }
        // Project onto one exact chart only when every complete underlying
        // surface has an algebraic common-carrier certificate. The later
        // noding/occupancy phases operate on the original analytic edges.
        let chart = surface.try_into_homogeneous_bspline_surface()?;
        for region in source {
            let other = region.surface().try_into_homogeneous_bspline_surface()?;
            if !common_ruled_extrusion_strip(&chart, &other)
                && !common_circular_extrusion_cylinder(&chart, &other)
            {
                return None;
            }
        }
    }
    let mut samples = HashMap::default();
    let first = source
        .first()?
        .absolute_boundaries()
        .first()?
        .front_vertex()?
        .point();
    let mut cylinder = CylinderChart::new(&surface, first, tol);
    if let Some(chart) = cylinder.as_mut() {
        let mut world_samples = Vec::new();
        for region in source {
            for edge in region
                .absolute_boundaries()
                .iter()
                .flat_map(|wire| wire.iter())
            {
                let poly = samples.entry(edge.id()).or_insert_with(|| {
                    let curve = edge.curve();
                    PolylineCurve(
                        curve
                            .parameter_division(
                                curve.range_tuple(),
                                BooleanTolerance::new(tol).mesh_seed,
                            )
                            .1,
                    )
                });
                world_samples.extend(poly.iter().copied());
            }
        }
        chart.set_seam(world_samples)?;
    }
    let mut projected = Vec::with_capacity(source.len());
    for region in source {
        let mut region_loops = Vec::new();
        for wire in region.absolute_boundaries() {
            let mut loop_uv = if let Some(chart) = cylinder.as_ref() {
                project_cylinder_wire(chart, wire, &mut samples, tol)?
            } else {
                create_parameter_boundary_on_surface(&surface, wire, None, &mut samples, tol)?
            };
            if loop_uv.len() < 3 {
                return None;
            }
            let front = wire.front_vertex()?;
            if front.id() == wire.back_vertex()?.id() {
                reconcile_local_wire_uv_closure(
                    &surface,
                    cylinder.as_ref(),
                    front.point(),
                    &mut loop_uv.0,
                );
            }
            region_loops.push(loop_uv.0);
        }
        projected.push(region_loops);
    }
    reconcile_shared_vertex_uv(
        source,
        &surface,
        cylinder.as_ref(),
        &samples,
        &mut projected,
    )?;

    let mut polygons = Vec::new();
    let mut edges = HashMap::<EdgeId<AltCurve<C, S>>, ProjectedEdge<AltCurve<C, S>>>::default();
    let mut span_min = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut span_max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for (region, region_loops) in source.iter().zip(projected) {
        for (wire, poly) in region.absolute_boundaries().iter().zip(&region_loops) {
            for &point in poly {
                span_min.x = span_min.x.min(point.x);
                span_min.y = span_min.y.min(point.y);
                span_max.x = span_max.x.max(point.x);
                span_max.y = span_max.y.max(point.y);
            }
            // Every exact edge must use its parent wire's *normalized*
            // chart coordinates, rather than a separate inverse projection.
            let mut start = 0;
            for edge in wire {
                let segment_count = samples.get(&edge.id())?.len().checked_sub(1)?;
                let end = start + segment_count;
                let mut uv = poly.get(start..=end)?.to_vec();
                if !edge.orientation() {
                    uv.reverse();
                }
                edges.entry(edge.id()).or_insert_with(|| ProjectedEdge {
                    edge: edge.absolute_clone(),
                    uv,
                });
                start = end;
            }
            if start + 1 != poly.len() {
                return None;
            }
        }
        polygons.push(region_loops);
    }
    let span = (span_max.x - span_min.x).max(span_max.y - span_min.y);
    if !span.is_finite() || span <= 0.0 {
        return None;
    }
    // Edge IDs are allocation-specific; use a geometric traversal order
    // so intersection evidence never depends on hash-table iteration.
    let mut edges = edges.into_values().collect::<Vec<_>>();
    edges.sort_by(|a, b| {
        let (af, ab, bf, bb) = (
            a.edge.front().point(),
            a.edge.back().point(),
            b.edge.front().point(),
            b.edge.back().point(),
        );
        let curve_middle = |edge: &ProjectedEdge<AltCurve<C, S>>| {
            let curve = edge.edge.oriented_curve();
            let (t0, t1) = curve.range_tuple();
            curve.evaluate(f64::midpoint(t0, t1))
        };
        point_lex_cmp(&af, &bf)
            .then_with(|| point_lex_cmp(&ab, &bb))
            .then_with(|| point_lex_cmp(&curve_middle(a), &curve_middle(b)))
            .then_with(|| a.uv.len().cmp(&b.uv.len()))
    });
    Some(PlanarArrangement {
        surface,
        cylinder,
        polygons,
        edges,
        span,
    })
}

/// Probe both sides of each exact edge, retaining only the union boundary.
fn selected_boundary_with<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    arrangement: &PlanarArrangement<C, S>,
    membership: impl Fn(Point2) -> bool,
) -> Option<Vec<ProjectedEdge<AltCurve<C, S>>>> {
    let mut boundary = Vec::new();
    for item in &arrangement.edges {
        let middle = item.uv.len() / 2;
        let (start, end) = (item.uv[middle - 1], item.uv[middle]);
        let tangent = end - start;
        let length = tangent.magnitude();
        if length <= 0.0 {
            return None;
        }
        let left = Vector2::new(-tangent.y, tangent.x) / length;
        let center = start.midpoint(end);
        let mut edge_transition = None;
        for ratio in DEFAULT_NUMERICAL_TOLERANCES.boolean_planar_primary_probe_ratios {
            let offset = ratio * arrangement.span;
            let (left_covered, right_covered) = (
                membership(center + left * offset),
                membership(center - left * offset),
            );
            if left_covered != right_covered {
                edge_transition = Some((left_covered, right_covered));
                break;
            }
        }
        match edge_transition {
            Some((true, false)) => boundary.push(ProjectedEdge {
                edge: item.edge.clone(),
                uv: item.uv.clone(),
            }),
            Some((false, true)) => {
                boundary.push(ProjectedEdge {
                    edge: item.edge.inverse(),
                    uv: item.uv.iter().rev().copied().collect(),
                });
            }
            _ => {}
        }
    }
    Some(boundary)
}

fn selected_boundary<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    arrangement: &PlanarArrangement<C, S>,
) -> Option<Vec<ProjectedEdge<AltCurve<C, S>>>> {
    selected_boundary_with(arrangement, |point| covered(&arrangement.polygons, point))
        .filter(|boundary| !boundary.is_empty())
}

/// Reconstruct exact-edge loops by choosing clockwise exits at junctions.
fn trace_boundary<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    mut boundary: Vec<ProjectedEdge<AltCurve<C, S>>>,
) -> Option<Vec<OverlayLoop<C, S>>> {
    boundary.sort_by(|lhs, rhs| {
        let lhs_uv = &lhs.uv;
        let rhs_uv = &rhs.uv;
        lhs_uv[0]
            .x
            .total_cmp(&rhs_uv[0].x)
            .then_with(|| lhs_uv[0].y.total_cmp(&rhs_uv[0].y))
            .then_with(|| {
                lhs_uv[lhs_uv.len() - 1]
                    .x
                    .total_cmp(&rhs_uv[rhs_uv.len() - 1].x)
            })
            .then_with(|| {
                lhs_uv[lhs_uv.len() - 1]
                    .y
                    .total_cmp(&rhs_uv[rhs_uv.len() - 1].y)
            })
    });
    let mut outgoing = HashMap::<VertexId<Point3>, Vec<usize>>::default();
    for (index, item) in boundary.iter().enumerate() {
        outgoing
            .entry(item.edge.front().id())
            .or_default()
            .push(index);
    }
    let mut unseen = (0..boundary.len()).collect::<BTreeSet<_>>();
    let mut loops = Vec::new();
    while let Some(&first) = unseen.first() {
        let mut index = first;
        let start = boundary[first].edge.front().id();
        let (mut chain, mut points) = (Vec::new(), Vec::new());
        loop {
            if !unseen.remove(&index) {
                return None;
            }
            let item = &boundary[index];
            chain.push(item.edge.clone());
            points.extend(item.uv.iter().take(item.uv.len() - 1).copied());
            let back = item.edge.back().id();
            if back == start {
                break;
            }
            let candidates = outgoing.get(&back)?;
            let previous = item.uv[item.uv.len() - 1] - item.uv[item.uv.len() - 2];
            let reverse_angle = (-previous.y).atan2(-previous.x);
            let next = candidates
                .iter()
                .filter(|&&candidate| unseen.contains(&candidate))
                .min_by(|&&lhs, &&rhs| {
                    let angle = |candidate: usize| {
                        let uv = &boundary[candidate].uv;
                        let tangent = uv[1] - uv[0];
                        (reverse_angle - tangent.y.atan2(tangent.x))
                            .rem_euclid(std::f64::consts::TAU)
                    };
                    angle(lhs).total_cmp(&angle(rhs))
                })?;
            index = *next;
        }
        loops.push((Wire::from(chain), points));
    }
    Some(loops)
}

fn validates_occupancy_with<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    arrangement: &PlanarArrangement<C, S>,
    loops: &[OverlayLoop<C, S>],
    membership: impl Fn(Point2) -> bool,
) -> bool {
    let traced = |point| {
        loops
            .iter()
            .map(|(_, uv)| winding(uv, point))
            .sum::<isize>()
            != 0
    };
    arrangement.edges.iter().all(|item| {
        let middle = item.uv.len() / 2;
        let start = item.uv[middle - 1];
        let end = item.uv[middle];
        let tangent = end - start;
        let length = tangent.magnitude();
        if length <= 0.0 {
            return false;
        }
        let left = Vector2::new(-tangent.y, tangent.x) / length;
        let center = start.midpoint(end);
        DEFAULT_NUMERICAL_TOLERANCES
            .boolean_planar_secondary_probe_ratios
            .into_iter()
            .all(|ratio| {
                let offset = arrangement.span * ratio;
                let left_probe = center + left * offset;
                let right_probe = center - left * offset;
                membership(left_probe) == traced(left_probe)
                    && membership(right_probe) == traced(right_probe)
            })
    })
}

fn validates_occupancy<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    arrangement: &PlanarArrangement<C, S>,
    loops: &[OverlayLoop<C, S>],
) -> bool {
    validates_occupancy_with(arrangement, loops, |point| {
        covered(&arrangement.polygons, point)
    })
}

fn assemble_regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    loops: Vec<OverlayLoop<C, S>>,
    reference: &AltRegion<C, S>,
    surface: &S,
) -> Option<Vec<AltRegion<C, S>>> {
    let mut outer = Vec::<(Vec<Wire<Point3, AltCurve<C, S>>>, Vec<Point2>)>::new();
    let mut holes = Vec::new();
    for (wire, poly) in loops {
        if signed_area(&poly) > 0.0 {
            outer.push((vec![wire], poly));
        } else {
            holes.push((wire, poly));
        }
    }
    for (hole, uv) in holes {
        let point = *uv.first()?;
        let parent = outer
            .iter()
            .enumerate()
            .filter(|(_, (_, poly))| winding(poly, point) != 0)
            .min_by(|(_, lhs), (_, rhs)| signed_area(&lhs.1).total_cmp(&signed_area(&rhs.1)))
            .map(|(index, _)| index)?;
        outer[parent].0.push(hole);
    }
    Some(
        outer
            .into_iter()
            .enumerate()
            .map(|(key, (wires, _))| {
                BooleanRegion::new(
                    wires,
                    surface.clone(),
                    reference.orientation(),
                    key,
                    reference.source_face(),
                )
            })
            .collect(),
    )
}

/// Union one source face's selected patches without changing exact edge curves.
fn overlay_one_source<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    source: &[AltRegion<C, S>],
    tol: f64,
) -> Option<Vec<AltRegion<C, S>>> {
    let reference = source.first()?;
    // Separate STEP faces may encode coincident curves with independent IDs.
    // Canonicalize their exact geometric edge identities before computing
    // occupancy; otherwise two independently encoded copies would be emitted
    // as two exterior boundaries rather than one filled face.
    let canonical_regions = node_selected_regions(source, &[], BooleanTolerance::new(tol))?;
    let initial = project_source(&canonical_regions, tol, false)?;
    let arrangement = if needs_exact_noding(&initial.edges, initial.span, tol) {
        let points = exact_node_points(&initial, tol)?;
        let noded = node_selected_regions(&canonical_regions, &points, BooleanTolerance::new(tol))?;
        project_source(&noded, tol, false)?
    } else {
        initial
    };
    let boundary = selected_boundary(&arrangement)?;
    // Crossings between edges entirely inside the union do not affect the
    // extracted boundary. Reject only unresolved noding on surviving edges.
    if needs_exact_noding(&boundary, arrangement.span, tol) {
        return None;
    }
    let loops = trace_boundary(boundary)?;
    if !validates_occupancy(&arrangement, &loops) {
        return None;
    }
    assemble_regions(loops, reference, &arrangement.surface)
}

/// Verify coplanarity over the actual trimmed geometry, not just at a plane
/// origin. Nearly parallel supports may diverge by macroscopic distances
/// across a large face even when their normals look identical.
fn region_lies_on_plane<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    region: &AltRegion<C, S>,
    plane: &Plane,
    tol: f64,
) -> bool {
    let normal = plane.normal();
    let origin = plane.origin();
    region
        .absolute_boundaries()
        .iter()
        .flat_map(|wire| wire.iter())
        .all(|edge| {
            let curve = edge.oriented_curve();
            let (start, end) = curve.range_tuple();
            [start, f64::midpoint(start, end), end]
                .into_iter()
                .all(|parameter| {
                    let point = curve.evaluate(parameter);
                    (point - origin).dot(normal).abs() <= tol
                })
        })
}

/// Coplanar faces with the same oriented world-space normal share material
/// ownership even if STEP supplied distinct source faces or parameter frames.
///
/// Comparing carrier geometry instead of STEP entity IDs avoids retaining
/// multiple coincident caps; comparing the *oriented* normal prevents merging
/// the two opposite sides of a thin zero-thickness interface.
fn same_oriented_support<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    lhs: &AltRegion<C, S>,
    rhs: &AltRegion<C, S>,
    tol: f64,
) -> bool {
    let (Some(AnalyticSurfaceKind::Plane(lhs_plane)), Some(AnalyticSurfaceKind::Plane(rhs_plane))) = (
        lhs.surface().try_into_analytic_surface_kind(),
        rhs.surface().try_into_analytic_surface_kind(),
    ) else {
        return false;
    };
    let lhs_normal = if lhs.orientation() {
        lhs_plane.normal()
    } else {
        -lhs_plane.normal()
    };
    let rhs_normal = if rhs.orientation() {
        rhs_plane.normal()
    } else {
        -rhs_plane.normal()
    };
    lhs_normal.dot(rhs_normal)
        >= 1.0 - DEFAULT_NUMERICAL_TOLERANCES.boolean_planar_normal_cosine_slack
        && region_lies_on_plane(lhs, &rhs_plane, tol)
        && region_lies_on_plane(rhs, &lhs_plane, tol)
}

/// Certified common ruled carrier with the SAME original material-facing
/// normal. The carrier identity is an algebraic rational control-net/circle
/// certificate, while the evaluated normal is used ONLY to distinguish its
/// two possible orientations at a shared geometric point.
fn same_oriented_difference_support<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    lhs: &AltRegion<C, S>,
    rhs: &AltRegion<C, S>,
    tol: f64,
) -> bool {
    if same_oriented_support(lhs, rhs, tol) {
        return true;
    }
    let (Some(a), Some(b)) = (
        lhs.surface().try_into_homogeneous_bspline_surface(),
        rhs.surface().try_into_homogeneous_bspline_surface(),
    ) else {
        return false;
    };
    let cylindrical = common_circular_extrusion_cylinder(&a, &b);
    if !common_ruled_extrusion_strip(&a, &b) && !cylindrical {
        return false;
    }
    if cylindrical {
        let (Some((ca, aa, _, _)), Some((cb, ab, _, _))) = (
            certified_circular_cylinder(&a),
            certified_circular_cylinder(&b),
        ) else {
            return false;
        };
        let Some(lhs_sign) = oriented_cylinder_radial_sign(lhs, ca, aa, tol) else {
            return false;
        };
        let Some(rhs_sign) = oriented_cylinder_radial_sign(rhs, cb, ab, tol) else {
            return false;
        };
        return lhs_sign * rhs_sign > 0.0;
    }
    let Some(point) = rhs
        .absolute_boundaries()
        .first()
        .and_then(Wire::front_vertex)
        .map(|vertex| vertex.point())
    else {
        return false;
    };
    let mut lhs_surface = lhs.surface();
    let mut rhs_surface = rhs.surface();
    if !lhs.orientation() {
        lhs_surface.invert();
    }
    if !rhs.orientation() {
        rhs_surface.invert();
    }
    let Some(lhs_uv) = lhs_surface
        .search_parameter(point, None, 100)
        .or_else(|| lhs_surface.search_nearest_parameter(point, None, 100))
    else {
        return false;
    };
    let Some(rhs_uv) = rhs_surface
        .search_parameter(point, None, 100)
        .or_else(|| rhs_surface.search_nearest_parameter(point, None, 100))
    else {
        return false;
    };
    if lhs_surface.evaluate(lhs_uv.0, lhs_uv.1).distance(point) > tol
        || rhs_surface.evaluate(rhs_uv.0, rhs_uv.1).distance(point) > tol
    {
        return false;
    }
    lhs_surface
        .normal(lhs_uv.0, lhs_uv.1)
        .dot(rhs_surface.normal(rhs_uv.0, rhs_uv.1))
        > 1.0 - DEFAULT_NUMERICAL_TOLERANCES.boolean_planar_tangent_cosine_slack
}

fn oriented_cylinder_radial_sign<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    region: &AltRegion<C, S>,
    center: Point3,
    axis: Vector3,
    tol: f64,
) -> Option<f64> {
    let point = region
        .absolute_boundaries()
        .first()?
        .front_vertex()?
        .point();
    let surface = region.surface();
    let uv = surface
        .search_parameter(point, None, 100)
        .or_else(|| surface.search_nearest_parameter(point, None, 100))?;
    if surface.evaluate(uv.0, uv.1).distance(point) > tol {
        return None;
    }
    let radial = point - center - axis * (point - center).dot(axis);
    let normal = surface.normal(uv.0, uv.1);
    let length = normal.magnitude() * radial.magnitude();
    if !length.is_finite() || length <= 0.0 {
        return None;
    }
    let sign = normal.dot(radial) / length;
    if sign.abs() < DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_normal_alignment_cosine {
        return None;
    }
    Some(if region.orientation() { sign } else { -sign })
}

/// Union selected cells sharing oriented planar support, including distinct source faces.
pub(super) fn overlay_source_regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: Vec<AltRegion<C, S>>,
    tol: f64,
) -> Vec<AltRegion<C, S>> {
    let mut groups = Vec::<Vec<AltRegion<C, S>>>::new();
    for region in regions {
        let group = groups.iter().position(|group| {
            same_oriented_support(&group[0], &region, tol)
                || group[0].source_face() == region.source_face()
                    && group[0].orientation() == region.orientation()
        });
        if let Some(index) = group {
            groups[index].push(region);
        } else {
            groups.push(vec![region]);
        }
    }
    groups
        .into_iter()
        .flat_map(|group| {
            if group.len() > 1
                && let Some(union) = overlay_one_source(&group, tol)
            {
                return union;
            }
            group
        })
        .collect()
}

/// Exact positive-minus-cutter occupancy on a certified common surface
/// chart (planar or algebraically shared rational ruled). Operands can be independently parameterized or have
/// partial edge overlaps. Intersections become shared exact curve nodes.
fn overlay_planar_subtract<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    positive: &[AltRegion<C, S>],
    cutter: &[AltRegion<C, S>],
    tol: f64,
) -> Option<Vec<AltRegion<C, S>>> {
    let joined = positive.iter().chain(cutter).cloned().collect::<Vec<_>>();
    let canonical = node_selected_regions(&joined, &[], BooleanTolerance::new(tol))?;
    let initial = project_source(&canonical, tol, true)?;
    let arrangement = if needs_exact_noding(&initial.edges, initial.span, tol) {
        let points = exact_node_points(&initial, tol)?;
        let noded = node_selected_regions(&canonical, &points, BooleanTolerance::new(tol))?;
        project_source(&noded, tol, true)?
    } else {
        initial
    };
    let positive_count = positive.len();
    let membership = |point| {
        covered(&arrangement.polygons[..positive_count], point)
            && !covered(&arrangement.polygons[positive_count..], point)
    };
    let boundary = selected_boundary_with(&arrangement, membership)?;
    if needs_exact_noding(&boundary, arrangement.span, tol) {
        return None;
    }
    if boundary.is_empty() {
        return validates_occupancy_with(&arrangement, &[], membership).then_some(Vec::new());
    }
    let loops = trace_boundary(boundary)?;
    if !validates_occupancy_with(&arrangement, &loops, membership) {
        return None;
    }
    assemble_regions(loops, &positive[0], &arrangement.surface)
}

/// Remove selected cutter patches that share a plane and *original outward
/// facing* with selected positive patches. Without this overlay the two
/// operands each own a coincident planar wall, creating triply-used edges at
/// the neighboring curved wall. Merely removing those faces would erase the
/// positive wall instead of carving its correctly trimmed remainder.
///
/// The inverted cutter regions have the opposite final orientation to their
/// original outside-facing surfaces, so compare positive regions against
/// the *uninverted* cutter to identify the same-facing subtraction interface.
pub(super) fn overlay_difference_regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: Vec<AltRegion<C, S>>,
    positive_count: usize,
    original_positive: &[AltRegion<C, S>],
    original_cutter: &[AltRegion<C, S>],
    tol: f64,
) -> Option<Vec<AltRegion<C, S>>> {
    if positive_count > regions.len() {
        return None;
    }
    // The face divider may produce a non-simple planar walk, e.g. an end
    // cap with retraced zero-area digons. Its boundary still describes a
    // well-defined selected planar region. Canonicalize its actual material
    // occupancy with exact analytic edges *before* comparing with other
    // planar supports. Keep the operand-side partition through expansion.
    let mut normalized = Vec::new();
    let mut normalized_positive_count = 0;
    for (index, region) in regions.into_iter().enumerate() {
        let fragments = if region.try_face().is_none()
            && matches!(
                region.surface().try_into_analytic_surface_kind(),
                Some(AnalyticSurfaceKind::Plane(_))
            ) {
            overlay_one_source(std::slice::from_ref(&region), tol)?
        } else {
            vec![region]
        };
        if index < positive_count {
            normalized_positive_count += fragments.len();
        }
        normalized.extend(fragments);
    }
    let regions = normalized;
    let positive_count = normalized_positive_count;
    // A shared *infinite* support is not a common trimmed NURBS patch.
    // Every original source face owns its finite parameter domain. Combining
    // fragments from distinct source faces here and assigning the union to
    // the first fragment's surface can create a topologically closed shell
    // whose boundary curves are unrepresentable by that finite patch.
    //
    // Clip each original positive face independently. A complete original
    // cutter footprint is a geometric subtractor even when its own inverted
    // face was excluded by volume classification. Only the separately
    // selected inverted cutter fragments can survive as inward walls.
    let mut taken_positive = vec![false; positive_count];
    let mut result = Vec::new();
    for first in 0..positive_count {
        if taken_positive[first] {
            continue;
        }
        let positive_indices = (first..positive_count)
            .filter(|&index| {
                !taken_positive[index]
                    && regions[index].source_face() == regions[first].source_face()
                    && regions[index].orientation() == regions[first].orientation()
            })
            .collect::<Vec<_>>();
        let positive = positive_indices
            .iter()
            .map(|&index| regions[index].clone())
            .collect::<Vec<_>>();
        let original_clipping = original_cutter
            .iter()
            .filter(|region| same_oriented_difference_support(&regions[first], region, tol))
            .cloned()
            .collect::<Vec<_>>();
        if original_clipping.is_empty() {
            result.extend(positive);
        } else {
            result.extend(overlay_planar_subtract(&positive, &original_clipping, tol)?);
        }
        for index in positive_indices {
            taken_positive[index] = true;
        }
    }

    // An inverted cutter patch is kept only where the COMPLETE source
    // positive did not already own the same outward-facing surface.
    // This stage must not depend on preclassified positive remnants.
    // Distinct cutter source regions are never assigned to each other's
    // finite surface parameter domains.
    for cutter in regions.into_iter().skip(positive_count) {
        let mut original = cutter.clone();
        original.invert();
        let complete_positive = original_positive
            .iter()
            .filter(|region| same_oriented_difference_support(&original, region, tol))
            .cloned()
            .collect::<Vec<_>>();
        if complete_positive.is_empty() {
            result.push(cutter);
        } else {
            result.extend(overlay_planar_subtract(
                std::slice::from_ref(&cutter),
                &complete_positive,
                tol,
            )?);
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::super::Alternative;
    use super::*;
    use monstertruck_geometry::prelude::{Line, Plane};
    use monstertruck_modeling::{Curve, Surface};
    use monstertruck_topology::Vertex;

    type TestRegion = AltRegion<Curve, Surface>;
    type TestEdge = Edge<Point3, AltCurve<Curve, Surface>>;

    struct PlanarFixture {
        vertices: Vec<Vertex<Point3>>,
        edges: HashMap<(usize, usize), TestEdge>,
    }

    impl PlanarFixture {
        fn new(points: &[[f64; 2]]) -> Self {
            Self {
                vertices: points
                    .iter()
                    .map(|p| Vertex::new(Point3::new(p[0], p[1], 0.0)))
                    .collect(),
                edges: HashMap::default(),
            }
        }

        fn edge(&mut self, a: usize, b: usize) -> TestEdge {
            if let Some(edge) = self.edges.get(&(a, b)) {
                return edge.clone();
            }
            if let Some(edge) = self.edges.get(&(b, a)) {
                return edge.inverse();
            }
            let (start, end) = (&self.vertices[a], &self.vertices[b]);
            let curve = Alternative::FirstType(Curve::from(Line(start.point(), end.point())));
            let edge = Edge::new_unchecked(start, end, curve);
            self.edges.insert((a, b), edge.clone());
            edge
        }

        fn region(&mut self, path: &[usize]) -> TestRegion {
            self.region_on(path, Plane::xy(), true, 0)
        }

        fn region_on(
            &mut self,
            path: &[usize],
            plane: Plane,
            orientation: bool,
            source_face: usize,
        ) -> TestRegion {
            let wire = (0..path.len() - 1)
                .map(|i| self.edge(path[i], path[i + 1]))
                .collect::<Vec<_>>();
            BooleanRegion::new(
                vec![Wire::from(wire)],
                Surface::Plane(plane),
                orientation,
                0,
                source_face,
            )
        }
    }

    /// An analytic quadratic profile swept in x: every edge is an original
    /// line or whole rational spline, not a chord sampled from a mesh.
    fn ruled_difference_region(
        axial: (f64, f64),
        orientation: bool,
        source_face: usize,
    ) -> TestRegion {
        use monstertruck_geometry::prelude::{
            BsplineCurve, BsplineSurface, KnotVector, NurbsCurve, NurbsSurface,
        };
        let profile = [(0.0, 0.0), (0.5, 0.4), (1.0, 0.0)];
        let net = profile
            .into_iter()
            .map(|(y, z)| {
                vec![
                    Vector4::new(axial.0, y, z, 1.0),
                    Vector4::new(axial.1, y, z, 1.0),
                ]
            })
            .collect::<Vec<_>>();
        let vertices = [
            Vertex::new(Point3::new(axial.0, 0.0, 0.0)),
            Vertex::new(Point3::new(axial.1, 0.0, 0.0)),
            Vertex::new(Point3::new(axial.1, 1.0, 0.0)),
            Vertex::new(Point3::new(axial.0, 1.0, 0.0)),
        ];
        let line = |a: usize, b: usize| {
            Edge::new_unchecked(
                &vertices[a],
                &vertices[b],
                Alternative::FirstType(Curve::from(Line(vertices[a].point(), vertices[b].point()))),
            )
        };
        let profile_edge = |column: usize, a: usize, b: usize| {
            Edge::new_unchecked(
                &vertices[a],
                &vertices[b],
                Alternative::FirstType(Curve::NurbsCurve(NurbsCurve::new(
                    BsplineCurve::new_unchecked(
                        KnotVector::bezier_knot(2),
                        net.iter().map(|row| row[column]).collect(),
                    ),
                ))),
            )
        };
        let wire = Wire::from(vec![
            line(0, 1),
            profile_edge(1, 1, 2),
            line(2, 3),
            profile_edge(0, 0, 3).inverse(),
        ]);
        BooleanRegion::new(
            vec![wire],
            Surface::NurbsSurface(NurbsSurface::new(BsplineSurface::new_unchecked(
                (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
                net,
            ))),
            orientation,
            0,
            source_face,
        )
    }

    fn circular_arc_region(
        axial: (f64, f64),
        angles: (f64, f64),
        radius: f64,
        source_face: usize,
    ) -> TestRegion {
        let mid = f64::midpoint(angles.0, angles.1);
        let w = ((angles.1 - angles.0) * 0.5).cos();
        let pts = [
            (radius * angles.0.cos(), radius * angles.0.sin(), 1.0),
            (radius * mid.cos() / w, radius * mid.sin() / w, w),
            (radius * angles.1.cos(), radius * angles.1.sin(), 1.0),
        ];
        let net = pts
            .into_iter()
            .map(|(y, z, weight)| {
                vec![
                    Vector4::new(axial.0 * weight, y * weight, z * weight, weight),
                    Vector4::new(axial.1 * weight, y * weight, z * weight, weight),
                ]
            })
            .collect::<Vec<_>>();
        let point = |x: f64, angle: f64| Point3::new(x, radius * angle.cos(), radius * angle.sin());
        let vertices = [
            Vertex::new(point(axial.0, angles.0)),
            Vertex::new(point(axial.1, angles.0)),
            Vertex::new(point(axial.1, angles.1)),
            Vertex::new(point(axial.0, angles.1)),
        ];
        let line = |a: usize, b: usize| {
            Edge::new_unchecked(
                &vertices[a],
                &vertices[b],
                Alternative::FirstType(Curve::from(Line(vertices[a].point(), vertices[b].point()))),
            )
        };
        let arc = |column: usize, a: usize, b: usize| {
            Edge::new_unchecked(
                &vertices[a],
                &vertices[b],
                Alternative::FirstType(Curve::NurbsCurve(NurbsCurve::new(
                    BsplineCurve::new_unchecked(
                        KnotVector::bezier_knot(2),
                        net.iter().map(|row| row[column]).collect(),
                    ),
                ))),
            )
        };
        BooleanRegion::new(
            vec![Wire::from(vec![
                line(0, 1),
                arc(1, 1, 2),
                line(2, 3),
                arc(0, 0, 3).inverse(),
            ])],
            Surface::NurbsSurface(NurbsSurface::new(BsplineSurface::new_unchecked(
                (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
                net,
            ))),
            true,
            0,
            source_face,
        )
    }

    #[test]
    fn adjacent_circular_source_faces_keep_distinct_finite_surface_ownership() {
        let r = 0.6731;
        let a = circular_arc_region((0.0, 2.0), (0.0, std::f64::consts::FRAC_PI_4), r, 101);
        let b = circular_arc_region(
            (0.0, 2.0),
            (std::f64::consts::FRAC_PI_4, std::f64::consts::FRAC_PI_2),
            r,
            202,
        );
        // The cutter touches A only on the common arc boundary; its actual
        // material overlap is a notch in B. Both positives share the exact
        // infinite cylinder but NOT the finite source-arc parameter domain.
        let cutter = circular_arc_region(
            (0.5, 1.0),
            (std::f64::consts::FRAC_PI_4, std::f64::consts::PI / 3.0),
            r,
            303,
        );
        assert!(same_oriented_difference_support(&a, &b, 1.0e-6));
        let result = overlay_difference_regions(
            vec![a.clone(), b.clone()],
            2,
            &[a.clone(), b.clone()],
            &[cutter],
            1.0e-6,
        )
        .expect("each original positive face must retain its own finite cylindrical support");
        assert_eq!(result.len(), 2);
        assert!(result.iter().all(|region| region.try_face().is_some()));
        let area = |region: &TestRegion| {
            let arrangement = project_source(std::slice::from_ref(region), 1.0e-6, true)
                .expect("a source face remains a coherent circular parameter region");
            arrangement.polygons[0]
                .iter()
                .map(|wire| signed_area(wire))
                .sum::<f64>()
                .abs()
        };
        let face_a = result
            .iter()
            .find(|region| region.source_face() == 101)
            .expect("untouched adjacent cylinder arc retains its own source face");
        let face_b = result
            .iter()
            .find(|region| region.source_face() == 202)
            .expect("notched source face is not merged onto its neighboring arc");
        assert!(
            (area(face_a) - area(&a)).abs() < 2.0e-5,
            "a touching cutter must not remove finite surface coverage from the neighbor"
        );
        let removed = area(&b) - area(face_b);
        let expected = r * 0.5 * (std::f64::consts::PI / 3.0 - std::f64::consts::FRAC_PI_4);
        assert!(
            (removed - expected).abs() < 2.0e-5,
            "source B lost wrong area: removed={removed}, expected={expected}"
        );
    }

    #[test]
    fn rational_arcs_on_one_cylinder_use_a_shared_metric_chart() {
        let radius = 0.6731;
        let outer = circular_arc_region((0.0, 2.0), (0.0, std::f64::consts::FRAC_PI_2), radius, 19);
        let inner = circular_arc_region((0.5, 1.0), (0.0, std::f64::consts::FRAC_PI_4), radius, 22);
        assert!(same_oriented_difference_support(&outer, &inner, 1.0e-6));
        let source = [outer.clone(), inner.clone()];
        let chart = project_source(&source, 1.0e-6, true)
            .expect("different rational arc charts share the exact cylinder");
        assert!(chart.cylinder.is_some());
        let result = overlay_planar_subtract(&[outer], &[inner], 1.0e-6)
            .expect("circular common carrier has one exact clipped wire");
        assert_eq!(result.len(), 1);
        assert!(result[0].try_face().is_some());
        assert_eq!(result[0].absolute_boundaries().len(), 1);
        let regions = project_source(&result, 1.0e-6, true).unwrap();
        let actual = signed_area(&regions.polygons[0][0]).abs();
        let expected = radius * 7.0 * std::f64::consts::PI / 8.0;
        assert!(
            (actual - expected).abs() < 2.0e-5,
            "actual={actual},expected={expected}"
        );
    }

    #[test]
    fn projected_uv_t_junction_at_shared_world_vertex_does_not_request_noding() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]);
        let horizontal = fixture.edge(0, 1);
        let vertical = fixture.edge(0, 2);
        let shared = Point2::new(0.0, 0.62);
        let edges = [
            ProjectedEdge {
                edge: horizontal,
                uv: vec![shared, Point2::new(1.0, 0.62)],
            },
            ProjectedEdge {
                edge: vertical,
                uv: vec![
                    Point2::new(0.0, 0.6200002),
                    Point2::new(0.0, 0.6199998),
                    Point2::new(0.0, 0.5),
                ],
            },
        ];
        assert!(
            !needs_exact_noding(&edges, 1.0, 1.0e-6),
            "inverse-mapping drift at the SAME analytic 3D endpoint is not an interior cut"
        );
    }

    #[test]
    fn projected_uv_t_junction_at_distinct_world_vertex_still_requests_noding() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [2.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
        let edges = [
            ProjectedEdge {
                edge: fixture.edge(0, 1),
                uv: vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)],
            },
            ProjectedEdge {
                edge: fixture.edge(2, 3),
                uv: vec![Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)],
            },
        ];
        assert!(
            needs_exact_noding(&edges, 2.0, 1.0e-6),
            "a distinct 3D edge endpoint on an exact edge interior needs an analytic split"
        );
    }

    #[test]
    fn cylinder_chart_radial_validation_respects_operation_precision() {
        let region =
            circular_arc_region((0.0, 1.0), (0.0, std::f64::consts::FRAC_PI_2), 0.6731, 101);
        let surface = region.surface();
        let point = region
            .absolute_boundaries()
            .first()
            .and_then(Wire::front_vertex)
            .expect("source cylinder must have an initial vertex")
            .point();
        let displaced = point + Vector3::new(0.0, 5.0e-7, 0.0);
        assert!(
            CylinderChart::new(&surface, displaced, 1.0e-8).is_none(),
            "a 500 nm radial error is not a 10 nm certified chart point"
        );
        assert!(
            CylinderChart::new(&surface, displaced, 1.0e-6).is_some(),
            "an explicitly coarser tolerance may qualify the same chart point"
        );
    }

    #[test]
    fn projected_uv_endpoint_drift_requires_operation_local_world_agreement() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [1.0, 0.0], [5.0e-7, 0.0], [0.0, 1.0]]);
        let edges = [
            ProjectedEdge {
                edge: fixture.edge(0, 1),
                uv: vec![Point2::new(0.0, 0.62), Point2::new(1.0, 0.62)],
            },
            ProjectedEdge {
                edge: fixture.edge(2, 3),
                uv: vec![
                    Point2::new(0.0, 0.6200002),
                    Point2::new(0.0, 0.6199998),
                    Point2::new(0.0, 0.5),
                ],
            },
        ];
        assert!(
            !needs_exact_noding(&edges, 1.0, 1.0e-6),
            "a coarsely admitted shared endpoint does not need a second node"
        );
        assert!(
            needs_exact_noding(&edges, 1.0, 1.0e-8),
            "two endpoints 500 nm apart cannot share identity at 10 nm precision"
        );
    }

    #[test]
    fn closed_wire_projection_reconciles_local_inverse_drift_without_erasing_winding() {
        let surface = Plane::xy();
        let vertex = Point3::new(0.0, 0.62, 0.0);
        let mut locally_drifting = vec![
            Point2::new(2.8e-7, 0.62),
            Point2::new(0.25, 0.62),
            Point2::new(0.25, 1.0),
            Point2::new(0.0, 1.0),
            Point2::new(0.0, 0.62),
        ];
        reconcile_local_wire_uv_closure(&surface, None, vertex, &mut locally_drifting);
        assert_eq!(locally_drifting[0], locally_drifting[4]);
        assert_eq!(locally_drifting[0], Point2::new(0.0, 0.62));

        // A closed 3D wire wrapping an intrinsic periodic parameter remains
        // open in the universal-cover UV chart by its actual winding period.
        let mut winding = vec![
            Point2::new(0.0, 0.62),
            Point2::new(0.25, 0.62),
            Point2::new(0.5, 0.62),
            Point2::new(0.75, 0.62),
            Point2::new(1.0, 0.62),
        ];
        reconcile_local_wire_uv_closure(&surface, None, vertex, &mut winding);
        assert_eq!(winding[0], Point2::new(0.0, 0.62));
        assert_eq!(winding[4], Point2::new(1.0, 0.62));
    }

    #[test]
    fn equivalent_welded_vertices_share_a_coherent_chart_across_regions() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]]);
        let a = fixture.region(&[0, 1, 2, 3, 0]);
        let b = fixture.region(&[0, 1, 2, 3, 0]);
        let source = [a, b];
        let mut polylines = HashMap::default();
        let mut projected = source
            .iter()
            .map(|region| {
                region
                    .absolute_boundaries()
                    .iter()
                    .map(|wire| {
                        create_parameter_boundary_on_surface(
                            &region.surface(),
                            wire,
                            None,
                            &mut polylines,
                            1.0e-6,
                        )
                        .unwrap()
                        .0
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        // The same topological first vertex appears in both source regions,
        // but one independent inverse mapping wandered in its UV chart.
        projected[1][0][0].x = 2.8e-7;
        projected[1][0][4].x = 2.8e-7;
        reconcile_shared_vertex_uv(
            &source,
            &Surface::Plane(Plane::xy()),
            None,
            &polylines,
            &mut projected,
        )
        .expect("welded source vertices admit chart reconciliation");
        assert_eq!(projected[1][0][0], projected[0][0][0]);
        assert_eq!(projected[1][0][4], projected[0][0][4]);

        // Different lifts of one periodic world vertex must not be merged
        // merely because their welded topological ID agrees.
        let mut winding = projected.clone();
        for uv in &mut winding[1][0] {
            uv.x += 2.0;
        }
        reconcile_shared_vertex_uv(
            &source,
            &Surface::Plane(Plane::xy()),
            None,
            &polylines,
            &mut winding,
        )
        .expect("separate chart lifts are retained");
        assert_eq!(winding[1][0][0].x, projected[1][0][0].x + 2.0);
        assert_eq!(winding[1][0][4].x, projected[1][0][4].x + 2.0);
    }

    #[test]
    fn projected_edge_charts_follow_their_parent_wire_walk() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        // Arrange that the last wire edge is inverse to its stored edge.
        fixture.edge(0, 3);
        let planar = fixture.region(&[0, 1, 2, 3, 0]);
        let curved = ruled_difference_region((0.0, 2.0), true, 4);
        for region in [planar, curved] {
            let arrangement = project_source(std::slice::from_ref(&region), 1.0e-6, true)
                .expect("wire and edge curves must share one UV projection");
            let poly = &arrangement.polygons[0][0];
            let mut start = 0;
            for edge in &region.absolute_boundaries()[0] {
                let item = arrangement
                    .edges
                    .iter()
                    .find(|item| item.edge.id() == edge.id())
                    .expect("projected edge retained original identity");
                let mut uv = item.uv.clone();
                if !edge.orientation() {
                    uv.reverse();
                }
                assert_eq!(
                    &poly[start..start + uv.len()],
                    uv,
                    "every projected edge must be a contiguous part of its source wire"
                );
                start += uv.len() - 1;
            }
            assert_eq!(start + 1, poly.len());
        }
    }

    #[test]
    fn ruled_difference_masks_cutter_against_original_not_selected_remnant() {
        let original = ruled_difference_region((0.0, 2.0), true, 10);
        let selected = ruled_difference_region((0.0, 0.4), true, 10);
        let inverted_cutter = ruled_difference_region((0.7, 1.2), false, 7);
        let mut unflipped = inverted_cutter.clone();
        unflipped.invert();
        assert!(same_oriented_difference_support(
            &original, &unflipped, 1.0e-6
        ));
        let result = overlay_difference_regions(
            vec![selected, inverted_cutter],
            1,
            &[original],
            &[unflipped],
            1.0e-6,
        )
        .expect("exact common ruled patch must support original ownership");
        assert_eq!(
            result.len(),
            1,
            "the cutter inside original source cannot be resurrected"
        );
        assert!(
            result[0].try_face().is_some(),
            "the remaining curved face is valid"
        );
    }

    #[test]
    fn original_same_facing_cutter_trims_positive_even_when_not_selected() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [1.0, 3.0],
            [3.0, 3.0],
            [3.0, 4.0],
            [1.0, 4.0],
        ]);
        let host = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 17);
        let original_cutter = fixture.region_on(&[4, 5, 6, 7, 4], Plane::xy(), true, 3);
        let result =
            overlay_difference_regions(vec![host.clone()], 1, &[host], &[original_cutter], 1.0e-6)
                .expect("a coincident unselected cutter still clips the source positive face");
        assert_eq!(result.len(), 1);
        let face = result[0].try_face().expect("notched face remains valid");
        assert_eq!(face.absolute_boundaries().len(), 1);
        assert_eq!(face.absolute_boundaries()[0].len(), 8);
        let projected = project_source(&result, 1.0e-6, false).unwrap();
        assert!(
            (signed_area(&projected.polygons[0][0]) - 14.0).abs() < 1.0e-7,
            "the whole original coincident cutter footprint must carve a 2-square-unit notch"
        );
    }

    #[test]
    fn disjoint_unselected_original_cutter_does_not_remove_positive_material() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [4.0, 0.0],
            [6.0, 0.0],
            [6.0, 2.0],
            [4.0, 2.0],
        ]);
        let host = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 1);
        let original_cutter = fixture.region_on(&[4, 5, 6, 7, 4], Plane::xy(), true, 2);
        let result =
            overlay_difference_regions(vec![host.clone()], 1, &[host], &[original_cutter], 1.0e-6)
                .expect("source cutter on same carrier but away from material is harmless");
        assert_eq!(result.len(), 1);
        let poly = project_source(&result, 1.0e-6, false).unwrap();
        assert!((signed_area(&poly.polygons[0][0]) - 4.0).abs() < 1.0e-7);
    }

    #[test]
    fn planar_difference_normalizes_single_selected_non_simple_cap() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [4.0, 0.0], [4.0, 3.0], [0.0, 3.0]]);
        // The same exact horizontal edge is traversed twice with opposite
        // incidence: this is an arrangement walk, not a valid B-rep face.
        let cap = fixture.region_on(&[0, 1, 2, 3, 0, 1, 0], Plane::xy(), true, 0);
        assert!(
            cap.try_face().is_none(),
            "raw retraced cap must be non-simple"
        );
        let pieces = overlay_difference_regions(vec![cap.clone()], 1, &[cap], &[], 1.0e-6)
            .expect("exact planar material occupancy must normalize a cap walk");
        assert_eq!(pieces.len(), 1);
        let face = pieces[0]
            .try_face()
            .expect("normalized cap must be a valid face");
        assert_eq!(face.absolute_boundaries().len(), 1);
        assert_eq!(face.absolute_boundaries()[0].len(), 4);
        let mut samples = HashMap::default();
        let poly = create_parameter_boundary_on_surface(
            &pieces[0].surface(),
            &pieces[0].absolute_boundaries()[0],
            None,
            &mut samples,
            1.0e-6,
        )
        .expect("normalized cap retained source analytic edges");
        assert!((signed_area(&poly.0) - 12.0).abs() < 1.0e-6);
    }

    #[test]
    fn planar_difference_masks_cutter_by_original_positive_not_preclassified_remnant() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 3.0],
            [0.0, 3.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [1.0, 3.0],
            [3.0, 3.0],
            [3.0, 4.0],
            [1.0, 4.0],
        ]);
        let original = fixture.region_on(&[0, 1, 4, 5, 0], Plane::xy(), true, 1);
        // The classifier already removed y=3..4 from the selected
        // positive: comparing the cutter against this selected fragment
        // would incorrectly resurrect the discarded shared source wall.
        let selected = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 1);
        let cutter = fixture.region_on(&[6, 7, 8, 9, 6], Plane::xy(), false, 2);
        let mut original_cutter = cutter.clone();
        original_cutter.invert();
        let result = overlay_difference_regions(
            vec![selected, cutter],
            1,
            &[original],
            &[original_cutter],
            1e-6,
        )
        .expect("complete original support is needed for cutter ownership");
        assert_eq!(
            result.len(),
            1,
            "same-facing cutter inside original must cancel"
        );
        let face = result[0]
            .try_face()
            .expect("selected positive remnant remains intact");
        assert_eq!(face.absolute_boundaries().len(), 1);
        let mut samples = HashMap::default();
        let poly = create_parameter_boundary_on_surface(
            &result[0].surface(),
            &result[0].absolute_boundaries()[0],
            None,
            &mut samples,
            1e-6,
        )
        .expect("remnant retains exact source boundary");
        assert!((signed_area(&poly.0) - 12.0).abs() < 1e-6);
    }

    #[test]
    fn planar_difference_carves_partial_same_facing_wall_without_losing_remainder() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [1.0, 3.0],
            [3.0, 3.0],
            [3.0, 4.0],
            [1.0, 4.0],
        ]);
        let host = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 17);
        // A same-original-facing cutter has opposite orientation after
        // inversion to become a Boolean difference face.
        let void_inverted = fixture.region_on(&[4, 5, 6, 7, 4], Plane::xy(), false, 3);
        let mut original_cutter = void_inverted.clone();
        original_cutter.invert();
        let remainder = overlay_difference_regions(
            vec![host.clone(), void_inverted],
            1,
            &[host],
            &[original_cutter],
            1.0e-6,
        )
        .expect("shared planar wall must support exact subtraction overlay");
        assert_eq!(remainder.len(), 1, "partial notch leaves one host face");
        let face = remainder[0]
            .try_face()
            .expect("notched wall remains a face");
        assert_eq!(face.absolute_boundaries().len(), 1);
        assert_eq!(face.absolute_boundaries()[0].len(), 8);
        let mut samples = HashMap::default();
        let projected = create_parameter_boundary_on_surface(
            &remainder[0].surface(),
            &remainder[0].absolute_boundaries()[0],
            None,
            &mut samples,
            1.0e-6,
        )
        .expect("analytic boundary must project");
        assert!((signed_area(&projected.0) - 14.0).abs() < 1.0e-7);
    }

    #[test]
    fn planar_difference_removes_fully_covered_common_wall() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let host = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 1);
        let void_inverted = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), false, 2);
        let mut original_cutter = void_inverted.clone();
        original_cutter.invert();
        let remainder = overlay_difference_regions(
            vec![host.clone(), void_inverted],
            1,
            &[host],
            &[original_cutter],
            1.0e-6,
        )
        .expect("fully coincident source wall should be consumed exactly");
        assert!(
            remainder.is_empty(),
            "identical same-facing contact patches cancel"
        );
    }

    #[test]
    fn planar_difference_does_not_remove_unrelated_or_thin_separate_wall() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [4.0, 0.0],
            [6.0, 0.0],
            [6.0, 2.0],
            [4.0, 2.0],
        ]);
        let host = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 1);
        let void_inverted = fixture.region_on(&[4, 5, 6, 7, 4], Plane::xy(), false, 2);
        let mut original_cutter = void_inverted.clone();
        original_cutter.invert();
        let result = overlay_difference_regions(
            vec![host.clone(), void_inverted],
            1,
            &[host],
            &[original_cutter],
            1.0e-6,
        )
        .expect("disjoint surface footprints remain usable");
        assert_eq!(result.len(), 2, "do not drop a disjoint cutter face");
    }

    #[test]
    fn different_source_faces_on_same_plane_share_union_boundary() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, -1.0],
            [3.0, -1.0],
            [3.0, 1.0],
            [1.0, 1.0],
        ]);
        let first = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 11);
        let transformed = Plane::new(
            Point3::new(3.0, 7.0, 0.0),
            Point3::new(3.0, 8.0, 0.0),
            Point3::new(2.0, 7.0, 0.0),
        );
        let second = fixture.region_on(&[4, 5, 6, 7, 4], transformed, true, 22);
        let union = overlay_source_regions(vec![first, second], 1.0e-7);
        assert_eq!(union.len(), 1);
        assert_eq!(union[0].absolute_boundaries()[0].len(), 8);
        let mut samples = HashMap::default();
        let projected = create_parameter_boundary_on_surface(
            &union[0].surface(),
            &union[0].absolute_boundaries()[0],
            None,
            &mut samples,
            1.0e-7,
        )
        .unwrap();
        assert!((signed_area(&projected.0) - 7.0).abs() < 1.0e-8);
    }

    #[test]
    fn reversed_parameter_orientation_with_flipped_face_normal_still_unifies() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, -1.0],
            [3.0, -1.0],
            [3.0, 1.0],
            [1.0, 1.0],
        ]);
        let first = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 11);
        let flipped = Plane::new(
            Point3::new(-2.0, 4.0, 0.0),
            Point3::new(-1.0, 4.0, 0.0),
            Point3::new(-2.0, 3.0, 0.0),
        );
        let second = fixture.region_on(&[4, 5, 6, 7, 4], flipped, false, 22);
        let union = overlay_source_regions(vec![first, second], 1.0e-7);
        assert_eq!(union.len(), 1);
        assert_eq!(union[0].absolute_boundaries()[0].len(), 8);
        let mut samples = HashMap::default();
        let projected = create_parameter_boundary_on_surface(
            &union[0].surface(),
            &union[0].absolute_boundaries()[0],
            None,
            &mut samples,
            1.0e-7,
        )
        .unwrap();
        assert!((signed_area(&projected.0) - 7.0).abs() < 1.0e-8);
    }

    #[test]
    fn opposite_oriented_plane_faces_do_not_merge() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let first = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 11);
        let second = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), false, 22);
        let union = overlay_source_regions(vec![first, second], 1.0e-7);
        assert_eq!(union.len(), 2);
    }

    #[test]
    fn duplicate_independent_face_topologies_unify_once() {
        // Two source shells often encode coincident faces with entirely
        // independent vertices/edges. Equality of STEP IDs is irrelevant.
        let points = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let mut first_fixture = PlanarFixture::new(&points);
        let mut second_fixture = PlanarFixture::new(&points);
        let first = first_fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 11);
        let second = second_fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 22);
        let union = overlay_source_regions(vec![first, second], 1.0e-7);
        assert_eq!(
            union.len(),
            1,
            "two coincident independent caps share coverage"
        );
        assert_eq!(union[0].absolute_boundaries().len(), 1);
        assert_eq!(union[0].absolute_boundaries()[0].len(), 4);
    }

    #[test]
    fn nearly_parallel_planes_are_not_coplanar_across_large_faces() {
        let points = [[0.0, 0.0], [10_000.0, 0.0], [10_000.0, 2.0], [0.0, 2.0]];
        let mut flat_fixture = PlanarFixture::new(&points);
        let mut tilted_fixture = PlanarFixture::new(&points);
        // Inclination of 1 µrad passes a naive normal-dot threshold and
        // coincident origins, but produces 10 mm separation at the far end.
        tilted_fixture.vertices = points
            .iter()
            .map(|p| Vertex::new(Point3::new(p[0], p[1], p[0] * 1.0e-6)))
            .collect();
        let flat = flat_fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 11);
        let tilted_plane = Plane::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 1.0e-6),
            Point3::new(0.0, 1.0, 0.0),
        );
        let tilted = tilted_fixture.region_on(&[0, 1, 2, 3, 0], tilted_plane, true, 22);
        assert!(!same_oriented_support(&flat, &tilted, 1.0e-6));
        assert_eq!(overlay_source_regions(vec![flat, tilted], 1.0e-6).len(), 2);
    }

    #[test]
    fn parallel_offset_planes_do_not_merge() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]);
        let first = fixture.region_on(&[0, 1, 2, 3, 0], Plane::xy(), true, 11);
        let offset = Plane::new(
            Point3::new(0.0, 0.0, 0.01),
            Point3::new(1.0, 0.0, 0.01),
            Point3::new(0.0, 1.0, 0.01),
        );
        let second = fixture.region_on(&[0, 1, 2, 3, 0], offset, true, 22);
        let union = overlay_source_regions(vec![first, second], 1.0e-7);
        assert_eq!(union.len(), 2);
    }

    #[test]
    fn overlapping_same_source_regions_are_a_nonzero_union() {
        // An outer six-edge square and its small corner triangle are present
        // as both a repeated subwalk and a separately selected positive cell.
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [0.0, 1.0],
        ]);
        let repeated = fixture.region(&[0, 1, 5, 0, 1, 2, 3, 4, 5, 0]);
        let triangle = fixture.region(&[0, 1, 5, 0]);
        let union = overlay_source_regions(vec![repeated, triangle], 1.0e-6);
        assert_eq!(union.len(), 1);
        assert_eq!(union[0].absolute_boundaries().len(), 1);
        let wire = &union[0].absolute_boundaries()[0];
        assert_eq!(wire.len(), 6);
        assert!(!wire.iter().any(|edge| edge.id() == fixture.edge(1, 5).id()));
        let mut polys = HashMap::default();
        let projected = create_parameter_boundary_on_surface(
            &union[0].surface(),
            wire,
            None,
            &mut polys,
            1.0e-6,
        )
        .unwrap();
        assert!((signed_area(&projected.0) - 16.0).abs() < 1.0e-8);
    }

    #[test]
    fn adjacent_coplanar_cells_drop_shared_internal_edge() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ]);
        let left = fixture.region(&[0, 1, 2, 3, 0]);
        let right = fixture.region(&[1, 4, 5, 2, 1]);
        let shared = fixture.edge(1, 2).id();
        let union = overlay_source_regions(vec![left, right], 1.0e-6);
        assert_eq!(union.len(), 1);
        let wire = &union[0].absolute_boundaries()[0];
        assert_eq!(wire.len(), 6);
        assert!(!wire.iter().any(|edge| edge.id() == shared));
    }

    #[test]
    fn crossing_rectangles_are_noded_and_union_exactly() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, -1.0],
            [3.0, -1.0],
            [3.0, 1.0],
            [1.0, 1.0],
        ]);
        let first = fixture.region(&[0, 1, 2, 3, 0]);
        let second = fixture.region(&[4, 5, 6, 7, 4]);
        let union = overlay_source_regions(vec![first, second], 1.0e-7);
        assert_eq!(
            union.len(),
            1,
            "proper crossings must be split into exact edges"
        );
        assert_eq!(union[0].absolute_boundaries().len(), 1);
        let wire = &union[0].absolute_boundaries()[0];
        assert_eq!(wire.len(), 8);
        let mut samples = HashMap::default();
        let poly = create_parameter_boundary_on_surface(
            &union[0].surface(),
            wire,
            None,
            &mut samples,
            1.0e-7,
        )
        .unwrap();
        assert!((signed_area(&poly.0) - 7.0).abs() < 1.0e-8);
    }

    #[test]
    fn endpoint_interior_t_junction_is_noded_without_facet_edges() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 2.0],
            [0.0, 2.0],
            [3.0, 0.5],
            [4.0, 0.5],
            [4.0, 1.5],
            [3.0, 1.5],
        ]);
        let main = fixture.region(&[0, 1, 2, 3, 0]);
        let protrusion = fixture.region(&[4, 5, 6, 7, 4]);
        let union = overlay_source_regions(vec![main, protrusion], 1.0e-7);
        assert_eq!(
            union.len(),
            1,
            "exact T junctions must share split vertices"
        );
        assert_eq!(union[0].absolute_boundaries().len(), 1);
        let wire = &union[0].absolute_boundaries()[0];
        assert_eq!(wire.len(), 8);
        let mut samples = HashMap::default();
        let poly = create_parameter_boundary_on_surface(
            &union[0].surface(),
            wire,
            None,
            &mut samples,
            1.0e-7,
        )
        .unwrap();
        assert!((signed_area(&poly.0) - 7.0).abs() < 1.0e-8);
    }

    #[test]
    fn duplicate_positive_regions_do_not_cancel_by_parity() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let a = fixture.region(&[0, 1, 2, 3, 0]);
        let b = fixture.region(&[0, 1, 2, 3, 0]);
        let union = overlay_source_regions(vec![a, b], 1.0e-6);
        assert_eq!(union.len(), 1);
        assert_eq!(union[0].absolute_boundaries()[0].len(), 4);
    }
}
