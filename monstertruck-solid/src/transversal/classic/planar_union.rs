//! Planar Boolean overlay retaining original analytic edges.
//!
//! This operates on noded exact-edge arrangements. Repeated boundary walks
//! and overlapping selected cells are unified by nonzero material occupancy,
//! never by parity or by discarding positional face indices.

use super::divide_face::create_parameter_boundary_on_surface;
use super::faces_classification::BooleanRegion;
use super::weld::node_selected_regions;
use super::{AltCurve, AltRegion, ShapeOpsCurve, ShapeOpsSurface};
use monstertruck_geometry::prelude::*;
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
    polygons: Vec<Vec<Vec<Point2>>>,
    edges: Vec<ProjectedEdge<AltCurve<C, S>>>,
    span: f64,
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
fn needs_exact_noding<C>(edges: &[ProjectedEdge<C>], span: f64) -> bool {
    let epsilon = span * 1.0e-10;
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
        position > 1.0e-6
            && position < 1.0 - 1.0e-6
            && (point - (start + direction * position)).magnitude2() < epsilon * epsilon
    };
    for (index, lhs) in edges.iter().enumerate() {
        for rhs in edges.iter().skip(index + 1) {
            let lhs_uv = &lhs.uv;
            let rhs_uv = &rhs.uv;
            for &(point, other) in &[
                (lhs_uv[0], rhs_uv),
                (lhs_uv[lhs_uv.len() - 1], rhs_uv),
                (rhs_uv[0], lhs_uv),
                (rhs_uv[rhs_uv.len() - 1], lhs_uv),
            ] {
                if other.windows(2).any(|w| on_interior(point, w[0], w[1])) {
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
                    let approximate = arrangement.surface.evaluate(uv.x, uv.y);
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

fn project_source<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    source: &[AltRegion<C, S>],
    tol: f64,
) -> Option<PlanarArrangement<C, S>> {
    let surface = source.first()?.surface();
    if !matches!(
        surface.try_into_analytic_surface_kind(),
        Some(AnalyticSurfaceKind::Plane(_))
    ) {
        return None;
    }
    let mut samples = HashMap::default();
    let mut polygons = Vec::new();
    let mut edges = HashMap::<EdgeId<AltCurve<C, S>>, ProjectedEdge<AltCurve<C, S>>>::default();
    let mut span_min = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut span_max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for region in source {
        let mut region_loops = Vec::new();
        for wire in region.absolute_boundaries() {
            let projected =
                create_parameter_boundary_on_surface(&surface, wire, None, &mut samples, tol)?;
            if projected.len() < 3 {
                return None;
            }
            for &point in &projected.0 {
                span_min.x = span_min.x.min(point.x);
                span_min.y = span_min.y.min(point.y);
                span_max.x = span_max.x.max(point.x);
                span_max.y = span_max.y.max(point.y);
            }
            region_loops.push(projected.0);
            for edge in wire {
                if edges.contains_key(&edge.id()) {
                    continue;
                }
                let absolute = edge.absolute_clone();
                let one = Wire::from(vec![absolute.clone()]);
                let projected =
                    create_parameter_boundary_on_surface(&surface, &one, None, &mut samples, tol)?;
                if projected.len() < 2 {
                    return None;
                }
                edges.insert(
                    edge.id(),
                    ProjectedEdge {
                        edge: absolute,
                        uv: projected.0,
                    },
                );
            }
        }
        polygons.push(region_loops);
    }
    let span = (span_max.x - span_min.x).max(span_max.y - span_min.y);
    if !span.is_finite() || span <= 0.0 {
        return None;
    }
    Some(PlanarArrangement {
        surface,
        polygons,
        edges: edges.into_values().collect(),
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
        for ratio in [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5] {
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
        [1.0e-7, 1.0e-6].into_iter().all(|ratio| {
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
    let canonical_regions = node_selected_regions(source, &[])?;
    let initial = project_source(&canonical_regions, tol)?;
    let arrangement = if needs_exact_noding(&initial.edges, initial.span) {
        let points = exact_node_points(&initial, tol)?;
        let noded = node_selected_regions(&canonical_regions, &points)?;
        project_source(&noded, tol)?
    } else {
        initial
    };
    let boundary = selected_boundary(&arrangement)?;
    // Crossings between edges entirely inside the union do not affect the
    // extracted boundary. Reject only unresolved noding on surviving edges.
    if needs_exact_noding(&boundary, arrangement.span) {
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
    lhs_normal.dot(rhs_normal) >= 1.0 - 1.0e-10
        && region_lies_on_plane(lhs, &rhs_plane, tol)
        && region_lies_on_plane(rhs, &lhs_plane, tol)
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

/// Exact planar positive-minus-cutter occupancy on an already-identified
/// same-facing support. Operands can be independently parameterized or have
/// partial edge overlaps. Intersections become shared exact curve nodes.
fn overlay_planar_subtract<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    positive: &[AltRegion<C, S>],
    cutter: &[AltRegion<C, S>],
    tol: f64,
) -> Option<Vec<AltRegion<C, S>>> {
    let joined = positive.iter().chain(cutter).cloned().collect::<Vec<_>>();
    let canonical = node_selected_regions(&joined, &[])?;
    let initial = project_source(&canonical, tol)?;
    let arrangement = if needs_exact_noding(&initial.edges, initial.span) {
        let points = exact_node_points(&initial, tol)?;
        let noded = node_selected_regions(&canonical, &points)?;
        project_source(&noded, tol)?
    } else {
        initial
    };
    let positive_count = positive.len();
    let membership = |point| {
        covered(&arrangement.polygons[..positive_count], point)
            && !covered(&arrangement.polygons[positive_count..], point)
    };
    let boundary = selected_boundary_with(&arrangement, membership)?;
    if needs_exact_noding(&boundary, arrangement.span) {
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
    tol: f64,
) -> Option<Vec<AltRegion<C, S>>> {
    if positive_count > regions.len() {
        return None;
    }
    let mut taken = vec![false; regions.len()];
    let mut additions = HashMap::<usize, Vec<AltRegion<C, S>>>::default();
    for first in 0..positive_count {
        if taken[first] {
            continue;
        }
        let positive_indices = (first..positive_count)
            .filter(|&i| !taken[i] && same_oriented_support(&regions[first], &regions[i], tol))
            .collect::<Vec<_>>();
        let cutter_indices = (positive_count..regions.len())
            .filter(|&i| {
                if taken[i] {
                    return false;
                }
                let mut original = regions[i].clone();
                original.invert();
                same_oriented_support(&regions[first], &original, tol)
            })
            .collect::<Vec<_>>();
        if cutter_indices.is_empty() {
            continue;
        }
        let positive = positive_indices
            .iter()
            .map(|&i| regions[i].clone())
            .collect::<Vec<_>>();
        let cutter = cutter_indices
            .iter()
            .map(|&i| regions[i].clone())
            .collect::<Vec<_>>();
        // Common same-facing source/cutter patches cancel, but a cutter
        // patch outside this particular positive source face can still be
        // legitimately exposed by the cut. Preserve both non-overlapping
        // remainders with their respective original face orientations.
        let mut remainder = overlay_planar_subtract(&positive, &cutter, tol)?;
        remainder.extend(overlay_planar_subtract(&cutter, &positive, tol)?);
        additions.insert(first, remainder);
        for i in positive_indices.into_iter().chain(cutter_indices) {
            taken[i] = true;
        }
    }
    let mut result = Vec::new();
    for (i, region) in regions.into_iter().enumerate() {
        if let Some(mut replacement) = additions.remove(&i) {
            result.append(&mut replacement);
        } else if !taken[i] {
            result.push(region);
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
        let remainder = overlay_difference_regions(vec![host, void_inverted], 1, 1.0e-6)
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
        let remainder = overlay_difference_regions(vec![host, void_inverted], 1, 1.0e-6)
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
        let result = overlay_difference_regions(vec![host, void_inverted], 1, 1.0e-6)
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
