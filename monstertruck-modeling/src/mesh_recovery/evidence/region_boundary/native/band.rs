//! Periodic isoparametric bands, independent of the surface primitive.
//! Two noncontractible parallel boundary cycles define an annular domain.
//! Cut at two parameter values, globally subdivide their shared ring edges,
//! and use the same two seam sequences in opposite chart-face orientations.

use super::super::{
    MeshEvidence, OrientedTrim, RegionTrimLoops, SharedTrimGeometry, SharedTrimGraph,
    SurfaceCarrier,
};
use super::normal_basis;
use super::{
    NativeTrimError, NativeTrimOptions, chart_point, native_circle_edge, periodic_coordinate,
};
use crate::{
    Edge, Face, IncludeCurve, InnerSpace, Point3, Surface, Vector3, Vertex, Wire, builder,
};
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

mod variable_rim;

#[derive(Clone, Copy)]
struct Ring {
    axis: usize,
    fixed: f64,
    positive: bool,
}

#[derive(Clone, Copy)]
struct Circle {
    center: Point3,
    normal: Vector3,
    radius: f64,
}

pub(super) struct BandPlan {
    vertices: [[Vertex; 2]; 2],
    seams: [Vec<Edge>; 2],
}

impl BandPlan {
    pub(super) fn faces(
        &self,
        wires: &[Wire],
        surface: &Surface,
    ) -> Result<Vec<Face>, NativeTrimError> {
        let [first, rest @ ..] = wires else {
            return Err(NativeTrimError::InvalidBoundary);
        };
        let [a0, a1] = &self.vertices[0];
        let [b0, b1] = &self.vertices[1];
        let [first_a, first_b] = split_wire(first, a0, a1).ok_or(NativeTrimError::FaceTopology)?;
        let [second_a, second_b] = match rest {
            [second] => split_wire(second, b1, b0).ok_or(NativeTrimError::FaceTopology)?,
            [] if b0 == b1 => [Vec::new(), Vec::new()],
            _ => return Err(NativeTrimError::InvalidBoundary),
        };
        let wire_a = first_a
            .into_iter()
            .chain(self.seams[1].iter().cloned())
            .chain(second_a)
            .chain(self.seams[0].iter().rev().map(Edge::inverse))
            .collect::<Wire>();
        let wire_b = first_b
            .into_iter()
            .chain(self.seams[0].iter().cloned())
            .chain(second_b)
            .chain(self.seams[1].iter().rev().map(Edge::inverse))
            .collect::<Wire>();
        [wire_a, wire_b]
            .into_iter()
            .map(|wire| {
                let face = Face::new(vec![wire], surface.clone())
                    .map_err(|_| NativeTrimError::FaceTopology)?;
                if !face.is_geometric_consistent() {
                    return Err(NativeTrimError::FaceGeometry);
                }
                Ok(face)
            })
            .collect()
    }
}

pub(super) struct BandContext<'a> {
    pub carrier: SurfaceCarrier,
    pub surface: &'a Surface,
    pub outward: bool,
    pub evidence: &'a MeshEvidence,
    pub supports: &'a [usize],
    pub positions: &'a [Point3],
    pub graph: &'a SharedTrimGraph,
    pub options: NativeTrimOptions,
}

impl BandContext<'_> {
    pub(super) fn plan(
        &self,
        boundary: &RegionTrimLoops,
        trims: &mut [Result<Vec<Edge>, NativeTrimError>],
    ) -> Option<BandPlan> {
        let [first, second] = boundary.loops.as_slice() else {
            return self.cap(boundary, trims);
        };
        let rings = match (self.ring(first), self.ring(second)) {
            (Some(a), Some(b)) => [a, b],
            (Some(ring), None) => return self.variable_band([first, second], 0, ring, trims),
            (None, Some(ring)) => return self.variable_band([first, second], 1, ring, trims),
            (None, None) => return None,
        };
        if rings[0].axis != rings[1].axis || rings[0].positive == rings[1].positive {
            return None;
        }
        let axis = rings[0].axis;
        let other = 1 - axis;
        // In an outward (azimuth,meridian/axial) chart, a positive
        // azimuth winding bounds the lower meridian coordinate; a positive
        // meridian winding bounds the upper azimuth coordinate.
        let first_lower = (rings[0].positive == (axis == 0)) == self.outward;
        let (lower, upper) = if first_lower {
            (rings[0].fixed, rings[1].fixed)
        } else {
            (rings[1].fixed, rings[0].fixed)
        };
        let span = if periodic_coordinate(self.carrier, other) {
            (upper - lower).rem_euclid(TAU)
        } else {
            upper - lower
        };
        if span <= self.coordinate_tolerance(other) {
            return None;
        }
        if !self.covers_band(other, lower, span) {
            return None;
        }
        let mut changed = copy_rings(first.iter().chain(second), trims)?;
        let mut endpoints = Vec::new();
        for (cycle, ring) in [first, second].into_iter().zip(rings) {
            endpoints.push(self.ring_vertices(cycle, ring, &mut changed)?);
        }
        let vertices: [[Vertex; 2]; 2] = endpoints.try_into().ok()?;
        let signed_span = if first_lower { span } else { -span };
        let seams = [0.0, PI]
            .into_iter()
            .enumerate()
            .map(|(id, angle)| {
                self.seam(
                    axis,
                    angle,
                    rings[0].fixed,
                    signed_span,
                    &vertices[0][id],
                    &vertices[1][id],
                )
            })
            .collect::<Option<Vec<_>>>()?
            .try_into()
            .ok()?;
        commit_rings(changed, trims);
        Some(BandPlan { vertices, seams })
    }

    /// A closed, regular, degree-one source chart covers the complete torus.
    /// Build four bounded faces with eight shared half-circle seams, not a
    /// boundaryless provisional face or a synthetic material solid.
    pub(super) fn closed_torus(&self) -> Option<Vec<Face>> {
        let SurfaceCarrier::CircularRevolution {
            radial_center,
            profile_radius,
            ..
        } = self.carrier
        else {
            return None;
        };
        if radial_center - profile_radius <= self.options.intersection_tolerance
            || profile_radius <= self.options.intersection_tolerance
            || !self.has_unit_torus_degree()
        {
            return None;
        }
        let angles = [0.0, PI];
        let mut vertices = Vec::with_capacity(4);
        for u in angles {
            for v in angles {
                vertices.push(builder::vertex(world_point(self.carrier, [u, v])?));
            }
        }
        let mut u_edges = Vec::with_capacity(4);
        let mut v_edges = Vec::with_capacity(4);
        for (i, u) in angles.into_iter().enumerate() {
            for (j, v) in angles.into_iter().enumerate() {
                let start = &vertices[2 * i + j];
                for (varying, fixed, end, edges) in [
                    (0, v, &vertices[2 * (1 - i) + j], &mut u_edges),
                    (1, u, &vertices[2 * i + 1 - j], &mut v_edges),
                ] {
                    let circle = parameter_circle(self.carrier, varying, fixed)?;
                    let edge = native_circle_edge(
                        start,
                        end,
                        circle.center,
                        circle.normal,
                        circle.radius,
                        PI,
                        self.options.intersection_tolerance,
                    )
                    .ok()?;
                    if !edge.is_geometric_consistent() || !self.surface.include(&edge.curve()) {
                        return None;
                    }
                    edges.push(edge);
                }
            }
        }
        let mut faces = Vec::with_capacity(4);
        for i in 0..2 {
            for j in 0..2 {
                let wire = Wire::from(vec![
                    u_edges[2 * i + j].clone(),
                    v_edges[2 * (1 - i) + j].clone(),
                    u_edges[2 * i + 1 - j].inverse(),
                    v_edges[2 * i + j].inverse(),
                ]);
                let wire = if self.outward { wire } else { wire.inverse() };
                let face = Face::new(vec![wire], self.surface.clone()).ok()?;
                if !face.is_geometric_consistent() {
                    return None;
                }
                faces.push(face);
            }
        }
        let shell = faces.iter().cloned().collect::<crate::Shell>();
        (shell.shell_condition() == monstertruck_topology::shell::ShellCondition::Closed
            && shell.singular_vertices().is_empty())
        .then_some(faces)
    }

    // The source polygons form an oriented closed 2-chain in the periodic
    // parameter domain. Check opposite incidence of every edge, unique
    // local lifts, and consistent nonnegative chart orientation. The total
    // signed area divided by (2*pi)^2 is its integer mapping degree.
    // Degree one excludes duplicated covers and a cover with a missing patch.
    fn has_unit_torus_degree(&self) -> bool {
        let mut incidence = BTreeMap::<[usize; 2], (usize, bool)>::new();
        let mut total = 0.0;
        let mut compensation = 0.0;
        let mut polygon_count = 0_usize;
        for &support in self.supports {
            for &polygon in &self.evidence.planar_supports[support].source_face_indices {
                let Some(ids) = self.evidence.source_polygon_vertices(polygon) else {
                    return false;
                };
                let Some(area) = self.oriented_chart_polygon_area(ids) else {
                    return false;
                };
                for (&a, &b) in ids.iter().zip(ids.iter().cycle().skip(1)) {
                    if a == b {
                        return false;
                    }
                    let key = [a.min(b), a.max(b)];
                    let entry = incidence.entry(key).or_insert((0, a < b));
                    if entry.0 >= 2 || (entry.0 == 1 && entry.1 == (a < b)) {
                        return false;
                    }
                    entry.0 += 1;
                }
                let corrected = area - compensation;
                let next = total + corrected;
                compensation = (next - total) - corrected;
                total = next;
                polygon_count += 1;
            }
        }
        // An integer-degree chain must have an area of exactly one period
        // square. Account only for accumulated floating point error.
        let roundoff = f64::EPSILON * 128.0 * crate::usize_to_f64(polygon_count) * TAU * TAU;
        polygon_count != 0
            && incidence.values().all(|&(count, _)| count == 2)
            && (total - TAU * TAU).abs() <= roundoff
    }

    fn oriented_chart_polygon_area(&self, ids: &[usize]) -> Option<f64> {
        if ids.len() < 3 {
            return None;
        }
        let coords = ids
            .iter()
            .map(|&id| {
                chart_point(
                    self.carrier,
                    *self.positions.get(id)?,
                    self.options.angular_resolution,
                )
            })
            .collect::<Option<Vec<_>>>()?;
        let origin = *coords.first()?;
        let mut previous = origin;
        let mut area = 0.0;
        let sign = if self.outward { 1.0 } else { -1.0 };
        for raw in coords.iter().skip(1).chain(coords.first()) {
            let mut next = *raw;
            for axis in 0..2 {
                let delta = coordinate_difference(self.carrier, axis, raw[axis], previous[axis]);
                if periodic_coordinate(self.carrier, axis)
                    && delta.abs() >= PI - self.options.angular_resolution
                {
                    return None;
                }
                next[axis] = previous[axis] + delta;
            }
            let a = [previous[0] - origin[0], previous[1] - origin[1]];
            let b = [next[0] - origin[0], next[1] - origin[1]];
            let twice = sign * a[0].mul_add(b[1], -a[1] * b[0]);
            let scale = a[0]
                .abs()
                .max(a[1].abs())
                .max(b[0].abs())
                .max(b[1].abs())
                .max(1.0);
            let roundoff = 32.0 * f64::EPSILON * scale * scale;
            if twice < -roundoff {
                return None;
            }
            area += twice;
            previous = next;
        }
        ((previous[0] - origin[0]).abs() <= 32.0 * f64::EPSILON * TAU
            && (previous[1] - origin[1]).abs() <= 32.0 * f64::EPSILON * origin[1].abs().max(TAU)
            && area > 0.0)
            .then_some(area * 0.5)
    }

    // A cap is the limiting annular chart whose second ring collapses
    // to a singular point. Both chart seams share one topological tip.
    fn cap(
        &self,
        boundary: &RegionTrimLoops,
        trims: &mut [Result<Vec<Edge>, NativeTrimError>],
    ) -> Option<BandPlan> {
        let [cycle] = boundary.loops.as_slice() else {
            return None;
        };
        let ring = self.ring(cycle)?;
        if ring.axis != 0 {
            return None;
        }
        let positive = ring.positive == self.outward;
        let (tip_coordinate, tip) = match self.carrier {
            SurfaceCarrier::CircularRevolution {
                axial_center,
                axis,
                radial_center,
                profile_radius,
            } if radial_center.abs() <= crate::KERNEL_EPSILON => {
                let direction = if positive { 1.0 } else { -1.0 };
                (
                    direction * PI * 0.5,
                    axial_center + axis * (direction * profile_radius),
                )
            }
            SurfaceCarrier::LinearRevolution {
                origin,
                axis,
                radius_at_origin,
                slope,
            } if slope != 0.0 => {
                let z = -radius_at_origin / slope;
                let tip = origin + axis * z;
                // A finite cone with one rim must actually reach its
                // singular end in the observed source, not extrapolate
                // across an unobserved portion of the infinite carrier.
                let witnessed = self.supports.iter().any(|&support| {
                    self.evidence.planar_supports[support]
                        .source_face_indices
                        .iter()
                        .any(|&polygon| {
                            self.evidence
                                .source_polygon_vertices(polygon)
                                .is_some_and(|ids| {
                                    ids.iter().any(|&id| {
                                        (self.positions[id] - tip).magnitude()
                                            <= self.options.max_vertex_displacement
                                    })
                                })
                        })
                });
                if !witnessed || (z > ring.fixed) != positive {
                    return None;
                }
                (z, tip)
            }
            _ => return None,
        };
        let sweep = tip_coordinate - ring.fixed;
        if sweep.abs() <= self.coordinate_tolerance(1)
            || !self.covers_band(1, tip_coordinate.min(ring.fixed), sweep.abs())
        {
            return None;
        }
        let mut changed = copy_rings(cycle.iter(), trims)?;
        let endpoints = self.ring_vertices(cycle, ring, &mut changed)?;
        let pole = builder::vertex(tip);
        let seams = [0.0, PI]
            .into_iter()
            .enumerate()
            .map(|(id, angle)| self.seam(0, angle, ring.fixed, sweep, &endpoints[id], &pole))
            .collect::<Option<Vec<_>>>()?
            .try_into()
            .ok()?;
        commit_rings(changed, trims);
        Some(BandPlan {
            vertices: [endpoints, [pole.clone(), pole]],
            seams,
        })
    }

    fn ring_vertices(
        &self,
        cycle: &[OrientedTrim],
        ring: Ring,
        changed: &mut BTreeMap<usize, Vec<Edge>>,
    ) -> Option<[Vertex; 2]> {
        let vertices = [0.0, PI].map(|angle| {
            let mut uv = [0.0; 2];
            uv[ring.axis] = angle;
            uv[1 - ring.axis] = ring.fixed;
            let point = world_point(self.carrier, uv)?;
            cut_ring(cycle, point, changed, self.options.intersection_tolerance)
        });
        let [a, b] = vertices;
        Some([a?, b?])
    }

    fn ring(&self, cycle: &[OrientedTrim]) -> Option<Ring> {
        let mut coords = Vec::new();
        for reference in cycle {
            let trim = self.graph.trims.get(reference.trim)?;
            if !matches!(trim.geometry, SharedTrimGeometry::Circle(_)) {
                return None;
            }
            if reference.reversed {
                for edge in trim.source_edges.iter().rev() {
                    coords.push(chart_point(
                        self.carrier,
                        self.positions[edge.vertices[1]],
                        self.options.angular_resolution,
                    )?);
                }
            } else {
                for edge in &trim.source_edges {
                    coords.push(chart_point(
                        self.carrier,
                        self.positions[edge.vertices[0]],
                        self.options.angular_resolution,
                    )?);
                }
            }
        }
        let first = *coords.first()?;
        for axis in 0..2 {
            if !periodic_coordinate(self.carrier, axis) {
                continue;
            }
            let other = 1 - axis;
            let fixed = first[other];
            if !coords.iter().all(|uv| {
                coordinate_difference(self.carrier, other, uv[other], fixed).abs()
                    <= self.coordinate_tolerance(other)
            }) {
                continue;
            }
            let sweep = coords
                .iter()
                .zip(coords.iter().cycle().skip(1))
                .map(|(a, b)| (b[axis] - a[axis] + PI).rem_euclid(TAU) - PI)
                .sum::<f64>();
            if (sweep.abs() - TAU).abs() > self.options.angular_resolution {
                continue;
            }
            let expected = parameter_circle(self.carrier, axis, fixed)?;
            // All constituent analytic trims must be this same complete
            // isoparametric circle, not merely interpolate ring samples.
            let valid = cycle.iter().all(|reference| {
                let SharedTrimGeometry::Circle(circle) = self.graph.trims[reference.trim].geometry
                else {
                    return false;
                };
                let tol = self.options.intersection_tolerance;
                (circle.center - expected.center).magnitude() <= tol
                    && (circle.radius - expected.radius).abs() <= tol
                    && circle.normal.cross(expected.normal).magnitude() * expected.radius <= tol
            });
            if valid {
                return Some(Ring {
                    axis,
                    fixed,
                    positive: sweep > 0.0,
                });
            }
        }
        None
    }

    const fn coordinate_tolerance(&self, axis: usize) -> f64 {
        if axis == 1 && matches!(self.carrier, SurfaceCarrier::LinearRevolution { .. }) {
            self.options.intersection_tolerance
        } else {
            self.options.angular_resolution
        }
    }

    fn covers_band(&self, other: usize, lower: f64, span: f64) -> bool {
        self.supports.iter().all(|&id| {
            self.evidence.planar_supports[id]
                .source_face_indices
                .iter()
                .all(|&polygon| {
                    let Some(ids) = self.evidence.source_polygon_vertices(polygon) else {
                        return false;
                    };
                    ids.iter().all(|&id| {
                        let point = self.positions[id];
                        let value =
                            chart_point(self.carrier, point, self.options.angular_resolution)
                                .map(|uv| uv[other])
                                .or_else(|| {
                                    if other != 1 {
                                        return None;
                                    }
                                    match self.carrier {
                                        SurfaceCarrier::CircularRevolution {
                                            axial_center,
                                            axis,
                                            radial_center,
                                            ..
                                        } if radial_center.abs() <= crate::KERNEL_EPSILON => Some(
                                            (point - axial_center).dot(axis).signum() * PI * 0.5,
                                        ),
                                        SurfaceCarrier::LinearRevolution {
                                            origin,
                                            axis,
                                            radius_at_origin,
                                            slope,
                                        } if slope != 0.0 => {
                                            let apex = -radius_at_origin / slope;
                                            let tip = origin + axis * apex;
                                            ((point - tip).magnitude()
                                                <= self.options.max_vertex_displacement)
                                                .then_some(apex)
                                        }
                                        _ => None,
                                    }
                                });
                        let Some(value) = value else {
                            return false;
                        };
                        let mut coordinate = value - lower;
                        if periodic_coordinate(self.carrier, other) {
                            coordinate = coordinate.rem_euclid(TAU);
                            if TAU - coordinate <= self.options.angular_resolution {
                                coordinate = 0.0;
                            }
                        }
                        coordinate >= -self.coordinate_tolerance(other)
                            && coordinate <= span + self.coordinate_tolerance(other)
                    })
                })
        })
    }

    fn seam(
        &self,
        axis: usize,
        angle: f64,
        start: f64,
        sweep: f64,
        front: &Vertex,
        back: &Vertex,
    ) -> Option<Vec<Edge>> {
        let other = 1 - axis;
        let edges = if let Some(circle) = parameter_circle(self.carrier, other, angle) {
            // Isoparametric curve can itself be periodic. Split long
            // connector arcs so each native edge has distinct endpoints.
            let count = if sweep.abs() > PI { 2 } else { 1 };
            let mut edges = Vec::new();
            let mut vertex = front.clone();
            for i in 0..count {
                let next = if i + 1 == count {
                    back.clone()
                } else {
                    let mut uv = [0.0; 2];
                    uv[axis] = angle;
                    uv[other] = sweep.mul_add(0.5, start);
                    builder::vertex(world_point(self.carrier, uv)?)
                };
                edges.push(
                    native_circle_edge(
                        &vertex,
                        &next,
                        circle.center,
                        circle.normal,
                        circle.radius,
                        sweep / f64::from(count),
                        self.options.intersection_tolerance,
                    )
                    .ok()?,
                );
                vertex = next;
            }
            edges
        } else {
            vec![builder::line(front, back)]
        };
        edges
            .iter()
            .all(|e| e.is_geometric_consistent() && self.surface.include(&e.curve()))
            .then_some(edges)
    }
}

fn coordinate_difference(carrier: SurfaceCarrier, axis: usize, a: f64, b: f64) -> f64 {
    if periodic_coordinate(carrier, axis) {
        (a - b + PI).rem_euclid(TAU) - PI
    } else {
        a - b
    }
}

fn world_point(carrier: SurfaceCarrier, [theta, profile]: [f64; 2]) -> Option<Point3> {
    match carrier {
        SurfaceCarrier::Plane { .. } => None,
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } => {
            let (u, v) = normal_basis(axis);
            let (sin, cos) = theta.sin_cos();
            Some(
                origin
                    + axis * profile
                    + (u * cos + v * sin) * slope.mul_add(profile, radius_at_origin),
            )
        }
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } => {
            let (u, v) = normal_basis(axis);
            let (sin, cos) = theta.sin_cos();
            let (z, rho) = profile.sin_cos();
            Some(
                axial_center
                    + axis * (profile_radius * z)
                    + (u * cos + v * sin) * profile_radius.mul_add(rho, radial_center),
            )
        }
    }
}

fn parameter_circle(carrier: SurfaceCarrier, varying: usize, fixed: f64) -> Option<Circle> {
    match carrier {
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } if varying == 0 => Some(Circle {
            center: origin + axis * fixed,
            normal: axis,
            radius: slope.mul_add(fixed, radius_at_origin),
        }),
        SurfaceCarrier::Plane { .. } | SurfaceCarrier::LinearRevolution { .. } => None,
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } => {
            if varying == 0 {
                let (sin, cos) = fixed.sin_cos();
                Some(Circle {
                    center: axial_center + axis * (profile_radius * sin),
                    normal: axis,
                    radius: profile_radius.mul_add(cos, radial_center),
                })
            } else {
                let (u, v) = normal_basis(axis);
                let (sin, cos) = fixed.sin_cos();
                let radial = u * cos + v * sin;
                Some(Circle {
                    center: axial_center + radial * radial_center,
                    normal: radial.cross(axis),
                    radius: profile_radius,
                })
            }
        }
    }
    .filter(|circle| circle.radius > 0.0)
}

fn copy_rings<'a>(
    references: impl Iterator<Item = &'a OrientedTrim>,
    trims: &[Result<Vec<Edge>, NativeTrimError>],
) -> Option<BTreeMap<usize, Vec<Edge>>> {
    let mut changed = BTreeMap::new();
    for reference in references {
        if let std::collections::btree_map::Entry::Vacant(entry) = changed.entry(reference.trim) {
            entry.insert(trims.get(reference.trim)?.as_ref().ok()?.clone());
        }
    }
    Some(changed)
}

fn commit_rings(
    changed: BTreeMap<usize, Vec<Edge>>,
    trims: &mut [Result<Vec<Edge>, NativeTrimError>],
) {
    for (id, edges) in changed {
        trims[id] = Ok(edges);
    }
}

fn cut_ring(
    cycle: &[OrientedTrim],
    point: Point3,
    trims: &mut BTreeMap<usize, Vec<Edge>>,
    tolerance: f64,
) -> Option<Vertex> {
    for reference in cycle {
        let edges = trims.get(&reference.trim)?;
        for edge in edges {
            for vertex in [edge.front(), edge.back()] {
                if (vertex.point() - point).magnitude() <= tolerance {
                    return Some(vertex.clone());
                }
            }
        }
    }
    let vertex = builder::vertex(point);
    let mut cuts = Vec::new();
    for reference in cycle {
        let edges = trims.get(&reference.trim)?;
        for (index, edge) in edges.iter().enumerate() {
            if let Some((left, right)) = edge.cut(&vertex)
                && left.is_geometric_consistent()
                && right.is_geometric_consistent()
            {
                cuts.push((reference.trim, index, left, right));
            }
        }
    }
    let [(id, index, left, right)] = cuts.as_slice() else {
        return None;
    };
    trims
        .get_mut(id)?
        .splice(*index..=*index, [left.clone(), right.clone()]);
    Some(vertex)
}

fn split_wire(wire: &Wire, first: &Vertex, second: &Vertex) -> Option<[Vec<Edge>; 2]> {
    let edges = wire.edge_iter().cloned().collect::<Vec<_>>();
    let start = edges.iter().position(|edge| edge.front() == first)?;
    let end = edges.iter().position(|edge| edge.front() == second)?;
    if start == end {
        return None;
    }
    let distance = (end + edges.len() - start) % edges.len();
    Some([
        (0..distance)
            .map(|i| edges[(start + i) % edges.len()].clone())
            .collect(),
        (distance..edges.len())
            .map(|i| edges[(start + i) % edges.len()].clone())
            .collect(),
    ])
}
