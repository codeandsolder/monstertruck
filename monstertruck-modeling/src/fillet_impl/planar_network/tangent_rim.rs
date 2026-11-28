//! Orthogonal tangent cap-rim rolling ball: exact straight cylinders and
//! circular tori, with one shared meridian at each original wire vertex.
//! This is a local analytic construction; full solid interference and source
//! equivalence are separate admission obligations.
use super::{ANGULAR, Incidence, certify_generated_face, line_axis, orient};
use crate::{Curve, Edge, EdgeId, Face, Shell, Surface, Vertex, VertexId, Wire, builder};
use monstertruck_geometry::{
    analytic_carrier::AnalyticSupport,
    analytic_curve_incidence::{CircularArcFrame, CurveIncidence, certify_circular_arc_coverage},
    prelude::*,
};
use rustc_hash::{FxHashMap as Map, FxHashSet as Set};

#[derive(Clone, Copy)]
enum Side {
    Planar {
        outward: Vector3,
    },
    Circular {
        center: Point3,
        radius: f64,
        sense: f64,
    },
}

impl Side {
    fn outward(self, point: Point3) -> Vector3 {
        match self {
            Self::Planar { outward } => outward,
            Self::Circular {
                center,
                radius,
                sense,
            } => (point - center) * (sense / radius),
        }
    }

    fn cap_curve(self, edge: &Edge, radius: f64) -> Curve {
        let transform = match self {
            Self::Planar { outward } => Matrix4::new_translation(&(-outward * radius)),
            Self::Circular {
                center,
                radius: original,
                sense,
            } => {
                Matrix4::new_translation(&center.position_vector())
                    * Matrix4::new_scaling((sense.mul_add(-radius, original)) / original)
                    * Matrix4::new_translation(&-center.position_vector())
            }
        };
        edge.oriented_curve().transformed(transform)
    }

    fn surface(self, cap: &Edge, normal: Vector3, offset: f64) -> Option<Surface> {
        match self {
            Self::Planar { .. } => {
                let start = cap.front().point();
                let end = cap.back().point();
                let axis = (end - start).normalize();
                let mut surface =
                    Surface::Revolution(Processor::new(RevolutionSurface::by_revolution(
                        Curve::Line(Line(start, end)),
                        start - normal * offset,
                        axis,
                    )));
                if surface.normal(0.5, 0.0).dot(&normal) < 0.0 {
                    surface.invert();
                }
                Some(surface)
            }
            Self::Circular {
                center,
                radius: original,
                sense,
            } => {
                let radius = offset.abs();
                let major = sense.mul_add(-offset, original);
                // A ring torus has a regular, globally oriented support.
                if major <= radius + KERNEL_EPSILON {
                    return None;
                }
                let x = (cap.front().point() - center).normalize();
                let y = normal.cross(&x).normalize();
                let transform = Matrix4::from_columns(&[
                    x.push(0.0),
                    y.push(0.0),
                    normal.push(0.0),
                    (center - normal * offset).to_homogeneous(),
                ]);
                let torus = Torus::new(Point3::origin(), major, radius).ok()?;
                let mut surface = Surface::Torus(Processor::with_transform(torus, transform));
                if offset < 0.0 {
                    surface.invert();
                }
                Some(surface)
            }
        }
    }
}

struct RimEdge {
    edge: Edge,
    side_face: usize,
    side: Side,
    convex: bool,
}

/// Determine whether the complete directed cap edge bounds a convex rim or
/// a reentrant shoulder, using the certified line/circle carrier tangent.
fn rim_convexity(edge: &Edge, side: Side, normal: Vector3) -> Option<bool> {
    let curve = edge.oriented_curve();
    let (start, end) = curve.range_tuple();
    // The directed cap wire's outside conormal equals the side's
    // physical outward normal at convex rims and its negative at
    // reentrant shoulders. Both must agree over the complete edge.
    let initial = curve.derivative(start);
    let magnitude = initial.norm();
    if !magnitude.is_finite() || magnitude <= KERNEL_EPSILON {
        return None;
    }
    let convex = initial
        .normalize()
        .cross(&normal)
        .dot(&side.outward(edge.front().point()))
        > 0.0;
    let sense = if convex { 1.0 } else { -1.0 };
    for (t, vertex) in [(start, edge.front()), (end, edge.back())] {
        let tangent = curve.derivative(t);
        let magnitude = tangent.norm();
        let residual =
            (tangent.normalize().cross(&normal) - side.outward(vertex.point()) * sense).norm();
        if !magnitude.is_finite()
            || magnitude <= KERNEL_EPSILON
            || !residual.is_finite()
            || residual > ANGULAR
        {
            return None;
        }
    }
    Some(convex)
}

impl RimEdge {
    fn new(
        shell: &Shell,
        incidence: &Incidence,
        edge: Edge,
        cap: usize,
        normal: Vector3,
        radius: f64,
    ) -> Option<Self> {
        let source = incidence.edges.get(&edge.id())?;
        let [first, second] = source.uses.as_slice() else {
            return None;
        };
        if first.forward == second.forward || first.face == second.face {
            return None;
        }
        let side_face = if first.face == cap {
            second.face
        } else if second.face == cap {
            first.face
        } else {
            return None;
        };
        let surface = shell.get(side_face)?.oriented_surface();
        let carrier = crate::geometry::analytic_carrier(&surface)?;
        if crate::certify_analytic_curve_on_surface(&edge.curve(), &surface, KERNEL_EPSILON)
            != CurveIncidence::Proven
            || crate::certify_analytic_curve_on_surface(
                &edge.curve(),
                &shell[cap].surface(),
                KERNEL_EPSILON,
            ) != CurveIncidence::Proven
        {
            return None;
        }
        let sign = if carrier.positive_orientation() {
            1.0
        } else {
            -1.0
        };
        let point = edge.front().point();
        let side = match carrier.support() {
            AnalyticSupport::Plane {
                normal: outward, ..
            } => {
                line_axis(&edge)?;
                let outward = outward * sign;
                if normal.dot(&outward).abs() > ANGULAR {
                    return None;
                }
                Side::Planar { outward }
            }
            AnalyticSupport::Cylinder {
                origin,
                axis,
                radius: original,
            } => {
                if normal.cross(&axis).norm() > ANGULAR {
                    return None;
                }
                let center = origin + axis * (point - origin).dot(&axis);
                let frame = CircularArcFrame::new(center, normal, original, point, KERNEL_EPSILON)?;
                let rational =
                    NurbsCurve::new(edge.oriented_curve().try_into_homogeneous_bspline_curve()?);
                certify_circular_arc_coverage(&rational, frame)?;
                Side::Circular {
                    center,
                    radius: original,
                    sense: sign,
                }
            }
            _ => return None,
        };
        let convex = rim_convexity(&edge, side, normal)?;
        let sense = if convex { 1.0 } else { -1.0 };
        if let Side::Circular {
            radius: original,
            sense: radial,
            ..
        } = side
            && (radial * sense).mul_add(-radius, original) <= radius + KERNEL_EPSILON
        {
            return None;
        }
        Some(Self {
            edge,
            side_face,
            side,
            convex,
        })
    }
}

struct Junction {
    cap: Vertex,
    side: Vertex,
    meridian: Edge,
}

fn junctions(rim: &[RimEdge], normal: Vector3, radius: f64) -> Option<Vec<Junction>> {
    let mut seen = Set::default();
    rim.iter()
        .enumerate()
        .map(|(i, edge)| {
            let previous = &rim[(i + rim.len() - 1) % rim.len()];
            let vertex = edge.edge.front();
            if !seen.insert(vertex.id()) || previous.edge.back() != vertex {
                return None;
            }
            let point = vertex.point();
            let outward = edge.side.outward(point);
            if (previous.side.outward(point) - outward).norm() > ANGULAR {
                return None;
            }
            let cap = Vertex::new(point - outward * radius);
            let side = Vertex::new(point - normal * radius);
            let center = point - (outward + normal) * radius;
            let middle = center + (outward + normal) * (radius * std::f64::consts::FRAC_1_SQRT_2);
            let meridian = builder::try_circle_arc(&cap, &side, middle).ok()?;
            Some(Junction {
                cap,
                side,
                meridian,
            })
        })
        .collect()
}

struct Patch {
    cap: Edge,
    side: Edge,
    face: Face,
}

impl Patch {
    fn new(
        rim: &RimEdge,
        start: &Junction,
        end: &Junction,
        normal: Vector3,
        radius: f64,
    ) -> Option<Self> {
        let cap = Edge::new(&start.cap, &end.cap, rim.side.cap_curve(&rim.edge, radius)).ok()?;
        let side = Edge::new(
            &start.side,
            &end.side,
            rim.edge
                .oriented_curve()
                .transformed(Matrix4::new_translation(&(-normal * radius))),
        )
        .ok()?;
        let surface = rim.side.surface(&cap, normal, radius)?;
        let wire = Wire::from(vec![
            cap.inverse(),
            start.meridian.clone(),
            side.clone(),
            end.meridian.inverse(),
        ]);
        let face = Face::new(vec![wire], surface).ok()?;
        certify_generated_face(&face).then_some(Self { cap, side, face })
    }
}

/// The old side generators must each be a complete straight segment whose
/// new terminal vertex lies strictly inside its original finite extent.
fn trim_generators(
    incidence: &Incidence,
    rim: &[RimEdge],
    junctions: &[Junction],
    selected: &Set<EdgeId>,
    normal: Vector3,
    offset: f64,
) -> Option<Map<EdgeId, Edge>> {
    let mut replacements = Map::default();
    let vertices = rim
        .iter()
        .zip(junctions)
        .map(|(edge, junction)| (edge.edge.front().id(), &junction.side))
        .collect::<Map<VertexId, _>>();
    for (i, edge) in rim.iter().enumerate() {
        let vertex = edge.edge.front();
        let prev = &rim[(i + rim.len() - 1) % rim.len()];
        for id in &incidence.at_vertex[&vertex.id()] {
            if selected.contains(id) {
                continue;
            }
            let original = &incidence.edges[id];
            if replacements.contains_key(id)
                || original.uses.len() != 2
                || original
                    .uses
                    .iter()
                    .any(|u| u.face != edge.side_face && u.face != prev.side_face)
            {
                return None;
            }
            let source = &original.edge;
            line_axis(source)?;
            let (other, at_front) = if source.absolute_front() == vertex {
                (source.absolute_back(), true)
            } else {
                (source.absolute_front(), false)
            };
            if vertices.contains_key(&other.id()) {
                return None;
            }
            let delta = other.point() - vertex.point();
            if delta.dot(&-normal) * offset.signum() <= offset.abs() + KERNEL_EPSILON
                || delta.cross(&normal).norm() > KERNEL_EPSILON
            {
                return None;
            }
            let replacement = if at_front {
                builder::line(vertices[&vertex.id()], other)
            } else {
                builder::line(other, vertices[&vertex.id()])
            };
            replacements.insert(*id, replacement);
        }
    }
    Some(replacements)
}

fn replace_faces(
    shell: &Shell,
    cap: usize,
    rim: &[RimEdge],
    patches: &[Patch],
    retained: &Map<EdgeId, Edge>,
) -> Option<Shell> {
    let by_id = rim
        .iter()
        .zip(patches)
        .map(|(rim, patch)| (rim.edge.id(), (rim, patch)))
        .collect::<Map<_, _>>();
    let mut result = Vec::with_capacity(shell.len() + patches.len());
    for (index, face) in shell.iter().enumerate() {
        let mut changed = false;
        let mut boundaries = Vec::new();
        for wire in face.boundaries() {
            let edges = wire
                .iter()
                .map(|edge| {
                    let replacement = if let Some((rim, patch)) = by_id.get(&edge.id()) {
                        let rail = if index == cap {
                            &patch.cap
                        } else if index == rim.side_face {
                            &patch.side
                        } else {
                            return None;
                        };
                        orient(rail, edge.orientation() == rim.edge.orientation())
                    } else if let Some(replacement) = retained.get(&edge.id()) {
                        orient(replacement, edge.orientation())
                    } else {
                        return Some(edge.clone());
                    };
                    changed = true;
                    Some(replacement)
                })
                .collect::<Option<Vec<_>>>()?;
            boundaries.push(Wire::from(edges));
        }
        if changed {
            let mut updated = Face::new(boundaries, face.oriented_surface()).ok()?;
            updated.set_stable_id(face.stable_id());
            if !updated.is_geometric_consistent() {
                return None;
            }
            result.push(updated);
        } else {
            result.push(face.clone());
        }
    }
    result.extend(patches.iter().map(|patch| patch.face.clone()));
    let result = Shell::from(result);
    result.check_solid_boundary().ok()?;
    Some(result)
}

pub(super) fn try_round_tangent_rim(shell: &Shell, edges: &[Edge], radii: &[f64]) -> Option<Shell> {
    let radius = *radii.first()?;
    if !radius.is_finite()
        || radius <= KERNEL_EPSILON
        || edges.len() != radii.len()
        || radii
            .iter()
            .any(|r| !r.is_finite() || (r - radius).abs() > KERNEL_EPSILON)
    {
        return None;
    }
    let selected = edges.iter().map(Edge::id).collect::<Set<_>>();
    if selected.len() != edges.len() {
        return None;
    }
    let mut candidates = shell.iter().enumerate().filter_map(|(index, face)| {
        let Surface::Plane(plane) = face.oriented_surface() else {
            return None;
        };
        let mut wires = face.boundaries().into_iter().filter(|wire| {
            wire.len() == selected.len() && wire.iter().all(|edge| selected.contains(&edge.id()))
        });
        let wire = wires.next()?;
        (wires.next().is_none() && wire.is_closed()).then_some((index, plane.normal(), wire))
    });
    let (cap, normal, wire) = candidates.next()?;
    if candidates.next().is_some() {
        return None;
    }
    let incidence = Incidence::new(shell);
    let rim = wire
        .iter()
        .map(|edge| RimEdge::new(shell, &incidence, edge.clone(), cap, normal, radius))
        .collect::<Option<Vec<_>>>()?;
    let convex = rim.first()?.convex;
    if rim.iter().any(|edge| edge.convex != convex) {
        return None;
    }
    let offset = if convex { radius } else { -radius };
    let junctions = junctions(&rim, normal, offset)?;
    let retained = trim_generators(&incidence, &rim, &junctions, &selected, normal, offset)?;
    let patches = rim
        .iter()
        .enumerate()
        .map(|(i, edge)| {
            Patch::new(
                edge,
                &junctions[i],
                &junctions[(i + 1) % junctions.len()],
                normal,
                offset,
            )
        })
        .collect::<Option<Vec<_>>>()?;
    replace_faces(shell, cap, &rim, &patches, &retained)
}
