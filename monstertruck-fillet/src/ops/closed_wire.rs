//! Closed-wire reconstruction when neighboring fillet sections meet at the
//! same sharp junction. Vertices and seams are owned by the ORIGINAL wire
//! incidence; no geometric weld of unrelated shell topology is performed.

use super::{Edge, Face, FaceBoundaryEdgeIndex, Shell, Vertex, Wire};
use crate::types::{EdgeId, VertexId};
use monstertruck_geometry::prelude::*;
use monstertruck_topology::shell::ShellCondition;
use std::collections::{HashMap, HashSet};

struct Junction {
    cap: Vertex,
    side: Vertex,
    seam: Edge,
}

struct RailPair {
    source_forward: bool,
    side_face: usize,
    cap: Edge,
    side: Edge,
}

fn make_junctions(wire: &Wire, surfaces: &[NurbsSurface<Vector4>]) -> Option<Vec<Junction>> {
    let count = surfaces.len();
    let mut corners = HashSet::new();
    let mut result = Vec::with_capacity(count);
    for (i, edge) in wire.edge_iter().enumerate() {
        if !corners.insert(edge.front().id()) {
            return None;
        }
        let previous = &surfaces[(i + count - 1) % count];
        let surface = &surfaces[i];
        if !previous.is_clamped() || !surface.is_clamped() {
            return None;
        }
        let seam = surface.curve_u(0);
        let preceding = previous.curve_u(previous.control_points()[0].len() - 1);
        // The reconciled sections must have the same knot frame and
        // homogeneous controls, not merely spatially close terminal points.
        if seam.knot_vector() != preceding.knot_vector()
            || seam.control_points().len() != preceding.control_points().len()
            || !seam
                .control_points()
                .iter()
                .zip(preceding.control_points())
                .all(|(a, b)| a.near(b))
        {
            return None;
        }
        let cap = Vertex::new(seam.front());
        let side = Vertex::new(seam.back());
        if cap.point().near(&side.point()) {
            return None;
        }
        let seam = Edge::new(&cap, &side, seam.into()).ok()?;
        result.push(Junction { cap, side, seam });
    }
    Some(result)
}

fn make_rails(
    wire: &Wire,
    adjacent: &[FaceBoundaryEdgeIndex],
    surfaces: &[NurbsSurface<Vector4>],
    junctions: &[Junction],
) -> Option<HashMap<EdgeId, RailPair>> {
    wire.edge_iter()
        .enumerate()
        .map(|(i, original)| {
            let a = &junctions[i];
            let b = &junctions[(i + 1) % junctions.len()];
            let surface = &surfaces[i];
            let cap = Edge::new(&a.cap, &b.cap, surface.curve_v(0).into()).ok()?;
            let side = Edge::new(
                &a.side,
                &b.side,
                surface.curve_v(surface.control_points().len() - 1).into(),
            )
            .ok()?;
            if !cap.is_geometric_consistent() || !side.is_geometric_consistent() {
                return None;
            }
            Some((
                original.id(),
                RailPair {
                    source_forward: original.orientation(),
                    side_face: adjacent[i].face_index,
                    cap,
                    side,
                },
            ))
        })
        .collect()
}

/// Trim a surviving original edge to the new side-rail vertices. The original
/// curve remains authoritative: moving an endpoint off its support is not an
/// allowed topological stitch.
fn trim_retained(edge: &Edge, corners: &HashMap<VertexId, &Vertex>) -> Option<Edge> {
    let front = edge.absolute_front();
    let back = edge.absolute_back();
    let new_front = corners.get(&front.id()).copied().unwrap_or(front);
    let new_back = corners.get(&back.id()).copied().unwrap_or(back);
    let mut curve = edge.curve();
    let (start, end) = curve.range_tuple();
    let parameter = |vertex: &Vertex| {
        if vertex.point().near(&curve.front()) {
            return Some(start);
        }
        if vertex.point().near(&curve.back()) {
            return Some(end);
        }
        let t = curve
            .search_parameter(vertex.point(), None, 100)
            .or_else(|| curve.search_nearest_parameter(vertex.point(), None, 100))?;
        (t.is_finite() && t >= start && t <= end && curve.evaluate(t).near(&vertex.point()))
            .then_some(t)
    };
    let a = parameter(new_front)?;
    let b = parameter(new_back)?;
    if a >= b || new_front.point().near(&new_back.point()) {
        return None;
    }
    if a > start {
        curve = curve.cut(a);
    }
    if b < end {
        let _ = curve.cut(b);
    }
    let replacement = Edge::new(new_front, new_back, curve).ok()?;
    replacement.is_geometric_consistent().then_some(replacement)
}

fn shortened_edges(
    shell: &Shell,
    wire: &Wire,
    junctions: &[Junction],
    rails: &HashMap<EdgeId, RailPair>,
) -> Option<HashMap<EdgeId, Edge>> {
    let corners = wire
        .edge_iter()
        .zip(junctions)
        .map(|(edge, j)| (edge.front().id(), &j.side))
        .collect::<HashMap<_, _>>();
    let mut result = HashMap::new();
    for edge in shell.edge_iter() {
        if rails.contains_key(&edge.id()) || result.contains_key(&edge.id()) {
            continue;
        }
        if corners.contains_key(&edge.absolute_front().id())
            || corners.contains_key(&edge.absolute_back().id())
        {
            result.insert(edge.id(), trim_retained(&edge, &corners)?);
        }
    }
    Some(result)
}

fn rewrite_face(
    face: &Face,
    index: usize,
    cap_face: usize,
    rails: &HashMap<EdgeId, RailPair>,
    shortened: &HashMap<EdgeId, Edge>,
) -> Option<Face> {
    let mut changed = false;
    let mut boundaries = Vec::new();
    for mut wire in face.boundaries() {
        for edge in wire.iter_mut() {
            let replacement = if let Some(pair) = rails.get(&edge.id()) {
                let rail = if index == cap_face {
                    &pair.cap
                } else if index == pair.side_face {
                    &pair.side
                } else {
                    return None;
                };
                Some(if edge.orientation() == pair.source_forward {
                    rail.clone()
                } else {
                    rail.inverse()
                })
            } else {
                shortened.get(&edge.id()).map(|trimmed| {
                    if edge.orientation() {
                        trimmed.clone()
                    } else {
                        trimmed.inverse()
                    }
                })
            };
            if let Some(replacement) = replacement {
                changed = true;
                *edge = replacement;
            }
        }
        boundaries.push(wire);
    }
    if changed {
        Face::new(boundaries, face.oriented_surface()).ok()
    } else {
        Some(face.clone())
    }
}

/// Build a detached shell using one topological object for each retained cut,
/// rail, and terminal seam. Applicable to any clamped generated fillet
/// surfaces whose junctions actually coincide with retained source supports,
/// including side faces spanning multiple selected edges. Non-matching
/// junctions need a miter/corner intersection rather than a forced weld.
pub(super) fn reconstruct(
    shell: &Shell,
    wire: &Wire,
    cap_face: usize,
    adjacent: &[FaceBoundaryEdgeIndex],
    surfaces: &[NurbsSurface<Vector4>],
) -> Option<Shell> {
    if wire.len() != surfaces.len() || adjacent.len() != surfaces.len() || !wire.is_closed() {
        return None;
    }
    let junctions = make_junctions(wire, surfaces)?;
    let rails = make_rails(wire, adjacent, surfaces, &junctions)?;
    let shortened = shortened_edges(shell, wire, &junctions, &rails)?;
    let mut result = shell
        .iter()
        .enumerate()
        .map(|(i, face)| rewrite_face(face, i, cap_face, &rails, &shortened))
        .collect::<Option<Shell>>()?;
    for (i, original) in wire.edge_iter().enumerate() {
        let rail = &rails[&original.id()];
        let seam_start = &junctions[i].seam;
        let seam_end = &junctions[(i + 1) % junctions.len()].seam;
        let boundary = [
            rail.cap.inverse(),
            seam_start.clone(),
            rail.side.clone(),
            seam_end.inverse(),
        ]
        .into();
        result.push(Face::new(vec![boundary], surfaces[i].clone()).ok()?);
    }
    if shell.is_connected() && !result.is_connected() {
        return None;
    }
    if shell.shell_condition() == ShellCondition::Closed
        && result.shell_condition() != ShellCondition::Closed
    {
        return None;
    }
    Some(result)
}
