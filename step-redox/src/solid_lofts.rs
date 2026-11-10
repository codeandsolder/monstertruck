//! Exact recognition of planar, sectioned polyhedra as editable ruled lofts.
//!
//! A section is inferred from actual source `EDGE_CURVE` vertices. No faceting
//! or interpolation is used: every claimed loft wall is a matching source
//! planar quadrilateral with identical topology and geometric vertices.

use crate::brep::{CurveSupport, OrientedEdgeUse, SurfaceSupport};
use crate::cad_ir::{CadModel, CadNode, LoftSection, Profile2d, RigidTransform};
use crate::cad_recovery::{CadFragment, CadFragmentSource};
use crate::math3::{cross, distance, dot, normalize, sub};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use ruststep::ast::EntityInstance;
use std::collections::{BTreeMap, BTreeSet, HashMap};

const TOL_MM: f64 = 1.0e-7;
const DIR_TOL: f64 = 1.0e-10;

#[derive(Clone, Copy)]
struct Vertex {
    id: u64,
    point: [f64; 3],
}

struct Section {
    level: f64,
    corners: [Vertex; 4],
    profile: Profile2d,
}

fn plane_quad(face: &ReducedFace) -> bool {
    matches!(face.surface, SurfaceSupport::Plane(_))
        && face.loops.len() == 1
        && face.loops[0].outer
        && face.loops[0].edges.len() == 4
        && face.loops[0]
            .edges
            .iter()
            .all(|e| matches!(e.support, CurveSupport::Line(_)))
}

fn points(solid: &ReducedSolid) -> Option<BTreeMap<u64, [f64; 3]>> {
    let mut points = BTreeMap::new();
    for face in solid.faces() {
        if !plane_quad(face) {
            return None;
        }
        let SurfaceSupport::Plane(plane) = &face.surface else {
            return None;
        };
        for edge in face.loop_edges() {
            for (vertex, point) in [
                (edge.start_vertex, edge.start_mm),
                (edge.end_vertex, edge.end_mm),
            ] {
                if !point.iter().all(|v| v.is_finite())
                    || dot(sub(point, plane.origin_mm), plane.normal).abs() > TOL_MM
                {
                    return None;
                }
                if points
                    .insert(vertex, point)
                    .is_some_and(|prior| distance(prior, point) > TOL_MM)
                {
                    return None;
                }
            }
        }
    }
    Some(points)
}

fn section(
    level: f64,
    vertices: &[Vertex],
    origin: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
) -> Option<Section> {
    if vertices.len() != 4 {
        return None;
    }
    let coords = vertices
        .iter()
        .map(|v| (dot(sub(v.point, origin), x), dot(sub(v.point, origin), y)))
        .collect::<Vec<_>>();
    let xmin = coords.iter().map(|v| v.0).min_by(f64::total_cmp)?;
    let xmax = coords.iter().map(|v| v.0).max_by(f64::total_cmp)?;
    let ymin = coords.iter().map(|v| v.1).min_by(f64::total_cmp)?;
    let ymax = coords.iter().map(|v| v.1).max_by(f64::total_cmp)?;
    if xmax - xmin <= TOL_MM || ymax - ymin <= TOL_MM {
        return None;
    }
    let corners_2d = [(xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax)];
    let mut matched = [None; 4];
    for (v, (cx, cy)) in vertices.iter().zip(coords) {
        let index = corners_2d
            .iter()
            .position(|(x, y)| (x - cx).abs() <= TOL_MM && (y - cy).abs() <= TOL_MM)?;
        if matched[index].replace(*v).is_some() {
            return None;
        }
    }
    let corners = matched
        .into_iter()
        .collect::<Option<Vec<_>>>()?
        .try_into()
        .ok()?;
    let profile = Profile2d::polygon(corners_2d.map(Into::into).to_vec()).ok()?;
    Some(Section {
        level,
        corners,
        profile,
    })
}

const fn edge_pair(a: u64, b: u64) -> (u64, u64) { if a < b { (a, b) } else { (b, a) } }

fn source_face_has_edges(face: &ReducedFace, expected: &[(u64, u64); 4]) -> bool {
    face.loops[0]
        .edges
        .iter()
        .map(|edge| edge_pair(edge.start_vertex, edge.end_vertex))
        .collect::<BTreeSet<_>>()
        == expected.iter().copied().collect::<BTreeSet<_>>()
}

fn face_partition_proven(
    solid: &ReducedSolid,
    sections: &[Section],
    axis: [f64; 3],
) -> Option<[u64; 2]> {
    let mut vertex_section = BTreeMap::<u64, (usize, usize)>::new();
    for (stage, section) in sections.iter().enumerate() {
        for (corner, vertex) in section.corners.iter().enumerate() {
            vertex_section.insert(vertex.id, (stage, corner));
        }
    }
    let mut caps = [None; 2];
    let mut sides = BTreeSet::new();
    for face in solid.faces() {
        let verts = face.loops[0]
            .edges
            .iter()
            .map(|edge| edge.start_vertex)
            .collect::<Vec<_>>();
        if verts.iter().copied().collect::<BTreeSet<_>>().len() != 4 {
            return None;
        }
        let mut indices = verts
            .iter()
            .map(|v| vertex_section.get(v).copied())
            .collect::<Option<Vec<_>>>()?;
        indices.sort_unstable();
        let (stage0, _) = indices[0];
        if indices.iter().all(|(stage, _)| *stage == stage0) {
            let which = if stage0 == 0 {
                0
            } else if stage0 == sections.len() - 1 {
                1
            } else {
                return None;
            };
            let SurfaceSupport::Plane(plane) = face.surface else {
                return None;
            };
            if dot(normalize(plane.normal, DIR_TOL)?, axis).abs() < 1.0 - DIR_TOL {
                return None;
            }
            let c = &sections[stage0].corners;
            let edges = std::array::from_fn(|i| edge_pair(c[i].id, c[(i + 1) % 4].id));
            if !source_face_has_edges(face, &edges) || caps[which].replace(face.id).is_some() {
                return None;
            }
        } else {
            let next = stage0 + 1;
            if next >= sections.len()
                || indices.iter().filter(|(s, _)| *s == stage0).count() != 2
                || indices.iter().filter(|(s, _)| *s == next).count() != 2
                || indices.iter().any(|(s, _)| *s != stage0 && *s != next)
            {
                return None;
            }
            let bottom = indices
                .iter()
                .filter_map(|(s, i)| (*s == stage0).then_some(*i))
                .collect::<Vec<_>>();
            let top = indices
                .iter()
                .filter_map(|(s, i)| (*s == next).then_some(*i))
                .collect::<Vec<_>>();
            let side = (0..4).find(|&i| {
                let target = [i, (i + 1) % 4];
                bottom.iter().all(|v| target.contains(v)) && top.iter().all(|v| target.contains(v))
            })?;
            if !sides.insert((stage0, side)) {
                return None;
            }
            let a = &sections[stage0].corners;
            let b = &sections[next].corners;
            let other = (side + 1) % 4;
            let expected = [
                edge_pair(a[side].id, a[other].id),
                edge_pair(a[other].id, b[other].id),
                edge_pair(b[other].id, b[side].id),
                edge_pair(b[side].id, a[side].id),
            ];
            if !source_face_has_edges(face, &expected) {
                return None;
            }
        }
    }
    if sides.len() != 4 * (sections.len() - 1) {
        return None;
    }
    Some([caps[0]?, caps[1]?])
}

fn candidate(
    solid: &ReducedSolid,
    vertices: &BTreeMap<u64, [f64; 3]>,
    cap: &ReducedFace,
) -> Option<(f64, CadFragment)> {
    let SurfaceSupport::Plane(plane) = cap.surface else {
        return None;
    };
    let axis = normalize(plane.normal, DIR_TOL)?;
    let origin = cap.loops[0].edges[0].start_mm;
    let x = normalize(sub(cap.loops[0].edges[0].end_mm, origin), DIR_TOL)?;
    if dot(x, axis).abs() > DIR_TOL {
        return None;
    }
    let y = normalize(cross(axis, x), DIR_TOL)?;
    let mut ranked = vertices
        .iter()
        .map(|(&id, &point)| (dot(sub(point, origin), axis), Vertex { id, point }))
        .collect::<Vec<_>>();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut groups = Vec::<(f64, Vec<Vertex>)>::new();
    for (level, vertex) in ranked {
        if let Some((at, vs)) = groups.last_mut()
            && (level - *at).abs() <= TOL_MM
        {
            vs.push(vertex);
        } else {
            groups.push((level, vec![vertex]));
        }
    }
    if groups.len() < 3 || groups.iter().any(|(_, group)| group.len() != 4) {
        return None;
    }
    let sections = groups
        .into_iter()
        .map(|(level, group)| section(level, &group, origin, x, y))
        .collect::<Option<Vec<_>>>()?;
    let cap_ids = face_partition_proven(solid, &sections, axis)?;
    let span = sections.last()?.level - sections.first()?.level;
    if span <= TOL_MM {
        return None;
    }
    let transforms = sections
        .iter()
        .map(|section| LoftSection {
            profile: section.profile.clone(),
            transform: RigidTransform {
                matrix: [
                    [
                        x[0],
                        y[0],
                        axis[0],
                        axis[0].mul_add(section.level, origin[0]),
                    ],
                    [
                        x[1],
                        y[1],
                        axis[1],
                        axis[1].mul_add(section.level, origin[1]),
                    ],
                    [
                        x[2],
                        y[2],
                        axis[2],
                        axis[2].mul_add(section.level, origin[2]),
                    ],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        })
        .collect();
    let mut model = CadModel::new();
    let root = model.add_node(CadNode::Loft {
        sections: transforms,
    });
    model.add_root(root).ok()?;
    model.validate().ok()?;
    let source = CadFragmentSource::SolidLoft {
        solid_id: solid.solid_id(),
        cap_face_ids: cap_ids,
        side_face_ids: solid
            .faces()
            .filter_map(|f| (!cap_ids.contains(&f.id)).then_some(f.id))
            .collect(),
    };
    Some((
        span,
        CadFragment {
            source,
            model,
            root,
        },
    ))
}

/// Recognize any number of exact rectangular stages and bevel transitions,
/// without imposing source axis, position, scale, face IDs or section count.
///
/// Candidate faces establish a local frame, but the widest valid loft axis
/// is selected so a long chamfered bar is represented along its long axis.
pub fn detect_sectioned_planar_loft(solid: &ReducedSolid) -> Option<CadFragment> {
    solid.validate_closed_topology()?;
    let vertices = points(solid)?;
    let mut candidates = solid
        .faces()
        .filter_map(|cap| candidate(solid, &vertices, cap))
        .collect::<Vec<_>>();
    candidates.sort_by(|a, b| b.0.total_cmp(&a.0));
    candidates.into_iter().next().map(|(_, fragment)| fragment)
}

/// An exact planar end treatment, whose terminal face and one or more
/// four-sided ruled wall rings are real source topology. The shoulder section is an implied interface,
/// never an invented source face.
pub struct RuledEndExtension {
    pub fragment: CadFragment,
    pub tip_face_id: u64,
    pub axis: [f64; 3],
    /// Section plane at the connection to preceding source material.
    pub shoulder_point_mm: [f64; 3],
    /// A point on the actual terminal source cap plane.
    pub tip_point_mm: [f64; 3],
}

/// Discover terminal four-sided ruled caps. Unlike the whole-solid loft
/// recognizer, this works with a mixed analytic source: curved elbow walls
/// remain untouched while the exact planar tip treatment is recovered.
pub fn detect_ruled_end_extensions(solid: &ReducedSolid) -> Vec<RuledEndExtension> {
    solid
        .faces()
        .filter_map(|cap| terminal_ruled_extension(solid, cap))
        .collect()
}

/// Discover source-witnessed ruled end sections, including multi-stage
/// planar strips and analytic walls continuing into an earlier core.
pub fn detect_extended_ruled_end_extensions(solid: &ReducedSolid) -> Vec<RuledEndExtension> {
    solid
        .faces()
        .filter_map(|cap| {
            terminal_sectioned_ruled_extension(solid, cap)
                .or_else(|| terminal_ruled_extension(solid, cap))
                .or_else(|| terminal_extended_planar_ruled_extension(solid, cap))
        })
        .collect()
}

fn matching_ruled_end_shoulder(
    solid: &ReducedSolid,
    cap: &ReducedFace,
    cap_vertices: &BTreeMap<u64, [f64; 3]>,
) -> Option<(BTreeMap<u64, Vertex>, [u64; 4])> {
    let cap_edges = &cap.loops[0].edges;
    let mut matching_shoulder = BTreeMap::<u64, Vertex>::new();
    let mut walls = Vec::new();
    let mut shoulder_edges = BTreeSet::new();
    for edge in cap_edges {
        let adjacent = solid
            .faces()
            .filter(|f| f.id != cap.id && f.loop_edges().any(|use_| use_.edge_id == edge.edge_id))
            .collect::<Vec<_>>();
        let [wall] = adjacent.as_slice() else {
            return None;
        };
        if !plane_quad(wall) || walls.contains(&wall.id) {
            return None;
        }
        walls.push(wall.id);
        let opposite = wall
            .loop_edges()
            .filter(|other| {
                !cap_vertices.contains_key(&other.start_vertex)
                    && !cap_vertices.contains_key(&other.end_vertex)
            })
            .collect::<Vec<_>>();
        let [shoulder_edge] = opposite.as_slice() else {
            return None;
        };
        shoulder_edges.insert(shoulder_edge.edge_id);
        for end_id in [edge.start_vertex, edge.end_vertex] {
            let connectors = wall
                .loop_edges()
                .filter_map(|other| {
                    if other.start_vertex == end_id && !cap_vertices.contains_key(&other.end_vertex)
                    {
                        Some(Vertex {
                            id: other.end_vertex,
                            point: other.end_mm,
                        })
                    } else if other.end_vertex == end_id
                        && !cap_vertices.contains_key(&other.start_vertex)
                    {
                        Some(Vertex {
                            id: other.start_vertex,
                            point: other.start_mm,
                        })
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            let [shoulder_vertex] = connectors.as_slice() else {
                return None;
            };
            if matching_shoulder
                .insert(end_id, *shoulder_vertex)
                .is_some_and(|prior| {
                    prior.id != shoulder_vertex.id
                        || distance(prior.point, shoulder_vertex.point) > TOL_MM
                })
            {
                return None;
            }
        }
    }
    if walls.len() != 4
        || matching_shoulder.len() != 4
        || matching_shoulder
            .values()
            .map(|v| v.id)
            .collect::<BTreeSet<_>>()
            .len()
            != 4
        || shoulder_edges.len() != 4
    {
        return None;
    }
    // If the shoulder itself is another source cap, this is already a
    // complete simple prism/loft, not a terminal added material feature.
    if solid.faces().any(|f| {
        plane_quad(f)
            && f.loop_edges().map(|e| e.edge_id).collect::<BTreeSet<_>>() == shoulder_edges
    }) {
        return None;
    }

    Some((matching_shoulder, walls.try_into().ok()?))
}

fn terminal_ruled_extension(solid: &ReducedSolid, cap: &ReducedFace) -> Option<RuledEndExtension> {
    if !plane_quad(cap) {
        return None;
    }
    let SurfaceSupport::Plane(plane) = &cap.surface else {
        return None;
    };
    let cap_edges = &cap.loops[0].edges;
    let cap_vertices = cap_edges
        .iter()
        .map(|e| (e.start_vertex, e.start_mm))
        .collect::<BTreeMap<_, _>>();
    if cap_vertices.len() != 4 {
        return None;
    }

    let (matching_shoulder, wall_face_ids) =
        matching_ruled_end_shoulder(solid, cap, &cap_vertices)?;
    let axis = normalize(plane.normal, DIR_TOL)?;
    let origin = cap_edges[0].start_mm;
    let x = normalize(sub(cap_edges[0].end_mm, origin), DIR_TOL)?;
    if dot(axis, x).abs() > DIR_TOL {
        return None;
    }
    let y = normalize(cross(axis, x), DIR_TOL)?;
    let tip_vertices = cap_edges
        .iter()
        .map(|edge| Vertex {
            id: edge.start_vertex,
            point: edge.start_mm,
        })
        .collect::<Vec<_>>();
    let shoulder_vertices = cap_edges
        .iter()
        .map(|edge| matching_shoulder.get(&edge.start_vertex).copied())
        .collect::<Option<Vec<_>>>()?;
    let depth = dot(sub(shoulder_vertices[0].point, origin), axis);
    if depth.abs() <= TOL_MM
        || shoulder_vertices
            .iter()
            .any(|v| (dot(sub(v.point, origin), axis) - depth).abs() > TOL_MM)
    {
        return None;
    }
    let tip = section(0.0, &tip_vertices, origin, x, y)?;
    let shoulder = section(depth, &shoulder_vertices, origin, x, y)?;
    // The STEP tip face's outward normal points *away* from its shoulder.
    // Builder::try_skin_wires and its start/end cap conventions assume the
    // sections progress in the positive local axis direction. Reversing a
    // negative-depth pair here is essential: both solids can be closed while
    // the loft's oriented material is inverted, so a Boolean would then keep
    // the entire shared interface as an external face.
    let ordered = if depth < 0.0 {
        [shoulder, tip]
    } else {
        [tip, shoulder]
    };
    let sections = ordered
        .into_iter()
        .map(|section| LoftSection {
            profile: section.profile,
            transform: RigidTransform {
                matrix: [
                    [
                        x[0],
                        y[0],
                        axis[0],
                        axis[0].mul_add(section.level, origin[0]),
                    ],
                    [
                        x[1],
                        y[1],
                        axis[1],
                        axis[1].mul_add(section.level, origin[1]),
                    ],
                    [
                        x[2],
                        y[2],
                        axis[2],
                        axis[2].mul_add(section.level, origin[2]),
                    ],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        })
        .collect::<Vec<_>>();
    let mut model = CadModel::new();
    let root = model.add_node(CadNode::Loft { sections });
    model.add_root(root).ok()?;
    model.validate().ok()?;
    Some(RuledEndExtension {
        fragment: CadFragment {
            source: CadFragmentSource::RuledEndExtension {
                solid_id: solid.solid_id(),
                tip_face_id: cap.id,
                wall_face_ids: wall_face_ids.to_vec(),
            },
            model,
            root,
        },
        tip_face_id: cap.id,
        axis,
        shoulder_point_mm: crate::math3::add(origin, crate::math3::mul(axis, depth)),
        tip_point_mm: origin,
    })
}

/// Follow a source four-sided planar wall ring, preserving the actual STEP
/// vertex and edge identities instead of extrapolating an idealized shoulder.
/// Previous walls are excluded; a cap or a non-quad continuation ends the
/// strip rather than being interpreted as another section.
fn next_quad_wall_ring(
    solid: &ReducedSolid,
    ring: &[OrientedEdgeUse],
    consumed: &BTreeSet<u64>,
) -> Option<(Vec<OrientedEdgeUse>, Vec<u64>)> {
    let vertices = ring
        .iter()
        .flat_map(|edge| [edge.start_vertex, edge.end_vertex])
        .collect::<BTreeSet<_>>();
    if ring.len() != 4 || vertices.len() != 4 {
        return None;
    }
    let mut next = Vec::with_capacity(4);
    let mut walls = Vec::with_capacity(4);
    for edge in ring {
        let adjacent = solid
            .faces()
            .filter(|face| {
                !consumed.contains(&face.id)
                    && face.loop_edges().any(|used| used.edge_id == edge.edge_id)
            })
            .collect::<Vec<_>>();
        let [wall] = adjacent.as_slice() else {
            return None;
        };
        if !plane_quad(wall) || walls.contains(&wall.id) {
            return None;
        }
        let opposite = wall
            .loop_edges()
            .filter(|used| {
                !vertices.contains(&used.start_vertex) && !vertices.contains(&used.end_vertex)
            })
            .collect::<Vec<_>>();
        let [far_edge] = opposite.as_slice() else {
            return None;
        };
        walls.push(wall.id);
        next.push((*far_edge).clone());
    }
    let mut degrees = BTreeMap::<u64, usize>::new();
    for edge in &next {
        *degrees.entry(edge.start_vertex).or_default() += 1;
        *degrees.entry(edge.end_vertex).or_default() += 1;
    }
    if degrees.len() != 4
        || degrees.values().any(|&degree| degree != 2)
        || next
            .iter()
            .map(|edge| edge.edge_id)
            .collect::<BTreeSet<_>>()
            .len()
            != 4
    {
        return None;
    }
    Some((next, walls))
}

fn section_from_wall_ring(
    ring: &[OrientedEdgeUse],
    origin: [f64; 3],
    axis: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
) -> Option<Section> {
    let mut seen = BTreeMap::<u64, [f64; 3]>::new();
    for edge in ring {
        for (id, point) in [
            (edge.start_vertex, edge.start_mm),
            (edge.end_vertex, edge.end_mm),
        ] {
            if seen
                .insert(id, point)
                .is_some_and(|old| distance(old, point) > TOL_MM)
            {
                return None;
            }
        }
    }
    let vertices = seen
        .into_iter()
        .map(|(id, point)| Vertex { id, point })
        .collect::<Vec<_>>();
    let level = dot(sub(vertices.first()?.point, origin), axis);
    if vertices
        .iter()
        .any(|vertex| (dot(sub(vertex.point, origin), axis) - level).abs() > TOL_MM)
    {
        return None;
    }
    section(level, &vertices, origin, x, y)
}

/// A source terminal cap may lead into several real four-quad rings before
/// disappearing into a larger solid. Each ring yields an *editable* loft
/// section; no fictitious shoulder/cap is introduced into source provenance.
fn terminal_sectioned_ruled_extension(
    solid: &ReducedSolid,
    cap: &ReducedFace,
) -> Option<RuledEndExtension> {
    if !plane_quad(cap) {
        return None;
    }
    let SurfaceSupport::Plane(plane) = cap.surface else {
        return None;
    };
    let origin = cap.loops.first()?.edges.first()?.start_mm;
    let axis = normalize(plane.normal, DIR_TOL)?;
    let first_edge = cap.loops.first()?.edges.first()?;
    let x = normalize(sub(first_edge.end_mm, first_edge.start_mm), DIR_TOL)?;
    if dot(x, axis).abs() > DIR_TOL {
        return None;
    }
    let y = normalize(cross(axis, x), DIR_TOL)?;
    let mut current = cap.loops.first()?.edges.clone();
    let mut sections = vec![section_from_wall_ring(&current, origin, axis, x, y)?];
    let mut consumed = BTreeSet::from([cap.id]);
    let mut wall_ids = Vec::new();
    let mut direction = 0.0_f64;
    loop {
        let Some((next, walls)) = next_quad_wall_ring(solid, &current, &consumed) else {
            break;
        };
        let Some(stage) = section_from_wall_ring(&next, origin, axis, x, y) else {
            break;
        };
        let previous = sections.last()?.level;
        let delta = stage.level - previous;
        if delta.abs() <= TOL_MM || (direction != 0.0 && delta.signum() != direction.signum()) {
            break;
        }
        direction = delta;
        wall_ids.extend(&walls);
        consumed.extend(walls);
        current = next;
        sections.push(stage);
    }
    // The ordinary terminal extension already owns a single transition.
    // This recognizer is for actual multi-stage walls only.
    if sections.len() < 3 {
        return None;
    }
    let depth = sections.last()?.level;
    if depth < 0.0 {
        sections.reverse();
    }
    let sections = sections
        .into_iter()
        .map(|s| LoftSection {
            profile: s.profile,
            transform: RigidTransform {
                matrix: [
                    [x[0], y[0], axis[0], axis[0].mul_add(s.level, origin[0])],
                    [x[1], y[1], axis[1], axis[1].mul_add(s.level, origin[1])],
                    [x[2], y[2], axis[2], axis[2].mul_add(s.level, origin[2])],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        })
        .collect();
    let mut model = CadModel::new();
    let root = model.add_node(CadNode::Loft { sections });
    model.add_root(root).ok()?;
    model.validate().ok()?;
    Some(RuledEndExtension {
        fragment: CadFragment {
            source: CadFragmentSource::RuledEndExtension {
                solid_id: solid.solid_id(),
                tip_face_id: cap.id,
                wall_face_ids: wall_ids,
            },
            model,
            root,
        },
        tip_face_id: cap.id,
        axis,
        shoulder_point_mm: crate::math3::add(origin, crate::math3::mul(axis, depth)),
        tip_point_mm: origin,
    })
}

/// Nearest non-tip axial station witnessed by source vertices on actual
/// incident wall faces, not by extrapolation of the planar carriers.
fn source_witnessed_terminal_depth(
    solid: &ReducedSolid,
    cap_edges: &[OrientedEdgeUse],
    wall_face_ids: &[u64],
    origin: [f64; 3],
    axis: [f64; 3],
) -> Option<f64> {
    let cap_vertices = cap_edges
        .iter()
        .map(|edge| edge.start_vertex)
        .collect::<BTreeSet<_>>();
    let depth = solid
        .faces()
        .filter(|wall| wall_face_ids.contains(&wall.id))
        .flat_map(ReducedFace::loop_edges)
        .flat_map(|edge| {
            [
                (edge.start_vertex, edge.start_mm),
                (edge.end_vertex, edge.end_mm),
            ]
        })
        .filter(|(id, _)| !cap_vertices.contains(id))
        .filter_map(|(_, point)| {
            let distance = dot(sub(point, origin), axis);
            (distance.abs() > TOL_MM).then_some(distance)
        })
        .min_by(|a, b| a.abs().total_cmp(&b.abs()))?;
    (depth.abs() > TOL_MM).then_some(depth)
}

// A terminal loft may share its side carriers with a preceding extrusion.
// In that case two side faces extend past the shoulder and are not four-edge
// quads. Their analytic planes still determine the exact bounded shoulder.
fn terminal_extended_planar_ruled_extension(
    solid: &ReducedSolid,
    cap: &ReducedFace,
) -> Option<RuledEndExtension> {
    if !plane_quad(cap) {
        return None;
    }
    let SurfaceSupport::Plane(cap_plane) = cap.surface else {
        return None;
    };
    let cap_edges = &cap.loops[0].edges;
    let mut wall_planes = Vec::with_capacity(cap_edges.len());
    let mut wall_face_ids = Vec::with_capacity(cap_edges.len());
    for edge in cap_edges {
        let mut touching = solid.faces().filter(|face| {
            face.id != cap.id && face.loop_edges().any(|used| used.edge_id == edge.edge_id)
        });
        let wall = touching.next()?;
        if touching.next().is_some() || wall_face_ids.contains(&wall.id) {
            return None;
        }
        let SurfaceSupport::Plane(plane) = wall.surface else {
            return None;
        };
        wall_face_ids.push(wall.id);
        wall_planes.push(plane);
    }

    let axis = normalize(cap_plane.normal, DIR_TOL)?;
    let origin = cap_edges.first()?.start_mm;
    let x = normalize(sub(cap_edges.first()?.end_mm, origin), DIR_TOL)?;
    let y = normalize(cross(axis, x), DIR_TOL)?;
    let depth = source_witnessed_terminal_depth(solid, cap_edges, &wall_face_ids, origin, axis)?;
    let shoulder_origin = crate::math3::add(origin, crate::math3::mul(axis, depth));
    let mut shoulder = Vec::with_capacity(4);
    for i in 0..4 {
        let prior = wall_planes[(i + 3) % 4];
        let next = wall_planes[i];
        let point = intersect_three_planes(
            (axis, shoulder_origin),
            (prior.normal, prior.origin_mm),
            (next.normal, next.origin_mm),
        )?;
        if !point.iter().all(|v| v.is_finite()) {
            return None;
        }
        shoulder.push(point);
    }

    let profile = |points: &[[f64; 3]]| {
        Profile2d::polygon(
            points
                .iter()
                .map(|&p| [dot(sub(p, origin), x), dot(sub(p, origin), y)])
                .collect(),
        )
        .ok()
    };
    let tip = cap_edges
        .iter()
        .map(|edge| edge.start_mm)
        .collect::<Vec<_>>();
    let tip_profile = profile(&tip)?;
    let shoulder_profile = profile(&shoulder)?;
    let frame = |level| RigidTransform {
        matrix: [
            [x[0], y[0], axis[0], axis[0].mul_add(level, origin[0])],
            [x[1], y[1], axis[1], axis[1].mul_add(level, origin[1])],
            [x[2], y[2], axis[2], axis[2].mul_add(level, origin[2])],
            [0.0, 0.0, 0.0, 1.0],
        ],
    };
    let tip = LoftSection {
        profile: tip_profile,
        transform: frame(0.0),
    };
    let shoulder = LoftSection {
        profile: shoulder_profile,
        transform: frame(depth),
    };
    let sections = if depth < 0.0 {
        vec![shoulder, tip]
    } else {
        vec![tip, shoulder]
    };
    let mut model = CadModel::new();
    let root = model.add_node(CadNode::Loft { sections });
    model.add_root(root).ok()?;
    model.validate().ok()?;
    Some(RuledEndExtension {
        fragment: CadFragment {
            source: CadFragmentSource::RuledEndExtension {
                solid_id: solid.solid_id(),
                tip_face_id: cap.id,
                wall_face_ids,
            },
            model,
            root,
        },
        tip_face_id: cap.id,
        axis,
        shoulder_point_mm: shoulder_origin,
        tip_point_mm: origin,
    })
}

fn intersect_three_planes(
    a: ([f64; 3], [f64; 3]),
    b: ([f64; 3], [f64; 3]),
    c: ([f64; 3], [f64; 3]),
) -> Option<[f64; 3]> {
    let ab = cross(a.0, b.0);
    let bc = cross(b.0, c.0);
    let ca = cross(c.0, a.0);
    let det = dot(a.0, bc);
    if det.abs() <= DIR_TOL {
        return None;
    }
    let out = crate::math3::add(
        crate::math3::add(
            crate::math3::mul(bc, dot(a.0, a.1)),
            crate::math3::mul(ca, dot(b.0, b.1)),
        ),
        crate::math3::mul(ab, dot(c.0, c.1)),
    );
    Some(crate::math3::mul(out, 1.0 / det))
}

/// Read-only parser witness for any original STEP solid; the same immutable
/// reduced-B-rep recognizer also consumes the kernel-owned native B-rep.
pub fn detect_sectioned_planar_lofts(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Vec<CadFragment> {
    entities
        .iter()
        .filter_map(|e| {
            let record = crate::step_graph::simple_record(e)?;
            (record.name == "MANIFOLD_SOLID_BREP").then_some(crate::step_graph::entity_id(e))
        })
        .filter_map(|id| {
            let reduced = ReducedSolid::from_step(id, entities, index)?;
            detect_sectioned_planar_loft(&reduced)
        })
        .collect()
}

#[cfg(test)]
mod terminal_extension_tests {
    use super::*;
    use crate::step_graph::build_index;
    use monstertruck_io::step::load::ParsedStepSource;

    #[test]
    fn four_section_beveled_bar_is_an_exact_editable_loft() -> anyhow::Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_four_section_beveled_bar.step");
        let source = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &source.data[0].entities;
        let reduced = ReducedSolid::from_step(24, entities, &build_index(entities))
            .ok_or_else(|| anyhow::anyhow!("no source B-rep"))?;
        let recovered = detect_sectioned_planar_loft(&reduced)
            .ok_or_else(|| anyhow::anyhow!("real four-section beveled bar not recovered"))?;
        assert!(matches!(
            recovered.source,
            CadFragmentSource::SolidLoft { solid_id: 24, .. }
        ));
        recovered.model.validate()?;
        #[cfg(feature = "cad-kernel-monstertruck")]
        {
            use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&recovered.model, recovered.root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            let reconstructed = kernel.to_step(&evaluated)?;
            let output = ParsedStepSource::parse(reconstructed.as_bytes())?.exchange;
            let out_entities = &output.data[0].entities;
            let solid_ids = out_entities
                .iter()
                .filter_map(|e| {
                    crate::step_graph::simple_record(e)
                        .filter(|record| record.name == "MANIFOLD_SOLID_BREP")
                        .map(|_| crate::step_graph::entity_id(e))
                })
                .collect::<Vec<_>>();
            let [output_id] = solid_ids.as_slice() else {
                anyhow::bail!("ruled loft did not produce exactly one STEP solid");
            };
            let boundary = crate::boundary_proof::compare_analytic_boundaries(
                24,
                entities,
                *output_id,
                out_entities,
            );
            assert!(
                boundary.equivalent,
                "source vs ruled loft boundary differs: {boundary:?}"
            );
            let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(
                report.unrecovered_solid_count, 0,
                "{:?}",
                report.kernel_recovery_failures
            );
            let [chosen] = report.solid_lofts.as_slice() else {
                anyhow::bail!("complete IR did not select the source-exact four-section loft");
            };
            assert_eq!(chosen.fragment.source, recovered.source);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn source_4a6356_exposes_exact_multistage_terminal_loft() -> anyhow::Result<()> {
        use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
        let bytes = include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");
        let exchange = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &exchange.data[0].entities;
        let solid = ReducedSolid::from_step(171, entities, &build_index(entities))
            .ok_or_else(|| anyhow::anyhow!("source mixed carrier B-rep is missing"))?;
        let candidates = detect_extended_ruled_end_extensions(&solid);
        let end =
            candidates
                .iter()
                .find(|extension| {
                    extension.fragment.model.nodes.iter().any(
                        |node| matches!(node, CadNode::Loft { sections } if sections.len() == 4),
                    )
                })
                .ok_or_else(|| anyhow::anyhow!("failed to find source-derived 4-section lobe"))?;
        let CadFragmentSource::RuledEndExtension { wall_face_ids, .. } = &end.fragment.source
        else {
            anyhow::bail!("unexpected extension provenance");
        };
        assert_eq!(wall_face_ids.len(), 12);
        let kernel = MonstertruckKernel::project_standard();
        let result = kernel.evaluate(&end.fragment.model, end.fragment.root)?;
        assert!(kernel.summarize(&result).is_single_solid());
        if let Ok(output) = std::env::var("STEP_PARTIAL_LOFT_OUTPUT") {
            std::fs::write(output, kernel.to_step(&result)?)?;
        }
        Ok(())
    }

    #[test]
    fn discovers_exact_planar_beveled_ends_of_mixed_analytic_elbow() -> anyhow::Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_beveled_l_hook.step");
        let source = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &source.data[0].entities;
        let reduced = ReducedSolid::from_step(84, entities, &build_index(entities))
            .ok_or_else(|| anyhow::anyhow!("source B-rep is missing"))?;
        let mut ends = detect_ruled_end_extensions(&reduced)
            .into_iter()
            .filter_map(|extension| match extension.fragment.source {
                CadFragmentSource::RuledEndExtension { wall_face_ids, .. } => {
                    Some((extension.tip_face_id, wall_face_ids))
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        ends.sort_unstable_by_key(|entry| entry.0);
        assert_eq!(
            ends.len(),
            2,
            "both source end treatments must be reconstructed"
        );
        assert_eq!(ends[0].0, 330);
        assert_eq!(ends[1].0, 355);
        Ok(())
    }
}
