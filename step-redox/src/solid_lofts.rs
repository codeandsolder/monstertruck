//! Exact recognition of planar, sectioned polyhedra as editable ruled lofts.
//!
//! A section is inferred from actual source `EDGE_CURVE` vertices. No faceting
//! or interpolation is used: every claimed loft wall is a matching source
//! planar quadrilateral with identical topology and geometric vertices.

use crate::brep::{CurveSupport, SurfaceSupport};
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

/// An exact planar end treatment, whose terminal face and four ruled walls
/// are real source topology. The shoulder section is an implied interface,
/// never an invented source face.
pub struct RuledEndExtension {
    pub fragment: CadFragment,
    pub tip_face_id: u64,
    pub axis: [f64; 3],
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
                wall_face_ids,
            },
            model,
            root,
        },
        tip_face_id: cap.id,
        axis,
    })
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
