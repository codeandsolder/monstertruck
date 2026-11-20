//! Editable decomposition of a piecewise-prismatic planar source solid.
//!
//! Between transverse source stations, all remaining support planes contain
//! the extrusion direction. Intersect their finite trims with an interior
//! slice and join by shared STEP edge identity. No inferred sharp corners,
//! faceting or arbitrary cap-pair subset enumeration is required.

use super::{
    CadFragment, CadFragmentSource, add_bounded_union_tree, rebuilt_kernel_boundary_matches_source,
    solid::local_frame_transform,
};
use crate::brep::{CurveSupport, SurfaceSupport};
use crate::cad_ir::{CadModel, CadNode, Profile2d, ProofStatus, Provenance};
use crate::cad_kernel::CadKernel;
use crate::math2::{cross as cross2, distance as distance2, sub as sub2};
use crate::math3::{
    add, canonical_unit_direction, cross, distance, dot, mul, norm, normalize,
    point_to_unit_line_distance, sub,
};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use anyhow::Result;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::{BTreeMap, HashMap};
use step_p21::ast::EntityInstance;

const LENGTH: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIRECTION: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

#[derive(Clone, Copy)]
struct Frame {
    origin: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    axis: [f64; 3],
}

impl Frame {
    fn new(origin: [f64; 3], axis: [f64; 3]) -> Option<Self> {
        let seed_index = (0..3).min_by(|&a, &b| axis[a].abs().total_cmp(&axis[b].abs()))?;
        let mut seed = [0.0; 3];
        seed[seed_index] = 1.0;
        let x = normalize(cross(seed, axis), DIRECTION)?;
        Some(Self {
            origin,
            x,
            y: cross(axis, x),
            axis,
        })
    }

    fn axial(self, point: [f64; 3]) -> f64 { dot(sub(point, self.origin), self.axis) }

    fn project(self, point: [f64; 3]) -> [f64; 2] {
        let local = sub(point, self.origin);
        [dot(local, self.x), dot(local, self.y)]
    }
}

fn planar_source(solid: &ReducedSolid) -> Option<Vec<[f64; 3]>> {
    solid.validate_closed_topology()?;
    let mut normals = Vec::new();
    for face in solid.faces() {
        let SurfaceSupport::Plane(plane) = face.surface else {
            return None;
        };
        let normal = canonical_unit_direction(plane.normal, DIRECTION)?;
        if !normals.iter().any(|&n| norm(cross(n, normal)) <= DIRECTION) {
            normals.push(normal);
        }
        for edge in face.loop_edges() {
            let CurveSupport::Line(line) = edge.support else {
                return None;
            };
            let direction = normalize(line.direction, DIRECTION)?;
            if distance(edge.start_mm, edge.end_mm) <= LENGTH {
                return None;
            }
            for point in [edge.start_mm, edge.end_mm] {
                if !point.iter().all(|v| v.is_finite())
                    || dot(sub(point, plane.origin_mm), normal).abs() > LENGTH
                    || point_to_unit_line_distance(point, line.origin_mm, direction) > LENGTH
                {
                    return None;
                }
            }
        }
    }
    Some(normals)
}

fn stations(solid: &ReducedSolid, frame: Frame, normals: &[[f64; 3]]) -> Option<Vec<f64>> {
    if normals.iter().any(|&normal| {
        dot(normal, frame.axis).abs() > DIRECTION && norm(cross(normal, frame.axis)) > DIRECTION
    }) {
        return None;
    }
    let mut levels = solid
        .faces()
        .flat_map(ReducedFace::loop_edges)
        .flat_map(|edge| [frame.axial(edge.start_mm), frame.axial(edge.end_mm)])
        .collect::<Vec<_>>();
    levels.sort_by(f64::total_cmp);
    levels.dedup_by(|a, b| (*a - *b).abs() <= LENGTH);
    (levels.len() >= 3
        && levels
            .windows(2)
            .all(|pair| pair[1] - pair[0] > 2.0 * LENGTH))
    .then_some(levels)
}

#[derive(Default)]
struct SliceGraph {
    points: BTreeMap<u64, [f64; 2]>,
    next: BTreeMap<u64, u64>,
}

impl SliceGraph {
    fn insert(&mut self, id: u64, at: [f64; 2]) -> Option<()> {
        if self
            .points
            .insert(id, at)
            .is_some_and(|old| distance2(old, at) > LENGTH)
        {
            return None;
        }
        Some(())
    }

    fn add_face(&mut self, face: &ReducedFace, frame: Frame, level: f64) -> Option<()> {
        let SurfaceSupport::Plane(plane) = face.surface else {
            return None;
        };
        let normal = normalize(plane.normal, DIRECTION)?;
        if norm(cross(normal, frame.axis)) <= DIRECTION {
            return Some(());
        }
        let tangent_3d = cross(
            frame.axis,
            mul(normal, if face.surface_forward { 1.0 } else { -1.0 }),
        );
        let tangent = [dot(tangent_3d, frame.x), dot(tangent_3d, frame.y)];
        let mut hits = Vec::new();
        for edge in face.loop_edges() {
            let a = frame.axial(edge.start_mm);
            let b = frame.axial(edge.end_mm);
            if (a - level).abs() <= LENGTH || (b - level).abs() <= LENGTH {
                return None;
            }
            // The slice is strictly between source vertex stations, so no endpoint or coplanar-edge cases occur.
            if (a < level) == (b < level) {
                continue;
            }
            let t = (level - a) / (b - a);
            let at = frame.project(add(edge.start_mm, mul(sub(edge.end_mm, edge.start_mm), t)));
            self.insert(edge.edge_id, at)?;
            let order = crate::math2::dot(at, tangent);
            hits.push((order, edge.edge_id));
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (pairs, remainder) = hits.as_chunks::<2>();
        for &[(a, start), (b, end)] in pairs {
            // Alternating intersections along this plane's finite trim recover all its material intervals, including holes.
            if b - a <= LENGTH || self.next.insert(start, end).is_some() {
                return None;
            }
        }
        remainder.is_empty().then_some(())
    }

    fn rings(mut self) -> Option<Vec<Vec<[f64; 2]>>> {
        if self.points.len() != self.next.len() {
            return None;
        }
        let mut rings = Vec::new();
        while let Some((&start, _)) = self.next.first_key_value() {
            let mut ring = Vec::new();
            let mut at = start;
            loop {
                ring.push(*self.points.get(&at)?);
                at = self.next.remove(&at)?;
                if at == start {
                    break;
                }
            }
            if ring.len() < 3 {
                return None;
            }
            rings.push(ring);
        }
        Some(rings)
    }
}

fn area(ring: &[[f64; 2]]) -> f64 {
    let first = ring[0];
    ring[1..]
        .windows(2)
        .map(|p| cross2(sub2(p[0], first), sub2(p[1], first)))
        .sum::<f64>()
        * 0.5
}

fn contains(ring: &[[f64; 2]], point: [f64; 2]) -> bool {
    let mut inside = false;
    for (&start, &end) in ring.iter().zip(ring.iter().cycle().skip(1)) {
        let [start_x, start_y] = sub2(start, point);
        let [end_x, end_y] = sub2(end, point);
        if (start_y > 0.0) != (end_y > 0.0)
            && (end_x - start_x).mul_add(-start_y / (end_y - start_y), start_x) > 0.0
        {
            inside = !inside;
        }
    }
    inside
}

fn section_profiles(rings: &[Vec<[f64; 2]>]) -> Option<Vec<Profile2d>> {
    let areas = rings.iter().map(|ring| area(ring)).collect::<Vec<_>>();
    if areas
        .iter()
        .any(|a| !a.is_finite() || a.abs() <= LENGTH * LENGTH)
    {
        return None;
    }
    let mut outers = BTreeMap::new();
    for (i, ring) in rings.iter().enumerate() {
        if areas[i] > 0.0 {
            outers.insert(i, Profile2d::polygon(ring.clone()).ok()?);
        }
    }
    for (i, ring) in rings.iter().enumerate() {
        if areas[i] > 0.0 {
            continue;
        }
        // Assign a hole to its smallest containing positive ring. Nested
        // islands are independent positive profiles, not holes of holes.
        let parent = outers
            .keys()
            .copied()
            .filter(|&j| areas[j] > -areas[i] && contains(&rings[j], ring[0]))
            .min_by(|&a, &b| areas[a].total_cmp(&areas[b]))?;
        outers
            .get_mut(&parent)?
            .loops
            .extend(Profile2d::polygon(ring.clone()).ok()?.loops);
    }
    (!outers.is_empty()).then(|| outers.into_values().collect())
}

fn slab_profiles(solid: &ReducedSolid, frame: Frame, level: f64) -> Option<Vec<Profile2d>> {
    let mut graph = SliceGraph::default();
    for face in solid.faces() {
        graph.add_face(face, frame, level)?;
    }
    section_profiles(&graph.rings()?)
}

fn sectional_model(
    solid: &ReducedSolid,
    frame: Frame,
    levels: &[f64],
) -> Result<Option<(CadModel, crate::cad_ir::NodeId)>> {
    let mut model = CadModel::new();
    let mut children = Vec::new();
    for limits in levels.windows(2) {
        let [a, b] = [limits[0], limits[1]];
        let Some(profiles) = slab_profiles(solid, frame, a.midpoint(b)) else {
            return Ok(None);
        };
        let origin = add(frame.origin, mul(frame.axis, a));
        for profile in profiles {
            let body = model.add_node(CadNode::Extrude {
                profile,
                vector_mm: [0.0, 0.0, b - a],
            });
            children.push(model.add_node(CadNode::Transform {
                transform: local_frame_transform(origin, frame.x, frame.y, frame.axis),
                child: body,
            }));
        }
    }
    if children.len() < 2 {
        return Ok(None);
    }
    let root = add_bounded_union_tree(&mut model, children);
    model.add_root(root)?;
    model.validate()?;
    Ok(Some((model, root)))
}

pub(super) fn recover<K: CadKernel>(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    kernel: &K,
) -> Result<Option<CadFragment>> {
    let Some(solid) = ReducedSolid::from_step(solid_id, entities, index) else {
        return Ok(None);
    };
    let Some(normals) = planar_source(&solid) else {
        return Ok(None);
    };
    let Some(origin) = solid
        .faces()
        .flat_map(ReducedFace::loop_edges)
        .next()
        .map(|edge| edge.start_mm)
    else {
        return Ok(None);
    };
    let mut candidates = Vec::new();
    for &axis in &normals {
        let Some(frame) = Frame::new(origin, axis) else {
            continue;
        };
        let Some(levels) = stations(&solid, frame, &normals) else {
            continue;
        };
        if let Some((model, root)) = sectional_model(&solid, frame, &levels)? {
            candidates.push((model.complexity_score(root)?, model, root, levels.len() - 1));
        }
    }
    candidates.sort_by_key(|c| c.0);
    for (_, mut model, root, slabs) in candidates {
        let Ok(shape) = kernel.evaluate(&model, root) else {
            continue;
        };
        if !kernel.summarize(&shape).is_single_solid()
            || !rebuilt_kernel_boundary_matches_source(kernel, &shape, solid_id, entities)
                .unwrap_or(false)
        {
            continue;
        }
        let face_ids = solid.faces().map(|face| face.id).collect::<Vec<_>>();
        let mut source_entity_ids = vec![solid_id];
        source_entity_ids.extend(&face_ids);
        model.set_provenance(
            root,
            Provenance {
                source_entity_ids,
                proof: ProofStatus::WithinTolerance,
                max_residual_mm: Some(crate::boundary_proof::ANALYTIC_BOUNDARY_PROOF_TOL_MM),
            },
        )?;
        return Ok(Some(CadFragment {
            source: CadFragmentSource::SectionalExtrusionUnion {
                solid_id,
                face_ids,
                slabs,
            },
            model,
            root,
        }));
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rectangle(x: f64, y: f64, w: f64, h: f64) -> Vec<[f64; 2]> {
        vec![[x, y], [x + w, y], [x + w, y + h], [x, y + h]]
    }

    #[test]
    fn sectional_profiles_keep_disjoint_islands_and_nested_voids() {
        let outer = rectangle(0.0, 0.0, 10.0, 10.0);
        let mut hole = rectangle(1.0, 1.0, 8.0, 8.0);
        hole.reverse();
        let island = rectangle(2.0, 2.0, 6.0, 6.0);
        let mut inner_hole = rectangle(3.0, 3.0, 4.0, 4.0);
        inner_hole.reverse();
        let disjoint = rectangle(20.0, 20.0, 2.0, 2.0);
        let profiles = section_profiles(&[hole, inner_hole, disjoint, island, outer]).unwrap();
        assert_eq!(profiles.len(), 3);
        let counts = profiles.iter().map(|p| p.loops.len()).collect::<Vec<_>>();
        assert_eq!(counts, [1, 2, 2]);
        let mut uncontained = rectangle(0.0, 0.0, 1.0, 1.0);
        uncontained.reverse();
        assert!(section_profiles(&[uncontained]).is_none());
    }

    #[test]
    fn sectional_slice_graph_refuses_open_or_wrongly_oriented_edges() {
        let mut graph = SliceGraph::default();
        for (i, p) in [(1, [0.0, 0.0]), (2, [1.0, 0.0]), (3, [0.0, 1.0])] {
            graph.insert(i, p).unwrap();
        }
        graph.next.extend([(1, 2), (2, 3), (3, 1)]);
        assert_eq!(graph.rings().unwrap().len(), 1);
        let mut graph = SliceGraph::default();
        graph
            .points
            .extend([(1, [0.0, 0.0]), (2, [1.0, 0.0]), (3, [0.0, 1.0])]);
        graph.next.extend([(1, 2), (2, 1), (3, 2)]);
        assert!(graph.rings().is_none());
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn sectional_real_sources_are_automatically_source_equivalent() -> Result<()> {
        let cases: &[(&str, &[u8], u64, usize)] = &[
            (
                "recess 6cba",
                include_bytes!("../../validation/fixtures/native_sectional_recess_6cbaef14.step"),
                55,
                2,
            ),
            (
                "recess 27af",
                include_bytes!("../../validation/fixtures/native_sectional_recess_27afac28.step"),
                49,
                2,
            ),
            (
                "stepped body 6558",
                include_bytes!("../../validation/fixtures/native_sectional_stepped_655821f7.step"),
                99,
                5,
            ),
        ];
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for &(label, bytes, expected_id, expected_slabs) in cases {
            let source = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
            let entities = &source.data[0].entities;
            let index = crate::step_graph::build_index(entities);
            // The generic sectional path must independently produce the correct
            // body even when an earlier, simpler cap reconstruction now wins.
            let sectioned = recover(expected_id, entities, &index, &kernel)?.ok_or_else(|| {
                anyhow::anyhow!("{label}: no independently proven sectional union")
            })?;
            assert!(
                matches!(sectioned.source,CadFragmentSource::SectionalExtrusionUnion {solid_id,slabs,..} if solid_id==expected_id && slabs==expected_slabs),
                "{label}: {:?}",
                sectioned.source
            );
            let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1, "{label}");
            assert_eq!(
                report.constructively_recovered_solid_count, 1,
                "{label}: {:?}",
                report.kernel_recovery_failures
            );
            let [selected] = report.extrusion_core_unions.as_slice() else {
                anyhow::bail!("{label}: expected one source-derived cross-section union");
            };
            let fragment = &selected.fragment;
            // The two recessed prisms may select a simpler earlier proof.
            // The stepped source is the genuine new sectional admission.
            if expected_id == 99 {
                assert!(
                    matches!(
                        fragment.source,
                        CadFragmentSource::SectionalExtrusionUnion { solid_id: 99, .. }
                    ),
                    "{label}: {:?}",
                    fragment.source
                );
            }
            fragment.model.validate()?;
            let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid(), "{label}");
            assert!(
                rebuilt_kernel_boundary_matches_source(&kernel, &evaluated, expected_id, entities)?,
                "{label}"
            );
        }
        Ok(())
    }
}
