//! Editable decomposition of a piecewise-prismatic analytic source solid.
//!
//! Between transverse source stations, side planes and cylinders contain
//! the extrusion direction. Intersect their finite trims with an interior
//! slice and join by shared STEP edge identity. No inferred sharp corners,
//! faceting or arbitrary cap-pair subset enumeration is required.

use super::{
    CadFragment, CadFragmentSource, add_bounded_union_tree, rebuilt_kernel_boundary_matches_source,
    solid::{local_frame_transform, recovered_profile_curve_to_ir},
};
use crate::brep::{CurveSupport, OrientedEdgeUse, SurfaceSupport};
use crate::cad_ir::{CadModel, CadNode, Profile2d, ProfileLoop, ProofStatus, Provenance};
use crate::cad_kernel::CadKernel;
use crate::math2::{distance as distance2, sub as sub2};
use crate::math3::{
    add, canonical_unit_direction, cross, distance, dot, mul, norm, normalize,
    point_to_unit_line_distance, sub,
};
use crate::profile_curves::RecoveredProfileCurve;
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use crate::solid_extrusions::{signed_area_curves, strict_profile_contains};
use anyhow::Result;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;
use step_p21::ast::EntityInstance;

mod circular;
mod transitions;

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

fn source_axes(solid: &ReducedSolid) -> Option<Vec<[f64; 3]>> {
    solid.validate_closed_topology()?;
    let mut directions = Vec::new();
    for face in solid.faces() {
        let direction = match face.surface {
            SurfaceSupport::Plane(plane) => plane.normal,
            SurfaceSupport::Cylinder(cylinder) => cylinder.axis,
            _ => return None,
        };
        let axis = canonical_unit_direction(direction, DIRECTION)?;
        if !directions
            .iter()
            .any(|&other| norm(cross(other, axis)) <= DIRECTION)
        {
            directions.push(axis);
        }
    }
    Some(directions)
}

fn valid_line(edge: &OrientedEdgeUse) -> bool {
    let CurveSupport::Line(line) = edge.support else {
        return false;
    };
    let Some(direction) = normalize(line.direction, DIRECTION) else {
        return false;
    };
    distance(edge.start_mm, edge.end_mm) > LENGTH
        && [edge.start_mm, edge.end_mm].into_iter().all(|point| {
            point.iter().all(|v| v.is_finite())
                && point_to_unit_line_distance(point, line.origin_mm, direction) <= LENGTH
        })
}

fn valid_planar_face(face: &ReducedFace, frame: Frame) -> Option<()> {
    let SurfaceSupport::Plane(plane) = face.surface else {
        return None;
    };
    let normal = normalize(plane.normal, DIRECTION)?;
    let transverse = norm(cross(normal, frame.axis)) <= DIRECTION;
    if !transverse && dot(normal, frame.axis).abs() > DIRECTION {
        return None;
    }
    for edge in face.loop_edges() {
        if ![edge.start_mm, edge.end_mm]
            .into_iter()
            .all(|point| dot(sub(point, plane.origin_mm), normal).abs() <= LENGTH)
        {
            return None;
        }
        match edge.support {
            CurveSupport::Line(_) if valid_line(edge) => {}
            CurveSupport::Circle(circle) if transverse => {
                circular::transverse_arc(edge, frame)?;
                if dot(sub(circle.center_mm, plane.origin_mm), normal).abs() > LENGTH {
                    return None;
                }
            }
            _ => return None,
        }
    }
    Some(())
}

fn stations(solid: &ReducedSolid, frame: Frame) -> Option<Vec<f64>> {
    for face in solid.faces() {
        match face.surface {
            SurfaceSupport::Plane(_) => valid_planar_face(face, frame)?,
            SurfaceSupport::Cylinder(_) => circular::CylinderTrim::new(face, frame).map(|_| ())?,
            _ => return None,
        }
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
    next: BTreeMap<u64, (u64, RecoveredProfileCurve)>,
    closed: Vec<Vec<RecoveredProfileCurve>>,
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

    fn connect(&mut self, start: u64, end: u64, curve: RecoveredProfileCurve) -> Option<()> {
        if distance2(*self.points.get(&start)?, curve.start_point()?) > LENGTH
            || distance2(*self.points.get(&end)?, curve.end_point()?) > LENGTH
        {
            return None;
        }
        self.next
            .insert(start, (end, curve))
            .is_none()
            .then_some(())
    }

    fn add_plane(&mut self, face: &ReducedFace, frame: Frame, level: f64) -> Option<()> {
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
            if b - a <= LENGTH {
                return None;
            }
            let curve = RecoveredProfileCurve::Line {
                source_edge_ids: vec![start, end],
                start_mm: *self.points.get(&start)?,
                end_mm: *self.points.get(&end)?,
            };
            self.connect(start, end, curve)?;
        }
        remainder.is_empty().then_some(())
    }

    fn rings(mut self) -> Option<Vec<Vec<RecoveredProfileCurve>>> {
        if self.points.len() != self.next.len() {
            return None;
        }
        while let Some((&start, _)) = self.next.first_key_value() {
            let mut ring = Vec::new();
            let mut at = start;
            loop {
                let (end, curve) = self.next.remove(&at)?;
                ring.push(curve);
                at = end;
                if at == start {
                    break;
                }
            }
            if ring.len() < 3
                && ring
                    .iter()
                    .all(|curve| matches!(curve, RecoveredProfileCurve::Line { .. }))
            {
                return None;
            }
            self.closed.push(ring);
        }
        Some(self.closed)
    }
}

/// Pick a ray separated from all endpoints and circular tangencies, then use
/// the same analytic strict inclusion test as the extrusion-core recognizer.
/// No perturbed sample point or polygonal approximation is involved.
fn contains(ring: &[RecoveredProfileCurve], point: [f64; 2]) -> Option<bool> {
    let mut forbidden = Vec::new();
    for curve in ring {
        for vertex in [curve.start_point()?, curve.end_point()?] {
            let relative = sub2(vertex, point);
            if relative[0].hypot(relative[1]) <= LENGTH {
                return None;
            }
            forbidden.push(relative[1].atan2(relative[0]).rem_euclid(PI));
        }
        if let RecoveredProfileCurve::CircleArc {
            center_mm,
            radius_mm,
            ..
        } = curve
        {
            let relative = sub2(*center_mm, point);
            let distance = relative[0].hypot(relative[1]);
            if distance >= *radius_mm {
                let bearing = relative[1].atan2(relative[0]);
                let tangent = (radius_mm / distance).asin();
                forbidden.extend([
                    (bearing - tangent).rem_euclid(PI),
                    (bearing + tangent).rem_euclid(PI),
                ]);
            }
        }
    }
    forbidden.sort_by(f64::total_cmp);
    let first = *forbidden.first()?;
    forbidden.push(first + PI);
    let interval = forbidden
        .windows(2)
        .max_by(|a, b| (a[1] - a[0]).total_cmp(&(b[1] - b[0])))?;
    let angle = interval[0].midpoint(interval[1]);
    let (sin, cos) = angle.sin_cos();
    let rotate = |p| {
        let [x, y] = sub2(p, point);
        [x.mul_add(cos, y * sin), y.mul_add(cos, -x * sin)]
    };
    let rotated = ring
        .iter()
        .map(|curve| match curve {
            RecoveredProfileCurve::Line {
                source_edge_ids,
                start_mm,
                end_mm,
            } => Some(RecoveredProfileCurve::Line {
                source_edge_ids: source_edge_ids.clone(),
                start_mm: rotate(*start_mm),
                end_mm: rotate(*end_mm),
            }),
            RecoveredProfileCurve::CircleArc {
                source_edge_ids,
                center_mm,
                radius_mm,
                start_angle_rad,
                end_angle_rad,
            } => Some(RecoveredProfileCurve::CircleArc {
                source_edge_ids: source_edge_ids.clone(),
                center_mm: rotate(*center_mm),
                radius_mm: *radius_mm,
                start_angle_rad: start_angle_rad - angle,
                end_angle_rad: end_angle_rad - angle,
            }),
            _ => None,
        })
        .collect::<Option<Vec<_>>>()?;
    strict_profile_contains(&rotated, [0.0; 2], LENGTH)
}

fn section_profiles(rings: &[Vec<RecoveredProfileCurve>]) -> Option<Vec<Profile2d>> {
    let areas = rings
        .iter()
        .map(|ring| signed_area_curves(ring))
        .collect::<Vec<_>>();
    if areas
        .iter()
        .any(|a| !a.is_finite() || a.abs() <= LENGTH * LENGTH)
    {
        return None;
    }
    let mut outers = BTreeMap::new();
    for (i, ring) in rings.iter().enumerate() {
        if areas[i] > 0.0 {
            outers.insert(
                i,
                Profile2d {
                    loops: vec![ProfileLoop {
                        curves: ring.iter().map(recovered_profile_curve_to_ir).collect(),
                    }],
                },
            );
        }
    }
    for (i, ring) in rings.iter().enumerate() {
        if areas[i] > 0.0 {
            continue;
        }
        // Assign a hole to its smallest containing positive ring. Nested
        // islands are independent positive profiles, not holes of holes.
        let mut parents = Vec::new();
        let witness = ring.first()?.start_point()?;
        for &j in outers.keys() {
            if areas[j] > -areas[i] && contains(&rings[j], witness)? {
                parents.push(j);
            }
        }
        let parent = parents
            .into_iter()
            .min_by(|&a, &b| areas[a].total_cmp(&areas[b]))?;
        outers.get_mut(&parent)?.loops.push(ProfileLoop {
            curves: ring.iter().map(recovered_profile_curve_to_ir).collect(),
        });
    }
    (!outers.is_empty()).then(|| outers.into_values().collect())
}

fn slab_profiles(solid: &ReducedSolid, frame: Frame, level: f64) -> Option<Vec<Profile2d>> {
    let mut graph = SliceGraph::default();
    for face in solid.faces() {
        match face.surface {
            SurfaceSupport::Plane(_) => graph.add_plane(face, frame, level)?,
            SurfaceSupport::Cylinder(_) => {
                circular::CylinderTrim::new(face, frame)?.add_slice(&mut graph, level)?;
            }
            _ => return None,
        }
    }
    section_profiles(&graph.rings()?)
}

fn extrusion_node(
    model: &mut CadModel,
    frame: Frame,
    a: f64,
    b: f64,
    profile: Profile2d,
) -> crate::cad_ir::NodeId {
    let origin = add(frame.origin, mul(frame.axis, a));
    let body = model.add_node(CadNode::Extrude {
        profile,
        vector_mm: [0.0, 0.0, b - a],
    });
    model.add_node(CadNode::Transform {
        transform: local_frame_transform(origin, frame.x, frame.y, frame.axis),
        child: body,
    })
}

fn sectional_model(
    frame: Frame,
    levels: &[f64],
    sections: &[Vec<Profile2d>],
) -> Result<Option<(CadModel, crate::cad_ir::NodeId)>> {
    let mut model = CadModel::new();
    let mut children = Vec::new();
    for (limits, profiles) in levels.windows(2).zip(sections) {
        let [a, b] = [limits[0], limits[1]];
        let parts = profiles
            .iter()
            .map(|profile| extrusion_node(&mut model, frame, a, b, profile.clone()))
            .collect::<Vec<_>>();
        if !parts.is_empty() {
            children.push(add_bounded_union_tree(&mut model, parts));
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

pub fn recover<K: CadKernel>(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    kernel: &K,
) -> Result<Option<CadFragment>> {
    let Some(solid) = ReducedSolid::from_step(solid_id, entities, index) else {
        return Ok(None);
    };
    let Some(axes) = source_axes(&solid) else {
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
    for axis in axes {
        let Some(frame) = Frame::new(origin, axis) else {
            continue;
        };
        let Some(levels) = stations(&solid, frame) else {
            continue;
        };
        let Some(sections) = levels
            .windows(2)
            .map(|pair| slab_profiles(&solid, frame, pair[0].midpoint(pair[1])))
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };
        // Both representations derive from the same proven analytic sections.
        // Try the simpler editable construction first, including continuous
        // host extrusions with source-witnessed recesses.
        let slabs = levels.len() - 1;
        if let Some((model, root)) =
            transitions::model(&solid, frame, &levels, &sections[0], entities, index)?
        {
            candidates.push((model.complexity_score(root)?, model, root, slabs, true));
        }
        if let Some((model, root)) = sectional_model(frame, &levels, &sections)? {
            candidates.push((model.complexity_score(root)?, model, root, slabs, false));
        }
    }
    candidates.sort_by_key(|c| c.0);
    for (_, mut model, root, slabs, transitions) in candidates {
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
            source: if transitions {
                CadFragmentSource::SectionalExtrusionTransitions {
                    solid_id,
                    face_ids,
                    slabs,
                }
            } else {
                CadFragmentSource::SectionalExtrusionUnion {
                    solid_id,
                    face_ids,
                    slabs,
                }
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

    fn rectangle(x: f64, y: f64, w: f64, h: f64) -> Vec<RecoveredProfileCurve> {
        let points = [[x, y], [x + w, y], [x + w, y + h], [x, y + h]];
        points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .map(|(&start_mm, &end_mm)| RecoveredProfileCurve::Line {
                source_edge_ids: vec![],
                start_mm,
                end_mm,
            })
            .collect()
    }
    fn reverse(ring: &mut [RecoveredProfileCurve]) {
        ring.reverse();
        for curve in ring {
            *curve = curve.reversed();
        }
    }

    #[test]
    fn sectional_profiles_keep_disjoint_islands_and_nested_voids() {
        let outer = rectangle(0.0, 0.0, 10.0, 10.0);
        let mut hole = rectangle(1.0, 1.0, 8.0, 8.0);
        reverse(&mut hole);
        let island = rectangle(2.0, 2.0, 6.0, 6.0);
        let mut inner_hole = rectangle(3.0, 3.0, 4.0, 4.0);
        reverse(&mut inner_hole);
        let disjoint = rectangle(20.0, 20.0, 2.0, 2.0);
        let profiles = section_profiles(&[hole, inner_hole, disjoint, island, outer]).unwrap();
        assert_eq!(profiles.len(), 3);
        let counts = profiles.iter().map(|p| p.loops.len()).collect::<Vec<_>>();
        assert_eq!(counts, [1, 2, 2]);
        let mut uncontained = rectangle(0.0, 0.0, 1.0, 1.0);
        reverse(&mut uncontained);
        assert!(section_profiles(&[uncontained]).is_none());
    }

    #[test]
    fn sectional_slice_graph_refuses_open_or_wrongly_oriented_edges() {
        let mut graph = SliceGraph::default();
        for (i, p) in [(1, [0.0, 0.0]), (2, [1.0, 0.0]), (3, [0.0, 1.0])] {
            graph.insert(i, p).unwrap();
        }
        for (a, b) in [(1, 2), (2, 3), (3, 1)] {
            graph
                .connect(
                    a,
                    b,
                    RecoveredProfileCurve::Line {
                        source_edge_ids: vec![],
                        start_mm: graph.points[&a],
                        end_mm: graph.points[&b],
                    },
                )
                .unwrap();
        }
        assert_eq!(graph.rings().unwrap().len(), 1);
        let mut graph = SliceGraph::default();
        graph
            .points
            .extend([(1, [0.0, 0.0]), (2, [1.0, 0.0]), (3, [0.0, 1.0])]);
        for (a, b) in [(1, 2), (2, 1), (3, 2)] {
            graph
                .connect(
                    a,
                    b,
                    RecoveredProfileCurve::Line {
                        source_edge_ids: vec![],
                        start_mm: graph.points[&a],
                        end_mm: graph.points[&b],
                    },
                )
                .unwrap();
        }
        assert!(graph.rings().is_none());
    }

    #[test]
    fn sectional_circular_holes_remain_exact_curves() {
        let outer = rectangle(-4.0, -4.0, 8.0, 8.0);
        let circle = vec![RecoveredProfileCurve::CircleArc {
            source_edge_ids: vec![],
            center_mm: [0.0, 0.0],
            radius_mm: 2.0,
            start_angle_rad: 0.0,
            end_angle_rad: -std::f64::consts::TAU,
        }];
        let island = rectangle(-0.5, -0.5, 1.0, 1.0);
        assert_eq!(contains(&circle, [0.0, 0.0]), Some(true));
        assert_eq!(contains(&circle, [3.0, 0.0]), Some(false));
        assert_eq!(contains(&circle, [2.0, 0.0]), None);
        let profiles = section_profiles(&[outer, circle, island]).unwrap();
        assert_eq!(profiles.len(), 2);
        assert_eq!(profiles[0].loops.len(), 2);
        assert!(matches!(
            profiles[0].loops[1].curves.as_slice(),
            [crate::cad_ir::Curve2d::CircleArc { .. }]
        ));
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn sectional_rounded_real_sources_are_automatically_source_equivalent() -> Result<()> {
        use crate::cad_ir::Curve2d;
        let cases: &[(&str, &[u8], u64)] = &[
            (
                "20e5",
                include_bytes!("../../validation/fixtures/native_sectional_rounded_20e5585e.step"),
                138,
            ),
            (
                "9a9d",
                include_bytes!("../../validation/fixtures/native_sectional_rounded_9a9d454d.step"),
                152,
            ),
        ];
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        for &(label, bytes, id) in cases {
            let source = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
            let entities = &source.data[0].entities;
            let index = crate::step_graph::build_index(entities);
            let fragment = recover(id, entities, &index, &kernel)?
                .ok_or_else(|| anyhow::anyhow!("{label}: no sectional reconstruction"))?;
            assert!(fragment.model.nodes.iter().any(|node|matches!(node,CadNode::Extrude {profile,..} if profile.loops.iter().flat_map(|l|&l.curves).any(|c|matches!(c,Curve2d::CircleArc {..})))));
            let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(
                report.constructively_recovered_solid_count, 1,
                "{label}: {:?}",
                report.kernel_recovery_failures
            );
            let [selected] = report.extrusion_core_unions.as_slice() else {
                anyhow::bail!("{label}: expected one recovered union");
            };
            let fragment = &selected.fragment;
            let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            assert!(rebuilt_kernel_boundary_matches_source(
                &kernel, &evaluated, id, entities
            )?);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn sectional_real_sources_are_automatically_source_equivalent() -> Result<()> {
        let cases: &[(&str, &[u8], u64, usize)] = &[
            (
                "multiple planar transitions b906",
                include_bytes!(
                    "../../validation/fixtures/native_sectional_transitions_b906e9ab.step"
                ),
                146,
                5,
            ),
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
                matches!(sectioned.source,CadFragmentSource::SectionalExtrusionUnion {solid_id,slabs,..} | CadFragmentSource::SectionalExtrusionTransitions {solid_id,slabs,..} if solid_id==expected_id && slabs==expected_slabs),
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
                            | CadFragmentSource::SectionalExtrusionTransitions { solid_id: 99, .. }
                    ),
                    "{label}: {:?}",
                    fragment.source
                );
            }
            if expected_id == 146 {
                assert!(
                    matches!(
                        fragment.source,
                        CadFragmentSource::SectionalExtrusionTransitions { solid_id: 146, .. }
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
