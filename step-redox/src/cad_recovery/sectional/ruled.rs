//! Piecewise ruled reconstruction of planar polyhedra.
//!
//! Between consecutive vertex stations true source edges cross the whole slab.
//! Their identities give a one-to-one correspondence of section vertices; the
//! same planar source face bounds each ruled edge at both slice endpoints.

use super::{
    DIRECTION, Frame, LENGTH, add_bounded_union_tree, extrusion_node, local_frame_transform,
    section_groups, slice_rings, station_levels, valid_planar_face,
};
use crate::brep::{OrientedEdgeUse, SurfaceSupport};
use crate::cad_ir::{BooleanOp, CadModel, CadNode, LoftSection, NodeId, Profile2d};
use crate::math2::distance;
use crate::math3::{add, cross, dot, mul, norm, normalize, sub};
use crate::profile_curves::RecoveredProfileCurve;
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use anyhow::Result;
use std::collections::BTreeMap;

fn point_at(edge: &OrientedEdgeUse, frame: Frame, level: f64) -> Option<[f64; 2]> {
    let a = frame.axial(edge.start_mm);
    let b = frame.axial(edge.end_mm);
    if level < a.min(b) - LENGTH || level > a.max(b) + LENGTH || (a - b).abs() <= LENGTH {
        return None;
    }
    let t = (level - a) / (b - a);
    Some(frame.project(add(edge.start_mm, mul(sub(edge.end_mm, edge.start_mm), t))))
}

fn section_points(
    ring: &[RecoveredProfileCurve],
    edges: &BTreeMap<u64, &OrientedEdgeUse>,
    frame: Frame,
    level: f64,
) -> Option<Vec<[f64; 2]>> {
    let points = ring
        .iter()
        .map(|curve| {
            let RecoveredProfileCurve::Line {
                source_edge_ids, ..
            } = curve
            else {
                return None;
            };
            let [start, _] = source_edge_ids.as_slice() else {
                return None;
            };
            point_at(edges.get(start)?, frame, level)
        })
        .collect::<Option<Vec<_>>>()?;
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .all(|(&a, &b)| distance(a, b) > LENGTH)
        .then_some(points)
}

struct RuledProfile {
    bottom: Vec<[f64; 2]>,
    top: Vec<[f64; 2]>,
}

impl RuledProfile {
    fn node(self, model: &mut CadModel, frame: Frame, limits: [f64; 2]) -> Result<NodeId> {
        let [a, b] = limits;
        if self
            .bottom
            .iter()
            .zip(&self.top)
            .all(|(&lo, &hi)| distance(lo, hi) <= LENGTH)
        {
            return Ok(extrusion_node(
                model,
                frame,
                a,
                b,
                Profile2d::polygon(self.bottom)?,
            ));
        }
        let sections = [self.bottom, self.top]
            .into_iter()
            .zip(limits)
            .map(|(points, level)| {
                Ok(LoftSection {
                    profile: Profile2d::polygon(points)?,
                    transform: local_frame_transform(
                        add(frame.origin, mul(frame.axis, level)),
                        frame.x,
                        frame.y,
                        frame.axis,
                    ),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(model.add_node(CadNode::Loft { sections }))
    }
}

fn slab(
    solid: &ReducedSolid,
    edges: &BTreeMap<u64, &OrientedEdgeUse>,
    frame: Frame,
    limits: [f64; 2],
) -> Option<Vec<Vec<RuledProfile>>> {
    let rings = slice_rings(solid, frame, limits[0].midpoint(limits[1]))?;
    section_groups(&rings)?
        .into_iter()
        .map(|group| {
            group
                .into_iter()
                .enumerate()
                .map(|(i, index)| {
                    let mut bottom = section_points(&rings[index], edges, frame, limits[0])?;
                    let mut top = section_points(&rings[index], edges, frame, limits[1])?;
                    // Hole rings become positive cutter profiles, with matching vertex order.
                    if i != 0 {
                        bottom.reverse();
                        top.reverse();
                    }
                    Some(RuledProfile { bottom, top })
                })
                .collect()
        })
        .collect()
}

pub(super) fn model(
    solid: &ReducedSolid,
    frame: Frame,
) -> Result<Option<(CadModel, NodeId, usize)>> {
    let mut inclined = false;
    let mut transverse = 0;
    for face in solid.faces() {
        let SurfaceSupport::Plane(plane) = face.surface else {
            return Ok(None);
        };
        let Some(normal) = normalize(plane.normal, DIRECTION) else {
            return Ok(None);
        };
        if norm(cross(normal, frame.axis)) <= DIRECTION {
            transverse += 1;
        } else if dot(normal, frame.axis).abs() > DIRECTION {
            inclined = true;
        }
        if valid_planar_face(face, frame, true).is_none() {
            return Ok(None);
        }
        // The ruled construction has only straight generators and section edges.
        if face.loop_edges().any(|edge| !super::valid_line(edge)) {
            return Ok(None);
        }
    }
    if !inclined || transverse < 2 {
        return Ok(None);
    }
    let Some(levels) = station_levels(solid, frame) else {
        return Ok(None);
    };
    let edges = solid
        .faces()
        .flat_map(ReducedFace::loop_edges)
        .map(|edge| (edge.edge_id, edge))
        .collect::<BTreeMap<_, _>>();
    let mut model = CadModel::new();
    let mut children = Vec::new();
    for pair in levels.windows(2) {
        let limits = [pair[0], pair[1]];
        let Some(groups) = slab(solid, &edges, frame, limits) else {
            return Ok(None);
        };
        let parts = groups
            .into_iter()
            .map(|group| {
                let mut nodes = group
                    .into_iter()
                    .map(|profile| profile.node(&mut model, frame, limits))
                    .collect::<Result<Vec<_>>>()?;
                Ok(if nodes.len() == 1 {
                    nodes.remove(0)
                } else {
                    model.add_node(CadNode::Boolean {
                        op: BooleanOp::Difference,
                        children: nodes,
                    })
                })
            })
            .collect::<Result<Vec<_>>>()?;
        children.push(add_bounded_union_tree(&mut model, parts));
    }
    let root = add_bounded_union_tree(&mut model, children);
    model.add_root(root)?;
    model.validate()?;
    Ok(Some((model, root, levels.len() - 1)))
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
    use crate::complete_ir::recover_complete_ir_with_kernel_bytes;
    use anyhow::Result;

    #[test]
    fn inclined_relief_host_recovery_is_source_verified() -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let bytes =
            include_bytes!("../../../validation/fixtures/native_relief_oblique_8ddfec04.step");
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(report.unrecovered_solid_count, 0);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        let fragments = report.selected_fragments().collect::<Vec<_>>();
        assert_eq!(fragments.len(), 1);
        let fragment = &fragments[0].fragment;
        let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
        let summary = kernel.summarize(&evaluated);
        assert!(summary.geometrically_consistent);
        assert_eq!(summary.shells, 1);
        Ok(())
    }
}
