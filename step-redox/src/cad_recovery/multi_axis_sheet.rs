//! Recover a formed sheet from maximal parallel bend groups and the remaining
//! planar patches. A group can have branching flanges, multiple coaxial-radius pairs,
//! and partitioned end caps; a complete common transverse cap is an exact
//! analytic extrusion proposal for the group. All proposals remain subject to
//! the independent whole-solid reconstruction proof.
use super::{
    CadFragment, CadFragmentSource, compose_union_models, rebuilt_kernel_boundary_matches_source,
    solid::{local_frame_transform, recover_source_cap_extrusion_model},
};
use crate::brep::{CurveSupport, SurfaceSupport};
use crate::cad_ir::{CadModel, CadNode, NodeId};
use crate::cad_kernel::CadKernel;
use crate::formed_sheet::{self, FormedSheetEvidence, SheetCylinderPair};
use crate::math3::{add, cross, dot, mul, norm, normalize, sub};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use crate::solid_extrusions::{
    RecoveredSourceCapExtrusion, cap_contains, signed_area_curves,
    source_cap_extrusion_proposal_at_station,
};
use anyhow::Result;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::{HashMap, HashSet};
use step_p21::ast::EntityInstance;

const LENGTH: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIRECTION: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

struct BendGroup<'a> {
    axis: [f64; 3],
    pairs: Vec<&'a SheetCylinderPair>,
}

fn groups(evidence: &FormedSheetEvidence) -> Option<Vec<BendGroup<'_>>> {
    let mut groups: Vec<BendGroup<'_>> = Vec::new();
    for pair in &evidence.cylinder_pairs {
        let axis = normalize(pair.axis, DIRECTION)?;
        if let Some(group) = groups
            .iter_mut()
            .find(|g| norm(cross(g.axis, axis)) <= DIRECTION)
        {
            group.pairs.push(pair);
        } else {
            groups.push(BendGroup {
                axis,
                pairs: vec![pair],
            });
        }
    }
    (groups.len() > 1).then_some(groups)
}

fn face_station(face: &ReducedFace, axis: [f64; 3]) -> Option<[f64; 2]> {
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for edge in face.loop_edges() {
        for point in [edge.start_mm, edge.end_mm] {
            let at = dot(point, axis);
            if !at.is_finite() {
                return None;
            }
            low = low.min(at);
            high = high.max(at);
        }
    }
    (high >= low).then_some([low, high])
}

fn group_cap(
    solid: &ReducedSolid,
    group: &BendGroup<'_>,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<RecoveredSourceCapExtrusion> {
    let axis = group.axis;
    let face_ids = group
        .pairs
        .iter()
        .flat_map(|pair| [pair.inner_face_id, pair.outer_face_id])
        .collect::<Vec<_>>();
    let mut edge_sets = Vec::new();
    let mut common: Option<[f64; 2]> = None;
    for face_id in face_ids {
        let face = solid.faces().find(|f| f.id == face_id)?;
        let [low, high] = face_station(face, axis)?;
        if high - low <= LENGTH {
            return None;
        }
        if let Some([a, b]) = common {
            if (low - a).abs() > LENGTH || (high - b).abs() > LENGTH {
                return None;
            }
        } else {
            common = Some([low, high]);
        }
        // A complete cap must witness the circular trim of every member,
        // not simply touch one outer radius at a shared vertex.
        let edges = face
            .loop_edges()
            .filter(|e| matches!(e.support, CurveSupport::Circle(_)))
            .map(|e| e.edge_id)
            .collect::<HashSet<_>>();
        if edges.is_empty() {
            return None;
        }
        edge_sets.push(edges);
    }
    let [low, high] = common?;
    let mut candidates = Vec::new();
    for face in solid.faces() {
        let SurfaceSupport::Plane(plane) = face.surface else {
            continue;
        };
        let normal = normalize(plane.normal, DIRECTION)?;
        if norm(cross(normal, axis)) > DIRECTION {
            continue;
        }
        let station = dot(plane.origin_mm, axis);
        let inward = if (station - low).abs() <= LENGTH {
            axis
        } else if (station - high).abs() <= LENGTH {
            mul(axis, -1.0)
        } else {
            continue;
        };
        let outward = mul(normal, if face.surface_forward { 1.0 } else { -1.0 });
        if dot(outward, inward) > -1.0 + DIRECTION {
            continue;
        }
        let cap_edges = face.loop_edges().map(|e| e.edge_id).collect::<HashSet<_>>();
        if !edge_sets.iter().all(|edges| !edges.is_disjoint(&cap_edges)) {
            continue;
        }
        let Some(cap) = source_cap_extrusion_proposal_at_station(
            solid.solid_id(),
            face.id,
            inward,
            high - low,
            entities,
            index,
        ) else {
            continue;
        };
        let area = signed_area_curves(&cap.profile_curves).abs()
            - cap
                .inner_profile_loops
                .iter()
                .map(|ring| signed_area_curves(ring).abs())
                .sum::<f64>();
        if area > LENGTH * LENGTH {
            candidates.push((area, cap));
        }
    }
    candidates
        .into_iter()
        .max_by(|a, b| {
            a.0.total_cmp(&b.0)
                .then_with(|| b.1.cap_face_id.cmp(&a.1.cap_face_id))
        })
        .map(|(_, cap)| cap)
}

fn remaining_plane_caps(
    solid: &ReducedSolid,
    evidence: &FormedSheetEvidence,
    bend_caps: &[RecoveredSourceCapExtrusion],
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<Vec<RecoveredSourceCapExtrusion>> {
    let mut caps = Vec::new();
    for pair in &evidence.plane_pairs {
        let negative = solid.faces().find(|f| f.id == pair.negative_face_id)?;
        let positive = solid.faces().find(|f| f.id == pair.positive_face_id)?;
        // Containment of corner witnesses is a proposal-selection heuristic,
        // not a material proof. A missing patch cannot pass full boundary equality.
        if [negative, positive]
            .into_iter()
            .flat_map(ReducedFace::loop_edges)
            .flat_map(|e| [e.start_mm, e.end_mm])
            .all(|p| bend_caps.iter().any(|cap| cap_contains(cap, p)))
        {
            continue;
        }
        let inward = normalize(pair.normal, DIRECTION)?;
        let cap = source_cap_extrusion_proposal_at_station(
            solid.solid_id(),
            negative.id,
            inward,
            pair.separation_mm,
            entities,
            index,
        )?;
        caps.push(cap);
    }
    Some(caps)
}

fn recover_caps(
    evidence: &FormedSheetEvidence,
    solid: &ReducedSolid,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<Vec<RecoveredSourceCapExtrusion>> {
    let mut groups = groups(evidence)?;
    if evidence.paired_cylindrical_faces != evidence.cylindrical_faces
        || !evidence.patch_graph_connected
    {
        return None;
    }
    // Start on a simple bend group; attach its planar flange before more
    // complicated branched groups. Prefer contact continuity over face IDs.
    groups.sort_by_key(|g| g.pairs.len());
    let mut bend_caps = groups
        .iter()
        .map(|g| group_cap(solid, g, entities, index))
        .collect::<Option<Vec<_>>>()?;
    let planar = remaining_plane_caps(solid, evidence, &bend_caps, entities, index)?;
    let remaining = bend_caps.split_off(1);
    bend_caps.extend(planar);
    bend_caps.extend(remaining);
    Some(bend_caps)
}

#[derive(Clone, Copy, Debug)]
struct AssemblyFrame {
    origin: [f64; 3],
    axes: [[f64; 3]; 3],
}

impl AssemblyFrame {
    fn direction(self, v: [f64; 3]) -> [f64; 3] { self.axes.map(|axis| dot(v, axis)) }

    fn cap(self, cap: &RecoveredSourceCapExtrusion) -> RecoveredSourceCapExtrusion {
        let mut cap = cap.clone();
        cap.origin_mm = self.direction(sub(cap.origin_mm, self.origin));
        cap.x_axis = self.direction(cap.x_axis);
        cap.y_axis = self.direction(cap.y_axis);
        cap.z_axis = self.direction(cap.z_axis);
        cap
    }
}

/// A common assembly frame determined by the terminal and branch directions,
/// not by world coordinate signs. Keep the original exact local sketches.
fn assembly_frame(
    caps: &[RecoveredSourceCapExtrusion],
    solid: &ReducedSolid,
) -> Option<AssemblyFrame> {
    let first = caps.first()?;
    let last = caps.last()?;
    let transverse = cross(first.z_axis, last.z_axis);
    let mut y = normalize(transverse, DIRECTION)?;
    // The first uncovered flat flange supplies the free-end direction.
    let terminal = caps.get(1)?;
    let center = |cap: &RecoveredSourceCapExtrusion| {
        add(cap.origin_mm, mul(cap.z_axis, cap.height_mm * 0.5))
    };
    if dot(sub(center(terminal), center(last)), y) < 0.0 {
        y = mul(y, -1.0);
    }
    let mut z = last.z_axis;
    if dot(sub(center(last), center(first)), z) < 0.0 {
        z = mul(z, -1.0);
    }
    let x = normalize(cross(y, z), DIRECTION)?;
    let y = cross(z, x);
    let axes = [x, y, z];
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for edge in solid.faces().flat_map(ReducedFace::loop_edges) {
        for point in [edge.start_mm, edge.end_mm] {
            for (i, axis) in axes.iter().enumerate() {
                let value = dot(point, *axis);
                min[i] = min[i].min(value);
                max[i] = max[i].max(value);
            }
        }
    }
    let origin = axes
        .into_iter()
        .zip(min.into_iter().zip(max))
        .fold([0.0; 3], |a, (axis, (lo, hi))| {
            add(a, mul(axis, lo.midpoint(hi)))
        });
    Some(AssemblyFrame { origin, axes })
}

fn compose_caps(
    caps: &[RecoveredSourceCapExtrusion],
    frame: Option<AssemblyFrame>,
) -> Result<(CadModel, NodeId)> {
    let models = caps
        .iter()
        .map(|cap| {
            frame.map_or_else(
                || recover_source_cap_extrusion_model(cap),
                |frame| recover_source_cap_extrusion_model(&frame.cap(cap)),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let parts = models
        .iter()
        .map(|(model, root)| (model, *root))
        .collect::<Vec<_>>();
    let (mut model, mut root) = compose_union_models(&parts)?;
    if let Some(frame) = frame {
        let [x, y, z] = frame.axes;
        root = model.add_node(CadNode::Transform {
            child: root,
            transform: local_frame_transform(frame.origin, x, y, z),
        });
        model.roots.clear();
        model.add_root(root)?;
    }
    model.validate()?;
    Ok((model, root))
}

pub(super) fn recover<K: CadKernel>(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    kernel: &K,
) -> Result<Option<CadFragment>> {
    let Some(evidence) =
        formed_sheet::detect_formed_sheet_evidence_with_min_pairs(entities, index, 2, 1)
            .into_iter()
            .find(|e| e.solid_id == solid_id)
    else {
        return Ok(None);
    };
    let Some(solid) = ReducedSolid::from_step(solid_id, entities, index) else {
        return Ok(None);
    };
    if solid.validate_closed_topology().is_none() {
        return Ok(None);
    }
    let Some(caps) = recover_caps(&evidence, &solid, entities, index) else {
        return Ok(None);
    };
    let frame = assembly_frame(&caps, &solid);
    let mut selected = None;
    for frame in std::iter::once(None).chain(frame.map(Some)) {
        let (model, root) = compose_caps(&caps, frame)?;
        let Ok(shape) = kernel.evaluate(&model, root) else {
            continue;
        };
        if kernel.summarize(&shape).is_single_solid()
            && rebuilt_kernel_boundary_matches_source(kernel, &shape, solid_id, entities)
                .unwrap_or(false)
        {
            selected = Some((model, root));
            break;
        }
    }
    let Some((mut model, root)) = selected else {
        return Ok(None);
    };
    let face_ids = super::set_proven_source_boundary(&mut model, root, &solid)?;
    Ok(Some(CadFragment {
        source: CadFragmentSource::FormedSheetExtrusionUnion {
            solid_id,
            face_ids,
            cap_face_ids: caps.iter().map(|cap| cap.cap_face_id).collect(),
        },
        model,
        root,
    }))
}
#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;
    #[test]
    fn multi_axis_sheet_source_recovery() -> Result<()> {
        let kernel = crate::cad_kernel::monstertruck::MonstertruckKernel::project_standard();
        let cases: &[(&[u8], u64)] = &[
            (
                include_bytes!("../../validation/fixtures/native_multiaxis_sheet_67a9092b.step"),
                217,
            ),
            (
                include_bytes!("../../validation/fixtures/native_multiaxis_sheet_844ae71b.step"),
                190,
            ),
        ];
        for &(bytes, id) in cases {
            let source = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
            let entities = &source.data[0].entities;
            let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(
                report.constructively_recovered_solid_count, 1,
                "source #{id}: {:?}",
                report.kernel_recovery_failures
            );
            let [selected] = report.extrusion_core_unions.as_slice() else {
                anyhow::bail!("expected one grouped-bend extrusion union for source #{id}")
            };
            let fragment = &selected.fragment;
            assert!(
                matches!(&fragment.source,CadFragmentSource::FormedSheetExtrusionUnion{solid_id,cap_face_ids,..} if *solid_id==id && cap_face_ids.len()==3),
                "{:?}",
                fragment.source
            );
            fragment.model.validate()?;
            assert_eq!(
                fragment
                    .model
                    .nodes
                    .iter()
                    .filter(|n| matches!(n, CadNode::Extrude { .. }))
                    .count(),
                3
            );
            let evaluated = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            assert!(rebuilt_kernel_boundary_matches_source(
                &kernel, &evaluated, id, entities
            )?);
        }
        Ok(())
    }
}
