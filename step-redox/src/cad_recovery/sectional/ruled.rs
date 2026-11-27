//! Analytic swept slabs with linearly changing planar edges and fixed arcs. Within consecutive
//! vertex stations, each finite face intersects a slice in the same ordered
//! intervals. The mid-slice graph supplies correspondence; the actual
//! original crossing edges supply both end profiles, including collapsed
//! corner transitions. Axis-parallel cylinders have invariant section arcs. Independent full-source proof remains mandatory.
use super::{
    DIRECTION, Frame, LENGTH, circular::CylinderTrim, extrusion_node, section_groups, slice_rings,
    station_levels, valid_planar_face,
};
use crate::{
    brep::{OrientedEdgeUse, SurfaceSupport},
    cad_ir::{CadModel, CadNode, Curve2d, LoftSection, NodeId, Profile2d, ProfileLoop},
    cad_recovery::{
        add_bounded_union_tree,
        solid::{local_frame_transform, recovered_profile_curve_to_ir},
    },
    math2::distance as distance2,
    math3::{add, cross, dot, mul, norm, normalize, sub},
    profile_curves::RecoveredProfileCurve,
    reduced_brep::{ReducedFace, ReducedSolid},
};
use anyhow::Result;
use std::collections::BTreeMap;

struct RuledSweep<'a> {
    edges: BTreeMap<u64, &'a OrientedEdgeUse>,
    levels: Vec<f64>,
}

impl<'a> RuledSweep<'a> {
    fn new(solid: &'a ReducedSolid, frame: Frame) -> Option<Self> {
        let mut inclined = false;
        for face in solid.faces() {
            match face.surface {
                SurfaceSupport::Plane(plane) => {
                    let normal = normalize(plane.normal, DIRECTION)?;
                    inclined |= norm(cross(normal, frame.axis)) > DIRECTION
                        && dot(normal, frame.axis).abs() > DIRECTION;
                    valid_planar_face(face, frame, true)?;
                }
                SurfaceSupport::Cylinder(_) => {
                    CylinderTrim::new(face, frame)?;
                }
                _ => return None,
            }
        }
        if !inclined {
            return None;
        }
        let levels = station_levels(solid, frame)?;
        let edges = solid
            .faces()
            .flat_map(ReducedFace::loop_edges)
            .map(|edge| (edge.edge_id, edge))
            .collect();
        Some(Self { edges, levels })
    }

    fn point(&self, id: u64, frame: Frame, level: f64) -> Option<[f64; 2]> {
        let edge = self.edges.get(&id)?;
        let [a, b] = [edge.start_mm, edge.end_mm].map(|p| frame.axial(p));
        if (a - b).abs() <= LENGTH || level < a.min(b) - LENGTH || level > a.max(b) + LENGTH {
            return None;
        }
        let point = if (a - level).abs() <= LENGTH {
            edge.start_mm
        } else if (b - level).abs() <= LENGTH {
            edge.end_mm
        } else {
            add(
                edge.start_mm,
                mul(sub(edge.end_mm, edge.start_mm), (level - a) / (b - a)),
            )
        };
        Some(frame.project(point))
    }

    fn profile(
        &self,
        rings: &[Vec<RecoveredProfileCurve>],
        group: &[usize],
        frame: Frame,
        level: f64,
    ) -> Option<Profile2d> {
        let loops = group
            .iter()
            .map(|&index| {
                let ring = &rings[index];
                let mut points = ring
                    .iter()
                    .map(|curve| match curve {
                        RecoveredProfileCurve::Line {
                            source_edge_ids, ..
                        } => {
                            let [start, _] = source_edge_ids.as_slice() else {
                                return None;
                            };
                            self.point(*start, frame, level)
                        }
                        RecoveredProfileCurve::CircleArc { .. } => curve.start_point(),
                        _ => None,
                    })
                    .collect::<Option<Vec<_>>>()?;
                if points.len() >= 3 {
                    canonicalize_collapsed_corners(&mut points)?;
                }
                let curves = ring
                    .iter()
                    .zip(points.iter().zip(points.iter().cycle().skip(1)))
                    .map(|(curve, (&start_mm, &end_mm))| match curve {
                        RecoveredProfileCurve::Line {
                            source_edge_ids, ..
                        } => {
                            let [_, end] = source_edge_ids.as_slice() else {
                                return None;
                            };
                            (distance2(end_mm, self.point(*end, frame, level)?) <= LENGTH)
                                .then_some(Curve2d::Line { start_mm, end_mm })
                        }
                        RecoveredProfileCurve::CircleArc { .. } => {
                            (distance2(start_mm, curve.start_point()?) <= LENGTH
                                && distance2(end_mm, curve.end_point()?) <= LENGTH)
                                .then(|| recovered_profile_curve_to_ir(curve))
                        }
                        _ => None,
                    })
                    .collect::<Option<Vec<_>>>()?;
                Some(ProfileLoop { curves })
            })
            .collect::<Option<Vec<_>>>()?;
        Some(Profile2d { loops })
    }
}

/// Canonicalize contiguous limiting corners, using the original point
/// distances rather than the progressively modified coordinates. Every
/// displacement must fit the existing length budget, even along a chain
/// of individually close points.
fn canonicalize_collapsed_corners(points: &mut [[f64; 2]]) -> Option<()> {
    let count = points.len();
    if count < 3 {
        return None;
    }
    let offset = (0..count)
        .find(|&i| distance2(points[i], points[(i + 1) % count]) > LENGTH)
        .map_or(0, |i| (i + 1) % count);
    points.rotate_left(offset);
    let mut start = 0;
    while start < count {
        let mut end = start + 1;
        while end < count && distance2(points[end - 1], points[end]) <= LENGTH {
            end += 1;
        }
        let group = &mut points[start..end];
        let canonical = *group
            .iter()
            .min_by(|a, b| a[0].total_cmp(&b[0]).then_with(|| a[1].total_cmp(&b[1])))?;
        if group.iter().any(|&p| distance2(p, canonical) > LENGTH) {
            return None;
        }
        group.fill(canonical);
        start = end;
    }
    points.rotate_right(offset);
    Some(())
}

fn slab_nodes(
    model: &mut CadModel,
    sweep: &RuledSweep<'_>,
    solid: &ReducedSolid,
    frame: Frame,
    lo: f64,
    hi: f64,
) -> Option<Vec<NodeId>> {
    let rings = slice_rings(solid, frame, lo.midpoint(hi))?;
    if rings.is_empty() {
        return Some(Vec::new());
    }
    let groups = section_groups(&rings)?;
    groups
        .iter()
        .map(|group| {
            let start = sweep.profile(&rings, group, frame, lo)?;
            let end = sweep.profile(&rings, group, frame, hi)?;
            // Constant slabs stay ordinary compact extrusions. Preserve
            // authored cyclic correspondence even when vertex ordering differs
            // from the section's geometric bounding box.
            if same_profiles(&start, &end) {
                return Some(extrusion_node(model, frame, lo, hi, start));
            }
            Some(
                model.add_node(CadNode::Loft {
                    sections: [(lo, start), (hi, end)]
                        .into_iter()
                        .map(|(at, profile)| LoftSection {
                            profile,
                            transform: local_frame_transform(
                                add(frame.origin, mul(frame.axis, at)),
                                frame.x,
                                frame.y,
                                frame.axis,
                            ),
                        })
                        .collect(),
                }),
            )
        })
        .collect()
}

fn same_profiles(a: &Profile2d, b: &Profile2d) -> bool {
    a.loops.len() == b.loops.len()
        && a.loops.iter().zip(&b.loops).all(|(a, b)| {
            a.curves.len() == b.curves.len()
                && a.curves.iter().zip(&b.curves).all(|(a, b)| match (a, b) {
                    (
                        Curve2d::Line {
                            start_mm: a,
                            end_mm: ae,
                        },
                        Curve2d::Line {
                            start_mm: b,
                            end_mm: be,
                        },
                    ) => {
                        distance2(*a, *b) <= LENGTH
                            && distance2(*ae, *be) <= LENGTH
                            && distance2(*a, *ae) > LENGTH
                    }
                    (Curve2d::CircleArc { .. }, Curve2d::CircleArc { .. }) => a == b,
                    _ => false,
                })
        })
}

pub(super) fn model(
    solid: &ReducedSolid,
    frame: Frame,
) -> Result<Option<(CadModel, NodeId, usize)>> {
    let Some(sweep) = RuledSweep::new(solid, frame) else {
        return Ok(None);
    };
    let mut model = CadModel::new();
    let mut children = Vec::new();
    for ab in sweep.levels.windows(2) {
        let Some(parts) = slab_nodes(&mut model, &sweep, solid, frame, ab[0], ab[1]) else {
            return Ok(None);
        };
        if !parts.is_empty() {
            children.push(add_bounded_union_tree(&mut model, parts));
        }
    }
    if children.is_empty() {
        return Ok(None);
    }
    let root = add_bounded_union_tree(&mut model, children);
    model.add_root(root)?;
    model.validate()?;
    Ok(Some((model, root, sweep.levels.len() - 1)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collapsed_corners_preserve_cyclic_correspondence() -> Result<()> {
        let mut points = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [LENGTH * 0.2, 0.0]];
        anyhow::Context::context(
            canonicalize_collapsed_corners(&mut points),
            "valid cyclic collapse",
        )?;
        assert!(distance2(points[0], points[3]) <= f64::EPSILON);
        assert!(distance2(points[1], [2.0, 0.0]) <= f64::EPSILON);
        assert!(distance2(points[2], [2.0, 2.0]) <= f64::EPSILON);
        Ok(())
    }

    #[test]
    fn collapsed_corners_reject_transitive_displacement() {
        let mut points = [
            [0.0, 0.0],
            [0.75 * LENGTH, 0.0],
            [1.5 * LENGTH, 0.0],
            [2.0, 2.0],
        ];
        assert!(canonicalize_collapsed_corners(&mut points).is_none());
    }
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod native_tests {
    use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
    use crate::complete_ir::recover_complete_ir_with_kernel_bytes;
    use anyhow::Result;

    #[test]
    fn partitioned_spline_relief_retains_raw_vertex_incidence() -> Result<()> {
        let kernel = MonstertruckKernel::project_standard();
        let bytes =
            include_bytes!("../../../validation/fixtures/native_relief_raw_trim_52f96f03.step");
        let report = recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
        assert_eq!(report.body_count, 1);
        assert_eq!(report.constructively_recovered_solid_count, 1);
        assert_eq!(
            report.source_boundary_verified_solids,
            [report.body_scan.solids[0].source]
        );
        Ok(())
    }

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
