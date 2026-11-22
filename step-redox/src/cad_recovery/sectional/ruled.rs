//! Planar swept slabs with linearly changing sections. Within consecutive
//! vertex stations, each finite face intersects a slice in the same ordered
//! intervals. The mid-slice graph supplies correspondence; the actual
//! original crossing edges supply both end profiles, including collapsed
//! corner transitions. Independent full-source proof remains mandatory.
use super::{
    DIRECTION, Frame, LENGTH, SliceGraph, extrusion_node, section_ring_groups, valid_line,
};
use crate::{
    brep::{OrientedEdgeUse, SurfaceSupport},
    cad_ir::{CadModel, CadNode, Curve2d, LoftSection, NodeId, Profile2d, ProfileLoop},
    cad_recovery::{add_bounded_union_tree, solid::local_frame_transform},
    math2::distance as distance2,
    math3::{add, cross, dot, mul, norm, normalize, sub},
    profile_curves::RecoveredProfileCurve,
    reduced_brep::ReducedSolid,
};
use anyhow::Result;
use std::collections::BTreeMap;

struct PlanarSweep<'a> {
    edges: BTreeMap<u64, &'a OrientedEdgeUse>,
    levels: Vec<f64>,
}

impl<'a> PlanarSweep<'a> {
    fn new(solid: &'a ReducedSolid, frame: Frame) -> Option<Self> {
        let mut edges = BTreeMap::new();
        let mut levels = Vec::new();
        let mut tilted = false;
        for face in solid.faces() {
            let SurfaceSupport::Plane(plane) = face.surface else {
                return None;
            };
            let normal = normalize(plane.normal, DIRECTION)?;
            tilted |= norm(cross(normal, frame.axis)) > DIRECTION
                && dot(normal, frame.axis).abs() > DIRECTION;
            for edge in face.loop_edges() {
                if !valid_line(edge)
                    || [edge.start_mm, edge.end_mm]
                        .into_iter()
                        .any(|point| dot(sub(point, plane.origin_mm), normal).abs() > LENGTH)
                {
                    return None;
                }
                edges.entry(edge.edge_id).or_insert(edge);
                levels.extend([frame.axial(edge.start_mm), frame.axial(edge.end_mm)]);
            }
        }
        levels.sort_by(f64::total_cmp);
        levels.dedup_by(|a, b| (*a - *b).abs() <= LENGTH);
        (tilted && levels.len() >= 2 && levels.windows(2).all(|ab| ab[1] - ab[0] > 2. * LENGTH))
            .then_some(Self { edges, levels })
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
                let mut points = rings[index]
                    .iter()
                    .map(|curve| {
                        let RecoveredProfileCurve::Line {
                            source_edge_ids, ..
                        } = curve
                        else {
                            return None;
                        };
                        self.point(*source_edge_ids.first()?, frame, level)
                    })
                    .collect::<Option<Vec<_>>>()?;
                canonicalize_collapsed_corners(&mut points)?;
                Some(ProfileLoop {
                    curves: points
                        .iter()
                        .zip(points.iter().cycle().skip(1))
                        .map(|(&start_mm, &end_mm)| Curve2d::Line { start_mm, end_mm })
                        .collect(),
                })
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
    sweep: &PlanarSweep<'_>,
    solid: &ReducedSolid,
    frame: Frame,
    lo: f64,
    hi: f64,
) -> Option<Vec<NodeId>> {
    let mut graph = SliceGraph::default();
    for face in solid.faces() {
        graph.add_plane(face, frame, lo.midpoint(hi))?;
    }
    let rings = graph.rings()?;
    if rings.is_empty() {
        return Some(Vec::new());
    }
    let groups = section_ring_groups(&rings)?;
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
    a.loops.iter().zip(&b.loops).all(|(a, b)| {
        a.curves.iter().zip(&b.curves).all(|(a, b)| {
            let (
                Curve2d::Line {
                    start_mm: a,
                    end_mm: ae,
                },
                Curve2d::Line {
                    start_mm: b,
                    end_mm: be,
                },
            ) = (a, b)
            else {
                return false;
            };
            distance2(*a, *b) <= LENGTH
                && distance2(*ae, *be) <= LENGTH
                && distance2(*a, *ae) > LENGTH
        })
    })
}

pub(super) fn model(
    solid: &ReducedSolid,
    frame: Frame,
) -> Result<Option<(CadModel, NodeId, usize)>> {
    let Some(sweep) = PlanarSweep::new(solid, frame) else {
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
    fn collapsed_corners_preserve_cyclic_correspondence() {
        let mut points = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [LENGTH * 0.2, 0.0]];
        canonicalize_collapsed_corners(&mut points).unwrap();
        assert!(distance2(points[0], points[3]) <= f64::EPSILON);
        assert!(distance2(points[1], [2.0, 0.0]) <= f64::EPSILON);
        assert!(distance2(points[2], [2.0, 2.0]) <= f64::EPSILON);
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
