//! Reconstruct a tangent flat flange continuing a coaxial annular bend.
//! The actual cap outlines, tangent cylinders and taper-to-host line edges
//! fix the editable two-prism intersection. No source fixture dimensions
//! are part of this recognizer.

use super::{
    channel_round_ends::coaxial, solid::extrude_profile_in_frame, source_incidence::FaceIncidence,
};
use crate::{
    brep::{CurveSupport, CylinderSupport, SurfaceSupport},
    cad_ir::{BooleanOp, CadModel, CadNode, Curve2d, NodeId, ProofStatus, Provenance},
    math2,
    math3::{add, cross, dot, mul, norm, normalize, sub},
    profile_curves::RecoveredProfileCurve,
    reduced_brep::{ReducedFace, ReducedSolid},
    solid_extrusions::{CapInterval, RecoveredSolidExtrusion, SignedCapPairProposal},
};
use anyhow::Result;
use monstertruck_core::tolerance::{DEFAULT_STEP_RECOVERY_TOLERANCES, NATIVE_STEP_RECOVERY_MM};
use std::collections::BTreeSet;

const LENGTH: f64 = NATIVE_STEP_RECOVERY_MM;
const DIRECTION: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

#[derive(Debug, Clone, Copy)]
struct Frame {
    origin: [f64; 3],
    width: [f64; 3],
    side: [f64; 3],
    radial: [f64; 3],
}

impl Frame {
    fn local(self, point: [f64; 3]) -> [f64; 3] {
        let delta = sub(point, self.origin);
        [
            dot(delta, self.width),
            dot(delta, self.side),
            dot(delta, self.radial),
        ]
    }
}

#[derive(Debug, Clone, Copy)]
struct Band {
    inner: f64,
    outer: f64,
}

impl Band {
    fn radial_at(side: f64, radius: f64) -> Option<f64> {
        let gap = radius - side.abs();
        (gap > LENGTH && gap.is_finite()).then(|| (gap * (radius + side.abs())).sqrt())
    }
}

#[derive(Debug, Clone)]
pub(super) struct LipProposal {
    frame: Frame,
    band: Band,
    width_mm: f64,
    back: f64,
    front: f64,
    polygon: Vec<[f64; 2]>,
    source_face_ids: Vec<u64>,
    solid_id: u64,
}

impl LipProposal {
    pub(super) fn source_face_ids(&self) -> &[u64] { &self.source_face_ids }

    fn annular_profile(&self, back: f64, front: f64) -> Option<Vec<Curve2d>> {
        let inner = self.band.inner;
        let outer = self.band.outer;
        let inner_back = [back, Band::radial_at(back, inner)?];
        let inner_tangent = [0.0, inner];
        let inner_front = [front, inner];
        let outer_front = [front, outer];
        let outer_tangent = [0.0, outer];
        let outer_back = [back, Band::radial_at(back, outer)?];
        let line = |start_mm, end_mm| Curve2d::Line { start_mm, end_mm };
        let arc = |radius_mm, start: [f64; 2], end: [f64; 2]| Curve2d::CircleArc {
            center_mm: [0.0, 0.0],
            radius_mm,
            start_angle_rad: start[1].atan2(start[0]),
            end_angle_rad: end[1].atan2(end[0]),
        };
        Some(if front > 0.0 {
            vec![
                arc(inner, inner_back, inner_tangent),
                line(inner_tangent, inner_front),
                line(inner_front, outer_front),
                line(outer_front, outer_tangent),
                arc(outer, outer_tangent, outer_back),
                line(outer_back, inner_back),
            ]
        } else {
            vec![
                line(inner_back, outer_back),
                arc(outer, outer_back, outer_tangent),
                line(outer_tangent, outer_front),
                line(outer_front, inner_front),
                line(inner_front, inner_tangent),
                arc(inner, inner_tangent, inner_back),
            ]
        })
    }

    pub(super) fn model(&self) -> Result<(CadModel, NodeId)> {
        let sign = self.front.signum();
        // All extra boundaries are outside the source-witnessed taper
        // polygon. The annular root remains separated from tangency.
        let clearance = 0.25 * (self.band.outer - self.band.inner);
        let back_clearance = clearance.min(0.25 * (self.band.inner - self.back.abs()));
        let back = sign.mul_add(-back_clearance, self.back);
        let front = sign.mul_add(clearance, self.front);
        let mut model = CadModel::new();
        let profile = self
            .annular_profile(back, front)
            .ok_or_else(|| anyhow::anyhow!("annular lip profile reaches tangency"))?;
        let f = self.frame;
        let curved = extrude_profile_in_frame(
            &mut model,
            profile,
            2.0_f64.mul_add(clearance, self.width_mm),
            add(f.origin, mul(f.width, -clearance)),
            [f.side, f.radial, f.width],
        );
        let lower = Band::radial_at(back, self.band.inner)
            .ok_or_else(|| anyhow::anyhow!("annular lip lower station"))?
            - clearance;
        let upper = self.band.outer + clearance;
        let polygon = self
            .polygon
            .iter()
            .zip(self.polygon.iter().cycle().skip(1))
            .map(|(&start_mm, &end_mm)| Curve2d::Line { start_mm, end_mm })
            .collect();
        let taper = extrude_profile_in_frame(
            &mut model,
            polygon,
            upper - lower,
            add(f.origin, mul(f.radial, lower)),
            [f.width, f.side, f.radial],
        );
        let root = model.add_node(CadNode::Boolean {
            op: BooleanOp::Intersection,
            children: vec![curved, taper],
        });
        model.set_provenance(
            root,
            Provenance {
                source_entity_ids: std::iter::once(self.solid_id)
                    .chain(self.source_face_ids.iter().copied())
                    .collect(),
                proof: ProofStatus::HeuristicCandidate,
                max_residual_mm: None,
            },
        )?;
        model.add_root(root)?;
        model.validate()?;
        Ok((model, root))
    }
}

struct TangentFace {
    cylinder: CylinderSupport,
    face_id: u64,
    outward: bool,
    sign: f64,
}

fn tangent_face(
    graph: &FaceIncidence<'_>,
    cap: &ReducedFace,
    width: [f64; 3],
    radial: [f64; 3],
) -> Option<TangentFace> {
    let SurfaceSupport::Plane(plane) = cap.surface else {
        return None;
    };
    let mut matches = cap.loop_edges().filter_map(|edge| {
        let CurveSupport::Line(line) = edge.support else {
            return None;
        };
        if norm(cross(normalize(line.direction, DIRECTION)?, width)) > DIRECTION {
            return None;
        }
        let neighbor = graph.neighbor(cap.id, edge.edge_id)?;
        let SurfaceSupport::Cylinder(cylinder) = neighbor.surface else {
            return None;
        };
        if norm(cross(normalize(cylinder.axis, DIRECTION)?, width)) > DIRECTION {
            return None;
        }
        let signed = dot(sub(plane.origin_mm, cylinder.axis_origin_mm), radial);
        if (signed.abs() - cylinder.radius_mm).abs() > LENGTH || cylinder.radius_mm <= LENGTH {
            return None;
        }
        // The actual straight boundary must be the tangent generator.
        if [edge.start_mm, edge.end_mm].into_iter().any(|point| {
            let delta = sub(point, cylinder.axis_origin_mm);
            norm(sub(
                sub(delta, mul(width, dot(delta, width))),
                mul(radial, signed),
            )) > LENGTH
        }) {
            return None;
        }
        Some(TangentFace {
            cylinder,
            face_id: neighbor.id,
            outward: neighbor.surface_forward,
            sign: signed.signum(),
        })
    });
    let result = matches.next()?;
    // Independent segmentation can produce multiple coedges to the same
    // tangent support. A competing tangent carrier is ambiguous.
    matches
        .all(|m| m.face_id == result.face_id)
        .then_some(result)
}

fn side_faces<'a>(
    graph: &FaceIncidence<'a>,
    caps: [&ReducedFace; 2],
    radial: [f64; 3],
) -> Vec<&'a ReducedFace> {
    let sides = |cap: &ReducedFace| {
        cap.loop_edges()
            .filter_map(|edge| {
                let face = graph.neighbor(cap.id, edge.edge_id)?;
                let SurfaceSupport::Plane(plane) = face.surface else {
                    return None;
                };
                (dot(plane.normal, radial).abs() <= DIRECTION).then_some(face.id)
            })
            .collect::<BTreeSet<_>>()
    };
    sides(caps[0])
        .intersection(&sides(caps[1]))
        .filter_map(|&id| graph.face(id))
        .collect()
}

fn back_station(
    graph: &FaceIncidence<'_>,
    sides: &[&ReducedFace],
    host: &RecoveredSolidExtrusion,
    frame: Frame,
    band: Band,
) -> Option<f64> {
    let mut found = [None::<f64>; 2];
    for face in sides {
        for edge in face.loop_edges() {
            let Some(neighbor) = graph.neighbor(face.id, edge.edge_id) else {
                continue;
            };
            let Some(cap) = host.cap_face_ids.iter().position(|&id| id == neighbor.id) else {
                continue;
            };
            let CurveSupport::Line(line) = edge.support else {
                return None;
            };
            if norm(cross(normalize(line.direction, DIRECTION)?, frame.radial)) > DIRECTION {
                return None;
            }
            let a = frame.local(edge.start_mm);
            let b = frame.local(edge.end_mm);
            let width = if cap == 0 { 0.0 } else { host.height_mm };
            if (a[0] - width).abs() > LENGTH
                || (b[0] - width).abs() > LENGTH
                || (a[1] - b[1]).abs() > LENGTH
            {
                return None;
            }
            let station = a[1].midpoint(b[1]);
            let inner = Band::radial_at(station, band.inner)?;
            let outer = Band::radial_at(station, band.outer)?;
            if (a[2].min(b[2]) - inner).abs() > LENGTH || (a[2].max(b[2]) - outer).abs() > LENGTH {
                return None;
            }
            if found[cap].is_some_and(|previous| (previous - station).abs() > LENGTH) {
                return None;
            }
            found[cap] = Some(station);
        }
    }
    let [Some(a), Some(b)] = found else {
        return None;
    };
    ((a - b).abs() <= LENGTH).then_some(a.midpoint(b))
}

fn convex_polygon(points: &[[f64; 2]]) -> bool {
    points.len() >= 3
        && points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .all(|(&a, &b)| {
                let delta = math2::sub(b, a);
                let Some(direction) = math2::normalize(delta, LENGTH) else {
                    return false;
                };
                points
                    .iter()
                    .all(|&p| math2::cross(direction, math2::sub(p, a)) >= -LENGTH)
            })
}

fn projected_polygon(lip: &RecoveredSolidExtrusion, frame: Frame) -> Option<Vec<[f64; 2]>> {
    if !lip.inner_profile_loops.is_empty() {
        return None;
    }
    let mut points = lip
        .profile_curves
        .iter()
        .map(|curve| {
            let RecoveredProfileCurve::Line { start_mm, .. } = curve else {
                return None;
            };
            let world = add(
                add(lip.origin_mm, mul(lip.x_axis, start_mm[0])),
                mul(lip.y_axis, start_mm[1]),
            );
            let [x, y, _] = frame.local(world);
            (x.is_finite() && y.is_finite()).then_some([x, y])
        })
        .collect::<Option<Vec<_>>>()?;
    let area = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(&a, &b)| math2::cross(a, b))
        .sum::<f64>();
    if area < 0.0 {
        points.reverse();
    }
    convex_polygon(&points).then_some(points)
}

// Use the authored side-plane equations, not the tiny differences
// between separately evaluated cap vertices, to define the taper.
fn extend_rear(
    points: &mut [[f64; 2]],
    back: f64,
    width: f64,
    frame: Frame,
    sides: &[&ReducedFace],
) -> Option<f64> {
    let mut rear = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .enumerate()
        .filter(|(_, pair)| pair.0[1].abs() <= LENGTH && pair.1[1].abs() <= LENGTH)
        .map(|(i, _)| i);
    let index = rear.next()?;
    if rear.next().is_some() {
        return None;
    }
    let front = points
        .iter()
        .map(|p| p[1])
        .max_by(|a, b| a.abs().total_cmp(&b.abs()))?;
    if front.abs() <= LENGTH
        || front * back >= 0.0
        || points.iter().any(|p| p[1] * front.signum() < -LENGTH)
    {
        return None;
    }
    let equations = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .enumerate()
        .map(|(i, (&a, &b))| {
            if i == index {
                return Some(([0.0, 1.0], back));
            }
            let mut candidates = sides.iter().filter_map(|face| {
                let SurfaceSupport::Plane(plane) = face.surface else {
                    return None;
                };
                let normal = [
                    dot(plane.normal, frame.width),
                    dot(plane.normal, frame.side),
                ];
                let offset = dot(sub(plane.origin_mm, frame.origin), plane.normal);
                let evaluate = |p: [f64; 2]| normal[0].mul_add(p[0], normal[1] * p[1]) - offset;
                (evaluate(a).abs() <= LENGTH && evaluate(b).abs() <= LENGTH)
                    .then_some((normal, offset))
            });
            let equation = candidates.next()?;
            candidates.next().is_none().then_some(equation)
        })
        .collect::<Option<Vec<_>>>()?;
    let intersection = |(n, a): ([f64; 2], f64), (m, b): ([f64; 2], f64)| {
        let det = math2::cross(n, m);
        (det.abs() > DIRECTION).then(|| {
            [
                a.mul_add(m[1], -b * n[1]) / det,
                n[0].mul_add(b, -m[0] * a) / det,
            ]
        })
    };
    let count = points.len();
    let next = (index + 1) % count;
    for i in 0..count {
        let p = intersection(equations[(i + count - 1) % count], equations[i])?;
        if !p.iter().all(|v| v.is_finite())
            || (i != index && i != next && math2::norm(math2::sub(p, points[i])) > LENGTH)
        {
            return None;
        }
        points[i] = p;
    }
    let a = points[index];
    let b = points[next];
    if a[0].min(b[0]).abs() > LENGTH || (a[0].max(b[0]) - width).abs() > LENGTH {
        return None;
    }
    convex_polygon(points).then_some(front)
}

fn proposal(
    graph: &FaceIncidence<'_>,
    host: &RecoveredSolidExtrusion,
    lip: &RecoveredSolidExtrusion,
) -> Option<LipProposal> {
    let width = normalize(host.z_axis, DIRECTION)?;
    let axis = normalize(lip.z_axis, DIRECTION)?;
    if dot(width, axis).abs() > DIRECTION {
        return None;
    }
    let caps = [
        graph.face(lip.cap_face_ids[0])?,
        graph.face(lip.cap_face_ids[1])?,
    ];
    let first = tangent_face(graph, caps[0], width, axis)?;
    let second = tangent_face(graph, caps[1], width, axis)?;
    if first.outward == second.outward
        || (first.sign - second.sign).abs() > DIRECTION
        || !coaxial(first.cylinder, second.cylinder)
    {
        return None;
    }
    let (inner, outer) = if first.outward {
        (second.cylinder, first.cylinder)
    } else {
        (first.cylinder, second.cylinder)
    };
    let band = Band {
        inner: inner.radius_mm,
        outer: outer.radius_mm,
    };
    if band.outer - band.inner <= LENGTH || (band.outer - band.inner - lip.height_mm).abs() > LENGTH
    {
        return None;
    }
    let radial = mul(axis, first.sign);
    let origin = add(
        inner.axis_origin_mm,
        mul(width, dot(sub(host.origin_mm, inner.axis_origin_mm), width)),
    );
    let frame = Frame {
        origin,
        width,
        radial,
        side: normalize(cross(radial, width), DIRECTION)?,
    };
    let sides = side_faces(graph, caps, radial);
    let back = back_station(graph, &sides, host, frame, band)?;
    let mut polygon = projected_polygon(lip, frame)?;
    let front = extend_rear(&mut polygon, back, host.height_mm, frame, &sides)?;
    let source_face_ids = host
        .cap_face_ids
        .into_iter()
        .chain(lip.cap_face_ids)
        .chain([first.face_id, second.face_id])
        .chain(sides.iter().map(|f| f.id))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    Some(LipProposal {
        frame,
        band,
        width_mm: host.height_mm,
        back,
        front,
        polygon,
        source_face_ids,
        solid_id: host.solid_id,
    })
}

pub(super) fn detect(
    solid: &ReducedSolid,
    host: &RecoveredSolidExtrusion,
    proposals: &[SignedCapPairProposal],
) -> Vec<LipProposal> {
    let graph = FaceIncidence::new(solid);
    proposals
        .iter()
        .filter(|p| p.interval == CapInterval::Material)
        .filter_map(|p| proposal(&graph, host, &p.prism))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{solid_extrusions::detect_signed_cap_pairs, step_graph::build_index};
    use anyhow::{Context, ensure};

    fn mixed_lip() -> Result<LipProposal> {
        let exchange = step_p21::parser::parse(include_str!(
            "../../validation/fixtures/native_mixed_treatment_4a6356.step"
        ))?;
        let entities = &exchange.data[0].entities;
        let proposals = detect_signed_cap_pairs(entities);
        let host = proposals
            .iter()
            .find(|p| p.interval == CapInterval::Material && p.prism.cap_face_ids == [620, 1159])
            .context("host cap pair")?;
        let solid = ReducedSolid::from_step(host.prism.solid_id, entities, &build_index(entities))
            .context("source topology")?;
        let mut lips = detect(&solid, &host.prism, &proposals);
        ensure!(
            lips.len() == 1,
            "expected one source lip, found {}",
            lips.len()
        );
        lips.pop().context("source lip")
    }

    #[test]
    fn mixed_tapered_lip_recovers_from_tangent_caps_and_actual_back_edges() -> Result<()> {
        let lip = mixed_lip()?;
        assert!(lip.source_face_ids().contains(&175) && lip.source_face_ids().contains(&299));
        assert!((lip.band.inner - 0.3429).abs() < LENGTH);
        assert!((lip.band.outer - 0.6731).abs() < LENGTH);
        assert!((lip.back.abs() - 0.0127).abs() < LENGTH);
        assert!((lip.front.abs() - 0.4953).abs() < LENGTH);
        let (model, root) = lip.model()?;
        assert_eq!(model.roots, [root]);
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn source_tapered_lip_replays_reference_boundary() -> Result<()> {
        use crate::{
            boundary_proof::compare_analytic_boundaries,
            cad_kernel::{CadKernel, monstertruck::MonstertruckKernel},
            cad_recovery::{evaluated_step_solid, mixed_fixture_models},
        };
        let kernel = MonstertruckKernel::standard_mm(1e-5);
        let (model, root) = mixed_lip()?.model()?;
        let (reference, reference_root) = mixed_fixture_models::curved_lip()?;
        let actual = kernel.evaluate(&model, root)?;
        let expected = kernel.evaluate(&reference, reference_root)?;
        let (a, ae) =
            evaluated_step_solid(&kernel, &actual)?.context("source-derived lip solid")?;
        let (b, be) = evaluated_step_solid(&kernel, &expected)?.context("reference lip solid")?;
        let comparison = compare_analytic_boundaries(a, &ae, b, &be);
        ensure!(
            comparison.equivalent,
            "tapered lip differs: {comparison:#?}"
        );
        Ok(())
    }
}
