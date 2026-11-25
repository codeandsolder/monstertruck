//! Source-guided semicylindrical terminations of a paired-cap empty channel.
//!
//! The channel spacing supplies the diameter. Shared source edges supply
//! the two generators and the actual end carriers. This produces editable
//! cutter hypotheses only; admitting any composed body still requires the
//! independent oriented whole-source boundary proof.

use super::source_incidence::FaceIncidence;
use crate::{
    brep::{CurveSupport, CylinderSupport, OrientedEdgeUse, SurfaceSupport},
    math3::{canonical_unit_direction, cross, distance, dot, norm, normalize, sub},
    reduced_brep::{ReducedFace, ReducedSolid},
    solid_extrusions::{RecoveredSolidExtrusion, source_boundary_edge_witness_points},
};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::BTreeSet;

mod model;

const LENGTH: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIRECTION: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

#[derive(Debug, Clone, Copy)]
struct RoundFrame {
    origin: [f64; 3],
    width_axis: [f64; 3],
    side_axis: [f64; 3],
    axis: [f64; 3],
    radius: f64,
    positive_side: bool,
}

#[derive(Debug, Clone, Copy)]
enum EndCarriers {
    Planar {
        lower: f64,
        upper: f64,
    },
    Coaxial {
        center_side: f64,
        center_axial: f64,
        inner_radius: f64,
        outer_radius: f64,
        positive_root: bool,
    },
}

#[derive(Debug, Clone)]
pub(super) struct RoundEndProposal {
    frame: RoundFrame,
    ends: EndCarriers,
    source_entity_ids: Vec<u64>,
    round_face_ids: Vec<u64>,
}

impl RoundEndProposal {
    pub(super) fn round_face_ids(&self) -> &[u64] { &self.round_face_ids }
    pub(super) fn source_entity_ids(&self) -> &[u64] { &self.source_entity_ids }
}

impl RoundFrame {
    fn from_channel(cylinder: CylinderSupport, channel: &RecoveredSolidExtrusion) -> Option<Self> {
        let axis = canonical_unit_direction(cylinder.axis, DIRECTION)?;
        let width_axis = normalize(channel.z_axis, DIRECTION)?;
        let radius = cylinder.radius_mm;
        if !radius.is_finite()
            || radius <= LENGTH
            || dot(axis, width_axis).abs() > DIRECTION
            || (2.0 * radius - channel.height_mm).abs() > LENGTH
            || (dot(sub(cylinder.axis_origin_mm, channel.origin_mm), width_axis) - radius).abs()
                > LENGTH
        {
            return None;
        }
        Some(Self {
            origin: cylinder.axis_origin_mm,
            width_axis,
            side_axis: normalize(cross(axis, width_axis), DIRECTION)?,
            axis,
            radius,
            positive_side: false,
        })
    }

    fn local(self, point: [f64; 3]) -> [f64; 3] {
        let delta = sub(point, self.origin);
        [
            dot(delta, self.width_axis),
            dot(delta, self.side_axis),
            dot(delta, self.axis),
        ]
    }

    fn diametral_generator(self, edge: &OrientedEdgeUse, positive_width: bool) -> bool {
        let CurveSupport::Line(line) = edge.support else {
            return false;
        };
        let Some(direction) = normalize(line.direction, DIRECTION) else {
            return false;
        };
        if norm(cross(direction, self.axis)) > DIRECTION
            || distance(edge.start_mm, edge.end_mm) <= LENGTH
        {
            return false;
        }
        let expected = if positive_width {
            self.radius
        } else {
            -self.radius
        };
        [edge.start_mm, edge.end_mm].into_iter().all(|point| {
            let [width, side, _] = self.local(point);
            (width - expected).abs() <= LENGTH && side.abs() <= LENGTH
        })
    }
}

pub(super) fn coaxial(first: CylinderSupport, second: CylinderSupport) -> bool {
    let (Some(a), Some(b)) = (
        normalize(first.axis, DIRECTION),
        normalize(second.axis, DIRECTION),
    ) else {
        return false;
    };
    norm(cross(a, b)) <= DIRECTION
        && norm(cross(sub(second.axis_origin_mm, first.axis_origin_mm), a)) <= LENGTH
}

fn same_cylinder(first: CylinderSupport, second: CylinderSupport) -> bool {
    coaxial(first, second) && (first.radius_mm - second.radius_mm).abs() <= LENGTH
}

fn group<'a>(
    incidence: &FaceIncidence<'a>,
    first: &'a ReducedFace,
    cylinder: CylinderSupport,
) -> Option<Vec<&'a ReducedFace>> {
    let mut seen = BTreeSet::from([first.id]);
    let mut faces = vec![first];
    let mut index = 0;
    while let Some(face) = faces.get(index).copied() {
        for edge in face.loop_edges() {
            // An authored periodic seam has two uses on this same face.
            // It is not an external neighbor.
            if incidence
                .owners(edge.edge_id)?
                .iter()
                .all(|f| f.id == face.id)
            {
                continue;
            }
            let neighbor = incidence.neighbor(face.id, edge.edge_id)?;
            if !neighbor.surface_forward
                && let SurfaceSupport::Cylinder(other) = neighbor.surface
                && same_cylinder(cylinder, other)
                && seen.insert(neighbor.id)
            {
                faces.push(neighbor);
            }
        }
        index += 1;
    }
    Some(faces)
}

struct EndBoundary<'a> {
    face: &'a ReducedFace,
    edge: &'a OrientedEdgeUse,
}

fn external_boundaries<'a>(
    incidence: &FaceIncidence<'a>,
    group: &[&'a ReducedFace],
    channel: &RecoveredSolidExtrusion,
    frame: RoundFrame,
) -> Option<Vec<EndBoundary<'a>>> {
    let mut cap_seen = [false; 2];
    let mut ends = Vec::new();
    let ids = group.iter().map(|f| f.id).collect::<BTreeSet<_>>();
    for face in group {
        for edge in face.loop_edges() {
            if incidence
                .owners(edge.edge_id)?
                .iter()
                .all(|f| ids.contains(&f.id))
            {
                continue;
            }
            let other = incidence.neighbor(face.id, edge.edge_id)?;
            if let Some(cap) = channel.cap_face_ids.iter().position(|&id| id == other.id) {
                if !frame.diametral_generator(edge, cap == 1) {
                    return None;
                }
                cap_seen[cap] = true;
            } else {
                ends.push(EndBoundary { face: other, edge });
            }
        }
    }
    (cap_seen.into_iter().all(|seen| seen) && !ends.is_empty()).then_some(ends)
}

fn side(frame: RoundFrame, boundaries: &[EndBoundary<'_>]) -> Option<bool> {
    let mut lo = 0.0_f64;
    let mut hi = 0.0_f64;
    for boundary in boundaries {
        let edge = boundary.edge;
        for point in [edge.start_mm, edge.end_mm]
            .into_iter()
            .chain(source_boundary_edge_witness_points(edge))
        {
            let [_, y, _] = frame.local(point);
            if !y.is_finite() {
                return None;
            }
            lo = lo.min(y);
            hi = hi.max(y);
        }
    }
    match (lo < -LENGTH, hi > LENGTH) {
        (true, false) => Some(false),
        (false, true) => Some(true),
        _ => None,
    }
}

fn planar_ends(frame: RoundFrame, boundaries: &[EndBoundary<'_>]) -> Option<EndCarriers> {
    let mut stations = Vec::<f64>::new();
    for boundary in boundaries {
        let SurfaceSupport::Plane(plane) = boundary.face.surface else {
            return None;
        };
        let CurveSupport::Circle(circle) = boundary.edge.support else {
            return None;
        };
        if norm(cross(plane.normal, frame.axis)) > DIRECTION
            || norm(cross(circle.normal, frame.axis)) > DIRECTION
            || (circle.radius_mm - frame.radius).abs() > LENGTH
        {
            return None;
        }
        let station = frame.local(plane.origin_mm)[2];
        let [width, side, center_station] = frame.local(circle.center_mm);
        if width.hypot(side) > LENGTH
            || (center_station - station).abs() > LENGTH
            || [boundary.edge.start_mm, boundary.edge.end_mm]
                .into_iter()
                .any(|point| (frame.local(point)[2] - station).abs() > LENGTH)
        {
            return None;
        }
        if !stations.iter().any(|&s| (s - station).abs() <= LENGTH) {
            stations.push(station);
        }
    }
    stations.sort_by(f64::total_cmp);
    let [lower, upper] = stations.as_slice() else {
        return None;
    };
    (*upper - *lower > LENGTH).then_some(EndCarriers::Planar {
        lower: *lower,
        upper: *upper,
    })
}

fn coaxial_ends(frame: RoundFrame, boundaries: &[EndBoundary<'_>]) -> Option<EndCarriers> {
    let mut inner = None::<CylinderSupport>;
    let mut outer = None::<CylinderSupport>;
    let mut root = None;
    for boundary in boundaries {
        let SurfaceSupport::Cylinder(cylinder) = boundary.face.surface else {
            return None;
        };
        if norm(cross(cylinder.axis, frame.width_axis)) > DIRECTION {
            return None;
        }
        let slot = if boundary.face.surface_forward {
            &mut outer
        } else {
            &mut inner
        };
        if slot.is_some_and(|prior| !same_cylinder(prior, cylinder)) {
            return None;
        }
        *slot = Some(cylinder);
        let center = frame.local(cylinder.axis_origin_mm)[2];
        for point in [boundary.edge.start_mm, boundary.edge.end_mm]
            .into_iter()
            .chain(source_boundary_edge_witness_points(boundary.edge))
        {
            let delta = frame.local(point)[2] - center;
            if !delta.is_finite() || delta.abs() <= LENGTH {
                return None;
            }
            let positive = delta > 0.0;
            if root.is_some_and(|old| old != positive) {
                return None;
            }
            root = Some(positive);
        }
    }
    let (inner, outer) = (inner?, outer?);
    if !coaxial(inner, outer) || outer.radius_mm - inner.radius_mm <= LENGTH {
        return None;
    }
    let [_, center_side, center_axial] = frame.local(inner.axis_origin_mm);
    // The complete semicircular footprint must fit a separated real root
    // of each limiting cylindrical carrier. Tangent or branch-switching
    // cuts require a different constructive grammar.
    let tip = if frame.positive_side {
        frame.radius
    } else {
        -frame.radius
    };
    if center_side.abs().max((tip - center_side).abs()) >= inner.radius_mm - LENGTH {
        return None;
    }
    Some(EndCarriers::Coaxial {
        center_side,
        center_axial,
        inner_radius: inner.radius_mm,
        outer_radius: outer.radius_mm,
        positive_root: root?,
    })
}

fn proposal(
    solid_id: u64,
    incidence: &FaceIncidence<'_>,
    group: &[&ReducedFace],
    channel: &RecoveredSolidExtrusion,
    mut frame: RoundFrame,
) -> Option<RoundEndProposal> {
    let boundaries = external_boundaries(incidence, group, channel, frame)?;
    frame.positive_side = side(frame, &boundaries)?;
    let ends = planar_ends(frame, &boundaries).or_else(|| coaxial_ends(frame, &boundaries))?;
    let source_entity_ids = std::iter::once(solid_id)
        .chain(channel.cap_face_ids)
        .chain(group.iter().map(|f| f.id))
        .chain(boundaries.iter().map(|b| b.face.id))
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    Some(RoundEndProposal {
        frame,
        ends,
        source_entity_ids,
        round_face_ids: group.iter().map(|f| f.id).collect(),
    })
}

pub(super) fn detect(
    solid: &ReducedSolid,
    channel: &RecoveredSolidExtrusion,
) -> Vec<RoundEndProposal> {
    let incidence = FaceIncidence::new(solid);
    let mut seen = BTreeSet::new();
    let mut proposals = Vec::new();
    for face in solid.faces() {
        let SurfaceSupport::Cylinder(cylinder) = face.surface else {
            continue;
        };
        if face.surface_forward || seen.contains(&face.id) {
            continue;
        }
        let Some(frame) = RoundFrame::from_channel(cylinder, channel) else {
            continue;
        };
        let Some(group) = group(&incidence, face, cylinder) else {
            continue;
        };
        seen.extend(group.iter().map(|f| f.id));
        if let Some(proposal) = proposal(solid.solid_id(), &incidence, &group, channel, frame) {
            proposals.push(proposal);
        }
    }
    proposals
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        solid_extrusions::{CapInterval, detect_signed_cap_pairs},
        step_graph::build_index,
    };
    use anyhow::{Context, Result};

    fn mixed_rounds() -> Result<Vec<RoundEndProposal>> {
        let source = step_p21::parser::parse(include_str!(
            "../../validation/fixtures/native_mixed_treatment_4a6356.step"
        ))?;
        let entities = &source.data[0].entities;
        let channel = detect_signed_cap_pairs(entities)
            .into_iter()
            .find(|p| p.interval == CapInterval::Void && p.prism.cap_face_ids == [391, 429])
            .context("source channel absent")?;
        let solid =
            ReducedSolid::from_step(channel.prism.solid_id, entities, &build_index(entities))
                .context("source topology absent")?;
        Ok(detect(&solid, &channel.prism))
    }

    #[test]
    fn mixed_source_round_ends_are_derived_from_channel_and_end_carriers() -> Result<()> {
        let rounds = mixed_rounds()?;
        assert_eq!(rounds.len(), 2);
        let planar = rounds
            .iter()
            .find(|p| matches!(p.ends, EndCarriers::Planar { .. }))
            .context("partitioned lower round absent")?;
        assert!(planar.source_entity_ids.contains(&68) && planar.source_entity_ids.contains(&1079));
        let EndCarriers::Planar { lower, upper } = planar.ends else {
            anyhow::bail!("planar ends");
        };
        assert!((upper - lower - 0.3302).abs() < LENGTH);
        let curved = rounds
            .iter()
            .find(|p| matches!(p.ends, EndCarriers::Coaxial { .. }))
            .context("crossing-cylinder upper round absent")?;
        assert!(curved.source_entity_ids.contains(&680));
        let EndCarriers::Coaxial {
            inner_radius,
            outer_radius,
            ..
        } = curved.ends
        else {
            anyhow::bail!("coaxial ends");
        };
        assert!((inner_radius - 0.3429).abs() < LENGTH);
        assert!((outer_radius - 0.6731).abs() < LENGTH);
        for round in rounds {
            assert!((round.frame.radius - 0.1524).abs() < LENGTH);
            let (model, root) = round.model()?;
            model.validate()?;
            assert_eq!(model.roots, [root]);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn source_round_end_models_reproduce_bounded_cutter_geometry() -> Result<()> {
        use crate::{
            boundary_proof::compare_analytic_boundaries,
            cad_kernel::{CadKernel, monstertruck::MonstertruckKernel},
            cad_recovery::{evaluated_step_solid, mixed_fixture_models},
        };
        let kernel = MonstertruckKernel::standard_mm(1e-5);
        let (reference, upper, lower) = mixed_fixture_models::rounding_caps()?;
        for proposal in mixed_rounds()? {
            let reference_root = match proposal.ends {
                EndCarriers::Planar { .. } => lower,
                EndCarriers::Coaxial { .. } => upper,
            };
            let (model, root) = proposal.model()?;
            let expected = kernel.evaluate(&reference, reference_root)?;
            let actual = kernel.evaluate(&model, root)?;
            assert!(kernel.summarize(&actual).is_single_solid());
            let (first, first_entities) =
                evaluated_step_solid(&kernel, &expected)?.context("reference solid")?;
            let (second, second_entities) =
                evaluated_step_solid(&kernel, &actual)?.context("recovered cutter solid")?;
            let proof =
                compare_analytic_boundaries(first, &first_entities, second, &second_entities);
            anyhow::ensure!(proof.equivalent, "source cutter proof failed: {proof:#?}");
        }
        Ok(())
    }
}
