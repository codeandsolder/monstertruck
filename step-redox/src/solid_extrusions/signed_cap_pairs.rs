//! Paired planar caps of material prisms and empty channels.
//!
//! The two cap profiles are compared in one physical frame after directed
//! line/arc chain simplification, rather than matching the source edge
//! counts or topology IDs. Face normals distinguish an occupied interval
//! hypothesis from a cutter hypothesis. Neither implies complete side-wall
//! coverage or source-solid equivalence.

use super::{
    CanonicalProfile, DIR_TOL, FaceInfo, GEOM_TOL_MM, RecoveredProfileCurve,
    RecoveredSolidExtrusion, canonical_direction, canonical_profile, canonical_profile_loop,
    cap_loop_roles, distance2, dot, face_info, face_lies_on_plane, mul, parallel,
};
use crate::{
    brep::{self, SurfaceSupport},
    step_graph::build_index,
};
use serde::Serialize;
use std::collections::HashMap;
use step_p21::ast::EntityInstance;

/// Orientation of two translated source caps relative to their interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CapInterval {
    Material,
    Void,
}

/// A paired-cap prism hypothesis, not a recovered whole source solid.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SignedCapPairProposal {
    pub interval: CapInterval,
    pub prism: RecoveredSolidExtrusion,
}

fn control_residual(first: &[[f64; 2]], second: &[[f64; 2]]) -> Option<f64> {
    if first.len() != second.len() || first.is_empty() {
        return None;
    }
    first
        .iter()
        .zip(second)
        .try_fold(0.0_f64, |maximum, (&a, &b)| {
            let residual = distance2(a, b);
            residual.is_finite().then_some(maximum.max(residual))
        })
}

fn same_parameters(first: &[f64], second: &[f64]) -> bool {
    first.len() == second.len()
        && first
            .iter()
            .zip(second)
            .all(|(a, b)| a.is_finite() && b.is_finite() && a.total_cmp(b).is_eq())
}

fn curve_residual(first: &RecoveredProfileCurve, second: &RecoveredProfileCurve) -> Option<f64> {
    use RecoveredProfileCurve::{BSpline, Bezier, CircleArc, Line};
    let endpoints = distance2(first.start_point()?, second.start_point()?)
        .max(distance2(first.end_point()?, second.end_point()?));
    let support = match (first, second) {
        (Line { .. }, Line { .. }) => 0.0,
        (
            CircleArc {
                center_mm: a,
                radius_mm: ar,
                start_angle_rad: as_,
                end_angle_rad: ae,
                ..
            },
            CircleArc {
                center_mm: b,
                radius_mm: br,
                start_angle_rad: bs,
                end_angle_rad: be,
                ..
            },
        ) => {
            let sweep = (ae - as_) - (be - bs);
            distance2(*a, *b) + (ar - br).abs() + sweep.abs() * ar.max(*br)
        }
        (
            Bezier {
                control_points_mm: a,
                ..
            },
            Bezier {
                control_points_mm: b,
                ..
            },
        ) => control_residual(a, b)?,
        (
            BSpline {
                degree: ad,
                control_points_mm: a,
                knots: ak,
                weights: aw,
                ..
            },
            BSpline {
                degree: bd,
                control_points_mm: b,
                knots: bk,
                weights: bw,
                ..
            },
        ) if ad == bd
            && same_parameters(ak, bk)
            && match (aw, bw) {
                (None, None) => true,
                (Some(a), Some(b)) => same_parameters(a, b) && a.iter().all(|&w| w > 0.0),
                _ => false,
            } =>
        {
            control_residual(a, b)?
        }
        _ => return None,
    };
    let residual = endpoints.max(support);
    (residual.is_finite() && residual <= GEOM_TOL_MM).then_some(residual)
}

fn directed_loop_residual(
    first: &[RecoveredProfileCurve],
    second: &[RecoveredProfileCurve],
) -> Option<f64> {
    if first.len() != second.len() || first.is_empty() {
        return None;
    }
    (0..second.len())
        .filter_map(|shift| {
            first
                .iter()
                .zip(second.iter().cycle().skip(shift))
                .try_fold(0.0_f64, |max, (a, b)| Some(max.max(curve_residual(a, b)?)))
        })
        .min_by(f64::total_cmp)
}

fn loop_residual(first: &[RecoveredProfileCurve], second: &[RecoveredProfileCurve]) -> Option<f64> {
    directed_loop_residual(first, second).or_else(|| {
        let reversed = second
            .iter()
            .rev()
            .map(RecoveredProfileCurve::reversed)
            .collect::<Vec<_>>();
        directed_loop_residual(first, &reversed)
    })
}

fn profile_residual(profile: &CanonicalProfile, top: &FaceInfo, z_axis: [f64; 3]) -> Option<f64> {
    let roles = cap_loop_roles(top, z_axis)?;
    let mut top_loops = top
        .loops
        .iter()
        .zip(roles)
        .map(|(wire, outer)| {
            let curves = canonical_profile_loop(
                &wire.edges,
                profile.origin_mm,
                profile.x_axis,
                profile.y_axis,
                z_axis,
                outer,
            )?;
            Some((outer, curves))
        })
        .collect::<Option<Vec<_>>>()?;
    if top_loops.len() != profile.inner_loops.len() + 1 {
        return None;
    }
    let mut maximum = 0.0_f64;
    for (outer, curves) in std::iter::once((true, profile.curves.as_slice())).chain(
        profile
            .inner_loops
            .iter()
            .map(|ring| (false, ring.as_slice())),
    ) {
        let mut matches = top_loops
            .iter()
            .enumerate()
            .filter_map(|(i, (role, other))| {
                (*role == outer)
                    .then(|| loop_residual(curves, other))
                    .flatten()
                    .map(|r| (i, r))
            });
        let (index, residual) = matches.next()?;
        if matches.next().is_some() {
            return None;
        }
        top_loops.swap_remove(index);
        maximum = maximum.max(residual);
    }
    Some(maximum)
}

fn pair(solid_id: u64, first: &FaceInfo, second: &FaceInfo) -> Option<SignedCapPairProposal> {
    let (SurfaceSupport::Plane(a), SurfaceSupport::Plane(b)) = (&first.surface, &second.surface)
    else {
        return None;
    };
    if !parallel(a.normal, b.normal)
        || first.loops.len() != second.loops.len()
        || !face_lies_on_plane(first, *a)
        || !face_lies_on_plane(second, *b)
    {
        return None;
    }
    let axis = canonical_direction(a.normal);
    let distance = dot(axis, b.origin_mm) - dot(axis, a.origin_mm);
    if !distance.is_finite() || distance.abs() <= GEOM_TOL_MM {
        return None;
    }
    let (bottom, top, bottom_plane, top_plane) = if distance > 0.0 {
        (first, second, a, b)
    } else {
        (second, first, b, a)
    };
    let orientation = |face: &FaceInfo, normal: [f64; 3]| {
        dot(
            axis,
            mul(normal, if face.surface_forward { 1.0 } else { -1.0 }),
        )
    };
    let (bottom_sign, top_sign) = (
        orientation(bottom, bottom_plane.normal),
        orientation(top, top_plane.normal),
    );
    let interval = if bottom_sign <= -1.0 + DIR_TOL && top_sign >= 1.0 - DIR_TOL {
        CapInterval::Material
    } else if bottom_sign >= 1.0 - DIR_TOL && top_sign <= -1.0 + DIR_TOL {
        CapInterval::Void
    } else {
        return None;
    };
    let profile = canonical_profile(bottom, axis, &cap_loop_roles(bottom, axis)?)?;
    let residual = profile_residual(&profile, top, axis)?
        .max(a.max_residual_mm)
        .max(b.max_residual_mm);
    Some(SignedCapPairProposal {
        interval,
        prism: RecoveredSolidExtrusion {
            solid_id,
            cap_face_ids: [bottom.id, top.id],
            side_face_ids: Vec::new(),
            profile_curves: profile.curves,
            inner_profile_loops: profile.inner_loops,
            origin_mm: profile.origin_mm,
            x_axis: profile.x_axis,
            y_axis: profile.y_axis,
            z_axis: axis,
            height_mm: distance.abs(),
            max_residual_mm: residual,
        },
    })
}

pub fn detect_signed_cap_pairs_for_solid_with_index(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Vec<SignedCapPairProposal> {
    let Some(ids) = brep::solid_face_ids(solid_id, entities, index) else {
        return Vec::new();
    };
    let faces = ids
        .into_iter()
        .filter_map(|id| face_info(id, entities, index))
        .filter(|f| matches!(f.surface, SurfaceSupport::Plane(_)))
        .collect::<Vec<_>>();
    let mut out = Vec::new();
    for (i, first) in faces.iter().enumerate() {
        for second in &faces[i + 1..] {
            if let Some(proposal) = pair(solid_id, first, second) {
                out.push(proposal);
            }
        }
    }
    out.sort_by_key(|p| p.prism.cap_face_ids);
    out
}

/// Find translated planar cap-pair proposals, including outward-facing
/// channels whose caps have independently partitioned analytic boundaries.
#[must_use = "inspect signed cap proposals before constructive replay"]
pub fn detect_signed_cap_pairs(entities: &[EntityInstance]) -> Vec<SignedCapPairProposal> {
    let index = build_index(entities);
    brep::manifold_solid_ids(entities)
        .into_iter()
        .flat_map(|solid_id| {
            detect_signed_cap_pairs_for_solid_with_index(solid_id, entities, &index)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Context, Result};

    #[test]
    fn mixed_source_void_caps_are_recognized_despite_arc_partition() -> Result<()> {
        let bytes = include_str!("../../validation/fixtures/native_mixed_treatment_4a6356.step");
        let source = step_p21::parser::parse(bytes)?;
        let proposals = detect_signed_cap_pairs(&source.data[0].entities);
        let void = proposals
            .iter()
            .find(|p| p.prism.cap_face_ids == [391, 429])
            .context("translated source channel not recognized")?;
        assert_eq!(void.interval, CapInterval::Void);
        assert!((void.prism.height_mm - 0.3048).abs() < GEOM_TOL_MM);
        let material = proposals
            .iter()
            .find(|p| p.prism.cap_face_ids == [620, 1159])
            .context("translated source host caps not recognized")?;
        assert_eq!(material.interval, CapInterval::Material);
        assert!((material.prism.height_mm - 1.27).abs() < GEOM_TOL_MM);
        Ok(())
    }
}
