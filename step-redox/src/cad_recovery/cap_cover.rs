//! Find a constructive union of material-only source cap extrusions.
//!
//! Finite boundary witnesses are used solely to select proposals, never to
//! establish equivalence. Only a native whole-source boundary proof admits a
//! recovered solid.
use super::{
    CadFragment, CadFragmentSource, SourceRecovery,
    cap_assembly::{AssemblyFrame, proven_cap_assembly},
};
use crate::cad_kernel::CadKernel;
use crate::math3::{dot, sub};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use crate::solid_extrusions::{
    RecoveredSourceCapExtrusion, cap_contains, partial_source_cap_extrusion_candidates_with_index,
    profile_contains_with_rotated_ray, signed_area_curves, source_boundary_edge_witness_points,
};
use anyhow::Result;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::HashSet;

const LENGTH: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;

fn strict_contains(cap: &RecoveredSourceCapExtrusion, point: [f64; 3]) -> bool {
    let delta = sub(point, cap.origin_mm);
    let depth = dot(delta, cap.z_axis);
    let margin = 4.0 * LENGTH;
    if depth <= margin || depth >= cap.height_mm - margin {
        return false;
    }
    let point = [dot(delta, cap.x_axis), dot(delta, cap.y_axis)];
    profile_contains_with_rotated_ray(&cap.profile_curves, point, margin) == Some(true)
        && cap
            .inner_profile_loops
            .iter()
            .all(|ring| profile_contains_with_rotated_ray(ring, point, margin) == Some(false))
}

fn volume(cap: &RecoveredSourceCapExtrusion) -> Option<f64> {
    if cap
        .profile_curves
        .iter()
        .chain(cap.inner_profile_loops.iter().flatten())
        .any(crate::profile_curves::RecoveredProfileCurve::is_spline)
    {
        return None;
    }
    let volume = (signed_area_curves(&cap.profile_curves).abs()
        - cap
            .inner_profile_loops
            .iter()
            .map(|r| signed_area_curves(r).abs())
            .sum::<f64>())
        * cap.height_mm;
    (volume.is_finite() && volume > 0.0).then_some(volume)
}

struct Candidate {
    cap: RecoveredSourceCapExtrusion,
    volume: f64,
    witnesses: Vec<usize>,
}

fn select_caps(
    caps: Vec<RecoveredSourceCapExtrusion>,
    witnesses: &[[f64; 3]],
) -> Option<Vec<RecoveredSourceCapExtrusion>> {
    if witnesses.is_empty() {
        return None;
    }
    let mut candidates = caps
        .into_iter()
        .filter_map(|cap| {
            let volume = volume(&cap)?;
            if witnesses.iter().any(|&p| strict_contains(&cap, p)) {
                return None;
            }
            let covered = witnesses
                .iter()
                .enumerate()
                .filter_map(|(i, &p)| cap_contains(&cap, p).then_some(i))
                .collect::<Vec<_>>();
            (!covered.is_empty()).then_some(Candidate {
                cap,
                volume,
                witnesses: covered,
            })
        })
        .collect::<Vec<_>>();
    let mut uncovered = vec![true; witnesses.len()];
    let mut remaining = witnesses.len();
    let mut selected = Vec::new();
    while remaining > 0 {
        let (index, gain) = candidates
            .iter()
            .enumerate()
            .map(|(i, c)| (i, c.witnesses.iter().filter(|&&j| uncovered[j]).count()))
            .max_by(|(i, a), (j, b)| {
                a.cmp(b)
                    .then_with(|| candidates[*i].volume.total_cmp(&candidates[*j].volume))
                    .then_with(|| {
                        candidates[*i]
                            .cap
                            .cap_face_id
                            .cmp(&candidates[*j].cap.cap_face_id)
                    })
            })?;
        if gain == 0 {
            return None;
        }
        let candidate = candidates.swap_remove(index);
        for i in candidate.witnesses {
            uncovered[i] = false;
        }
        remaining -= gain;
        selected.push(candidate.cap);
    }
    (selected.len() >= 2).then_some(selected)
}

pub(super) fn recover<K: CadKernel>(
    source: &SourceRecovery<'_, K>,
    ceiling: Option<u64>,
) -> Result<Option<CadFragment>> {
    let SourceRecovery {
        solid_id,
        entities,
        index,
        ..
    } = *source;
    let Some(solid) = ReducedSolid::from_step(solid_id, entities, index) else {
        return Ok(None);
    };
    if solid.validate_closed_topology().is_none() {
        return Ok(None);
    }
    let mut seen = HashSet::new();
    let witnesses = solid
        .faces()
        .flat_map(ReducedFace::loop_edges)
        .filter(|edge| seen.insert(edge.edge_id))
        .flat_map(|edge| {
            [edge.start_mm, edge.end_mm]
                .into_iter()
                .chain(source_boundary_edge_witness_points(edge))
        })
        .collect::<Vec<_>>();
    let caps = partial_source_cap_extrusion_candidates_with_index(solid_id, entities, index);
    let Some(caps) = select_caps(caps, &witnesses) else {
        return Ok(None);
    };
    // For more than two components, assemble smaller features before the
    // large host, in a shared local frame. This avoids repeatedly classifying
    // the largest curved shell against small, nearly coincident attachments.
    let mut compact = caps.clone();
    compact.sort_by(|a, b| {
        volume(a)
            .unwrap_or(0.0)
            .total_cmp(&volume(b).unwrap_or(0.0))
            .then_with(|| a.cap_face_id.cmp(&b.cap_face_id))
    });
    let centered = AssemblyFrame::centered(
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        witnesses.iter().copied(),
    );
    let mut proposals = Vec::new();
    if caps.len() > 2
        && let Some(frame) = centered
    {
        proposals.push((compact.as_slice(), Some(frame)));
    }
    proposals.push((caps.as_slice(), None));
    if caps.len() == 2
        && let Some(frame) = centered
    {
        proposals.push((compact.as_slice(), Some(frame)));
    }
    let mut admitted = None;
    for (parts, frame) in proposals {
        if let Some((model, root)) = proven_cap_assembly(parts, frame, source, ceiling)? {
            admitted = Some((model, root, parts));
            break;
        }
    }
    let Some((mut model, root, caps)) = admitted else {
        return Ok(None);
    };
    let face_ids = super::set_proven_source_boundary(&mut model, root, &solid)?;
    Ok(Some(CadFragment {
        source: CadFragmentSource::SourceCapExtrusionUnion {
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
    use crate::cad_ir::CadNode;
    use crate::cad_kernel::monstertruck::MonstertruckKernel;
    use crate::cad_recovery::rebuilt_kernel_boundary_matches_source;

    #[test]
    fn source_cap_cover_recovers_bent_sheets() -> Result<()> {
        for (bytes, id, count) in [
            (
                include_bytes!("../../validation/fixtures/native_cap_cover_0712fbbd.step")
                    .as_slice(),
                101,
                2,
            ),
            (
                include_bytes!("../../validation/fixtures/native_centered_cap_1968bb53.step")
                    .as_slice(),
                115,
                3,
            ),
        ] {
            let kernel = MonstertruckKernel::project_standard();
            let parsed = monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange;
            let entities = &parsed.data[0].entities;
            let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(
                report.constructively_recovered_solid_count, 1,
                "{:?}",
                report.kernel_recovery_failures
            );
            let [selected] = report.extrusion_core_unions.as_slice() else {
                anyhow::bail!("expected one source-cap extrusion union");
            };
            let fragment = &selected.fragment;
            assert!(
                matches!(&fragment.source, CadFragmentSource::SourceCapExtrusionUnion { solid_id, cap_face_ids, .. } if *solid_id == id && cap_face_ids.len() == count),
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
                count
            );
            let shape = kernel.evaluate(&fragment.model, fragment.root)?;
            assert!(kernel.summarize(&shape).is_single_solid());
            assert!(rebuilt_kernel_boundary_matches_source(
                &kernel, &shape, id, entities
            )?);
        }
        Ok(())
    }
}
