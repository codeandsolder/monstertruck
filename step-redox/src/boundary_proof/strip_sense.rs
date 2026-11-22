//! Explicit proof-only reconciliation of inconsistent spline-strip face sense.
//!
//! The closed shell's oriented boundary traversal is redundant with the
//! surface chart and face sense on a completely trimmed embedded strip.
//! Only a proof of the opposite complete traversal permits a metadata
//! correction. Nothing in the authored geometry or topology may be edited.
use super::{
    AnalyticBoundaryDiagnostic, CellWitnesses, analytic_boundary_diagnostic, brep,
    oriented_support_key, spline_strip, strict_closed_shell_faces,
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use step_p21::ast::{EntityInstance, Parameter};

/// One contradictory face-sense flag; IDs are local to the original DATA section.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct StripFaceSenseRepair {
    pub face_id: u64,
    pub surface_id: u64,
    pub original_same_sense: bool,
    pub reconciled_same_sense: bool,
}

pub(crate) struct ReconciledStripBoundary {
    pub entities: Vec<EntityInstance>,
    pub diagnostic: AnalyticBoundaryDiagnostic,
    pub repairs: Vec<StripFaceSenseRepair>,
}

/// Build a distinct proof view. Every unsupported face must have an exact,
/// uniquely reversed spline-strip trim certificate; every other face retains
/// its original semantics. The final complete boundary proof is mandatory.
/// An empty repair plan is not an orientation-reconciled certificate.
pub(crate) fn reconcile_strip_face_senses(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<ReconciledStripBoundary> {
    let faces = strict_closed_shell_faces(solid_id, entities, index).ok()?;
    let mut repairs = Vec::new();
    let mut cells = CellWitnesses::default();
    for face_id in faces {
        if oriented_support_key(face_id, entities, index, &mut cells).is_some() {
            continue;
        }
        let surface_id = brep::face_surface(face_id, entities, index)?;
        let forward = match brep::face_sense(face_id, entities, index)?.as_str() {
            "T" => true,
            "F" => false,
            _ => return None,
        };
        if !spline_strip::reversed_face_sense_is_proven(
            face_id, surface_id, forward, entities, index,
        ) {
            return None;
        }
        repairs.push(StripFaceSenseRepair {
            face_id,
            surface_id,
            original_same_sense: forward,
            reconciled_same_sense: !forward,
        });
    }
    if repairs.is_empty() {
        return None;
    }
    repairs.sort_unstable_by_key(|r| r.face_id);
    let mut reconciled = entities.to_vec();
    for repair in &repairs {
        let EntityInstance::Simple { id, record } =
            reconciled.get_mut(*index.get(&repair.face_id)?)?
        else {
            return None;
        };
        if *id != repair.face_id
            || !matches!(record.name.as_str(), "ADVANCED_FACE" | "FACE_SURFACE")
        {
            return None;
        }
        let Parameter::List(params) = &mut record.parameter else {
            return None;
        };
        let Parameter::Enumeration(sense) = params.get_mut(3)? else {
            return None;
        };
        if sense != if repair.original_same_sense { "T" } else { "F" } {
            return None;
        }
        *sense = if repair.reconciled_same_sense {
            "T"
        } else {
            "F"
        }
        .to_owned();
    }
    let diagnostic = analytic_boundary_diagnostic(solid_id, &reconciled, index);
    diagnostic.supported.then_some(ReconciledStripBoundary {
        entities: reconciled,
        diagnostic,
        repairs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_graph::{build_index, simple_record_mut};
    use anyhow::{Context, Result};
    use monstertruck_io::step::load::ParsedStepSource;

    const HOUSING: &[u8] =
        include_bytes!("../../validation/fixtures/native_relief_housing_0914dda4.step");
    const CONSISTENT: &[u8] =
        include_bytes!("../../validation/fixtures/native_relief_box_0b5fc9e7.step");

    #[test]
    fn strip_sense_reconciliation_only_changes_the_proven_flags() -> Result<()> {
        let source = ParsedStepSource::parse(HOUSING)?.exchange;
        let entities = &source.data[0].entities;
        let index = build_index(entities);
        assert!(!analytic_boundary_diagnostic(4138, entities, &index).supported);
        let reconciled = reconcile_strip_face_senses(4138, entities, &index)
            .context("complete opposite trim proof")?;
        assert!(reconciled.diagnostic.supported);
        assert_eq!(reconciled.repairs.len(), 10);
        assert!(reconcile_strip_face_senses(4138, &reconciled.entities, &index).is_none());
        let mut restored = source.clone();
        restored.data[0].entities = reconciled.entities;
        for repair in reconciled.repairs {
            let original = &entities[index[&repair.face_id]];
            let revised = &mut restored.data[0].entities[index[&repair.face_id]];
            let original_record = super::super::simple_record(original).context("raw face")?;
            let revised_record = simple_record_mut(revised).context("reconciled face")?;
            let (Parameter::List(original_params), Parameter::List(revised_params)) =
                (&original_record.parameter, &mut revised_record.parameter)
            else {
                anyhow::bail!("invalid source face parameter list");
            };
            assert_ne!(original_params[3], revised_params[3]);
            revised_params[3] = original_params[3].clone();
        }
        assert_eq!(
            crate::write_exchange(&source)?,
            crate::write_exchange(&restored)?
        );
        let consistent = ParsedStepSource::parse(CONSISTENT)?.exchange;
        let entities = &consistent.data[0].entities;
        assert!(reconcile_strip_face_senses(1202, entities, &build_index(entities)).is_none());
        Ok(())
    }

    #[test]
    fn strip_sense_reconciliation_rejects_changed_geometry_or_broken_loop() -> Result<()> {
        let original = ParsedStepSource::parse(HOUSING)?.exchange;
        let entities = &original.data[0].entities;
        let index = build_index(entities);
        let mut changed = entities.clone();
        // Change an interior control of spline #3009, without changing its
        // surface's controls or endpoints. Reversing a flag cannot heal this.
        let record = simple_record_mut(&mut changed[index[&8140]]).context("cap control")?;
        let Parameter::List(params) = &mut record.parameter else {
            anyhow::bail!("control parameters");
        };
        let Parameter::List(coords) = &mut params[1] else {
            anyhow::bail!("control coordinates");
        };
        let Parameter::Real(x) = &mut coords[0] else {
            anyhow::bail!("control x-coordinate");
        };
        *x += 1.0e-4;
        assert!(reconcile_strip_face_senses(4138, &changed, &index).is_none());

        let mut disconnected = entities.clone();
        // Break a distinct already-supported planar face's ordered loop.
        // The repaired strips alone remain valid, but the shell must refuse.
        let faces =
            strict_closed_shell_faces(4138, entities, &index).map_err(anyhow::Error::msg)?;
        let planar = faces
            .into_iter()
            .find(|&face| {
                matches!(
                    oriented_support_key(face, entities, &index, &mut CellWitnesses::default()),
                    Some(super::super::OrientedSupportKey::Plane { .. })
                )
            })
            .context("planar face")?;
        let edge = brep::face_loops(planar, entities, &index).context("planar loops")?[0].edges[0]
            .oriented_edge_id;
        let record = simple_record_mut(&mut disconnected[index[&edge]]).context("planar coedge")?;
        let Parameter::List(params) = &mut record.parameter else {
            anyhow::bail!("coedge parameters");
        };
        let Parameter::Enumeration(sense) = params.last_mut().context("missing coedge sense")?
        else {
            anyhow::bail!("coedge sense");
        };
        *sense = if sense == "T" { "F" } else { "T" }.to_owned();
        assert!(reconcile_strip_face_senses(4138, &disconnected, &index).is_none());
        Ok(())
    }
}
