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
