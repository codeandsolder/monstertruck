//! Stable kernel identity assigned while importing STEP topology.
//!
//! STEP entity numbers and kernel [`StableId`]s are deliberately separate
//! namespaces.  These records bridge the two after a successful conversion
//! without making source IDs part of the kernel's identity scheme.

use super::ShellLoadReport;
use monstertruck_topology::StableId;
use serde::{Deserialize, Serialize};

/// One kept STEP entity and the kernel-local stable ID assigned to its imported
/// topology occurrence.
///
/// A source entity can appear more than once in a malformed or deliberately
/// shared input.  In that case the journal can contain repeated `source_id`
/// entries carrying the same `stable_id`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct StepEntityProvenance {
    /// STEP entity number in the imported DATA section.
    pub source_id: u64,
    /// Kernel-local identity assigned during solid import.
    pub stable_id: StableId,
}

/// Source identity retained for one imported face occurrence.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct StepFaceSource {
    /// The underlying `FACE_SURFACE` / `ADVANCED_FACE` entity when referenced.
    ///
    /// This is the source identity corresponding to the kernel face.  It is
    /// intentionally distinct from an `ORIENTED_FACE` wrapper.
    pub face_surface_id: Option<u64>,
    /// The entity referenced directly by the shell's face list.
    ///
    /// For an `ORIENTED_FACE` this is the wrapper ID; for a direct face
    /// reference it equals `face_surface_id`; owned holders have no entity ID.
    pub shell_member_id: Option<u64>,
}

/// One kept face occurrence and the kernel identity assigned to it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct StepFaceProvenance {
    /// STEP source identity retained while resolving the shell member.
    pub source: StepFaceSource,
    /// Kernel-local identity assigned during solid import.
    pub stable_id: StableId,
}

/// Source-to-kernel identity journal for topology kept from one imported solid.
///
/// Entries are in compressed-topology occurrence order across boundary shells.
/// Repeated STEP source entities therefore remain visible rather than being
/// collapsed into a map.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct StepTopologyProvenance {
    /// Imported `VERTEX_POINT` occurrences.
    pub vertices: Vec<StepEntityProvenance>,
    /// Imported `EDGE_CURVE` occurrences.
    pub edges: Vec<StepEntityProvenance>,
    /// Imported face occurrences, retaining both face and shell-member source IDs.
    pub faces: Vec<StepFaceProvenance>,
}

/// Identity-bearing result of importing one STEP solid.
///
/// `shell_reports` remains positional with `solid.boundaries`; provenance is a
/// separate concern describing what WAS kept, while the reports explain what
/// was lost during conversion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StepSolidImport<S> {
    /// Imported compressed solid carrying stable-ID sidecars and allocator state.
    pub solid: S,
    /// One load-loss report per boundary shell.
    pub shell_reports: Vec<ShellLoadReport>,
    /// Source-to-kernel topology identity journal.
    pub provenance: StepTopologyProvenance,
}
