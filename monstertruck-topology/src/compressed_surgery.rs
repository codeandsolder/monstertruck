//! Transactional StableId-addressed surgery for exact-trim compressed topology.
//!
//! Geometry proof belongs above this layer. Callers provide a fully normalized patch:
//! exact replacement curves, explicit endpoint topology, and face-boundary splices. This module
//! resolves every retained target through StableIds, allocates identity only for genuinely new
//! topology, regenerates exact face-local trims for every changed edge use, compacts dead storage,
//! cleans dangling topology attributes, validates the resulting runtime topology, and commits only
//! if the complete transaction succeeds.

use crate::compress::{
    CompressedEdge, CompressedEdgeUse, CompressedTrimmedShell, CompressedTrimmedSolid,
    StableIdEntity, StableIdLocation,
};
use crate::trimmed::TrimmedSolid;
use crate::{ExactParameterBoundary2D, StableId, StableIdExhausted};
use rustc_hash::FxHashSet as HashSet;
use std::collections::BTreeSet;
use thiserror::Error;

/// Vertex reference used by a normalized topology patch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchVertexRef {
    /// Reuse one uniquely-addressable existing vertex.
    Existing(StableId),
    /// Reference a vertex from [`CompressedTrimmedTopologyPatch::new_vertices`].
    New(usize),
}

/// Edge reference used when rebuilding one face boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchEdgeRef {
    /// Reuse one uniquely-addressable existing edge.
    Existing(StableId),
    /// Reference an edge from [`CompressedTrimmedTopologyPatch::new_edges`].
    New(usize),
}

/// One oriented edge use inserted by a face-boundary splice.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PatchEdgeUse {
    /// Shared edge to use.
    pub edge: PatchEdgeRef,
    /// Local face orientation of the edge.
    pub orientation: bool,
}

/// New shared edge created by a normalized patch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NewCompressedEdge<C> {
    /// Endpoints of the new edge.
    pub vertices: (PatchVertexRef, PatchVertexRef),
    /// Exact shared 3D edge curve.
    pub curve: C,
}

/// Geometry/topology replacement for one retained shared edge.
///
/// The `StableId` of the edge itself is preserved. Callers must provide the complete replacement
/// curve rather than relying on this topology layer to extend or otherwise reinterpret arbitrary
/// curve types.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RetainedCompressedEdge<C> {
    /// Existing edge to rewrite.
    pub edge: StableId,
    /// Replacement endpoints.
    pub vertices: (PatchVertexRef, PatchVertexRef),
    /// Exact replacement 3D curve.
    pub curve: C,
}

/// Replaces one contiguous cyclic run of edges in one face boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedFaceSplice {
    /// Face whose boundary is rebuilt.
    pub face: StableId,
    /// Existing edge `StableId`s forming exactly one contiguous run in one wire.
    pub remove_edges: Vec<StableId>,
    /// Ordered replacement edge uses. May be empty when the neighboring retained edges are
    /// rewritten to meet directly.
    pub replacement: Vec<PatchEdgeUse>,
}

/// Proof-independent normalized patch for one compressed trimmed boundary shell.
///
/// Every existing entity is addressed by `StableId`. The patch may create vertices/edges, rewrite
/// retained shared edges, splice arbitrary affected faces, remove faces, and require selected old
/// edges to become completely unused before canonical compaction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedTrimmedTopologyPatch<P, C> {
    /// Faces deleted by the operation.
    pub remove_faces: Vec<StableId>,
    /// New derived vertices. Fresh `StableId`s are allocated transactionally.
    pub new_vertices: Vec<P>,
    /// New derived shared edges. Fresh `StableId`s are allocated transactionally.
    pub new_edges: Vec<NewCompressedEdge<C>>,
    /// Existing edges whose identity survives but whose exact geometry/endpoints change.
    pub rewrite_edges: Vec<RetainedCompressedEdge<C>>,
    /// Face-boundary edits, including carrier and terminal/end-face repairs.
    pub face_splices: Vec<CompressedFaceSplice>,
    /// Existing edges that must have no surviving uses after face edits/removals.
    ///
    /// This is the fail-closed guard against an incomplete inverse-operation plan: if a treatment
    /// interface or terminal edge is still referenced by any surviving face, the transaction is
    /// rejected rather than leaving a non-manifold partial heal.
    pub require_unused_edges: Vec<StableId>,
}

/// Fresh identities allocated by a successful topology patch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedTrimmedTopologyPatchResult {
    /// `StableId`s assigned to `new_vertices`, in the same order.
    pub new_vertex_ids: Vec<StableId>,
    /// `StableId`s assigned to `new_edges`, in the same order.
    pub new_edge_ids: Vec<StableId>,
}

/// Failure of a transactional compressed-topology patch.
#[derive(Debug, Error)]
pub enum CompressedTopologyPatchError {
    /// Existing compressed topology was malformed or the final runtime topology was invalid.
    #[error(transparent)]
    Topology(#[from] crate::errors::Error),
    /// No fresh `StableId` can be allocated.
    #[error(transparent)]
    StableIdExhausted(#[from] StableIdExhausted),
    /// A requested `StableId` is absent or unassigned.
    #[error("topology patch target {id} is missing")]
    MissingTarget { id: StableId },
    /// A requested `StableId` occurs more than once and cannot identify one mutation target.
    #[error("topology patch target {id} is ambiguous")]
    AmbiguousTarget { id: StableId },
    /// A `StableId` resolved to the wrong topology entity class.
    #[error("topology patch target {id} is {actual:?}, expected {expected:?}")]
    WrongEntity {
        id: StableId,
        expected: StableIdEntity,
        actual: StableIdEntity,
    },
    /// One patch attempted to cross compressed boundary-shell storage domains.
    #[error("topology patch spans boundary shells {expected} and {actual}")]
    CrossBoundary { expected: usize, actual: usize },
    /// A patch contains no existing target from which its boundary shell can be resolved.
    #[error("topology patch must address at least one existing StableId target")]
    NoBoundaryTarget,
    /// A new-entity reference points beyond the corresponding patch vector.
    #[error("topology patch {entity} reference {index} is out of bounds for {len} new entities")]
    NewEntityOutOfBounds {
        entity: &'static str,
        index: usize,
        len: usize,
    },
    /// The same retained edge was rewritten more than once.
    #[error("topology patch rewrites retained edge {id} more than once")]
    DuplicateEdgeRewrite { id: StableId },
    /// A face scheduled for deletion was also edited.
    #[error("topology patch both removes and splices face {id}")]
    RemovedFaceSpliced { id: StableId },
    /// One face splice did not describe one unique contiguous cyclic run.
    #[error("face {face} does not contain removed edges as one unique contiguous boundary run")]
    NonContiguousFaceSplice { face: StableId },
    /// A face splice would leave an empty wire.
    #[error("face {face} splice would leave an empty boundary wire")]
    EmptyFaceBoundary { face: StableId },
    /// An edge promised to disappear is still used after all patch edits.
    #[error("edge {edge} is still referenced after normalized topology edits")]
    RequiredEdgeStillUsed { edge: StableId },
    /// Exact face-local trim regeneration failed for changed geometry.
    #[error("exact face-local trim is unavailable for changed edge {edge} on face {face}")]
    MissingExactTrim { edge: StableId, face: StableId },
}

#[derive(Clone, Copy, Debug)]
struct ResolvedLocation {
    boundary: usize,
    index: usize,
}

#[derive(Clone, Copy, Debug)]
enum ResolvedVertexRef {
    Existing(usize),
    New(usize),
}

#[derive(Clone, Copy, Debug)]
enum ResolvedEdgeRef {
    Existing(usize),
    New(usize),
}

#[derive(Clone, Debug)]
struct ResolvedEdgeRewrite<C> {
    edge_index: usize,
    vertices: (ResolvedVertexRef, ResolvedVertexRef),
    curve: C,
}

#[derive(Clone, Debug)]
struct ResolvedFaceSplice {
    face_id: StableId,
    face_index: usize,
    remove_edges: Vec<StableId>,
    replacement: Vec<(ResolvedEdgeRef, bool)>,
}

#[derive(Clone, Debug)]
struct ResolvedPatch<C> {
    boundary: usize,
    remove_faces: Vec<(StableId, usize)>,
    new_edges: Vec<(ResolvedVertexRef, ResolvedVertexRef, C)>,
    rewrite_edges: Vec<ResolvedEdgeRewrite<C>>,
    face_splices: Vec<ResolvedFaceSplice>,
    require_unused_edges: Vec<(StableId, usize)>,
}

#[derive(Default)]
struct TopologyIds {
    vertices: BTreeSet<StableId>,
    edges: BTreeSet<StableId>,
    faces: BTreeSet<StableId>,
}

impl<P, C, S, T> CompressedTrimmedSolid<P, C, S, T>
where
    P: Clone,
    C: Clone + ExactParameterBoundary2D<S, BoundaryCurve = T>,
    S: Clone,
    T: Clone,
{
    /// Applies a fully normalized topology patch transactionally.
    ///
    /// Existing targets are resolved against the unmodified solid and must be unique. Every change
    /// is made on a clone, changed edge uses receive newly generated exact face-local trims, dead
    /// topology is compacted canonically, and the complete runtime topology is validated before
    /// `self` is replaced.
    ///
    /// # Errors
    /// Returns [`CompressedTopologyPatchError`] without modifying `self` if target resolution,
    /// normalized-plan structure, exact trim regeneration, compaction, `StableId` allocation, or final
    /// topology validation fails.
    pub fn apply_topology_patch(
        &mut self,
        patch: &CompressedTrimmedTopologyPatch<P, C>,
    ) -> Result<CompressedTrimmedTopologyPatchResult, CompressedTopologyPatchError> {
        let mut candidate = self.clone();
        validate_compressed_references(&candidate)?;
        let before_ids = topology_ids(&candidate);
        let resolved = resolve_patch(&candidate, patch)?;
        let result = apply_resolved_patch(&mut candidate, patch, resolved)?;
        let after_ids = topology_ids(&candidate);
        cleanup_removed_attributes(&mut candidate, &before_ids, &after_ids);
        let _: TrimmedSolid<P, C, S, T> = candidate.clone().try_into()?;
        *self = candidate;
        Ok(result)
    }
}

fn validate_compressed_references<P: Clone, C: Clone, S: Clone, T: Clone>(
    solid: &CompressedTrimmedSolid<P, C, S, T>,
) -> Result<(), crate::errors::Error> {
    for shell in &solid.boundaries {
        let mut clone = shell.clone();
        clone.compact_unused_topology()?;
    }
    Ok(())
}

fn resolve_target<P, C, S, T>(
    solid: &CompressedTrimmedSolid<P, C, S, T>,
    id: StableId,
    expected: StableIdEntity,
) -> Result<ResolvedLocation, CompressedTopologyPatchError> {
    match solid.location_by_stable_id(id)? {
        StableIdLocation::Missing => Err(CompressedTopologyPatchError::MissingTarget { id }),
        StableIdLocation::Ambiguous => Err(CompressedTopologyPatchError::AmbiguousTarget { id }),
        StableIdLocation::Unique {
            entity,
            boundary,
            index,
        } if entity == expected => Ok(ResolvedLocation { boundary, index }),
        StableIdLocation::Unique { entity, .. } => Err(CompressedTopologyPatchError::WrongEntity {
            id,
            expected,
            actual: entity,
        }),
    }
}

const fn accept_boundary(
    selected: &mut Option<usize>,
    location: ResolvedLocation,
) -> Result<(), CompressedTopologyPatchError> {
    match *selected {
        Some(expected) if expected != location.boundary => {
            Err(CompressedTopologyPatchError::CrossBoundary {
                expected,
                actual: location.boundary,
            })
        }
        Some(_) => Ok(()),
        None => {
            *selected = Some(location.boundary);
            Ok(())
        }
    }
}

fn resolve_vertex_ref<P, C, S, T>(
    solid: &CompressedTrimmedSolid<P, C, S, T>,
    reference: PatchVertexRef,
    new_len: usize,
    boundary: &mut Option<usize>,
) -> Result<ResolvedVertexRef, CompressedTopologyPatchError> {
    match reference {
        PatchVertexRef::Existing(id) => {
            let location = resolve_target(solid, id, StableIdEntity::Vertex)?;
            accept_boundary(boundary, location)?;
            Ok(ResolvedVertexRef::Existing(location.index))
        }
        PatchVertexRef::New(index) if index < new_len => Ok(ResolvedVertexRef::New(index)),
        PatchVertexRef::New(index) => Err(CompressedTopologyPatchError::NewEntityOutOfBounds {
            entity: "vertex",
            index,
            len: new_len,
        }),
    }
}

fn resolve_edge_ref<P, C, S, T>(
    solid: &CompressedTrimmedSolid<P, C, S, T>,
    reference: PatchEdgeRef,
    new_len: usize,
    boundary: &mut Option<usize>,
) -> Result<ResolvedEdgeRef, CompressedTopologyPatchError> {
    match reference {
        PatchEdgeRef::Existing(id) => {
            let location = resolve_target(solid, id, StableIdEntity::Edge)?;
            accept_boundary(boundary, location)?;
            Ok(ResolvedEdgeRef::Existing(location.index))
        }
        PatchEdgeRef::New(index) if index < new_len => Ok(ResolvedEdgeRef::New(index)),
        PatchEdgeRef::New(index) => Err(CompressedTopologyPatchError::NewEntityOutOfBounds {
            entity: "edge",
            index,
            len: new_len,
        }),
    }
}

fn resolve_patch<P, C: Clone, S, T>(
    solid: &CompressedTrimmedSolid<P, C, S, T>,
    patch: &CompressedTrimmedTopologyPatch<P, C>,
) -> Result<ResolvedPatch<C>, CompressedTopologyPatchError> {
    let mut boundary = None;
    let mut remove_faces = Vec::with_capacity(patch.remove_faces.len());
    let mut removed_face_ids = BTreeSet::new();
    for &id in &patch.remove_faces {
        let location = resolve_target(solid, id, StableIdEntity::Face)?;
        accept_boundary(&mut boundary, location)?;
        if removed_face_ids.insert(id) {
            remove_faces.push((id, location.index));
        }
    }

    let mut new_edges = Vec::with_capacity(patch.new_edges.len());
    for edge in &patch.new_edges {
        let front = resolve_vertex_ref(
            solid,
            edge.vertices.0,
            patch.new_vertices.len(),
            &mut boundary,
        )?;
        let back = resolve_vertex_ref(
            solid,
            edge.vertices.1,
            patch.new_vertices.len(),
            &mut boundary,
        )?;
        new_edges.push((front, back, edge.curve.clone()));
    }

    let mut rewritten_ids = BTreeSet::new();
    let mut rewrite_edges = Vec::with_capacity(patch.rewrite_edges.len());
    for edge in &patch.rewrite_edges {
        if !rewritten_ids.insert(edge.edge) {
            return Err(CompressedTopologyPatchError::DuplicateEdgeRewrite { id: edge.edge });
        }
        let location = resolve_target(solid, edge.edge, StableIdEntity::Edge)?;
        accept_boundary(&mut boundary, location)?;
        let front = resolve_vertex_ref(
            solid,
            edge.vertices.0,
            patch.new_vertices.len(),
            &mut boundary,
        )?;
        let back = resolve_vertex_ref(
            solid,
            edge.vertices.1,
            patch.new_vertices.len(),
            &mut boundary,
        )?;
        rewrite_edges.push(ResolvedEdgeRewrite {
            edge_index: location.index,
            vertices: (front, back),
            curve: edge.curve.clone(),
        });
    }

    let mut face_splices = Vec::with_capacity(patch.face_splices.len());
    for splice in &patch.face_splices {
        if removed_face_ids.contains(&splice.face) {
            return Err(CompressedTopologyPatchError::RemovedFaceSpliced { id: splice.face });
        }
        let face = resolve_target(solid, splice.face, StableIdEntity::Face)?;
        accept_boundary(&mut boundary, face)?;
        for &edge in &splice.remove_edges {
            let location = resolve_target(solid, edge, StableIdEntity::Edge)?;
            accept_boundary(&mut boundary, location)?;
        }
        let replacement = splice
            .replacement
            .iter()
            .map(|edge_use| {
                Ok((
                    resolve_edge_ref(solid, edge_use.edge, patch.new_edges.len(), &mut boundary)?,
                    edge_use.orientation,
                ))
            })
            .collect::<Result<Vec<_>, CompressedTopologyPatchError>>()?;
        face_splices.push(ResolvedFaceSplice {
            face_id: splice.face,
            face_index: face.index,
            remove_edges: splice.remove_edges.clone(),
            replacement,
        });
    }

    let mut require_unused_edges = Vec::with_capacity(patch.require_unused_edges.len());
    let mut required_ids = BTreeSet::new();
    for &id in &patch.require_unused_edges {
        if !required_ids.insert(id) {
            continue;
        }
        let location = resolve_target(solid, id, StableIdEntity::Edge)?;
        accept_boundary(&mut boundary, location)?;
        require_unused_edges.push((id, location.index));
    }

    let boundary = boundary.ok_or(CompressedTopologyPatchError::NoBoundaryTarget)?;
    Ok(ResolvedPatch {
        boundary,
        remove_faces,
        new_edges,
        rewrite_edges,
        face_splices,
        require_unused_edges,
    })
}

fn max_assigned_id<P, C, S, T>(solid: &CompressedTrimmedSolid<P, C, S, T>) -> u64 {
    solid
        .boundaries
        .iter()
        .flat_map(|shell| {
            shell
                .vertex_stable_ids
                .iter()
                .flatten()
                .chain(shell.edge_stable_ids.iter().flatten())
                .chain(shell.face_stable_ids.iter().flatten())
        })
        .filter(|id| id.is_assigned())
        .map(|id| id.raw())
        .max()
        .unwrap_or(0)
}

fn allocate_fresh_ids<P, C, S, T>(
    solid: &mut CompressedTrimmedSolid<P, C, S, T>,
    count: usize,
) -> Result<Vec<StableId>, CompressedTopologyPatchError> {
    if count == 0 {
        return Ok(Vec::new());
    }
    let max_id = max_assigned_id(solid);
    let allocator = solid.id_allocator.get_or_insert_with(Default::default);
    allocator.ensure_above(max_id)?;
    allocator.ensure_capacity(count)?;
    (0..count)
        .map(|_| {
            allocator
                .allocate()
                .map_err(CompressedTopologyPatchError::from)
        })
        .collect()
}

const fn resolved_vertex_index(reference: ResolvedVertexRef, new_vertex_base: usize) -> usize {
    match reference {
        ResolvedVertexRef::Existing(index) => index,
        ResolvedVertexRef::New(index) => new_vertex_base + index,
    }
}

const fn resolved_edge_index(reference: ResolvedEdgeRef, new_edge_base: usize) -> usize {
    match reference {
        ResolvedEdgeRef::Existing(index) => index,
        ResolvedEdgeRef::New(index) => new_edge_base + index,
    }
}

fn apply_resolved_patch<P, C, S, T>(
    solid: &mut CompressedTrimmedSolid<P, C, S, T>,
    patch: &CompressedTrimmedTopologyPatch<P, C>,
    resolved: ResolvedPatch<C>,
) -> Result<CompressedTrimmedTopologyPatchResult, CompressedTopologyPatchError>
where
    P: Clone,
    C: Clone + ExactParameterBoundary2D<S, BoundaryCurve = T>,
    S: Clone,
    T: Clone,
{
    let allocation_count = patch
        .new_vertices
        .len()
        .checked_add(patch.new_edges.len())
        .ok_or(StableIdExhausted)?;
    let allocated = allocate_fresh_ids(solid, allocation_count)?;
    let (new_vertex_ids, new_edge_ids) = allocated.split_at(patch.new_vertices.len());
    let new_vertex_ids = new_vertex_ids.to_vec();
    let new_edge_ids = new_edge_ids.to_vec();

    let shell = &mut solid.boundaries[resolved.boundary];
    let new_vertex_base = shell.vertices.len();
    let vertex_ids = shell
        .vertex_stable_ids
        .get_or_insert_with(|| vec![StableId::UNASSIGNED; shell.vertices.len()]);
    shell.vertices.extend(patch.new_vertices.iter().cloned());
    vertex_ids.extend(new_vertex_ids.iter().copied());

    let new_edge_base = shell.edges.len();
    let edge_ids = shell
        .edge_stable_ids
        .get_or_insert_with(|| vec![StableId::UNASSIGNED; shell.edges.len()]);
    for ((front, back, curve), stable_id) in resolved.new_edges.iter().zip(&new_edge_ids) {
        shell.edges.push(CompressedEdge {
            vertices: (
                resolved_vertex_index(*front, new_vertex_base),
                resolved_vertex_index(*back, new_vertex_base),
            ),
            curve: curve.clone(),
        });
        edge_ids.push(*stable_id);
    }

    let mut changed_edges = HashSet::default();
    for rewrite in &resolved.rewrite_edges {
        let edge = &mut shell.edges[rewrite.edge_index];
        edge.vertices = (
            resolved_vertex_index(rewrite.vertices.0, new_vertex_base),
            resolved_vertex_index(rewrite.vertices.1, new_vertex_base),
        );
        edge.curve = rewrite.curve.clone();
        changed_edges.insert(rewrite.edge_index);
    }
    changed_edges.extend(new_edge_base..shell.edges.len());

    for splice in &resolved.face_splices {
        let replacement = splice
            .replacement
            .iter()
            .map(|&(edge, orientation)| {
                let index = resolved_edge_index(edge, new_edge_base);
                changed_edges.insert(index);
                CompressedEdgeUse {
                    index,
                    orientation,
                    trim_curve: None,
                }
            })
            .collect::<Vec<_>>();
        splice_face_boundary_run(
            shell,
            splice.face_index,
            splice.face_id,
            &splice.remove_edges,
            replacement,
        )?;
    }

    let mut remove_faces = resolved.remove_faces;
    remove_faces.sort_unstable_by_key(|(_, index)| std::cmp::Reverse(*index));
    for (_, index) in remove_faces {
        shell.remove_face(index);
    }

    for &(stable_id, edge_index) in &resolved.require_unused_edges {
        if shell
            .faces
            .iter()
            .flat_map(|face| face.boundaries.iter())
            .flatten()
            .any(|edge_use| edge_use.index == edge_index)
        {
            return Err(CompressedTopologyPatchError::RequiredEdgeStillUsed { edge: stable_id });
        }
    }

    regenerate_exact_trims(shell, &changed_edges)?;
    shell.compact_unused_topology()?;

    Ok(CompressedTrimmedTopologyPatchResult {
        new_vertex_ids,
        new_edge_ids,
    })
}

fn splice_face_boundary_run<P, C, S, T: Clone>(
    shell: &mut CompressedTrimmedShell<P, C, S, T>,
    face_index: usize,
    face_id: StableId,
    remove_edges: &[StableId],
    replacement: Vec<CompressedEdgeUse<T>>,
) -> Result<(), CompressedTopologyPatchError> {
    if remove_edges.is_empty() {
        return Err(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id });
    }
    let Some(edge_ids) = shell.edge_stable_ids.as_deref() else {
        return Err(CompressedTopologyPatchError::MissingTarget {
            id: remove_edges[0],
        });
    };
    let remove = remove_edges.iter().copied().collect::<BTreeSet<_>>();
    if remove.len() != remove_edges.len() {
        return Err(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id });
    }

    let face = &mut shell.faces[face_index];
    let matching_wires = face
        .boundaries
        .iter()
        .enumerate()
        .filter_map(|(wire_index, wire)| {
            let count = wire
                .iter()
                .filter(|edge_use| {
                    edge_ids
                        .get(edge_use.index)
                        .is_some_and(|id| remove.contains(id))
                })
                .count();
            (count != 0).then_some((wire_index, count))
        })
        .collect::<Vec<_>>();
    let [(wire_index, count)] = matching_wires.as_slice() else {
        return Err(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id });
    };
    if *count != remove.len() {
        return Err(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id });
    }

    let wire = &face.boundaries[*wire_index];
    let edge_count = wire.len();
    if edge_count == 0 || remove.len() > edge_count {
        return Err(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id });
    }
    let removed_at = |index: usize| {
        edge_ids
            .get(wire[index].index)
            .is_some_and(|id| remove.contains(id))
    };
    let run_start = if remove.len() == edge_count {
        0
    } else {
        (0..edge_count)
            .find(|&index| removed_at(index) && !removed_at((index + edge_count - 1) % edge_count))
            .ok_or(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id })?
    };
    let mut run_len = 0;
    while run_len < edge_count && removed_at((run_start + run_len) % edge_count) {
        run_len += 1;
    }
    if run_len != remove.len() {
        return Err(CompressedTopologyPatchError::NonContiguousFaceSplice { face: face_id });
    }

    let remaining_len = edge_count - run_len;
    if remaining_len == 0 && replacement.is_empty() {
        return Err(CompressedTopologyPatchError::EmptyFaceBoundary { face: face_id });
    }
    let mut rebuilt = Vec::with_capacity(remaining_len + replacement.len());
    rebuilt.extend(replacement);
    let first_remaining = (run_start + run_len) % edge_count;
    for offset in 0..remaining_len {
        rebuilt.push(wire[(first_remaining + offset) % edge_count].clone());
    }
    face.boundaries[*wire_index] = rebuilt;
    Ok(())
}

fn regenerate_exact_trims<P, C, S, T>(
    shell: &mut CompressedTrimmedShell<P, C, S, T>,
    changed_edges: &HashSet<usize>,
) -> Result<(), CompressedTopologyPatchError>
where
    C: ExactParameterBoundary2D<S, BoundaryCurve = T>,
{
    let edge_ids = shell.edge_stable_ids.as_deref();
    let face_ids = shell.face_stable_ids.as_deref();
    let edges = &shell.edges;
    for (face_index, face) in shell.faces.iter_mut().enumerate() {
        for edge_use in face.boundaries.iter_mut().flatten() {
            if !changed_edges.contains(&edge_use.index) {
                continue;
            }
            let edge = &edges[edge_use.index];
            let trim_curve = edge
                .curve
                .exact_parameter_boundary_2d(&face.surface)
                .ok_or_else(|| CompressedTopologyPatchError::MissingExactTrim {
                    edge: edge_ids
                        .and_then(|ids| ids.get(edge_use.index))
                        .copied()
                        .unwrap_or(StableId::UNASSIGNED),
                    face: face_ids
                        .and_then(|ids| ids.get(face_index))
                        .copied()
                        .unwrap_or(StableId::UNASSIGNED),
                })?;
            edge_use.trim_curve = Some(trim_curve);
        }
    }
    Ok(())
}

fn topology_ids<P, C, S, T>(solid: &CompressedTrimmedSolid<P, C, S, T>) -> TopologyIds {
    let mut ids = TopologyIds::default();
    for shell in &solid.boundaries {
        if let Some(sidecar) = &shell.vertex_stable_ids {
            ids.vertices
                .extend(sidecar.iter().copied().filter(|id| id.is_assigned()));
        }
        if let Some(sidecar) = &shell.edge_stable_ids {
            ids.edges
                .extend(sidecar.iter().copied().filter(|id| id.is_assigned()));
        }
        if let Some(sidecar) = &shell.face_stable_ids {
            ids.faces
                .extend(sidecar.iter().copied().filter(|id| id.is_assigned()));
        }
    }
    ids
}

fn cleanup_removed_attributes<P, C, S, T>(
    solid: &mut CompressedTrimmedSolid<P, C, S, T>,
    before: &TopologyIds,
    after: &TopologyIds,
) {
    let Some(attributes) = solid.attributes.as_mut() else {
        return;
    };
    let removed_vertices = before
        .vertices
        .difference(&after.vertices)
        .copied()
        .collect::<BTreeSet<_>>();
    let removed_edges = before
        .edges
        .difference(&after.edges)
        .copied()
        .collect::<BTreeSet<_>>();
    let removed_faces = before
        .faces
        .difference(&after.faces)
        .copied()
        .collect::<BTreeSet<_>>();
    attributes.remove_topology_ids(&removed_vertices, &removed_edges, &removed_faces);
}

#[cfg(test)]
mod tests;
