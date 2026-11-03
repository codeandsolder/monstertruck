//! Serialized topological data exchange format
//!
//! Topological data structures in monstertruck are subject to editing and has complex reference relationships.
//! They are not suitable for direct serialization and must be converted to lighter and simpler data structures.
//! These structures, prefixed with `Compressed`, are a group of structures that are easy to serialize,
//! but not suitable for real-time shape editing.
//!
//! They directly reflect the results of parsing data from json or STEP, and all member variables are public.
//! Boundary connectivity and closure are checked when converting to proprietary data structures, `Vertex`, `Edge`, and so on.

#[cfg(test)]
use crate::*;
use monstertruck_core::{
    ContentHasher, DeterministicContentHash,
    bounding_box::{Bounded, BoundingBox},
    cgmath64::Point3,
};
use rustc_hash::FxHashMap as HashMap;
use serde::{Deserialize, Serialize};

use crate::{
    Debug, Edge, EdgeId, ExactParameterBoundary2D, Face, Hash, Hasher, KERNEL_EPSILON,
    ParametricCurve, ParametricSurface, Result, Shell, Solid, SolidAttributes, StableId,
    StableIdAllocator, Vertex, VertexId, Wire, errors::Error,
};

/// Serialized compressed edge
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompressedEdge<C> {
    /// vertices of the edge
    pub vertices: (usize, usize),
    /// curve geometry of the edge
    pub curve: C,
}

impl<C> CompressedEdge<C> {
    #[inline]
    fn create_edge<P>(self, v: &[Vertex<P>]) -> Result<Edge<P, C>> {
        let front = &v[self.vertices.0];
        let back = &v[self.vertices.1];
        Edge::new(front, back, self.curve)
    }
}

/// The index of an edge in `CompressedShell`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CompressedEdgeIndex {
    /// the index of the edge
    pub index: usize,
    /// the orientation of the edge
    pub orientation: bool,
}

impl From<(usize, bool)> for CompressedEdgeIndex {
    fn from((index, orientation): (usize, bool)) -> Self { Self { index, orientation } }
}

/// Face-local use of an edge in a compressed boundary.
///
/// This is the serialized equivalent of a coedge or edge-use: it refers to a shared
/// edge by `index`, stores the local `orientation`, and may also carry an exact trim
/// curve for this face-use.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompressedEdgeUse<T> {
    /// The index of the shared edge in the parent `CompressedShell`.
    pub index: usize,
    /// The local orientation of the edge-use.
    pub orientation: bool,
    /// The exact face-local trim carried by this edge-use, if available.
    pub trim_curve: Option<T>,
}

impl<T> From<CompressedEdgeIndex> for CompressedEdgeUse<T> {
    fn from(CompressedEdgeIndex { index, orientation }: CompressedEdgeIndex) -> Self {
        Self {
            index,
            orientation,
            trim_curve: None,
        }
    }
}

impl<T> From<(usize, bool, Option<T>)> for CompressedEdgeUse<T> {
    fn from((index, orientation, trim_curve): (usize, bool, Option<T>)) -> Self {
        Self {
            index,
            orientation,
            trim_curve,
        }
    }
}

/// Serialized compressed face
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompressedFace<S> {
    /// Boundaries of the face
    pub boundaries: Vec<Vec<CompressedEdgeIndex>>,
    /// orientation of the face
    pub orientation: bool,
    /// surface geometry of the face
    pub surface: S,
}

/// Serialized compressed face with exact face-local trim storage.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompressedTrimmedFace<S, T> {
    /// Boundaries of the face as oriented edge-uses.
    pub boundaries: Vec<Vec<CompressedEdgeUse<T>>>,
    /// orientation of the face
    pub orientation: bool,
    /// surface geometry of the face
    pub surface: S,
}

impl<S> CompressedFace<S> {
    fn create_face<P, C>(self, edges: &[Edge<P, C>]) -> Result<Face<P, C, S>> {
        let wires: Vec<Wire<P, C>> = self
            .boundaries
            .into_iter()
            .map(|wire| {
                wire.into_iter()
                    .map(|CompressedEdgeIndex { index, orientation }| {
                        if orientation {
                            edges[index].clone()
                        } else {
                            edges[index].inverse()
                        }
                    })
                    .collect()
            })
            .collect();
        let mut face = Face::new(wires, self.surface)?;
        if !self.orientation {
            face.invert();
        }
        Ok(face)
    }
}

/// Serialized compressed shell
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompressedShell<P, C, S> {
    /// all geometries of vertices
    pub vertices: Vec<P>,
    /// all geometries and end vertices of edges
    pub edges: Vec<CompressedEdge<C>>,
    /// all geometries and boundaries of faces
    pub faces: Vec<CompressedFace<S>>,
    /// Stable IDs for vertices (same order as `vertices`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vertex_stable_ids: Option<Vec<StableId>>,
    /// Stable IDs for edges (same order as `edges`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edge_stable_ids: Option<Vec<StableId>>,
    /// Stable IDs for faces (same order as `faces`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub face_stable_ids: Option<Vec<StableId>>,
}

/// Serialized compressed shell with exact face-local trim storage.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CompressedTrimmedShell<P, C, S, T> {
    /// all geometries of vertices
    pub vertices: Vec<P>,
    /// all geometries and end vertices of edges
    pub edges: Vec<CompressedEdge<C>>,
    /// all geometries and boundaries of faces, including optional face-local trims
    pub faces: Vec<CompressedTrimmedFace<S, T>>,
    /// Stable IDs for vertices (same order as `vertices`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vertex_stable_ids: Option<Vec<StableId>>,
    /// Stable IDs for edges (same order as `edges`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub edge_stable_ids: Option<Vec<StableId>>,
    /// Stable IDs for faces (same order as `faces`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub face_stable_ids: Option<Vec<StableId>>,
}

/// Serialized compressed solid
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CompressedSolid<P, C, S> {
    /// all boundaries of solid
    pub boundaries: Vec<CompressedShell<P, C, S>>,
    /// The stable ID allocator state, preserved across serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id_allocator: Option<StableIdAllocator>,
    /// Per-element attributes, preserved across serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub attributes: Option<SolidAttributes>,
}

/// Topology entity class addressed by stable-ID lookup.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StableIdEntity {
    /// Compressed vertex.
    Vertex,
    /// Compressed edge.
    Edge,
    /// Compressed face.
    Face,
}

/// Result of resolving one assigned [`StableId`] against compressed topology.
///
/// Import provenance may intentionally reuse a stable ID when one source entity
/// is referenced more than once. Mutators must therefore distinguish an absent
/// identity from a unique target and an ambiguous repeated occurrence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StableIdLocation {
    /// The ID is unassigned or does not occur in the compressed solid.
    Missing,
    /// Exactly one topology entity carries the ID.
    Unique {
        /// Topology entity class carrying the ID.
        entity: StableIdEntity,
        /// Boundary-shell index in the compressed solid.
        boundary: usize,
        /// Entity index inside that shell's parallel topology array.
        index: usize,
    },
    /// More than one topology entity carries the ID.
    Ambiguous,
}

/// Serialized compressed solid with exact face-local trim storage.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CompressedTrimmedSolid<P, C, S, T> {
    /// all boundaries of solid
    pub boundaries: Vec<CompressedTrimmedShell<P, C, S, T>>,
    /// The stable ID allocator state, preserved across serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id_allocator: Option<StableIdAllocator>,
    /// Per-element attributes, preserved across serialization.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub attributes: Option<SolidAttributes>,
}

trait StableIdIndexedShell {
    fn validate_stable_ids(&self) -> Result<()>;
    fn stable_ids(&self, entity: StableIdEntity) -> Option<&[StableId]>;
}

impl<P, C, S> StableIdIndexedShell for CompressedShell<P, C, S> {
    fn validate_stable_ids(&self) -> Result<()> {
        validate_stable_id_counts(
            self.vertex_stable_ids.as_deref(),
            self.vertices.len(),
            self.edge_stable_ids.as_deref(),
            self.edges.len(),
            self.face_stable_ids.as_deref(),
            self.faces.len(),
        )
    }

    fn stable_ids(&self, entity: StableIdEntity) -> Option<&[StableId]> {
        match entity {
            StableIdEntity::Vertex => self.vertex_stable_ids.as_deref(),
            StableIdEntity::Edge => self.edge_stable_ids.as_deref(),
            StableIdEntity::Face => self.face_stable_ids.as_deref(),
        }
    }
}

impl<P, C, S, T> StableIdIndexedShell for CompressedTrimmedShell<P, C, S, T> {
    fn validate_stable_ids(&self) -> Result<()> {
        validate_stable_id_counts(
            self.vertex_stable_ids.as_deref(),
            self.vertices.len(),
            self.edge_stable_ids.as_deref(),
            self.edges.len(),
            self.face_stable_ids.as_deref(),
            self.faces.len(),
        )
    }

    fn stable_ids(&self, entity: StableIdEntity) -> Option<&[StableId]> {
        match entity {
            StableIdEntity::Vertex => self.vertex_stable_ids.as_deref(),
            StableIdEntity::Edge => self.edge_stable_ids.as_deref(),
            StableIdEntity::Face => self.face_stable_ids.as_deref(),
        }
    }
}

fn find_stable_id_location<B: StableIdIndexedShell>(
    boundaries: &[B],
    target: StableId,
) -> Result<StableIdLocation> {
    if !target.is_assigned() {
        return Ok(StableIdLocation::Missing);
    }
    let mut found = None;
    for (boundary, shell) in boundaries.iter().enumerate() {
        shell.validate_stable_ids()?;
        for entity in [
            StableIdEntity::Vertex,
            StableIdEntity::Edge,
            StableIdEntity::Face,
        ] {
            let Some(ids) = shell.stable_ids(entity) else {
                continue;
            };
            for (index, &stable_id) in ids.iter().enumerate() {
                if stable_id != target {
                    continue;
                }
                if found.is_some() {
                    return Ok(StableIdLocation::Ambiguous);
                }
                found = Some(StableIdLocation::Unique {
                    entity,
                    boundary,
                    index,
                });
            }
        }
    }
    Ok(found.unwrap_or(StableIdLocation::Missing))
}

fn compact_items_with_stable_ids<T>(
    items: &mut Vec<T>,
    stable_ids: &mut Option<Vec<StableId>>,
    used: &[bool],
) -> Vec<usize> {
    debug_assert_eq!(items.len(), used.len());
    let old_items = std::mem::take(items);
    let old_ids = stable_ids.take();
    let mut retained = Vec::with_capacity(used.iter().filter(|&&keep| keep).count());
    let mut retained_ids = old_ids
        .as_ref()
        .map(|_| Vec::with_capacity(retained.capacity()));
    let mut remap = vec![usize::MAX; old_items.len()];
    for (old_index, item) in old_items.into_iter().enumerate() {
        if !used[old_index] {
            continue;
        }
        remap[old_index] = retained.len();
        retained.push(item);
        if let (Some(old_ids), Some(new_ids)) = (&old_ids, &mut retained_ids) {
            new_ids.push(old_ids[old_index]);
        }
    }
    *items = retained;
    *stable_ids = retained_ids;
    remap
}

fn compact_edge_vertex_storage<P, C>(
    vertices: &mut Vec<P>,
    edges: &mut Vec<CompressedEdge<C>>,
    vertex_stable_ids: &mut Option<Vec<StableId>>,
    edge_stable_ids: &mut Option<Vec<StableId>>,
    used_edges: &[bool],
) -> Result<Vec<usize>> {
    debug_assert_eq!(edges.len(), used_edges.len());
    let mut used_vertices = vec![false; vertices.len()];
    for (edge_index, edge) in edges.iter().enumerate() {
        if !used_edges[edge_index] {
            continue;
        }
        for vertex in <[usize; 2]>::from(edge.vertices) {
            let Some(used) = used_vertices.get_mut(vertex) else {
                return Err(Error::CompressedReferenceOutOfBounds {
                    entity: "vertex",
                    index: vertex,
                    len: vertices.len(),
                });
            };
            *used = true;
        }
    }

    let edge_remap = compact_items_with_stable_ids(edges, edge_stable_ids, used_edges);
    let vertex_remap = compact_items_with_stable_ids(vertices, vertex_stable_ids, &used_vertices);
    for edge in edges {
        edge.vertices.0 = vertex_remap[edge.vertices.0];
        edge.vertices.1 = vertex_remap[edge.vertices.1];
    }
    Ok(edge_remap)
}

fn retain_with_stable_ids<T>(
    items: &mut Vec<T>,
    stable_ids: &mut Option<Vec<StableId>>,
    mut keep: impl FnMut(&T) -> bool,
) {
    if let Some(ids) = stable_ids.take() {
        let mut ids = ids.into_iter();
        let mut retained_ids = Vec::with_capacity(items.len());
        items.retain(|item| {
            let stable_id = ids.next().unwrap_or(StableId::UNASSIGNED);
            if keep(item) {
                retained_ids.push(stable_id);
                true
            } else {
                false
            }
        });
        *stable_ids = Some(retained_ids);
    } else {
        items.retain(keep);
    }
}

fn remove_with_stable_id<T>(
    items: &mut Vec<T>,
    stable_ids: &mut Option<Vec<StableId>>,
    index: usize,
) -> T {
    if let Some(ids) = stable_ids {
        ids.resize(items.len(), StableId::UNASSIGNED);
        ids.remove(index);
    }
    items.remove(index)
}

pub(crate) fn validate_stable_id_counts(
    vertex_ids: Option<&[StableId]>,
    vertex_count: usize,
    edge_ids: Option<&[StableId]>,
    edge_count: usize,
    face_ids: Option<&[StableId]>,
    face_count: usize,
) -> Result<()> {
    for (entity, ids, expected) in [
        ("vertices", vertex_ids, vertex_count),
        ("edges", edge_ids, edge_count),
        ("faces", face_ids, face_count),
    ] {
        if let Some(ids) = ids
            && ids.len() != expected
        {
            return Err(Error::StableIdCountMismatch {
                entity,
                expected,
                actual: ids.len(),
            });
        }
    }
    Ok(())
}

struct CompressDirector<P, C> {
    vmap: HashMap<VertexId<P>, (usize, P)>,
    emap: HashMap<EdgeId<C>, (usize, CompressedEdge<C>)>,
    vertex_stable_ids: HashMap<VertexId<P>, (usize, StableId)>,
    edge_stable_ids: HashMap<EdgeId<C>, (usize, StableId)>,
}

impl<P: Clone, C: Clone> CompressDirector<P, C> {
    #[inline]
    fn new() -> Self {
        Self {
            vmap: HashMap::default(),
            emap: HashMap::default(),
            vertex_stable_ids: HashMap::default(),
            edge_stable_ids: HashMap::default(),
        }
    }
    #[inline]
    fn get_vid(&mut self, vertex: &Vertex<P>) -> usize {
        let id = self.vmap.len();
        let vid = vertex.id();
        self.vertex_stable_ids
            .entry(vid)
            .or_insert_with(|| (id, vertex.stable_id()));
        self.vmap
            .entry(vid)
            .or_insert_with(|| (id, vertex.point()))
            .0
    }

    #[inline]
    fn get_eid(&mut self, edge: &Edge<P, C>) -> CompressedEdgeIndex {
        if let Some(got) = self.emap.get(&edge.id()) {
            (got.0, edge.orientation()).into()
        } else {
            let id = self.emap.len();
            let front_id = self.get_vid(edge.absolute_front());
            let back_id = self.get_vid(edge.absolute_back());
            let curve = edge.curve();
            let cedge = CompressedEdge {
                vertices: (front_id, back_id),
                curve,
            };
            self.edge_stable_ids
                .insert(edge.id(), (id, edge.stable_id()));
            self.emap.insert(edge.id(), (id, cedge));
            (id, edge.orientation()).into()
        }
    }

    #[inline]
    fn create_boundary(&mut self, boundary: &Wire<P, C>) -> Vec<CompressedEdgeIndex> {
        boundary.iter().map(|edge| self.get_eid(edge)).collect()
    }

    #[inline]
    fn create_trimmed_boundary<S, T, F>(
        &mut self,
        boundary: &Wire<P, C>,
        surface: &S,
        trim_curve: &mut F,
    ) -> Vec<CompressedEdgeUse<T>>
    where
        F: FnMut(&Edge<P, C>, &S) -> Option<T>,
    {
        boundary
            .iter()
            .map(|edge| {
                let CompressedEdgeIndex { index, orientation } = self.get_eid(edge);
                CompressedEdgeUse {
                    index,
                    orientation,
                    trim_curve: trim_curve(edge, surface),
                }
            })
            .collect()
    }

    #[inline]
    fn create_cface<S: Clone>(&mut self, face: &Face<P, C, S>) -> CompressedFace<S> {
        CompressedFace {
            boundaries: face
                .boundaries
                .iter()
                .map(|wire| self.create_boundary(wire))
                .collect(),
            orientation: face.orientation(),
            surface: face.surface(),
        }
    }

    #[inline]
    fn create_trimmed_cface<S: Clone, T, F>(
        &mut self,
        face: &Face<P, C, S>,
        trim_curve: &mut F,
    ) -> CompressedTrimmedFace<S, T>
    where
        F: FnMut(&Edge<P, C>, &S) -> Option<T>,
    {
        let surface = face.surface();
        CompressedTrimmedFace {
            boundaries: face
                .boundaries
                .iter()
                .map(|wire| self.create_trimmed_boundary(wire, &surface, trim_curve))
                .collect(),
            orientation: face.orientation(),
            surface,
        }
    }

    #[inline]
    fn map2vec<K, T>(map: HashMap<K, (usize, T)>) -> Vec<T> {
        let mut vec: Vec<_> = map.into_iter().map(|entry| entry.1).collect();
        vec.sort_by_key(|x| x.0);
        vec.into_iter().map(|x| x.1).collect()
    }

    #[inline]
    fn stable_ids(self) -> (Vec<P>, Vec<CompressedEdge<C>>, Vec<StableId>, Vec<StableId>) {
        let vertices = Self::map2vec(self.vmap);
        let edges = Self::map2vec(self.emap);
        let vsids = Self::map2vec_simple(self.vertex_stable_ids);
        let esids = Self::map2vec_simple(self.edge_stable_ids);
        (vertices, edges, vsids, esids)
    }

    #[inline]
    fn map2vec_simple<K, T>(map: HashMap<K, (usize, T)>) -> Vec<T> {
        let mut vec: Vec<_> = map.into_iter().map(|entry| entry.1).collect();
        vec.sort_by_key(|x| x.0);
        vec.into_iter().map(|x| x.1).collect()
    }
}

impl<P: Clone, C: Clone, S: Clone> Shell<P, C, S> {
    /// Compresses the shell into the serialized compressed shell.
    #[must_use]
    pub fn compress(&self) -> CompressedShell<P, C, S> {
        let mut director = CompressDirector::new();
        let face_stable_ids: Vec<StableId> = self.iter().map(super::Face::stable_id).collect();
        let mut face_closure = |face: &Face<P, C, S>| director.create_cface(face);
        let faces: Vec<_> = self.iter().map(&mut face_closure).collect();
        let (vertices, edges, vsids, esids) = director.stable_ids();
        let has_any_stable_id = vsids.iter().any(|id| id.is_assigned())
            || esids.iter().any(|id| id.is_assigned())
            || face_stable_ids.iter().any(|id| id.is_assigned());
        CompressedShell {
            vertices,
            edges,
            faces,
            vertex_stable_ids: if has_any_stable_id { Some(vsids) } else { None },
            edge_stable_ids: if has_any_stable_id { Some(esids) } else { None },
            face_stable_ids: if has_any_stable_id {
                Some(face_stable_ids)
            } else {
                None
            },
        }
    }

    /// Compresses the shell into serialized compressed topology with face-local trims.
    pub fn compress_with_face_trims<T, F>(
        &self,
        mut trim_curve: F,
    ) -> CompressedTrimmedShell<P, C, S, T>
    where
        F: FnMut(&Edge<P, C>, &S) -> Option<T>,
    {
        let mut director = CompressDirector::new();
        let face_stable_ids: Vec<StableId> = self.iter().map(super::Face::stable_id).collect();
        let mut face_closure =
            |face: &Face<P, C, S>| director.create_trimmed_cface(face, &mut trim_curve);
        let faces = self.iter().map(&mut face_closure).collect();
        let (vertices, edges, vsids, esids) = director.stable_ids();
        let has_any_stable_id = vsids.iter().any(|id| id.is_assigned())
            || esids.iter().any(|id| id.is_assigned())
            || face_stable_ids.iter().any(|id| id.is_assigned());
        CompressedTrimmedShell {
            vertices,
            edges,
            faces,
            vertex_stable_ids: if has_any_stable_id { Some(vsids) } else { None },
            edge_stable_ids: if has_any_stable_id { Some(esids) } else { None },
            face_stable_ids: if has_any_stable_id {
                Some(face_stable_ids)
            } else {
                None
            },
        }
    }

    /// Compresses the shell while preserving exact face-local trims carried by the edge curves.
    #[must_use]
    pub fn compress_with_exact_face_trims(
        &self,
    ) -> CompressedTrimmedShell<P, C, S, <C as ExactParameterBoundary2D<S>>::BoundaryCurve>
    where C: ExactParameterBoundary2D<S> {
        self.compress_with_face_trims(|edge, surface| {
            edge.with_curve(|curve| curve.exact_parameter_boundary_2d(surface))
        })
    }
}

impl<P, C, S> Shell<P, C, S> {
    /// Extracts the serialized compressed shell into the shell.
    ///
    /// # Errors
    /// Returns [`Error::NotSimpleWire`](crate::errors::Error::NotSimpleWire) if any boundary wire has repeated vertices.
    /// STEP files from other CAD systems often produce such wires. To handle this,
    /// apply `SplitClosedEdgesAndFaces` or `SplitClosedEdgesAndFaces`
    /// (from `monstertruck-healing`) to the `CompressedShell` before calling `extract`.
    /// The convenience function `monstertruck_healing::extract_healed` does both steps.
    pub fn extract(cshell: CompressedShell<P, C, S>) -> Result<Self> {
        validate_stable_id_counts(
            cshell.vertex_stable_ids.as_deref(),
            cshell.vertices.len(),
            cshell.edge_stable_ids.as_deref(),
            cshell.edges.len(),
            cshell.face_stable_ids.as_deref(),
            cshell.faces.len(),
        )?;
        let CompressedShell {
            vertices,
            edges,
            faces,
            vertex_stable_ids,
            edge_stable_ids,
            face_stable_ids,
        } = cshell;
        let mut vertices: Vec<_> = vertices.into_iter().map(Vertex::new).collect();
        if let Some(vsids) = vertex_stable_ids {
            for (v, sid) in vertices.iter_mut().zip(vsids) {
                v.set_stable_id(sid);
            }
        }
        let mut edges = edges
            .into_iter()
            .map(move |edge| edge.create_edge(&vertices))
            .collect::<Result<Vec<_>>>()?;
        if let Some(esids) = edge_stable_ids {
            for (e, sid) in edges.iter_mut().zip(esids) {
                e.set_stable_id(sid);
            }
        }
        let shell: Self = faces
            .into_iter()
            .map(move |face| face.create_face(&edges))
            .collect::<Result<Self>>()?;
        if let Some(fsids) = face_stable_ids {
            let mut face_list: Vec<Face<P, C, S>> = shell.into_iter().collect();
            for (f, sid) in face_list.iter_mut().zip(fsids) {
                f.set_stable_id(sid);
            }
            Ok(Self::from(face_list))
        } else {
            Ok(shell)
        }
    }
}

impl<S, T> CompressedTrimmedFace<S, T> {
    /// Drops face-local trim data and returns a plain `CompressedFace`.
    pub fn erase_trims(self) -> CompressedFace<S> {
        CompressedFace {
            boundaries: self
                .boundaries
                .into_iter()
                .map(|wire| {
                    wire.into_iter()
                        .map(
                            |CompressedEdgeUse {
                                 index, orientation, ..
                             }| {
                                CompressedEdgeIndex { index, orientation }
                            },
                        )
                        .collect()
                })
                .collect(),
            orientation: self.orientation,
            surface: self.surface,
        }
    }

    /// Clones into a plain `CompressedFace` without cloning trim data.
    pub fn cloned_without_trims(&self) -> CompressedFace<S>
    where S: Clone {
        CompressedFace {
            boundaries: self
                .boundaries
                .iter()
                .map(|wire| {
                    wire.iter()
                        .map(
                            |CompressedEdgeUse {
                                 index, orientation, ..
                             }| CompressedEdgeIndex {
                                index: *index,
                                orientation: *orientation,
                            },
                        )
                        .collect()
                })
                .collect(),
            orientation: self.orientation,
            surface: self.surface.clone(),
        }
    }
}

/// Samples per parameter direction taken from each edge curve when computing
/// a shell's bounding box.
const CURVE_BOUNDS_SAMPLES: u32 = 16;

/// Samples per parameter direction taken from each face surface when a shell
/// carries no vertices or edges at all. Only a fallback, so a coarse grid is
/// plenty.
const SURFACE_BOUNDS_SAMPLES: u32 = 8;

impl<P, C, S, T> CompressedTrimmedShell<P, C, S, T>
where
    P: Bounded,
    C: ParametricCurve<Point = P>,
    S: ParametricSurface<Point = P>,
{
    /// Bounding box of the shell's geometry.
    ///
    /// Vertices and sampled edge curves define the box. Unbounded curves are
    /// skipped rather than panicking.
    #[must_use]
    pub fn bounding_box(&self) -> BoundingBox<P> { self.bounding_box_with(|point| point) }

    /// Bounding box of the shell's geometry, with every sampled point passed
    /// through `map` first.
    ///
    /// Use this when the shell's geometry is not already in the space you
    /// want the box in -- an assembly placement, say. Mapping each sample is
    /// exact, where transforming the corners of an axis-aligned box computed
    /// beforehand inflates it by up to `sqrt(3)` under rotation.
    ///
    /// Points come from vertices and sampled edge curves. Unbounded curves
    /// are skipped rather than panicking. If neither yields anything the box
    /// falls back to sampling the faces' surfaces; see [`Self::bounding_box`].
    pub fn bounding_box_with(&self, mut map: impl FnMut(P) -> P) -> BoundingBox<P> {
        let mut bdd_box = BoundingBox::new();
        for point in &self.vertices {
            bdd_box.push(map(*point));
        }
        for edge in &self.edges {
            let Some((t0, t1)) = edge.curve.try_range_tuple() else {
                continue;
            };
            for i in 0..=CURVE_BOUNDS_SAMPLES {
                let t = t0 + (t1 - t0) * f64::from(i) / f64::from(CURVE_BOUNDS_SAMPLES);
                bdd_box.push(map(edge.curve.evaluate(t)));
            }
        }
        if !bdd_box.is_empty() {
            return bdd_box;
        }

        // No vertices and no edges at all: every face is untrimmed. Exporters
        // write such faces as a surface bounded by a degenerate vertex loop,
        // which carries no cartesian geometry, so the pass above finds
        // nothing and callers deriving a tolerance from the box would get
        // zero. With no trim curves the surface's parameter rectangle *is*
        // the face, so sampling it is correct. Reached only when the topology
        // pass came up empty, so trimmed shells are never inflated to their
        // surfaces' parameter rectangles.
        let samples = f64::from(SURFACE_BOUNDS_SAMPLES);
        for face in &self.faces {
            let (Some((u0, u1)), Some((v0, v1))) = face.surface.try_range_tuple() else {
                continue;
            };
            for i in 0..=SURFACE_BOUNDS_SAMPLES {
                let u = u0 + (u1 - u0) * f64::from(i) / samples;
                for j in 0..=SURFACE_BOUNDS_SAMPLES {
                    let v = v0 + (v1 - v0) * f64::from(j) / samples;
                    bdd_box.push(map(face.surface.evaluate(u, v)));
                }
            }
        }

        bdd_box
    }
}

impl<C, S, T> CompressedTrimmedShell<Point3, C, S, T>
where
    C: ParametricCurve<Point = Point3>,
    S: ParametricSurface<Point = Point3>,
{
    /// Absolute tessellation tolerance as a fraction of the shell's own size.
    ///
    /// Tessellation tolerances are absolute lengths, so a fixed value means
    /// something different for a 1 mm part than for a 10 m one. `factor` is a
    /// fraction of the shell's bounding-box diameter -- `1.0e-3` is a good
    /// default -- and the result is clamped up to [`KERNEL_EPSILON`] so a
    /// degenerate or single-point shell can never ask for an infinitely fine
    /// mesh.
    ///
    /// ```no_run
    /// # use monstertruck_topology::compress::CompressedTrimmedShell;
    /// # fn demo<C, S, T>(shell: &CompressedTrimmedShell<monstertruck_core::cgmath64::Point3, C, S, T>)
    /// # where C: monstertruck_traits::ParametricCurve<Point = monstertruck_core::cgmath64::Point3>,
    /// #       S: monstertruck_traits::ParametricSurface<Point = monstertruck_core::cgmath64::Point3> {
    /// let tolerance = shell.relative_tolerance(1.0e-3);
    /// # let _ = tolerance;
    /// # }
    /// ```
    #[must_use]
    pub fn relative_tolerance(&self, factor: f64) -> f64 {
        let tolerance = self.bounding_box().diameter() * factor;
        if tolerance.is_finite() {
            f64::max(tolerance, KERNEL_EPSILON)
        } else {
            KERNEL_EPSILON
        }
    }
}

impl<P, C, S, T> CompressedTrimmedShell<P, C, S, T> {
    /// Removes edges and vertices unused by any face boundary, remapping indices and stable IDs.
    ///
    /// The operation validates all stable-ID sidecar cardinalities and compressed references before
    /// mutating the shell. Face order and face stable IDs are unchanged.
    ///
    /// # Errors
    /// Returns a typed topology error when a stable-ID sidecar is malformed or a face/edge
    /// references an out-of-range edge/vertex index.
    pub fn compact_unused_topology(&mut self) -> Result<()> {
        validate_stable_id_counts(
            self.vertex_stable_ids.as_deref(),
            self.vertices.len(),
            self.edge_stable_ids.as_deref(),
            self.edges.len(),
            self.face_stable_ids.as_deref(),
            self.faces.len(),
        )?;
        let mut used_edges = vec![false; self.edges.len()];
        for edge_use in self
            .faces
            .iter()
            .flat_map(|face| face.boundaries.iter())
            .flatten()
        {
            let Some(used) = used_edges.get_mut(edge_use.index) else {
                return Err(Error::CompressedReferenceOutOfBounds {
                    entity: "edge",
                    index: edge_use.index,
                    len: self.edges.len(),
                });
            };
            *used = true;
        }
        let edge_remap = compact_edge_vertex_storage(
            &mut self.vertices,
            &mut self.edges,
            &mut self.vertex_stable_ids,
            &mut self.edge_stable_ids,
            &used_edges,
        )?;
        for edge_use in self
            .faces
            .iter_mut()
            .flat_map(|face| face.boundaries.iter_mut())
            .flatten()
        {
            edge_use.index = edge_remap[edge_use.index];
        }
        Ok(())
    }

    /// Discards all stable topology identity carried by this shell.
    ///
    /// Use this before arbitrary topology surgery that cannot preserve or remap the
    /// existing identity domain. Geometry and connectivity are unchanged.
    pub fn discard_topology_identity(&mut self) {
        self.vertex_stable_ids = None;
        self.edge_stable_ids = None;
        self.face_stable_ids = None;
    }

    /// Retains only faces accepted by `keep`, compacting the parallel face stable-ID sidecar.
    pub fn retain_faces(&mut self, keep: impl FnMut(&CompressedTrimmedFace<S, T>) -> bool) {
        retain_with_stable_ids(&mut self.faces, &mut self.face_stable_ids, keep);
    }

    /// Removes one face and the matching face stable-ID entry when present.
    ///
    /// # Panics
    /// Panics when `index` is outside the face array, matching [`Vec::remove`].
    pub fn remove_face(&mut self, index: usize) -> CompressedTrimmedFace<S, T> {
        remove_with_stable_id(&mut self.faces, &mut self.face_stable_ids, index)
    }

    /// Drops all face-local trim data and returns a plain `CompressedShell`.
    pub fn erase_trims(self) -> CompressedShell<P, C, S> {
        CompressedShell {
            vertices: self.vertices,
            edges: self.edges,
            faces: self
                .faces
                .into_iter()
                .map(CompressedTrimmedFace::erase_trims)
                .collect(),
            vertex_stable_ids: self.vertex_stable_ids,
            edge_stable_ids: self.edge_stable_ids,
            face_stable_ids: self.face_stable_ids,
        }
    }

    /// Clones into a plain `CompressedShell` without cloning trim data.
    pub fn cloned_without_trims(&self) -> CompressedShell<P, C, S>
    where
        P: Clone,
        C: Clone,
        S: Clone, {
        CompressedShell {
            vertices: self.vertices.clone(),
            edges: self.edges.clone(),
            faces: self
                .faces
                .iter()
                .map(CompressedTrimmedFace::cloned_without_trims)
                .collect(),
            vertex_stable_ids: self.vertex_stable_ids.clone(),
            edge_stable_ids: self.edge_stable_ids.clone(),
            face_stable_ids: self.face_stable_ids.clone(),
        }
    }
}

impl<P, C, S, T> CompressedTrimmedSolid<P, C, S, T> {
    /// Resolve an assigned stable topology ID to one unique compressed location.
    ///
    /// # Errors
    /// Returns [`Error::StableIdCountMismatch`] when the requested stable-ID sidecar is malformed.
    pub fn location_by_stable_id(&self, target: StableId) -> Result<StableIdLocation> {
        find_stable_id_location(&self.boundaries, target)
    }

    /// Discards the solid's entire stable topology identity domain.
    ///
    /// This is the correct operation before arbitrary topology composition or
    /// surgery where the existing sidecars cannot be proven valid. The allocator
    /// and attributes are discarded too because both are keyed by that identity domain.
    pub fn discard_topology_identity(&mut self) {
        for shell in &mut self.boundaries {
            shell.discard_topology_identity();
        }
        self.id_allocator = None;
        self.attributes = None;
    }

    /// Drops all face-local trim data and returns a plain `CompressedSolid`.
    pub fn erase_trims(self) -> CompressedSolid<P, C, S> {
        CompressedSolid {
            boundaries: self
                .boundaries
                .into_iter()
                .map(CompressedTrimmedShell::erase_trims)
                .collect(),
            id_allocator: self.id_allocator,
            attributes: self.attributes,
        }
    }

    /// Clones into a plain `CompressedSolid` without cloning trim data.
    pub fn cloned_without_trims(&self) -> CompressedSolid<P, C, S>
    where
        P: Clone,
        C: Clone,
        S: Clone, {
        CompressedSolid {
            boundaries: self
                .boundaries
                .iter()
                .map(CompressedTrimmedShell::cloned_without_trims)
                .collect(),
            id_allocator: self.id_allocator.clone(),
            attributes: self.attributes.clone(),
        }
    }
}

impl<P, C, S> CompressedShell<P, C, S> {
    /// Removes edges and vertices unused by any face boundary, remapping indices and stable IDs.
    ///
    /// The operation validates all stable-ID sidecar cardinalities and compressed references before
    /// mutating the shell. Face order and face stable IDs are unchanged.
    ///
    /// # Errors
    /// Returns a typed topology error when a stable-ID sidecar is malformed or a face/edge
    /// references an out-of-range edge/vertex index.
    pub fn compact_unused_topology(&mut self) -> Result<()> {
        validate_stable_id_counts(
            self.vertex_stable_ids.as_deref(),
            self.vertices.len(),
            self.edge_stable_ids.as_deref(),
            self.edges.len(),
            self.face_stable_ids.as_deref(),
            self.faces.len(),
        )?;
        let mut used_edges = vec![false; self.edges.len()];
        for edge_use in self
            .faces
            .iter()
            .flat_map(|face| face.boundaries.iter())
            .flatten()
        {
            let Some(used) = used_edges.get_mut(edge_use.index) else {
                return Err(Error::CompressedReferenceOutOfBounds {
                    entity: "edge",
                    index: edge_use.index,
                    len: self.edges.len(),
                });
            };
            *used = true;
        }
        let edge_remap = compact_edge_vertex_storage(
            &mut self.vertices,
            &mut self.edges,
            &mut self.vertex_stable_ids,
            &mut self.edge_stable_ids,
            &used_edges,
        )?;
        for edge_use in self
            .faces
            .iter_mut()
            .flat_map(|face| face.boundaries.iter_mut())
            .flatten()
        {
            edge_use.index = edge_remap[edge_use.index];
        }
        Ok(())
    }

    /// Discards all stable topology identity carried by this shell.
    pub fn discard_topology_identity(&mut self) {
        self.vertex_stable_ids = None;
        self.edge_stable_ids = None;
        self.face_stable_ids = None;
    }

    /// Retains only faces accepted by `keep`, compacting the parallel face stable-ID sidecar.
    pub fn retain_faces(&mut self, keep: impl FnMut(&CompressedFace<S>) -> bool) {
        retain_with_stable_ids(&mut self.faces, &mut self.face_stable_ids, keep);
    }

    /// Removes one face and the matching face stable-ID entry when present.
    ///
    /// # Panics
    /// Panics when `index` is outside the face array, matching [`Vec::remove`].
    pub fn remove_face(&mut self, index: usize) -> CompressedFace<S> {
        remove_with_stable_id(&mut self.faces, &mut self.face_stable_ids, index)
    }
}

impl<P: Clone, C: Clone, S: Clone> CompressedShell<P, C, S> {
    /// Maps the shared edge curves into a new [`CompressedShell`].
    pub fn map_curves<D, F>(&self, mut map_curve: F) -> CompressedShell<P, D, S>
    where F: FnMut(&C) -> D {
        CompressedShell {
            vertices: self.vertices.clone(),
            edges: self
                .edges
                .iter()
                .map(|edge| CompressedEdge {
                    vertices: edge.vertices,
                    curve: map_curve(&edge.curve),
                })
                .collect(),
            faces: self.faces.clone(),
            vertex_stable_ids: self.vertex_stable_ids.clone(),
            edge_stable_ids: self.edge_stable_ids.clone(),
            face_stable_ids: self.face_stable_ids.clone(),
        }
    }

    /// Attaches face-local trims to each edge-use in this compressed shell.
    pub fn with_face_trims<T, F>(&self, mut trim_curve: F) -> CompressedTrimmedShell<P, C, S, T>
    where F: FnMut(&C, &S) -> Option<T> {
        CompressedTrimmedShell {
            vertices: self.vertices.clone(),
            edges: self.edges.clone(),
            faces: self
                .faces
                .iter()
                .map(|face| CompressedTrimmedFace {
                    boundaries: face
                        .boundaries
                        .iter()
                        .map(|wire| {
                            wire.iter()
                                .map(|CompressedEdgeIndex { index, orientation }| {
                                    CompressedEdgeUse {
                                        index: *index,
                                        orientation: *orientation,
                                        trim_curve: trim_curve(
                                            &self.edges[*index].curve,
                                            &face.surface,
                                        ),
                                    }
                                })
                                .collect()
                        })
                        .collect(),
                    orientation: face.orientation,
                    surface: face.surface.clone(),
                })
                .collect(),
            vertex_stable_ids: self.vertex_stable_ids.clone(),
            edge_stable_ids: self.edge_stable_ids.clone(),
            face_stable_ids: self.face_stable_ids.clone(),
        }
    }

    /// Attaches exact face-local trims to each edge-use in this compressed shell.
    #[must_use]
    pub fn with_exact_face_trims(
        &self,
    ) -> CompressedTrimmedShell<P, C, S, <C as ExactParameterBoundary2D<S>>::BoundaryCurve>
    where C: ExactParameterBoundary2D<S> {
        self.with_face_trims(
            monstertruck_traits::ExactParameterBoundary2D::exact_parameter_boundary_2d,
        )
    }
}

impl<P, C, S> CompressedSolid<P, C, S> {
    /// Resolve an assigned stable topology ID to one unique compressed location.
    ///
    /// # Errors
    /// Returns [`Error::StableIdCountMismatch`] when the requested stable-ID sidecar is malformed.
    pub fn location_by_stable_id(&self, target: StableId) -> Result<StableIdLocation> {
        find_stable_id_location(&self.boundaries, target)
    }

    /// Discards the solid's entire stable topology identity domain.
    ///
    /// The allocator and attributes are discarded with the shell sidecars because
    /// they are meaningful only within that same stable-ID namespace.
    pub fn discard_topology_identity(&mut self) {
        for shell in &mut self.boundaries {
            shell.discard_topology_identity();
        }
        self.id_allocator = None;
        self.attributes = None;
    }
}

impl<P: Clone, C: Clone, S: Clone> CompressedSolid<P, C, S> {
    /// Maps the shared edge curves into a new [`CompressedSolid`].
    pub fn map_curves<D, F>(&self, mut map_curve: F) -> CompressedSolid<P, D, S>
    where F: FnMut(&C) -> D {
        CompressedSolid {
            boundaries: self
                .boundaries
                .iter()
                .map(|shell| shell.map_curves(&mut map_curve))
                .collect(),
            id_allocator: self.id_allocator.clone(),
            attributes: self.attributes.clone(),
        }
    }

    /// Attaches face-local trims to each edge-use in this compressed solid.
    pub fn with_face_trims<T, F>(&self, mut trim_curve: F) -> CompressedTrimmedSolid<P, C, S, T>
    where F: FnMut(&C, &S) -> Option<T> {
        CompressedTrimmedSolid {
            boundaries: self
                .boundaries
                .iter()
                .map(|shell| shell.with_face_trims(&mut trim_curve))
                .collect(),
            id_allocator: self.id_allocator.clone(),
            attributes: self.attributes.clone(),
        }
    }

    /// Attaches exact face-local trims to each edge-use in this compressed solid.
    #[must_use]
    pub fn with_exact_face_trims(
        &self,
    ) -> CompressedTrimmedSolid<P, C, S, <C as ExactParameterBoundary2D<S>>::BoundaryCurve>
    where C: ExactParameterBoundary2D<S> {
        self.with_face_trims(
            monstertruck_traits::ExactParameterBoundary2D::exact_parameter_boundary_2d,
        )
    }
}

impl<P: Clone, C: Clone, S: Clone> Solid<P, C, S> {
    /// Compresses the solid into the serialized compressed solid.
    #[must_use]
    pub fn compress(&self) -> CompressedSolid<P, C, S> {
        CompressedSolid {
            boundaries: self
                .boundaries()
                .iter()
                .map(super::Shell::compress)
                .collect(),
            id_allocator: Some(self.id_allocator.clone()),
            attributes: if self.attributes.is_empty() {
                None
            } else {
                Some(self.attributes.clone())
            },
        }
    }

    /// Compresses the solid into serialized compressed topology with face-local trims.
    pub fn compress_with_face_trims<T, F>(
        &self,
        mut trim_curve: F,
    ) -> CompressedTrimmedSolid<P, C, S, T>
    where
        F: FnMut(&Edge<P, C>, &S) -> Option<T>,
    {
        CompressedTrimmedSolid {
            boundaries: self
                .boundaries()
                .iter()
                .map(|shell| shell.compress_with_face_trims(&mut trim_curve))
                .collect(),
            id_allocator: Some(self.id_allocator.clone()),
            attributes: if self.attributes.is_empty() {
                None
            } else {
                Some(self.attributes.clone())
            },
        }
    }

    /// Compresses the solid while preserving exact face-local trims carried by the edge curves.
    #[must_use]
    pub fn compress_with_exact_face_trims(
        &self,
    ) -> CompressedTrimmedSolid<P, C, S, <C as ExactParameterBoundary2D<S>>::BoundaryCurve>
    where C: ExactParameterBoundary2D<S> {
        self.compress_with_face_trims(|edge, surface| {
            edge.curve().exact_parameter_boundary_2d(surface)
        })
    }

    /// Extracts the serialized compressed solid into the solid.
    ///
    /// # Errors
    /// Returns an error when a compressed shell is invalid or the extracted
    /// shells do not form valid solid boundaries.
    pub fn extract(csolid: CompressedSolid<P, C, S>) -> Result<Self> {
        let id_allocator = csolid.id_allocator.unwrap_or_default();
        let attributes = csolid.attributes.unwrap_or_default();
        let shells: Result<Vec<Shell<P, C, S>>> =
            csolid.boundaries.into_iter().map(Shell::extract).collect();
        let mut solid = Self::new(shells?)?;
        solid.id_allocator = id_allocator;
        solid.attributes = attributes;
        Ok(solid)
    }
}

// ---------------------------------------------------------------------------
// Deterministic topology hashing.
// ---------------------------------------------------------------------------

/// Options controlling which data is included in a [`Solid`] hash.
#[derive(Clone, Copy, Debug)]
pub struct SolidHashOptions {
    /// Include analytic geometry content (point, curve, surface data).
    pub include_geometry: bool,
    /// Include sparse per-element attributes.
    pub include_attributes: bool,
    /// Include stable IDs in the hash.
    pub include_stable_ids: bool,
}

impl Default for SolidHashOptions {
    fn default() -> Self {
        Self {
            include_geometry: false,
            include_attributes: false,
            include_stable_ids: true,
        }
    }
}

/// Hash the structural connectivity of a compressed shell, ignoring geometry.
fn hash_shell_topology<P, C, S>(shell: &CompressedShell<P, C, S>, hasher: &mut ContentHasher) {
    // Vertex count.
    hasher.write_usize(shell.vertices.len());
    // Edge connectivity (vertex indices only).
    hasher.write_usize(shell.edges.len());
    shell.edges.iter().for_each(|e| {
        e.vertices.0.content_hash(hasher);
        e.vertices.1.content_hash(hasher);
    });
    // Face boundary structure.
    hasher.write_usize(shell.faces.len());
    shell.faces.iter().for_each(|f| {
        f.orientation.content_hash(hasher);
        hasher.write_usize(f.boundaries.len());
        for wire in &f.boundaries {
            hasher.write_usize(wire.len());
            for edge_use in wire {
                edge_use.index.content_hash(hasher);
                edge_use.orientation.content_hash(hasher);
            }
        }
    });
}

/// Hash the stable IDs carried by a compressed shell.
fn hash_shell_stable_ids<P, C, S>(shell: &CompressedShell<P, C, S>, hasher: &mut ContentHasher) {
    shell.vertex_stable_ids.content_hash(hasher);
    shell.edge_stable_ids.content_hash(hasher);
    shell.face_stable_ids.content_hash(hasher);
}

/// Hash the geometry content of a compressed shell.
fn hash_shell_geometry<P, C, S>(shell: &CompressedShell<P, C, S>, hasher: &mut ContentHasher)
where
    P: DeterministicContentHash,
    C: DeterministicContentHash,
    S: DeterministicContentHash, {
    shell.vertices.content_hash(hasher);
    shell
        .edges
        .iter()
        .for_each(|e| e.curve.content_hash(hasher));
    shell
        .faces
        .iter()
        .for_each(|f| f.surface.content_hash(hasher));
}

impl<P: Clone, C: Clone, S: Clone> Solid<P, C, S> {
    /// Compute a deterministic hash of the topology structure and stable IDs.
    ///
    /// This hash is independent of pointer identity, attribute edits, and
    /// geometry content. Two solids with the same connectivity and stable IDs
    /// produce the same hash.
    #[must_use]
    pub fn topology_hash(&self) -> u64 {
        let compressed = self.compress();
        let mut hasher = ContentHasher::new();
        hasher.write_usize(compressed.boundaries.len());
        compressed.boundaries.iter().for_each(|shell| {
            hash_shell_topology(shell, &mut hasher);
            hash_shell_stable_ids(shell, &mut hasher);
        });
        // Include allocator state (semantically relevant).
        if let Some(ref alloc) = compressed.id_allocator {
            hasher.write_u8(1);
            alloc.peek().content_hash(&mut hasher);
        } else {
            hasher.write_u8(0);
        }
        hasher.finish()
    }

    /// Compute a deterministic hash of topology, stable IDs, and sparse attributes.
    ///
    /// Unlike [`topology_hash`](Self::topology_hash), this changes when per-element
    /// attributes (e.g. selection stamps) are modified.
    #[must_use]
    pub fn topology_attribute_hash(&self) -> u64 {
        let compressed = self.compress();
        let mut hasher = ContentHasher::new();
        hasher.write_usize(compressed.boundaries.len());
        compressed.boundaries.iter().for_each(|shell| {
            hash_shell_topology(shell, &mut hasher);
            hash_shell_stable_ids(shell, &mut hasher);
        });
        if let Some(ref alloc) = compressed.id_allocator {
            hasher.write_u8(1);
            alloc.peek().content_hash(&mut hasher);
        } else {
            hasher.write_u8(0);
        }
        // Include attributes.
        self.attributes.content_hash(&mut hasher);
        hasher.finish()
    }
}

impl<P, C, S> Solid<P, C, S>
where
    P: Clone + DeterministicContentHash,
    C: Clone + DeterministicContentHash,
    S: Clone + DeterministicContentHash,
{
    /// Compute a full deterministic content hash.
    ///
    /// Includes canonical topology, stable IDs, sparse attributes, and
    /// analytic geometry. Remains deterministic across runs and independent
    /// of pointer identity.
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        let compressed = self.compress();
        let mut hasher = ContentHasher::new();
        hasher.write_usize(compressed.boundaries.len());
        compressed.boundaries.iter().for_each(|shell| {
            hash_shell_topology(shell, &mut hasher);
            hash_shell_stable_ids(shell, &mut hasher);
            hash_shell_geometry(shell, &mut hasher);
        });
        if let Some(ref alloc) = compressed.id_allocator {
            hasher.write_u8(1);
            alloc.peek().content_hash(&mut hasher);
        } else {
            hasher.write_u8(0);
        }
        self.attributes.content_hash(&mut hasher);
        hasher.finish()
    }
}

// -------------------------- test --------------------------

#[test]
fn compress_extract() {
    let cube = solid::cube();
    let shell0 = &cube.boundaries()[0];
    let shell1 = Shell::extract(shell0.compress()).unwrap();
    assert!(same_topology(shell0, &shell1));
}

#[cfg(test)]
fn vmap_subroutin<P, Q>(
    v0: &Vertex<P>,
    v1: &Vertex<Q>,
    vmap: &mut HashMap<VertexId<P>, VertexId<Q>>,
) -> bool {
    let existing = vmap.get(&v0.id()).copied();
    existing.map_or_else(
        || {
            vmap.insert(v0.id(), v1.id());
            true
        },
        |got| got == v1.id(),
    )
}

#[cfg(test)]
fn emap_subroutin<P, Q, C, D>(
    edge0: &Edge<P, C>,
    edge1: &Edge<Q, D>,
    vmap: &mut HashMap<VertexId<P>, VertexId<Q>>,
    emap: &mut HashMap<EdgeId<C>, EdgeId<D>>,
) -> bool {
    let existing = emap.get(&edge0.id()).copied();
    existing.map_or_else(
        || {
            emap.insert(edge0.id(), edge1.id());
            vmap_subroutin(edge0.front(), edge1.front(), vmap)
                && vmap_subroutin(edge0.back(), edge1.back(), vmap)
        },
        |got| got == edge1.id(),
    )
}

#[cfg(test)]
fn same_topology<P, C, S, Q, D, T>(one: &Shell<P, C, S>, other: &Shell<Q, D, T>) -> bool {
    let mut vmap = HashMap::<VertexId<P>, VertexId<Q>>::default();
    let mut emap = HashMap::<EdgeId<C>, EdgeId<D>>::default();
    if one.len() != other.len() {
        return false;
    }
    for (face0, face1) in one.iter().zip(other.iter()) {
        let biters0 = face0.boundary_iters();
        let biters1 = face1.boundary_iters();
        if biters0.len() != biters1.len() {
            return false;
        }
        for (biter0, biter1) in biters0.into_iter().zip(biters1) {
            if biter0.len() != biter1.len() {
                return false;
            }
            for (edge0, edge1) in biter0.zip(biter1) {
                if !emap_subroutin(&edge0, &edge1, &mut vmap, &mut emap) {
                    return false;
                }
            }
        }
    }
    true
}

#[test]
fn extract_rejects_misaligned_stable_id_sidecars() {
    let cube = solid::cube();
    let shell = &cube.boundaries()[0];
    let compressed = shell.compress();

    let mut bad_vertices = compressed.clone();
    bad_vertices.vertex_stable_ids =
        Some(vec![StableId::UNASSIGNED; bad_vertices.vertices.len() - 1]);
    assert!(matches!(
        Shell::extract(bad_vertices),
        Err(Error::StableIdCountMismatch {
            entity: "vertices",
            ..
        })
    ));

    let mut bad_edges = compressed.clone();
    bad_edges.edge_stable_ids = Some(vec![StableId::UNASSIGNED; bad_edges.edges.len() - 1]);
    assert!(matches!(
        Shell::extract(bad_edges),
        Err(Error::StableIdCountMismatch {
            entity: "edges",
            ..
        })
    ));

    let mut bad_faces = compressed;
    bad_faces.face_stable_ids = Some(vec![StableId::UNASSIGNED; bad_faces.faces.len() - 1]);
    assert!(matches!(
        Shell::extract(bad_faces),
        Err(Error::StableIdCountMismatch {
            entity: "faces",
            ..
        })
    ));
}

#[test]
fn stable_ids_survive_compress_extract() -> std::result::Result<(), StableIdExhausted> {
    // Build a cube, take it apart, assign stable IDs to faces, rebuild.
    let cube = solid::cube();
    let mut shells = cube.into_boundaries();
    let shell = &mut shells[0];
    let mut alloc = StableIdAllocator::new();
    for face in shell.face_iter_mut() {
        face.set_stable_id(alloc.allocate()?);
    }

    let orig_sids: Vec<_> = shell.face_iter().map(|f| f.stable_id()).collect();
    let compressed = shell.compress();

    // All face stable IDs should appear in the compressed form.
    let csids = compressed.face_stable_ids.as_ref().unwrap();
    assert_eq!(&orig_sids, csids);

    // Extract and verify.
    let extracted_shell = Shell::extract(compressed).unwrap();
    let ext_sids: Vec<_> = extracted_shell.face_iter().map(|f| f.stable_id()).collect();
    assert_eq!(orig_sids, ext_sids);
    Ok(())
}

#[test]
fn trimmed_compression_preserves_topology_stable_ids()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut cube = solid::cube();
    cube.ensure_topology_stable_ids()?;
    let plain = cube.boundaries()[0].compress();
    let trimmed = plain.with_face_trims(|_, _| None::<()>);

    assert_eq!(trimmed.vertex_stable_ids, plain.vertex_stable_ids);
    assert_eq!(trimmed.edge_stable_ids, plain.edge_stable_ids);
    assert_eq!(trimmed.face_stable_ids, plain.face_stable_ids);
    assert_eq!(trimmed.clone().erase_trims(), plain);
    assert_eq!(trimmed.cloned_without_trims(), plain);
    Ok(())
}

#[test]
fn runtime_trimmed_roundtrip_preserves_topology_stable_ids()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut cube = solid::cube();
    cube.ensure_topology_stable_ids()?;
    let source = cube.boundaries()[0].compress_with_face_trims(|_, _| None::<()>);
    let runtime: crate::trimmed::TrimmedShell<_, _, _, ()> = source.clone().try_into()?;
    let roundtrip = CompressedTrimmedShell::from(&runtime);

    assert_eq!(roundtrip.vertex_stable_ids, source.vertex_stable_ids);
    assert_eq!(roundtrip.edge_stable_ids, source.edge_stable_ids);
    assert_eq!(roundtrip.face_stable_ids, source.face_stable_ids);
    Ok(())
}

#[test]
fn trimmed_extract_rejects_misaligned_stable_id_sidecars()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut cube = solid::cube();
    cube.ensure_topology_stable_ids()?;
    let mut trimmed = cube.boundaries()[0].compress_with_face_trims(|_, _| None::<()>);
    if let Some(ids) = &mut trimmed.vertex_stable_ids {
        ids.pop();
    }

    assert!(matches!(
        crate::trimmed::TrimmedShell::try_from(trimmed),
        Err(Error::StableIdCountMismatch {
            entity: "vertices",
            ..
        })
    ));
    Ok(())
}

#[test]
fn runtime_trimmed_solid_assigns_stable_ids_without_erasing_trims()
-> std::result::Result<(), StableIdExhausted> {
    let cube = solid::cube();
    let mut trimmed = cube.to_trimmed_with_face_trims(|_, _| None::<()>);
    let before = CompressedTrimmedSolid::from(&trimmed);
    assert!(before.boundaries.iter().all(|shell| {
        shell.vertex_stable_ids.is_none()
            && shell.edge_stable_ids.is_none()
            && shell.face_stable_ids.is_none()
    }));

    trimmed.ensure_topology_stable_ids()?;

    let after = CompressedTrimmedSolid::from(&trimmed);
    assert!(after.boundaries.iter().all(|shell| {
        shell
            .vertex_stable_ids
            .as_ref()
            .is_some_and(|ids| ids.iter().all(|id| id.is_assigned()))
            && shell
                .edge_stable_ids
                .as_ref()
                .is_some_and(|ids| ids.iter().all(|id| id.is_assigned()))
            && shell
                .face_stable_ids
                .as_ref()
                .is_some_and(|ids| ids.iter().all(|id| id.is_assigned()))
    }));
    Ok(())
}

#[test]
fn compressed_solid_trim_roundtrip_preserves_allocator_and_attributes()
-> std::result::Result<(), StableIdExhausted> {
    use crate::attributes::AttributeValue;

    let mut compressed = solid::cube().compress();
    let mut allocator = StableIdAllocator::new();
    let face_id = allocator.allocate()?;
    let mut attributes = SolidAttributes::new();
    attributes.faces.set(
        "marker",
        face_id,
        AttributeValue::String("trimmed".to_owned()),
    );
    compressed.id_allocator = Some(allocator);
    compressed.attributes = Some(attributes);

    let trimmed = compressed.with_face_trims(|_, _| None::<()>);
    assert_eq!(trimmed.id_allocator, compressed.id_allocator);
    assert_eq!(trimmed.attributes, compressed.attributes);
    assert_eq!(trimmed.erase_trims(), compressed);
    Ok(())
}

#[test]
fn runtime_trimmed_solid_preserves_allocator_and_attributes()
-> std::result::Result<(), Box<dyn std::error::Error>> {
    use crate::attributes::AttributeValue;

    let mut solid = solid::cube();
    solid.ensure_topology_stable_ids()?;
    let face_id = solid.alloc_id()?;
    solid.face_attributes_mut().set(
        "marker",
        face_id,
        AttributeValue::String("runtime-trimmed".to_owned()),
    );
    let allocator_peek = solid.id_allocator().peek();
    let attributes = solid.attributes().clone();

    let trimmed = solid.to_trimmed_with_face_trims(|_, _| None::<()>);
    assert_eq!(trimmed.id_allocator().peek(), allocator_peek);
    assert_eq!(trimmed.attributes(), &attributes);

    let compressed = CompressedTrimmedSolid::from(&trimmed);
    let runtime: crate::trimmed::TrimmedSolid<_, _, _, ()> = compressed.try_into()?;
    assert_eq!(runtime.id_allocator().peek(), allocator_peek);
    assert_eq!(runtime.attributes(), &attributes);

    let plain = runtime.erase_trims();
    assert_eq!(plain.id_allocator().peek(), allocator_peek);
    assert_eq!(plain.attributes(), &attributes);
    Ok(())
}

#[test]
fn attributes_survive_compress_extract() -> std::result::Result<(), StableIdExhausted> {
    use crate::attributes::AttributeValue;
    let cube = solid::cube();
    let mut alloc = StableIdAllocator::new();
    let face_id = alloc.allocate()?;

    // Build compressed solid manually to set attributes.
    let mut compressed = cube.compress();
    let mut attrs = SolidAttributes::new();
    attrs.faces.set(
        "color",
        face_id,
        AttributeValue::Color([1.0, 0.0, 0.0, 1.0]),
    );
    compressed.attributes = Some(attrs);
    compressed.id_allocator = Some(alloc);

    let extracted = Solid::extract(compressed).unwrap();
    assert_eq!(
        extracted.face_attributes().get("color", face_id),
        Some(&AttributeValue::Color([1.0, 0.0, 0.0, 1.0]))
    );
    Ok(())
}

#[test]
fn id_allocator_survives_compress_extract() {
    let mut cube = solid::cube();
    // Advance the allocator.
    let _a = cube.alloc_id();
    let _b = cube.alloc_id();
    let peek_before = cube.id_allocator().peek();

    let compressed = cube.compress();
    let extracted = Solid::extract(compressed).unwrap();
    assert_eq!(extracted.id_allocator().peek(), peek_before);
}

impl<P, C, S> Serialize for Shell<P, C, S>
where
    P: Clone + Serialize,
    C: Clone + Serialize,
    S: Clone + Serialize,
{
    fn serialize<Serializer>(
        &self,
        serializer: Serializer,
    ) -> std::result::Result<Serializer::Ok, Serializer::Error>
    where
        Serializer: serde::Serializer,
    {
        self.compress().serialize(serializer)
    }
}

impl<'de, P, C, S> Deserialize<'de> for Shell<P, C, S>
where
    P: Clone + Deserialize<'de>,
    C: Clone + Deserialize<'de>,
    S: Clone + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        use serde::de::Error;
        let compressed = CompressedShell::<P, C, S>::deserialize(deserializer)?;
        Self::extract(compressed).map_err(D::Error::custom)
    }
}

impl<P, C, S> Serialize for Solid<P, C, S>
where
    P: Clone + Serialize,
    C: Clone + Serialize,
    S: Clone + Serialize,
{
    fn serialize<Serializer>(
        &self,
        serializer: Serializer,
    ) -> std::result::Result<Serializer::Ok, Serializer::Error>
    where
        Serializer: serde::Serializer,
    {
        self.compress().serialize(serializer)
    }
}

impl<'de, P, C, S> Deserialize<'de> for Solid<P, C, S>
where
    P: Clone + Deserialize<'de>,
    C: Clone + Deserialize<'de>,
    S: Clone + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        use serde::de::Error;
        let compressed = CompressedSolid::<P, C, S>::deserialize(deserializer)?;
        Self::extract(compressed).map_err(D::Error::custom)
    }
}

impl<P, C, S> Serialize for Face<P, C, S>
where
    P: Clone + Serialize,
    C: Clone + Serialize,
    S: Clone + Serialize,
{
    fn serialize<Serializer>(
        &self,
        serializer: Serializer,
    ) -> std::result::Result<Serializer::Ok, Serializer::Error>
    where
        Serializer: serde::Serializer,
    {
        Shell::from(vec![self.clone()]).serialize(serializer)
    }
}

impl<'de, P, C, S> Deserialize<'de> for Face<P, C, S>
where
    P: Clone + Deserialize<'de>,
    C: Clone + Deserialize<'de>,
    S: Clone + Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        Shell::deserialize(deserializer).and_then(|mut shell| {
            shell.pop().ok_or_else(|| {
                <D::Error as serde::de::Error>::custom("serialized face contained no face")
            })
        })
    }
}

#[cfg(test)]
mod tests;
