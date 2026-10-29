#[cfg(test)]
use crate::shell::ShellCondition;
use crate::{
    BoundedCurve, Concat, CurveParameter, Cut, Debug, DebugDisplay, Edge, EdgeId,
    ElementAttributes, Face, FaceId, Formatter, IncludeCurve, Invertible, ParameterTransform,
    Result, SearchParameter, Shell, Solid, SolidAttributes, SolidDisplayFormat, StableId,
    StableIdAllocator, StableIdExhausted, Tolerance, Vertex, VertexId,
};
use rustc_hash::FxHashMap as HashMap;
use std::{hash::Hash, vec::Vec};

fn unassigned_shared_count<K>(uses: impl IntoIterator<Item = (K, StableId)>) -> usize
where K: Eq + Hash {
    let mut assignments = HashMap::<K, bool>::default();
    for (key, stable_id) in uses {
        let assigned = stable_id.is_assigned();
        assignments
            .entry(key)
            .and_modify(|any_assigned| *any_assigned |= assigned)
            .or_insert(assigned);
    }
    assignments.values().filter(|&&assigned| !assigned).count()
}

fn assign_vertex_stable_id<P>(
    vertex: &mut Vertex<P>,
    ids: &mut HashMap<VertexId<P>, StableId>,
    id_allocator: &mut StableIdAllocator,
) -> std::result::Result<(), StableIdExhausted> {
    let vertex_id = vertex.id();
    let stable_id = if let Some(stable_id) = ids.get(&vertex_id).copied() {
        stable_id
    } else {
        let stable_id = id_allocator.allocate()?;
        ids.insert(vertex_id, stable_id);
        stable_id
    };
    vertex.set_stable_id(stable_id);
    Ok(())
}

fn collect_assigned_vertex_stable_id<P>(
    vertex: &Vertex<P>,
    ids: &mut HashMap<VertexId<P>, StableId>,
) {
    if vertex.stable_id().is_assigned() && !ids.contains_key(&vertex.id()) {
        ids.insert(vertex.id(), vertex.stable_id());
    }
}

fn collect_assigned_edge_stable_id<P, C>(
    edge: &Edge<P, C>,
    ids: &mut HashMap<EdgeId<C>, StableId>,
) {
    if edge.stable_id().is_assigned() && !ids.contains_key(&edge.id()) {
        ids.insert(edge.id(), edge.stable_id());
    }
}

impl<P, C, S> Solid<P, C, S> {
    /// Creates a solid whose boundaries are non-empty, connected, closed manifolds.
    ///
    /// # Errors
    /// Returns the topology validation error reported by [`Shell::check_solid_boundary`]
    /// for the first invalid boundary shell.
    #[inline]
    pub fn new(boundaries: Vec<Shell<P, C, S>>) -> Result<Self> {
        for shell in &boundaries {
            shell.check_solid_boundary()?;
        }
        Ok(Self::new_unchecked(boundaries))
    }

    /// Create without validation (for performance-critical code).
    #[inline]
    #[must_use]
    pub fn new_unchecked(boundaries: Vec<Shell<P, C, S>>) -> Self {
        Self {
            boundaries,
            id_allocator: StableIdAllocator::new(),
            attributes: SolidAttributes::new(),
        }
    }

    /// Allocates a fresh [`StableId`] from this solid's allocator.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] if the identifier space is exhausted.
    #[inline]
    pub fn alloc_id(&mut self) -> std::result::Result<StableId, StableIdExhausted> {
        self.id_allocator.allocate()
    }

    /// Returns a reference to this solid's [`StableIdAllocator`].
    #[inline]
    #[must_use]
    pub const fn id_allocator(&self) -> &StableIdAllocator { &self.id_allocator }

    /// Returns a mutable reference to this solid's [`StableIdAllocator`].
    #[inline]
    pub const fn id_allocator_mut(&mut self) -> &mut StableIdAllocator { &mut self.id_allocator }

    /// Returns a reference to this solid's [`SolidAttributes`].
    #[inline]
    #[must_use]
    pub const fn attributes(&self) -> &SolidAttributes { &self.attributes }

    /// Returns a mutable reference to this solid's [`SolidAttributes`].
    #[inline]
    pub const fn attributes_mut(&mut self) -> &mut SolidAttributes { &mut self.attributes }

    /// Returns a reference to the face [`ElementAttributes`].
    #[inline]
    #[must_use]
    pub const fn face_attributes(&self) -> &ElementAttributes { &self.attributes.faces }

    /// Returns a mutable reference to the face [`ElementAttributes`].
    #[inline]
    pub const fn face_attributes_mut(&mut self) -> &mut ElementAttributes {
        &mut self.attributes.faces
    }

    /// Returns a reference to the edge [`ElementAttributes`].
    #[inline]
    #[must_use]
    pub const fn edge_attributes(&self) -> &ElementAttributes { &self.attributes.edges }

    /// Returns a mutable reference to the edge [`ElementAttributes`].
    #[inline]
    pub const fn edge_attributes_mut(&mut self) -> &mut ElementAttributes {
        &mut self.attributes.edges
    }

    /// Returns a reference to the vertex [`ElementAttributes`].
    #[inline]
    #[must_use]
    pub const fn vertex_attributes(&self) -> &ElementAttributes { &self.attributes.vertices }

    /// Returns a mutable reference to the vertex [`ElementAttributes`].
    #[inline]
    pub const fn vertex_attributes_mut(&mut self) -> &mut ElementAttributes {
        &mut self.attributes.vertices
    }

    /// Returns the reference of boundary shells
    #[inline]
    #[must_use]
    pub const fn boundaries(&self) -> &Vec<Shell<P, C, S>> { &self.boundaries }
    /// Returns the boundary shells
    #[inline]
    #[must_use]
    pub fn into_boundaries(self) -> Vec<Shell<P, C, S>> { self.boundaries }

    /// The valid EMPTY solid: zero boundary shells, hence zero faces and zero
    /// volume -- the geometric empty set as a first-class B-rep value.
    ///
    /// This is a well-formed [`Solid`]: [`Solid::new`] accepts an empty
    /// boundary list (the per-shell validity loop is vacuously satisfied), so
    /// the empty set is representable without any `unchecked` escape. Boolean
    /// operations return it for results that are geometrically, provably empty
    /// (e.g. the intersection of disjoint solids), so a consumer receives a
    /// closed zero-volume value rather than an error signal.
    #[inline]
    #[must_use]
    pub fn empty() -> Self { Self::new_unchecked(Vec::new()) }

    /// True when this solid is the EMPTY set: it has no boundary shells (hence
    /// no faces and zero volume). See [`Solid::empty`].
    #[inline]
    #[must_use]
    pub const fn is_empty(&self) -> bool { self.boundaries.is_empty() }

    /// Returns an iterator over the faces.
    #[inline]
    pub fn face_iter(&self) -> impl Iterator<Item = &Face<P, C, S>> {
        self.boundaries.iter().flatten()
    }

    /// Returns an iterator over the edges.
    #[inline]
    pub fn edge_iter(&self) -> impl Iterator<Item = Edge<P, C>> + '_ {
        self.face_iter().flat_map(Face::boundaries).flatten()
    }

    /// Returns an iterator over the vertices.
    #[inline]
    pub fn vertex_iter(&self) -> impl Iterator<Item = Vertex<P>> + '_ {
        self.edge_iter().map(|edge| edge.front().clone())
    }

    /// Returns every stored vertex use, including both endpoints of every edge use.
    ///
    /// Unlike [`Solid::vertex_iter`], this deliberately retains duplicate uses so `StableId`
    /// reconciliation can observe metadata attached to any clone of a shared vertex.
    fn vertex_use_iter(&self) -> impl Iterator<Item = Vertex<P>> + '_ {
        self.edge_iter().flat_map(|edge| {
            let (front, back) = edge.absolute_ends();
            [front.clone(), back.clone()]
        })
    }

    /// invert all faces
    #[inline]
    pub fn not(&mut self) {
        self.boundaries
            .iter_mut()
            .flat_map(|shell| shell.face_iter_mut())
            .for_each(|face| {
                face.invert();
            });
    }

    fn unassigned_face_count(&self) -> usize {
        self.face_iter()
            .filter(|face| !face.stable_id().is_assigned())
            .count()
    }

    fn unassigned_vertex_count(&self) -> usize {
        unassigned_shared_count(
            self.vertex_use_iter()
                .map(|vertex| (vertex.id(), vertex.stable_id())),
        )
    }

    fn unassigned_edge_count(&self) -> usize {
        unassigned_shared_count(self.edge_iter().map(|edge| (edge.id(), edge.stable_id())))
    }

    /// Assigns a fresh [`StableId`] to every face that does not have one yet,
    /// leaving already-assigned faces untouched.
    ///
    /// Faces of a freshly built solid are [`StableId::UNASSIGNED`] until ids
    /// are explicitly assigned, so per-face attributes (e.g. a colour) cannot
    /// be keyed on them. Calling this first gives every face a stable id to key
    /// attributes on.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] before mutation if there are not enough identifiers left.
    pub fn ensure_face_stable_ids(&mut self) -> std::result::Result<(), StableIdExhausted> {
        self.ensure_allocator_above_existing_topology_ids()?;
        self.id_allocator
            .ensure_capacity(self.unassigned_face_count())?;
        let id_allocator = &mut self.id_allocator;
        for shell in &mut self.boundaries {
            for face in shell.face_iter_mut() {
                if !face.stable_id().is_assigned() {
                    face.set_stable_id(id_allocator.allocate()?);
                }
            }
        }
        Ok(())
    }

    /// Assigns a [`StableId`] to every vertex use that does not have one yet.
    ///
    /// Shared vertex uses receive the same id, preferring an existing assigned
    /// id for that topology key before allocating a fresh one.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] before mutation if there are not enough identifiers left.
    pub fn ensure_vertex_stable_ids(&mut self) -> std::result::Result<(), StableIdExhausted> {
        self.ensure_allocator_above_existing_topology_ids()?;
        self.id_allocator
            .ensure_capacity(self.unassigned_vertex_count())?;
        let mut ids = HashMap::<VertexId<P>, StableId>::default();
        for vertex in self.vertex_use_iter() {
            collect_assigned_vertex_stable_id(&vertex, &mut ids);
        }
        let id_allocator = &mut self.id_allocator;
        for shell in &mut self.boundaries {
            for face in shell.face_iter_mut() {
                for wire in &mut face.boundaries {
                    for edge in wire.edge_iter_mut() {
                        assign_vertex_stable_id(&mut edge.vertices.0, &mut ids, id_allocator)?;
                        assign_vertex_stable_id(&mut edge.vertices.1, &mut ids, id_allocator)?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Assigns a [`StableId`] to every edge use that does not have one yet.
    ///
    /// Shared edge uses receive the same id, preferring an existing assigned id
    /// for that topology key before allocating a fresh one.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] before mutation if there are not enough identifiers left.
    pub fn ensure_edge_stable_ids(&mut self) -> std::result::Result<(), StableIdExhausted> {
        self.ensure_allocator_above_existing_topology_ids()?;
        self.id_allocator
            .ensure_capacity(self.unassigned_edge_count())?;
        let mut ids = HashMap::<EdgeId<C>, StableId>::default();
        for edge in self.edge_iter() {
            collect_assigned_edge_stable_id(&edge, &mut ids);
        }
        let id_allocator = &mut self.id_allocator;
        for shell in &mut self.boundaries {
            for face in shell.face_iter_mut() {
                for wire in &mut face.boundaries {
                    for edge in wire.edge_iter_mut() {
                        let edge_id = edge.id();
                        let stable_id = if let Some(stable_id) = ids.get(&edge_id).copied() {
                            stable_id
                        } else {
                            let stable_id = id_allocator.allocate()?;
                            ids.insert(edge_id, stable_id);
                            stable_id
                        };
                        edge.set_stable_id(stable_id);
                    }
                }
            }
        }
        Ok(())
    }

    /// Assigns stable ids to every topology element that does not have one yet.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] if the identifier space is exhausted.
    pub fn ensure_topology_stable_ids(&mut self) -> std::result::Result<(), StableIdExhausted> {
        self.ensure_allocator_above_existing_topology_ids()?;
        let required = self
            .unassigned_vertex_count()
            .checked_add(self.unassigned_edge_count())
            .and_then(|count| count.checked_add(self.unassigned_face_count()))
            .ok_or(StableIdExhausted)?;
        self.id_allocator.ensure_capacity(required)?;
        self.ensure_vertex_stable_ids()?;
        self.ensure_edge_stable_ids()?;
        self.ensure_face_stable_ids()
    }

    /// Replaces every topology element's [`StableId`] with a caller-supplied
    /// value, so result ids derive from a caller-computed canonical ordering
    /// rather than construction/traversal order.
    ///
    /// `vertex_ids` and `edge_ids` map each shared vertex/edge topology key
    /// ([`Vertex::id`]/[`Edge::id`]) to the id every use of that element must
    /// carry, so all uses of one shared vertex or edge receive the same id.
    /// Faces are not shared, so `face_ids_in_order` supplies the id of each
    /// face in [`face_iter`](Self::face_iter) order. Any element absent from its
    /// map (or beyond `face_ids_in_order`) keeps its current id. The allocator
    /// is reset and then advanced past `next_free_id - 1`, so ids handed out
    /// afterwards never collide with the assigned range.
    ///
    /// Unlike [`ensure_topology_stable_ids`](Self::ensure_topology_stable_ids),
    /// which allocates in traversal order and skips already-assigned elements,
    /// this reassigns every mapped element unconditionally. The caller owns the
    /// canonical order; this method only applies it.
    ///
    /// # Errors
    /// Returns [`StableIdExhausted`] if advancing the allocator past the
    /// supplied stable-ID range would exhaust the ID space.
    pub fn apply_stable_ids(
        &mut self,
        vertex_ids: &HashMap<VertexId<P>, StableId>,
        edge_ids: &HashMap<EdgeId<C>, StableId>,
        face_ids_in_order: &[StableId],
        next_free_id: u64,
    ) -> std::result::Result<(), StableIdExhausted> {
        let mut next_allocator = StableIdAllocator::new();
        next_allocator.ensure_above(next_free_id.saturating_sub(1))?;
        let mut face_ids = face_ids_in_order.iter().copied();
        for shell in &mut self.boundaries {
            for face in shell.face_iter_mut() {
                if let Some(id) = face_ids.next() {
                    face.set_stable_id(id);
                }
                for wire in &mut face.boundaries {
                    for edge in wire.edge_iter_mut() {
                        if let Some(id) = edge_ids.get(&edge.id()).copied() {
                            edge.set_stable_id(id);
                        }
                        if let Some(id) = vertex_ids.get(&edge.vertices.0.id()).copied() {
                            edge.vertices.0.set_stable_id(id);
                        }
                        if let Some(id) = vertex_ids.get(&edge.vertices.1.id()).copied() {
                            edge.vertices.1.set_stable_id(id);
                        }
                    }
                }
            }
        }
        self.id_allocator = next_allocator;
        Ok(())
    }

    fn ensure_allocator_above_existing_topology_ids(
        &mut self,
    ) -> std::result::Result<(), StableIdExhausted> {
        let max_id = self
            .boundaries
            .iter()
            .flat_map(|shell| shell.face_iter())
            .fold(0_u64, |max_id, face| {
                let max_id = max_id.max(face.stable_id().raw());
                face.boundaries
                    .iter()
                    .flat_map(|wire| wire.iter())
                    .fold(max_id, |max_id, edge| {
                        max_id
                            .max(edge.stable_id().raw())
                            .max(edge.vertices.0.stable_id().raw())
                            .max(edge.vertices.1.stable_id().raw())
                    })
            });
        self.id_allocator.ensure_above(max_id)
    }

    /// Returns a new solid whose surfaces are mapped by `surface_mapping`,
    /// curves are mapped by `curve_mapping` and points are mapped by `point_mapping`.
    /// # Remarks
    /// Accessing geometry elements directly in the closure will result in a deadlock.
    /// So, this method does not appear to the document.
    #[doc(hidden)]
    #[inline]
    pub fn try_mapped<Q, D, T>(
        &self,
        mut point_mapping: impl FnMut(&P) -> Option<Q>,
        mut curve_mapping: impl FnMut(&C) -> Option<D>,
        mut surface_mapping: impl FnMut(&S) -> Option<T>,
    ) -> Option<Solid<Q, D, T>> {
        // Spec 012 U4 / ledger C11 instance 2. This is the exact site the 011
        // backtrace landed on: `extract_healed_trimmed_solid` -> `erase_trims`
        // -> here, with a shell that is short a face. It used to abort the
        // process in debug and hand back `Some(invalid solid)` in release.
        // `.ok()` makes both profiles answer the same shape of thing -- a
        // `None` this function's signature already promised.
        Solid::new(
            self.boundaries()
                .iter()
                .map(move |shell| {
                    shell.try_mapped(&mut point_mapping, &mut curve_mapping, &mut surface_mapping)
                })
                .collect::<Option<Vec<_>>>()?,
        )
        .ok()
    }

    /// Returns a new solid whose surfaces are mapped by `surface_mapping`,
    /// curves are mapped by `curve_mapping` and points are mapped by `point_mapping`.
    /// # Remarks
    /// Accessing geometry elements directly in the closure will result in a deadlock.
    /// So, this method does not appear to the document.
    #[doc(hidden)]
    #[inline]
    pub fn mapped<Q, D, T>(
        &self,
        mut point_mapping: impl FnMut(&P) -> Q,
        mut curve_mapping: impl FnMut(&C) -> D,
        mut surface_mapping: impl FnMut(&S) -> T,
    ) -> Solid<Q, D, T> {
        // Infallible signature, so there is no channel to report a violation
        // on -- and `new_unchecked` is stated here rather than hidden behind
        // the former profile-dependent constructor, because the alternative was a panic in
        // debug and this exact call in release (spec 012 U4).
        //
        // It costs nothing, and that is provable rather than hoped for.
        // `Solid::new`'s condition is `Shell::check_solid_boundary`, which
        // is PURELY COMBINATORIAL -- non-empty, connected, closed, no singular
        // vertex -- and `mapped` rebuilds the identical face/edge/vertex
        // incidence structure with only the geometry replaced. So the result
        // satisfies the condition exactly when the RECEIVER does, and a
        // violation here is always a violation the receiver arrived with
        // (`Solid::new_unchecked` / `TrimmedSolid::erase_trims` upstream). This
        // is not the place that can diagnose it: see `try_mapped` above, which
        // has an `Option` and now uses it.
        Solid::new_unchecked(
            self.boundaries()
                .iter()
                .map(move |shell| {
                    shell.mapped(&mut point_mapping, &mut curve_mapping, &mut surface_mapping)
                })
                .collect(),
        )
    }

    /// Returns the consistence of the geometry of end vertices
    /// and the geometry of edge.
    #[inline]
    #[must_use]
    pub fn is_geometric_consistent(&self) -> bool
    where
        P: Tolerance,
        C: BoundedCurve<Point = P>,
        S: IncludeCurve<C>, {
        self.boundaries()
            .iter()
            .all(super::Shell::is_geometric_consistent)
    }

    /// Cuts one edge into two edges at vertex.
    #[inline]
    pub fn cut_edge(
        &mut self,
        edge_id: EdgeId<C>,
        vertex: &Vertex<P>,
    ) -> Option<(Edge<P, C>, Edge<P, C>)>
    where
        P: Clone,
        C: Cut<Point = P> + SearchParameter<CurveParameter, Point = P>,
    {
        let res = self
            .boundaries
            .iter_mut()
            .find_map(|shell| shell.cut_edge(edge_id, vertex));
        #[cfg(debug_assertions)]
        let _ = Self::new_unchecked(self.boundaries.clone());
        res
    }
    /// Removes `vertex` from `self` by concat two edges on both sides.
    #[inline]
    pub fn remove_vertex_by_concat_edges(&mut self, vertex_id: VertexId<P>) -> Option<Edge<P, C>>
    where
        P: Debug,
        C: Concat<C, Point = P, Output = C> + Invertible + ParameterTransform, {
        let res = self
            .boundaries
            .iter_mut()
            .find_map(|shell| shell.remove_vertex_by_concat_edges(vertex_id));
        #[cfg(debug_assertions)]
        let _ = Self::new_unchecked(self.boundaries.clone());
        res
    }

    /// Cut a face with `face_id` by edge.
    #[inline]
    pub fn cut_face_by_edge(&mut self, face_id: FaceId<S>, edge: Edge<P, C>) -> bool
    where S: Clone {
        let tuple = self.boundaries.iter_mut().find_map(|shell| {
            let find_res = shell
                .face_iter_mut()
                .enumerate()
                .find(move |(_, face)| face.id() == face_id)
                .map(move |(i, _)| i);
            find_res.map(move |i| (shell, i))
        });
        if let Some((shell, i)) = tuple
            && let Some((face0, face1)) = shell[i].cut_by_edge(edge)
        {
            shell[i] = face0;
            shell.push(face1);
            return true;
        }
        false
    }

    /// Creates display struct for debugging the solid.
    #[inline]
    #[must_use]
    pub const fn display(
        &self,
        format: SolidDisplayFormat,
    ) -> DebugDisplay<'_, Self, SolidDisplayFormat> {
        DebugDisplay {
            entity: self,
            format,
        }
    }
}

impl<P: Clone, C: Clone, S: Clone> Solid<P, C, Option<S>> {
    /// Returns the value with the Option removed if there is no `None` in the surfaces of the faces.
    #[inline]
    pub fn collect_option(&self) -> Option<Solid<P, C, S>> {
        self.try_mapped(|x| Some(x.clone()), |x| Some(x.clone()), Clone::clone)
    }
}

impl<P, C, S> PartialEq for Solid<P, C, S> {
    fn eq(&self, other: &Self) -> bool { self.boundaries == other.boundaries }
}

impl<P, C, S> Eq for Solid<P, C, S> {}

impl<P: Debug, C: Debug, S: Debug> Debug for DebugDisplay<'_, Solid<P, C, S>, SolidDisplayFormat> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self.format {
            SolidDisplayFormat::ShellsList { shell_format } => f
                .debug_list()
                .entries(
                    self.entity
                        .boundaries
                        .iter()
                        .map(|shell| shell.display(shell_format)),
                )
                .finish(),
            SolidDisplayFormat::ShellsListTuple { shell_format } => f
                .debug_tuple("Solid")
                .field(&DebugDisplay {
                    entity: self.entity,
                    format: SolidDisplayFormat::ShellsList { shell_format },
                })
                .finish(),
            SolidDisplayFormat::Struct { shell_format } => f
                .debug_struct("Solid")
                .field(
                    "boundaries",
                    &DebugDisplay {
                        entity: self.entity,
                        format: SolidDisplayFormat::ShellsList { shell_format },
                    },
                )
                .finish(),
        }
    }
}

#[cfg(test)]
pub(super) fn cube() -> Solid<(), (), ()> {
    use crate::*;
    let v = Vertex::from_points([(); 8]);
    let edge = [
        Edge::new_unchecked(&v[0], &v[1], ()), // 0
        Edge::new_unchecked(&v[1], &v[2], ()), // 1
        Edge::new_unchecked(&v[2], &v[3], ()), // 2
        Edge::new_unchecked(&v[3], &v[0], ()), // 3
        Edge::new_unchecked(&v[0], &v[4], ()), // 4
        Edge::new_unchecked(&v[1], &v[5], ()), // 5
        Edge::new_unchecked(&v[2], &v[6], ()), // 6
        Edge::new_unchecked(&v[3], &v[7], ()), // 7
        Edge::new_unchecked(&v[4], &v[5], ()), // 8
        Edge::new_unchecked(&v[5], &v[6], ()), // 9
        Edge::new_unchecked(&v[6], &v[7], ()), // 10
        Edge::new_unchecked(&v[7], &v[4], ()), // 11
    ];

    let wire0 = wire![&edge[0], &edge[1], &edge[2], &edge[3]];
    let face0 = Face::new_unchecked(vec![wire0], ());

    let wire1 = wire![&edge[4], &edge[8], &edge[5].inverse(), &edge[0].inverse(),];
    let face1 = Face::new_unchecked(vec![wire1], ());

    let wire2 = wire![&edge[5], &edge[9], &edge[6].inverse(), &edge[1].inverse(),];
    let face2 = Face::new_unchecked(vec![wire2], ());

    let wire3 = wire![&edge[6], &edge[10], &edge[7].inverse(), &edge[2].inverse(),];
    let face3 = Face::new_unchecked(vec![wire3], ());
    let wire4 = wire![&edge[7], &edge[11], &edge[4].inverse(), &edge[3].inverse(),];
    let face4 = Face::new_unchecked(vec![wire4], ());
    let wire5 = wire![
        &edge[11].inverse(),
        &edge[10].inverse(),
        &edge[9].inverse(),
        &edge[8].inverse(),
    ];
    let face5 = Face::new_unchecked(vec![wire5], ());

    let mut shell = Shell::new();
    shell.push(face0);
    shell.push(face5);
    assert!(!shell.is_connected());
    shell.push(face1);
    assert_eq!(shell.shell_condition(), ShellCondition::Oriented);
    assert!(shell.is_connected());
    shell.push(face2);
    shell.push(face3);
    shell.push(face4);

    Solid::new_unchecked(vec![shell])
}

#[test]
fn cube_test() { cube(); }

#[test]
fn ensure_face_stable_ids_assigns_unassigned_faces() -> std::result::Result<(), StableIdExhausted> {
    use monstertruck_core::StableId;
    let mut solid = cube();
    // A freshly built cube has unassigned face ids.
    assert!(solid.face_iter().all(|f| !f.stable_id().is_assigned()));

    solid.ensure_face_stable_ids()?;
    let ids: Vec<StableId> = solid.face_iter().map(|f| f.stable_id()).collect();
    assert_eq!(ids.len(), 6);
    assert!(ids.iter().all(|id| id.is_assigned()), "all faces assigned");
    let mut unique = ids.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(unique.len(), 6, "face ids are unique");

    // Idempotent: a second call leaves already-assigned faces untouched.
    solid.ensure_face_stable_ids()?;
    let ids_again: Vec<StableId> = solid.face_iter().map(|f| f.stable_id()).collect();
    assert_eq!(ids, ids_again, "second call is a no-op for assigned faces");
    Ok(())
}

#[test]
fn ensure_vertex_stable_ids_assigns_each_topological_vertex_once()
-> std::result::Result<(), StableIdExhausted> {
    use monstertruck_core::StableId;
    let mut solid = cube();
    assert!(
        solid
            .vertex_use_iter()
            .all(|vertex| !vertex.stable_id().is_assigned())
    );

    solid.ensure_vertex_stable_ids()?;
    let ids: Vec<StableId> = solid
        .vertex_use_iter()
        .map(|vertex| vertex.stable_id())
        .collect();
    assert!(
        ids.iter().all(|id| id.is_assigned()),
        "all vertex uses assigned"
    );
    let mut unique = ids.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(unique.len(), 8, "shared vertex uses receive one id");

    solid.ensure_vertex_stable_ids()?;
    let ids_again: Vec<StableId> = solid
        .vertex_use_iter()
        .map(|vertex| vertex.stable_id())
        .collect();
    assert_eq!(
        ids, ids_again,
        "second call is a no-op for assigned vertices"
    );
    Ok(())
}

#[test]
fn ensure_vertex_stable_ids_preserves_id_found_only_on_back_endpoint()
-> std::result::Result<(), StableIdExhausted> {
    use monstertruck_core::StableId;
    let mut solid = cube();
    let target_id = solid
        .edge_iter()
        .next()
        .expect("cube has edges")
        .absolute_back()
        .id();
    let existing_id = StableId::new(100);
    let mut assigned_one_back_use = false;

    'outer: for shell in &mut solid.boundaries {
        for face in shell.face_iter_mut() {
            for wire in &mut face.boundaries {
                for edge in wire.edge_iter_mut() {
                    if edge.vertices.1.id() == target_id {
                        edge.vertices.1.set_stable_id(existing_id);
                        assigned_one_back_use = true;
                        break 'outer;
                    }
                }
            }
        }
    }
    assert!(assigned_one_back_use);

    solid.ensure_vertex_stable_ids()?;
    let target_ids: Vec<StableId> = solid
        .vertex_use_iter()
        .filter(|vertex| vertex.id() == target_id)
        .map(|vertex| vertex.stable_id())
        .collect();
    assert!(!target_ids.is_empty());
    assert!(target_ids.iter().all(|&id| id == existing_id));
    assert!(
        solid
            .vertex_use_iter()
            .all(|vertex| vertex.stable_id().is_assigned())
    );
    assert!(solid.id_allocator().peek() > existing_id.raw());
    Ok(())
}

#[test]
fn ensure_edge_stable_ids_assigns_each_topological_edge_once()
-> std::result::Result<(), StableIdExhausted> {
    use monstertruck_core::StableId;
    let mut solid = cube();
    assert!(solid.edge_iter().all(|e| !e.stable_id().is_assigned()));

    solid.ensure_edge_stable_ids()?;
    let ids: Vec<StableId> = solid.edge_iter().map(|e| e.stable_id()).collect();
    assert_eq!(ids.len(), 24);
    assert!(
        ids.iter().all(|id| id.is_assigned()),
        "all edge uses assigned"
    );
    let mut unique = ids.clone();
    unique.sort();
    unique.dedup();
    assert_eq!(unique.len(), 12, "shared edge uses receive one id");

    solid.ensure_edge_stable_ids()?;
    let ids_again: Vec<StableId> = solid.edge_iter().map(|e| e.stable_id()).collect();
    assert_eq!(ids, ids_again, "second call is a no-op for assigned edges");
    Ok(())
}

#[test]
fn ensure_edge_stable_ids_preserves_existing_shared_edge_id()
-> std::result::Result<(), StableIdExhausted> {
    use monstertruck_core::StableId;
    let mut solid = cube();
    let target_id = solid.edge_iter().next().expect("cube has edges").id();
    let existing_id = StableId::new(100);
    let mut seen = 0;

    for shell in &mut solid.boundaries {
        for face in shell.face_iter_mut() {
            for wire in &mut face.boundaries {
                for edge in wire.edge_iter_mut() {
                    if edge.id() == target_id {
                        seen += 1;
                        if seen == 2 {
                            edge.set_stable_id(existing_id);
                        }
                    }
                }
            }
        }
    }
    assert_eq!(seen, 2, "cube edge should have two edge uses");

    solid.ensure_edge_stable_ids()?;
    let target_ids: Vec<StableId> = solid
        .edge_iter()
        .filter(|edge| edge.id() == target_id)
        .map(|edge| edge.stable_id())
        .collect();
    assert_eq!(target_ids, vec![existing_id, existing_id]);
    assert!(solid.edge_iter().all(|edge| edge.stable_id().is_assigned()));
    assert!(solid.id_allocator().peek() > existing_id.raw());
    Ok(())
}
