use super::{Arc, Error, Result, RwLock, TnurccAcwPointIter, TnurccControlPoint, TnurccEdge};

impl<P> TnurccControlPoint<P> {
    /// Creates a new `TnurccControlPoint` instance with index `index`, cartesian point `point`,
    /// valence `0`, and no incoming edge (`None`).
    pub const fn new(index: usize, point: P) -> Self {
        Self {
            index,
            valence: 0,
            point,
            incoming_edge: None,
        }
    }

    /// Returns a vector containing all edges connected to `p`. Returns an empty vector if `p` does
    /// not have an `incoming_edge`. Collects in an anti-clockwise fashion starting from `p`'s
    /// reference `incoming_edge`.
    ///
    /// # Panics
    /// Panics if `p`'s `incoming_edge` does not have `p` as an end.
    pub fn radial_edges(p: &Arc<RwLock<Self>>) -> Result<Vec<Arc<RwLock<TnurccEdge<P>>>>> {
        let edge = p
            .read()
            .incoming_edge
            .as_ref()
            .map(Arc::clone)
            .ok_or(Error::TnurccMalformedFace)?;
        let point_end = edge.read().point_end(p).ok_or(Error::TnurccMalformedFace)?;
        Ok(TnurccAcwPointIter::from_edge(edge, point_end).collect())
    }

    /// Returns the edge connecting two points, if it exists.
    ///
    /// # Returns
    /// - `None` if `center` has no `incoming_edge`, if that `incoming_edge` is incorrectly assigned to `center`.
    ///   or if no edge connects `center` to `op`.
    ///
    /// - `Some(edge)` if the edge is found.
    ///
    /// # Borrows
    /// Immutably borrows `center`, `op`, and any edges connected to `center`.
    ///
    /// # Panics
    /// Panics if any borrow fails
    pub fn edge_from_opposing_point(
        center: &Arc<RwLock<Self>>,
        op: &Arc<RwLock<Self>>,
    ) -> Option<Arc<RwLock<TnurccEdge<P>>>> {
        let incoming_edge = center.read().incoming_edge.as_ref().map(Arc::clone);
        let in_edge = incoming_edge?;
        let endpoint = in_edge.read().point_end(center)?;
        TnurccAcwPointIter::from_edge(Arc::clone(&in_edge), endpoint)
            .find(|edge| edge.read().point_end(op).is_some())
    }
}

impl<P> Drop for TnurccControlPoint<P> {
    fn drop(&mut self) { self.incoming_edge = None; }
}
