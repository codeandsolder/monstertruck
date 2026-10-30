use super::{
    Arc, RwLock, TnurccAcwFaceIter, TnurccAcwPointIter, TnurccEdge, TnurccFace, TnurccFaceSide,
    TnurccVertexEnd,
};

impl<P> TnurccAcwPointIter<P> {
    /// Creates a new `TnurccAcwPointIter` which iterates the edges around `e`'s vertex `end` in an anti-clockwise manner.
    /// Returns `None` without making a full rotation, that is, will return `e` as the first element, but not the last.
    ///
    /// # Borrows
    /// Immutably borrows every edge connected to the point at `e`'s `end` when calling `next`.
    pub fn from_edge(e: Arc<RwLock<TnurccEdge<P>>>, end: TnurccVertexEnd) -> Self {
        let point = Arc::clone(&e.read().point_at_end(end));
        let start = Arc::clone(&e);
        Self {
            point,
            start,
            cur: Some(e),
        }
    }
}

impl<P> TnurccAcwFaceIter<P> {
    /// Creates a new `TnurccAcwFaceIter` which iterates the edges around `e`'s face `side` in an anti-clockwise manner.
    /// `next` returns `None` without making a full rotation, that is, will return `e` as the first element, but not the last.
    ///
    /// # Returns
    /// - `None` if `e` does not have a face on `side`.
    ///
    /// - `Some(iter)` otherwise.
    pub fn try_from_edge(e: Arc<RwLock<TnurccEdge<P>>>, side: TnurccFaceSide) -> Option<Self> {
        let face = e.read().face_from_side(side)?;
        let start = Arc::clone(&e);
        Some(Self {
            face,
            start,
            cur: Some(e),
        })
    }

    /// Creates a new `TnurccAcwFaceIter` which iterates the edges around `f` in an anti-clockwise manner.
    /// `next` returns `None` without making a full rotation, that is, will not return the first edge twice.
    ///
    /// # Returns
    /// - `None` if `f` does not have a reference edge.
    ///
    /// - `Some(iter)` otherwise.
    #[expect(
        dead_code,
        reason = "face-seeded traversal constructor complements the edge-seeded iterator API"
    )]
    pub fn try_from_face(f: Arc<RwLock<TnurccFace<P>>>) -> Option<Self> {
        let edge = f.read().edge.as_ref().map(Arc::clone)?;
        let start = Arc::clone(&edge);
        Some(Self {
            face: f,
            start,
            cur: Some(edge),
        })
    }
}

impl<P> Iterator for TnurccAcwPointIter<P> {
    type Item = Arc<RwLock<TnurccEdge<P>>>;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = Arc::clone(self.cur.as_ref()?);
        let ret = Some(Arc::clone(&edge));

        // Is point the origin or dest?
        let end = edge.read().point_end(&self.point)?;

        // Get the next ACW edge for point
        let new_edge = edge.read().acw_edge_from_end(end)?;

        // If the new edge is the starting edge, stop the iterator by setting cur to none
        // Otherwise, keep going
        if std::ptr::eq(self.start.as_ref(), new_edge.as_ref()) {
            self.cur = None;
        } else {
            self.cur = Some(new_edge);
        }

        ret
    }
}

impl<P> Iterator for TnurccAcwFaceIter<P> {
    type Item = Arc<RwLock<TnurccEdge<P>>>;

    fn next(&mut self) -> Option<Self::Item> {
        let edge = Arc::clone(self.cur.as_ref()?);
        let ret = Some(Arc::clone(&edge));

        // Is point the origin or dest?
        let side = edge.read().face_side(&self.face)?;

        // Get the next ACW edge for point
        let new_edge = edge.read().acw_edge_from_side(side)?;

        // If the new edge is the starting edge, stop the iterator by setting cur to none
        // Otherwise, keep going
        if std::ptr::eq(self.start.as_ref(), new_edge.as_ref()) {
            self.cur = None;
        } else {
            self.cur = Some(new_edge);
        }

        ret
    }
}
