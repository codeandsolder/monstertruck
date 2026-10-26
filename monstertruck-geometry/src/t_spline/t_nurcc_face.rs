use super::*;

impl<P> TnurccFace<P> {
    /// Returns all boundary vertices in anti-clockwise order.
    pub fn boundary_vertices(
        f: Arc<RwLock<Self>>,
    ) -> Result<Vec<Arc<RwLock<TnurccControlPoint<P>>>>> {
        let edge = f
            .read()
            .edge
            .as_ref()
            .map(Arc::clone)
            .ok_or(Error::TnurccMalformedFace)?;
        let face_side = edge
            .read()
            .face_side(Arc::clone(&f))
            .ok_or(Error::TnurccMalformedFace)?;
        let iter = TnurccAcwFaceIter::try_from_edge(Arc::clone(&edge), face_side)
            .ok_or(Error::TnurccMalformedFace)?;

        iter.map(|edge| {
            let side = edge
                .read()
                .face_side(Arc::clone(&f))
                .ok_or(Error::TnurccMalformedFace)?;
            Ok(match side {
                TnurccFaceSide::Left => Arc::clone(&edge.read().origin),
                TnurccFaceSide::Right => Arc::clone(&edge.read().dest),
            })
        })
        .collect()
    }

    /// Returns all boundary edges in anti-clockwise order.
    pub fn border_edges(f: Arc<RwLock<Self>>) -> Result<Vec<Arc<RwLock<TnurccEdge<P>>>>> {
        let edge = f
            .read()
            .edge
            .as_ref()
            .map(Arc::clone)
            .ok_or(Error::TnurccMalformedFace)?;
        let face_side = edge
            .read()
            .face_side(Arc::clone(&f))
            .ok_or(Error::TnurccMalformedFace)?;
        let iter =
            TnurccAcwFaceIter::try_from_edge(edge, face_side).ok_or(Error::TnurccMalformedFace)?;
        Ok(iter.collect())
    }
}

impl<P> Drop for TnurccFace<P> {
    fn drop(&mut self) {
        for i in 0..self.corners.len() {
            self.corners[i] = None;
        }
        self.edge = None;
    }
}
