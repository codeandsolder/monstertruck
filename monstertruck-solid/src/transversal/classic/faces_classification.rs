//! Classic (0.3.2) per-face and/or/unknown classification.
//!
//! Ported verbatim from the published 0.3.2 crate's
//! `transversal::faces_classification`. Kept private to the classic subtree so
//! it binds to the classic [`ShapesOpStatus`] (the crate-level
//! `faces_classification` binds to the loops-store status of the same
//! name, a distinct type).

use super::loops_store::ShapesOpStatus;
use monstertruck_topology::*;
use rustc_hash::FxHashMap as HashMap;

#[derive(Clone, Debug)]
pub(super) struct FacesClassification<P, C, S> {
    shell: Shell<P, C, S>,
    status: HashMap<FaceId<S>, ShapesOpStatus>,
}

impl<P, C, S> Default for FacesClassification<P, C, S> {
    fn default() -> Self {
        Self {
            shell: Default::default(),
            status: HashMap::default(),
        }
    }
}

impl<P, C, S> FacesClassification<P, C, S> {
    pub(super) fn push(&mut self, face: Face<P, C, S>, status: ShapesOpStatus) {
        self.status.insert(face.id(), status);
        self.shell.push(face);
    }

    pub(super) fn and_or_unknown(&self) -> [Shell<P, C, S>; 3] {
        let [mut and, mut or, mut unknown] = <[Shell<P, C, S>; 3]>::default();
        for face in &self.shell {
            // SAFETY: `push()` inserts every face id into `status`.
            match self.status.get(&face.id()).unwrap() {
                ShapesOpStatus::And => and.push(face.clone()),
                ShapesOpStatus::Or => or.push(face.clone()),
                ShapesOpStatus::Unknown => unknown.push(face.clone()),
            }
        }
        [and, or, unknown]
    }

    pub(super) fn integrate_by_component(&mut self) {
        let [and, or, unknown] = self.and_or_unknown();
        let and_boundary = and.extract_boundaries();
        let or_boundary = or.extract_boundaries();
        let components = unknown.connected_components();
        for comp in components {
            let boundary = comp.extract_boundaries();
            let touches_and = boundary.iter().flatten().any(|boundary_edge| {
                and_boundary
                    .iter()
                    .flatten()
                    .any(|edge| edge.id() == boundary_edge.id())
            });
            let touches_or = boundary.iter().flatten().any(|boundary_edge| {
                or_boundary
                    .iter()
                    .flatten()
                    .any(|edge| edge.id() == boundary_edge.id())
            });
            let status = match (touches_and, touches_or) {
                (true, false) => ShapesOpStatus::And,
                (false, true) => ShapesOpStatus::Or,
                // Mixed evidence means this connected unknown region bridges
                // both classified sides. A boundary-less component has no
                // topological evidence at all. In both cases leave it unknown
                // for geometric classification instead of letting edge order
                // pick a status for the whole component.
                _ => continue,
            };
            comp.iter().for_each(|face| {
                // SAFETY: face originated from `self.shell` via `and_or_unknown()`.
                *self.status.get_mut(&face.id()).unwrap() = status;
            });
        }
    }
}
