//! Section-local adjacency of actual source face and edge identities.
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use std::collections::BTreeMap;

pub(super) struct FaceIncidence<'a> {
    faces: BTreeMap<u64, &'a ReducedFace>,
    edge_faces: BTreeMap<u64, Vec<&'a ReducedFace>>,
}

impl<'a> FaceIncidence<'a> {
    pub(super) fn new(solid: &'a ReducedSolid) -> Self {
        let mut edge_faces = BTreeMap::<_, Vec<_>>::new();
        let mut faces = BTreeMap::new();
        for face in solid.faces() {
            faces.insert(face.id, face);
            for edge in face.loop_edges() {
                edge_faces.entry(edge.edge_id).or_default().push(face);
            }
        }
        Self { faces, edge_faces }
    }

    pub(super) fn face(&self, id: u64) -> Option<&'a ReducedFace> { self.faces.get(&id).copied() }

    pub(super) fn owners(&self, edge: u64) -> Option<&[&'a ReducedFace]> {
        Some(self.edge_faces.get(&edge)?.as_slice())
    }

    pub(super) fn neighbor(&self, face: u64, edge: u64) -> Option<&'a ReducedFace> {
        let [first, second] = self.owners(edge)? else {
            return None;
        };
        match (first.id == face, second.id == face) {
            (true, false) => Some(second),
            (false, true) => Some(first),
            _ => None,
        }
    }
}
