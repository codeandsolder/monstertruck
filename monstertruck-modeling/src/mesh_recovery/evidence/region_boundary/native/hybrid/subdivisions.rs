//! Chart seams may subdivide a lifted source edge after its approximation
//! was certified. Keep the original prescribed edge for the facet proof,
//! then propagate the topology-only subdivision to its patch face.

use super::{HybridRecoveryError, NativeRegionFace, NativeTrimError};
use crate::{Curve, Edge, Face, Wire};
use monstertruck_topology::EdgeId;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::collections::BTreeMap;

type SubdivisionMap = HashMap<EdgeId<Curve>, Vec<Edge>>;

fn subdivision_map(
    original: &BTreeMap<usize, Vec<Edge>>,
    final_trims: &[Result<Vec<Edge>, NativeTrimError>],
) -> Result<SubdivisionMap, HybridRecoveryError> {
    let mut result = HashMap::default();
    for (&id, old_edges) in original {
        let mut remaining = final_trims
            .get(id)
            .ok_or(HybridRecoveryError::InvalidInput)?
            .as_ref()
            .map_err(|&error| error)?
            .as_slice();
        for old in old_edges {
            let last = remaining
                .iter()
                .position(|edge| edge.back() == old.back())
                .ok_or(HybridRecoveryError::InvalidInput)?;
            let (pieces, tail) = remaining.split_at(last + 1);
            remaining = tail;
            if pieces[0].front() != old.front()
                || pieces.windows(2).any(|p| p[0].back() != p[1].front())
            {
                return Err(HybridRecoveryError::InvalidInput);
            }
            if pieces.len() == 1 {
                if pieces[0].id() != old.id() {
                    return Err(HybridRecoveryError::InvalidInput);
                }
            } else if result.insert(old.id(), pieces.to_vec()).is_some() {
                return Err(HybridRecoveryError::InvalidInput);
            }
        }
        if !remaining.is_empty() {
            return Err(HybridRecoveryError::InvalidInput);
        }
    }
    Ok(result)
}

fn replace_face(
    face: &Face,
    subdivisions: &SubdivisionMap,
    used: &mut HashSet<EdgeId<Curve>>,
) -> Result<Option<Face>, HybridRecoveryError> {
    if !face
        .absolute_boundaries()
        .iter()
        .flat_map(|wire| wire.edge_iter())
        .any(|edge| subdivisions.contains_key(&edge.id()))
    {
        return Ok(None);
    }
    let mut wires = Vec::new();
    for wire in face.boundaries() {
        let mut boundary = Vec::new();
        for edge in wire.edge_iter() {
            let Some(pieces) = subdivisions.get(&edge.id()) else {
                boundary.push(edge.clone());
                continue;
            };
            // A prescribed edge is on the outside of exactly one patch.
            // Every subdivision must therefore be propagated exactly once.
            if !used.insert(edge.id()) {
                return Err(HybridRecoveryError::InvalidInput);
            }
            let first = &pieces[0];
            let last = pieces.last().ok_or(HybridRecoveryError::InvalidInput)?;
            if edge.front() == first.front() && edge.back() == last.back() {
                boundary.extend(pieces.iter().cloned());
            } else if edge.front() == last.back() && edge.back() == first.front() {
                boundary.extend(pieces.iter().rev().map(Edge::inverse));
            } else {
                return Err(HybridRecoveryError::InvalidInput);
            }
        }
        wires.push(Wire::from(boundary));
    }
    let mut replacement =
        Face::new(wires, face.oriented_surface()).map_err(|_| NativeTrimError::FaceTopology)?;
    replacement.set_stable_id(face.stable_id());
    if !replacement.is_geometric_consistent() {
        return Err(NativeTrimError::FaceGeometry.into());
    }
    Ok(Some(replacement))
}

/// Only chart planning can modify these edges, and it does so exclusively
/// by cutting the existing native curve. This is a subdivision of the
/// certified boundary locus, not another approximation or carrier fit.
pub(super) fn propagate(
    regions: &mut [NativeRegionFace],
    original: &BTreeMap<usize, Vec<Edge>>,
    final_trims: &[Result<Vec<Edge>, NativeTrimError>],
) -> Result<(), HybridRecoveryError> {
    let subdivisions = subdivision_map(original, final_trims)?;
    if subdivisions.is_empty() {
        return Ok(());
    }
    let mut used = HashSet::default();
    for region in regions {
        for face in region.faces.as_mut().map_err(|&mut error| error)? {
            if let Some(replacement) = replace_face(face, &subdivisions, &mut used)? {
                *face = replacement;
            }
        }
    }
    if used.len() != subdivisions.len() {
        return Err(HybridRecoveryError::InvalidInput);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Plane, Point3, Surface, builder};

    #[test]
    fn lifted_subdivision_reuses_both_orientations_without_changing_face_locus()
    -> Result<(), Box<dyn std::error::Error>> {
        let a = builder::vertex(Point3::new(0.0, 0.0, 0.0));
        let b = builder::vertex(Point3::new(1.0, 0.0, 0.0));
        let c = builder::vertex(Point3::new(0.0, 1.0, 0.0));
        let middle = builder::vertex(Point3::new(0.3, 0.0, 0.0));
        let edge = builder::line(&a, &b);
        let (left, right) = edge.cut(&middle).ok_or("cut")?;
        let original = BTreeMap::from([(0, vec![edge.clone()])]);
        let final_trims = vec![Ok(vec![left.clone(), right.clone()])];
        let face = Face::new(
            vec![Wire::from(vec![
                edge,
                builder::line(&b, &c),
                builder::line(&c, &a),
            ])],
            Surface::Plane(Plane::new(a.point(), b.point(), c.point())),
        )?;
        for source in [face.clone(), face.inverse()] {
            let mut regions = vec![NativeRegionFace {
                region: 0,
                faces: Ok(vec![source]),
                requires_chart_cut: false,
                chart_seams: 0,
            }];
            propagate(&mut regions, &original, &final_trims)?;
            let faces = regions[0].faces.as_ref().map_err(|&error| error)?;
            assert_eq!(faces.len(), 1);
            assert!(faces[0].is_geometric_consistent());
            let edges = faces[0].boundaries();
            assert_eq!(edges[0].len(), 4);
            for piece in [&left, &right] {
                assert_eq!(
                    edges[0]
                        .edge_iter()
                        .filter(|e| e.id() == piece.id())
                        .count(),
                    1
                );
            }
        }
        Ok(())
    }
}
