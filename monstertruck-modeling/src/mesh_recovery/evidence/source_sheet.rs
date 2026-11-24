//! Optional cancellation of coincident, oppositely oriented planar sheets.
//! This preserves the oriented surface chain, not the original CAD topology:
//! an undersampled feature may have collapsed to a zero-thickness ribbon.
//! Never infer its missing width, profile or volume from cancellation.

use super::super::{
    EdgeIncidence, EdgeKey, InputPolygon, MeshRecoveryError, edge_key, input_polygon_from_ids,
    source_edge_incidence,
};
use super::{EvidenceOptions, MeshEvidence, analyze_welded_evidence, validate_evidence_options};
use crate::{InnerSpace, KERNEL_EPSILON};
use std::collections::{BTreeMap, BTreeSet};

/// An oppositely oriented pair of source planar sheets with the same
/// geometric support and identical reverse directed boundary.
#[derive(Clone, Debug)]
pub struct CanceledSourceSheet {
    /// Source polygon IDs on each side before normalization.
    pub sides: [Vec<usize>; 2],
    /// Ambiguous original edges whose incidence this cancellation resolves.
    pub resolved_edges: Vec<[usize; 2]>,
    /// Largest source vertex residual from the common plane.
    pub max_plane_residual: f64,
}

/// A separate evidence table preserving a map to the unmodified source.
#[derive(Debug)]
pub struct NormalizedSourceEvidence {
    /// Evidence on retained source polygons. Canonical vertex IDs unchanged.
    pub evidence: MeshEvidence,
    /// Active polygon index to original input polygon index.
    pub original_polygon_ids: Vec<usize>,
    /// Every canceled source sheet pair and its certificate.
    pub canceled_sheets: Vec<CanceledSourceSheet>,
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct SheetKey {
    boundary: Vec<EdgeKey>,
}
struct Sheet {
    support: usize,
    directions: Vec<bool>,
    vertices: Vec<usize>,
}

impl MeshEvidence {
    /// Cancel uniquely paired opposite planar sheets that have identical
    /// polygonal support domains and cause ambiguous source incidence.
    /// The input remains immutable. No vertices move, and no new open
    /// edges or ambiguous edges may be introduced. This is an opt-in
    /// oriented-chain simplification, not a certificate of the intended
    /// solid topology. An actual missing feature is not reconstructed.
    ///
    /// # Errors
    /// Invalid options or inconsistent source polygon references.
    pub fn normalize_coincident_sheets(
        &self,
        options: EvidenceOptions,
    ) -> Result<Option<NormalizedSourceEvidence>, MeshRecoveryError> {
        validate_evidence_options(options)?;
        if self.non_manifold_edge_count == 0 {
            return Ok(None);
        }
        let polygons = (0..self.source_planar_support.len())
            .map(|id| {
                let ids = self
                    .source_polygon_vertices(id)
                    .ok_or(MeshRecoveryError::DegeneratePolygon(id))?;
                input_polygon_from_ids(&self.positions, ids.to_vec(), id, KERNEL_EPSILON)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
        let mut groups = BTreeMap::<SheetKey, Vec<Sheet>>::new();
        let incident_supports = incidence
            .values()
            .filter(|uses| uses.len() > 2)
            .flat_map(|uses| {
                uses.iter()
                    .map(|&(polygon, _)| self.source_planar_support[polygon])
            })
            .collect::<BTreeSet<_>>();
        for support in incident_supports {
            if let Some((key, sheet)) = sheet_key(self, support, &polygons, &incidence) {
                groups.entry(key).or_default().push(sheet);
            }
        }
        let mut canceled_sheets = Vec::new();
        let mut removed = BTreeSet::new();
        for sheets in groups.into_values() {
            let [a, b] = sheets.as_slice() else {
                continue;
            };
            if !a.directions.iter().zip(&b.directions).all(|(a, b)| a != b) {
                continue;
            }
            let first = &self.planar_supports[a.support];
            let second = &self.planar_supports[b.support];
            if first.normal.dot(second.normal) >= 0.0 {
                continue;
            }
            // Exact common boundary samples anchor the plane. Every
            // interior vertex of either triangulation must also lie on
            // it within numerical resolution; nearby displaced sheets
            // cannot be paired by this operation.
            let residual = a
                .vertices
                .iter()
                .chain(&b.vertices)
                .map(|&id| (self.positions[id] - first.origin).dot(first.normal).abs())
                .fold(0.0_f64, f64::max);
            if residual > KERNEL_EPSILON.min(options.planar_support.distance_tolerance) {
                continue;
            }
            let pair = first
                .source_face_indices
                .iter()
                .chain(&second.source_face_indices)
                .copied()
                .collect::<BTreeSet<_>>();
            let Some(resolved_edges) = cancellation_edges(&polygons, &incidence, &pair, &removed)
            else {
                continue;
            };
            removed.extend(pair);
            canceled_sheets.push(CanceledSourceSheet {
                sides: [
                    first.source_face_indices.clone(),
                    second.source_face_indices.clone(),
                ],
                resolved_edges,
                max_plane_residual: residual,
            });
        }
        if canceled_sheets.is_empty() {
            return Ok(None);
        }
        let (original_polygon_ids, retained): (Vec<_>, Vec<_>) = polygons
            .into_iter()
            .enumerate()
            .filter(|(id, _)| !removed.contains(id))
            .unzip();
        let evidence = analyze_welded_evidence(self.positions.clone(), retained, options);
        Ok(Some(NormalizedSourceEvidence {
            evidence,
            original_polygon_ids,
            canceled_sheets,
        }))
    }
}

// For a compact oriented planar chain its directed boundary determines
// its signed coverage. Require a consistently oriented manifold interior,
// the exact reverse boundary, and coplanarity of every vertex; this
// also handles differing triangulations of the same planar sheet.
fn sheet_key(
    evidence: &MeshEvidence,
    support: usize,
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
) -> Option<(SheetKey, Sheet)> {
    let region = &evidence.planar_supports[support];
    let mut vertices = BTreeSet::new();
    let mut edges = BTreeMap::<EdgeKey, Vec<bool>>::new();
    for &id in &region.source_face_indices {
        let ids = &polygons[id].vertices;
        vertices.extend(ids.iter().copied());
        for (&a, &b) in ids.iter().zip(ids.iter().cycle().skip(1)) {
            let (key, direction) = edge_key(a, b);
            edges.entry(key).or_default().push(direction);
        }
    }
    let mut boundary = Vec::new();
    let mut directions = Vec::new();
    for (key, uses) in edges {
        match uses.as_slice() {
            [forward] => {
                boundary.push(key);
                directions.push(*forward);
            }
            [a, b] if a != b => {}
            _ => return None,
        }
    }
    if boundary.is_empty()
        || !boundary
            .iter()
            .any(|key| incidence.get(key).is_some_and(|u| u.len() > 2))
    {
        return None;
    }
    Some((
        SheetKey { boundary },
        Sheet {
            support,
            directions,
            vertices: vertices.into_iter().collect(),
        },
    ))
}

// Test the entire affected incidence transaction. No additional geometric
// boundaries or nonmanifold contacts may be introduced by cancellation.
fn cancellation_edges(
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    pair: &BTreeSet<usize>,
    removed: &BTreeSet<usize>,
) -> Option<Vec<EdgeKey>> {
    let keys = pair
        .iter()
        .flat_map(|&id| {
            let ids = &polygons[id].vertices;
            ids.iter()
                .zip(ids.iter().cycle().skip(1))
                .map(|(&a, &b)| edge_key(a, b).0)
        })
        .collect::<BTreeSet<_>>();
    let mut resolved = Vec::new();
    for key in keys {
        let before = incidence
            .get(&key)?
            .iter()
            .filter(|(id, _)| !removed.contains(id))
            .copied()
            .collect::<Vec<_>>();
        let after = before
            .iter()
            .filter(|(id, _)| !pair.contains(id))
            .copied()
            .collect::<Vec<_>>();
        match after.as_slice() {
            [] if before.len() >= 2 => {}
            [single] if before.as_slice() == [*single] => {}
            [a, b] if a.1 != b.1 => {
                if before.len() > 2 {
                    resolved.push(key);
                }
            }
            _ => return None,
        }
        if after.is_empty() && before.len() > 2 {
            resolved.push(key);
        }
    }
    (!resolved.is_empty()).then_some(resolved)
}

#[cfg(test)]
mod tests {
    use super::super::analyze_mesh_evidence;
    use super::*;
    use crate::Point3;
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    #[test]
    fn cancels_differently_triangulated_opposite_sheet_attached_to_a_plane()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, -1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                    Point3::new(1.0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([
                [0, 1, 2],
                [1, 0, 3],
                [0, 1, 4],
                [1, 5, 4],
                [1, 0, 5],
                [0, 4, 5],
            ]),
        )?;
        let options = EvidenceOptions::default();
        let source = analyze_mesh_evidence(&mesh, options)?;
        assert_eq!(source.non_manifold_edge_count, 1);
        let result = source
            .normalize_coincident_sheets(options)?
            .ok_or("cancel pair")?;
        assert_eq!(result.original_polygon_ids, vec![0, 1]);
        assert_eq!(result.evidence.non_manifold_edge_count, 0);
        assert_eq!(result.evidence.planar_supports.len(), 1);
        assert_eq!(result.canceled_sheets.len(), 1);
        assert_eq!(source.source_planar_support.len(), 6);
        Ok(())
    }

    #[test]
    fn retains_ambiguous_noncoincident_faces_and_unattached_reverse_sheets()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, -1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                    Point3::new(0.0, 0.0, -1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [1, 0, 3], [0, 1, 4], [1, 0, 5]]),
        )?;
        let options = EvidenceOptions::default();
        let source = analyze_mesh_evidence(&mesh, options)?;
        assert!(source.normalize_coincident_sheets(options)?.is_none());
        Ok(())
    }
}
