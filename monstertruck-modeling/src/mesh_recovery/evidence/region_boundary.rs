//! Winding-preserving boundaries of arbitrary unions of original supports.
//! These polygonal loops carry the provenance needed to fit shared trims.

use super::super::{
    EdgeIncidence, edge_key, source_edge_incidence, trace_oriented_boundary_cycles,
};
use super::MeshEvidence;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

mod contacts;
mod containment;
mod curve_evidence;
pub use contacts::{RegionContactGraph, RegionContactPath};
pub use containment::{CarrierResidual, SurfaceCarrier, SurfaceCurveAgreement};
pub use curve_evidence::{ContactCircleEvidence, ContactCurveEvidence, ContactLineEvidence};

/// One directed source edge on a candidate region's boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RegionBoundaryEdge {
    /// Original canonical position IDs, in the interior polygon's winding.
    pub vertices: [usize; 2],
    /// Original source polygon on the selected side.
    pub interior_polygon: usize,
    /// Original source polygon across the edge, or none on an open mesh.
    pub exterior_polygon: Option<usize>,
}

/// Exact polygonal boundary of a union of source planar supports.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionBoundary {
    /// Oriented loops, starting at their lexicographically smallest directed
    /// canonical source edge.
    /// No outer/inner ordering is inferred for a nonplanar surface.
    /// Empty for a closed region; not a certificate of manifoldness.
    /// Distinct sectors can share a canonical source vertex; see
    /// `Self::pinched_vertices`.
    pub loops: Vec<Vec<RegionBoundaryEdge>>,
    /// Number of distinct selected original source polygons.
    pub source_polygon_count: usize,
    /// Canonical source vertex IDs shared by distinct boundary sectors.
    /// A trim crossing or touching itself at a vertex is not a simple
    /// analytic face boundary, even if all its directed edge cycles close.
    pub pinched_vertices: Vec<usize>,
    /// Canonically sorted, unique original support IDs in this region.
    pub support_indices: Vec<usize>,
}

/// One bounded exterior polygon fan between selected sectors at a pinch.
/// Adding the entire support union is only a proposal for geometric testing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionBoundaryGap {
    /// Original vertex shared by the two selected sectors.
    pub vertex: usize,
    /// Original exterior polygons encountered in fan order.
    pub source_polygons: Vec<usize>,
    /// Unique, sorted IDs of all source planar supports containing the fan.
    pub added_support_indices: Vec<usize>,
}

/// Reason a candidate union cannot produce an unambiguous oriented boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum RegionBoundaryError {
    /// At least one support must be selected.
    #[error("no source supports selected")]
    EmptyRegion,
    /// A selected support ID is not present in the source evidence.
    #[error("unknown source support {0}")]
    InvalidSupport(usize),
    /// Source evidence contains inconsistent polygon or edge references.
    #[error("inconsistent source evidence")]
    InvalidEvidence,
    /// A partition must assign a region to every original source polygon.
    #[error("source partition must assign every original polygon")]
    InvalidPartition,
    /// A selected polygon touches a nonmanifold or winding-conflicting edge.
    #[error("ambiguous source edge {0:?}")]
    AmbiguousEdge([usize; 2]),
    /// Source face sectors fail to form closed directed boundary cycles.
    #[error("region boundary sectors do not form closed oriented cycles")]
    NonCyclicBoundary,
}

/// Reusable full source incidence for arbitrary overlapping region hypotheses.
/// Includes suppressed tessellation diagonals and ambiguous source contacts.
pub struct RegionBoundaryIndex<'a> {
    evidence: &'a MeshEvidence,
    incidence: EdgeIncidence,
}

impl<'a> RegionBoundaryIndex<'a> {
    /// Index original polygon incidence once for any number of region queries.
    ///
    /// # Errors
    /// Returns an error for inconsistent original polygon references.
    pub fn new(evidence: &'a MeshEvidence) -> Result<Self, RegionBoundaryError> {
        let polygons = (0..evidence.source_planar_support.len())
            .map(|id| {
                evidence
                    .source_polygon_vertices(id)
                    .ok_or(RegionBoundaryError::InvalidEvidence)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if polygons.iter().any(|p| p.len() < 3) {
            return Err(RegionBoundaryError::InvalidEvidence);
        }
        Ok(Self {
            evidence,
            incidence: source_edge_incidence(polygons),
        })
    }

    /// Extract the exact oriented source boundary of the selected support union.
    /// Selection order and duplicate support IDs do not affect the result.
    ///
    /// # Errors
    /// Returns an error for invalid selection, ambiguous source incidence, or
    /// a boundary whose original polygon sectors fail to close. Pinches
    /// are retained explicitly rather than choosing a geometric successor.
    pub fn boundary(&self, supports: &[usize]) -> Result<RegionBoundary, RegionBoundaryError> {
        if supports.is_empty() {
            return Err(RegionBoundaryError::EmptyRegion);
        }
        let mut selected = supports.to_vec();
        selected.sort_unstable();
        selected.dedup();
        let mut boundary = Vec::new();
        let mut source_polygon_count = 0;
        for &support in &selected {
            let source = self
                .evidence
                .planar_supports
                .get(support)
                .ok_or(RegionBoundaryError::InvalidSupport(support))?;
            for &polygon in &source.source_face_indices {
                if self.evidence.source_planar_support.get(polygon) != Some(&support) {
                    return Err(RegionBoundaryError::InvalidEvidence);
                }
                source_polygon_count += 1;
                let vertices = self
                    .evidence
                    .source_polygon_vertices(polygon)
                    .ok_or(RegionBoundaryError::InvalidEvidence)?;
                for (&a, &b) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
                    let exterior_polygon = self.opposite_polygon(polygon, [a, b])?;
                    if let Some(other) = exterior_polygon {
                        let other_support = self
                            .evidence
                            .source_planar_support
                            .get(other)
                            .ok_or(RegionBoundaryError::InvalidEvidence)?;
                        if selected.binary_search(other_support).is_ok() {
                            continue;
                        }
                    }
                    let edge = RegionBoundaryEdge {
                        vertices: [a, b],
                        interior_polygon: polygon,
                        exterior_polygon,
                    };
                    boundary.push(edge);
                }
            }
        }
        self.trace_boundary(boundary, source_polygon_count, selected)
    }

    fn trace_boundary(
        &self,
        mut boundary: Vec<RegionBoundaryEdge>,
        source_polygon_count: usize,
        support_indices: Vec<usize>,
    ) -> Result<RegionBoundary, RegionBoundaryError> {
        boundary.sort_unstable_by_key(|edge| edge.vertices);
        let edge_ids = boundary
            .iter()
            .enumerate()
            .map(|(id, edge)| (edge.vertices, id))
            .collect::<HashMap<_, _>>();
        if edge_ids.len() != boundary.len() {
            return Err(RegionBoundaryError::NonCyclicBoundary);
        }
        let mut outgoing = HashMap::default();
        let mut pinches = HashSet::default();
        for (id, edge) in boundary.iter().enumerate() {
            let start = edge.vertices[0];
            if outgoing.insert(start, id).is_some() {
                pinches.insert(start);
            }
        }

        // A unique outgoing edge is already the sector successor. At a
        // pinch, follow the original incident polygon fan; proximity or
        // angular turn would incorrectly connect unrelated source sectors.
        let next = boundary
            .iter()
            .enumerate()
            .map(|(id, edge)| {
                let vertex = edge.vertices[1];
                let successor = if pinches.contains(&vertex) {
                    self.sector_successor(edge, &boundary, &edge_ids)?
                } else {
                    *outgoing
                        .get(&vertex)
                        .ok_or(RegionBoundaryError::NonCyclicBoundary)?
                };
                Ok((id, successor))
            })
            .collect::<Result<HashMap<_, _>, RegionBoundaryError>>()?;
        let mut pinched_vertices = pinches.into_iter().collect::<Vec<_>>();
        pinched_vertices.sort_unstable();
        let cycles =
            trace_oriented_boundary_cycles(next).ok_or(RegionBoundaryError::NonCyclicBoundary)?;
        let loops = cycles
            .into_iter()
            .map(|cycle| cycle.into_iter().map(|id| boundary[id]).collect())
            .collect();
        Ok(RegionBoundary {
            loops,
            source_polygon_count,
            pinched_vertices,
            support_indices,
        })
    }

    fn sector_successor(
        &self,
        edge: &RegionBoundaryEdge,
        boundary: &[RegionBoundaryEdge],
        edge_ids: &HashMap<[usize; 2], usize>,
    ) -> Result<usize, RegionBoundaryError> {
        let vertex = edge.vertices[1];
        let mut polygon = edge.interior_polygon;
        for _ in 0..self.evidence.source_planar_support.len() {
            let directed = self.outgoing_edge(polygon, vertex)?;
            if let Some(&id) = edge_ids.get(&directed) {
                if boundary[id].interior_polygon != polygon {
                    return Err(RegionBoundaryError::InvalidEvidence);
                }
                return Ok(id);
            }
            polygon = self
                .opposite_polygon(polygon, directed)?
                .ok_or(RegionBoundaryError::NonCyclicBoundary)?;
        }
        Err(RegionBoundaryError::NonCyclicBoundary)
    }

    fn outgoing_edge(
        &self,
        polygon: usize,
        vertex: usize,
    ) -> Result<[usize; 2], RegionBoundaryError> {
        let vertices = self
            .evidence
            .source_polygon_vertices(polygon)
            .ok_or(RegionBoundaryError::InvalidEvidence)?;
        let i = vertices
            .iter()
            .position(|&v| v == vertex)
            .ok_or(RegionBoundaryError::InvalidEvidence)?;
        Ok([vertex, vertices[(i + 1) % vertices.len()]])
    }

    fn opposite_polygon(
        &self,
        polygon: usize,
        directed: [usize; 2],
    ) -> Result<Option<usize>, RegionBoundaryError> {
        let (key, forward) = edge_key(directed[0], directed[1]);
        let uses = self
            .incidence
            .get(&key)
            .ok_or(RegionBoundaryError::InvalidEvidence)?;
        match uses.as_slice() {
            [only] if *only == (polygon, forward) => Ok(None),
            [a, b] if *a == (polygon, forward) && b.1 != forward => Ok(Some(b.0)),
            [a, b] if *b == (polygon, forward) && a.1 != forward => Ok(Some(a.0)),
            _ => Err(RegionBoundaryError::AmbiguousEdge(key)),
        }
    }

    /// Enumerate bounded exterior source fans at pinched vertices.
    ///
    /// This proposes adding complete observed supports, never creating source
    /// vertices or changing mesh topology. Open or ambiguous exterior fans do
    /// not establish bounded gaps. A proposed union requires independent
    /// whole-region geometric and topological validation.
    #[must_use = "gap proposals require independent validation"]
    pub fn pinch_gaps(&self, boundary: &RegionBoundary) -> Vec<RegionBoundaryGap> {
        boundary
            .loops
            .iter()
            .flatten()
            .filter(|edge| {
                boundary
                    .pinched_vertices
                    .binary_search(&edge.vertices[0])
                    .is_ok()
            })
            .filter_map(|edge| self.bounded_gap(edge, &boundary.support_indices))
            .collect()
    }

    fn bounded_gap(
        &self,
        edge: &RegionBoundaryEdge,
        selected: &[usize],
    ) -> Option<RegionBoundaryGap> {
        let vertex = edge.vertices[0];
        let mut polygon = edge.exterior_polygon?;
        let mut source_polygons = Vec::new();
        for _ in 0..self.evidence.source_planar_support.len() {
            let support = *self.evidence.source_planar_support.get(polygon)?;
            if selected.binary_search(&support).is_ok() {
                if source_polygons.is_empty() || polygon == edge.interior_polygon {
                    return None;
                }
                let mut added_support_indices = source_polygons
                    .iter()
                    .map(|&id| self.evidence.source_planar_support[id])
                    .collect::<Vec<_>>();
                added_support_indices.sort_unstable();
                added_support_indices.dedup();
                return Some(RegionBoundaryGap {
                    vertex,
                    source_polygons,
                    added_support_indices,
                });
            }
            source_polygons.push(polygon);
            let directed = self.outgoing_edge(polygon, vertex).ok()?;
            polygon = self.opposite_polygon(polygon, directed).ok()??;
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Point3,
        mesh_recovery::{EvidenceOptions, analyze_mesh_evidence},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    fn two_faces() -> Result<MeshEvidence, Box<dyn std::error::Error>> {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [1, 0, 3]]),
        )?;
        Ok(analyze_mesh_evidence(&mesh, EvidenceOptions::default())?)
    }

    #[test]
    fn adjacent_regions_share_reversed_boundary_with_polygon_provenance()
    -> Result<(), Box<dyn std::error::Error>> {
        let evidence = two_faces()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let a = index.boundary(&[0])?;
        let b = index.boundary(&[1])?;
        assert_eq!(a.loops.len(), 1);
        assert_eq!(b.loops.len(), 1);
        let seam_a = a
            .loops
            .iter()
            .flatten()
            .find(|edge| edge.exterior_polygon.is_some())
            .ok_or("missing seam a")?;
        let seam_b = b
            .loops
            .iter()
            .flatten()
            .find(|edge| edge.exterior_polygon.is_some())
            .ok_or("missing seam b")?;
        assert_eq!(seam_a.vertices, [seam_b.vertices[1], seam_b.vertices[0]]);
        assert_eq!(Some(seam_a.interior_polygon), seam_b.exterior_polygon);
        assert_eq!(Some(seam_b.interior_polygon), seam_a.exterior_polygon);
        let joined = index.boundary(&[1, 0, 1])?;
        assert_eq!(joined, index.boundary(&[0, 1])?);
        assert_eq!(joined.source_polygon_count, 2);
        assert_eq!(joined.loops.len(), 1);
        assert_eq!(joined.loops[0].len(), 4);
        assert!(
            joined
                .loops
                .iter()
                .flatten()
                .all(|edge| edge.exterior_polygon.is_none())
        );
        Ok(())
    }

    #[test]
    fn exactly_welded_vertex_contact_preserves_distinct_source_sectors()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(-1.0, 0.0, 0.0),
                    Point3::new(0.0, -1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 3, 4]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let region = index.boundary(&[0, 1])?;
        assert!(index.pinch_gaps(&region).is_empty());
        assert_eq!(region.source_polygon_count, 2);
        assert_eq!(region.loops.len(), 2);
        assert_eq!(region.pinched_vertices.len(), 1);
        for cycle in &region.loops {
            assert_eq!(cycle.len(), 3);
            assert!(
                cycle
                    .iter()
                    .all(|edge| edge.interior_polygon == cycle[0].interior_polygon)
            );
        }
        Ok(())
    }

    #[test]
    fn bounded_polygon_fans_offer_complete_source_support_unions()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 1.0),
                    Point3::new(0.5, 1.0, 1.0),
                    Point3::new(-0.5, 1.0, 1.0),
                    Point3::new(-1.0, 0.0, 1.0),
                    Point3::new(-0.5, -1.0, 1.0),
                    Point3::new(0.5, -1.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([
                [0, 1, 2],
                [0, 2, 3],
                [0, 3, 4],
                [0, 4, 5],
                [0, 5, 6],
                [0, 6, 1],
            ]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let selected = [
            evidence.source_planar_support[0],
            evidence.source_planar_support[3],
        ];
        let region = index.boundary(&selected)?;
        let gaps = index.pinch_gaps(&region);
        assert_eq!(gaps.len(), 2);
        assert_eq!(region.pinched_vertices.len(), 1);
        let mut actual_polygons = gaps
            .iter()
            .map(|gap| {
                let mut ids = gap.source_polygons.clone();
                ids.sort_unstable();
                ids
            })
            .collect::<Vec<_>>();
        actual_polygons.sort_unstable();
        assert_eq!(actual_polygons, vec![vec![1, 2], vec![4, 5]]);
        for gap in gaps {
            let mut union = selected.to_vec();
            union.extend(gap.added_support_indices);
            let completed = index.boundary(&union)?;
            assert!(completed.pinched_vertices.is_empty());
            assert_eq!(completed.source_polygon_count, 4);
            assert_eq!(completed.loops.len(), 1);
        }
        Ok(())
    }

    #[test]
    fn invalid_selection_is_explicit() -> Result<(), Box<dyn std::error::Error>> {
        let evidence = two_faces()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        assert_eq!(index.boundary(&[]), Err(RegionBoundaryError::EmptyRegion));
        assert_eq!(
            index.boundary(&[2]),
            Err(RegionBoundaryError::InvalidSupport(2))
        );
        Ok(())
    }
}
