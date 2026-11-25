//! Assembly of source facet approximations with one shared edge per
//! source incidence. External analytic trims and their vertices are
//! retained by identity, not copied into independent approximations.

use super::{
    BoundaryFacetError, BoundaryFacetInput, SurfaceRegionBudget, edge_key, source_edge_incidence,
};
use crate::{Edge, InnerSpace, Point3, Shell, Vertex, builder};
mod coalesce;

use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    ops::Range,
};

/// A native shared edge associated with an oriented source segment.
#[derive(Debug, Clone)]
pub struct PrescribedSourceEdge {
    /// Source vertex indices, in the edge's orientation.
    pub vertices: [usize; 2],
    /// Exact edge whose topological identity is retained.
    pub edge: Edge,
}

/// Bounded approximation of a complete connected source facet patch.
#[derive(Debug)]
pub struct BoundaryPatchApproximation {
    /// Oriented native shell, open or closed according to source incidence.
    pub shell: Shell,
    /// Native face index range for each input polygon, in source order.
    /// Polygons certified together as one cell share the same range.
    pub source_face_ranges: Vec<Range<usize>>,
    /// Number of reconstructed cells covering multiple original polygons.
    pub coalesced_cells: usize,
    /// Largest displacement of an original source vertex.
    pub vertex_displacement: f64,
    /// Maximum continuous two-sided deviation of a reconstructed facet
    /// from its source centroid fan.
    pub deviation_upper_bound: f64,
    /// Maximum bounded radial seam discrepancy over all reconstructed facets.
    pub seam_deviation_upper_bound: f64,
}

/// Numerical and local source-cell search policy for a patch.
#[derive(Clone, Copy, Debug)]
pub struct BoundaryPatchOptions<'a> {
    /// Independent vertex displacement and continuous chord bounds.
    pub budget: SurfaceRegionBudget,
    /// Optional original region ownership: a cell may not cross owners.
    pub owner_by_polygon: Option<&'a [usize]>,
    /// Maximum connected-cell union hypotheses for the entire patch.
    /// Zero disables local coalescence.
    pub max_coalescence_trials: usize,
}

impl BoundaryPatchOptions<'_> {
    /// Standard bounded search for an unpartitioned source patch.
    #[must_use]
    pub const fn new(budget: SurfaceRegionBudget) -> Self {
        Self {
            budget,
            owner_by_polygon: None,
            max_coalescence_trials: 1024,
        }
    }
}

/// The patch could not be constructed with the requested topology and budgets.
#[derive(Debug, Clone, thiserror::Error)]
pub enum BoundaryPatchError {
    /// Source or prescribed edge incidence is inconsistent.
    #[error("invalid boundary patch incidence")]
    InvalidIncidence,
    /// Conflicting topological identities or excessive source vertex movement.
    #[error("incompatible boundary patch junction")]
    InvalidJunction,
    /// The requested local facet approximation cannot be established.
    #[error("source facet {polygon} approximation failed: {error}")]
    Facet {
        /// Index of the failing input polygon.
        polygon: usize,
        /// Local geometric rejection.
        error: BoundaryFacetError,
        /// Original source vertex IDs in polygon boundary order.
        source_vertex_ids: Vec<usize>,
        /// Exact geometric input needed to replay this failed certificate.
        input: Box<BoundaryFacetInput>,
    },
    /// Native output failed connected, oriented manifold validation.
    #[error("invalid boundary patch shell")]
    InvalidShell,
}

fn bind_vertex(
    vertices: &mut BTreeMap<usize, Vertex>,
    index: usize,
    vertex: &Vertex,
) -> Result<(), BoundaryPatchError> {
    use std::collections::btree_map::Entry;
    match vertices.entry(index) {
        Entry::Vacant(entry) => {
            entry.insert(vertex.clone());
        }
        Entry::Occupied(entry) => {
            if entry.get() != vertex {
                return Err(BoundaryPatchError::InvalidJunction);
            }
        }
    }
    Ok(())
}

fn source_vertices(
    source: &[Point3],
    polygons: &[&[usize]],
) -> Result<BTreeSet<usize>, BoundaryPatchError> {
    let mut used = BTreeSet::new();
    for polygon in polygons {
        let mut local = BTreeSet::new();
        if polygon.len() < 3
            || polygon.iter().any(|&index| {
                source
                    .get(index)
                    .is_none_or(|p| ![p.x, p.y, p.z].into_iter().all(f64::is_finite))
                    || !local.insert(index)
            })
        {
            return Err(BoundaryPatchError::InvalidIncidence);
        }
        used.extend(local);
    }
    if used.is_empty() {
        return Err(BoundaryPatchError::InvalidIncidence);
    }
    Ok(used)
}

/// Construct a connected native approximation from source polygons and
/// exact prescribed source edges.
///
/// The original source topology supplies
/// all unspecified straight edges, including interior edges. A source
/// vertex must have one native identity throughout the whole patch.
///
/// The input facets must admit the certified projected convex
/// fan used by [`super::approximate_boundary_facet`]. This does not infer the
/// original CAD control net. If a facet cannot be certified, a bounded
/// local union search can remove unprescribed internal diagonals, while
/// retaining every source vertex and region owner. Deviation is measured
/// against the original polygonal source union. Each prescribed edge must correspond to one actual
/// source segment, in either orientation. Vertex and facet approximation
/// budgets are independent.
///
/// # Errors
/// Rejects non-manifold, disconnected or inconsistently oriented incidence,
/// duplicate/conflicting edge and vertex identities, inadmissible local
/// approximations and geometrically inconsistent native output.
pub fn approximate_boundary_patch(
    source: &[Point3],
    polygons: &[&[usize]],
    prescribed: &[PrescribedSourceEdge],
    options: BoundaryPatchOptions<'_>,
) -> Result<BoundaryPatchApproximation, BoundaryPatchError> {
    use BoundaryPatchError::{InvalidIncidence, InvalidJunction};
    let budget = options.budget;
    if options
        .owner_by_polygon
        .is_some_and(|owners| owners.len() != polygons.len())
    {
        return Err(InvalidIncidence);
    }
    if !budget.vertex.is_finite()
        || budget.vertex < 0.0
        || !budget.chord.is_finite()
        || budget.chord <= 0.0
    {
        return Err(InvalidJunction);
    }
    let used = source_vertices(source, polygons)?;
    let incidence = source_edge_incidence(polygons.iter().copied());
    if incidence
        .values()
        .any(|uses| uses.len() > 2 || (uses.len() == 2 && uses[0].1 == uses[1].1))
    {
        return Err(InvalidIncidence);
    }
    let expected_closed = incidence.values().all(|uses| uses.len() == 2);
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    let mut protected = BTreeSet::new();
    for prescribed in prescribed {
        let [a, b] = prescribed.vertices;
        let (key, forward) = edge_key(a, b);
        if !incidence.contains_key(&key)
            || !prescribed.edge.is_geometric_consistent()
            || edges.contains_key(&key)
        {
            return Err(InvalidIncidence);
        }
        protected.insert(key);
        bind_vertex(&mut vertices, a, prescribed.edge.front())?;
        bind_vertex(&mut vertices, b, prescribed.edge.back())?;
        edges.insert(
            key,
            if forward {
                prescribed.edge.clone()
            } else {
                prescribed.edge.inverse()
            },
        );
    }
    let mut vertex_displacement = 0.0_f64;
    let mut identities = HashSet::new();
    for index in used {
        let vertex = vertices
            .entry(index)
            .or_insert_with(|| builder::vertex(source[index]));
        let movement = (vertex.point() - source[index]).magnitude();
        if !movement.is_finite() || movement > budget.vertex || !identities.insert(vertex.id()) {
            return Err(InvalidJunction);
        }
        vertex_displacement = vertex_displacement.max(movement);
    }
    for key in incidence.keys().copied() {
        edges
            .entry(key)
            .or_insert_with(|| builder::line(&vertices[&key[0]], &vertices[&key[1]]));
    }
    coalesce::AssemblyInput {
        source,
        polygons,
        edges: &edges,
        incidence: &incidence,
        protected: &protected,
        options,
        vertex_displacement,
        expected_closed,
    }
    .assemble()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BsplineCurve, Curve, KnotVector, NurbsCurve, Vector4};

    #[test]
    fn sliver_boundary_recovers_with_certified_convex_source_union()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.005, 0.0, 0.0),
            Point3::new(1.327, 0.00024, 0.0),
            Point3::new(1.327, 0.1, 0.0),
        ];
        let [a, b, c, _d] = points.map(builder::vertex);
        let curve = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 1.0),
                Vector4::new(0.0025, 0.0001, 0.0, 1.0),
                Vector4::new(0.005, 0.0, 0.0, 1.0),
            ],
        )?));
        let edge = Edge::new(&a, &b, curve)?;
        let prescribed = PrescribedSourceEdge {
            vertices: [0, 1],
            edge: edge.clone(),
        };
        let polygons: [&[usize]; 2] = [&[0, 1, 2], &[0, 2, 3]];
        let budget = SurfaceRegionBudget {
            vertex: 0.0,
            chord: 0.001,
        };
        let options = BoundaryPatchOptions::new(budget);
        assert!(matches!(
            approximate_boundary_patch(
                &points,
                &polygons,
                std::slice::from_ref(&prescribed),
                BoundaryPatchOptions {
                    max_coalescence_trials: 0,
                    ..options
                }
            ),
            Err(BoundaryPatchError::Facet { polygon: 0, .. })
        ));
        let merged = approximate_boundary_patch(
            &points,
            &polygons,
            std::slice::from_ref(&prescribed),
            options,
        )?;
        assert_eq!(merged.coalesced_cells, 1);
        assert_eq!(merged.source_face_ranges, vec![0..4, 0..4]);
        assert!(merged.deviation_upper_bound <= budget.chord);
        assert_eq!(merged.vertex_displacement, 0.0);
        assert_eq!(
            merged
                .shell
                .edge_iter()
                .filter(|e| e.id() == edge.id())
                .count(),
            1
        );
        assert!(merged.shell.is_geometric_consistent());
        assert!(matches!(
            approximate_boundary_patch(
                &points,
                &polygons,
                std::slice::from_ref(&prescribed),
                BoundaryPatchOptions {
                    owner_by_polygon: Some(&[0, 1]),
                    ..options
                }
            ),
            Err(BoundaryPatchError::Facet { polygon: 0, .. })
        ));
        let diagonal = PrescribedSourceEdge {
            vertices: [0, 2],
            edge: builder::line(&a, &c),
        };
        assert!(matches!(
            approximate_boundary_patch(&points, &polygons, &[prescribed, diagonal], options),
            Err(BoundaryPatchError::Facet { polygon: 0, .. })
        ));
        Ok(())
    }

    #[test]
    fn connected_patch_preserves_shared_curve_and_vertex_identity()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let a = builder::vertex(points[0]);
        let c = builder::vertex(points[2]);
        let diagonal = Edge::new(
            &a,
            &c,
            Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
                KnotVector::bezier_knot(2),
                vec![
                    Vector4::new(0.0, 0.0, 0.0, 1.0),
                    Vector4::new(1.0, 1.0, 0.1, 1.0),
                    Vector4::new(2.0, 2.0, 0.0, 1.0),
                ],
            )?)),
        )?;
        let polygons: [&[usize]; 2] = [&[0, 1, 2], &[0, 2, 3]];
        let prescribed = [PrescribedSourceEdge {
            vertices: [0, 2],
            edge: diagonal.clone(),
        }];
        let budget = SurfaceRegionBudget {
            vertex: 0.0,
            chord: 0.15,
        };
        let recovered = approximate_boundary_patch(
            &points,
            &polygons,
            &prescribed,
            BoundaryPatchOptions::new(budget),
        )?;
        assert_eq!(recovered.source_face_ranges, vec![0..3, 3..6]);
        assert_eq!(
            recovered
                .shell
                .edge_iter()
                .filter(|edge| edge.id() == diagonal.id())
                .count(),
            2
        );
        assert!(recovered.deviation_upper_bound <= budget.chord);
        assert_eq!(recovered.vertex_displacement, 0.0);
        assert!(
            approximate_boundary_patch(
                &points,
                &polygons,
                &prescribed,
                BoundaryPatchOptions::new(SurfaceRegionBudget {
                    vertex: 0.0,
                    chord: 0.001
                })
            )
            .is_err()
        );
        assert!(
            approximate_boundary_patch(
                &points,
                &[&[0, 1, 2], &[0, 3, 2]],
                &prescribed,
                BoundaryPatchOptions::new(budget)
            )
            .is_err()
        );
        let conflicting = PrescribedSourceEdge {
            vertices: [0, 1],
            edge: builder::line(&builder::vertex(points[0]), &builder::vertex(points[1])),
        };
        assert!(matches!(
            approximate_boundary_patch(
                &points,
                &polygons,
                &[prescribed[0].clone(), conflicting],
                BoundaryPatchOptions::new(budget)
            ),
            Err(BoundaryPatchError::InvalidJunction)
        ));
        Ok(())
    }
}
