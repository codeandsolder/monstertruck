//! Primitive-independent geometric evidence from polygon mesh incidence.
//!
//! Source incidence and collinear chains are topology/coordinate facts.
//! Dihedral selection and parallel groups are hypotheses with explicit
//! thresholds; they must not alter the source vertices or edge graph.

use super::{
    EdgeIncidence, InputPolygon, MeshRecoveryError, PlanarRecoveryOptions, input_polygon_from_ids,
    positions_and_polygons, seams, segment_patches, signed_plane_distance, source_edge_incidence,
};
use crate::{InnerSpace, KERNEL_EPSILON, Point3, Vector3};
use monstertruck_mesh::PolygonMesh;
use std::collections::{HashMap, HashSet};

mod circle_axes;
mod curvature;
mod growth;
mod join;
mod meridian;
mod region_boundary;
mod source_sheet;
mod source_topology;
mod spatial;
mod station_graph;
mod support;
pub use circle_axes::{
    CircleAxisDiscovery, CircleAxisEvidence, CircleAxisOptions, CircleAxisStation,
    SourceCircleObservation, discover_circle_axes,
};
pub use curvature::{CurvatureTrack, CurvatureTrackOptions, trace_curvature_tracks};
pub use growth::{
    CurvatureGrowthHypothesis, CurvatureGrowthIndex, CurvatureGrowthOptions, CurvatureGrowthStep,
};
pub use join::{SpatialTrackJoinEvidence, measure_curvature_spatial_join};
pub use meridian::{CircularMeridianEvidence, measure_circular_meridian};
pub use region_boundary::{
    CarrierCoincidenceEvidence, CarrierResidual, ContactCircleEvidence, ContactCurveEvidence,
    ContactLineEvidence, NativeRegionFace, NativeTrimError, NativeTrimGraph, NativeTrimOptions,
    OrientedTrim, ReconciledPartition, ReconciledPlane, ReconciliationOptions, RegionBoundary,
    RegionBoundaryEdge, RegionBoundaryError, RegionBoundaryGap, RegionBoundaryIndex,
    RegionContactGraph, RegionContactPath, RegionTrimBoundary, RegionTrimLoops, SharedAnalyticTrim,
    SharedTrimGeometry, SharedTrimGraph, SharedTrimOptions, SharedTrimPiece,
    SharedTrimSegmentation, SurfaceCarrier, SurfaceCurveAgreement, SurfaceRegionBudget,
    SurfaceRegionPartition, SurfaceRegionProposal,
};
pub use source_sheet::{CanceledSourceSheet, NormalizedSourceEvidence};
pub use source_topology::SeparatedSourceVertex;
pub use spatial::{
    AxialRadialStation, SpatialTrackEvidence, measure_curvature_spatial,
    measure_spatial_support_set,
};
pub use station_graph::{
    MeridianGraphOptions, MeridianStation, MeridianStationGraph, MeridianTransition,
    StationLinearCandidate, StationMeridianCandidate, StationMeridianDiscovery,
};
use support::build_support_interfaces;
pub use support::{SupportBoundaryTopology, SupportInterface};

/// Thresholds only for feature discovery and direction clustering.
#[derive(Clone, Copy, Debug)]
pub struct EvidenceOptions {
    /// Source edge is retained when its normal discontinuity exceeds this
    /// angle in radians. Coplanar triangulation edges are suppressed.
    pub min_dihedral_angle_radians: f64,
    /// Undirected straight chains differing by at most this angle in radians
    /// may join a candidate parallel direction group.
    pub max_parallel_angle_radians: f64,
    /// Seed-anchored angle and position bounds for connected planar supports.
    /// A one-triangle planar support is still evidence, not a CAD plane proof.
    pub planar_support: PlanarRecoveryOptions,
}

impl Default for EvidenceOptions {
    fn default() -> Self {
        Self {
            min_dihedral_angle_radians: 1.0e-6,
            max_parallel_angle_radians: 1.0e-6,
            planar_support: PlanarRecoveryOptions::default(),
        }
    }
}

/// Reason a source segment is geometrically interesting.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FeatureEdgeKind {
    /// Source surface boundary (one incident polygon).
    Boundary,
    /// Two oppositely oriented polygon edge uses have different normals.
    /// This may simply be a triangulated smooth surface.
    Dihedral,
    /// Two source polygon uses traverse the same edge in the same direction.
    /// The evidence remains observable but is not manifold B-rep topology.
    WindingConflict,
    /// More than two polygon uses coincide on a source edge. Frequently
    /// occurs where independent CAD bodies share geometrical contacts.
    NonManifoldJunction,
}

/// Source edge with canonical positions and original polygon provenance.
#[derive(Clone, Debug)]
pub struct FeatureEdge {
    /// Undirected canonical source position indices.
    pub vertices: [usize; 2],
    /// Original source polygon indices, retaining multi-use junctions.
    pub source_face_indices: Vec<usize>,
    /// Whether each source polygon traverses this canonical edge forward.
    pub source_forward_winding: Vec<bool>,
    /// Boundary, dihedral, or ambiguous source incidence.
    pub kind: FeatureEdgeKind,
    /// Dihedral angle in radians; none for source boundaries.
    pub dihedral_angle_radians: Option<f64>,
}

/// Maximal exact straight-through feature chain.
#[derive(Clone, Debug)]
pub struct CollinearChain {
    /// Ordered canonical vertices, retaining all original intermediate points.
    pub vertices: Vec<usize>,
    /// Indices in the feature-edge array between successive chain vertices.
    pub edge_indices: Vec<usize>,
    /// Category shared across the chain.
    pub kind: FeatureEdgeKind,
}

/// Relative line placement in one candidate parallel bundle's seed frame.
#[derive(Clone, Debug)]
pub struct ParallelPlacement {
    /// Index in [`MeshEvidence::chains`].
    pub chain_index: usize,
    /// Perpendicular displacement from seed line to this chain's first point.
    /// For approximately parallel directions this is a pointwise offset,
    /// not a globally constant separation.
    pub transverse_offset: Vector3,
    /// Magnitude of `transverse_offset` in mesh coordinate units.
    pub transverse_distance: f64,
    /// Projected start/end extent in the seed line's local frame.
    pub longitudinal_interval: [f64; 2],
    /// Angular difference from the *fixed seed direction*, in radians.
    pub angular_residual_radians: f64,
}

/// Angular hypothesis with spatial evidence for each straight-line member.
///
/// Every candidate is checked against the seed, never against an arbitrary
/// neighboring member; tolerance chains cannot drift through orientation
/// space and incorrectly merge nonparallel carriers.
#[derive(Clone, Debug)]
pub struct ParallelBundle {
    /// Fixed, canonical undirected unit seed axis.
    pub direction: Vector3,
    /// First point on the seed chain, providing a local reference frame.
    pub origin: Point3,
    /// Source chains and their relative placement measurements.
    pub placements: Vec<ParallelPlacement>,
}

impl ParallelBundle {
    /// Transverse separation of two members in the bundle's seed plane.
    /// For nearly parallel candidates, it is only a local measurement.
    #[must_use]
    pub fn transverse_separation(&self, first: usize, second: usize) -> Option<f64> {
        let a = self.placements.get(first)?.transverse_offset;
        let b = self.placements.get(second)?.transverse_offset;
        Some((a - b).magnitude())
    }
}

/// Relationship between two straight chain support lines, independent of a
/// high-level surface or CAD feature hypothesis.
#[derive(Clone, Debug)]
pub struct ChainPairEvidence {
    /// Undirected angle in radians, from parallel (zero) to perpendicular.
    pub angle_radians: f64,
    /// Adaptive exact 3D orientation predicate: all four stored source
    /// endpoints belong to a common supporting plane.
    pub exactly_coplanar: bool,
    /// Both source endpoints of the second chain lie on the first infinite
    /// line, according to adaptive exact 2D projected orientation predicates.
    pub exactly_coincident_carriers: bool,
    /// Closest distance between infinite support lines (not finite segments).
    /// For nearly parallel lines this may occur far outside either segment.
    pub supporting_line_distance: f64,
    /// Separation of the second line's first point from the first line;
    /// preserves a useful local spatial measurement for near-parallel lines.
    pub reference_transverse_distance: f64,
}

/// Connected seed-anchored approximately planar source facets.
///
/// This is an observed local polygonal support, not proof that the design
/// intended an analytic plane: a faceted cylinder also contains small planar
/// chord strips. Preserve small/isolated supports for later global fitting.
#[derive(Clone, Debug)]
pub struct PlanarSupport {
    /// Indices into the original source polygon sequence.
    pub source_face_indices: Vec<usize>,
    /// Oriented normal of the fixed seed polygon.
    pub normal: Vector3,
    /// Seed origin in the source coordinate frame.
    pub origin: Point3,
    /// Total original polygonal area in mesh coordinate units squared.
    pub area: f64,
    /// Maximum point-to-seed-plane residual among member vertices.
    pub max_vertex_plane_residual: f64,
    /// Maximum angular deviation from the fixed seed normal in radians.
    pub max_normal_angle_radians: f64,
}

/// Three observed planar supports whose central support touches both
/// neighbors along two source boundary components.
///
/// These measurements are useful for finding a smoothly rotating normal
/// sequence (for example a polygonal chord approximation of a curved surface)
/// without assuming a cylinder, torus or other source primitive.
#[derive(Clone, Debug)]
pub struct SupportNormalRotationEvidence {
    /// Middle support and its two neighbors, in the requested interface order.
    pub center_support: usize,
    pub neighbor_supports: [usize; 2],
    /// The two interface component IDs used by this local observation.
    pub interface_indices: [usize; 2],
    /// Canonical reference direction from the first interface's longest edge.
    pub reference_direction: Vector3,
    /// Difference between the two interfaces' undirected seed directions.
    pub interface_axis_residual_radians: f64,
    /// Actual source normal-turn angles across both interfaces.
    pub turn_radians: [f64; 2],
    /// Absolute mismatch between successive source-normal turn magnitudes.
    pub turn_magnitude_mismatch_radians: f64,
    /// Both normal rotations have the same signed sense around the reference
    /// axis, independent of the order in which the two interfaces are supplied.
    pub coherent_rotation: bool,
    /// Largest absolute dot product between reference axis and the normals
    /// of the three support seed planes; zero means all are transverse.
    pub max_normal_axis_dot: f64,
    /// Perpendicular distance from the second interface's first point to
    /// the first seed line; not the global distance for nonparallel axes.
    pub transverse_interface_separation: f64,
    /// Largest absolute departure from the seed line across both boundaries.
    pub max_interface_seed_line_offset: f64,
    /// Largest undirected source segment angular spread from either seed.
    pub max_interface_segment_spread_radians: f64,
    /// Candidate normal-space rotation axis inferred from the differences
    /// of all three consecutive source support normals. Unlike an interface
    /// tangent this can identify conical as well as cylindrical normal cones.
    /// None means the differences are rank-deficient: the axis is not
    /// identifiable from this triplet.
    pub normal_space_axis: Option<Vector3>,
    /// Sine of the angle between the two successive normal-difference vectors.
    /// Near zero means an ill-conditioned normal-space axis.
    pub normal_space_axis_strength: f64,
    /// Maximum undirected angular divergence between the inferred axis and
    /// the two observed boundary segment seed directions; None if the
    /// normal-space axis is unidentifiable.
    pub normal_space_axis_to_boundaries_radians: Option<f64>,
}

/// Reusable primitive-free geometric evidence for a source mesh.
#[derive(Clone, Debug)]
pub struct MeshEvidence {
    /// Exactly welded source positions, except separately certified closed
    /// component vertices that retain distinct topological IDs at contacts.
    pub positions: Vec<Point3>,
    /// Coordinate-identical copies created to preserve closed components.
    pub separated_vertices: Vec<SeparatedSourceVertex>,
    /// Boundary and dihedral source edges in canonical index order.
    pub edges: Vec<FeatureEdge>,
    /// Maximal exact collinear straight-through chains.
    pub chains: Vec<CollinearChain>,
    /// Approximate direction groups with measured spatial relationships.
    pub parallel_bundles: Vec<ParallelBundle>,
    /// Connected seed-anchored planar facet supports, including singletons.
    pub planar_supports: Vec<PlanarSupport>,
    /// Index into the planar-support list for each original source polygon.
    pub source_planar_support: Vec<usize>,
    /// Compressed source polygon positions after exact welding and T-junction
    /// insertion. Complete polygon vertices remain available for spatial
    /// verification, not just vertices on neighboring support interfaces.
    /// Offsets have one entry more than the source polygon count.
    pub source_polygon_vertex_offsets: Vec<usize>,
    /// Flattened canonical source-position IDs, indexed by the offsets.
    pub source_polygon_vertex_indices: Vec<usize>,
    /// Connected interfaces between different planar-support regions.
    /// These are source topology facts and measurements, not designed seams.
    pub support_interfaces: Vec<SupportInterface>,
    /// Interface IDs incident to each region, retaining disconnected
    /// interfaces between the same pair as distinct boundary components.
    pub support_interface_indices: Vec<Vec<usize>>,
    /// Locally coplanar triangulation edges not included as features.
    pub suppressed_interior_edge_count: usize,
    /// Retained edges used by more than two source polygons.
    pub non_manifold_edge_count: usize,
    /// Retained two-face edges whose oriented source uses conflict.
    pub inconsistent_winding_edge_count: usize,
}

impl MeshEvidence {
    /// Signed polyhedral volume of a closed, consistently oriented source chain.
    /// Returns none for open, ambiguous, or empty source incidence. This is
    /// a numerical integral, not a non-self-intersection certificate.
    #[must_use]
    pub fn closed_source_volume(&self) -> Option<f64> {
        if self
            .edges
            .iter()
            .any(|edge| edge.kind != FeatureEdgeKind::Dihedral)
            || self.non_manifold_edge_count != 0
            || self.inconsistent_winding_edge_count != 0
        {
            return None;
        }
        let polygons = (0..self.source_planar_support.len())
            .map(|id| self.source_polygon_vertices(id))
            .collect::<Option<Vec<_>>>()?;
        super::signed_polygon_volume(&self.positions, polygons)
    }

    /// Exact canonical source vertex IDs of an input polygon, including any
    /// collinear points inserted to resolve original T-junctions.
    #[must_use]
    pub fn source_polygon_vertices(&self, source_polygon: usize) -> Option<&[usize]> {
        let begin = *self.source_polygon_vertex_offsets.get(source_polygon)?;
        let end = *self.source_polygon_vertex_offsets.get(source_polygon + 1)?;
        self.source_polygon_vertex_indices.get(begin..end)
    }

    /// Refine selected coplanar support groups into individual original
    /// polygons while preserving all canonical source vertices, polygon IDs,
    /// edge incidence, and feature observations. Rebuild only support-level
    /// membership and interfaces; the original grouping remains available as
    /// an independent planar explanation in the caller's baseline evidence.
    ///
    /// # Errors
    /// Returns an error when the refinement references invalid support or
    /// source polygon IDs, or the original polygon geometry is inconsistent.
    pub fn refine_planar_supports(&self, split: &[usize]) -> Result<Self, MeshRecoveryError> {
        let split = split.iter().copied().collect::<HashSet<_>>();
        if split.iter().any(|&id| id >= self.planar_supports.len()) {
            return Err(MeshRecoveryError::InvalidEvidenceOptions);
        }
        let polygons = (0..self.source_planar_support.len())
            .map(|id| {
                let vertices = self
                    .source_polygon_vertices(id)
                    .ok_or(MeshRecoveryError::DegeneratePolygon(id))?;
                input_polygon_from_ids(&self.positions, vertices.to_vec(), id, KERNEL_EPSILON)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut groups = Vec::new();
        for (id, support) in self.planar_supports.iter().enumerate() {
            if split.contains(&id) {
                groups.extend(support.source_face_indices.iter().map(|&id| vec![id]));
            } else {
                groups.push(support.source_face_indices.clone());
            }
        }
        let mut membership = vec![usize::MAX; polygons.len()];
        let mut supports = Vec::with_capacity(groups.len());
        for (id, group) in groups.into_iter().enumerate() {
            for &polygon in &group {
                let slot = membership
                    .get_mut(polygon)
                    .ok_or(MeshRecoveryError::DegeneratePolygon(polygon))?;
                if *slot != usize::MAX {
                    return Err(MeshRecoveryError::DegeneratePolygon(polygon));
                }
                *slot = id;
            }
            supports.push(measure_planar_support(&polygons, &self.positions, group));
        }
        if membership.contains(&usize::MAX) {
            return Err(MeshRecoveryError::InvalidEvidenceOptions);
        }
        let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
        let (support_interfaces, support_interface_indices) = build_support_interfaces(
            &self.positions,
            &polygons,
            &incidence,
            &membership,
            supports.len(),
        );
        let mut refined = self.clone();
        refined.planar_supports = supports;
        refined.source_planar_support = membership;
        refined.support_interfaces = support_interfaces;
        refined.support_interface_indices = support_interface_indices;
        Ok(refined)
    }

    /// Planar-support region IDs on the two sides of a regular manifold edge.
    ///
    /// This does not certify a sharp CAD seam: regions may be chord facets of
    /// a curved surface. Returns None for open or ambiguous source incidence.
    #[must_use]
    pub fn edge_planar_support_pair(&self, edge_index: usize) -> Option<[usize; 2]> {
        let edge = self.edges.get(edge_index)?;
        if edge.kind != FeatureEdgeKind::Dihedral {
            return None;
        }
        let &[first, second] = edge.source_face_indices.as_slice() else {
            return None;
        };
        Some([
            *self.source_planar_support.get(first)?,
            *self.source_planar_support.get(second)?,
        ])
    }

    /// Measure a potential smooth normal-rotation continuation through a
    /// central support, using two distinct adjacent source interfaces.
    ///
    /// This is a lazy local relation rather than a quadratic matrix of
    /// candidate surface tracks. A match is not evidence that the original
    /// designer used an analytic curved surface: it can also describe a
    /// deliberately faceted design.
    #[must_use]
    pub fn support_normal_rotation(
        &self,
        center: usize,
        first: usize,
        second: usize,
    ) -> Option<SupportNormalRotationEvidence> {
        if first == second || self.planar_supports.get(center).is_none() {
            return None;
        }
        let [a, b] = [
            self.support_interfaces.get(first)?,
            self.support_interfaces.get(second)?,
        ];
        let neighbor_of = |interface: &SupportInterface| {
            if interface.supports[0] == center {
                Some(interface.supports[1])
            } else if interface.supports[1] == center {
                Some(interface.supports[0])
            } else {
                None
            }
        };
        let [neighbor_a, neighbor_b] = [neighbor_of(a)?, neighbor_of(b)?];
        if neighbor_a == neighbor_b {
            return None;
        }
        let axis = a.seed_direction;
        let second_axis = b.seed_direction;
        let interface_axis_residual_radians = axis
            .cross(second_axis)
            .magnitude()
            .atan2(axis.dot(second_axis).abs());
        let [left, middle, right] = [
            self.planar_supports.get(neighbor_a)?.normal,
            self.planar_supports.get(center)?.normal,
            self.planar_supports.get(neighbor_b)?.normal,
        ];
        let normal_angle = |first: Vector3, second: Vector3| {
            first.cross(second).magnitude().atan2(first.dot(second))
        };
        let turn_radians = [normal_angle(left, middle), normal_angle(middle, right)];
        let signed_first = axis.dot(left.cross(middle));
        let signed_second = axis.dot(middle.cross(right));
        let first_difference = middle - left;
        let second_difference = right - middle;
        let normal_cross = first_difference.cross(second_difference);
        let normal_strength = normal_cross.magnitude();
        let normal_difference_scale = first_difference.magnitude() * second_difference.magnitude();
        let normal_space_axis_strength = if normal_difference_scale > 0.0 {
            (normal_strength / normal_difference_scale).min(1.0)
        } else {
            0.0
        };
        // Two collinear normal-difference vectors cannot locate a unique
        // normal-space rotation axis. Do not invent one by selecting a
        // coordinate direction or using the interface tangent as a fallback.
        let normal_space_axis =
            (normal_space_axis_strength > 1.0e-12).then(|| canonical_unit(normal_cross));
        let normal_space_axis_to_boundaries_radians = normal_space_axis.map(|candidate| {
            [axis, second_axis]
                .iter()
                .map(|&direction| {
                    candidate
                        .cross(direction)
                        .magnitude()
                        .atan2(candidate.dot(direction).abs())
                })
                .fold(0.0_f64, f64::max)
        });
        let reference_origin = self.positions[*a.source_edges.first()?.first()?];
        let second_origin = self.positions[*b.source_edges.first()?.first()?];
        Some(SupportNormalRotationEvidence {
            center_support: center,
            neighbor_supports: [neighbor_a, neighbor_b],
            interface_indices: [first, second],
            reference_direction: axis,
            interface_axis_residual_radians,
            turn_radians,
            turn_magnitude_mismatch_radians: (turn_radians[0] - turn_radians[1]).abs(),
            coherent_rotation: signed_first * signed_second > 0.0,
            max_normal_axis_dot: [left, middle, right]
                .into_iter()
                .map(|normal| axis.dot(normal).abs())
                .fold(0.0_f64, f64::max),
            transverse_interface_separation: ((second_origin - reference_origin).cross(axis))
                .magnitude(),
            max_interface_seed_line_offset: a.max_seed_line_offset.max(b.max_seed_line_offset),
            max_interface_segment_spread_radians: a
                .max_seed_direction_spread_radians
                .max(b.max_seed_direction_spread_radians),
            normal_space_axis,
            normal_space_axis_strength,
            normal_space_axis_to_boundaries_radians,
        })
    }

    /// Compare the infinite support lines of any two recovered feature
    /// chains, computing geometric relations lazily rather than materializing
    /// a quadratic matrix of all chain pairs.
    #[must_use]
    pub fn chain_pair(&self, first: usize, second: usize) -> Option<ChainPairEvidence> {
        let a = self.chains.get(first)?;
        let b = self.chains.get(second)?;
        let [a0, a1] = [
            self.positions[*a.vertices.first()?],
            self.positions[*a.vertices.last()?],
        ];
        let [b0, b1] = [
            self.positions[*b.vertices.first()?],
            self.positions[*b.vertices.last()?],
        ];
        let u = (a1 - a0).normalize();
        let v = (b1 - b0).normalize();
        let cross = u.cross(v);
        let sine = cross.magnitude();
        let angle = sine.atan2(u.dot(v).abs());
        let delta = b0 - a0;
        let reference_distance = delta.cross(u).magnitude();
        let supporting_distance = if sine == 0.0 {
            reference_distance
        } else {
            delta.dot(cross).abs() / sine
        };
        let coord = |p: Point3| robust::Coord3D {
            x: p.x,
            y: p.y,
            z: p.z,
        };
        Some(ChainPairEvidence {
            angle_radians: angle,
            exactly_coplanar: robust::orient3d(coord(a0), coord(a1), coord(b0), coord(b1)) == 0.0,
            exactly_coincident_carriers: seams::on_line(a0, a1, b0) && seams::on_line(a0, a1, b1),
            supporting_line_distance: supporting_distance,
            reference_transverse_distance: reference_distance,
        })
    }
}

/// Feature selection is descriptive rather than constructive: even a polygon
/// soup with overlapping bodies keeps its source edge incidence as evidence.
struct EdgeSelection {
    edges: Vec<FeatureEdge>,
    suppressed_coplanar: usize,
    non_manifold: usize,
    winding_conflicts: usize,
}

fn select_edges(
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    options: EvidenceOptions,
) -> EdgeSelection {
    let mut keys = incidence.keys().copied().collect::<Vec<_>>();
    keys.sort_unstable();
    let mut edges = Vec::new();
    let mut suppressed_coplanar = 0;
    let mut non_manifold = 0;
    let mut winding_conflicts = 0;
    for key in keys {
        let uses = &incidence[&key];
        let (kind, angle) = match uses.as_slice() {
            [(_face, _)] => (Some(FeatureEdgeKind::Boundary), None),
            [(left, l_forward), (right, r_forward)] => {
                let a = polygons[*left].normal;
                let b = polygons[*right].normal;
                let angle = a.cross(b).magnitude().atan2(a.dot(b));
                if l_forward == r_forward {
                    winding_conflicts += 1;
                    (Some(FeatureEdgeKind::WindingConflict), Some(angle))
                } else if angle > options.min_dihedral_angle_radians {
                    (Some(FeatureEdgeKind::Dihedral), Some(angle))
                } else {
                    suppressed_coplanar += 1;
                    (None, None)
                }
            }
            [_, _, _, ..] => {
                non_manifold += 1;
                (Some(FeatureEdgeKind::NonManifoldJunction), None)
            }
            [] => continue,
        };
        if let Some(kind) = kind {
            edges.push(FeatureEdge {
                vertices: key,
                source_face_indices: uses.iter().map(|(face, _)| *face).collect(),
                source_forward_winding: uses.iter().map(|(_, forward)| *forward).collect(),
                kind,
                dihedral_angle_radians: angle,
            });
        }
    }
    EdgeSelection {
        edges,
        suppressed_coplanar,
        non_manifold,
        winding_conflicts,
    }
}

const fn opposite_vertex(edge: &FeatureEdge, vertex: usize) -> usize {
    if edge.vertices[0] == vertex {
        edge.vertices[1]
    } else {
        edge.vertices[0]
    }
}

const fn endpoint_slot(edge: &FeatureEdge, vertex: usize) -> usize {
    if edge.vertices[0] == vertex { 0 } else { 1 }
}

/// Pair an incoming edge with its only exact straight-through continuation,
/// even when other unrelated source edges meet at that vertex. Both incoming
/// and outgoing rays must uniquely choose *each other*; overlapping or branched
/// collinear tracks remain separate rather than hiding a real ambiguity.
fn straight_continuations(positions: &[Point3], edges: &[FeatureEdge]) -> Vec<[Option<usize>; 2]> {
    let mut incidence = vec![Vec::new(); positions.len()];
    for (edge_id, edge) in edges.iter().enumerate() {
        for &vertex in &edge.vertices {
            incidence[vertex].push(edge_id);
        }
    }
    let mut continuation = vec![[None; 2]; edges.len()];
    for (vertex, adjacent) in incidence.iter().enumerate() {
        let mut candidates = vec![Vec::new(); adjacent.len()];
        for (i, &first_id) in adjacent.iter().enumerate() {
            let first = &edges[first_id];
            let opposite_first = positions[opposite_vertex(first, vertex)];
            let toward_first = opposite_first - positions[vertex];
            for (j, &second_id) in adjacent.iter().enumerate().skip(i + 1) {
                let second = &edges[second_id];
                if first.kind != second.kind {
                    continue;
                }
                let opposite_second = positions[opposite_vertex(second, vertex)];
                if toward_first.dot(opposite_second - positions[vertex]) > 0.0 {
                    continue;
                }
                if seams::on_segment(opposite_first, opposite_second, positions[vertex]).is_some() {
                    candidates[i].push(j);
                    candidates[j].push(i);
                }
            }
        }
        for (i, options) in candidates.iter().enumerate() {
            let [j] = options.as_slice() else {
                continue;
            };
            if candidates[*j].as_slice() == [i] {
                let edge_id = adjacent[i];
                continuation[edge_id][endpoint_slot(&edges[edge_id], vertex)] = Some(adjacent[*j]);
            }
        }
    }
    continuation
}

fn trace_chains(positions: &[Point3], edges: &[FeatureEdge]) -> Vec<CollinearChain> {
    let continuation = straight_continuations(positions, edges);
    let mut visited = vec![false; edges.len()];
    let mut chains = Vec::new();
    for (first_id, first_edge) in edges.iter().enumerate() {
        if visited[first_id] {
            continue;
        }
        // Rewind along the unique exact line to its first unvisited anchor.
        // The bounded loop is defensive against malformed cyclic source
        // topology (a geometrically valid straight finite cycle cannot exist).
        let mut start_edge = first_id;
        let mut start_vertex = first_edge.vertices[0];
        for _ in 0..edges.len() {
            let Some(previous) =
                continuation[start_edge][endpoint_slot(&edges[start_edge], start_vertex)]
            else {
                break;
            };
            if previous == first_id || visited[previous] {
                break;
            }
            start_vertex = opposite_vertex(&edges[previous], start_vertex);
            start_edge = previous;
        }
        let mut vertices = vec![start_vertex];
        let mut ids = Vec::new();
        let mut current = start_vertex;
        let mut next_id = start_edge;
        loop {
            if visited[next_id] {
                break;
            }
            visited[next_id] = true;
            ids.push(next_id);
            current = opposite_vertex(&edges[next_id], current);
            vertices.push(current);
            let Some(next) = continuation[next_id][endpoint_slot(&edges[next_id], current)] else {
                break;
            };
            next_id = next;
        }
        // Do not rely on merely pairwise straightness: prove every
        // intermediate vertex lies on the final full endpoint chord.
        let first = positions[vertices[0]];
        let last = positions[vertices[vertices.len() - 1]];
        let exact = first != last
            && vertices[1..vertices.len() - 1]
                .iter()
                .all(|&v| seams::on_segment(first, last, positions[v]).is_some());
        if exact {
            chains.push(CollinearChain {
                vertices,
                edge_indices: ids,
                kind: first_edge.kind,
            });
        } else {
            for id in ids {
                chains.push(CollinearChain {
                    vertices: edges[id].vertices.to_vec(),
                    edge_indices: vec![id],
                    kind: edges[id].kind,
                });
            }
        }
    }
    chains
}

pub(super) fn canonical_unit(vector: Vector3) -> Vector3 {
    let axis = (0..3)
        .max_by(|&a, &b| vector[a].abs().total_cmp(&vector[b].abs()))
        .unwrap_or(0);
    let directed = if vector[axis].is_sign_negative() {
        -vector
    } else {
        vector
    };
    directed.normalize()
}

// Angular-neighbor grid bounds the number of bundle seeds inspected per
// chain, even when a complex curved mesh has many unrelated directions.
// Register both +/- unit direction to avoid hemisphere-boundary artifacts.
fn direction_cell(direction: Vector3, width: f64) -> [f64; 3] {
    // Minimum cell width keeps indices below 2^53: integer cell values and
    // their immediate neighbors are exactly representable, without casts.
    [direction.x, direction.y, direction.z].map(|component| ((component + 1.0) / width).floor())
}

fn parallel_bundles(
    positions: &[Point3],
    chains: &[CollinearChain],
    options: EvidenceOptions,
) -> Vec<ParallelBundle> {
    let mut candidates = chains
        .iter()
        .enumerate()
        .map(|(id, chain)| {
            let delta =
                positions[chain.vertices[chain.vertices.len() - 1]] - positions[chain.vertices[0]];
            (id, delta.magnitude2(), canonical_unit(delta))
        })
        .collect::<Vec<_>>();
    candidates.sort_unstable_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    let width = (2.0 * (options.max_parallel_angle_radians * 0.5).sin()).max(f64::EPSILON * 8.0);
    let mut bins = HashMap::<[u64; 3], Vec<usize>>::new();
    let mut bundles: Vec<ParallelBundle> = Vec::new();
    for (chain_index, _, direction) in candidates {
        let start = positions[chains[chain_index].vertices[0]];
        let end = positions[chains[chain_index].vertices[chains[chain_index].vertices.len() - 1]];
        let mut possible = HashSet::new();
        for sign in [direction, -direction] {
            let cell = direction_cell(sign, width);
            for x in [-1.0, 0.0, 1.0] {
                for y in [-1.0, 0.0, 1.0] {
                    for z in [-1.0, 0.0, 1.0] {
                        let neighbor = [
                            (cell[0] + x).to_bits(),
                            (cell[1] + y).to_bits(),
                            (cell[2] + z).to_bits(),
                        ];
                        if let Some(ids) = bins.get(&neighbor) {
                            possible.extend(ids.iter().copied());
                        }
                    }
                }
            }
        }
        let closest = possible
            .into_iter()
            .filter_map(|id| {
                let dot = bundles[id].direction.dot(direction);
                let angle = bundles[id]
                    .direction
                    .cross(direction)
                    .magnitude()
                    .atan2(dot.abs());
                (angle <= options.max_parallel_angle_radians).then_some((id, angle))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        let index = if let Some((id, _)) = closest {
            id
        } else {
            let id = bundles.len();
            bundles.push(ParallelBundle {
                direction,
                origin: start,
                placements: Vec::new(),
            });
            for sign in [direction, -direction] {
                bins.entry(direction_cell(sign, width).map(f64::to_bits))
                    .or_default()
                    .push(id);
            }
            id
        };
        let group = &mut bundles[index];
        let delta = start - group.origin;
        let offset = delta - group.direction * delta.dot(group.direction);
        let begin = delta.dot(group.direction);
        let finish = (end - group.origin).dot(group.direction);
        let angular = group
            .direction
            .cross(direction)
            .magnitude()
            .atan2(group.direction.dot(direction).abs());
        group.placements.push(ParallelPlacement {
            chain_index,
            transverse_offset: offset,
            transverse_distance: offset.magnitude(),
            longitudinal_interval: [begin.min(finish), begin.max(finish)],
            angular_residual_radians: angular,
        });
    }
    bundles
}

/// Build approximate planar supports without assuming valid global manifold
/// topology: only oppositely oriented regular two-face edges connect regions.
/// This leaves multi-body contact and non-manifold junctions as boundaries.
fn planar_supports(
    polygons: &[InputPolygon],
    positions: &[Point3],
    incidence: &EdgeIncidence,
    options: PlanarRecoveryOptions,
) -> (Vec<usize>, Vec<PlanarSupport>) {
    let mut neighbors = vec![Vec::new(); polygons.len()];
    for uses in incidence.values() {
        if let [(first, first_forward), (second, second_forward)] = uses.as_slice()
            && first_forward != second_forward
        {
            neighbors[*first].push(*second);
            neighbors[*second].push(*first);
        }
    }
    for adjacent in &mut neighbors {
        adjacent.sort_unstable();
    }
    let (membership, groups) = segment_patches(polygons, positions, &neighbors, options);
    let supports = groups
        .into_iter()
        .map(|group| measure_planar_support(polygons, positions, group))
        .collect();
    (membership, supports)
}

fn measure_planar_support(
    polygons: &[InputPolygon],
    positions: &[Point3],
    source_face_indices: Vec<usize>,
) -> PlanarSupport {
    let seed = &polygons[source_face_indices[0]];
    let mut area = 0.0;
    let mut max_vertex_plane_residual = 0.0_f64;
    let mut max_normal_angle_radians = 0.0_f64;
    for &face_id in &source_face_indices {
        let face = &polygons[face_id];
        let reference = positions[face.vertices[0]];
        // Signed fan contributions cancel for concave polygons;
        // summing triangle magnitudes would overestimate their area.
        let mut signed_twice_area = 0.0;
        for vertex in face.vertices.windows(2).skip(1) {
            let first = positions[vertex[0]] - reference;
            let second = positions[vertex[1]] - reference;
            signed_twice_area += first.cross(second).dot(face.normal);
        }
        area = signed_twice_area.abs().mul_add(0.5, area);
        max_vertex_plane_residual = face
            .vertices
            .iter()
            .map(|&id| signed_plane_distance(positions[id], seed.origin, seed.normal).abs())
            .fold(max_vertex_plane_residual, f64::max);
        let angular = face
            .normal
            .cross(seed.normal)
            .magnitude()
            .atan2(face.normal.dot(seed.normal));
        max_normal_angle_radians = max_normal_angle_radians.max(angular);
    }
    PlanarSupport {
        source_face_indices,
        normal: seed.normal,
        origin: seed.origin,
        area,
        max_vertex_plane_residual,
        max_normal_angle_radians,
    }
}

/// Extract maximal exact straight feature chains and candidate parallel
/// spatial arrangements without presupposing any analytic CAD primitive.
///
/// Source coordinates are exactly welded initially. Independently closed,
/// nonzero-volume components retain distinct topological IDs at coincident
/// contacts; exact T-junction repair operates only on unresolved components.
/// Coplanar triangulation edges are suppressed.
/// Other source edges retain polygon adjacency and measured dihedrals.
/// Chain collinearity is exact in the stored coordinates, while direction
/// groups are only bounded-angular hypotheses.
///
/// # Errors
/// Returns [`MeshRecoveryError`] for invalid thresholds or invalid source
/// geometry. Non-manifold and inconsistent-winding edge uses remain
/// observable as explicitly flagged features rather than invalidating all
/// local geometric evidence.
pub fn analyze_mesh_evidence(
    mesh: &PolygonMesh,
    options: EvidenceOptions,
) -> Result<MeshEvidence, MeshRecoveryError> {
    validate_evidence_options(options)?;
    let (mut positions, mut polygons) = positions_and_polygons(mesh, KERNEL_EPSILON)?;
    let mut separated_vertices = Vec::new();
    let closed = source_topology::separate_closed_components(
        &mut positions,
        &mut polygons,
        &mut separated_vertices,
    );
    // A component whose original edge uses already form a closed chain
    // has no T-junctions to repair. Do not split it at another body's
    // coincident geometry, or use its vertices to attach unresolved faces.
    if closed.contains(&false) {
        seams::split_source_t_junctions(
            &positions,
            polygons
                .iter_mut()
                .zip(&closed)
                .filter_map(|(p, &closed)| (!closed).then_some(p)),
        );
        // Exact subdivision can reveal additional independently closed chains.
        source_topology::separate_closed_components(
            &mut positions,
            &mut polygons,
            &mut separated_vertices,
        );
    }
    Ok(analyze_welded_evidence(
        positions,
        &polygons,
        options,
        separated_vertices,
    ))
}

fn validate_evidence_options(options: EvidenceOptions) -> Result<(), MeshRecoveryError> {
    if !options.min_dihedral_angle_radians.is_finite()
        || !(0.0..std::f64::consts::PI).contains(&options.min_dihedral_angle_radians)
        || !options.max_parallel_angle_radians.is_finite()
        || !(0.0..std::f64::consts::FRAC_PI_2).contains(&options.max_parallel_angle_radians)
        || !options.planar_support.distance_tolerance.is_finite()
        || options.planar_support.distance_tolerance < 0.0
        || !options.planar_support.angular_tolerance_radians.is_finite()
        || !(0.0..std::f64::consts::FRAC_PI_2)
            .contains(&options.planar_support.angular_tolerance_radians)
    {
        return Err(MeshRecoveryError::InvalidEvidenceOptions);
    }
    Ok(())
}

fn analyze_welded_evidence(
    positions: Vec<Point3>,
    polygons: &[InputPolygon],
    options: EvidenceOptions,
    separated_vertices: Vec<SeparatedSourceVertex>,
) -> MeshEvidence {
    let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
    let EdgeSelection {
        edges,
        suppressed_coplanar,
        non_manifold,
        winding_conflicts,
    } = select_edges(polygons, &incidence, options);
    let (source_planar_support, planar_supports) =
        planar_supports(polygons, &positions, &incidence, options.planar_support);
    let (support_interfaces, support_interface_indices) = build_support_interfaces(
        &positions,
        polygons,
        &incidence,
        &source_planar_support,
        planar_supports.len(),
    );
    let chains = trace_chains(&positions, &edges);
    let parallel_bundles = parallel_bundles(&positions, &chains, options);
    let mut source_polygon_vertex_offsets = Vec::with_capacity(polygons.len() + 1);
    let mut source_polygon_vertex_indices = Vec::new();
    for polygon in polygons {
        source_polygon_vertex_offsets.push(source_polygon_vertex_indices.len());
        source_polygon_vertex_indices.extend_from_slice(&polygon.vertices);
    }
    source_polygon_vertex_offsets.push(source_polygon_vertex_indices.len());
    MeshEvidence {
        positions,
        separated_vertices,
        edges,
        chains,
        parallel_bundles,
        planar_supports,
        source_planar_support,
        source_polygon_vertex_offsets,
        source_polygon_vertex_indices,
        support_interfaces,
        support_interface_indices,
        suppressed_interior_edge_count: suppressed_coplanar,
        non_manifold_edge_count: non_manifold,
        inconsistent_winding_edge_count: winding_conflicts,
    }
}

#[cfg(test)]
mod tests;
