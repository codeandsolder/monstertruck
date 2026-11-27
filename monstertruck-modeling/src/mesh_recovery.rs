//! Recovery of analytic planar B-rep faces from polygon meshes.
//!
//! This is inverse meshing, not tessellation: connected, consistently oriented
//! coplanar polygons become one planar face with their internal edges removed.
//! Source face membership and unsimplified boundary edges are retained, so
//! subsequent recovery stages can share and certify seams with curved patches.

use crate::{
    Curve, Face, InnerSpace, KERNEL_EPSILON, Plane, Point3, Shell, Surface, Vector3, Vertex, Wire,
    builder,
};
use monstertruck_mesh::PolygonMesh;
use monstertruck_topology::{Edge, errors::Error as TopologyError};
use rustc_hash::FxHashMap as HashMap;
use std::collections::{HashSet, VecDeque};

/// Geometric limits used to group planar mesh facets.
///
/// Both limits are absolute. The distance is expressed in the input mesh's
/// coordinate units; the angle is in radians. No coordinate snapping or
/// tolerance-based vertex welding is performed.
#[derive(Clone, Copy, Debug)]
pub struct PlanarRecoveryOptions {
    /// Maximum point-to-seed-plane distance for every vertex in a patch.
    pub distance_tolerance: f64,
    /// Maximum angular deviation of each polygon normal from the patch seed.
    pub angular_tolerance_radians: f64,
}

impl Default for PlanarRecoveryOptions {
    fn default() -> Self {
        Self {
            distance_tolerance: KERNEL_EPSILON,
            angular_tolerance_radians: 1.0e-6,
        }
    }
}

mod analytic;
pub use analytic::{
    AnalyticDiscoveryStats, AnalyticRecoveryError, AnalyticRecoveryOptions,
    AnalyticRegionDiscovery, AxisDiscoverySummary, CircleAxisSummary, CoplanarSupportReassignment,
    JointFitAttempt, JointFitMethod, JointFitRejection, JointFitReport, RecoveredAnalyticMesh,
    SelectedSurfaceProposal, SupportRefinementPass, construct_analytic_recovery,
    discover_analytic_regions, recover_analytic_evidence, recover_analytic_mesh,
};
mod boundary_facet;
pub use boundary_facet::{
    BoundaryFacetApproximation, BoundaryFacetError, BoundaryFacetInput, approximate_boundary_facet,
};
mod boundary_patch;
pub use boundary_patch::{
    BoundaryPatchApproximation, BoundaryPatchError, BoundaryPatchOptions, PrescribedSourceEdge,
    approximate_boundary_patch,
};
mod cylinder;
mod through_hole;
pub use cylinder::{CylinderRecoveryOptions, RecoveredCylinder, recover_cylindrical_solid};
pub use through_hole::{RecoveredCylindricalThroughHole, recover_cylindrical_through_hole};
mod evidence;
mod planar_domain;
pub use planar_domain::{PlanarDomainCertificate, certify_projected_planar_domain};
mod planar_fans;
pub use planar_fans::{PlanarFanOptions, PlanarFanSimplification, simplify_planar_fan_shell};
mod sample_fit;
mod seams;
pub use evidence::{
    AxialRadialStation, CanceledSourceSheet, CarrierCoincidenceEvidence, CarrierResidual,
    ChainPairEvidence, CircleAxisDiscovery, CircleAxisEvidence, CircleAxisOptions,
    CircleAxisStation, CircularMeridianEvidence, CollapsedSourceFacet, CollinearChain,
    ContactCircleEvidence, ContactConicEvidence, ContactCurveEvidence, ContactLineEvidence,
    ContactLinearSectionEvidence, ContactToroidalSectionEvidence, CurvatureGrowthHypothesis,
    CurvatureGrowthIndex, CurvatureGrowthOptions, CurvatureGrowthStep, CurvatureTrack,
    CurvatureTrackOptions, EvidenceOptions, FeatureEdge, FeatureEdgeKind, HybridConnectorOptions,
    HybridConnectorRecovery, HybridConnectorRelease, HybridPatchEvidence,
    HybridPlanarSimplification, HybridPolygonFaces, HybridRecovery, HybridRecoveryError,
    HybridRefinementOptions, HybridRefinementSeed, HybridSupportSearchOptions,
    HybridSupportSearchRecovery, LiftedRegionPatch, LiftedSourceChord, MeridianGraphOptions,
    MeridianStation, MeridianStationGraph, MeridianTransition, MeshEvidence, NativeChartObstacle,
    NativeRegionFace, NativeTrimError, NativeTrimGraph, NativeTrimOptions, OrientedTrim,
    ParallelBundle, ParallelPlacement, PatchSupportConflict, PatchSupportLocus, PlanarSupport,
    PlanarTrimCrossing, ReconciledPartition, ReconciledPlane, ReconciliationOptions,
    RefinedHybridRecovery, RegionBoundary, RegionBoundaryEdge, RegionBoundaryError,
    RegionBoundaryGap, RegionBoundaryIndex, RegionContactGraph, RegionContactPath,
    RegionPatchError, RegionPatchOptions, RegionTrimBoundary, RegionTrimLoops,
    SeparatedSourceVertex, SharedAnalyticTrim, SharedTrimGeometry, SharedTrimGraph,
    SharedTrimOptions, SharedTrimPiece, SharedTrimSegmentation, SourceCircleObservation,
    SourceForcedPair, SourceMultiplicityMethod, SourceMultiplicityPair,
    SourceMultiplicityResolution, SourceProvenance, SourceSectorPair, SourceTransverseCrossing,
    SpatialTrackEvidence, SpatialTrackJoinEvidence, StationLinearCandidate,
    StationMeridianCandidate, StationMeridianDiscovery, SupportBoundaryTopology, SupportInterface,
    SupportNormalRotationEvidence, SurfaceCarrier, SurfaceCurveAgreement, SurfaceRegionBudget,
    SurfaceRegionPartition, SurfaceRegionProposal, analyze_mesh_evidence, discover_circle_axes,
    measure_circular_meridian, measure_curvature_spatial, measure_curvature_spatial_join,
    measure_spatial_support_set, trace_curvature_tracks,
};

/// Source evidence for one reconstructed planar face.
#[derive(Clone, Debug)]
pub struct PlanarPatch {
    /// Source polygon indices in `mesh.face_iter()` order.
    pub source_face_indices: Vec<usize>,
    /// Winding-preserving loops of welded position indices, outer first,
    /// followed by holes. Collinear seam segments are simplified globally.
    pub boundary_loops: Vec<Vec<usize>>,
    /// Oriented normal of the seed polygon.
    pub normal: Vector3,
    /// Maximum observed source-vertex distance to the seed plane.
    pub max_planar_residual: f64,
}

/// Reconstructed planar faces with common topological vertices and edges.
///
/// The shell may be open. Even a `Closed` shell condition alone does not
/// certify positive material volume or absence of self-intersections.
#[derive(Debug)]
pub struct RecoveredPlanarMesh {
    /// Exact-coordinate-welded positions referenced by patch boundary loops.
    pub canonical_positions: Vec<Point3>,
    /// One native analytic plane face per planar patch.
    pub shell: Shell,
    /// Patches in the same order as faces in shell.
    pub patches: Vec<PlanarPatch>,
}

/// Errors that prevent trustworthy reconstruction.
#[derive(Debug, thiserror::Error)]
pub enum MeshRecoveryError {
    /// An option was NaN, infinite, negative, or outside its legal range.
    #[error("recovery tolerances must be finite and nonnegative, with angle < pi/2")]
    InvalidTolerance,
    /// A source point is not finite.
    #[error("mesh position {0} is not finite")]
    NonFinitePosition(usize),
    /// A polygon contains a position index outside the input buffer.
    #[error("mesh polygon {polygon} references missing position {position}")]
    InvalidPositionIndex { polygon: usize, position: usize },
    /// A polygon has a zero-length edge, duplicated vertex, or zero area.
    #[error("mesh polygon {0} has degenerate geometry")]
    DegeneratePolygon(usize),
    /// A polygon is warped, so its single polygonal boundary is not planar.
    #[error("mesh polygon {0} is not planar within the requested tolerance")]
    NonPlanarPolygon(usize),
    /// More than two polygons use the same undirected edge.
    #[error("mesh edge {0:?} has more than two incident polygons")]
    NonManifoldEdge([usize; 2]),
    /// Both incident polygons traverse their common edge in the same direction.
    #[error("mesh edge {0:?} has inconsistent incident winding")]
    InconsistentWinding([usize; 2]),
    /// Patch boundary is pinched, branched, open, or has ambiguous orientation.
    #[error("planar patch {0} has ambiguous boundary loops")]
    InvalidPatchBoundary(usize),
    /// The native cylinder builder could not produce valid shared trims.
    #[error("native analytic cylinder construction failed")]
    CylinderConstruction,
    /// Unsupported or non-finite cylinder recovery tolerance.
    #[error("cylinder tolerances must be finite and nonnegative, with at least six ring samples")]
    InvalidCylinderOptions,
    /// Feature detection or direction clustering thresholds are invalid.
    #[error(
        "evidence thresholds must be finite, with dihedral angle in [0,pi) and parallel angle in [0,pi/2)"
    )]
    InvalidEvidenceOptions,
    /// A patch would require nontrivial projection to become a valid B-rep face.
    #[error(
        "planar patch {patch} has source-to-plane residual {residual}, exceeding the kernel tolerance"
    )]
    ExcessivePlanarResidual { patch: usize, residual: f64 },
    /// Native topological validation failed.
    #[error(transparent)]
    Topology(#[from] TopologyError),
}

#[derive(Clone)]
struct InputPolygon {
    vertices: Vec<usize>,
    normal: Vector3,
    origin: Point3,
}

type EdgeKey = [usize; 2];
type EdgeUse = (usize, bool);
type EdgeIncidence = HashMap<EdgeKey, Vec<EdgeUse>>;
type PolygonNeighbors = Vec<Vec<usize>>;

const fn edge_key(a: usize, b: usize) -> (EdgeKey, bool) {
    if a < b {
        ([a, b], true)
    } else {
        ([b, a], false)
    }
}

fn signed_plane_distance(point: Point3, origin: Point3, normal: Vector3) -> f64 {
    (point - origin).dot(normal)
}

/// Integrate signed polyhedral volume by a fan triangulation of each source
/// polygon. All coordinates are relative to one source vertex, so translating
/// a small mechanical part far from the world origin does not unnecessarily
/// destroy precision. This certifies global winding *sign* only for an already
/// validated closed mesh; it does not certify non-self-intersection.
///
/// Compensated summation protects thin solids with cancelling oriented faces.
fn signed_source_volume(positions: &[Point3], polygons: &[InputPolygon]) -> Option<f64> {
    signed_polygon_volume(
        positions,
        polygons.iter().map(|polygon| polygon.vertices.as_slice()),
    )
}

fn signed_polygon_volume<'a>(
    positions: &[Point3],
    polygons: impl IntoIterator<Item = &'a [usize]>,
) -> Option<f64> {
    let mut polygons = polygons.into_iter();
    let first = polygons.next()?;
    let origin = *positions.get(*first.first()?)?;
    let mut six_volume = 0.0_f64;
    let mut compensation = 0.0_f64;
    for polygon in std::iter::once(first).chain(polygons) {
        if polygon.len() < 3 {
            return None;
        }
        let a = *positions.get(polygon[0])? - origin;
        for pair in polygon[1..].windows(2) {
            let b = *positions.get(pair[0])? - origin;
            let c = *positions.get(pair[1])? - origin;
            let term = a.dot(b.cross(c));
            let corrected = term - compensation;
            let next = six_volume + corrected;
            compensation = (next - six_volume) - corrected;
            six_volume = next;
        }
    }
    let volume = six_volume / 6.0;
    volume.is_finite().then_some(volume)
}

/// Positive material-volume orientation for a previously proven closed source
/// mesh. Does not imply a global non-self-intersection certificate.
fn positive_source_volume(positions: &[Point3], polygons: &[InputPolygon]) -> Option<f64> {
    signed_source_volume(positions, polygons).filter(|volume| *volume > 0.0)
}

fn positions_and_polygons(
    mesh: &PolygonMesh,
    tolerance: f64,
) -> Result<(Vec<Point3>, Vec<InputPolygon>), MeshRecoveryError> {
    let (positions, source) = welded_source_polygons(mesh)?;
    let polygons = source
        .into_iter()
        .enumerate()
        .map(|(id, vertices)| input_polygon_from_ids(&positions, vertices, id, tolerance))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((positions, polygons))
}

fn welded_source_polygons(
    mesh: &PolygonMesh,
) -> Result<(Vec<Point3>, Vec<Vec<usize>>), MeshRecoveryError> {
    // Only exactly identical coordinates are welded. Distinct nearby points
    // can be separate intentional material boundaries; fuzzy welding is unsafe.
    let mut positions = Vec::new();
    let mut index = HashMap::<[u64; 3], usize>::default();
    let mut remap = Vec::with_capacity(mesh.positions().len());
    for (id, &point) in mesh.positions().iter().enumerate() {
        if ![point.x, point.y, point.z].into_iter().all(f64::is_finite) {
            return Err(MeshRecoveryError::NonFinitePosition(id));
        }
        // Canonicalize signed zero, which is the same geometric position.
        let key = [point.x, point.y, point.z].map(|value| {
            if value == 0.0 {
                0.0_f64.to_bits()
            } else {
                value.to_bits()
            }
        });
        let canonical = *index.entry(key).or_insert_with(|| {
            let next = positions.len();
            positions.push(point);
            next
        });
        remap.push(canonical);
    }

    let mut polygons = Vec::new();
    for (source, face) in mesh.face_iter().enumerate() {
        let vertices = face
            .iter()
            .map(|v| {
                remap
                    .get(v.pos)
                    .copied()
                    .ok_or(MeshRecoveryError::InvalidPositionIndex {
                        polygon: source,
                        position: v.pos,
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        polygons.push(vertices);
    }
    Ok((positions, polygons))
}

/// Measure an already welded original polygon; shared by initial ingestion
/// and subsequent support refinement without changing source coordinates.
fn input_polygon_from_ids(
    positions: &[Point3],
    vertices: Vec<usize>,
    source: usize,
    tolerance: f64,
) -> Result<InputPolygon, MeshRecoveryError> {
    let mut seen = HashSet::new();
    if vertices.len() < 3 || !vertices.iter().all(|&id| seen.insert(id)) {
        return Err(MeshRecoveryError::DegeneratePolygon(source));
    }
    let origin = positions[vertices[0]];
    let mut area_vector = Vector3::new(0.0, 0.0, 0.0);
    let mut compensation = Vector3::new(0.0, 0.0, 0.0);
    for pair in vertices[1..].windows(2) {
        let term = seams::oriented_triangle_area(origin, positions[pair[0]], positions[pair[1]]);
        let corrected = term - compensation;
        let next = area_vector + corrected;
        compensation = (next - area_vector) - corrected;
        area_vector = next;
    }
    if vertices
        .iter()
        .zip(vertices.iter().cycle().skip(1))
        .any(|(&a, &b)| positions[a] == positions[b])
    {
        return Err(MeshRecoveryError::DegeneratePolygon(source));
    }
    // Scale before normalization so a finite nonzero area is not lost by
    // squaring subnormal components or overflowing a large area vector.
    let scale = area_vector
        .x
        .abs()
        .max(area_vector.y.abs())
        .max(area_vector.z.abs());
    if ![area_vector.x, area_vector.y, area_vector.z]
        .into_iter()
        .all(f64::is_finite)
        || scale <= 0.0
    {
        return Err(MeshRecoveryError::DegeneratePolygon(source));
    }
    let normal = (area_vector / scale).normalize();
    if vertices
        .iter()
        .any(|&id| signed_plane_distance(positions[id], origin, normal).abs() > tolerance)
    {
        return Err(MeshRecoveryError::NonPlanarPolygon(source));
    }
    Ok(InputPolygon {
        vertices,
        normal,
        origin,
    })
}

/// Collect every original oriented polygon use without enforcing manifoldness.
/// The evidence frontend must preserve ambiguous connections and report them,
/// while B-rep construction calls the strict validator below.
fn source_edge_incidence<'a>(polygons: impl IntoIterator<Item = &'a [usize]>) -> EdgeIncidence {
    let mut incidence = HashMap::<EdgeKey, Vec<EdgeUse>>::default();
    for (polygon_id, vertices) in polygons.into_iter().enumerate() {
        for i in 0..vertices.len() {
            let a = vertices[i];
            let b = vertices[(i + 1) % vertices.len()];
            let (key, forward) = edge_key(a, b);
            incidence
                .entry(key)
                .or_default()
                .push((polygon_id, forward));
        }
    }
    incidence
}

fn edge_incidence(
    polygons: &[InputPolygon],
) -> Result<(EdgeIncidence, PolygonNeighbors), MeshRecoveryError> {
    let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
    let mut neighbors = vec![Vec::new(); polygons.len()];
    for (&key, uses) in &incidence {
        match uses.as_slice() {
            [a, b] => {
                if a.1 == b.1 {
                    return Err(MeshRecoveryError::InconsistentWinding(key));
                }
                neighbors[a.0].push(b.0);
                neighbors[b.0].push(a.0);
            }
            [_, _, ..] => return Err(MeshRecoveryError::NonManifoldEdge(key)),
            _ => {}
        }
    }
    Ok((incidence, neighbors))
}

fn segment_patches(
    polygons: &[InputPolygon],
    positions: &[Point3],
    neighbors: &[Vec<usize>],
    options: PlanarRecoveryOptions,
) -> (Vec<usize>, Vec<Vec<usize>>) {
    let mut membership = vec![usize::MAX; polygons.len()];
    let mut patches = Vec::new();
    let cos_angle = options.angular_tolerance_radians.cos();
    for seed_id in 0..polygons.len() {
        if membership[seed_id] != usize::MAX {
            continue;
        }
        let patch_id = patches.len();
        membership[seed_id] = patch_id;
        let seed = &polygons[seed_id];
        let mut queue = VecDeque::from([seed_id]);
        let mut member_faces = Vec::new();
        while let Some(current) = queue.pop_front() {
            member_faces.push(current);
            for &adjacent in &neighbors[current] {
                if membership[adjacent] != usize::MAX {
                    continue;
                }
                let polygon = &polygons[adjacent];
                if polygon.normal.dot(seed.normal) < cos_angle {
                    continue;
                }
                if polygon.vertices.iter().any(|&v| {
                    signed_plane_distance(positions[v], seed.origin, seed.normal).abs()
                        > options.distance_tolerance
                }) {
                    continue;
                }
                membership[adjacent] = patch_id;
                queue.push_back(adjacent);
            }
        }
        member_faces.sort_unstable();
        patches.push(member_faces);
    }
    (membership, patches)
}

/// Decompose a one-successor source edge graph into deterministic oriented
/// cycles. No planar geometry or outer/inner classification is assumed.
fn trace_oriented_boundary_cycles(mut next: HashMap<usize, usize>) -> Option<Vec<Vec<usize>>> {
    let mut starts = next.keys().copied().collect::<Vec<_>>();
    starts.sort_unstable();
    let mut loops = Vec::new();
    for start in starts {
        if !next.contains_key(&start) {
            continue;
        }
        let mut cycle = Vec::new();
        let mut current = start;
        loop {
            cycle.push(current);
            current = next.remove(&current)?;
            if current == start {
                break;
            }
        }
        if cycle.len() < 3 {
            return None;
        }
        loops.push(cycle);
    }
    Some(loops)
}

fn boundary_loops(
    patch_id: usize,
    member_faces: &[usize],
    membership: &[usize],
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    positions: &[Point3],
    normal: Vector3,
) -> Result<Vec<Vec<usize>>, MeshRecoveryError> {
    let mut next = HashMap::<usize, usize>::default();
    let mut incoming = HashMap::<usize, usize>::default();
    for &face_id in member_faces {
        let vertices = &polygons[face_id].vertices;
        for i in 0..vertices.len() {
            let a = vertices[i];
            let b = vertices[(i + 1) % vertices.len()];
            let (key, _) = edge_key(a, b);
            let uses = &incidence[&key];
            if uses.len() == 2 && membership[uses[0].0] == membership[uses[1].0] {
                continue;
            }
            if next.insert(a, b).is_some() || incoming.insert(b, a).is_some() {
                return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
            }
        }
    }
    if next.is_empty() || next.len() != incoming.len() {
        return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
    }
    let mut loops = trace_oriented_boundary_cycles(next)
        .ok_or(MeshRecoveryError::InvalidPatchBoundary(patch_id))?;
    // The winding of an orientable planar patch gives one positive outer loop
    // and zero or more negative inner loops.
    let signed_area = |cycle: &[usize]| -> f64 {
        let origin = positions[cycle[0]];
        cycle
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                let a = positions[v] - origin;
                let b = positions[cycle[(i + 1) % cycle.len()]] - origin;
                a.cross(b).dot(normal)
            })
            .sum::<f64>()
            / 2.0
    };
    let mut outer = None;
    for (i, cycle) in loops.iter().enumerate() {
        let area = signed_area(cycle);
        if !area.is_finite() || area == 0.0 {
            return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
        }
        if area > 0.0 && outer.replace(i).is_some() {
            return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
        }
    }
    let Some(outer_index) = outer else {
        return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
    };
    loops.swap(0, outer_index);
    Ok(loops)
}

fn build_face(
    loops: &[Vec<usize>],
    seed: &InputPolygon,
    positions: &[Point3],
    vertices: &[Vertex],
    edges: &mut HashMap<EdgeKey, Edge<Point3, Curve>>,
) -> Result<Face, MeshRecoveryError> {
    let mut wires = Vec::with_capacity(loops.len());
    for cycle in loops {
        let mut contour = Vec::with_capacity(cycle.len());
        for i in 0..cycle.len() {
            let a = cycle[i];
            let b = cycle[(i + 1) % cycle.len()];
            let (key, forward) = edge_key(a, b);
            let edge = edges
                .entry(key)
                .or_insert_with(|| builder::line(&vertices[key[0]], &vertices[key[1]]));
            contour.push(if forward {
                edge.clone()
            } else {
                edge.inverse()
            });
        }
        wires.push(Wire::from(contour));
    }
    // Use the source plane rather than inferring it from an often very thin
    // polygon boundary; vertices are certified against this plane above.
    let p0 = seed.origin;
    let p1 = positions[seed.vertices[1]];
    let tangent = p1 - p0;
    let p2 = p0 + seed.normal.cross(tangent);
    let plane = Plane::new(p0, p1, p2);
    Ok(Face::new(wires, Surface::Plane(plane))?)
}

/// Reconstruct planar mesh regions as native B-rep faces.
///
/// Edges are shared between recovered faces when the source polygons share the
/// same exact geometric edge. The function validates source orientation,
/// boundary manifoldness, and planar residuals, but deliberately does not
/// merge short collinear boundary edges: doing so independently on neighboring
/// patches would create T-junctions. Existing geometric T-junctions are split
/// before incidence is built, then globally collinear seam chains are simplified
/// on both sides at once. The caller can check `shell_condition()` and validate
/// solid geometry separately.
///
/// The `distance_tolerance` is a maximum fitting residual, not a request to
/// weld nearby vertices. Resulting analytic faces can approximate the input
/// within that residual; exact reproduction is only guaranteed for coplanar
/// source polygons.
///
/// # Errors
/// Rejects nonfinite coordinates/tolerances, degenerate or nonplanar source
/// polygons, inconsistent winding, non-manifold edges, and invalid patch loops.
pub fn recover_planar_mesh(
    mesh: &PolygonMesh,
    options: PlanarRecoveryOptions,
) -> Result<RecoveredPlanarMesh, MeshRecoveryError> {
    if !options.distance_tolerance.is_finite()
        || options.distance_tolerance < 0.0
        || !options.angular_tolerance_radians.is_finite()
        || !(0.0..std::f64::consts::FRAC_PI_2).contains(&options.angular_tolerance_radians)
    {
        return Err(MeshRecoveryError::InvalidTolerance);
    }
    let (positions, mut polygons) = positions_and_polygons(mesh, options.distance_tolerance)?;
    seams::split_source_t_junctions(&positions, &mut polygons);
    let (incidence, neighbors) = edge_incidence(&polygons)?;
    let (membership, regions) = segment_patches(&polygons, &positions, &neighbors, options);
    let mut patches = Vec::with_capacity(regions.len());

    for (patch_id, member_faces) in regions.iter().enumerate() {
        let seed = &polygons[member_faces[0]];
        let loops = boundary_loops(
            patch_id,
            member_faces,
            &membership,
            &polygons,
            &incidence,
            &positions,
            seed.normal,
        )?;
        let max_planar_residual = member_faces
            .iter()
            .flat_map(|&face| polygons[face].vertices.iter().copied())
            .map(|v| signed_plane_distance(positions[v], seed.origin, seed.normal).abs())
            .fold(0.0_f64, f64::max);
        if max_planar_residual > KERNEL_EPSILON {
            return Err(MeshRecoveryError::ExcessivePlanarResidual {
                patch: patch_id,
                residual: max_planar_residual,
            });
        }
        patches.push(PlanarPatch {
            source_face_indices: member_faces.clone(),
            boundary_loops: loops,
            normal: seed.normal,
            max_planar_residual,
        });
    }
    seams::simplify_shared_seams(&mut patches, &positions, &incidence, &membership)?;
    let vertices = Vertex::from_points(&positions);
    let mut edges = HashMap::default();
    let mut faces = Vec::with_capacity(patches.len());
    for (patch, member_faces) in patches.iter().zip(&regions) {
        let seed = &polygons[member_faces[0]];
        faces.push(build_face(
            &patch.boundary_loops,
            seed,
            &positions,
            &vertices,
            &mut edges,
        )?);
    }
    Ok(RecoveredPlanarMesh {
        canonical_positions: positions,
        shell: faces.into_iter().collect(),
        patches,
    })
}

#[cfg(test)]
mod tests;
