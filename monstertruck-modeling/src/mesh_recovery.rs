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

/// Source evidence for one reconstructed planar face.
#[derive(Clone, Debug)]
pub struct PlanarPatch {
    /// Source polygon indices in mesh.faces().face_iter() order.
    pub source_face_indices: Vec<usize>,
    /// Winding-preserving loops of welded position indices, outer first,
    /// followed by holes. Collinear boundary vertices are retained.
    pub boundary_loops: Vec<Vec<usize>>,
    /// Oriented normal of the seed polygon.
    pub normal: Vector3,
    /// Maximum observed source-vertex distance to the seed plane.
    pub max_planar_residual: f64,
}

/// Reconstructed planar faces with common topological vertices and edges.
///
/// The shell may be open. Even a Closed shell condition alone does not
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

fn edge_key(a: usize, b: usize) -> (EdgeKey, bool) {
    if a < b {
        ([a, b], true)
    } else {
        ([b, a], false)
    }
}

fn signed_plane_distance(point: Point3, origin: Point3, normal: Vector3) -> f64 {
    (point - origin).dot(normal)
}

fn positions_and_polygons(
    mesh: &PolygonMesh,
    tolerance: f64,
) -> Result<(Vec<Point3>, Vec<InputPolygon>), MeshRecoveryError> {
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
        let mut seen = HashSet::new();
        if vertices.len() < 3 || !vertices.iter().all(|&id| seen.insert(id)) {
            return Err(MeshRecoveryError::DegeneratePolygon(source));
        }
        let origin = positions[vertices[0]];
        let mut area_vector = Vector3::new(0.0, 0.0, 0.0);
        for i in 0..vertices.len() {
            let p = positions[vertices[i]] - origin;
            let q = positions[vertices[(i + 1) % vertices.len()]] - origin;
            area_vector += p.cross(q);
            if positions[vertices[i]] == positions[vertices[(i + 1) % vertices.len()]] {
                return Err(MeshRecoveryError::DegeneratePolygon(source));
            }
        }
        let length_sq = area_vector.magnitude2();
        if !length_sq.is_finite() || length_sq <= 0.0 {
            return Err(MeshRecoveryError::DegeneratePolygon(source));
        }
        let normal = area_vector / length_sq.sqrt();
        if vertices
            .iter()
            .any(|&id| signed_plane_distance(positions[id], origin, normal).abs() > tolerance)
        {
            return Err(MeshRecoveryError::NonPlanarPolygon(source));
        }
        polygons.push(InputPolygon {
            vertices,
            normal,
            origin,
        });
    }
    Ok((positions, polygons))
}

fn edge_incidence(
    polygons: &[InputPolygon],
) -> Result<(EdgeIncidence, PolygonNeighbors), MeshRecoveryError> {
    let mut incidence = HashMap::<EdgeKey, Vec<EdgeUse>>::default();
    for (polygon_id, polygon) in polygons.iter().enumerate() {
        for i in 0..polygon.vertices.len() {
            let a = polygon.vertices[i];
            let b = polygon.vertices[(i + 1) % polygon.vertices.len()];
            let (key, forward) = edge_key(a, b);
            incidence
                .entry(key)
                .or_default()
                .push((polygon_id, forward));
        }
    }
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
    let mut loops = Vec::new();
    while !next.is_empty() {
        let Some(&start) = next.keys().min() else {
            return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
        };
        let mut cycle = Vec::new();
        let mut current = start;
        loop {
            cycle.push(current);
            let Some(other) = next.remove(&current) else {
                return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
            };
            if other == start {
                break;
            }
            if cycle.len() > incoming.len() || !incoming.contains_key(&other) {
                return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
            }
            current = other;
        }
        if cycle.len() < 3 {
            return Err(MeshRecoveryError::InvalidPatchBoundary(patch_id));
        }
        loops.push(cycle);
    }
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
/// patches would create T-junctions. The caller can check shell_condition()
/// and validate solid geometry separately.
///
/// The distance_tolerance is a maximum fitting residual, not a request to
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
    let (positions, polygons) = positions_and_polygons(mesh, options.distance_tolerance)?;
    let (incidence, neighbors) = edge_incidence(&polygons)?;
    let (membership, regions) = segment_patches(&polygons, &positions, &neighbors, options);
    let vertices = Vertex::from_points(&positions);
    let mut edges = HashMap::default();
    let mut faces = Vec::with_capacity(regions.len());
    let mut patches = Vec::with_capacity(regions.len());

    for (patch_id, member_faces) in regions.into_iter().enumerate() {
        let seed = &polygons[member_faces[0]];
        let loops = boundary_loops(
            patch_id,
            &member_faces,
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
        faces.push(build_face(&loops, seed, &positions, &vertices, &mut edges)?);
        patches.push(PlanarPatch {
            source_face_indices: member_faces,
            boundary_loops: loops,
            normal: seed.normal,
            max_planar_residual,
        });
    }
    Ok(RecoveredPlanarMesh {
        canonical_positions: positions,
        shell: faces.into_iter().collect(),
        patches,
    })
}

#[cfg(test)]
mod tests;
