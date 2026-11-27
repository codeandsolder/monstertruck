//! Strict assembly using the shared native edge sample identities.
//!
//! Face-local UVs and normals remain distinct. Only known native vertex
//! and edge sample identities can share a global position. Coordinate
//! proximity is used to associate a face boundary with its own source
//! samples, never to merge unrelated source identities.

use super::{
    CompressedEdge, CompressedFace, CompressedShell, Faces, Point3, PolygonMesh, PolylineCurve,
    StandardAttributes, StandardVertex,
};
use monstertruck_mesh::MetricSpace;
use std::collections::{BTreeMap, BTreeSet};

/// Failure to match the emitted face mesh to its native source boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeshStitchError {
    /// Invalid budget, source vertex or sampled source edge.
    InvalidSource,
    /// One native face is absent or has no usable polygons.
    MissingFace(usize),
    /// A face boundary position has no unique matching native sample.
    UnmatchedBoundary {
        /// Native face index.
        face: usize,
        /// Face-local position index.
        position: usize,
    },
    /// Emitted face incidence differs from the native sample chain.
    BoundaryMismatch(usize),
    /// The final mesh has inconsistent edge incidence or a disconnected vertex link.
    InvalidShell,
}

impl std::fmt::Display for MeshStitchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::InvalidSource => write!(f, "invalid native sample source or stitching budget"),
            Self::MissingFace(face) => write!(f, "native face {face} has no usable mesh"),
            Self::UnmatchedBoundary { face, position } => write!(
                f,
                "native face {face} position {position} has no unique source sample"
            ),
            Self::BoundaryMismatch(face) => write!(
                f,
                "native face {face} does not consume its source edge sample chain"
            ),
            Self::InvalidShell => write!(
                f,
                "stitched mesh has inconsistent edge incidence or disconnected vertex link"
            ),
        }
    }
}
impl std::error::Error for MeshStitchError {}

type Incidence = BTreeMap<[usize; 2], [usize; 2]>;

fn add_use(uses: &mut Incidence, a: usize, b: usize) {
    if a == b {
        return;
    }
    let (key, dir) = if a < b { ([a, b], 0) } else { ([b, a], 1) };
    uses.entry(key).or_default()[dir] += 1;
}

fn add_polygon(uses: &mut Incidence, polygon: &[StandardVertex]) {
    for (a, b) in polygon.iter().zip(polygon.iter().cycle().skip(1)) {
        add_use(uses, a.pos, b.pos);
    }
}

// Resolve every local boundary vertex, plus explicitly represented
// native seam segments that already have two face-local uses. A mere
// coincident interior position is not evidence for a shared identity.
fn face_position_map(
    face: usize,
    mesh: &PolygonMesh,
    samples: &FaceSamples,
    positions: &[Point3],
    tolerance: f64,
) -> Result<BTreeMap<usize, usize>, MeshStitchError> {
    let mut incidence = Incidence::new();
    for polygon in mesh.face_iter() {
        add_polygon(&mut incidence, polygon);
    }
    let mut result = BTreeMap::new();
    let mut candidates = BTreeMap::new();
    for (&[a, b], counts) in &incidence {
        let mut resolved = [None; 2];
        for (i, id) in [a, b].into_iter().enumerate() {
            let point = *mesh
                .positions()
                .get(id)
                .filter(|&&p| finite(p))
                .ok_or(MeshStitchError::InvalidSource)?;
            let candidate = *candidates
                .entry(id)
                .or_insert_with(|| samples.resolve(point, positions, tolerance));
            if counts[0] != counts[1] {
                result.insert(
                    id,
                    candidate.ok_or(MeshStitchError::UnmatchedBoundary { face, position: id })?,
                );
            }
            resolved[i] = candidate;
        }
        if let [Some(x), Some(y)] = resolved
            && samples.expected.contains_key(&[x.min(y), x.max(y)])
        {
            result.insert(a, x);
            result.insert(b, y);
        }
    }
    Ok(result)
}

// The link of a manifold mesh vertex must be one connected cycle or
// path. Edge incidence alone accepts bow-tie vertices.
fn manifold_links(faces: &Faces) -> bool {
    let mut links = BTreeMap::<usize, Vec<[usize; 2]>>::new();
    for polygon in faces.face_iter() {
        for (i, v) in polygon.iter().enumerate() {
            let before = polygon[(i + polygon.len() - 1) % polygon.len()].pos;
            let after = polygon[(i + 1) % polygon.len()].pos;
            links.entry(v.pos).or_default().push([before, after]);
        }
    }
    links.values().all(|pairs| {
        let mut neighbors = BTreeMap::<usize, Vec<usize>>::new();
        for &[a, b] in pairs {
            neighbors.entry(a).or_default().push(b);
            neighbors.entry(b).or_default().push(a);
        }
        if neighbors.values().any(|row| row.len() > 2) {
            return false;
        }
        let Some((&first, _)) = neighbors.first_key_value() else {
            return false;
        };
        let mut visited = BTreeSet::from([first]);
        let mut pending = vec![first];
        while let Some(id) = pending.pop() {
            pending.extend(
                neighbors[&id]
                    .iter()
                    .copied()
                    .filter(|&neighbor| visited.insert(neighbor)),
            );
        }
        visited.len() == neighbors.len()
    })
}

fn regular(uses: &Incidence) -> bool { uses.values().all(|count| count[0] <= 1 && count[1] <= 1) }

fn matches_chain(expected: &Incidence, actual: &Incidence) -> bool {
    regular(expected)
        && regular(actual)
        && expected
            .iter()
            .all(|(edge, count)| actual.get(edge).is_some_and(|actual| actual == count))
        && actual
            .iter()
            .all(|(edge, count)| expected.contains_key(edge) || *count == [1, 1])
}

fn finite(p: Point3) -> bool { [p.x, p.y, p.z].iter().all(|v| v.is_finite()) }

struct SampleTable {
    positions: Vec<Point3>,
    edges: Vec<Vec<usize>>,
}
impl SampleTable {
    fn new(
        shell: &CompressedShell<Point3, PolylineCurve, Option<PolygonMesh>>,
        tolerance: f64,
    ) -> Result<Self, MeshStitchError> {
        if !tolerance.is_finite() || tolerance <= 0.0 || !shell.vertices.iter().copied().all(finite)
        {
            return Err(MeshStitchError::InvalidSource);
        }
        let mut positions = shell.vertices.clone();
        let edges = shell
            .edges
            .iter()
            .map(|edge| {
                let (a, b) = edge.vertices;
                let front = *shell
                    .vertices
                    .get(a)
                    .ok_or(MeshStitchError::InvalidSource)?;
                let back = *shell
                    .vertices
                    .get(b)
                    .ok_or(MeshStitchError::InvalidSource)?;
                let samples = &edge.curve.0;
                if samples.len() < 2
                    || !samples.iter().copied().all(finite)
                    || samples[0].distance(front) > tolerance
                    || samples[samples.len() - 1].distance(back) > tolerance
                {
                    return Err(MeshStitchError::InvalidSource);
                }
                let mut ids = vec![a];
                for &p in &samples[1..samples.len() - 1] {
                    ids.push(positions.len());
                    positions.push(p);
                }
                ids.push(b);
                Ok(ids)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { positions, edges })
    }
}

#[derive(Clone, Copy)]
struct Candidate {
    id: usize,
    point: Point3,
}

struct FaceSamples {
    candidates: Vec<Candidate>,
    expected: Incidence,
}

impl FaceSamples {
    fn new(
        face: &CompressedFace<Option<PolygonMesh>>,
        edges: &[CompressedEdge<PolylineCurve>],
        table: &[Vec<usize>],
        positions: &[Point3],
    ) -> Option<Self> {
        let mut expected = Incidence::new();
        let mut candidates = Vec::new();
        let mut seen = BTreeSet::new();
        for edge_use in face.boundaries.iter().flatten() {
            let ids = table.get(edge_use.index)?;
            let edge = edges.get(edge_use.index)?;
            for pair in ids.windows(2) {
                let (a, b) = if edge_use.orientation == face.orientation {
                    (pair[0], pair[1])
                } else {
                    (pair[1], pair[0])
                };
                add_use(&mut expected, a, b);
            }
            if seen.insert(edge_use.index) {
                for (&id, &point) in ids.iter().zip(&edge.curve.0) {
                    candidates.push(Candidate { id, point });
                    // Include the original vertex coordinate as well as
                    // the small fitted curve endpoint residual.
                    candidates.push(Candidate {
                        id,
                        point: positions[id],
                    });
                }
            }
        }
        candidates.sort_unstable_by(|a, b| a.point.x.total_cmp(&b.point.x));
        Some(Self {
            candidates,
            expected,
        })
    }

    fn resolve(&self, p: Point3, positions: &[Point3], tolerance: f64) -> Option<usize> {
        let start = self
            .candidates
            .partition_point(|c| c.point.x < p.x - tolerance);
        let end = self
            .candidates
            .partition_point(|c| c.point.x <= p.x + tolerance);
        let mut exact = BTreeSet::new();
        let mut near = BTreeSet::new();
        for candidate in &self.candidates[start..end] {
            if p.distance(candidate.point) > tolerance
                || p.distance(positions[candidate.id]) > tolerance
            {
                continue;
            }
            if p == candidate.point {
                exact.insert(candidate.id);
            }
            near.insert(candidate.id);
        }
        let candidates = if exact.is_empty() { near } else { exact };
        (candidates.len() == 1)
            .then(|| candidates.first().copied())
            .flatten()
    }
}

struct Assembly {
    attributes: StandardAttributes,
    faces: Faces,
    incidence: Incidence,
}

impl Assembly {
    fn append_face(
        &mut self,
        index: usize,
        face: &CompressedFace<Option<PolygonMesh>>,
        samples: &FaceSamples,
        tolerance: f64,
    ) -> Result<(), MeshStitchError> {
        let mesh = face
            .surface
            .as_ref()
            .filter(|mesh| !mesh.faces().is_empty())
            .ok_or(MeshStitchError::MissingFace(index))?;
        let mut mapped =
            face_position_map(index, mesh, samples, &self.attributes.positions, tolerance)?;
        let uv_offset = self.attributes.uv_coords.len();
        let normal_offset = self.attributes.normals.len();
        self.attributes
            .uv_coords
            .extend_from_slice(mesh.uv_coords());
        self.attributes.normals.extend(
            mesh.normals()
                .iter()
                .map(|&n| if face.orientation { n } else { -n }),
        );
        let mut actual = Incidence::new();
        let initial_count = self.faces.len();
        for polygon in mesh.face_iter() {
            let mut converted = Vec::with_capacity(polygon.len());
            for &vertex in polygon {
                let id = if let Some(&id) = mapped.get(&vertex.pos) {
                    id
                } else {
                    let point = *mesh
                        .positions()
                        .get(vertex.pos)
                        .filter(|&&p| finite(p))
                        .ok_or(MeshStitchError::InvalidSource)?;
                    let id = self.attributes.positions.len();
                    self.attributes.positions.push(point);
                    mapped.insert(vertex.pos, id);
                    id
                };
                converted.push(StandardVertex {
                    pos: id,
                    uv: vertex.uv.map(|u| u + uv_offset),
                    nor: vertex.nor.map(|n| n + normal_offset),
                });
            }
            if !face.orientation {
                converted.reverse();
            }
            // Certified shared identities can collapse parameter-domain
            // pole edges. Their zero-area polygon chain cancels exactly.
            converted.dedup_by_key(|v| v.pos);
            if converted
                .first()
                .zip(converted.last())
                .is_some_and(|(a, b)| a.pos == b.pos)
            {
                converted.pop();
            }
            if converted.len() < 3 {
                continue;
            }
            if converted
                .iter()
                .map(|v| v.pos)
                .collect::<BTreeSet<_>>()
                .len()
                != converted.len()
            {
                return Err(MeshStitchError::BoundaryMismatch(index));
            }
            add_polygon(&mut actual, &converted);
            self.faces.push(converted);
        }
        if self.faces.len() == initial_count {
            return Err(MeshStitchError::MissingFace(index));
        }
        if !matches_chain(&samples.expected, &actual) {
            return Err(MeshStitchError::BoundaryMismatch(index));
        }
        for (edge, counts) in actual {
            let target = self.incidence.entry(edge).or_default();
            target[0] += counts[0];
            target[1] += counts[1];
        }
        Ok(())
    }
}

/// Join an already tessellated compressed shell using its native sample IDs.
///
/// Faces must emit the same oriented sample chain as their native edge uses.
/// A boundary vertex can move by at most `tolerance` to its unique original
/// sample position (including reconciled native edge endpoints). Unrelated
/// source identities are never welded, even at identical coordinates.
/// Interior positions and face-local normal/UV indices remain independent.
/// An open source shell remains open. This verifies edge incidence and vertex links, not
/// global geometric embedding or material volume.
///
/// # Errors
/// Returns a typed refusal for missing faces, unmatched or ambiguous
/// boundary positions, inconsistent edge incidence or singular vertices. Independently
/// resampled chart boundaries without native sample identities must
/// supply that topology before they can be stitched by this route.
pub fn stitch_compressed_shell(
    shell: &CompressedShell<Point3, PolylineCurve, Option<PolygonMesh>>,
    tolerance: f64,
) -> Result<PolygonMesh, MeshStitchError> {
    let SampleTable { positions, edges } = SampleTable::new(shell, tolerance)?;
    let mut assembly = Assembly {
        attributes: StandardAttributes {
            positions,
            ..Default::default()
        },
        faces: Faces::default(),
        incidence: Incidence::new(),
    };
    for (index, face) in shell.faces.iter().enumerate() {
        let samples = FaceSamples::new(face, &shell.edges, &edges, &assembly.attributes.positions)
            .ok_or(MeshStitchError::InvalidSource)?;
        assembly.append_face(index, face, &samples, tolerance)?;
    }
    if !regular(&assembly.incidence) || !manifold_links(&assembly.faces) {
        return Err(MeshStitchError::InvalidShell);
    }
    PolygonMesh::new(assembly.attributes, assembly.faces)
        .map_err(|_| MeshStitchError::InvalidSource)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_topology::compress::CompressedEdgeIndex;

    fn two_triangles() -> Result<
        CompressedShell<Point3, PolylineCurve, Option<PolygonMesh>>,
        Box<dyn std::error::Error>,
    > {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let middle = Point3::new(0.5, 0.5, 0.0);
        let edges = [[0, 1], [1, 2], [2, 0], [1, 3], [3, 2]]
            .into_iter()
            .enumerate()
            .map(|(i, [a, b])| {
                let points = if i == 1 {
                    vec![vertices[a], middle, vertices[b]]
                } else {
                    vec![vertices[a], vertices[b]]
                };
                CompressedEdge {
                    vertices: (a, b),
                    curve: monstertruck_mesh::PolylineCurve(points),
                }
            })
            .collect();
        let create = |uses: &[(usize, bool)],
                      positions: Vec<Point3>,
                      polygons: &[[usize; 3]]|
         -> Result<_, Box<dyn std::error::Error>> {
            Ok(CompressedFace {
                boundaries: vec![
                    uses.iter()
                        .map(|&(index, orientation)| CompressedEdgeIndex { index, orientation })
                        .collect(),
                ],
                orientation: true,
                surface: Some(PolygonMesh::new(
                    StandardAttributes {
                        positions,
                        ..Default::default()
                    },
                    Faces::from_iter(polygons.iter().copied()),
                )?),
            })
        };
        Ok(CompressedShell {
            vertices: vertices.clone(),
            edges,
            faces: vec![
                create(
                    &[(0, true), (1, true), (2, true)],
                    vec![vertices[0], vertices[1], middle, vertices[2]],
                    &[[0, 1, 2], [0, 2, 3]],
                )?,
                create(
                    &[(3, true), (4, true), (1, false)],
                    vec![vertices[1], vertices[3], vertices[2], middle],
                    &[[0, 1, 3], [1, 2, 3]],
                )?,
            ],
            vertex_stable_ids: None,
            edge_stable_ids: None,
            face_stable_ids: None,
        })
    }

    #[test]
    fn shares_native_samples_and_refuses_missing_boundary_segment()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut shell = two_triangles()?;
        let mesh = stitch_compressed_shell(&shell, 1.0e-9)?;
        assert_eq!(mesh.positions().len(), 5);
        assert_eq!(mesh.faces().len(), 4);
        let mut uses = Incidence::new();
        for polygon in mesh.face_iter() {
            add_polygon(&mut uses, polygon);
        }
        assert_eq!(uses[&[1, 4]], [1, 1]);
        assert_eq!(uses[&[2, 4]], [1, 1]);
        shell.faces[1].surface = Some(PolygonMesh::new(
            StandardAttributes {
                positions: shell.vertices[1..].to_vec(),
                ..Default::default()
            },
            Faces::from_iter([[0, 2, 1]]),
        )?);
        assert_eq!(
            stitch_compressed_shell(&shell, 1.0e-9),
            Err(MeshStitchError::BoundaryMismatch(1))
        );
        Ok(())
    }

    #[test]
    fn rejects_disconnected_vertex_link() {
        let faces = Faces::from_iter([[0, 1, 2], [0, 3, 4]]);
        assert!(!manifold_links(&faces));
        let faces = Faces::from_iter([[0, 1, 2], [0, 2, 3]]);
        assert!(manifold_links(&faces));
    }

    #[test]
    fn never_welds_coincident_native_components() -> Result<(), Box<dyn std::error::Error>> {
        let mut shell = two_triangles()?;
        let second = shell.clone();
        let offset = shell.vertices.len();
        let edge_offset = shell.edges.len();
        shell.vertices.extend(second.vertices);
        shell.edges.extend(second.edges.into_iter().map(|mut e| {
            e.vertices = (e.vertices.0 + offset, e.vertices.1 + offset);
            e
        }));
        shell.faces.extend(second.faces.into_iter().map(|mut f| {
            for e in f.boundaries.iter_mut().flatten() {
                e.index += edge_offset;
            }
            f
        }));
        let mesh = stitch_compressed_shell(&shell, 1.0e-9)?;
        assert_eq!(mesh.positions().len(), 10);
        assert_eq!(mesh.faces().len(), 8);
        Ok(())
    }
}
