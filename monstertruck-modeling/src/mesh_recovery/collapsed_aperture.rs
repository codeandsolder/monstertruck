//! Reconstruct source-supported openings of shallow prismatic cavities whose
//! mouth was collapsed into a zero-area contact graph by tessellation.
//! This is an explicit topology hypothesis, not a proof of the original CAD
//! topology. Retain the original mesh and every original polygon identity.

use super::{
    AnalyticRecoveryError, EdgeIncidence, MeshRecoveryError, NativeComponentReport,
    NativeMeshRecovery, NativeRecoveryOptions, input_polygon_from_ids, partition_recovery_mesh,
    recover_native_mesh, signed_polygon_volume, source_edge_incidence, welded_source_polygons,
};
use crate::{InnerSpace, KERNEL_EPSILON, Point3, Vector3};
use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::{BTreeMap, BTreeSet};

mod planar;
mod prism;

/// Provenance of each polygon in a proposed aperture reconstruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AperturePolygonOrigin {
    /// Unchanged directed polygon from the original mesh.
    Source(usize),
    /// Retriangulation of the surrounding planar support.
    ReconstructedCap,
    /// Reconstructed wall extending the observed floor outline.
    ReconstructedWall,
}

/// Local geometric evidence for a candidate aperture reconstruction.
#[derive(Clone, Debug)]
pub struct CollapsedApertureEvidence {
    /// Original cap polygons replaced by a planar domain with an opening.
    pub cap_polygons: Vec<usize>,
    /// Original cavity wall polygons replaced by a prismatic strip.
    pub wall_polygons: Vec<usize>,
    /// Original floor polygons retained without modification.
    pub floor_polygons: Vec<usize>,
    /// Coincident contact edges in the exact-coordinate-welded source.
    pub contact_edges: Vec<[usize; 2]>,
    /// Ordered original floor boundary in the exact-welded source.
    pub floor_ring: Vec<usize>,
    /// Depth of the observed floor from the surrounding plane.
    pub depth: f64,
    /// Local bidirectional source/candidate geometric deviation allowance.
    /// This must be reserved before further approximate recovery. It is not
    /// a certificate of the original CAD's missing aperture geometry.
    pub deviation_upper_bound: f64,
}

/// A possible reconstruction; the input and alternative topology remain
/// available to callers. No reference CAD face IDs are required.
#[derive(Debug)]
pub struct CollapsedApertureCandidate {
    /// Candidate geometry with original unaffected polygon cycles retained.
    mesh: PolygonMesh,
    /// Candidate face-order provenance (including generated faces).
    polygon_origins: Vec<AperturePolygonOrigin>,
    /// Geometric support and local approximation allowance.
    evidence: CollapsedApertureEvidence,
}

impl CollapsedApertureCandidate {
    /// Candidate geometry, with the original input still owned by the caller.
    #[must_use = "this is the source of the alternative recovery"]
    pub const fn mesh(&self) -> &PolygonMesh { &self.mesh }

    /// Candidate-face provenance in the candidate mesh's polygon order.
    #[must_use = "use the mapping to compose component-local provenance"]
    pub fn polygon_origins(&self) -> &[AperturePolygonOrigin] { &self.polygon_origins }

    /// Local geometric proof and explicit approximation allowance.
    #[must_use = "the allowance is required for subsequent recovery"]
    pub const fn evidence(&self) -> &CollapsedApertureEvidence { &self.evidence }

    fn remaining_options(
        &self,
        options: &NativeRecoveryOptions,
    ) -> Result<NativeRecoveryOptions, AnalyticRecoveryError> {
        let analytic = options.analytic_policy()?;
        let original = analytic
            .source_chord_budget
            .ok_or(AnalyticRecoveryError::InvalidOptions)?;
        let repair = self.evidence.deviation_upper_bound;
        let remaining = (original - repair).next_down();
        if remaining <= 0.0 {
            return Err(AnalyticRecoveryError::InvalidOptions);
        }
        let mut policy = options.clone();
        policy.analytic = analytic;
        policy.analytic.source_chord_budget = Some(remaining);
        if let Some(hybrid) = &mut policy.hybrid {
            if hybrid.budget.vertex < repair {
                return Err(AnalyticRecoveryError::InvalidOptions);
            }
            hybrid.budget.chord = remaining;
            hybrid.budget.vertex = (hybrid.budget.vertex - repair).next_down().max(0.0);
        }
        Ok(policy)
    }

    /// Recover this explicit alternative while reserving the local repair
    /// allowance from both subsequent chord and vertex budgets. Observer
    /// polygon IDs refer to the candidate mesh; compose them with
    /// `polygon_origins` for original-source provenance.
    ///
    /// # Errors
    /// Missing, inconsistent or insufficient geometric budgets, or a
    /// native recovery ingestion/accounting failure.
    pub fn recover_native(
        &self,
        options: &NativeRecoveryOptions,
        observe: impl FnMut(NativeComponentReport<'_>),
    ) -> Result<NativeMeshRecovery, AnalyticRecoveryError> {
        recover_native_mesh(&self.mesh, &self.remaining_options(options)?, observe)
    }
}

#[derive(Clone, Copy)]
struct Frame {
    origin: Point3,
    normal: Vector3,
    u: Vector3,
    v: Vector3,
    orientation: f64,
}

impl Frame {
    fn new(source: &Source, pair: &Contact, seed: usize) -> Option<Self> {
        let [a, b] = pair.edges[0];
        let normal = input_polygon_from_ids(
            &source.positions,
            source.polygons[seed].clone(),
            seed,
            KERNEL_EPSILON,
        )
        .ok()?
        .normal;
        let edge = source.positions[b] - source.positions[a];
        let tangent = edge - normal * edge.dot(normal);
        if tangent.magnitude() <= KERNEL_EPSILON {
            return None;
        }
        let u = tangent.normalize();
        Some(Self {
            origin: source.positions[a],
            normal,
            u,
            v: normal.cross(u),
            orientation: source.volumes[pair.outer]?.signum(),
        })
    }

    fn coordinates(self, p: Point3) -> [f64; 3] {
        let delta = p - self.origin;
        [
            delta.dot(self.u),
            delta.dot(self.v),
            -delta.dot(self.normal) * self.orientation,
        ]
    }

    fn lift(self, p: Point3) -> Point3 {
        p + self.normal * (self.coordinates(p)[2] * self.orientation)
    }
}

struct Source {
    positions: Vec<Point3>,
    polygons: Vec<Vec<usize>>,
    incidence: EdgeIncidence,
    components: Vec<Vec<usize>>,
    owner: Vec<usize>,
    volumes: Vec<Option<f64>>,
}

struct Contact {
    outer: usize,
    inner: usize,
    edges: Vec<[usize; 2]>,
}

impl Source {
    fn new(mesh: &PolygonMesh) -> Result<Self, MeshRecoveryError> {
        let (positions, polygons) = welded_source_polygons(mesh)?;
        let incidence = source_edge_incidence(polygons.iter().map(Vec::as_slice));
        let parts = partition_recovery_mesh(mesh)?;
        let mut components = Vec::with_capacity(parts.len());
        let mut volumes = Vec::with_capacity(parts.len());
        let mut owner = vec![usize::MAX; polygons.len()];
        for (id, part) in parts.into_iter().enumerate() {
            for &polygon in &part.original_polygons {
                owner[polygon] = id;
            }
            volumes.push(
                part.independently_closed
                    .then(|| {
                        signed_polygon_volume(
                            &positions,
                            part.original_polygons
                                .iter()
                                .map(|&i| polygons[i].as_slice()),
                        )
                    })
                    .flatten(),
            );
            components.push(part.original_polygons);
        }
        Ok(Self {
            positions,
            polygons,
            incidence,
            components,
            owner,
            volumes,
        })
    }

    fn contacts(&self) -> Vec<Contact> {
        let mut groups = BTreeMap::<(usize, usize), Vec<[usize; 2]>>::new();
        for (&edge, uses) in &self.incidence {
            let [a, b, c, d] = uses.as_slice() else {
                continue;
            };
            let mut owners = [
                (self.owner[a.0], a.1),
                (self.owner[b.0], b.1),
                (self.owner[c.0], c.1),
                (self.owner[d.0], d.1),
            ];
            owners.sort_unstable();
            let [(a, false), (b, true), (c, false), (d, true)] = owners else {
                continue;
            };
            if a != b || c != d || a == c {
                continue;
            }
            let (Some(va), Some(vc)) = (self.volumes[a], self.volumes[c]) else {
                continue;
            };
            if va.signum() == vc.signum() || va == 0.0 || vc == 0.0 {
                continue;
            }
            let key = if va.abs() > vc.abs() {
                (a, c)
            } else if vc.abs() > va.abs() {
                (c, a)
            } else {
                continue;
            };
            groups.entry(key).or_default().push(edge);
        }
        groups
            .into_iter()
            .map(|((outer, inner), mut edges)| {
                edges.sort_unstable();
                Contact {
                    outer,
                    inner,
                    edges,
                }
            })
            .collect()
    }

    fn cap(
        &self,
        pair: &Contact,
        coords: &[[f64; 3]],
        projected: &[[f64; 2]],
    ) -> Option<Vec<usize>> {
        let seeds = pair
            .edges
            .iter()
            .flat_map(|edge| {
                self.incidence[edge]
                    .iter()
                    .filter_map(|&(id, _)| (self.owner[id] == pair.outer).then_some(id))
            })
            .collect::<BTreeSet<_>>();
        let coplanar = |id: usize| {
            let face = &self.polygons[id];
            face.iter().all(|&v| coords[v][2].abs() <= KERNEL_EPSILON)
                && planar::area_twice(projected, face) > 0.0
                && planar::simple(projected, face)
        };
        if !seeds.iter().all(|&id| coplanar(id)) {
            return None;
        }
        let seed = *seeds.first()?;
        let mut cap = FxHashSet::default();
        cap.insert(seed);
        let mut queue = vec![seed];
        let mut cursor = 0;
        while let Some(&id) = queue.get(cursor) {
            for (a, b) in planar::edges(&self.polygons[id]) {
                let key = [a.min(b), a.max(b)];
                for &(other, _) in &self.incidence[&key] {
                    if self.owner[other] == pair.outer && !cap.contains(&other) && coplanar(other) {
                        cap.insert(other);
                        queue.push(other);
                    }
                }
            }
            cursor += 1;
        }
        if !seeds.iter().all(|id| cap.contains(id)) {
            return None;
        }
        queue.sort_unstable();
        Some(queue)
    }
}

struct Footprint {
    depth: f64,
    floor: Vec<usize>,
    wall: Vec<usize>,
    ring: Vec<usize>,
    contacts: BTreeSet<usize>,
}

impl Footprint {
    fn new(
        source: &Source,
        pair: &Contact,
        coords: &[[f64; 3]],
        projected: &[[f64; 2]],
    ) -> Option<Self> {
        let vertices = source.components[pair.inner]
            .iter()
            .flat_map(|&id| source.polygons[id].iter().copied())
            .collect::<BTreeSet<_>>();
        let depth = vertices
            .iter()
            .map(|&id| coords[id][2])
            .max_by(f64::total_cmp)?;
        if depth <= KERNEL_EPSILON || !depth.is_finite() {
            return None;
        }
        let contacts = pair
            .edges
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>();
        let top = vertices
            .iter()
            .copied()
            .filter(|&id| coords[id][2].abs() <= KERNEL_EPSILON)
            .collect::<BTreeSet<_>>();
        if top != contacts
            || vertices
                .iter()
                .any(|&id| !top.contains(&id) && (coords[id][2] - depth).abs() > KERNEL_EPSILON)
        {
            return None;
        }
        let (floor, wall): (Vec<_>, Vec<_>) = source.components[pair.inner]
            .iter()
            .copied()
            .partition(|&id| {
                source.polygons[id]
                    .iter()
                    .all(|&v| (coords[v][2] - depth).abs() <= KERNEL_EPSILON)
            });
        if floor.is_empty()
            || wall.is_empty()
            || !floor.iter().all(|&id| {
                planar::area_twice(projected, &source.polygons[id]) > 0.0
                    && planar::simple(projected, &source.polygons[id])
            })
        {
            return None;
        }
        let loops = planar::boundary(floor.iter().map(|&id| source.polygons[id].as_slice()))?;
        let [ring]: [Vec<usize>; 1] = loops.try_into().ok()?;
        if !planar::simple(projected, &ring) || planar::area_twice(projected, &ring) <= 0.0 {
            return None;
        }
        let sides = planar::edges(&ring).collect::<Vec<_>>();
        // Convex floor coverage allows a direct local Hausdorff bound via
        // projection onto the unchanged floor. Collinear samples survive.
        if sides.iter().any(|&(a, b)| {
            ring.iter()
                .any(|&id| planar::turn(projected[a], projected[b], projected[id]) < 0.0)
        }) {
            return None;
        }
        let boundary = ring.iter().copied().collect::<FxHashSet<_>>();
        if wall.iter().any(|&id| {
            let p = &source.polygons[id];
            p.iter().all(|v| top.contains(v))
                || p.iter().any(|v| !top.contains(v) && !boundary.contains(v))
        }) {
            return None;
        }
        // All observed floor facets must lie inside the convex domain.
        // Consistent oriented incidence and one simple boundary then prove
        // the complete filled planar chain independently of triangulation.
        if floor
            .iter()
            .flat_map(|&id| &source.polygons[id])
            .any(|&id| {
                sides
                    .iter()
                    .any(|&(a, b)| planar::turn(projected[a], projected[b], projected[id]) < 0.0)
            })
        {
            return None;
        }
        Some(Self {
            depth,
            floor,
            wall,
            ring,
            contacts,
        })
    }

    fn mouth(
        &self,
        source: &Source,
        frame: Frame,
        positions: &mut Vec<Point3>,
    ) -> Option<Vec<usize>> {
        let mut used = BTreeSet::new();
        let mut mouth = Vec::with_capacity(self.ring.len());
        for &id in &self.ring {
            let p = frame.lift(source.positions[id]);
            let mut matches = self
                .contacts
                .iter()
                .copied()
                .filter(|&v| (source.positions[v] - p).magnitude() <= KERNEL_EPSILON);
            let vertex = if let Some(v) = matches.next() {
                if matches.next().is_some() || !used.insert(v) {
                    return None;
                }
                v
            } else {
                let v = positions.len();
                positions.push(p);
                v
            };
            mouth.push(vertex);
        }
        (used == self.contacts).then_some(mouth)
    }
}

impl Source {
    fn candidate(&self, pair: &Contact, budget: f64) -> Option<CollapsedApertureCandidate> {
        let seed = self.incidence[&pair.edges[0]]
            .iter()
            .find(|&&(id, _)| self.owner[id] == pair.outer)?
            .0;
        let frame = Frame::new(self, pair, seed)?;
        let coords = self
            .positions
            .iter()
            .map(|&p| frame.coordinates(p))
            .collect::<Vec<_>>();
        let mut projected = coords.iter().map(|&[x, y, _]| [x, y]).collect::<Vec<_>>();
        let footprint = Footprint::new(self, pair, &coords, &projected)?;
        // Both old and new wall surfaces lie over the unchanged filled
        // floor. Removed cap coverage is at most the floor depth away.
        // Charge the plane/projection/vertex reconciliation allowances too.
        let deviation = (footprint.depth + 4.0 * KERNEL_EPSILON).next_up();
        if deviation >= budget {
            return None;
        }
        let cap = self.cap(pair, &coords, &projected)?;
        let mut loops = planar::boundary(cap.iter().map(|&id| self.polygons[id].as_slice()))?;
        if !planar::valid_domain(&projected, &loops) {
            return None;
        }
        let mut positions = self.positions.clone();
        let mouth = footprint.mouth(self, frame, &mut positions)?;
        projected.extend(positions[self.positions.len()..].iter().map(|&p| {
            let [x, y, _] = frame.coordinates(p);
            [x, y]
        }));
        if !mouth
            .iter()
            .all(|&id| planar::inside_domain(&projected, &loops, projected[id]))
            || loops
                .iter()
                .any(|ring| !planar::disjoint(&projected, ring, &mouth))
        {
            return None;
        }
        // Refuse to modify an aperture if an unrelated surface intersects
        // the repaired prism, even if its own shell is otherwise closed.
        let cap_set = cap.iter().copied().collect::<FxHashSet<_>>();
        if self.polygons.iter().enumerate().any(|(id, polygon)| {
            self.owner[id] != pair.inner
                && !cap_set.contains(&id)
                && prism::intersects(
                    &coords,
                    polygon,
                    &projected,
                    &footprint.ring,
                    footprint.depth,
                )
        }) {
            return None;
        }
        let mut hole = mouth.clone();
        hole.reverse();
        loops.push(hole);
        let cap_triangles = planar::triangulate(&projected, &loops)?;
        let mut walls = Vec::with_capacity(mouth.len() * 2);
        for ((a, b), (c, d)) in planar::edges(&footprint.ring).zip(planar::edges(&mouth)) {
            walls.extend([[b, a, c], [b, c, d]]);
        }
        let replaced = cap
            .iter()
            .chain(&footprint.wall)
            .copied()
            .collect::<FxHashSet<_>>();
        let mut faces = self
            .polygons
            .iter()
            .enumerate()
            .filter(|(id, _)| !replaced.contains(id))
            .map(|(id, p)| (p.clone(), AperturePolygonOrigin::Source(id)))
            .collect::<Vec<_>>();
        faces.extend(
            cap_triangles
                .into_iter()
                .map(|p| (p.to_vec(), AperturePolygonOrigin::ReconstructedCap)),
        );
        faces.extend(
            walls
                .into_iter()
                .map(|p| (p.to_vec(), AperturePolygonOrigin::ReconstructedWall)),
        );
        // Verify closed directed incidence of the reconstructed pair itself,
        // without using unrelated components as repairs.
        let repaired_boundary = planar::boundary(
            faces
                .iter()
                .filter(|(_, origin)| match origin {
                    AperturePolygonOrigin::Source(id) => {
                        matches!(self.owner[*id],o if o==pair.outer || o==pair.inner)
                    }
                    _ => true,
                })
                .map(|(p, _)| p.as_slice()),
        )?;
        if !repaired_boundary.is_empty() {
            return None;
        }
        let (mesh, polygon_origins) = assemble(&positions, faces);
        Some(CollapsedApertureCandidate {
            mesh,
            polygon_origins,
            evidence: CollapsedApertureEvidence {
                cap_polygons: cap,
                wall_polygons: footprint.wall,
                floor_polygons: footprint.floor,
                contact_edges: pair.edges.clone(),
                floor_ring: footprint.ring,
                depth: footprint.depth,
                deviation_upper_bound: deviation,
            },
        })
    }
}

fn assemble(
    positions: &[Point3],
    mut faces: Vec<(Vec<usize>, AperturePolygonOrigin)>,
) -> (PolygonMesh, Vec<AperturePolygonOrigin>) {
    // Faces stores triangles, quads and general polygons separately. Match
    // that order explicitly before collecting provenance.
    faces.sort_by_key(|(p, _)| p.len().min(5));
    let members = faces
        .iter()
        .flat_map(|(p, _)| p.iter().copied())
        .collect::<BTreeSet<_>>();
    let remap = members
        .iter()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect::<FxHashMap<_, _>>();
    let points = members.iter().map(|&id| positions[id]).collect();
    let (polygons, origins): (Vec<_>, Vec<_>) = faces
        .into_iter()
        .map(|(p, origin)| {
            (
                p.into_iter().map(|id| remap[&id]).collect::<Vec<_>>(),
                origin,
            )
        })
        .unzip();
    (
        PolygonMesh::new_unchecked(
            StandardAttributes {
                positions: points,
                ..Default::default()
            },
            Faces::from_iter(polygons),
        ),
        origins,
    )
}

/// Generate explicit aperture-repair hypotheses, one per certified
/// opposite-oriented contact pair. A candidate is not committed to the
/// original source. Reserve its deviation before later analytic recovery.
///
/// This recognizes a fully observed convex planar floor with two-level
/// sidewall evidence and a collapsed mouth, not arbitrary missing patches.
/// The source chord allowance must be finite and positive.
///
/// # Errors
/// Invalid budget or invalid source geometry/indexing.
pub fn collapsed_aperture_candidates(
    mesh: &PolygonMesh,
    source_chord_budget: f64,
) -> Result<Vec<CollapsedApertureCandidate>, MeshRecoveryError> {
    if !source_chord_budget.is_finite() || source_chord_budget <= 0.0 {
        return Err(MeshRecoveryError::InvalidTolerance);
    }
    let source = Source::new(mesh)?;
    Ok(source
        .contacts()
        .iter()
        .filter_map(|pair| source.candidate(pair, source_chord_budget))
        .collect())
}

#[cfg(test)]
mod tests;
