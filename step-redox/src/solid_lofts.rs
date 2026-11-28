//! Recover planar ruled lofts by following the source edge correspondence.
//! Whole solids and partial terminal features use the same wall-ring proof.
//! A section is not required to be rectangular or aligned to world axes.

use crate::brep::{CurveSupport, OrientedEdgeUse, PlaneSupport, SurfaceSupport};
use crate::cad_ir::{CadModel, CadNode, LoftSection, Profile2d, RigidTransform};
use crate::cad_recovery::{CadFragment, CadFragmentSource};
use crate::line_edges::line_chain;
use crate::math3::{add, cross, distance, dot, mul, normalize, sub};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use step_p21::ast::EntityInstance;

const TOL_MM: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIR_TOL: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.loft_direction;

#[derive(Clone, Copy)]
struct Vertex {
    id: u64,
    point: [f64; 3],
}

struct Ring {
    vertices: Vec<Vertex>,
    edge_ids: Vec<Vec<u64>>,
    level: f64,
}

#[derive(Clone, Copy)]
struct Frame {
    origin: [f64; 3],
    axis: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
}

impl Frame {
    fn from_cap(cap: &ReducedFace, planes: &BTreeMap<u64, Vec<PlaneSupport>>) -> Option<Self> {
        let SurfaceSupport::Plane(plane) = cap.surface else {
            return None;
        };
        let first = planar_loop(cap)?.first()?;
        let origin = planar_corner(planes.get(&first.start_vertex)?, first.start_mm)?;
        let axis = normalize(plane.normal, DIR_TOL)?;
        let CurveSupport::Line(line) = first.support else {
            return None;
        };
        let direction = normalize(line.direction, DIR_TOL)?;
        if dot(axis, direction).abs() > DIR_TOL {
            return None;
        }
        let y = normalize(cross(axis, direction), DIR_TOL)?;
        let x = normalize(cross(y, axis), DIR_TOL)?;
        Some(Self { origin, axis, x, y })
    }

    fn level(self, point: [f64; 3]) -> f64 { dot(sub(point, self.origin), self.axis) }

    fn section(self, ring: &Ring) -> Option<LoftSection> {
        let points = ring
            .vertices
            .iter()
            .map(|vertex| {
                let local = sub(vertex.point, self.origin);
                [dot(local, self.x), dot(local, self.y)]
            })
            .collect::<Vec<_>>();
        let origin = add(self.origin, mul(self.axis, ring.level));
        Some(LoftSection {
            profile: Profile2d::polygon(points).ok()?,
            transform: RigidTransform {
                matrix: [
                    [self.x[0], self.y[0], self.axis[0], origin[0]],
                    [self.x[1], self.y[1], self.axis[1], origin[1]],
                    [self.x[2], self.y[2], self.axis[2], origin[2]],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        })
    }
}

/// A planar polygon is a geometric condition, not a particular cap arity or
/// STEP outer-bound spelling. A wall slot is a quad, unless a source edge
/// has been split or a carrier continues beyond the recovered section.
fn planar_loop(face: &ReducedFace) -> Option<&[OrientedEdgeUse]> {
    let SurfaceSupport::Plane(plane) = face.surface else {
        return None;
    };
    let [loop_] = face.loops.as_slice() else {
        return None;
    };
    let edges = &loop_.edges;
    (edges.len() >= 3
        && edges.iter().all(|edge| {
            matches!(edge.support, CurveSupport::Line(_))
                && [edge.start_mm, edge.end_mm].into_iter().all(|point| {
                    point.iter().all(|value| value.is_finite())
                        && dot(sub(point, plane.origin_mm), plane.normal).abs() <= TOL_MM
                })
        }))
    .then_some(edges)
}

impl Ring {
    fn from_cap(cap: &ReducedFace, frame: Frame, topology: &Topology<'_>) -> Option<Self> {
        let edges = planar_loop(cap)?;
        let consumed = BTreeSet::from([cap.id]);
        let neighbors = edges
            .iter()
            .map(|edge| Some(topology.neighbor(edge.edge_id, &consumed)?.id))
            .collect::<Option<Vec<_>>>()?;
        let n = edges.len();
        // Begin at a real wall transition so wraparound subdivisions are
        // assigned to the same geometric chain.
        let first = (0..n).find(|&i| neighbors[i] != neighbors[(i + n - 1) % n])?;
        let mut chains = Vec::<Vec<&OrientedEdgeUse>>::new();
        let mut previous = None;
        for i in (0..n).map(|offset| (first + offset) % n) {
            if previous == Some(neighbors[i]) {
                chains.last_mut()?.push(&edges[i]);
            } else {
                chains.push(vec![&edges[i]]);
            }
            previous = Some(neighbors[i]);
        }
        let vertices = chains
            .iter()
            .map(|chain| {
                line_chain(chain)?;
                let edge = chain.first()?;
                Some(Vertex {
                    id: edge.start_vertex,
                    point: planar_corner(topology.planes.get(&edge.start_vertex)?, edge.start_mm)?,
                })
            })
            .collect::<Option<Vec<_>>>()?;
        Self::new(
            vertices,
            chains
                .iter()
                .map(|chain| chain.iter().map(|edge| edge.edge_id).collect())
                .collect(),
            frame,
        )
    }

    fn new(vertices: Vec<Vertex>, edge_ids: Vec<Vec<u64>>, frame: Frame) -> Option<Self> {
        let level = frame.level(vertices.first()?.point);
        if vertices.len() < 3
            || vertices.iter().map(|v| v.id).collect::<BTreeSet<_>>().len() != vertices.len()
            || edge_ids.len() != vertices.len()
            || edge_ids.iter().any(Vec::is_empty)
            || edge_ids
                .iter()
                .flatten()
                .copied()
                .collect::<BTreeSet<_>>()
                .len()
                != edge_ids.iter().map(Vec::len).sum::<usize>()
            || vertices
                .iter()
                .any(|v| (frame.level(v.point) - level).abs() > TOL_MM)
        {
            return None;
        }
        Some(Self {
            vertices,
            edge_ids,
            level,
        })
    }
}

/// Index the source adjacency once. Looking up an incident wall never scans
/// every face, and ambiguity is rejected rather than resolved by face order.
struct Topology<'a> {
    solid: &'a ReducedSolid,
    edge_faces: HashMap<u64, Vec<&'a ReducedFace>>,
    planes: BTreeMap<u64, Vec<PlaneSupport>>,
}

impl<'a> Topology<'a> {
    fn new(solid: &'a ReducedSolid) -> Self {
        let mut edge_faces = HashMap::<u64, Vec<&ReducedFace>>::new();
        for face in solid.faces() {
            for edge in face.loop_edges() {
                edge_faces.entry(edge.edge_id).or_default().push(face);
            }
        }
        Self {
            solid,
            edge_faces,
            planes: incident_planar_supports(solid),
        }
    }

    fn neighbor(&self, edge: u64, consumed: &BTreeSet<u64>) -> Option<&'a ReducedFace> {
        let mut faces = self
            .edge_faces
            .get(&edge)?
            .iter()
            .copied()
            .filter(|face| !consumed.contains(&face.id));
        let face = faces.next()?;
        faces.next().is_none().then_some(face)
    }

    fn closing_cap(&self, ring: &Ring, consumed: &BTreeSet<u64>, frame: Frame) -> Option<u64> {
        let face = self.neighbor(*ring.edge_ids.first()?.first()?, consumed)?;
        let SurfaceSupport::Plane(plane) = face.surface else {
            return None;
        };
        let edges = planar_loop(face)?;
        (dot(normalize(plane.normal, DIR_TOL)?, frame.axis).abs() >= 1.0 - DIR_TOL
            && edges
                .iter()
                .map(|edge| edge.edge_id)
                .collect::<BTreeSet<_>>()
                == ring
                    .edge_ids
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<BTreeSet<_>>())
        .then_some(face.id)
    }

    fn walls(&self, ring: &Ring, consumed: &BTreeSet<u64>) -> Option<Vec<&'a ReducedFace>> {
        let walls = ring
            .edge_ids
            .iter()
            .map(|chain| {
                let first = self.neighbor(*chain.first()?, consumed)?;
                chain
                    .iter()
                    .all(|&id| {
                        self.neighbor(id, consumed)
                            .is_some_and(|face| face.id == first.id)
                    })
                    .then_some(first)
            })
            .collect::<Option<Vec<_>>>()?;
        if walls
            .iter()
            .map(|face| face.id)
            .collect::<BTreeSet<_>>()
            .len()
            != walls.len()
        {
            return None;
        }
        Some(walls)
    }
}

/// One source wall supplies the opposite edge and the two corresponding
/// generator endpoints. No geometric sorting or nearest-point matching.
struct WallSlot {
    far_edges: Vec<u64>,
    start: Vertex,
    end: Vertex,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum WallRegion {
    Near,
    PreviousGenerator,
    NextGenerator,
    Far,
}

fn wall_region(
    edge: u64,
    near: &BTreeSet<u64>,
    previous: &BTreeSet<u64>,
    next: &BTreeSet<u64>,
) -> Option<WallRegion> {
    match (
        near.contains(&edge),
        previous.contains(&edge),
        next.contains(&edge),
    ) {
        (true, false, false) => Some(WallRegion::Near),
        (false, true, false) => Some(WallRegion::PreviousGenerator),
        (false, false, true) => Some(WallRegion::NextGenerator),
        (false, false, false) => Some(WallRegion::Far),
        _ => None,
    }
}

/// Partition an authored wall loop by actual incidence with the near ring
/// and adjacent walls. Each of the four regions is a complete straight
/// chain, not necessarily a single edge.
fn wall_slot(
    wall: &ReducedFace,
    previous_wall: &ReducedFace,
    next_wall: &ReducedFace,
    near_edges: &[u64],
    start: Vertex,
    end: Vertex,
) -> Option<WallSlot> {
    let edges = planar_loop(wall)?;
    let near = near_edges.iter().copied().collect::<BTreeSet<_>>();
    let previous = previous_wall
        .loop_edges()
        .map(|edge| edge.edge_id)
        .collect::<BTreeSet<_>>();
    let next = next_wall
        .loop_edges()
        .map(|edge| edge.edge_id)
        .collect::<BTreeSet<_>>();
    let regions = edges
        .iter()
        .map(|edge| wall_region(edge.edge_id, &near, &previous, &next))
        .collect::<Option<Vec<_>>>()?;
    let n = edges.len();
    let first = (0..n).find(|&i| {
        regions[i] == WallRegion::Near && regions[(i + n - 1) % n] != WallRegion::Near
    })?;
    let mut chains = Vec::<(WallRegion, Vec<&OrientedEdgeUse>)>::new();
    for i in (0..n).map(|offset| (first + offset) % n) {
        if let Some((_, chain)) = chains
            .last_mut()
            .filter(|(region, _)| *region == regions[i])
        {
            chain.push(&edges[i]);
        } else {
            chains.push((regions[i], vec![&edges[i]]));
        }
    }
    let [
        (WallRegion::Near, near_chain),
        (first_generator, _),
        (WallRegion::Far, far_chain),
        (last_generator, _),
    ] = chains.as_slice()
    else {
        return None;
    };
    if near_chain.len() != near.len() {
        return None;
    }
    for (_, chain) in &chains {
        line_chain(chain)?;
    }
    let near_first = near_chain.first()?;
    let near_last = near_chain.last()?;
    let far_first = far_chain.first()?;
    let far_last = far_chain.last()?;
    let far_start = Vertex {
        id: far_first.start_vertex,
        point: far_first.start_mm,
    };
    let far_end = Vertex {
        id: far_last.end_vertex,
        point: far_last.end_mm,
    };
    let (start, end) = if near_first.start_vertex == start.id
        && near_last.end_vertex == end.id
        && *first_generator == WallRegion::NextGenerator
        && *last_generator == WallRegion::PreviousGenerator
    {
        (far_end, far_start)
    } else if near_first.start_vertex == end.id
        && near_last.end_vertex == start.id
        && *first_generator == WallRegion::PreviousGenerator
        && *last_generator == WallRegion::NextGenerator
    {
        (far_start, far_end)
    } else {
        return None;
    };
    Some(WallSlot {
        far_edges: far_chain.iter().map(|edge| edge.edge_id).collect(),
        start,
        end,
    })
}

/// Refine section corners using the actual incident wall equations, measured
/// axial station and previous exact corner. This is the general form of
/// propagating drafted wall slopes; it avoids amplifying STEP vertex noise
/// into a short wall's angle. A collinear split has rank two, so retain its
/// source witness instead of inventing an additional constraint.
fn propagate_corner(
    previous: [f64; 3],
    witness: [f64; 3],
    before: PlaneSupport,
    after: PlaneSupport,
    frame: Frame,
    delta: f64,
) -> Option<[f64; 3]> {
    let displacement = intersect_three_planes(
        (frame.axis, mul(frame.axis, delta)),
        (before.normal, [0.0; 3]),
        (after.normal, [0.0; 3]),
    );
    let point = displacement.map_or(witness, |offset| add(previous, offset));
    (distance(point, witness) <= TOL_MM
        && [before, after].iter().all(|plane| {
            [previous, point]
                .into_iter()
                .all(|p| dot(sub(p, plane.origin_mm), plane.normal).abs() <= TOL_MM)
        }))
    .then_some(point)
}

fn wall_planes(walls: &[&ReducedFace]) -> Option<Vec<PlaneSupport>> {
    walls
        .iter()
        .map(|face| match face.surface {
            SurfaceSupport::Plane(plane) => Some(plane),
            _ => None,
        })
        .collect()
}

fn next_ring(
    current: &Ring,
    walls: &[&ReducedFace],
    topology: &Topology<'_>,
    frame: Frame,
) -> Option<Ring> {
    let n = current.vertices.len();
    let slots = walls
        .iter()
        .enumerate()
        .map(|(i, wall)| {
            wall_slot(
                wall,
                walls[(i + n - 1) % n],
                walls[(i + 1) % n],
                &current.edge_ids[i],
                current.vertices[i],
                current.vertices[(i + 1) % n],
            )
        })
        .collect::<Option<Vec<_>>>()?;
    let witnesses = slots
        .iter()
        .enumerate()
        .map(|(i, slot)| {
            let previous = &slots[(i + n - 1) % n].end;
            let v = slot.start;
            if previous.id != v.id || distance(previous.point, v.point) > TOL_MM {
                return None;
            }
            Some(Vertex {
                id: v.id,
                point: planar_corner(topology.planes.get(&v.id)?, v.point)?,
            })
        })
        .collect::<Option<Vec<_>>>()?;
    let level = frame.level(witnesses.first()?.point);
    let delta = level - current.level;
    if delta.abs() <= TOL_MM
        || witnesses
            .iter()
            .any(|v| (frame.level(v.point) - level).abs() > TOL_MM)
    {
        return None;
    }
    let planes = wall_planes(walls)?;
    let vertices = witnesses
        .into_iter()
        .enumerate()
        .map(|(i, vertex)| {
            let point = propagate_corner(
                current.vertices[i].point,
                vertex.point,
                planes[(i + n - 1) % n],
                planes[i],
                frame,
                delta,
            )?;
            Some(Vertex { point, ..vertex })
        })
        .collect::<Option<Vec<_>>>()?;
    Ring::new(
        vertices,
        slots.into_iter().map(|slot| slot.far_edges).collect(),
        frame,
    )
}

/// When a wall continues into the core, its authored plane and the nearest
/// actual incident non-ring vertex witness a candidate section. This is not
/// a source-face ownership proof; the caller must certify the resulting
/// Boolean against the original solid. It uses the same cyclic wall order
/// and section construction as an exact ring transition.
fn extended_ring(current: &Ring, walls: &[&ReducedFace], frame: Frame) -> Option<Ring> {
    let vertices = current
        .vertices
        .iter()
        .map(|v| v.id)
        .collect::<BTreeSet<_>>();
    let delta = walls
        .iter()
        .flat_map(|wall| wall.loop_edges())
        .flat_map(|edge| {
            [
                (edge.start_vertex, edge.start_mm),
                (edge.end_vertex, edge.end_mm),
            ]
        })
        .filter(|(id, _)| !vertices.contains(id))
        .map(|(_, p)| frame.level(p) - current.level)
        .filter(|depth| depth.abs() > TOL_MM && depth.is_finite())
        .min_by(|a, b| a.abs().total_cmp(&b.abs()))?;
    let level = current.level + delta;
    let planes = wall_planes(walls)?;
    let n = planes.len();
    let vertices = current
        .vertices
        .iter()
        .enumerate()
        .map(|(i, vertex)| {
            let prior = planes[(i + n - 1) % n];
            let next = planes[i];
            let point = add(
                vertex.point,
                intersect_three_planes(
                    (frame.axis, mul(frame.axis, delta)),
                    (prior.normal, [0.0; 3]),
                    (next.normal, [0.0; 3]),
                )?,
            );
            point
                .iter()
                .all(|v| v.is_finite())
                .then_some(Vertex { point, ..*vertex })
        })
        .collect::<Option<Vec<_>>>()?;
    Some(Ring {
        vertices,
        edge_ids: Vec::new(),
        level,
    })
}

struct LoftWalk {
    frame: Frame,
    sections: Vec<LoftSection>,
    wall_ids: Vec<u64>,
    closing_cap: Option<u64>,
    depth: f64,
}

impl LoftWalk {
    fn model(self, source: CadFragmentSource) -> Option<CadFragment> {
        let mut model = CadModel::new();
        let root = model.add_node(CadNode::Loft {
            sections: self.sections,
        });
        model.add_root(root).ok()?;
        model.validate().ok()?;
        Some(CadFragment {
            source,
            model,
            root,
        })
    }
}

/// Walk from a source cap through every witnessed wall ring. Each step
/// consumes previously unseen source faces, giving a topology-derived
/// termination invariant, not an arbitrary number of stages.
fn walk_loft(topology: &Topology<'_>, cap: &ReducedFace, allow_extended: bool) -> Option<LoftWalk> {
    let frame = Frame::from_cap(cap, &topology.planes)?;
    let mut ring = Ring::from_cap(cap, frame, topology)?;
    let mut sections = vec![frame.section(&ring)?];
    let mut consumed = BTreeSet::from([cap.id]);
    let mut wall_ids = Vec::new();
    let mut direction = 0.0_f64;
    let mut closing_cap = None;
    let mut extended = false;
    loop {
        if let Some(id) = topology.closing_cap(&ring, &consumed, frame) {
            closing_cap = Some(id);
            break;
        }
        let Some(walls) = topology.walls(&ring, &consumed) else {
            break;
        };
        let next = if let Some(next) = next_ring(&ring, &walls, topology, frame) {
            next
        } else if allow_extended && let Some(next) = extended_ring(&ring, &walls, frame) {
            extended = true;
            next
        } else {
            break;
        };
        let delta = next.level - ring.level;
        if delta.abs() <= TOL_MM || (direction != 0.0 && delta.signum() != direction.signum()) {
            break;
        }
        direction = delta;
        wall_ids.extend(walls.iter().map(|face| face.id));
        consumed.extend(walls.iter().map(|face| face.id));
        sections.push(frame.section(&next)?);
        ring = next;
        if extended {
            break;
        }
    }
    if sections.len() < 2 {
        return None;
    }
    let depth = ring.level;
    if depth < 0.0 {
        sections.reverse();
    }
    Some(LoftWalk {
        frame,
        sections,
        wall_ids,
        closing_cap,
        depth,
    })
}

/// Recognize a complete planar ruled loft using source cap-to-cap wall
/// correspondence. A general polygon and any number of stages are supported.
pub fn detect_sectioned_planar_loft(solid: &ReducedSolid) -> Option<CadFragment> {
    solid.validate_closed_topology()?;
    let topology = Topology::new(solid);
    solid
        .faces()
        .filter_map(|cap| {
            let walk = walk_loft(&topology, cap, false)?;
            let other = walk.closing_cap?;
            if walk.wall_ids.len() + 2 != solid.faces().count() {
                return None;
            }
            let span = walk.depth.abs();
            // Keep cap ordering canonical in the loft's positive axial direction.
            let cap_face_ids = if walk.depth < 0.0 {
                [other, cap.id]
            } else {
                [cap.id, other]
            };
            let side_face_ids = solid
                .faces()
                .filter_map(|f| (!cap_face_ids.contains(&f.id)).then_some(f.id))
                .collect();
            let source = CadFragmentSource::SolidLoft {
                solid_id: solid.solid_id(),
                cap_face_ids,
                side_face_ids,
            };
            Some((span, walk.model(source)?))
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(_, fragment)| fragment)
}

/// An exact or source-witnessed planar end treatment. The shoulder section
/// is an implied interface, not an invented source face.
pub struct RuledEndExtension {
    pub fragment: CadFragment,
    pub tip_face_id: u64,
    pub axis: [f64; 3],
    pub shoulder_point_mm: [f64; 3],
    pub tip_point_mm: [f64; 3],
}

fn terminal_extensions(solid: &ReducedSolid, allow_extended: bool) -> Vec<RuledEndExtension> {
    let topology = Topology::new(solid);
    solid
        .faces()
        .filter_map(|cap| {
            let walk = walk_loft(&topology, cap, allow_extended)?;
            if walk.closing_cap.is_some() {
                return None;
            }
            let frame = walk.frame;
            let shoulder_point_mm = add(frame.origin, mul(frame.axis, walk.depth));
            let source = CadFragmentSource::RuledEndExtension {
                solid_id: topology.solid.solid_id(),
                tip_face_id: cap.id,
                wall_face_ids: walk.wall_ids.clone(),
            };
            Some(RuledEndExtension {
                fragment: walk.model(source)?,
                tip_face_id: cap.id,
                axis: frame.axis,
                shoulder_point_mm,
                tip_point_mm: frame.origin,
            })
        })
        .collect()
}

/// Terminal features with fully witnessed source wall rings.
pub fn detect_ruled_end_extensions(solid: &ReducedSolid) -> Vec<RuledEndExtension> {
    terminal_extensions(solid, false)
}

/// Terminal features including an implied section on continuing source
/// carriers. The proposal still requires independent Boolean proof.
pub fn detect_extended_ruled_end_extensions(solid: &ReducedSolid) -> Vec<RuledEndExtension> {
    terminal_extensions(solid, true)
}

/// Index the actual analytic plane equations incident to each source vertex.
fn incident_planar_supports(solid: &ReducedSolid) -> BTreeMap<u64, Vec<PlaneSupport>> {
    let mut result = BTreeMap::<u64, Vec<PlaneSupport>>::new();
    for face in solid.faces() {
        let SurfaceSupport::Plane(plane) = face.surface else {
            continue;
        };
        let vertices = face
            .loop_edges()
            .flat_map(|edge| [edge.start_vertex, edge.end_vertex])
            .collect::<BTreeSet<_>>();
        for vertex in vertices {
            result.entry(vertex).or_default().push(plane);
        }
    }
    result
}

/// When three independent source planes determine a corner, intersect their
/// authored equations instead of amplifying rounded STEP vertex errors into
/// the normals of short loft walls. Solve in coordinates local to the witness.
fn planar_corner(planes: &[PlaneSupport], witness: [f64; 3]) -> Option<[f64; 3]> {
    let mut best = None;
    let mut best_det = DIR_TOL;
    for (i, first) in planes.iter().enumerate() {
        for (j, second) in planes.iter().enumerate().skip(i + 1) {
            for third in planes.iter().skip(j + 1) {
                let determinant = dot(first.normal, cross(second.normal, third.normal)).abs();
                if determinant > best_det {
                    best_det = determinant;
                    best = Some([first, second, third]);
                }
            }
        }
    }
    // Rank-deficient incidences do not determine a replacement point.
    let point = if let Some([first, second, third]) = best {
        let local = |plane: &PlaneSupport| (plane.normal, sub(plane.origin_mm, witness));
        crate::math3::add(
            witness,
            intersect_three_planes(local(first), local(second), local(third))?,
        )
    } else {
        witness
    };
    (distance(point, witness) <= TOL_MM
        && planes
            .iter()
            .all(|plane| dot(sub(point, plane.origin_mm), plane.normal).abs() <= TOL_MM))
    .then_some(point)
}

fn intersect_three_planes(
    a: ([f64; 3], [f64; 3]),
    b: ([f64; 3], [f64; 3]),
    c: ([f64; 3], [f64; 3]),
) -> Option<[f64; 3]> {
    let ab = cross(a.0, b.0);
    let bc = cross(b.0, c.0);
    let ca = cross(c.0, a.0);
    let det = dot(a.0, bc);
    if det.abs() <= DIR_TOL {
        return None;
    }
    let out = crate::math3::add(
        crate::math3::add(
            crate::math3::mul(bc, dot(a.0, a.1)),
            crate::math3::mul(ca, dot(b.0, b.1)),
        ),
        crate::math3::mul(ab, dot(c.0, c.1)),
    );
    Some(crate::math3::mul(out, 1.0 / det))
}

/// Read-only parser witness for any original STEP solid; the same immutable
/// reduced-B-rep recognizer also consumes the kernel-owned native B-rep.
pub fn detect_sectioned_planar_lofts(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    solid_ids: &[u64],
) -> Vec<CadFragment> {
    solid_ids
        .iter()
        .copied()
        .filter_map(|id| {
            let reduced = ReducedSolid::from_step(id, entities, index)?;
            detect_sectioned_planar_loft(&reduced)
        })
        .collect()
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
fn recover_proven_test_loft(
    step: &str,
    kernel: &crate::cad_kernel::monstertruck::MonstertruckKernel,
) -> anyhow::Result<CadFragment> {
    use crate::step_graph::build_index;
    use monstertruck_io::step::load::ParsedStepSource;
    let source = ParsedStepSource::parse(step.as_bytes())?.exchange;
    let entities = &source.data[0].entities;
    let index = build_index(entities);
    let ids = crate::brep::manifold_solid_ids(entities);
    let [id] = ids.as_slice() else {
        anyhow::bail!("expected a single source solid");
    };
    let reduced = ReducedSolid::from_step(*id, entities, &index)
        .ok_or_else(|| anyhow::anyhow!("source B-rep missing"))?;
    let recovered = detect_sectioned_planar_loft(&reduced)
        .ok_or_else(|| anyhow::anyhow!("source ring correspondence not recognized"))?;
    let proof = crate::cad_recovery::SourceRecovery::new(*id, entities, &index, kernel)
        .ok_or_else(|| anyhow::anyhow!("source boundary certificate unavailable"))?;
    assert!(proof.proves(&recovered.model, recovered.root, None)?);
    Ok(recovered)
}

#[cfg(test)]
mod terminal_extension_tests {
    use super::*;
    use crate::step_graph::build_index;
    use monstertruck_io::step::load::ParsedStepSource;

    #[test]
    fn authored_plane_corner_refines_noisy_vertex_without_inventing_rank() -> anyhow::Result<()> {
        let exact = [67.0, 5.0, -9.0];
        let witness = crate::math3::add(exact, [2e-9, -3e-9, 4e-9]);
        let plane = |normal, origin_mm| PlaneSupport {
            normal,
            origin_mm,
            max_residual_mm: 0.0,
        };
        let supports = [
            plane([1.0, 0.0, 0.0], exact),
            plane([0.0, 1.0, 0.0], exact),
            plane([0.0, 0.0, 1.0], exact),
            plane(
                [
                    std::f64::consts::FRAC_1_SQRT_2,
                    std::f64::consts::FRAC_1_SQRT_2,
                    0.0,
                ],
                exact,
            ),
        ];
        let refined = planar_corner(&supports, witness).ok_or_else(|| {
            anyhow::anyhow!("independent source planes must determine the corner")
        })?;
        assert!(distance(refined, exact) < 1e-14);
        let rank_deficient = planar_corner(&supports[..2], witness)
            .ok_or_else(|| anyhow::anyhow!("rank-deficient incidence should retain the witness"))?;
        assert!(distance(rank_deficient, witness) < f64::EPSILON);
        let inconsistent = [
            supports[0],
            supports[1],
            supports[2],
            plane(
                [1.0, 0.0, 0.0],
                crate::math3::add(exact, [10.0 * TOL_MM, 0.0, 0.0]),
            ),
        ];
        assert!(planar_corner(&inconsistent, witness).is_none());
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn two_section_drafted_rectangle_with_unlabelled_bounds_is_exact() -> anyhow::Result<()> {
        use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
        let kernel = MonstertruckKernel::project_standard();
        let mut model = CadModel::new();
        let root = model.add_node(CadNode::Loft {
            sections: vec![
                LoftSection {
                    profile: Profile2d::polygon(vec![
                        [-0.75, 0.0],
                        [0.75, 0.0],
                        [0.75, 0.9],
                        [-0.75, 0.9],
                    ])?,
                    transform: RigidTransform::identity(),
                },
                LoftSection {
                    profile: Profile2d::polygon(vec![
                        [-0.738, 0.012],
                        [0.738, 0.012],
                        [0.738, 0.888],
                        [-0.738, 0.888],
                    ])?,
                    transform: RigidTransform::translation_mm([0.0, 0.0, 1.35]),
                },
            ],
        });
        model.add_root(root)?;
        let authored = kernel.evaluate(&model, root)?;
        let step = kernel
            .to_step(&authored)?
            .replace("FACE_OUTER_BOUND(", "FACE_BOUND(");
        let recovered = recover_proven_test_loft(&step, &kernel)?;
        assert!(
            matches!(recovered.model.node(recovered.root)?, CadNode::Loft { sections } if sections.len() == 2)
        );
        Ok(())
    }

    #[test]
    fn four_section_beveled_bar_is_an_exact_editable_loft() -> anyhow::Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_four_section_beveled_bar.step");
        let source = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &source.data[0].entities;
        let reduced = ReducedSolid::from_step(24, entities, &build_index(entities))
            .ok_or_else(|| anyhow::anyhow!("no source B-rep"))?;
        let recovered = detect_sectioned_planar_loft(&reduced)
            .ok_or_else(|| anyhow::anyhow!("real four-section beveled bar not recovered"))?;
        assert!(matches!(
            recovered.source,
            CadFragmentSource::SolidLoft { solid_id: 24, .. }
        ));
        recovered.model.validate()?;
        #[cfg(feature = "cad-kernel-monstertruck")]
        {
            use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
            let kernel = MonstertruckKernel::project_standard();
            let evaluated = kernel.evaluate(&recovered.model, recovered.root)?;
            assert!(kernel.summarize(&evaluated).is_single_solid());
            let reconstructed = kernel.to_step(&evaluated)?;
            let output = ParsedStepSource::parse(reconstructed.as_bytes())?.exchange;
            let out_entities = &output.data[0].entities;
            let solid_ids = out_entities
                .iter()
                .filter_map(|e| {
                    crate::step_graph::simple_record(e)
                        .filter(|record| record.name == "MANIFOLD_SOLID_BREP")
                        .map(|_| crate::step_graph::entity_id(e))
                })
                .collect::<Vec<_>>();
            let [output_id] = solid_ids.as_slice() else {
                anyhow::bail!("ruled loft did not produce exactly one STEP solid");
            };
            let boundary = crate::boundary_proof::compare_analytic_boundaries(
                24,
                entities,
                *output_id,
                out_entities,
            );
            assert!(
                boundary.equivalent,
                "source vs ruled loft boundary differs: {boundary:?}"
            );
            let report = crate::complete_ir::recover_complete_ir_with_kernel_bytes(bytes, &kernel)?;
            assert_eq!(report.body_count, 1);
            assert_eq!(
                report.unrecovered_solid_count, 0,
                "{:?}",
                report.kernel_recovery_failures
            );
            let [chosen] = report.solid_lofts.as_slice() else {
                anyhow::bail!("complete IR did not select the source-exact four-section loft");
            };
            assert_eq!(chosen.fragment.source, recovered.source);
        }
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn source_4a6356_exposes_exact_multistage_terminal_loft() -> anyhow::Result<()> {
        use crate::cad_kernel::{CadKernel, monstertruck::MonstertruckKernel};
        let bytes = include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");
        let exchange = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &exchange.data[0].entities;
        let solid = ReducedSolid::from_step(171, entities, &build_index(entities))
            .ok_or_else(|| anyhow::anyhow!("source mixed carrier B-rep is missing"))?;
        let candidates = detect_extended_ruled_end_extensions(&solid);
        let end =
            candidates
                .iter()
                .find(|extension| {
                    extension.fragment.model.nodes.iter().any(
                        |node| matches!(node, CadNode::Loft { sections } if sections.len() == 4),
                    )
                })
                .ok_or_else(|| anyhow::anyhow!("failed to find source-derived 4-section lobe"))?;
        let CadFragmentSource::RuledEndExtension { wall_face_ids, .. } = &end.fragment.source
        else {
            anyhow::bail!("unexpected extension provenance");
        };
        assert_eq!(wall_face_ids.len(), 12);
        let kernel = MonstertruckKernel::project_standard();
        let result = kernel.evaluate(&end.fragment.model, end.fragment.root)?;
        assert!(kernel.summarize(&result).is_single_solid());
        if let Ok(output) = std::env::var("STEP_PARTIAL_LOFT_OUTPUT") {
            std::fs::write(output, kernel.to_step(&result)?)?;
        }
        Ok(())
    }

    #[test]
    fn discovers_exact_planar_beveled_ends_of_mixed_analytic_elbow() -> anyhow::Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_beveled_l_hook.step");
        let source = ParsedStepSource::parse(bytes)?.exchange;
        let entities = &source.data[0].entities;
        let reduced = ReducedSolid::from_step(84, entities, &build_index(entities))
            .ok_or_else(|| anyhow::anyhow!("source B-rep is missing"))?;
        let mut ends = detect_ruled_end_extensions(&reduced)
            .into_iter()
            .filter_map(|extension| match extension.fragment.source {
                CadFragmentSource::RuledEndExtension { wall_face_ids, .. } => {
                    Some((extension.tip_face_id, wall_face_ids))
                }
                _ => None,
            })
            .collect::<Vec<_>>();
        ends.sort_unstable_by_key(|entry| entry.0);
        assert_eq!(
            ends.len(),
            2,
            "both source end treatments must be reconstructed"
        );
        assert_eq!(ends[0].0, 330);
        assert_eq!(ends[1].0, 355);
        Ok(())
    }
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod polygon_tests;
