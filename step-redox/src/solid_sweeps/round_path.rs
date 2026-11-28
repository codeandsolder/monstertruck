//! One source-support graph for open and closed constant-round sweeps.
//! A graph edge is a complete, source-authored circular section. Its
//! incidence, rather than a guessed tangency of infinite supports, defines
//! the spine order and every segment endpoint.

use super::{
    CylinderRun, DIR_TOL, GEOM_TOL_MM, RecoveredRoundSweep, RecoveredSweepSegment, TorusBend,
    TorusBoundaryCircle, canonical_cylinder, canonical_torus, classify_torus_boundary_circle,
    common_value, lex_point, near, point_on_torus, prove_cylinder_between_points, same_cylinder,
    same_torus, torus_major_direction,
};
use crate::brep::{CircleSupport, CurveSupport, OrientedEdgeUse, PlaneSupport, SurfaceSupport};
use crate::circular_edges::circle_edge_sweep;
use crate::circular_edges::normalize_circular_edge;
use crate::math3::{cross, distance, dot, mul, norm, normalize, point_to_unit_line_distance, sub};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES as TOL;
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::f64::consts::TAU;
use step_p21::ast::EntityInstance;

#[derive(Clone)]
enum Carrier {
    Cylinder(CylinderRun),
    Torus(TorusBend),
    Cap(PlaneSupport),
}

impl Carrier {
    fn from_face(face: &ReducedFace) -> Option<Self> {
        match face.surface {
            SurfaceSupport::Cylinder(c) => Some(Self::Cylinder(canonical_cylinder(c)?)),
            SurfaceSupport::Torus(t) => Some(Self::Torus(canonical_torus(t)?)),
            SurfaceSupport::Plane(p) => Some(Self::Cap(p)),
            _ => None,
        }
    }

    fn same(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Cylinder(a), Self::Cylinder(b)) => same_cylinder(a, b),
            (Self::Torus(a), Self::Torus(b)) => same_torus(a, b),
            (Self::Cap(a), Self::Cap(b)) => {
                1.0 - dot(a.normal, b.normal).abs() <= DIR_TOL
                    && dot(sub(a.origin_mm, b.origin_mm), a.normal).abs() <= GEOM_TOL_MM
            }
            _ => false,
        }
    }

    fn set_faces(&mut self, faces: Vec<u64>) {
        match self {
            Self::Cylinder(run) => run.face_ids = faces,
            Self::Torus(bend) => bend.face_ids = faces,
            Self::Cap(_) => {}
        }
    }

    fn radius(&self) -> Option<f64> {
        match self {
            Self::Cylinder(run) => Some(run.radius_mm),
            Self::Torus(bend) => Some(bend.minor_radius_mm),
            Self::Cap(_) => None,
        }
    }

    fn is_cap(&self) -> bool { matches!(self, Self::Cap(_)) }

    fn contains_section(&self, circle: CircleSupport, radius: f64) -> bool {
        if !near(circle.radius_mm, radius) {
            return false;
        }
        match self {
            Self::Cylinder(run) => {
                near(run.radius_mm, radius)
                    && point_to_unit_line_distance(circle.center_mm, run.origin_mm, run.direction)
                        <= GEOM_TOL_MM
                    && 1.0 - dot(circle.normal, run.direction).abs() <= DIR_TOL
            }
            Self::Torus(bend) => matches!(
                classify_torus_boundary_circle(&circle, bend),
                Some(TorusBoundaryCircle::CrossSection { .. })
            ),
            Self::Cap(plane) => {
                dot(sub(circle.center_mm, plane.origin_mm), plane.normal).abs() <= GEOM_TOL_MM
                    && 1.0 - dot(circle.normal, plane.normal).abs() <= DIR_TOL
            }
        }
    }
}

type EdgeUses<'a> = BTreeMap<u64, Vec<(usize, &'a OrientedEdgeUse)>>;

struct Partition<'a> {
    pieces: Vec<Carrier>,
    face_groups: Vec<usize>,
    edges: EdgeUses<'a>,
}

impl<'a> Partition<'a> {
    fn from_solid(solid: &'a ReducedSolid) -> Option<Self> {
        solid.validate_closed_topology()?;
        let faces = solid.faces().collect::<Vec<_>>();
        let mut carriers = faces
            .iter()
            .map(|face| Carrier::from_face(face))
            .collect::<Option<Vec<_>>>()?;
        let mut edges = EdgeUses::new();
        for (i, face) in faces.iter().enumerate() {
            if face.loops.is_empty() {
                return None;
            }
            for edge in face.loop_edges() {
                edges.entry(edge.edge_id).or_default().push((i, edge));
            }
        }
        let mut neighbors = vec![BTreeSet::new(); faces.len()];
        for uses in edges.values() {
            let [(a, _), (b, _)] = uses.as_slice() else {
                return None;
            };
            if a != b {
                neighbors[*a].insert(*b);
                neighbors[*b].insert(*a);
            }
        }
        let mut pieces = Vec::new();
        let mut membership = vec![None; faces.len()];
        for seed in 0..faces.len() {
            if membership[seed].is_some() {
                continue;
            }
            let group = pieces.len();
            membership[seed] = Some(group);
            let mut queue = VecDeque::from([seed]);
            let mut face_ids = Vec::new();
            while let Some(face) = queue.pop_front() {
                face_ids.push(faces[face].id);
                for &neighbor in &neighbors[face] {
                    if membership[neighbor].is_none() && carriers[seed].same(&carriers[neighbor]) {
                        membership[neighbor] = Some(group);
                        queue.push_back(neighbor);
                    }
                }
            }
            face_ids.sort_unstable();
            carriers[seed].set_faces(face_ids);
            pieces.push(carriers[seed].clone());
        }
        Some(Self {
            pieces,
            face_groups: membership.into_iter().collect::<Option<_>>()?,
            edges,
        })
    }

    fn interfaces(&self, radius: f64) -> Option<Vec<Interface>> {
        let mut rings = BTreeMap::<[usize; 2], Vec<CircleRing>>::new();
        for uses in self.edges.values() {
            let [(a, edge), (b, _)] = uses.as_slice() else {
                return None;
            };
            let mut groups = [self.face_groups[*a], self.face_groups[*b]];
            if groups[0] == groups[1] {
                continue;
            }
            groups.sort_unstable();
            let edge = normalize_circular_edge(edge, None)?;
            let CurveSupport::Circle(circle) = edge.support else {
                return None;
            };
            if !groups
                .iter()
                .all(|&group| self.pieces[group].contains_section(circle, radius))
            {
                return None;
            }
            let candidates = rings.entry(groups).or_default();
            if let Some(ring) = candidates
                .iter_mut()
                .find(|ring| same_circle(ring.circle, circle))
            {
                ring.edges.push(edge);
            } else {
                candidates.push(CircleRing {
                    circle,
                    edges: vec![edge],
                });
            }
        }
        rings
            .into_iter()
            .flat_map(|(groups, rings)| rings.into_iter().map(move |ring| (groups, ring)))
            .map(|(groups, ring)| {
                ring.prove_complete()?;
                Some(Interface {
                    groups,
                    circle: ring.circle,
                })
            })
            .collect()
    }
}

fn same_circle(a: CircleSupport, b: CircleSupport) -> bool {
    near(a.radius_mm, b.radius_mm)
        && distance(a.center_mm, b.center_mm) <= GEOM_TOL_MM
        && 1.0 - dot(a.normal, b.normal).abs() <= DIR_TOL
}

struct CircleRing {
    circle: CircleSupport,
    edges: Vec<OrientedEdgeUse>,
}

fn angle_about(normal: [f64; 3], x: [f64; 3], radial: [f64; 3]) -> f64 {
    dot(normal, cross(x, radial)).atan2(dot(x, radial))
}

fn angular_tolerance(radius: f64) -> f64 { TOL.source_sweep_angle_rad.max(GEOM_TOL_MM / radius) }

fn split_interval(start: f64, length: f64) -> [(f64, f64); 2] {
    let end = start + length;
    [(start, end.min(TAU)), (0.0, (end - TAU).max(0.0))]
}

impl CircleRing {
    fn prove_complete(&self) -> Option<()> {
        let circle = self.circle;
        let normal = normalize(circle.normal, GEOM_TOL_MM)?;
        let x = normalize(circle.x_direction, GEOM_TOL_MM)?;
        let tolerance = angular_tolerance(circle.radius_mm);
        let mut degrees = HashMap::<u64, usize>::new();
        let mut intervals = Vec::new();
        for edge in &self.edges {
            let CurveSupport::Circle(authored) = edge.support else {
                return None;
            };
            if !same_circle(circle, authored) {
                return None;
            }
            *degrees.entry(edge.start_vertex).or_default() += 1;
            *degrees.entry(edge.end_vertex).or_default() += 1;
            let sweep = circle_edge_sweep(edge, authored)? * dot(authored.normal, normal).signum();
            let point = if sweep >= 0.0 {
                edge.start_mm
            } else {
                edge.end_mm
            };
            let angle = angle_about(normal, x, sub(point, circle.center_mm)).rem_euclid(TAU);
            if !sweep.is_finite() || sweep.abs() <= TOL.source_sweep_curve_parameter {
                return None;
            }
            intervals.extend(
                split_interval(angle, sweep.abs())
                    .into_iter()
                    .filter(|(a, b)| b > a),
            );
        }
        if degrees.is_empty() || degrees.values().any(|&count| count != 2) {
            return None;
        }
        intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut end = 0.0;
        let mut discrepancy = 0.0;
        for (start, stop) in intervals {
            discrepancy += (start - end).abs();
            end = stop.max(end);
        }
        (discrepancy + (TAU - end).abs() <= tolerance).then_some(())
    }
}

#[derive(Clone, Copy)]
struct Interface {
    groups: [usize; 2],
    circle: CircleSupport,
}

impl Interface {
    fn other(self, group: usize) -> Option<usize> {
        let [a, b] = self.groups;
        if group == a {
            Some(b)
        } else if group == b {
            Some(a)
        } else {
            None
        }
    }
}

struct Graph {
    pieces: Vec<Carrier>,
    interfaces: Vec<Interface>,
    adjacency: Vec<Vec<usize>>,
    radius: f64,
}

impl Graph {
    fn from_partition(partition: Partition<'_>) -> Option<Self> {
        let radii = partition
            .pieces
            .iter()
            .filter_map(Carrier::radius)
            .collect::<Vec<_>>();
        let radius = common_value(&radii)?;
        if radius <= GEOM_TOL_MM {
            return None;
        }
        let interfaces = partition.interfaces(radius)?;
        let mut adjacency = vec![Vec::new(); partition.pieces.len()];
        for (index, interface) in interfaces.iter().enumerate() {
            for group in interface.groups {
                adjacency[group].push(index);
            }
        }
        if partition
            .pieces
            .iter()
            .zip(&adjacency)
            .any(|(piece, edges)| edges.len() != if piece.is_cap() { 1 } else { 2 })
        {
            return None;
        }
        Some(Self {
            pieces: partition.pieces,
            interfaces,
            adjacency,
            radius,
        })
    }

    fn start_options(&self) -> Option<(bool, Vec<(usize, usize)>)> {
        let caps = self
            .pieces
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.is_cap().then_some(i))
            .collect::<Vec<_>>();
        match caps.as_slice() {
            [] => {
                let (index, interface) = self
                    .interfaces
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| lex_point(a.circle.center_mm, b.circle.center_mm))?;
                Some((
                    true,
                    interface
                        .groups
                        .into_iter()
                        .map(|group| (index, group))
                        .collect(),
                ))
            }
            [a, b] => {
                let ia = *self.adjacency[*a].first()?;
                let ib = *self.adjacency[*b].first()?;
                let (cap, index) = if lex_point(
                    self.interfaces[ia].circle.center_mm,
                    self.interfaces[ib].circle.center_mm,
                )
                .is_le()
                {
                    (*a, ia)
                } else {
                    (*b, ib)
                };
                Some((false, vec![(index, self.interfaces[index].other(cap)?)]))
            }
            _ => None,
        }
    }

    fn walk(
        &self,
        initial: usize,
        first: usize,
        closed: bool,
    ) -> Option<Vec<(usize, usize, usize)>> {
        let mut seen = BTreeSet::new();
        let mut path = Vec::new();
        let mut group = first;
        let mut entry = initial;
        loop {
            if self.pieces[group].is_cap() {
                if closed {
                    return None;
                }
                break;
            }
            if !seen.insert(group) {
                return None;
            }
            let exit = *self.adjacency[group].iter().find(|&&id| id != entry)?;
            path.push((group, entry, exit));
            group = self.interfaces[exit].other(group)?;
            entry = exit;
            if closed && entry == initial && group == first {
                break;
            }
        }
        let expected = self.pieces.iter().filter(|piece| !piece.is_cap()).count();
        (seen.len() == expected && !path.is_empty()).then_some(path)
    }

    fn segments(
        &self,
        path: &[(usize, usize, usize)],
        entities: &[EntityInstance],
        index: &HashMap<u64, usize>,
    ) -> Option<Vec<SegmentEvidence>> {
        path.iter()
            .map(|&(group, entry, exit)| {
                let start = self.interfaces[entry].circle.center_mm;
                let end = self.interfaces[exit].circle.center_mm;
                match &self.pieces[group] {
                    Carrier::Cylinder(run) => {
                        prove_cylinder_between_points(
                            run,
                            start,
                            end,
                            self.radius,
                            entities,
                            index,
                        )?;
                        line_evidence(start, end)
                    }
                    Carrier::Torus(bend) => prove_bend_between(bend, start, end, entities, index),
                    Carrier::Cap(_) => None,
                }
            })
            .collect()
    }
}

#[derive(Clone)]
pub(super) struct SegmentEvidence {
    pub(super) motion: RecoveredSweepSegment,
    pub(super) tangent_in: [f64; 3],
    pub(super) tangent_out: [f64; 3],
    pub(super) residual: f64,
}

pub(super) fn line_evidence(start: [f64; 3], end: [f64; 3]) -> Option<SegmentEvidence> {
    let vector_mm = sub(end, start);
    let tangent = normalize(vector_mm, GEOM_TOL_MM)?;
    Some(SegmentEvidence {
        motion: RecoveredSweepSegment::Translation { vector_mm },
        tangent_in: tangent,
        tangent_out: tangent,
        residual: 0.0,
    })
}

pub(super) fn bend_evidence(
    bend: &TorusBend,
    start: [f64; 3],
    end: [f64; 3],
    angle_rad: f64,
) -> Option<SegmentEvidence> {
    let radial_in = sub(start, bend.center_mm);
    let radial_out = sub(end, bend.center_mm);
    let residual = (norm(radial_in) - bend.major_radius_mm)
        .abs()
        .max((norm(radial_out) - bend.major_radius_mm).abs())
        .max(dot(radial_in, bend.axis).abs())
        .max(dot(radial_out, bend.axis).abs());
    if residual > GEOM_TOL_MM
        || !angle_rad.is_finite()
        || angle_rad.abs() <= TOL.source_sweep_angle_rad
        || angle_rad.abs() >= TAU
    {
        return None;
    }
    let tangent_in = normalize(
        mul(cross(bend.axis, radial_in), angle_rad.signum()),
        GEOM_TOL_MM,
    )?;
    let tangent_out = normalize(
        mul(cross(bend.axis, radial_out), angle_rad.signum()),
        GEOM_TOL_MM,
    )?;
    Some(SegmentEvidence {
        motion: RecoveredSweepSegment::Rotation {
            axis_origin_mm: bend.center_mm,
            axis_direction: bend.axis,
            angle_rad,
        },
        tangent_in,
        tangent_out,
        residual,
    })
}

struct BendChart<'a> {
    bend: &'a TorusBend,
    start_radial: [f64; 3],
    angle: f64,
    tolerance: f64,
}

impl BendChart<'_> {
    fn station(&self, radial: [f64; 3]) -> Option<f64> {
        let at = (self.angle.signum() * angle_about(self.bend.axis, self.start_radial, radial))
            .rem_euclid(TAU);
        let at = if TAU - at <= self.tolerance { 0.0 } else { at };
        (at <= self.angle.abs() + self.tolerance).then_some(at)
    }

    fn edge_interval(
        &self,
        edge: &OrientedEdgeUse,
        circle: CircleSupport,
    ) -> Option<Option<(f64, f64)>> {
        for point in [edge.start_mm, edge.end_mm] {
            if !point_on_torus(point, self.bend) {
                return None;
            }
            self.station(torus_major_direction(point, self.bend)?)?;
        }
        match classify_torus_boundary_circle(&circle, self.bend)? {
            TorusBoundaryCircle::CrossSection { radial_direction } => {
                self.station(radial_direction)?;
                Some(None)
            }
            TorusBoundaryCircle::MajorTrace => {
                let sweep = circle_edge_sweep(edge, circle)?;
                let along =
                    sweep * self.angle.signum() * dot(circle.normal, self.bend.axis).signum();
                let point = if along >= 0.0 {
                    edge.start_mm
                } else {
                    edge.end_mm
                };
                let start = self.station(torus_major_direction(point, self.bend)?)?;
                let end = start + sweep.abs();
                (sweep.abs() > TOL.source_sweep_curve_parameter
                    && end <= self.angle.abs() + self.tolerance)
                    .then_some(Some((start, end)))
            }
        }
    }

    fn proves(&self, entities: &[EntityInstance], index: &HashMap<u64, usize>) -> Option<()> {
        let mut intervals = Vec::new();
        for &face in &self.bend.face_ids {
            for loop_ in crate::brep::face_loops(face, entities, index)? {
                for edge in loop_.edges {
                    let edge = normalize_circular_edge(&edge, None)?;
                    let CurveSupport::Circle(circle) = edge.support else {
                        return None;
                    };
                    if let Some(interval) = self.edge_interval(&edge, circle)? {
                        intervals.push(interval)
                    }
                }
            }
        }
        // Major traces may overlap because different minor-angle seams
        // witness the same interval. Their union must cover the entire bend.
        intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut covered = 0.0_f64;
        for (start, end) in intervals {
            if start > covered + self.tolerance {
                return None;
            }
            covered = covered.max(end);
        }
        ((covered - self.angle.abs()).abs() <= self.tolerance).then_some(())
    }
}

/// Prove the directed angular interval from actual source trim traces. Both
/// the minor and major arc are considered; a carrier alone cannot choose.
pub(super) fn prove_bend_between(
    bend: &TorusBend,
    start: [f64; 3],
    end: [f64; 3],
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<SegmentEvidence> {
    let radial_in = normalize(sub(start, bend.center_mm), GEOM_TOL_MM)?;
    let radial_out = normalize(sub(end, bend.center_mm), GEOM_TOL_MM)?;
    let principal = angle_about(bend.axis, radial_in, radial_out);
    let mut candidates = [principal.rem_euclid(TAU), -(-principal).rem_euclid(TAU)]
        .into_iter()
        .filter_map(|angle| {
            let evidence = bend_evidence(bend, start, end, angle)?;
            BendChart {
                bend,
                start_radial: radial_in,
                angle,
                tolerance: angular_tolerance(bend.major_radius_mm),
            }
            .proves(entities, index)?;
            Some(evidence)
        });
    let candidate = candidates.next()?;
    candidates.next().is_none().then_some(candidate)
}

pub(super) fn assemble(
    solid_id: u64,
    face_ids: Vec<u64>,
    radius: f64,
    start: [f64; 3],
    closed: bool,
    segments: Vec<SegmentEvidence>,
) -> Option<RecoveredRoundSweep> {
    let first = segments.first()?;
    let last = segments.last()?;
    let tangent_matches =
        |a: [f64; 3], b: [f64; 3]| dot(a, b) >= 1.0 - TOL.source_sweep_tangent_cosine;
    if segments
        .windows(2)
        .any(|pair| !tangent_matches(pair[0].tangent_out, pair[1].tangent_in))
        || (closed && !tangent_matches(last.tangent_out, first.tangent_in))
    {
        return None;
    }
    let tangent = first.tangent_in;
    let x = segments
        .iter()
        .find_map(|e| match e.motion {
            RecoveredSweepSegment::Rotation { axis_direction, .. }
                if dot(axis_direction, tangent).abs() <= DIR_TOL =>
            {
                Some(axis_direction)
            }
            _ => None,
        })
        .or_else(|| normalize(cross(tangent, [1.0, 0.0, 0.0]), GEOM_TOL_MM))
        .or_else(|| normalize(cross(tangent, [0.0, 1.0, 0.0]), GEOM_TOL_MM))?;
    let y = normalize(cross(tangent, x), GEOM_TOL_MM)?;
    let max_residual_mm = segments.iter().map(|e| e.residual).fold(0.0, f64::max);
    Some(RecoveredRoundSweep {
        solid_id,
        source_face_ids: face_ids,
        profile_radius_mm: radius,
        path_start_mm: start,
        profile_x_axis: x,
        profile_y_axis: y,
        path_tangent: tangent,
        segments: segments.into_iter().map(|e| e.motion).collect(),
        closed,
        max_residual_mm,
    })
}

pub(super) fn detect(
    solid_id: u64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<RecoveredRoundSweep> {
    let solid = ReducedSolid::from_step(solid_id, entities, index)?;
    let partition = Partition::from_solid(&solid)?;
    let graph = Graph::from_partition(partition)?;
    let (closed, starts) = graph.start_options()?;
    let face_ids = solid.faces().map(|face| face.id).collect::<Vec<_>>();
    starts.into_iter().find_map(|(initial, first)| {
        let path = graph.walk(initial, first, closed)?;
        let segments = graph.segments(&path, entities, index)?;
        assemble(
            solid_id,
            face_ids.clone(),
            graph.radius,
            graph.interfaces[initial].circle.center_mm,
            closed,
            segments,
        )
    })
}

#[cfg(test)]
mod tests;
