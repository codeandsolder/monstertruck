//! Simultaneous analytic rolling-ball construction on trivalent planar edges.
//!
//! Each original corner owns its contact vertices. The selected edge owns
//! the two cylindrical rails and cross-sections, and a fully rounded corner
//! owns the spherical patch connecting the three cross-sections. This is
//! independent of source-edge iteration order and does not weld by proximity.

use crate::{Curve, Edge, EdgeId, Face, Shell, Surface, Vertex, VertexId, Wire, builder};
use monstertruck_geometry::{
    analytic_curve_incidence::CurveIncidence, analytic_line_segment::SegmentLocus, prelude::*,
};
use rustc_hash::{FxHashMap as Map, FxHashSet as Set};

mod tangent_rim;

const ANGULAR: f64 = monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES.angular;

#[derive(Clone, Copy)]
struct FaceUse {
    face: usize,
    forward: bool,
}

struct SourceEdge {
    edge: Edge,
    uses: Vec<FaceUse>,
}

struct Incidence {
    edges: Map<EdgeId, SourceEdge>,
    at_vertex: Map<VertexId, Vec<EdgeId>>,
}

impl Incidence {
    fn new(shell: &Shell) -> Self {
        let mut edges = Map::<EdgeId, SourceEdge>::default();
        for (face, original) in shell.iter().enumerate() {
            for edge in original.edge_iter() {
                edges
                    .entry(edge.id())
                    .or_insert_with(|| SourceEdge {
                        edge: edge.absolute_clone(),
                        uses: Vec::new(),
                    })
                    .uses
                    .push(FaceUse {
                        face,
                        forward: edge.orientation(),
                    });
            }
        }
        let mut at_vertex = Map::<VertexId, Vec<EdgeId>>::default();
        for (&id, edge) in &edges {
            for vertex in [edge.edge.absolute_front(), edge.edge.absolute_back()] {
                at_vertex.entry(vertex.id()).or_default().push(id);
            }
        }
        Self { edges, at_vertex }
    }
}

fn line_axis(edge: &Edge) -> Option<Vector3> {
    let start = edge.absolute_front().point();
    let end = edge.absolute_back().point();
    let delta = end - start;
    (delta.norm() > KERNEL_EPSILON
        && crate::certify_native_curve_segment_locus(&edge.curve(), start, end, KERNEL_EPSILON)
            == SegmentLocus::WithinTolerance)
        .then(|| delta.normalize())
}

fn normal(shell: &Shell, face: usize) -> Option<Vector3> {
    let Surface::Plane(p) = shell.get(face)?.oriented_surface() else {
        return None;
    };
    let normal = p.normal();
    normal.norm().near(&1.0).then_some(normal)
}

#[derive(Clone, Copy)]
enum PlanarTreatment {
    Round,
    Chamfer,
    CarrierChamfer {
        first_normal: Vector3,
        second_setback: f64,
    },
}

/// For planar carriers, equal transverse setbacks are related to the
/// auxiliary rolling-circle radius by setback = radius * tan(angle/2).
fn chamfer_radius(normals: [Vector3; 2], setback: f64) -> Option<f64> {
    let [a, b] = normals;
    let sine = a.cross(&b).norm();
    let cosine = a.dot(&b);
    if sine <= ANGULAR || !setback.is_finite() || setback <= KERNEL_EPSILON {
        return None;
    }
    let radius = setback * (1.0 + cosine) / sine;
    (radius.is_finite() && radius > KERNEL_EPSILON).then_some(radius)
}

struct SelectedEdge {
    id: EdgeId,
    radius: f64,
    axis: Vector3,
    faces: [FaceUse; 2],
    normals: [Vector3; 2],
    setbacks: Option<[f64; 2]>,
    convex: bool,
}

impl SelectedEdge {
    fn new(
        shell: &Shell,
        incidence: &Incidence,
        edge_id: EdgeId,
        size: f64,
        profile: PlanarTreatment,
    ) -> Option<Self> {
        if !size.is_finite() || size <= KERNEL_EPSILON {
            return None;
        }
        let original = incidence.edges.get(&edge_id)?;
        let axis = line_axis(&original.edge)?;
        let &[a, b] = original.uses.as_slice() else {
            return None;
        };
        if a.face == b.face || a.forward == b.forward {
            return None;
        }
        let na = normal(shell, a.face)?;
        let nb = normal(shell, b.face)?;
        let forward = if a.forward { axis } else { -axis };
        // The interior of face A lies to the left of its outward coedge.
        // Its relation to carrier B distinguishes convex and reentrant edges.
        let bend = nb.dot(&na.cross(&forward));
        if bend.abs() <= ANGULAR || na.dot(&axis).abs() > ANGULAR || nb.dot(&axis).abs() > ANGULAR {
            return None;
        }
        let normals = [na, nb];
        let (radius, setbacks) = match profile {
            PlanarTreatment::Round => (size, None),
            PlanarTreatment::Chamfer => (chamfer_radius(normals, size)?, Some([size; 2])),
            PlanarTreatment::CarrierChamfer {
                first_normal,
                second_setback,
            } => {
                if !second_setback.is_finite()
                    || second_setback <= KERNEL_EPSILON
                    || !first_normal.norm().near(&1.0)
                {
                    return None;
                }
                let distances = normals.map(|n| n.cross(&first_normal).norm());
                let setbacks = match distances.map(|d| d <= ANGULAR) {
                    [true, false] => [size, second_setback],
                    [false, true] => [second_setback, size],
                    _ => return None,
                };
                (size, Some(setbacks))
            }
        };
        Some(Self {
            id: edge_id,
            radius,
            axis,
            faces: [a, b],
            normals,
            setbacks,
            convex: bend < 0.0,
        })
    }
}

enum CornerKind {
    Sphere,
    Cap {
        face: usize,
        normal: Vector3,
        axis: Vector3,
    },
}

struct Corner {
    center: Point3,
    radius: f64,
    kind: CornerKind,
    contacts: Map<usize, Vertex>,
}

fn intersection_displacement(normals: [Vector3; 3], offsets: [f64; 3]) -> Option<Vector3> {
    let [a, b, c] = normals;
    let det = a.dot(&b.cross(&c));
    if det.abs() <= ANGULAR {
        return None;
    }
    let displacement =
        (b.cross(&c) * offsets[0] + c.cross(&a) * offsets[1] + a.cross(&b) * offsets[2]) / det;
    normals
        .iter()
        .zip(offsets)
        .all(|(n, d)| (n.dot(&displacement) - d).abs() <= KERNEL_EPSILON)
        .then_some(displacement)
}

fn projected_radial(radial: Vector3, kind: &CornerKind) -> Vector3 {
    match *kind {
        CornerKind::Sphere => radial,
        CornerKind::Cap { normal, axis, .. } => {
            radial - axis * (normal.dot(&radial) / normal.dot(&axis))
        }
    }
}

fn validate_corner_seams(
    incidence: &Incidence,
    vertex: &Vertex,
    faces: [usize; 3],
    normals: [Vector3; 3],
    offsets: [f64; 3],
    sense: f64,
) -> Option<()> {
    let incident = incidence.at_vertex.get(&vertex.id())?;
    let point = vertex.point();
    // A spherical junction requires a convex trihedron. For one treated
    // edge, the terminal plane can instead be a reentrant shoulder: the
    // bevel/round enlarges its opening or fills its notch. The two
    // retained edges must still point into the treated carrier halfspaces,
    // and the later segment-trimming proof must keep their new endpoints
    // inside their original finite spans.
    for id in incident {
        let source = &incidence.edges[id];
        let [a, b] = source.uses.as_slice() else {
            return None;
        };
        if a.face == b.face || a.forward == b.forward {
            return None;
        }
        let other = if source.edge.absolute_front().id() == vertex.id() {
            source.edge.absolute_back()
        } else {
            source.edge.absolute_front()
        };
        let missing = faces
            .iter()
            .position(|&f| source.uses.iter().all(|u| u.face != f))?;
        let advance = normals[missing].dot(&(other.point() - point));
        if offsets[missing] == 0.0 {
            if advance.abs() <= KERNEL_EPSILON {
                return None;
            }
        } else if sense * advance >= -KERNEL_EPSILON {
            return None;
        }
        line_axis(&source.edge)?;
    }
    Some(())
}

impl Corner {
    fn new(
        shell: &Shell,
        incidence: &Incidence,
        selected: &Map<EdgeId, SelectedEdge>,
        vertex: &Vertex,
    ) -> Option<Self> {
        let incident = incidence.at_vertex.get(&vertex.id())?;
        if incident.len() != 3 {
            return None;
        }
        let mut faces = incident
            .iter()
            .flat_map(|id| incidence.edges[id].uses.iter().map(|u| u.face))
            .collect::<Vec<_>>();
        faces.sort_unstable();
        faces.dedup();
        let [a, b, c] = faces.as_slice() else {
            return None;
        };
        let faces = [*a, *b, *c];
        let normals = [normal(shell, *a)?, normal(shell, *b)?, normal(shell, *c)?];
        let rounded = incident
            .iter()
            .filter_map(|id| selected.get(id))
            .collect::<Vec<_>>();
        let radius = rounded.first()?.radius;
        let sense = if rounded.first()?.convex { 1.0 } else { -1.0 };
        let mut offsets = [-sense * radius; 3];
        let kind = match rounded.as_slice() {
            [edge] => {
                let cap_index = faces
                    .iter()
                    .position(|&f| edge.faces.iter().all(|u| u.face != f))?;
                let cap_normal = normals[cap_index];
                if cap_normal.dot(&edge.axis).abs() <= ANGULAR {
                    return None;
                }
                offsets[cap_index] = 0.0;
                CornerKind::Cap {
                    face: faces[cap_index],
                    normal: cap_normal,
                    axis: edge.axis,
                }
            }
            [first, second, third]
                if rounded
                    .iter()
                    .all(|edge| edge.setbacks.is_none() && edge.convex)
                    && [second.radius, third.radius]
                        .iter()
                        .all(|r| r.near(&first.radius)) =>
            {
                CornerKind::Sphere
            }
            _ => return None,
        };
        let point = vertex.point();
        validate_corner_seams(incidence, vertex, faces, normals, offsets, sense)?;
        let center = point + intersection_displacement(normals, offsets)?;
        let mut contacts = Map::default();
        for (i, &face) in faces.iter().enumerate() {
            if offsets[i] == 0.0 {
                continue;
            }
            let contact = if let [edge] = rounded.as_slice()
                && let Some(setbacks) = edge.setbacks
            {
                let slot = edge.faces.iter().position(|u| u.face == face)?;
                let n = edge.normals[slot];
                let other = edge.normals[1 - slot];
                let direction = -sense * (other - n * n.dot(&other)).normalize();
                point + projected_radial(direction * setbacks[slot], &kind)
            } else {
                center + projected_radial(sense * radius * normals[i], &kind)
            };
            if normals.iter().zip(offsets).any(|(n, offset)| {
                let distance = n.dot(&(contact - point));
                if offset == 0.0 {
                    distance.abs() > KERNEL_EPSILON
                } else {
                    sense * distance > KERNEL_EPSILON
                }
            }) {
                return None;
            }
            contacts.insert(face, Vertex::new(contact));
        }
        Some(Self {
            center,
            radius,
            kind,
            contacts,
        })
    }
}

/// Exact (possibly elliptical) cross-section of a cylindrical round.
/// The cosine/sine vectors include axial projection into an oblique cap.
struct Section {
    center: Point3,
    cosine: Vector3,
    sine: Vector3,
    angle: f64,
}

impl Section {
    fn new(corner: &Corner, edge: &SelectedEdge) -> Option<Self> {
        let [a, b] = edge.normals;
        let cosine = a.dot(&b);
        let tangent = b - a * cosine;
        let angle = tangent.norm().atan2(cosine);
        if angle <= ANGULAR || std::f64::consts::PI - angle <= ANGULAR {
            return None;
        }
        let sense = if edge.convex { 1.0 } else { -1.0 };
        Some(Self {
            center: corner.center,
            cosine: projected_radial(sense * edge.radius * a, &corner.kind),
            sine: projected_radial(sense * edge.radius * tangent.normalize(), &corner.kind),
            angle,
        })
    }

    fn curve(&self) -> Curve {
        let basis = self.cosine.cross(&self.sine).normalize();
        let transform = Matrix4::from_columns(&[
            self.cosine.push(0.0),
            self.sine.push(0.0),
            basis.push(0.0),
            self.center.to_homogeneous(),
        ]);
        Processor::with_transform(
            TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, self.angle)),
            transform,
        )
        .to_same_geometry()
    }

    fn edge(&self, corner: &Corner, edge: &SelectedEdge) -> Option<Edge> {
        let a = corner.contacts.get(&edge.faces[0].face)?;
        let b = corner.contacts.get(&edge.faces[1].face)?;
        let seam = Edge::new(a, b, self.curve()).ok()?;
        seam.is_geometric_consistent().then_some(seam)
    }

    fn axial_coefficients(&self, axis: Vector3) -> [f64; 2] {
        [axis.dot(&self.cosine), axis.dot(&self.sine)]
    }
}

/// Exact extrema of the axial projection on a bounded circle parameter.
fn trig_bounds(a: f64, b: f64, angle: f64) -> [f64; 2] {
    let endpoint = a.mul_add(angle.cos(), b * angle.sin());
    let mut low = a.min(endpoint);
    let mut high = a.max(endpoint);
    let peak = b.atan2(a).rem_euclid(std::f64::consts::TAU);
    let minimum = (peak + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU);
    let amplitude = a.hypot(b);
    if peak <= angle {
        high = amplitude;
    }
    if minimum <= angle {
        low = -amplitude;
    }
    [low, high]
}

fn certify_generated_face(face: &Face) -> bool {
    let surface = face.surface();
    face.edge_iter().all(|edge| {
        edge.is_geometric_consistent()
            && crate::certify_analytic_curve_on_surface(&edge.curve(), &surface, KERNEL_EPSILON)
                == CurveIncidence::Proven
    })
}

struct TreatedEdge {
    rails: [Edge; 2],
    seams: [Edge; 2],
    face: Face,
}

fn exact_cylinder(start: &Section, end: &Section, edge: &SelectedEdge) -> Option<Surface> {
    let axis = edge.axis;
    let center_delta = end.center - start.center;
    let span = center_delta.dot(&axis);
    if (center_delta - axis * span).norm() > KERNEL_EPSILON {
        return None;
    }
    let [a, b] = start.axial_coefficients(axis);
    let [c, d] = end.axial_coefficients(axis);
    // Cross-sections cannot pass through one another anywhere in the
    // angular interval, even when the terminating caps are oblique.
    if span + trig_bounds(c - a, d - b, start.angle)[0] <= KERNEL_EPSILON {
        return None;
    }
    let low = trig_bounds(a, b, start.angle)[0];
    let high = span + trig_bounds(c, d, start.angle)[1];
    let sense = if edge.convex { 1.0 } else { -1.0 };
    let radial = edge.normals[0] * (sense * edge.radius);
    let generator = Curve::Line(Line(
        start.center + radial + axis * low,
        start.center + radial + axis * high,
    ));
    let orbit_axis = edge.normals[0].cross(&edge.normals[1]).normalize();
    let mut surface = Surface::Revolution(Processor::new(RevolutionSurface::by_revolution(
        generator,
        start.center,
        orbit_axis,
    )));
    if surface.normal(0.5, 0.0).dot(&edge.normals[0]) < 0.0 {
        surface.invert();
    }
    Some(surface)
}

impl TreatedEdge {
    fn new(
        edge: &SelectedEdge,
        original: &Edge,
        corners: &Map<VertexId, Corner>,
        profile: PlanarTreatment,
    ) -> Option<Self> {
        let a = corners.get(&original.absolute_front().id())?;
        let b = corners.get(&original.absolute_back().id())?;
        let rails = edge
            .faces
            .map(|u| builder::line(&a.contacts[&u.face], &b.contacts[&u.face]));
        let (seams, surface) = match profile {
            PlanarTreatment::Round => {
                let start = Section::new(a, edge)?;
                let end = Section::new(b, edge)?;
                (
                    [start.edge(a, edge)?, end.edge(b, edge)?],
                    Some(exact_cylinder(&start, &end, edge)?),
                )
            }
            PlanarTreatment::Chamfer | PlanarTreatment::CarrierChamfer { .. } => {
                // The bevel is a plane. Both rails must retain a positive
                // axial interval even with oblique termination planes.
                if rails.iter().any(|rail| {
                    (rail.back().point() - rail.front().point()).dot(&edge.axis) <= KERNEL_EPSILON
                }) {
                    return None;
                }
                let [first, second] = edge.faces.map(|u| u.face);
                (
                    [a, b].map(|corner| {
                        builder::line(&corner.contacts[&first], &corner.contacts[&second])
                    }),
                    None,
                )
            }
        };
        let wire = if edge.faces[0].forward {
            vec![
                rails[0].inverse(),
                seams[0].clone(),
                rails[1].clone(),
                seams[1].inverse(),
            ]
        } else {
            vec![
                rails[0].clone(),
                seams[1].clone(),
                rails[1].inverse(),
                seams[0].inverse(),
            ]
        };
        let wire = Wire::from(wire);
        let face = match surface {
            Some(surface) => Face::new(vec![wire], surface).ok()?,
            None => builder::try_attach_plane([wire]).ok()?,
        };
        certify_generated_face(&face).then_some(Self { rails, seams, face })
    }
}

fn orient(edge: &Edge, forward: bool) -> Edge {
    if forward {
        edge.clone()
    } else {
        edge.inverse()
    }
}

fn retained_endpoint(
    original: &SourceEdge,
    vertex: &Vertex,
    corners: &Map<VertexId, Corner>,
) -> Option<Vertex> {
    let Some(corner) = corners.get(&vertex.id()) else {
        return Some(vertex.clone());
    };
    let CornerKind::Cap { face, .. } = corner.kind else {
        return None;
    };
    if !original.uses.iter().any(|u| u.face == face) {
        return None;
    }
    let carrier = original.uses.iter().find(|u| u.face != face)?.face;
    corner.contacts.get(&carrier).cloned()
}

struct RetainedTrim {
    ends: [Vertex; 2],
    collapsed: bool,
}

/// Check the new endpoints against the finite original segment. A zero
/// limit is permitted only when both contact points actually coincide.
fn validate_trim(original: &Edge, ends: &[Vertex; 2]) -> Option<bool> {
    let old_a = original.absolute_front().point();
    let old_b = original.absolute_back().point();
    let axis = line_axis(original)?;
    let length = (old_b - old_a).dot(&axis);
    let [a, b] = ends.each_ref().map(Vertex::point);
    let ta = (a - old_a).dot(&axis);
    let tb = (b - old_a).dot(&axis);
    if [ta, tb]
        .into_iter()
        .any(|t| t < -KERNEL_EPSILON || t > length + KERNEL_EPSILON)
        || tb < ta - KERNEL_EPSILON
        || (old_a + axis * ta).distance(a) > KERNEL_EPSILON
        || (old_a + axis * tb).distance(b) > KERNEL_EPSILON
    {
        return None;
    }
    let collapsed = tb - ta <= KERNEL_EPSILON;
    (!collapsed || a.distance(b) <= KERNEL_EPSILON).then_some(collapsed)
}

fn retained_trims(
    incidence: &Incidence,
    selected: &Map<EdgeId, SelectedEdge>,
    corners: &Map<VertexId, Corner>,
) -> Option<Map<EdgeId, RetainedTrim>> {
    let mut trims = Map::default();
    for (&id, original) in &incidence.edges {
        if selected.contains_key(&id) {
            continue;
        }
        let old = [
            original.edge.absolute_front(),
            original.edge.absolute_back(),
        ];
        if !old.iter().any(|v| corners.contains_key(&v.id())) {
            continue;
        }
        let ends = [
            retained_endpoint(original, old[0], corners)?,
            retained_endpoint(original, old[1], corners)?,
        ];
        let collapsed = validate_trim(&original.edge, &ends)?;
        trims.insert(id, RetainedTrim { ends, collapsed });
    }
    Some(trims)
}

/// Contract connected components of proven zero-length retained segments.
/// Unlike spatial welding, only a previously shared topological edge can
/// identify its two proposed contact vertices. A transitive component must
/// still fit the original numerical coincidence bound.
fn collapsed_vertex_aliases(
    incidence: &Incidence,
    trims: &Map<EdgeId, RetainedTrim>,
) -> Option<Map<VertexId, Vertex>> {
    let mut adjacent = Map::<VertexId, Vec<Vertex>>::default();
    for trim in trims.values().filter(|trim| trim.collapsed) {
        let [a, b] = &trim.ends;
        adjacent.entry(a.id()).or_default().push(b.clone());
        adjacent.entry(b.id()).or_default().push(a.clone());
    }
    let mut aliases = Map::default();
    let mut seen = Set::default();
    for neighbors in adjacent.values() {
        let initial = neighbors.first()?;
        if !seen.insert(initial.id()) {
            continue;
        }
        let mut component = vec![initial.clone()];
        let mut cursor = 0;
        while cursor < component.len() {
            for neighbor in adjacent.get(&component[cursor].id())? {
                if seen.insert(neighbor.id()) {
                    component.push(neighbor.clone());
                }
            }
            cursor += 1;
        }
        // Preserve an unaffected source vertex when possible; otherwise
        // select a deterministic geometric representative.
        let canonical = component
            .iter()
            .min_by(|a, b| {
                let pa = a.point();
                let pb = b.point();
                incidence
                    .at_vertex
                    .contains_key(&b.id())
                    .cmp(&incidence.at_vertex.contains_key(&a.id()))
                    .then_with(|| pa.x.total_cmp(&pb.x))
                    .then_with(|| pa.y.total_cmp(&pb.y))
                    .then_with(|| pa.z.total_cmp(&pb.z))
                    .then_with(|| a.stable_id().cmp(&b.stable_id()))
            })?
            .clone();
        if component
            .iter()
            .any(|v| v.point().distance(canonical.point()) > KERNEL_EPSILON)
        {
            return None;
        }
        for vertex in component {
            aliases.insert(vertex.id(), canonical.clone());
        }
    }
    Some(aliases)
}

fn aliased_vertex(vertex: &Vertex, aliases: &Map<VertexId, Vertex>) -> Vertex {
    aliases
        .get(&vertex.id())
        .cloned()
        .unwrap_or_else(|| vertex.clone())
}

/// Keep a collapsed coedge as a vertex during face assembly, so its two
/// neighboring seam incidences still refer to their original junctions.
enum RebuiltEdge {
    Live(Edge),
    Collapsed(Vertex),
}

impl RebuiltEdge {
    fn front_id(&self) -> VertexId {
        match self {
            Self::Live(edge) => edge.front().id(),
            Self::Collapsed(vertex) => vertex.id(),
        }
    }

    fn back_id(&self) -> VertexId {
        match self {
            Self::Live(edge) => edge.back().id(),
            Self::Collapsed(vertex) => vertex.id(),
        }
    }

    fn oriented(&self, forward: bool) -> Self {
        match self {
            Self::Live(edge) => Self::Live(orient(edge, forward)),
            Self::Collapsed(vertex) => Self::Collapsed(vertex.clone()),
        }
    }
}

fn shortened_edges(
    incidence: &Incidence,
    selected: &Map<EdgeId, SelectedEdge>,
    trims: &Map<EdgeId, RetainedTrim>,
    aliases: &Map<VertexId, Vertex>,
) -> Option<Map<EdgeId, RebuiltEdge>> {
    let mut retained = Map::default();
    for (&id, original) in &incidence.edges {
        if selected.contains_key(&id) {
            continue;
        }
        let old = [
            original.edge.absolute_front(),
            original.edge.absolute_back(),
        ];
        let trim = trims.get(&id);
        let source = trim.map_or(old, |trim| [&trim.ends[0], &trim.ends[1]]);
        let ends = source.map(|v| aliased_vertex(v, aliases));
        let [a, b] = &ends;
        let replacement = if let Some(trim) = trim {
            let collapsed = validate_trim(&original.edge, &ends)?;
            if collapsed != trim.collapsed || (a.id() == b.id()) != collapsed {
                return None;
            }
            if collapsed {
                RebuiltEdge::Collapsed(a.clone())
            } else {
                RebuiltEdge::Live(builder::line(a, b))
            }
        } else if a.id() != old[0].id() || b.id() != old[1].id() {
            // Another edge already incident to an identified source vertex
            // keeps its geometry and follows the exact topology alias.
            if a.id() == b.id() {
                return None;
            }
            let mut edge = Edge::new(a, b, original.edge.curve()).ok()?;
            edge.set_stable_id(original.edge.stable_id());
            if !edge.is_geometric_consistent() {
                return None;
            }
            RebuiltEdge::Live(edge)
        } else {
            continue;
        };
        retained.insert(id, replacement);
    }
    Some(retained)
}

fn cap_seams(
    incidence: &Incidence,
    selected: &Map<EdgeId, SelectedEdge>,
    rounds: &Map<EdgeId, TreatedEdge>,
    corners: &Map<VertexId, Corner>,
) -> Option<Map<(VertexId, usize), Edge>> {
    let mut seams = Map::default();
    for (&id, corner) in corners {
        let CornerKind::Cap { face, .. } = corner.kind else {
            continue;
        };
        let edge = incidence.at_vertex[&id]
            .iter()
            .find(|id| selected.contains_key(id))?;
        let original = &incidence.edges[edge].edge;
        let index = usize::from(original.absolute_back().id() == id);
        let oriented = rounds[edge]
            .face
            .edge_iter()
            .find(|e| e.id() == rounds[edge].seams[index].id())?;
        seams.insert((id, face), oriented.inverse());
    }
    Some(seams)
}

struct ReplacementGeometry<'a> {
    selected: &'a Map<EdgeId, SelectedEdge>,
    rounds: &'a Map<EdgeId, TreatedEdge>,
    retained: &'a Map<EdgeId, RebuiltEdge>,
    cap_seams: &'a Map<(VertexId, usize), Edge>,
}

impl ReplacementGeometry<'_> {
    fn edge(&self, face: usize, edge: &Edge) -> Option<RebuiltEdge> {
        if let Some(round) = self.rounds.get(&edge.id()) {
            let slot = self.selected[&edge.id()]
                .faces
                .iter()
                .position(|u| u.face == face)?;
            return Some(RebuiltEdge::Live(orient(
                &round.rails[slot],
                edge.orientation(),
            )));
        }
        Some(self.retained.get(&edge.id()).map_or_else(
            || RebuiltEdge::Live(edge.clone()),
            |replacement| replacement.oriented(edge.orientation()),
        ))
    }

    fn face(&self, face_index: usize, face: &Face) -> Option<Face> {
        let mut changed = false;
        let mut boundaries = Vec::new();
        for wire in face.boundary_iters() {
            let source = wire.collect::<Vec<_>>();
            let edges = source
                .iter()
                .map(|edge| self.edge(face_index, edge))
                .collect::<Option<Vec<_>>>()?;
            let mut new_wire = Vec::new();
            for (i, edge) in edges.iter().enumerate() {
                match edge {
                    RebuiltEdge::Live(live) => {
                        changed |= live.id() != source[i].id();
                        new_wire.push(live.clone());
                    }
                    RebuiltEdge::Collapsed(_) => changed = true,
                }
                let next = &edges[(i + 1) % edges.len()];
                if edge.back_id() != next.front_id() {
                    let seam = self.cap_seams.get(&(source[i].back().id(), face_index))?;
                    if seam.front().id() != edge.back_id() || seam.back().id() != next.front_id() {
                        return None;
                    }
                    new_wire.push(seam.clone());
                    changed = true;
                }
            }
            if new_wire.is_empty() {
                return None;
            }
            boundaries.push(Wire::from(new_wire));
        }
        if !changed {
            return Some(face.clone());
        }
        let mut rebuilt = Face::new(boundaries, face.oriented_surface()).ok()?;
        rebuilt.set_stable_id(face.stable_id());
        rebuilt.is_geometric_consistent().then_some(rebuilt)
    }
}

fn spherical_face(corner: &Corner, original: &Vertex, mut seams: Vec<Edge>) -> Option<Face> {
    if seams.len() != 3 {
        return None;
    }
    let mut cycle = vec![seams.pop()?];
    while !seams.is_empty() {
        let last = cycle.last()?.back().id();
        let next = seams.iter().position(|edge| edge.front().id() == last)?;
        cycle.push(seams.swap_remove(next));
    }
    // All normals of the convex corner lie in the hemisphere pointing
    // towards the old sharp point. Place the spherical seam and poles in
    // the opposite hemisphere, away from the finite triangular trim.
    let x = (corner.center - original.point()).normalize();
    let seed = if x.x.abs() < x.y.abs() && x.x.abs() < x.z.abs() {
        Vector3::x()
    } else if x.y.abs() < x.z.abs() {
        Vector3::y()
    } else {
        Vector3::z()
    };
    let z = x.cross(&seed).normalize();
    let y = z.cross(&x);
    let transform = Matrix4::from_columns(&[
        x.push(0.0),
        y.push(0.0),
        z.push(0.0),
        corner.center.to_homogeneous(),
    ]);
    let surface = Surface::Sphere(Processor::with_transform(
        Sphere::new(Point3::origin(), corner.radius),
        transform,
    ));
    let face = Face::new(vec![Wire::from(cycle)], surface).ok()?;
    certify_generated_face(&face).then_some(face)
}

fn planar_network(
    shell: &Shell,
    edges: &[Edge],
    radii: &[f64],
    profile: PlanarTreatment,
) -> Option<Shell> {
    if edges.len() != radii.len() || edges.is_empty() {
        return None;
    }
    let incidence = Incidence::new(shell);
    let mut selected = Map::default();
    let mut vertices = Vec::<Vertex>::new();
    let mut seen = Set::default();
    for (edge, &radius) in edges.iter().zip(radii) {
        let id = edge.id();
        let selected_edge = SelectedEdge::new(shell, &incidence, id, radius, profile)?;
        if selected.insert(id, selected_edge).is_some() {
            return None;
        }
        for vertex in [edge.absolute_front(), edge.absolute_back()] {
            if seen.insert(vertex.id()) {
                vertices.push(vertex.clone());
            }
        }
    }
    let mut corners = vertices
        .iter()
        .map(|v| Some((v.id(), Corner::new(shell, &incidence, &selected, v)?)))
        .collect::<Option<Map<_, _>>>()?;
    let trims = retained_trims(&incidence, &selected, &corners)?;
    let aliases = collapsed_vertex_aliases(&incidence, &trims)?;
    for corner in corners.values_mut() {
        for contact in corner.contacts.values_mut() {
            *contact = aliased_vertex(contact, &aliases);
        }
    }
    // Both isolated straight rounds and trivalent rolling-ball networks use
    // exactly the same cap/rail assembly; no length heuristic is needed.
    let rounds = selected
        .values()
        .map(|edge| {
            Some((
                edge.id,
                TreatedEdge::new(edge, &incidence.edges[&edge.id].edge, &corners, profile)?,
            ))
        })
        .collect::<Option<Map<_, _>>>()?;
    let retained = shortened_edges(&incidence, &selected, &trims, &aliases)?;
    let cap_seams = cap_seams(&incidence, &selected, &rounds, &corners)?;
    let replacements = ReplacementGeometry {
        selected: &selected,
        rounds: &rounds,
        retained: &retained,
        cap_seams: &cap_seams,
    };
    let mut rebuilt = shell
        .iter()
        .enumerate()
        .map(|(i, face)| replacements.face(i, face))
        .collect::<Option<Vec<_>>>()?;
    for edge in edges {
        rebuilt.push(rounds[&edge.id()].face.clone());
    }
    for vertex in vertices {
        let corner = &corners[&vertex.id()];
        if !matches!(corner.kind, CornerKind::Sphere) {
            continue;
        }
        let mut seams = Vec::new();
        for id in &incidence.at_vertex[&vertex.id()] {
            let original = &incidence.edges[id].edge;
            let end = usize::from(original.absolute_back().id() == vertex.id());
            let round = &rounds[id];
            let edge = round
                .face
                .edge_iter()
                .find(|e| e.id() == round.seams[end].id())?;
            seams.push(edge.inverse());
        }
        rebuilt.push(spherical_face(corner, &vertex, seams)?);
    }
    let result = Shell::from(rebuilt);
    result.check_solid_boundary().ok()?;
    Some(result)
}

/// Construct an exact simultaneous convex planar three-way rolling-ball
/// network, including oblique planar terminal caps. The input is never
/// mutated and unsupported corner incidence is routed to the other solver.
pub(super) fn try_round_planar_network(
    shell: &Shell,
    edges: &[Edge],
    radii: &[f64],
) -> Option<Shell> {
    planar_network(shell, edges, radii, PlanarTreatment::Round)
        .or_else(|| tangent_rim::try_round_tangent_rim(shell, edges, radii))
}

pub(super) fn try_chamfer_planar_network(
    shell: &Shell,
    edges: &[Edge],
    setbacks: &[f64],
) -> Option<Shell> {
    planar_network(shell, edges, setbacks, PlanarTreatment::Chamfer)
}

pub(super) fn try_carrier_chamfer(
    shell: &Shell,
    edges: &[Edge],
    setbacks: &[f64],
    first_normal: Vector3,
    second_setback: f64,
) -> Option<Shell> {
    planar_network(
        shell,
        edges,
        setbacks,
        PlanarTreatment::CarrierChamfer {
            first_normal,
            second_setback,
        },
    )
}

pub(super) fn planar_chamfer_radii(
    shell: &Shell,
    edges: &[Edge],
    setbacks: &[f64],
) -> Option<Vec<f64>> {
    if edges.len() != setbacks.len() {
        return None;
    }
    let incidence = Incidence::new(shell);
    edges
        .iter()
        .zip(setbacks)
        .map(|(edge, &size)| {
            SelectedEdge::new(shell, &incidence, edge.id(), size, PlanarTreatment::Chamfer)
                .map(|selected| selected.radius)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Solid;

    fn oblique_prism() -> crate::Result<(Solid, Vec<Edge>)> {
        // An irregular trapezoidal profile in an oblique plane. Its top
        // has non-right-angle trihedral junctions and the lower round
        // termination curves are elliptical, not circular.
        let points = [[0.0, 0.0], [6.0, 0.0], [5.0, 4.0], [0.5, 4.0]]
            .map(|[x, y]| Point3::new(x, y, 0.2_f64.mul_add(x, 0.15 * y)));
        let vertices = points.map(Vertex::new);
        let wire = Wire::from(
            (0..4)
                .map(|i| builder::line(&vertices[i], &vertices[(i + 1) % 4]))
                .collect::<Vec<_>>(),
        );
        let cap: Face = builder::try_attach_plane([wire])?;
        let solid: Solid = builder::extrude(&cap, Vector3::new(0.0, 0.0, 4.0));
        let mut seen = Set::default();
        let edges = solid.boundaries()[0]
            .edge_iter()
            .filter(|edge| seen.insert(edge.id()))
            .filter(|edge| {
                [edge.front().point(), edge.back().point()]
                    .iter()
                    .any(|p| (p.z - 0.2_f64.mul_add(p.x, 0.15 * p.y) - 4.0).abs() <= KERNEL_EPSILON)
            })
            .collect();
        Ok((solid, edges))
    }

    #[test]
    fn oblique_round_contracts_exact_terminal_shoulder() -> crate::Result<()> {
        // The incline meets a vertical shoulder exactly at the 0.6 radius
        // tangency. At that limit one cap/side edge disappears while the
        // next retained edge and its original endpoint remain.
        let vertices = [
            [-4.3, 2.5],
            [-3.3, 2.5],
            [-3.3, 7.4],
            [-3.9, 8.0],
            [-4.3, 8.0],
        ]
        .map(|[y, z]| Vertex::new(Point3::new(0.0, y, z)));
        let wire = Wire::from(
            (0..5)
                .map(|i| builder::line(&vertices[i], &vertices[(i + 1) % 5]))
                .collect::<Vec<_>>(),
        );
        let cap: Face = builder::try_attach_plane([wire])?;
        let solid: Solid = builder::extrude(&cap, Vector3::new(6.65, 0.0, 0.0));
        let shell = &solid.boundaries()[0];
        let edge = shell
            .edge_iter()
            .find(|edge| {
                [edge.front().point(), edge.back().point()].iter().all(|p| {
                    (p.x - 6.65).abs() <= KERNEL_EPSILON && (p.z - 8.0).abs() <= KERNEL_EPSILON
                })
            })
            .expect("short sharp edge at the incline");
        for radius in [0.5, 0.6] {
            let rounded = planar_network(
                shell,
                std::slice::from_ref(&edge),
                &[radius],
                PlanarTreatment::Round,
            )
            .expect("analytic round fits up to and including the exact shoulder limit");
            rounded.check_solid_boundary()?;
            assert!(rounded.is_geometric_consistent());
            assert_eq!(rounded.len(), shell.len() + 1);
            assert!(rounded.iter().all(certify_generated_face));
        }
        assert!(planar_network(shell, &[edge], &[0.7], PlanarTreatment::Round).is_none());
        Ok(())
    }

    #[test]
    fn oblique_chamfers_use_physical_setback_and_compose() -> crate::Result<()> {
        let (solid, edges) = oblique_prism()?;
        let shell = &solid.boundaries()[0];
        let vertical_at = |edge: &Edge, x: f64, y: f64| {
            [edge.front().point(), edge.back().point()]
                .iter()
                .all(|p| (p.x - x).abs() < KERNEL_EPSILON && (p.y - y).abs() < KERNEL_EPSILON)
        };
        let edge = edges
            .iter()
            .find(|e| vertical_at(e, 5.0, 4.0))
            .expect("first corner");
        let setback = 0.25;
        let first = try_chamfer_planar_network(shell, std::slice::from_ref(edge), &[setback])
            .expect("native oblique chamfer");
        assert_eq!(first.len(), 7);
        first.check_solid_boundary()?;
        assert!(first.is_geometric_consistent());
        let bevel = first.last().expect("bevel face");
        assert!(matches!(bevel.surface(), Surface::Plane(_)));
        for vertex in bevel.vertex_iter() {
            let p = vertex.point();
            assert!(((p.x - 5.0).hypot(p.y - 4.0) - setback).abs() < KERNEL_EPSILON);
        }
        let second_edge = first
            .edge_iter()
            .find(|e| vertical_at(e, 0.0, 0.0))
            .expect("second corner");
        let second = try_chamfer_planar_network(&first, &[second_edge], &[setback])
            .expect("second native chamfer");
        assert_eq!(second.len(), 8);
        second.check_solid_boundary()?;
        assert!(second.is_geometric_consistent());
        assert!(
            second
                .edge_iter()
                .all(|e| matches!(e.curve(), Curve::Line(_)))
        );
        Ok(())
    }

    #[test]
    fn oblique_round_network_is_closed_and_order_independent() -> crate::Result<()> {
        let (solid, mut edges) = oblique_prism()?;
        let shell = &solid.boundaries()[0];
        assert_eq!(edges.len(), 8);
        let radii = [0.25; 8];
        let first =
            try_round_planar_network(shell, &edges, &radii).expect("analytic oblique network");
        assert_eq!(first.len(), 18);
        assert!(first.is_geometric_consistent());
        first.check_solid_boundary()?;
        edges.reverse();
        let second =
            try_round_planar_network(shell, &edges, &radii).expect("reverse selected order");
        second.check_solid_boundary()?;
        assert_eq!(second.len(), 18);
        let vertices = |shell: &Shell| {
            let mut seen = Set::default();
            shell
                .vertex_iter()
                .filter(|v| seen.insert(v.id()))
                .map(|v| v.point())
                .collect::<Vec<_>>()
        };
        let first_points = vertices(&first);
        let second_points = vertices(&second);
        assert_eq!(first_points.len(), second_points.len());
        assert!(
            first_points
                .iter()
                .all(|point| second_points.iter().any(|other| point.near(other)))
        );
        assert!(try_round_planar_network(shell, &edges, &[4.0; 8]).is_none());
        let mut incompatible = radii;
        incompatible[0] = 0.4;
        assert!(try_round_planar_network(shell, &edges, &incompatible).is_none());
        assert_eq!(shell.len(), 6);
        assert!(shell.is_geometric_consistent());
        Ok(())
    }
}
