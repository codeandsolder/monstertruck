//! Planar Boolean overlay retaining original analytic edges.
//!
//! This operates on noded exact-edge arrangements. In particular, repeated
//! boundary walks and overlapping selected cells are unified by nonzero
//! material occupancy, never by parity or discarding positional face indices.

use super::divide_face::create_parameter_boundary_on_surface;
use super::faces_classification::BooleanRegion;
use super::{AltCurve, AltRegion, ShapeOpsCurve, ShapeOpsSurface};
use monstertruck_geometry::prelude::*;
use monstertruck_topology::{Edge, EdgeId, VertexId, Wire};
use rustc_hash::FxHashMap as HashMap;
use std::collections::BTreeSet;

struct ProjectedEdge<C> {
    edge: Edge<Point3, C>,
    uv: Vec<Point2>,
}

fn signed_area(poly: &[Point2]) -> f64 {
    poly.iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
        .map(|(a, b)| a.x.mul_add(b.y, -b.x * a.y))
        .sum::<f64>()
        * 0.5
}

fn winding(poly: &[Point2], point: Point2) -> isize {
    poly.iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
        .fold(0_isize, |w, (a, b)| {
            let side = (b.x - a.x).mul_add(point.y - a.y, -(b.y - a.y) * (point.x - a.x));
            if a.y <= point.y && b.y > point.y && side > 0.0 {
                w + 1
            } else if b.y <= point.y && a.y > point.y && side < 0.0 {
                w - 1
            } else {
                w
            }
        })
}

/// Detect edges that would have to be cut before an exact-edge arrangement.
/// Proper segment crossings and a true endpoint-in-interior junction cannot
/// be represented by tracing the original unsplit edge identities.
fn needs_exact_noding<C>(edges: &[&ProjectedEdge<C>], span: f64) -> bool {
    let eps = span * 1.0e-10;
    let cross = |a: Point2, b: Point2, p: Point2| {
        let d = b - a;
        d.x.mul_add(p.y - a.y, -(d.y * (p.x - a.x)))
    };
    let on_interior = |p: Point2, a: Point2, b: Point2| {
        let d = b - a;
        let length2 = d.magnitude2();
        if length2 <= eps * eps {
            return false;
        }
        let t = (p - a).dot(d) / length2;
        t > 1.0e-6 && t < 1.0 - 1.0e-6 && (p - (a + d * t)).magnitude2() < eps * eps
    };
    for (i, lhs) in edges.iter().enumerate() {
        for rhs in edges.iter().skip(i + 1) {
            let a = &lhs.uv;
            let b = &rhs.uv;
            for &(point, ref other) in &[
                (a[0], b),
                (a[a.len() - 1], b),
                (b[0], a),
                (b[b.len() - 1], a),
            ] {
                if other.windows(2).any(|w| on_interior(point, w[0], w[1])) {
                    return true;
                }
            }
            for line_a in a.windows(2) {
                let (p, q) = (line_a[0], line_a[1]);
                for line_b in b.windows(2) {
                    let (r, s) = (line_b[0], line_b[1]);
                    let (c0, c1) = (cross(p, q, r), cross(p, q, s));
                    let (c2, c3) = (cross(r, s, p), cross(r, s, q));
                    let e0 = (q - p).magnitude() * eps;
                    let e1 = (s - r).magnitude() * eps;
                    if (c0 > e0 && c1 < -e0 || c0 < -e0 && c1 > e0)
                        && (c2 > e1 && c3 < -e1 || c2 < -e1 && c3 > e1)
                    {
                        return true;
                    }
                }
            }
        }
    }
    false
}

/// Planar union of regions already noded at their shared boundary vertices.
/// Aborts if original edge topology is insufficient to trace a closed union;
/// never replaces analytic edges with an approximate tessellation.
fn overlay_one_source<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    source: &[AltRegion<C, S>],
    tol: f64,
) -> Option<Vec<AltRegion<C, S>>> {
    let reference = source.first()?;
    let surface = reference.surface();
    if !matches!(
        surface.try_into_analytic_surface_kind(),
        Some(AnalyticSurfaceKind::Plane(_))
    ) {
        return None;
    }
    let mut samples = HashMap::default();
    let mut polygons = Vec::<Vec<Vec<Point2>>>::new();
    let mut edges = HashMap::<EdgeId<AltCurve<C, S>>, ProjectedEdge<AltCurve<C, S>>>::default();
    let mut span_min = Point2::new(f64::INFINITY, f64::INFINITY);
    let mut span_max = Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
    for region in source {
        let mut region_loops = Vec::new();
        for wire in region.absolute_boundaries() {
            let poly =
                create_parameter_boundary_on_surface(&surface, wire, None, &mut samples, tol)?;
            if poly.len() < 3 {
                return None;
            }
            for &p in &poly.0 {
                span_min.x = span_min.x.min(p.x);
                span_min.y = span_min.y.min(p.y);
                span_max.x = span_max.x.max(p.x);
                span_max.y = span_max.y.max(p.y);
            }
            region_loops.push(poly.0);
            for edge in wire {
                if edges.contains_key(&edge.id()) {
                    continue;
                }
                let absolute = edge.absolute_clone();
                let one = Wire::from(vec![absolute.clone()]);
                let projected =
                    create_parameter_boundary_on_surface(&surface, &one, None, &mut samples, tol)?;
                if projected.len() < 2 {
                    return None;
                }
                edges.insert(
                    edge.id(),
                    ProjectedEdge {
                        edge: absolute,
                        uv: projected.0,
                    },
                );
            }
        }
        polygons.push(region_loops);
    }
    let span = (span_max.x - span_min.x).max(span_max.y - span_min.y);
    if !span.is_finite() || span <= 0.0 {
        return None;
    }
    if needs_exact_noding(&edges.values().collect::<Vec<_>>(), span) {
        return None;
    }
    let covered = |point: Point2| {
        polygons.iter().any(|region| {
            region
                .iter()
                .map(|poly| winding(poly, point))
                .sum::<isize>()
                != 0
        })
    };

    // A boundary survives exactly when its sides differ in union occupancy.
    let mut boundary = Vec::<(Edge<Point3, AltCurve<C, S>>, Vec<Point2>)>::new();
    for item in edges.values() {
        let middle = item.uv.len() / 2;
        let (a, b) = (item.uv[middle - 1], item.uv[middle]);
        let tangent = b - a;
        let length = tangent.magnitude();
        if length <= 0.0 {
            return None;
        }
        let left = Vector2::new(-tangent.y, tangent.x) / length;
        let center = a.midpoint(b);
        let mut membership = None;
        for ratio in [1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5] {
            let offset = ratio * span;
            let (l, r) = (
                covered(center + left * offset),
                covered(center - left * offset),
            );
            if l != r {
                membership = Some((l, r));
                break;
            }
        }
        match membership {
            Some((true, false)) => boundary.push((item.edge.clone(), item.uv.clone())),
            Some((false, true)) => {
                boundary.push((item.edge.inverse(), item.uv.iter().rev().copied().collect()))
            }
            _ => {}
        }
    }
    if boundary.is_empty() {
        // Do not erase a filled source region merely because the probes were
        // unable to resolve a tiny or degenerate boundary.
        return None;
    }
    // Stable ordering matters at multi-edge junctions and also makes the
    // reconstructed face sequence independent of hash-table iteration.
    boundary.sort_by(|a, b| {
        let p = &a.1;
        let q = &b.1;
        p[0].x
            .total_cmp(&q[0].x)
            .then_with(|| p[0].y.total_cmp(&q[0].y))
            .then_with(|| p[p.len() - 1].x.total_cmp(&q[q.len() - 1].x))
            .then_with(|| p[p.len() - 1].y.total_cmp(&q[q.len() - 1].y))
    });

    let mut outgoing = HashMap::<VertexId<Point3>, Vec<usize>>::default();
    for (index, (edge, _)) in boundary.iter().enumerate() {
        outgoing.entry(edge.front().id()).or_default().push(index);
    }
    let mut unseen = (0..boundary.len()).collect::<BTreeSet<_>>();
    let mut loops = Vec::<(Wire<Point3, AltCurve<C, S>>, Vec<Point2>)>::new();
    while let Some(&first) = unseen.first() {
        let mut index = first;
        let start = boundary[first].0.front().id();
        let (mut chain, mut points) = (Vec::new(), Vec::new());
        loop {
            if !unseen.remove(&index) {
                return None;
            }
            let (edge, uv) = &boundary[index];
            chain.push(edge.clone());
            points.extend(uv.iter().take(uv.len() - 1).copied());
            let back = edge.back().id();
            if back == start {
                break;
            }
            let candidates = outgoing.get(&back)?;
            let last = uv[uv.len() - 1] - uv[uv.len() - 2];
            let reverse_angle = (-last.y).atan2(-last.x);
            let next = candidates
                .iter()
                .filter(|&&candidate| unseen.contains(&candidate))
                .min_by(|&&lhs, &&rhs| {
                    let angle = |candidate: usize| {
                        let u = &boundary[candidate].1;
                        let d = u[1] - u[0];
                        (reverse_angle - d.y.atan2(d.x)).rem_euclid(std::f64::consts::TAU)
                    };
                    angle(lhs).total_cmp(&angle(rhs))
                })?;
            index = *next;
        }
        loops.push((Wire::from(chain), points));
    }

    // Check that the traced boundary has the same material occupancy as the
    // original region union on both sides of every source edge. Detect graph
    // omissions and un-noded arrangements instead of emitting a plausible
    // but incorrect shell.
    let traced_coverage = |point: Point2| {
        loops
            .iter()
            .map(|(_, uv)| winding(uv, point))
            .sum::<isize>()
            != 0
    };
    for item in edges.values() {
        let middle = item.uv.len() / 2;
        let a = item.uv[middle - 1];
        let b = item.uv[middle];
        let segment = b - a;
        let length = segment.magnitude();
        if length <= 0.0 {
            return None;
        }
        let left = Vector2::new(-segment.y, segment.x) / length;
        let center = a.midpoint(b);
        for ratio in [1.0e-7, 1.0e-6] {
            let offset = span * ratio;
            let l = center + left * offset;
            let r = center - left * offset;
            if covered(l) != traced_coverage(l) || covered(r) != traced_coverage(r) {
                return None;
            }
        }
    }

    let mut outer = Vec::<(Vec<Wire<Point3, AltCurve<C, S>>>, Vec<Point2>)>::new();
    let mut holes = Vec::new();
    for (wire, poly) in loops {
        if signed_area(&poly) > 0.0 {
            outer.push((vec![wire], poly));
        } else {
            holes.push((wire, poly));
        }
    }
    for (hole, uv) in holes {
        let point = *uv.first()?;
        let parent = outer
            .iter()
            .enumerate()
            .filter(|(_, (_, poly))| winding(poly, point) != 0)
            .min_by(|(_, a), (_, b)| signed_area(&a.1).total_cmp(&signed_area(&b.1)))
            .map(|(i, _)| i)?;
        outer[parent].0.push(hole);
    }
    Some(
        outer
            .into_iter()
            .enumerate()
            .map(|(key, (wires, _))| {
                BooleanRegion::new(
                    wires,
                    surface.clone(),
                    reference.orientation(),
                    key,
                    reference.source_face(),
                )
            })
            .collect(),
    )
}

/// Unify coincident selected cells from each original planar source face.
pub(super) fn overlay_source_regions<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    regions: Vec<AltRegion<C, S>>,
    tol: f64,
) -> Vec<AltRegion<C, S>> {
    let mut groups = Vec::<Vec<AltRegion<C, S>>>::new();
    let mut group_indices = HashMap::<(usize, bool), usize>::default();
    for region in regions {
        let key = (region.source_face(), region.orientation());
        if let Some(&index) = group_indices.get(&key) {
            groups[index].push(region);
        } else {
            group_indices.insert(key, groups.len());
            groups.push(vec![region]);
        }
    }
    groups
        .into_iter()
        .flat_map(|group| {
            if group.len() > 1 {
                if let Some(union) = overlay_one_source(&group, tol) {
                    return union;
                }
            }
            group
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::super::Alternative;
    use super::*;
    use monstertruck_geometry::prelude::{Line, Plane};
    use monstertruck_modeling::{Curve, Surface};
    use monstertruck_topology::Vertex;

    type TestRegion = AltRegion<Curve, Surface>;
    type TestEdge = Edge<Point3, AltCurve<Curve, Surface>>;

    struct PlanarFixture {
        vertices: Vec<Vertex<Point3>>,
        edges: HashMap<(usize, usize), TestEdge>,
    }

    impl PlanarFixture {
        fn new(points: &[[f64; 2]]) -> Self {
            Self {
                vertices: points
                    .iter()
                    .map(|p| Vertex::new(Point3::new(p[0], p[1], 0.0)))
                    .collect(),
                edges: HashMap::default(),
            }
        }

        fn edge(&mut self, a: usize, b: usize) -> TestEdge {
            if let Some(edge) = self.edges.get(&(a, b)) {
                return edge.clone();
            }
            if let Some(edge) = self.edges.get(&(b, a)) {
                return edge.inverse();
            }
            let (start, end) = (&self.vertices[a], &self.vertices[b]);
            let curve = Alternative::FirstType(Curve::from(Line(start.point(), end.point())));
            let edge = Edge::new_unchecked(start, end, curve);
            self.edges.insert((a, b), edge.clone());
            edge
        }

        fn region(&mut self, path: &[usize]) -> TestRegion {
            let wire = (0..path.len() - 1)
                .map(|i| self.edge(path[i], path[i + 1]))
                .collect::<Vec<_>>();
            BooleanRegion::new(
                vec![Wire::from(wire)],
                Surface::Plane(Plane::xy()),
                true,
                0,
                0,
            )
        }
    }

    #[test]
    fn overlapping_same_source_regions_are_a_nonzero_union() {
        // An outer six-edge square and its small corner triangle are present
        // as both a repeated subwalk and a separately selected positive cell.
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
            [0.0, 1.0],
        ]);
        let repeated = fixture.region(&[0, 1, 5, 0, 1, 2, 3, 4, 5, 0]);
        let triangle = fixture.region(&[0, 1, 5, 0]);
        let union = overlay_source_regions(vec![repeated, triangle], 1.0e-6);
        assert_eq!(union.len(), 1);
        assert_eq!(union[0].absolute_boundaries().len(), 1);
        let wire = &union[0].absolute_boundaries()[0];
        assert_eq!(wire.len(), 6);
        assert!(!wire.iter().any(|edge| edge.id() == fixture.edge(1, 5).id()));
        let mut polys = HashMap::default();
        let projected = create_parameter_boundary_on_surface(
            &union[0].surface(),
            wire,
            None,
            &mut polys,
            1.0e-6,
        )
        .unwrap();
        assert!((signed_area(&projected.0) - 16.0).abs() < 1.0e-8);
    }

    #[test]
    fn adjacent_coplanar_cells_drop_shared_internal_edge() {
        let mut fixture = PlanarFixture::new(&[
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ]);
        let left = fixture.region(&[0, 1, 2, 3, 0]);
        let right = fixture.region(&[1, 4, 5, 2, 1]);
        let shared = fixture.edge(1, 2).id();
        let union = overlay_source_regions(vec![left, right], 1.0e-6);
        assert_eq!(union.len(), 1);
        let wire = &union[0].absolute_boundaries()[0];
        assert_eq!(wire.len(), 6);
        assert!(!wire.iter().any(|edge| edge.id() == shared));
    }

    #[test]
    fn duplicate_positive_regions_do_not_cancel_by_parity() {
        let mut fixture = PlanarFixture::new(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]);
        let a = fixture.region(&[0, 1, 2, 3, 0]);
        let b = fixture.region(&[0, 1, 2, 3, 0]);
        let union = overlay_source_regions(vec![a, b], 1.0e-6);
        assert_eq!(union.len(), 1);
        assert_eq!(union[0].absolute_boundaries()[0].len(), 4);
    }
}
