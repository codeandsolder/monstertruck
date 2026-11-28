//! Bulk planar normalization of certified collapsed rational fans.
//!
//! The source-backed recovery's ruled fan faces are recognizable by their
//! exact boundary control net and a collapsed linear surface direction.
//! Recheck each fan on a common plane, prove the oriented projected
//! coverage, cancel only opposite shared internal edges, and certify the
//! resulting outer/hole domain. Retain the original external edges.
//! All geometric movement is explicitly bounded and reported.

use super::{
    planar_domain::{certify_projected_planar_domain, certify_projected_radial_fan},
    sample_fit::normal_basis,
};
use crate::{
    Curve, Edge, EdgeId, Face, KERNEL_EPSILON, NurbsSurface, Plane, Point3, Shell, Surface,
    Vector3, Vector4, VertexId, Wire,
};
use monstertruck_geometry::prelude::TryIntoHomogeneousBsplineCurve;
use monstertruck_topology::shell::ShellCondition;
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::BTreeSet;

/// A bulk planar simplification and its face-level provenance.
#[derive(Debug)]
pub struct PlanarFanSimplification {
    /// Validated shell with the surviving native exterior edges unchanged.
    pub shell: Shell,
    /// New native face index for every original face.
    pub original_to_simplified: Vec<usize>,
    /// Number of accepted connected multi-face unions.
    pub merged_groups: usize,
    /// Number of connected candidate groups that could not be certified.
    pub rejected_groups: usize,
    /// Additional two-sided deviation bound for the endpoint-aligned
    /// planar reference domain, including plane and seam uncertainty.
    pub deviation_upper_bound: f64,
    /// Maximum endpoint-alignment uncertainty within accepted unions.
    pub seam_deviation_upper_bound: f64,
}

/// Explicit policy for a certified planar-fan simplification.
#[derive(Clone, Copy, Debug)]
pub struct PlanarFanOptions<'a> {
    /// Maximum additional geometric deviation; also limited by the kernel tolerance.
    pub maximum_additional_deviation: f64,
    /// Optional face ownership. Faces from different owners cannot be merged.
    pub owner_by_face: Option<&'a [usize]>,
    /// Native edges that must retain their topological identity.
    pub protected_edges: &'a [EdgeId],
}

impl PlanarFanOptions<'_> {
    /// An unpartitioned shell with no extra protected interior edges.
    #[must_use = "the constructed options carry the requested approximation budget"]
    pub const fn new(maximum_additional_deviation: f64) -> Self {
        Self {
            maximum_additional_deviation,
            owner_by_face: None,
            protected_edges: &[],
        }
    }
}

fn point(p: Vector4) -> Point3 { Point3::new(p.x / p.w, p.y / p.w, p.z / p.w) }

fn exact_radial(edge: &Edge) -> bool {
    matches!(edge.oriented_curve(), Curve::Line(line) if line.0==edge.front().point() && line.1==edge.back().point())
}

struct Fan {
    base: Edge,
    apex: Point3,
    points: Vec<Point3>,
    plane: Plane,
    margin: f64,
    collapse: f64,
}

impl Fan {
    fn from_face(face: &Face, limit: f64) -> Option<Self> {
        let Surface::NurbsSurface(mut surface) = face.surface() else {
            return None;
        };
        let wires = face.boundaries();
        let [wire] = wires.as_slice() else {
            return None;
        };
        if wire.len() != 3 || !wire.is_closed() || !wire.is_simple() {
            return None;
        }
        // Either parametric direction may be the collapsed linear one.
        for _ in 0..2 {
            if surface.degrees().1 == 1
                && surface.knot_vector_v().is_clamped(1)
                && surface.control_points().iter().all(|r| r.len() == 2)
            {
                for index in 0..3 {
                    let base = &wire[index];
                    let radial = &wire[(index + 1) % 3];
                    let back = &wire[(index + 2) % 3];
                    if !exact_radial(radial) || !exact_radial(back) {
                        continue;
                    }
                    let apex = radial.back().point();
                    if let Some(candidate) = Self::from_net(&surface, base, apex, limit) {
                        return Some(candidate);
                    }
                }
            }
            surface.swap_axes();
        }
        None
    }

    fn from_net(
        surface: &NurbsSurface<Vector4>,
        base: &Edge,
        apex: Point3,
        limit: f64,
    ) -> Option<Self> {
        let controls = surface.control_points();
        let [front, back] = [base.front().point(), base.back().point()];
        let orientation = (back - front).cross(&(apex - front));
        // The native spatial edge must be the same rational control net
        // as one full surface side, up to orientation. Carrier proximity
        // alone cannot establish that equality.
        let matches = |row: usize| {
            [base.clone(), base.inverse()].iter().any(|edge| {
                edge.oriented_curve()
                    .try_into_homogeneous_bspline_curve()
                    .is_some_and(|curve| {
                        curve.knot_vector() == surface.knot_vector_u()
                            && curve.control_points().len() == controls.len()
                            && curve
                                .control_points()
                                .iter()
                                .zip(controls)
                                .all(|(a, b)| *a == b[row])
                    })
            })
        };
        let row = (0..2).find(|&row| matches(row))?;
        if controls
            .iter()
            .flatten()
            .any(|p| p.w <= 0.0 || ![p.x, p.y, p.z, p.w].into_iter().all(f64::is_finite))
        {
            return None;
        }
        let points = controls
            .iter()
            .flatten()
            .copied()
            .map(point)
            .collect::<Vec<_>>();
        let world_scale = points
            .iter()
            .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
            .fold(1.0_f64, f64::max);
        let (min_weight, max_weight) = controls
            .iter()
            .flatten()
            .map(|p| p.w)
            .fold((f64::INFINITY, 0.0_f64), |(lo, hi), w| {
                (lo.min(w), hi.max(w))
            });
        let degree = crate::usize_to_f64(surface.degrees().0 + 2);
        let margin =
            256.0 * f64::EPSILON * world_scale * (max_weight / min_weight) * degree * degree;
        let collapse = controls
            .iter()
            .map(|r| (point(r[1 - row]) - apex).norm())
            .fold(0.0_f64, f64::max)
            + margin;
        if !collapse.is_finite() || collapse > limit {
            return None;
        }
        // Use the best-spanned normal of the whole positive-weight
        // net instead of amplifying a sliver's short chord cross product.
        let origin = *points.first()?;
        let axis = points
            .iter()
            .map(|&p| p - origin)
            .max_by(|a, b| a.norm_squared().total_cmp(&b.norm_squared()))?;
        let normal = points
            .iter()
            .map(|&p| axis.cross(&(p - origin)))
            .max_by(|a, b| a.norm_squared().total_cmp(&b.norm_squared()))?;
        let length = normal.norm();
        let sign = normal.dot(&orientation);
        if !length.is_finite() || length <= 0.0 || !sign.is_finite() || sign.abs() == 0.0 {
            return None;
        }
        let normal = normal * (sign.signum() / length);
        let (u, v) = normal_basis(normal);
        let plane = Plane::new(origin, origin + u, origin + v);
        let normal = plane.normal();
        if !normal.norm_squared().is_finite() || normal.dot(&orientation) <= 0.0 {
            return None;
        }
        let candidate = Self {
            base: base.clone(),
            apex,
            points,
            plane,
            margin,
            collapse,
        };
        (candidate.plane_bound(origin, normal) <= limit).then_some(candidate)
    }

    fn plane_bound(&self, origin: Point3, normal: Vector3) -> f64 {
        // The full rational control hull and native corners lie in this
        // slab; positive weights bound every surface point, not samples.
        let slab = self
            .points
            .iter()
            .copied()
            .chain([
                self.apex,
                self.base.front().point(),
                self.base.back().point(),
            ])
            .map(|p| normal.dot(&(p - origin)).abs())
            .fold(0.0_f64, f64::max)
            + self.margin;
        (128.0 * f64::EPSILON).mul_add(origin.x.abs().max(origin.y.abs()).max(origin.z.abs()), slab)
    }

    fn certify(&self, origin: Point3, normal: Vector3, limit: f64) -> Option<Bounds> {
        if self.plane.normal().dot(&normal) <= 0.0 {
            return None;
        }
        let plane = self.plane_bound(origin, normal);
        // The ideal radial surface differs from the original collapsed
        // row by at most collapse, with positive rational interpolation.
        let remaining = limit - plane - self.collapse;
        if !remaining.is_finite() || remaining <= 0.0 {
            return None;
        }
        let seam = certify_projected_radial_fan(&self.base, self.apex, origin, normal, remaining)?;
        let bound = Bounds {
            plane_and_collapse: plane + self.collapse,
            seam,
        };
        (bound.deviation() <= limit).then_some(bound)
    }
}

#[derive(Clone, Copy, Default)]
struct Bounds {
    plane_and_collapse: f64,
    seam: f64,
}

impl Bounds {
    fn deviation(self) -> f64 { self.plane_and_collapse + self.seam }

    const fn include(&mut self, other: Self) {
        self.plane_and_collapse = self.plane_and_collapse.max(other.plane_and_collapse);
        self.seam = self.seam.max(other.seam);
    }
}

type Incidence = FxHashMap<EdgeId, Vec<(usize, bool)>>;

struct UnionContext<'a> {
    shell: &'a Shell,
    incidence: Incidence,
    protected: FxHashSet<EdgeId>,
    adjacent: Vec<BTreeSet<usize>>,
    fans: Vec<Option<Fan>>,
    limit: f64,
}

impl<'a> UnionContext<'a> {
    fn new(shell: &'a Shell, options: PlanarFanOptions<'_>, limit: f64) -> Self {
        let protected = options
            .protected_edges
            .iter()
            .copied()
            .collect::<FxHashSet<_>>();
        let mut incidence = Incidence::default();
        for (index, face) in shell.iter().enumerate() {
            for edge in face.boundaries().iter().flatten() {
                incidence
                    .entry(edge.id())
                    .or_default()
                    .push((index, edge.orientation()));
            }
        }
        let mut adjacent = vec![BTreeSet::new(); shell.len()];
        for (id, uses) in &incidence {
            if let [(a, forward), (b, backward)] = uses.as_slice()
                && a != b
                && forward != backward
                && !protected.contains(id)
                && options
                    .owner_by_face
                    .is_none_or(|owner| owner[*a] == owner[*b])
            {
                adjacent[*a].insert(*b);
                adjacent[*b].insert(*a);
            }
        }
        let fans = shell
            .iter()
            .map(|face| Fan::from_face(face, limit))
            .collect();
        Self {
            shell,
            incidence,
            protected,
            adjacent,
            fans,
            limit,
        }
    }

    fn group(&self, seed: usize, unseen: &mut BTreeSet<usize>) -> Option<(Vec<usize>, Bounds)> {
        let first = self.fans[seed].as_ref()?;
        let mut bounds = first.certify(first.plane.origin(), first.plane.normal(), self.limit)?;
        let mut group = vec![seed];
        let mut index = 0;
        while index < group.len() {
            let current = group[index];
            for &other in &self.adjacent[current] {
                if !unseen.contains(&other) {
                    continue;
                }
                if let Some(candidate) = self.fans[other].as_ref()
                    && let Some(bound) =
                        candidate.certify(first.plane.origin(), first.plane.normal(), self.limit)
                {
                    unseen.remove(&other);
                    group.push(other);
                    bounds.include(bound);
                }
            }
            index += 1;
        }
        group.sort_unstable();
        Some((group, bounds))
    }

    fn wires(&self, group: &[usize]) -> Option<Vec<Wire>> {
        let mut boundary = Vec::new();
        for &index in group {
            for edge in self.shell[index].boundaries().into_iter().flatten() {
                let uses = self.incidence.get(&edge.id())?;
                match uses.as_slice() {
                    [(a, forward), (b, backward)] if a != b && forward != backward => {
                        if group.binary_search(a).is_ok() && group.binary_search(b).is_ok() {
                            if self.protected.contains(&edge.id()) {
                                return None;
                            }
                            continue;
                        }
                    }
                    [(_, _)] => {}
                    _ => return None,
                }
                boundary.push(edge);
            }
        }
        let mut outgoing = FxHashMap::<VertexId, usize>::default();
        let mut incoming = FxHashSet::<VertexId>::default();
        for (index, edge) in boundary.iter().enumerate() {
            if outgoing.insert(edge.front().id(), index).is_some()
                || !incoming.insert(edge.back().id())
            {
                return None;
            }
        }
        if outgoing.is_empty()
            || outgoing.len() != incoming.len()
            || !outgoing.keys().all(|id| incoming.contains(id))
        {
            return None;
        }
        let mut wires = Vec::new();
        // Use source order, not pointer/hash ordering, for reproducibility.
        for edge in &boundary {
            let start = edge.front().id();
            if !outgoing.contains_key(&start) {
                continue;
            }
            let mut next = start;
            let mut edges = Vec::new();
            while let Some(index) = outgoing.remove(&next) {
                let edge = &boundary[index];
                next = edge.back().id();
                edges.push(edge.clone());
            }
            if next != start {
                return None;
            }
            wires.push(Wire::from(edges));
        }
        Some(wires)
    }

    fn merge(&self, group: &[usize], mut bounds: Bounds) -> Option<(Face, Bounds)> {
        let first = self.fans[*group.first()?].as_ref()?;
        // The group was grown in the seed frame, and its seed is its
        // smallest still-unseen face. Recheck the surviving curved loops
        // positively in exactly that frame. A winding sum of the
        // positively oriented constituent domains cancels internal edges,
        // leaving the simple outer loop minus disjoint holes.
        let mut wires = self.wires(group)?;
        let cert = certify_projected_planar_domain(
            &wires,
            first.plane.origin(),
            first.plane.normal(),
            self.limit,
        )?;
        // The final loop proof and constituent proofs may subdivide at
        // different depths. They share the same aligned edge definition;
        // charge any larger arithmetic/seam bound to the result.
        bounds.seam = bounds.seam.max(cert.seam_displacement_upper_bound);
        if bounds.deviation() > self.limit {
            return None;
        }
        wires.swap(0, cert.outer_wire);
        let face = Face::new(wires, Surface::Plane(first.plane)).ok()?;
        face.is_geometric_consistent().then_some((face, bounds))
    }
}

/// Replace connected coplanar rational fan faces by certified simple
/// planar domains, without modifying their surviving native edges.
///
/// The caller supplies an additional geometric budget and must charge
/// the reported bound to any original approximation/vertex budgets.
/// Optional per-face ownership prevents cancellation across different
/// reconstruction regions. Non-fan faces and unproved groups are retained.
/// Positive certificates concern the endpoint-aligned domain, with
/// original boundary mismatch bounded separately. This does not prove
/// global solid embedding or material classification.
///
/// Returns None for malformed options or shell topology, or if the final
/// transaction cannot preserve its original topology.
#[must_use = "the returned shell contains the simplified faces and their deviation bound"]
pub fn simplify_planar_fan_shell(
    shell: &Shell,
    options: PlanarFanOptions<'_>,
) -> Option<PlanarFanSimplification> {
    let condition = shell.shell_condition();
    if !options.maximum_additional_deviation.is_finite()
        || options.maximum_additional_deviation < 0.0
        || options
            .owner_by_face
            .is_some_and(|owners| owners.len() != shell.len())
        || !matches!(condition, ShellCondition::Closed | ShellCondition::Oriented)
        || !shell.singular_vertices().is_empty()
        || !shell.is_geometric_consistent()
    {
        return None;
    }
    let limit = options.maximum_additional_deviation.min(KERNEL_EPSILON);
    let context = UnionContext::new(shell, options, limit);
    let mut unseen = context
        .fans
        .iter()
        .enumerate()
        .filter_map(|(i, fan)| fan.as_ref().map(|_| i))
        .collect::<BTreeSet<_>>();
    let mut membership = (0..shell.len()).collect::<Vec<_>>();
    let mut replacements = FxHashMap::default();
    let mut bounds = Bounds::default();
    let mut rejected_groups = 0;
    while let Some(seed) = unseen.pop_first() {
        let Some((group, group_bounds)) = context.group(seed, &mut unseen) else {
            continue;
        };
        if group.len() < 2 {
            continue;
        }
        if let Some((face, bound)) = context.merge(&group, group_bounds) {
            for index in group {
                membership[index] = seed;
            }
            replacements.insert(seed, face);
            bounds.include(bound);
        } else {
            rejected_groups += 1;
        }
    }
    let merged_groups = replacements.len();
    let mut faces = Vec::new();
    let mut original_to_simplified = vec![0; shell.len()];
    for (index, &seed) in membership.iter().enumerate() {
        if seed == index {
            original_to_simplified[index] = faces.len();
            faces.push(
                replacements
                    .remove(&index)
                    .unwrap_or_else(|| shell[index].clone()),
            );
        } else {
            original_to_simplified[index] = original_to_simplified[seed];
        }
    }
    let simplified = faces.into_iter().collect::<Shell>();
    if simplified.shell_condition() != condition
        || simplified.connected_components().len() != shell.connected_components().len()
        || !simplified.singular_vertices().is_empty()
        || !simplified.is_geometric_consistent()
    {
        return None;
    }
    Some(PlanarFanSimplification {
        shell: simplified,
        original_to_simplified,
        merged_groups,
        rejected_groups,
        deviation_upper_bound: bounds.deviation(),
        seam_deviation_upper_bound: bounds.seam,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        BsplineCurve, KnotVector, NurbsCurve, builder,
        mesh_recovery::{
            BoundaryPatchOptions, PrescribedSourceEdge, SurfaceRegionBudget,
            approximate_boundary_patch,
        },
    };

    #[test]
    fn coplanar_curved_fans_merge_with_original_external_edge_identity()
    -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(0., 0., 0.),
            Point3::new(2., 0., 0.),
            Point3::new(2., 2., 0.),
            Point3::new(0., 2., 0.),
        ];
        let vertices = points.map(builder::vertex);
        let curved = Edge::new(
            &vertices[0],
            &vertices[1],
            Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
                KnotVector::bezier_knot(2),
                vec![
                    Vector4::new(0., 0., 0., 1.),
                    Vector4::new(1., -0.1, 0., 1.),
                    Vector4::new(2., 0., 0., 1.),
                ],
            )?)),
        )?;
        let source: [&[usize]; 2] = [&[0, 1, 2], &[0, 2, 3]];
        let patch = approximate_boundary_patch(
            &points,
            &source,
            &[PrescribedSourceEdge {
                vertices: [0, 1],
                edge: curved.clone(),
            }],
            BoundaryPatchOptions::new(SurfaceRegionBudget {
                vertex: 1e-5,
                chord: 0.2,
            }),
        )?;
        assert_eq!(patch.shell.len(), 6);
        let merged = simplify_planar_fan_shell(&patch.shell, PlanarFanOptions::new(1e-8))
            .ok_or("simplification failed")?;
        assert_eq!(merged.shell.len(), 1);
        assert_eq!(merged.original_to_simplified, vec![0; 6]);
        assert_eq!(merged.merged_groups, 1);
        assert!(merged.deviation_upper_bound < 1e-8);
        assert!(matches!(merged.shell[0].surface(), Surface::Plane(_)));
        assert_eq!(merged.shell[0].boundaries()[0].len(), 4);
        assert_eq!(
            merged
                .shell
                .edge_iter()
                .filter(|e| e.id() == curved.id())
                .count(),
            1
        );
        let owned = simplify_planar_fan_shell(
            &patch.shell,
            PlanarFanOptions {
                owner_by_face: Some(&[0, 0, 0, 1, 1, 1]),
                ..PlanarFanOptions::new(1e-8)
            },
        )
        .ok_or("owned simplification failed")?;
        assert_eq!(owned.shell.len(), 2);
        assert_eq!(owned.original_to_simplified, vec![0, 0, 0, 1, 1, 1]);
        let disabled = simplify_planar_fan_shell(&patch.shell, PlanarFanOptions::new(0.))
            .ok_or("zero-budget simplification failed")?;
        assert_eq!(disabled.shell.len(), 6);
        Ok(())
    }

    #[test]
    fn planar_ring_coalesces_without_filling_its_hole() -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            (0., 0.),
            (4., 0.),
            (4., 4.),
            (0., 4.),
            (1., 1.),
            (3., 1.),
            (3., 3.),
            (1., 3.),
        ]
        .map(|(x, y)| Point3::new(x, y, 0.));
        let polygons = (0..4)
            .flat_map(|i| {
                let j = (i + 1) % 4;
                [[i, j, j + 4], [i, j + 4, i + 4]]
            })
            .collect::<Vec<_>>();
        let source = polygons.iter().map(|p| p.as_slice()).collect::<Vec<_>>();
        let patch = approximate_boundary_patch(
            &points,
            &source,
            &[],
            BoundaryPatchOptions::new(SurfaceRegionBudget {
                vertex: 1e-5,
                chord: 0.1,
            }),
        )?;
        let result = simplify_planar_fan_shell(&patch.shell, PlanarFanOptions::new(1e-8))
            .ok_or("ring simplification")?;
        assert_eq!(result.shell.len(), 1);
        assert_eq!(result.shell[0].boundaries().len(), 2);
        assert_eq!(result.original_to_simplified, vec![0; 24]);
        assert_eq!(result.shell[0].boundaries()[0].len(), 4);
        assert_eq!(result.shell[0].boundaries()[1].len(), 4);
        Ok(())
    }

    #[test]
    fn nonplanar_source_fan_preserves_bounded_geometry() -> Result<(), Box<dyn std::error::Error>> {
        let points = [
            Point3::new(0., 0., 0.),
            Point3::new(2., 0., 0.),
            Point3::new(0., 2., 0.),
        ];
        let vertices = points.map(builder::vertex);
        let curved = Edge::new(
            &vertices[0],
            &vertices[1],
            Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
                KnotVector::bezier_knot(2),
                vec![
                    Vector4::new(0., 0., 0., 1.),
                    Vector4::new(1., 0., 0.1, 1.),
                    Vector4::new(2., 0., 0., 1.),
                ],
            )?)),
        )?;
        let source: [&[usize]; 1] = [&[0, 1, 2]];
        let patch = approximate_boundary_patch(
            &points,
            &source,
            &[PrescribedSourceEdge {
                vertices: [0, 1],
                edge: curved.clone(),
            }],
            BoundaryPatchOptions::new(SurfaceRegionBudget {
                vertex: 1e-5,
                chord: 0.2,
            }),
        )?;
        let result = simplify_planar_fan_shell(&patch.shell, PlanarFanOptions::new(1e-8))
            .ok_or("simplification failed")?;
        assert_eq!(result.shell.len(), 2);
        assert_eq!(result.shell[0], patch.shell[0]);
        assert!(matches!(
            result.shell[0].surface(),
            Surface::NurbsSurface(_)
        ));
        assert!(matches!(result.shell[1].surface(), Surface::Plane(_)));
        Ok(())
    }
}
