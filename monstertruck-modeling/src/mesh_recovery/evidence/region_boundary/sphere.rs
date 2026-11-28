//! Axis-independent spherical region hypotheses from adjacent cyclic facets.
//! A single sphere through four points is interpolation; independent
//! quadratic sample rank distinguishes it from competing quadric carriers.

use super::super::super::sample_fit::{CircleCarrierFit, fit_circle_carrier, normal_basis};
use super::{
    MeshEvidence, RegionBoundaryIndex, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionProposal,
};
use crate::{InnerSpace, Point3, Vector3};
use std::collections::{BTreeSet, VecDeque};

#[derive(Clone, Copy)]
struct SphereSeed {
    center: Point3,
    radius: f64,
}

impl SphereSeed {
    fn from_circles(a: CircleCarrierFit, b: CircleCarrierFit, tolerance: f64) -> Option<Self> {
        let direction = a.normal.cross(b.normal);
        let denominator = direction.magnitude2();
        if denominator <= 64.0 * f64::EPSILON {
            return None;
        }
        let delta = b.center - a.center;
        let ta = delta.cross(b.normal).dot(direction) / denominator;
        let tb = delta.cross(a.normal).dot(direction) / denominator;
        let ca = a.center + a.normal * ta;
        let cb = b.center + b.normal * tb;
        if (ca - cb).magnitude() > tolerance {
            return None;
        }
        let center = ca + (cb - ca) * 0.5;
        let ra = (center - a.center).magnitude().hypot(a.radius);
        let rb = (center - b.center).magnitude().hypot(b.radius);
        let radius = ra.midpoint(rb);
        (radius.is_finite() && radius > 0.0 && (ra - rb).abs() <= tolerance)
            .then_some(Self { center, radius })
    }

    fn same(self, other: Self, tolerance: f64) -> bool {
        (self.center - other.center).magnitude() <= tolerance
            && (self.radius - other.radius).abs() <= tolerance
    }
}

// Fit the center in a source-local frame by solving the centered
// squared-distance equations directly with pivoted, reorthogonalized QR.
// Centering eliminates the constant coefficient without forming the
// normal equations (which would square the condition number).
fn fit_sphere(points: &[Point3], ids: &[usize]) -> Option<SphereSeed> {
    if ids.len() < 5 {
        return None;
    }
    let first = *points.get(*ids.first()?)?;
    let count = crate::usize_to_f64(ids.len());
    let mean = ids
        .iter()
        .map(|&id| points[id] - first)
        .fold(Vector3::new(0.0, 0.0, 0.0), |sum, p| sum + p)
        / count;
    let origin = first + mean;
    let local = ids
        .iter()
        .map(|&id| points[id] - origin)
        .collect::<Vec<_>>();
    let scale = local.iter().map(|p| p.magnitude()).fold(0.0_f64, f64::max);
    if scale <= 0.0 || !scale.is_finite() {
        return None;
    }
    let local = local.into_iter().map(|p| p / scale).collect::<Vec<_>>();
    let mean_squared = local.iter().map(|p| p.magnitude2()).sum::<f64>() / count;
    let rhs = local
        .iter()
        .map(|p| 0.5 * (p.magnitude2() - mean_squared))
        .collect::<Vec<_>>();
    let mut columns = [0, 1, 2].map(|axis| local.iter().map(|p| p[axis]).collect::<Vec<_>>());
    let dot = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .fold(0.0_f64, |sum, (&a, &b)| a.mul_add(b, sum))
    };
    let mut permutation = [0, 1, 2];
    let mut upper = [[0.0_f64; 3]; 3];
    let mut projected = [0.0_f64; 3];
    let initial_norm = columns
        .iter()
        .map(|c| dot(c, c))
        .fold(0.0_f64, f64::max)
        .sqrt();
    for i in 0..3 {
        let pivot = (i..3).max_by(|&a, &b| {
            dot(&columns[a], &columns[a]).total_cmp(&dot(&columns[b], &columns[b]))
        })?;
        columns.swap(i, pivot);
        permutation.swap(i, pivot);
        for row in upper.iter_mut().take(i) {
            row.swap(i, pivot);
        }
        let norm = dot(&columns[i], &columns[i]).sqrt();
        if norm <= initial_norm * 128.0 * f64::EPSILON || !norm.is_finite() {
            return None;
        }
        upper[i][i] = norm;
        for v in &mut columns[i] {
            *v /= norm;
        }
        projected[i] = dot(&columns[i], &rhs);
        for j in i + 1..3 {
            for _ in 0..2 {
                let projection = dot(&columns[i], &columns[j]);
                upper[i][j] += projection;
                for k in 0..ids.len() {
                    columns[j][k] = (-projection).mul_add(columns[i][k], columns[j][k]);
                }
            }
        }
    }
    let mut solution = [0.0_f64; 3];
    for i in (0..3).rev() {
        solution[i] = ((i + 1..3).fold(projected[i], |sum, j| {
            (-upper[i][j]).mul_add(solution[j], sum)
        })) / upper[i][i];
    }
    let mut offset = Vector3::new(0.0, 0.0, 0.0);
    for (i, axis) in permutation.into_iter().enumerate() {
        offset[axis] = solution[i] * scale;
    }
    let center = origin + offset;
    let radius = ids
        .iter()
        .map(|&id| (points[id] - center).magnitude())
        .sum::<f64>()
        / count;
    (radius > 0.0 && radius.is_finite()).then_some(SphereSeed { center, radius })
}

// Bound the residual of every point of an original polygon, not merely
// its vertices or perimeter. A planar disk with a co-spherical rim is not
// necessarily a good spherical approximation. Orthogonal projection onto
// the support plane preserves the polygon's filled domain; the support's
// measured point-to-plane envelope bounds displacement from that projection.
pub(super) fn facet_residual_bound(
    evidence: &MeshEvidence,
    support: usize,
    center: Point3,
    radius: f64,
) -> Option<f64> {
    let plane = evidence.planar_supports.get(support)?;
    let offset = (plane.origin - center).dot(plane.normal);
    let (u, v) = normal_basis(plane.normal);
    let mut bound = 0.0_f64;
    for &polygon in &plane.source_face_indices {
        let vertices = evidence.source_polygon_vertices(polygon)?;
        let mut maximum = 0.0_f64;
        for &id in vertices {
            maximum = maximum.max((evidence.positions.get(id)? - center).magnitude());
        }
        // Distance to the sphere is bounded by its maximum at the farthest
        // vertex and its minimum radial distance over the filled polygon.
        let axial_lower = (offset.abs() - plane.max_vertex_plane_residual).max(0.0);
        let transverse_lower = projected_polygon_distance(evidence, vertices, center, u, v)?;
        let minimum = axial_lower.hypot(transverse_lower);
        bound = bound
            .max((maximum - radius).abs())
            .max((radius - minimum).max(0.0));
    }
    bound.is_finite().then_some(bound)
}

// Distance from the projected center to a filled (possibly concave)
// polygon. Its closest point is either the center itself or on an edge.
// Use a half-open winding rule and adaptive orientation for crossings.
fn projected_polygon_distance(
    evidence: &MeshEvidence,
    ids: &[usize],
    center: Point3,
    u: Vector3,
    v: Vector3,
) -> Option<f64> {
    use robust::{Coord, orient2d};
    let project = |id: usize| {
        let delta = *evidence.positions.get(id)? - center;
        Some(Coord {
            x: delta.dot(u),
            y: delta.dot(v),
        })
    };
    let mut distance = f64::INFINITY;
    let mut winding = 0_i64;
    let zero = Coord { x: 0.0, y: 0.0 };
    for (&a, &b) in ids.iter().zip(ids.iter().cycle().skip(1)) {
        let start = project(a)?;
        let end = project(b)?;
        let dx = end.x - start.x;
        let dy = end.y - start.y;
        let length2 = dx.mul_add(dx, dy * dy);
        let t = if length2 > 0.0 {
            (-start.x.mul_add(dx, start.y * dy) / length2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        distance = distance.min(dx.mul_add(t, start.x).hypot(dy.mul_add(t, start.y)));
        let turn = orient2d(start, end, zero);
        if start.y <= 0.0 && end.y > 0.0 && turn > 0.0 {
            winding += 1;
        } else if start.y > 0.0 && end.y <= 0.0 && turn < 0.0 {
            winding -= 1;
        }
    }
    if winding != 0 {
        Some(0.0)
    } else {
        distance.is_finite().then_some(distance)
    }
}

struct SupportSamples {
    vertices: Vec<usize>,
    circle: Option<CircleCarrierFit>,
}

fn support_vertices(evidence: &MeshEvidence, support: usize) -> Option<Vec<usize>> {
    let mut ids = BTreeSet::new();
    for &polygon in &evidence.planar_supports.get(support)?.source_face_indices {
        ids.extend(evidence.source_polygon_vertices(polygon)?.iter().copied());
    }
    Some(ids.into_iter().collect())
}

/// Test the rank of the homogeneous quadratic sample design matrix.
/// A sphere is a quadric, so exact source samples have rank at most nine.
/// Rank nine identifies that quadratic carrier uniquely (up to scale).
/// Two rings of a cylinder also lie on a sphere but have additional
/// independent quadratic equations, and therefore do not pass this test.
pub(super) fn independently_tested(
    evidence: &MeshEvidence,
    supports: &[usize],
    residual: f64,
) -> bool {
    let mut ids = BTreeSet::new();
    for &support in supports {
        let Some(samples) = support_vertices(evidence, support) else {
            return false;
        };
        ids.extend(samples);
    }
    let Some(&first) = ids.first().and_then(|&id| evidence.positions.get(id)) else {
        return false;
    };
    if ids.len() < 9 {
        return false;
    }
    let span = ids
        .iter()
        .map(|&id| (evidence.positions[id] - first).magnitude())
        .fold(0.0_f64, f64::max);
    if !span.is_finite() || span <= 0.0 {
        return false;
    }
    let mut rows = ids
        .iter()
        .map(|&id| {
            let p = (evidence.positions[id] - first) / span;
            [
                1.0,
                p.x,
                p.y,
                p.z,
                p.x * p.x,
                p.y * p.y,
                p.z * p.z,
                p.x * p.y,
                p.x * p.z,
                p.y * p.z,
            ]
        })
        .collect::<Vec<_>>();
    let tolerance = 128.0_f64.mul_add(f64::EPSILON, 16.0 * residual / span);
    // Pivot by maximum remaining row norm; operate on the normalized
    // design matrix directly rather than squaring its condition number.
    let squared_norm = |row: &[f64; 10]| row.iter().map(|x| x * x).sum::<f64>();
    for _ in 0..9 {
        let Some((pivot, norm)) = rows
            .iter()
            .map(|r| (*r, squared_norm(r)))
            .max_by(|a, b| a.1.total_cmp(&b.1))
        else {
            return false;
        };
        if norm <= tolerance * tolerance {
            return false;
        }
        let basis = pivot.map(|x| x / norm.sqrt());
        for row in &mut rows {
            let dot = row.iter().zip(basis).map(|(a, b)| a * b).sum::<f64>();
            for (x, b) in row.iter_mut().zip(basis) {
                *x = (-dot).mul_add(b, *x);
            }
        }
    }
    true
}

struct SphereDiscovery<'a> {
    evidence: &'a MeshEvidence,
    samples: Vec<SupportSamples>,
    seen: Vec<Vec<usize>>,
    seeds: Vec<SphereSeed>,
    relative_tolerance: f64,
    absolute_tolerance: f64,
    chord_fraction: f64,
}

impl SphereDiscovery<'_> {
    fn strict_tolerance(&self, sphere: SphereSeed) -> f64 {
        self.absolute_tolerance.min(self.surface_tolerance(sphere))
    }

    fn surface_tolerance(&self, sphere: SphereSeed) -> f64 {
        self.relative_tolerance * sphere.radius
    }

    fn accepts(&self, support: usize, sphere: SphereSeed, outward: bool, tolerance: f64) -> bool {
        let sample = &self.samples[support];
        // A support need not itself be cyclic (e.g. a tessellated seam
        // with extra coplanar samples). Test its actual source vertices.
        let plane = &self.evidence.planar_supports[support];
        let offset = (plane.origin - sphere.center).dot(plane.normal);
        if offset.abs() <= tolerance || (offset > 0.0) != outward {
            return false;
        }
        sample.vertices.iter().all(|&id| {
            ((self.evidence.positions[id] - sphere.center).magnitude() - sphere.radius).abs()
                <= tolerance
        }) && facet_residual_bound(self.evidence, support, sphere.center, sphere.radius)
            .is_some_and(|error| error <= self.chord_fraction * sphere.radius)
    }

    fn grow(&self, seed: usize, sphere: SphereSeed, outward: bool, tolerance: f64) -> Vec<usize> {
        let mut tested = BTreeSet::from([seed]);
        let mut accepted = Vec::new();
        let mut queue = VecDeque::from([seed]);
        while let Some(id) = queue.pop_front() {
            if !self.accepts(id, sphere, outward, tolerance) {
                continue;
            }
            accepted.push(id);
            for &interface in &self.evidence.support_interface_indices[id] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let other = if a == id { b } else { a };
                if tested.insert(other) {
                    queue.push_back(other);
                }
            }
        }
        accepted.sort_unstable();
        accepted
    }
    // An initial local fit may have used perturbed trim samples.
    // The expanded candidate supplies interior vertices independently
    // of its boundary. Refit only those vertices, then independently
    // recover a precise connected core and revalidate the entire union.
    fn refit_interior(
        &self,
        index: &RegionBoundaryIndex<'_>,
        expanded: &[usize],
        seed: SphereSeed,
        outward: bool,
    ) -> Option<(SphereSeed, Vec<usize>, Vec<usize>)> {
        let boundary = index.boundary(expanded).ok()?;
        let boundary_vertices = boundary
            .loops
            .iter()
            .flatten()
            .flat_map(|edge| edge.vertices)
            .collect::<BTreeSet<_>>();
        let interior = expanded
            .iter()
            .flat_map(|&id| self.samples[id].vertices.iter().copied())
            .filter(|id| !boundary_vertices.contains(id))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let fitted = fit_sphere(&self.evidence.positions, &interior)?;
        let error = |sphere: SphereSeed| {
            interior
                .iter()
                .map(|&id| {
                    ((self.evidence.positions[id] - sphere.center).magnitude() - sphere.radius)
                        .abs()
                })
                .fold(0.0_f64, f64::max)
        };
        let tolerance = self.strict_tolerance(fitted);
        // A new carrier must be better supported on the interior, not
        // merely interpolate a different subset of seam samples.
        if error(fitted) > tolerance || error(fitted) >= error(seed) {
            return None;
        }
        let mut checked = BTreeSet::new();
        let mut best = Vec::new();
        for &id in expanded {
            if checked.contains(&id) {
                continue;
            }
            let core = self.grow(id, fitted, outward, tolerance);
            checked.insert(id);
            checked.extend(core.iter().copied());
            if core.len() > best.len() && independently_tested(self.evidence, &core, tolerance) {
                best = core;
            }
        }
        let &start = best.first()?;
        let grown = self.grow(start, fitted, outward, self.surface_tolerance(fitted));
        if !expanded.iter().all(|id| grown.binary_search(id).is_ok()) {
            return None;
        }
        Some((fitted, best, grown))
    }
}

impl<'a> RegionBoundaryIndex<'a> {
    /// Generate connected spherical alternatives without an assumed axis.
    /// Every candidate must pass the existing full source chord and topology
    /// checks, and the independent quadratic rank test.
    pub(in crate::mesh_recovery) fn discover_spherical_regions(
        &self,
        relative_tolerance: f64,
        absolute_tolerance: f64,
        chord_fraction: f64,
    ) -> Vec<SurfaceRegionProposal<'a>> {
        let samples = (0..self.evidence.planar_supports.len())
            .map(|id| {
                let vertices = support_vertices(self.evidence, id).unwrap_or_default();
                let circle = fit_circle_carrier(&self.evidence.positions, &vertices).filter(|c| {
                    c.max_plane_residual <= absolute_tolerance
                        && c.max_radial_residual <= absolute_tolerance
                });
                SupportSamples { vertices, circle }
            })
            .collect::<Vec<_>>();
        let mut discovery = SphereDiscovery {
            evidence: self.evidence,
            seen: vec![Vec::new(); samples.len()],
            samples,
            seeds: Vec::new(),
            relative_tolerance,
            absolute_tolerance,
            chord_fraction,
        };
        let mut proposals = Vec::new();
        for interface in &self.evidence.support_interfaces {
            let [a, b] = interface.supports;
            let (Some(first), Some(second)) =
                (discovery.samples[a].circle, discovery.samples[b].circle)
            else {
                continue;
            };
            let Some(sphere) = SphereSeed::from_circles(first, second, absolute_tolerance) else {
                continue;
            };
            let tolerance = discovery.strict_tolerance(sphere);
            if discovery.seen[a].iter().any(|&id| {
                discovery.seen[b].contains(&id) && discovery.seeds[id].same(sphere, tolerance)
            }) {
                continue;
            }
            let outward =
                (first.center - sphere.center).dot(self.evidence.planar_supports[a].normal) > 0.0;
            if !discovery.accepts(a, sphere, outward, tolerance)
                || !discovery.accepts(b, sphere, outward, tolerance)
            {
                continue;
            }
            let supports = discovery.grow(a, sphere, outward, tolerance);
            if supports.binary_search(&b).is_err() {
                continue;
            }
            let id = discovery.seeds.len();
            discovery.seeds.push(sphere);
            for &support in &supports {
                discovery.seen[support].push(id);
            }
            let expanded = discovery.grow(a, sphere, outward, discovery.surface_tolerance(sphere));
            let (sphere, supports, expanded) = discovery
                .refit_interior(self, &expanded, sphere, outward)
                .unwrap_or((sphere, supports, expanded));
            let tolerance = discovery.strict_tolerance(sphere);
            let carrier = SurfaceCarrier::CircularRevolution {
                axial_center: sphere.center,
                axis: Vector3::unit_z(),
                radial_center: 0.0,
                profile_radius: sphere.radius,
            };
            if !independently_tested(self.evidence, &supports, tolerance) {
                continue;
            }
            if expanded != supports
                && let Some(proposal) = self.measure_surface_region(
                    &expanded,
                    carrier,
                    SurfaceRegionBudget {
                        vertex: discovery.surface_tolerance(sphere),
                        chord: discovery.chord_fraction * sphere.radius,
                    },
                )
            {
                // The strict, connected core is a subset of the expanded
                // support union and independently tests the same carrier.
                proposals.push(proposal.with_independent_profile());
                // A second seed of this same sphere inside the expansion
                // cannot add a new maximal connected candidate.
                for support in expanded {
                    if !discovery.seen[support].contains(&id) {
                        discovery.seen[support].push(id);
                    }
                }
                continue;
            }
            if let Some(proposal) = self.measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: tolerance,
                    chord: discovery.chord_fraction * sphere.radius,
                },
            ) {
                proposals.push(proposal.with_independent_profile());
            }
        }
        proposals
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{EvidenceOptions, analyze_mesh_evidence};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use std::f64::consts::TAU;

    fn octant_mesh(
        steps: usize,
        center: Point3,
        radius: f64,
        boundary_offset: f64,
    ) -> Result<PolygonMesh, Box<dyn std::error::Error>> {
        // Barycentric subdivision projected to a sphere: no assumption
        // that samples lie on a set of latitudes around a preferred axis.
        let mut positions = Vec::new();
        let mut rows = Vec::new();
        for i in 0..=steps {
            let mut row = Vec::new();
            for j in 0..=steps - i {
                let v = Vector3::new(
                    crate::usize_to_f64(i),
                    crate::usize_to_f64(j),
                    crate::usize_to_f64(steps - i - j),
                )
                .normalize();
                // Cyclic permutation and arbitrary rigid placement.
                let rotated = Vector3::new(v.z, v.x, v.y);
                row.push(positions.len());
                let on_boundary = i == 0 || j == 0 || i + j == steps;
                let r = if on_boundary {
                    radius + boundary_offset
                } else {
                    radius
                };
                positions.push(center + rotated * r);
            }
            rows.push(row);
        }
        let mut triangles = Vec::new();
        for i in 0..steps {
            for j in 0..steps - i {
                let a = rows[i][j];
                let b = rows[i + 1][j];
                let c = rows[i][j + 1];
                triangles.push([a, b, c]);
                if j + 1 < steps - i {
                    triangles.push([b, rows[i + 1][j + 1], c]);
                }
            }
        }
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(triangles),
        )?)
    }

    #[test]
    fn independent_sphere_recovers_unstructured_spherical_patch()
    -> Result<(), Box<dyn std::error::Error>> {
        let steps = 6;
        let center = Point3::new(17.0, -3.0, 5.0);
        let radius = 2.5;
        let mesh = octant_mesh(steps, center, radius, 0.0)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let proposals = index.discover_spherical_regions(1.0e-4, 1.0e-8, 0.3);
        let proposal = proposals
            .iter()
            .find(|p| p.source_polygons() == steps * steps)
            .ok_or("sphere proposal")?;
        let SurfaceCarrier::CircularRevolution {
            axial_center,
            radial_center,
            profile_radius,
            ..
        } = proposal.carrier()
        else {
            return Err("unexpected carrier".into());
        };
        assert!((axial_center - center).magnitude() < 1.0e-10);
        assert!(radial_center.abs() < 1.0e-10 && (profile_radius - radius).abs() < 1.0e-10);
        // A numerically tiny fitted radial offset must not survive as
        // a toroidal implicit constraint at a spherical pole.
        let measured = index
            .measure_surface_region(
                proposal.supports(),
                SurfaceCarrier::CircularRevolution {
                    axial_center: center,
                    axis: Vector3::unit_z(),
                    radial_center: crate::KERNEL_EPSILON * 0.25,
                    profile_radius: radius,
                },
                SurfaceRegionBudget {
                    vertex: 1.0e-10,
                    chord: 0.5,
                },
            )
            .ok_or("canonical spherical proposal")?;
        assert!(matches!(
            measured.carrier(),
            SurfaceCarrier::CircularRevolution {
                radial_center: 0.0,
                ..
            }
        ));
        let pole = center + Vector3::unit_z() * radius;
        let (normal, residual) = measured
            .carrier()
            .signed_constraint(pole)
            .ok_or("sphere pole")?;
        assert_eq!(normal, Vector3::unit_z());
        assert_eq!(residual, 0.0);
        let recovered = crate::mesh_recovery::recover_analytic_mesh(&mesh, &Default::default())?;
        assert_eq!(recovered.native.regions.len(), 1);
        assert!(recovered.complete_shell().is_some());
        Ok(())
    }

    #[test]
    fn precise_core_extends_to_bounded_seam_samples() -> Result<(), Box<dyn std::error::Error>> {
        let steps = 12;
        let center = Point3::new(-1.0, 3.0, 17.0);
        let radius = 2.5;
        let seam_displacement = 2.0e-5;
        let mesh = octant_mesh(steps, center, radius, seam_displacement)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let proposals = index.discover_spherical_regions(1.0e-4, 1.0e-8, 0.3);
        let proposal = proposals
            .iter()
            .find(|p| p.source_polygons() == steps * steps)
            .ok_or("expanded spherical proposal")?;
        assert!(proposal.independently_tested_profile());
        let SurfaceCarrier::CircularRevolution {
            axial_center,
            radial_center,
            profile_radius,
            ..
        } = proposal.carrier()
        else {
            return Err("unexpected carrier".into());
        };
        assert!((axial_center - center).magnitude() < 1.0e-9);
        assert!(radial_center.abs() < 1.0e-9 && (profile_radius - radius).abs() < 1.0e-9);
        let max_error = proposal.relative_vertex_error() * radius * 1.0e-4;
        assert!((max_error - seam_displacement).abs() < 1.0e-9);
        // The identical geometry must not be admitted as a complete
        // region when the configured surface budget excludes the seam.
        let strict = index.discover_spherical_regions(1.0e-6, 1.0e-8, 0.3);
        assert!(strict.iter().all(|p| p.source_polygons() < steps * steps));
        Ok(())
    }

    #[test]
    fn interior_refit_removes_boundary_induced_carrier_bias()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(17.0, -3.0, 5.0);
        let radius = 0.04;
        let steps = 12;
        let mesh = octant_mesh(steps, center, radius, 2.0e-6)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let samples = (0..evidence.planar_supports.len())
            .map(|id| SupportSamples {
                vertices: support_vertices(&evidence, id).unwrap_or_default(),
                circle: None,
            })
            .collect::<Vec<_>>();
        let discovery = SphereDiscovery {
            evidence: &evidence,
            seen: vec![Vec::new(); samples.len()],
            samples,
            seeds: Vec::new(),
            relative_tolerance: 1.0e-4,
            absolute_tolerance: 1.0e-8,
            chord_fraction: 0.3,
        };
        let biased = SphereSeed {
            center: center + Vector3::new(1.0e-6, -0.5e-6, 0.7e-6),
            radius: radius + 0.5e-6,
        };
        let supports = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        let outward = (evidence.planar_supports[0].origin - center)
            .dot(evidence.planar_supports[0].normal)
            > 0.0;
        let (fitted, core, expanded) = discovery
            .refit_interior(&index, &supports, biased, outward)
            .ok_or("refit")?;
        assert!((fitted.center - center).magnitude() < 1.0e-12);
        assert!((fitted.radius - radius).abs() < 1.0e-12);
        assert_eq!(expanded, supports);
        assert!(core.len() < expanded.len());
        assert!(independently_tested(&evidence, &core, 1.0e-8));
        Ok(())
    }

    #[test]
    fn spherical_refit_rejects_coplanar_samples() {
        let points = (0..12)
            .map(|i| {
                let (s, c) = (TAU * f64::from(i) / 12.0).sin_cos();
                Point3::new(c, s, 2.0)
            })
            .collect::<Vec<_>>();
        assert!(fit_sphere(&points, &(0..points.len()).collect::<Vec<_>>()).is_none());
    }

    #[test]
    fn cospherical_planar_closure_does_not_poison_spherical_growth()
    -> Result<(), Box<dyn std::error::Error>> {
        let rings = 8_usize;
        let sides = 24_usize;
        let mut positions = vec![Point3::new(0.0, 0.0, 1.0)];
        let mut polygons = Vec::new();
        for ring in 1..=rings {
            // More than a hemisphere: the planar closure and the sphere
            // have the same oriented support-plane offset sign.
            let theta = TAU / 3.0 * crate::usize_to_f64(ring) / crate::usize_to_f64(rings);
            let (rho, z) = theta.sin_cos();
            for j in 0..sides {
                let (sin, cos) =
                    (TAU * crate::usize_to_f64(j) / crate::usize_to_f64(sides)).sin_cos();
                positions.push(Point3::new(rho * cos, rho * sin, z));
            }
        }
        for j in 0..sides {
            polygons.push(vec![0, 1 + j, 1 + (j + 1) % sides]);
        }
        for ring in 1..rings {
            for j in 0..sides {
                let a = 1 + (ring - 1) * sides;
                let b = a + sides;
                let k = (j + 1) % sides;
                polygons.push(vec![a + j, b + j, b + k, a + k]);
            }
        }
        let spherical_polygons = polygons.len();
        let cap_start = 1 + (rings - 1) * sides;
        let rim = (cap_start..cap_start + sides).rev().collect::<Vec<_>>();
        for triangulated in [false, true] {
            let mut faces = polygons.clone();
            if triangulated {
                faces.extend(
                    rim[1..]
                        .windows(2)
                        .map(|pair| vec![rim[0], pair[0], pair[1]]),
                );
            } else {
                faces.push(rim.clone());
            }
            let mesh = PolygonMesh::new(
                StandardAttributes {
                    positions: positions.clone(),
                    ..Default::default()
                },
                Faces::from_iter(faces),
            )?;
            let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
            let index = RegionBoundaryIndex::new(&evidence)?;
            let cap = evidence
                .planar_supports
                .iter()
                .position(|support| {
                    support.source_face_indices.iter().all(|&polygon| {
                        evidence
                            .source_polygon_vertices(polygon)
                            .is_some_and(|vertices| {
                                vertices
                                    .iter()
                                    .all(|&id| (evidence.positions[id].z + 0.5).abs() < 1.0e-10)
                            })
                    })
                })
                .ok_or("planar closure support")?;
            let center = Point3::new(0.0, 0.0, 0.0);
            let carrier = SurfaceCarrier::CircularRevolution {
                axial_center: center,
                axis: Vector3::unit_z(),
                radial_center: 0.0,
                profile_radius: 1.0,
            };
            let residual =
                facet_residual_bound(&evidence, cap, center, 1.0).ok_or("closure facet bound")?;
            assert!(
                (residual - 0.5).abs() < 1.0e-9,
                "closure residual {residual}"
            );
            if !triangulated {
                let perimeter = carrier
                    .measure_source_chords(&evidence, &[cap], 0.3)
                    .ok_or("closure perimeter")?;
                assert!(perimeter.upper_bound < 0.3);
            }
            let proposals = index.discover_spherical_regions(1.0e-4, 1.0e-8, 0.3);
            let found = proposals
                .iter()
                .find(|p| p.source_polygons() == spherical_polygons)
                .ok_or("sphere without its planar closure")?;
            assert!(!found.supports().contains(&cap));
            assert!(found.independently_tested_profile());
            let mut union = found.supports().to_vec();
            union.push(cap);
            union.sort_unstable();
            assert!(
                index
                    .measure_surface_region(
                        &union,
                        carrier,
                        SurfaceRegionBudget {
                            vertex: 1.0e-8,
                            chord: 0.3,
                        }
                    )
                    .is_none()
            );
        }
        Ok(())
    }

    #[test]
    fn two_ring_cylinder_does_not_certify_a_sphere() -> Result<(), Box<dyn std::error::Error>> {
        let ring = 12_usize;
        let mut positions = Vec::new();
        for z in [-0.25, 0.25] {
            for i in 0..ring {
                let (sin, cos) =
                    (TAU * crate::usize_to_f64(i) / crate::usize_to_f64(ring)).sin_cos();
                positions.push(Point3::new(cos, sin, z));
            }
        }
        let faces = (0..ring).map(|i| [i, (i + 1) % ring, (i + 1) % ring + ring, i + ring]);
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let supports = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        assert!(!independently_tested(&evidence, &supports, 1.0e-14));
        let index = RegionBoundaryIndex::new(&evidence)?;
        assert!(
            index
                .discover_spherical_regions(1.0e-4, 1.0e-8, 0.3)
                .is_empty()
        );
        Ok(())
    }
}
