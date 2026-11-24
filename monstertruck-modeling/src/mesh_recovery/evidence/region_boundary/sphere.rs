//! Axis-independent spherical region hypotheses from adjacent cyclic facets.
//! A single sphere through four points is interpolation; independent
//! quadratic sample rank distinguishes it from competing quadric carriers.

use super::super::super::sample_fit::{CircleCarrierFit, fit_circle_carrier};
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
    fn tolerance(&self, sphere: SphereSeed) -> f64 {
        self.absolute_tolerance
            .min(self.relative_tolerance * sphere.radius)
    }

    fn accepts(&self, support: usize, sphere: SphereSeed, outward: bool) -> bool {
        let sample = &self.samples[support];
        let Some(circle) = sample.circle else {
            return false;
        };
        let tolerance = self.tolerance(sphere);
        let normal = self.evidence.planar_supports[support].normal;
        let offset = (circle.center - sphere.center).dot(normal);
        if offset.abs() <= tolerance || (offset > 0.0) != outward {
            return false;
        }
        sample.vertices.iter().all(|&id| {
            ((self.evidence.positions[id] - sphere.center).magnitude() - sphere.radius).abs()
                <= tolerance
        })
    }

    fn grow(&self, seed: usize, sphere: SphereSeed, outward: bool) -> Vec<usize> {
        let mut tested = BTreeSet::from([seed]);
        let mut accepted = Vec::new();
        let mut queue = VecDeque::from([seed]);
        while let Some(id) = queue.pop_front() {
            if !self.accepts(id, sphere, outward) {
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
            let tolerance = discovery.tolerance(sphere);
            if discovery.seen[a].iter().any(|&id| {
                discovery.seen[b].contains(&id) && discovery.seeds[id].same(sphere, tolerance)
            }) {
                continue;
            }
            let outward =
                (first.center - sphere.center).dot(self.evidence.planar_supports[a].normal) > 0.0;
            if !discovery.accepts(a, sphere, outward) || !discovery.accepts(b, sphere, outward) {
                continue;
            }
            let supports = discovery.grow(a, sphere, outward);
            if supports.binary_search(&b).is_err() {
                continue;
            }
            let id = discovery.seeds.len();
            discovery.seeds.push(sphere);
            for &support in &supports {
                discovery.seen[support].push(id);
            }
            let carrier = SurfaceCarrier::CircularRevolution {
                axial_center: sphere.center,
                axis: Vector3::unit_z(),
                radial_center: 0.0,
                profile_radius: sphere.radius,
            };
            if let Some(proposal) = self.measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: tolerance,
                    chord: discovery.chord_fraction * sphere.radius,
                },
            ) && proposal.independently_tested_profile()
            {
                proposals.push(proposal);
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

    #[test]
    fn independent_sphere_recovers_unstructured_spherical_patch()
    -> Result<(), Box<dyn std::error::Error>> {
        let steps = 6_usize;
        let center = Point3::new(17.0, -3.0, 5.0);
        let radius = 2.5;
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
                positions.push(center + rotated * radius);
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
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(triangles),
        )?;
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
        let recovered = crate::mesh_recovery::recover_analytic_mesh(&mesh, Default::default())?;
        assert_eq!(recovered.native.regions.len(), 1);
        assert!(recovered.complete_shell().is_some());
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
