//! Extend independently measured circular-meridian carriers across
//! irregular source facets without allowing the new samples to refit or
//! independently establish the underlying profile.

use super::{
    MeshEvidence, RegionBoundaryIndex, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionProposal,
};
use crate::{InnerSpace, Vector3};
use std::collections::{BTreeSet, VecDeque};

fn orientation(
    evidence: &MeshEvidence,
    id: usize,
    carrier: SurfaceCarrier,
    budget: SurfaceRegionBudget,
    angular_tolerance: f64,
) -> Option<bool> {
    let support = evidence.planar_supports.get(id)?;
    let mut sign = None;
    for &polygon in &support.source_face_indices {
        let ids = evidence.source_polygon_vertices(polygon)?;
        let first = *evidence.positions.get(*ids.first()?)?;
        let mut center = Vector3::new(0.0, 0.0, 0.0);
        let mut normals = Vec::with_capacity(ids.len());
        for &vertex in ids {
            let point = *evidence.positions.get(vertex)?;
            let (normal, residual) = carrier.signed_constraint(point)?;
            if residual.abs() > budget.vertex {
                return None;
            }
            normals.push(normal);
            center += point - first;
        }
        let centroid = first + center / crate::usize_to_f64(ids.len());
        let (expected, _) = carrier.signed_constraint(centroid)?;
        let dot = expected.dot(support.normal);
        let deviation = expected.cross(support.normal).magnitude().atan2(dot.abs());
        let span = normals
            .into_iter()
            .map(|normal| {
                expected
                    .cross(normal)
                    .magnitude()
                    .atan2(expected.dot(normal))
            })
            .fold(0.0_f64, f64::max);
        if dot.abs() <= angular_tolerance.sin()
            || deviation > span + angular_tolerance + support.max_normal_angle_radians
            || sign.is_some_and(|previous| previous != (dot > 0.0))
        {
            return None;
        }
        sign = Some(dot > 0.0);
    }
    carrier
        .measure_circular_source_facets(evidence, &[id], budget.chord)
        .filter(|bound| bound.upper_bound <= budget.chord)?;
    sign
}

impl<'a> RegionBoundaryIndex<'a> {
    // Preserve the fixed placed core carrier and its independent evidence.
    // Each additional whole source support must meet both local vertex and
    // filled-facet budgets and the same oriented normal envelope.
    pub(in crate::mesh_recovery) fn expand_circular_region(
        &self,
        core: &SurfaceRegionProposal<'a>,
        budget: SurfaceRegionBudget,
        angular_tolerance: f64,
    ) -> Option<SurfaceRegionProposal<'a>> {
        let carrier = core.carrier();
        let SurfaceCarrier::CircularRevolution { radial_center, .. } = carrier else {
            return None;
        };
        // Spheres already have an independent core/growth path.
        if radial_center <= crate::KERNEL_EPSILON || !core.independently_tested_profile() {
            return None;
        }
        let first = *core.supports().first()?;
        let sign = orientation(self.evidence, first, carrier, budget, angular_tolerance)?;
        let mut visited = core.supports().iter().copied().collect::<BTreeSet<_>>();
        let mut supports = core.supports().to_vec();
        let mut pending = VecDeque::from(supports.clone());
        while let Some(id) = pending.pop_front() {
            for &interface in &self.evidence.support_interface_indices[id] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let other = if id == a { b } else { a };
                if visited.insert(other)
                    && orientation(self.evidence, other, carrier, budget, angular_tolerance)
                        == Some(sign)
                {
                    supports.push(other);
                    pending.push_back(other);
                }
            }
        }
        if supports.len() == core.supports().len() {
            return None;
        }
        supports.sort_unstable();
        self.measure_surface_region(&supports, carrier, budget)
            .map(SurfaceRegionProposal::with_independent_profile)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Point3,
        mesh_recovery::{EvidenceOptions, analyze_mesh_evidence},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    #[test]
    fn precise_toroid_core_extends_to_irregular_source_strip()
    -> Result<(), Box<dyn std::error::Error>> {
        let n = 16;
        let mut positions = Vec::new();
        for row in 0..5 {
            for j in 0..n {
                let a = std::f64::consts::TAU * crate::usize_to_f64(j) / crate::usize_to_f64(n);
                let b = 0.2 + 0.2 * crate::usize_to_f64(row);
                // The final source strip is not a repeated meridian ring
                // and has a small displacement inside the declared budget.
                let b = if row == 4 { b + 0.04 * a.sin() } else { b };
                let radius: f64 = if row == 4 { 0.500001 } else { 0.5 };
                let (sa, ca) = a.sin_cos();
                let (sb, cb) = b.sin_cos();
                let rho = radius.mul_add(cb, 2.0);
                positions.push(Point3::new(rho * ca, rho * sa, radius * sb));
            }
        }
        let faces = (0..4).flat_map(|row| {
            (0..n).flat_map(move |j| {
                let a = row * n + j;
                let b = row * n + (j + 1) % n;
                [[a, b, b + n], [a, b + n, a + n]]
            })
        });
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let mut supports = (0..6 * n)
            .map(|p| evidence.source_planar_support[p])
            .collect::<Vec<_>>();
        supports.sort_unstable();
        supports.dedup();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: 2.0,
            profile_radius: 0.5,
        };
        let budget = SurfaceRegionBudget {
            vertex: 1.0e-5,
            chord: 0.15,
        };
        let core = index
            .measure_surface_region(&supports, carrier, budget)
            .ok_or("core")?;
        assert!(core.independently_tested_profile());
        let expanded = index
            .expand_circular_region(&core, budget, 1.0e-5)
            .ok_or("expanded")?;
        assert_eq!(expanded.source_polygons(), 8 * n);
        assert!(
            index
                .expand_circular_region(
                    &core,
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        ..budget
                    },
                    1.0e-5
                )
                .is_none()
        );
        Ok(())
    }
}
