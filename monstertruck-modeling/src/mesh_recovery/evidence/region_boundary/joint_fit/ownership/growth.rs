//! Source evidence for adding a complete planar support at a failed curved
//! contact. The carrier is fixed and its endpoint errors do not establish
//! the uncertainty of the input.
use super::{MeshEvidence, Point3, SurfaceCarrier, chord_uncertainty};
use crate::{InnerProduct, KERNEL_EPSILON};

#[derive(Clone, Copy)]
pub(in crate::mesh_recovery) struct OwnershipGrowthOptions {
    pub maximum_adjustment: f64,
    pub declared_uncertainty: f64,
    pub angular_tolerance: f64,
    pub linear_chord_fraction: f64,
    pub circular_chord_fraction: f64,
}

fn polygon_points(evidence: &MeshEvidence, polygon: usize) -> Option<Vec<Point3>> {
    evidence
        .source_polygon_vertices(polygon)?
        .iter()
        .map(|&id| evidence.positions.get(id).copied())
        .collect()
}

fn center(points: &[Point3]) -> Option<Point3> {
    let first = *points.first()?;
    Some(
        first
            + points.iter().map(|&p| p - first).sum::<crate::Vector3>()
                / crate::usize_to_f64(points.len()),
    )
}

fn oriented_sign(evidence: &MeshEvidence, support: usize, carrier: SurfaceCarrier) -> Option<bool> {
    let plane = evidence.planar_supports.get(support)?;
    let points = polygon_points(evidence, *plane.source_face_indices.first()?)?;
    let normal = carrier.signed_constraint(center(&points)?)?.0;
    Some(normal.dot(&plane.normal) > 0.0)
}

fn vertex_error(
    carrier: SurfaceCarrier,
    points: &[Point3],
    policy: OwnershipGrowthOptions,
) -> Option<f64> {
    let mut maximum = 0.0_f64;
    for (i, &p) in points.iter().enumerate() {
        let before = points[(i + points.len() - 1) % points.len()];
        let after = points[(i + 1) % points.len()];
        let numerical =
            KERNEL_EPSILON.max(32.0 * f64::EPSILON * p.x.abs().max(p.y.abs()).max(p.z.abs()));
        let allowance = chord_uncertainty(carrier, p, before)
            .max(chord_uncertainty(carrier, p, after))
            .max(policy.declared_uncertainty)
            .max(numerical)
            .min(policy.maximum_adjustment);
        let residual = carrier.signed_constraint(p).map_or_else(
            || {
                let radial = carrier.measure_segment(p, p)?.witnessed_max;
                Some(match carrier {
                    SurfaceCarrier::LinearRevolution { slope, .. } => radial / 1.0_f64.hypot(slope),
                    _ => radial,
                })
            },
            |(_, value)| Some(value.abs()),
        )?;
        if residual > allowance {
            return None;
        }
        maximum = maximum.max(residual);
    }
    Some(maximum)
}

fn has_meridian_extent(carrier: SurfaceCarrier, points: &[Point3], angular: f64) -> Option<bool> {
    let (origin, axis, circular) = match carrier {
        SurfaceCarrier::Plane { .. } => return None,
        SurfaceCarrier::LinearRevolution { origin, axis, .. } => (origin, axis, false),
        SurfaceCarrier::CircularRevolution {
            axial_center, axis, ..
        } => (axial_center, axis, true),
    };
    let first = *points.first()?;
    let meridian = |p: Point3| {
        let d = p - origin;
        let z = d.dot(&axis);
        (z, (d - axis * z).norm())
    };
    let reference = meridian(first);
    let span = points
        .iter()
        .map(|&p| {
            let (z, r) = meridian(p);
            if circular {
                (z - reference.0).hypot(r - reference.1)
            } else {
                (z - reference.0).abs()
            }
        })
        .fold(0.0, f64::max);
    let source = points
        .iter()
        .map(|&p| (p - first).norm())
        .fold(0.0, f64::max);
    Some(span > angular * source)
}

fn normal_agrees(
    evidence: &MeshEvidence,
    support: usize,
    carrier: SurfaceCarrier,
    points: &[Point3],
    sign: bool,
    angular: f64,
) -> Option<bool> {
    let plane = evidence.planar_supports.get(support)?;
    let expected = carrier.signed_constraint(center(points)?)?.0;
    let dot = expected.dot(&plane.normal);
    let deviation = expected.cross(&plane.normal).norm().atan2(dot.abs());
    let span = points
        .iter()
        .filter_map(|&p| {
            carrier
                .signed_constraint(p)
                .map(|(n, _)| expected.cross(&n).norm().atan2(expected.dot(&n)))
        })
        .fold(0.0, f64::max);
    Some(
        dot.abs() > angular.sin()
            && (dot > 0.0) == sign
            && deviation <= span + angular + plane.max_normal_angle_radians,
    )
}

/// The smaller residual is the better-supported fixed carrier. Equal
/// residuals remain ambiguous, rather than granting an arbitrary owner.
pub(in crate::mesh_recovery) fn support_claim(
    evidence: &MeshEvidence,
    support: usize,
    curved_neighbor: usize,
    carrier: SurfaceCarrier,
    policy: OwnershipGrowthOptions,
) -> Option<f64> {
    let sign = oriented_sign(evidence, curved_neighbor, carrier)?;
    let mut maximum = 0.0_f64;
    for &polygon in &evidence.planar_supports.get(support)?.source_face_indices {
        let points = polygon_points(evidence, polygon)?;
        if !has_meridian_extent(carrier, &points, policy.angular_tolerance)?
            || !normal_agrees(
                evidence,
                support,
                carrier,
                &points,
                sign,
                policy.angular_tolerance,
            )?
        {
            return None;
        }
        maximum = maximum.max(vertex_error(carrier, &points, policy)?);
    }
    match carrier {
        SurfaceCarrier::Plane { .. } => return None,
        SurfaceCarrier::LinearRevolution { slope, .. } => {
            // The local radial filled-facet certificate uses unnormalized
            // meridian residual, unlike the signed distance witness.
            super::super::super::rim_linear::facet_bound(
                evidence,
                support,
                carrier,
                policy.maximum_adjustment * 1.0_f64.hypot(slope),
                policy.linear_chord_fraction,
            )?;
        }
        SurfaceCarrier::CircularRevolution { profile_radius, .. } => {
            let budget = policy.circular_chord_fraction * profile_radius;
            if carrier
                .measure_circular_source_facets(evidence, &[support], budget)?
                .upper_bound
                > budget
            {
                return None;
            }
        }
    }
    Some(maximum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Vector3,
        mesh_recovery::{EvidenceOptions, analyze_mesh_evidence},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    #[test]
    fn noisy_lateral_facet_passes_but_planar_cap_and_wrong_carrier_do_not()
    -> Result<(), Box<dyn std::error::Error>> {
        let angle = std::f64::consts::FRAC_PI_4;
        let (sin, cos) = angle.sin_cos();
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(cos, sin, 0.0),
                    Point3::new(cos, sin, 1.0),
                    Point3::new(1.0002, 0.0, 1.0),
                    Point3::new(0.0, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [0, 2, 3], [1, 0, 4]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let core = evidence.source_planar_support[0];
        let side = evidence.source_planar_support[1];
        let cap = evidence.source_planar_support[2];
        let carrier = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let policy = OwnershipGrowthOptions {
            maximum_adjustment: 0.01,
            declared_uncertainty: 0.0,
            angular_tolerance: 1.0e-5,
            linear_chord_fraction: 0.25,
            circular_chord_fraction: 0.25,
        };
        assert!(support_claim(&evidence, side, core, carrier, policy).is_some());
        assert!(support_claim(&evidence, cap, core, carrier, policy).is_none());
        assert!(
            support_claim(
                &evidence,
                side,
                core,
                carrier,
                OwnershipGrowthOptions {
                    maximum_adjustment: 1.0e-6,
                    ..policy
                }
            )
            .is_none()
        );
        let wrong = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius_at_origin: 0.8,
            slope: 0.0,
        };
        assert!(support_claim(&evidence, side, core, wrong, policy).is_none());
        Ok(())
    }
}
