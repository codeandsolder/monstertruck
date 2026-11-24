//! Distinguish common trim containment from common surface identity.
//! Compare carriers in a source-local domain, not at arbitrary world origins.

use super::SurfaceCarrier;
use crate::{InnerSpace, Point3, Vector3};

/// Difference between two same-family carrier residual functions over a
/// source-local ball.
///
/// A small bound establishes geometric compatibility, not
/// permission to remove topology (e.g. material or periodic seams).
#[derive(Clone, Copy, Debug)]
pub struct CarrierCoincidenceEvidence {
    /// Undirected normal/axis angle in radians.
    pub direction_angle_radians: f64,
    /// Plane offset difference, or distance between the closest points on
    /// both axes to the same reference point.
    pub reference_position_difference: f64,
    /// Profile difference at the reference: expected radius difference for
    /// linear profiles; meridian center displacement plus circle-radius
    /// difference for circular profiles. Zero for planes.
    pub reference_profile_difference: f64,
    /// Additional axial-frame variation per unit length: difference of
    /// radius gradient vectors for linear profiles; difference of aligned
    /// axial directions for circular profiles. Zero for planes.
    pub profile_gradient_difference: f64,
    /// Bound on the absolute difference of the two residual functions
    /// anywhere within the requested ball. Plane normals are sign-aligned;
    /// rotational residuals are signed distances to the meridian profile.
    /// This is a numerical analytic bound, not interval arithmetic.
    pub residual_difference_upper_bound: f64,
}

impl SurfaceCarrier {
    /// Compare independent carrier hypotheses over a ball enclosing the
    /// observed source domains. Infinite carrier identity should not depend
    /// on the arbitrary choice of axial origin or the sign of an axis.
    ///
    /// Returns None for differing carrier families or invalid input.
    /// For same-family carriers this measures a bound rather than deciding
    /// whether the domains, orientations, or topology permit a merge.
    #[must_use = "inspect geometric agreement before considering region unions"]
    pub fn compare_in_ball(
        self,
        other: Self,
        reference: Point3,
        extent: f64,
    ) -> Option<CarrierCoincidenceEvidence> {
        if !extent.is_finite() || extent < 0.0 || !finite_point(reference) {
            return None;
        }
        match (self, other) {
            (
                Self::Plane {
                    origin: first,
                    normal: first_normal,
                },
                Self::Plane {
                    origin: second,
                    normal: second_normal,
                },
            ) if super::containment::unit(first_normal)
                && super::containment::unit(second_normal)
                && finite_point(first)
                && finite_point(second) =>
            {
                let (angle, aligned_normal) = aligned_direction(first_normal, second_normal);
                let offset = ((reference - first).dot(first_normal)
                    - (reference - second).dot(aligned_normal))
                .abs();
                let normal_difference = (first_normal - aligned_normal).magnitude();
                Some(CarrierCoincidenceEvidence {
                    direction_angle_radians: angle,
                    reference_position_difference: offset,
                    reference_profile_difference: 0.0,
                    profile_gradient_difference: 0.0,
                    residual_difference_upper_bound: extent.mul_add(normal_difference, offset),
                })
            }
            (
                Self::LinearRevolution {
                    origin: first,
                    axis: first_axis,
                    radius_at_origin: first_radius,
                    slope: first_slope,
                },
                Self::LinearRevolution {
                    origin: second,
                    axis: second_axis,
                    radius_at_origin: second_radius,
                    slope: second_slope,
                },
            ) if super::containment::unit(first_axis)
                && super::containment::unit(second_axis)
                && finite_point(first)
                && finite_point(second)
                && [first_radius, second_radius, first_slope, second_slope]
                    .into_iter()
                    .all(f64::is_finite) =>
            {
                let (angle, _) = aligned_direction(first_axis, second_axis);
                let first_delta = reference - first;
                let second_delta = reference - second;
                let first_axial = first_delta.dot(first_axis);
                let second_axial = second_delta.dot(second_axis);
                let first_offset = first_delta - first_axis * first_axial;
                let second_offset = second_delta - second_axis * second_axial;
                let position = (first_offset - second_offset).magnitude();
                let radius = (first_slope.mul_add(first_axial, first_radius)
                    - second_slope.mul_add(second_axial, second_radius))
                .abs();
                let gradient = (first_axis * first_slope - second_axis * second_slope).magnitude();
                // ||P1-P2||_2 = sin(angle) for orthogonal projectors
                // transverse to unit axes. Apply the reverse triangle
                // inequality to their radial norms and add the affine
                // profile difference over the same domain.
                Some(CarrierCoincidenceEvidence {
                    direction_angle_radians: angle,
                    reference_position_difference: position,
                    reference_profile_difference: radius,
                    profile_gradient_difference: gradient,
                    residual_difference_upper_bound: extent
                        .mul_add(angle.sin() + gradient, position + radius),
                })
            }
            (
                Self::CircularRevolution {
                    axial_center: first,
                    axis: first_axis,
                    radial_center: first_radial,
                    profile_radius: first_radius,
                },
                Self::CircularRevolution {
                    axial_center: second,
                    axis: second_axis,
                    radial_center: second_radial,
                    profile_radius: second_radius,
                },
            ) => compare_circular(
                [first, second],
                [first_axis, second_axis],
                [first_radial, second_radial],
                [first_radius, second_radius],
                reference,
                extent,
            ),
            _ => None,
        }
        .filter(|fit| fit.residual_difference_upper_bound.is_finite())
    }
}

fn compare_circular(
    [first, second]: [Point3; 2],
    [first_axis, second_axis]: [Vector3; 2],
    [first_radial, second_radial]: [f64; 2],
    [first_radius, second_radius]: [f64; 2],
    reference: Point3,
    extent: f64,
) -> Option<CarrierCoincidenceEvidence> {
    if !super::containment::unit(first_axis)
        || !super::containment::unit(second_axis)
        || !finite_point(first)
        || !finite_point(second)
        || ![first_radial, second_radial, first_radius, second_radius]
            .into_iter()
            .all(f64::is_finite)
        || first_radius <= 0.0
        || second_radius <= 0.0
    {
        return None;
    }
    if first_radial == 0.0 && second_radial == 0.0 {
        // Spheres have no identifiable rotational axis.
        let position = (first - second).magnitude();
        let profile = (first_radius - second_radius).abs();
        return Some(CarrierCoincidenceEvidence {
            direction_angle_radians: 0.0,
            reference_position_difference: position,
            reference_profile_difference: profile,
            profile_gradient_difference: 0.0,
            residual_difference_upper_bound: position + profile,
        });
    }
    let (angle, aligned_axis) = aligned_direction(first_axis, second_axis);
    let first_delta = reference - first;
    let second_delta = reference - second;
    let first_axial = first_delta.dot(first_axis);
    let second_axial = second_delta.dot(aligned_axis);
    let first_offset = first_delta - first_axis * first_axial;
    let second_offset = second_delta - aligned_axis * second_axial;
    let position = (first_offset - second_offset).magnitude();
    let profile = (first_axial - second_axial).hypot(first_radial - second_radial)
        + (first_radius - second_radius).abs();
    let gradient = (first_axis - aligned_axis).magnitude();
    Some(CarrierCoincidenceEvidence {
        direction_angle_radians: angle,
        reference_position_difference: position,
        reference_profile_difference: profile,
        profile_gradient_difference: gradient,
        residual_difference_upper_bound: extent.mul_add(angle.sin() + gradient, position + profile),
    })
}

fn aligned_direction(first: Vector3, second: Vector3) -> (f64, Vector3) {
    let dot = first.dot(second);
    let angle = first.cross(second).magnitude().atan2(dot.abs());
    (angle, if dot < 0.0 { -second } else { second })
}

fn finite_point(point: Point3) -> bool {
    [point.x, point.y, point.z].into_iter().all(f64::is_finite)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cone_origin_translation_and_axis_reversal_preserve_identity()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(100.0, -200.0, 50.0);
        let axis = Vector3::unit_z();
        let first = SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: 3.0,
            slope: 0.25,
        };
        let second = SurfaceCarrier::LinearRevolution {
            origin: origin + axis * 8.0,
            axis: -axis,
            radius_at_origin: 5.0,
            slope: -0.25,
        };
        let fit = first
            .compare_in_ball(second, origin + Vector3::new(1.0, 2.0, 3.0), 10.0)
            .ok_or("missing comparison")?;
        assert!(fit.residual_difference_upper_bound < 1.0e-12);
        Ok(())
    }

    #[test]
    fn same_circle_does_not_make_different_cone_slopes_coincident()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let axis = Vector3::unit_z();
        let first = SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: 2.0,
            slope: 0.0,
        };
        let second = SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: 2.0,
            slope: 0.5,
        };
        let fit = first
            .compare_in_ball(second, origin, 4.0)
            .ok_or("missing comparison")?;
        assert!(fit.reference_profile_difference < 1.0e-12);
        assert!((fit.residual_difference_upper_bound - 2.0).abs() < 1.0e-12);
        let plane = SurfaceCarrier::Plane {
            origin,
            normal: axis,
        };
        assert!(first.compare_in_ball(plane, origin, 4.0).is_none());
        Ok(())
    }

    #[test]
    fn tilted_axis_comparison_bounds_both_residual_functions()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(2.0, 3.0, 4.0);
        let other_origin = origin + Vector3::new(0.2, 0.3, -1.0);
        let first_axis = Vector3::unit_z();
        let other_axis = Vector3::new(0.1, -0.2, 1.0).normalize();
        let first = SurfaceCarrier::LinearRevolution {
            origin,
            axis: first_axis,
            radius_at_origin: 3.0,
            slope: 0.25,
        };
        let second = SurfaceCarrier::LinearRevolution {
            origin: other_origin,
            axis: other_axis,
            radius_at_origin: 2.0,
            slope: 0.4,
        };
        let fit = first
            .compare_in_ball(second, origin, 5.0)
            .ok_or("missing comparison")?;
        for x in -2..=2 {
            for y in -2..=2 {
                for z in -2..=2 {
                    let point = origin + Vector3::new(f64::from(x), f64::from(y), f64::from(z));
                    let a = radial_error(point - origin, first_axis, 3.0, 0.25);
                    let b = radial_error(point - other_origin, other_axis, 2.0, 0.4);
                    assert!((a - b).abs() <= fit.residual_difference_upper_bound);
                }
            }
        }
        Ok(())
    }

    fn radial_error(delta: Vector3, axis: Vector3, radius: f64, slope: f64) -> f64 {
        let axial = delta.dot(axis);
        (delta - axis * axial).magnitude() - slope.mul_add(axial, radius)
    }
}
