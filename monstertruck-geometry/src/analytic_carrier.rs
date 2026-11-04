//! Kernel-native, parameterization-independent analytic surface identity.
//!
//! This module deliberately separates equality of a complete *carrier* from
//! equality of trimmed faces or occupied solids. A surface certificate does not
//! prove anything about its trims, its bounded domain or its adjacency graph.
//! Those claims belong to the later boundary/topology equivalence proof.

use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

/// Unoriented analytic support in 3D physical coordinates.
///
/// Axis directions are canonicalized. These values are only created by the
/// checked constructors below, not by trusting a sampled surface fit.
#[derive(Clone, Copy, Debug)]
pub enum AnalyticSupport {
    Plane {
        origin: Point3,
        normal: Vector3,
    },
    Sphere {
        center: Point3,
        radius: f64,
    },
    Cylinder {
        origin: Point3,
        axis: Vector3,
        radius: f64,
    },
    Cone {
        apex: Point3,
        axis: Vector3,
        slope: f64,
        /// A cone's two nappes have the same implicit equation but are
        /// different physical surface carriers.
        positive_nappe: bool,
    },
    Torus {
        center: Point3,
        axis: Vector3,
        major_radius: f64,
        minor_radius: f64,
    },
}

/// A structurally certified carrier and the global orientation of its chart.
#[derive(Clone, Copy, Debug)]
pub struct AnalyticCarrier {
    support: AnalyticSupport,
    /// Global chart orientation relative to the canonical unoriented support.
    orientation: bool,
}

/// The orientation relation between two identical analytic supports.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CarrierDomainRelation {
    SameOrientation,
    ReversedOrientation,
}

/// A finite, nonzero direction normalized with a deterministic sign.
#[must_use]
pub fn canonical_direction(mut vector: Vector3) -> Option<Vector3> {
    let magnitude = vector.magnitude();
    if !magnitude.is_finite() || magnitude <= KERNEL_EPSILON {
        return None;
    }
    vector /= magnitude;
    let flip = [vector.x, vector.y, vector.z]
        .into_iter()
        .find(|component| component.abs() > 64.0 * f64::EPSILON)
        .is_some_and(|component| component < 0.0);
    if flip {
        vector = -vector;
    }
    Some(vector)
}

/// Validate that a homogeneous transform is a finite uniform similarity.
///
/// Reject shear, perspective, singular transforms and non-finite coefficients.
#[must_use]
pub fn similarity_scale(transform: Matrix4) -> Option<f64> {
    let tolerance = DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity;
    if ![
        transform.x.x,
        transform.x.y,
        transform.x.z,
        transform.x.w,
        transform.y.x,
        transform.y.y,
        transform.y.z,
        transform.y.w,
        transform.z.x,
        transform.z.y,
        transform.z.z,
        transform.z.w,
        transform.w.x,
        transform.w.y,
        transform.w.z,
        transform.w.w,
    ]
    .into_iter()
    .all(f64::is_finite)
        || transform.x.w.abs() > tolerance
        || transform.y.w.abs() > tolerance
        || transform.z.w.abs() > tolerance
        || (transform.w.w - 1.0).abs() > tolerance
    {
        return None;
    }
    let columns = [
        transform.x.truncate(),
        transform.y.truncate(),
        transform.z.truncate(),
    ];
    let scale = columns[0].magnitude();
    if !scale.is_finite() || scale <= KERNEL_EPSILON {
        return None;
    }
    let length_tolerance = tolerance * scale.max(1.0);
    let dot_tolerance = tolerance * scale * scale;
    if columns[1..].iter().any(|column| {
        !column.magnitude().is_finite() || (column.magnitude() - scale).abs() > length_tolerance
    }) || [(0, 1), (0, 2), (1, 2)]
        .into_iter()
        .any(|(left, right)| columns[left].dot(columns[right]).abs() > dot_tolerance)
    {
        return None;
    }
    Some(scale)
}

fn processor_orientation(orientation: bool, transform: Matrix4) -> Option<f64> {
    similarity_scale(transform)?;
    let determinant = transform.determinant();
    if !determinant.is_finite() || determinant.abs() <= KERNEL_EPSILON {
        return None;
    }
    Some(if orientation {
        determinant.signum()
    } else {
        -determinant.signum()
    })
}

fn finite_point(point: Point3) -> bool {
    [point.x, point.y, point.z].into_iter().all(f64::is_finite)
}

impl AnalyticCarrier {
    /// The checked, unoriented geometric support (never a trim or solid proof).
    #[must_use]
    pub const fn support(self) -> AnalyticSupport { self.support }

    /// Certify a planar carrier with its parametric orientation.
    #[must_use]
    pub fn plane(plane: &Plane) -> Option<Self> {
        let normal = canonical_direction(plane.normal())?;
        let origin = plane.origin();
        if !finite_point(origin) {
            return None;
        }
        let orientation = plane.normal().dot(normal);
        if !orientation.is_finite() || orientation.abs() <= KERNEL_EPSILON {
            return None;
        }
        Some(Self {
            support: AnalyticSupport::Plane { origin, normal },
            orientation: orientation > 0.0,
        })
    }

    /// Certify a sphere under a true similarity transform.
    #[must_use]
    pub fn sphere(sphere: &Processor<Sphere, Matrix4>) -> Option<Self> {
        Self::sphere_pose(
            sphere.entity().center(),
            sphere.entity().radius(),
            *sphere.transform(),
            sphere.orientation(),
        )
    }

    /// Certify a sphere supplied by an exchange adapter without discarding its
    /// original chart transform or orientation.
    #[must_use]
    pub fn sphere_pose(
        center: Point3,
        radius: f64,
        transform: Matrix4,
        chart_orientation: bool,
    ) -> Option<Self> {
        let scale = similarity_scale(transform)?;
        let center = transform.transform_point(center);
        let radius = radius * scale;
        if !finite_point(center) || !radius.is_finite() || radius <= KERNEL_EPSILON {
            return None;
        }
        Some(Self {
            support: AnalyticSupport::Sphere { center, radius },
            orientation: processor_orientation(chart_orientation, transform)? > 0.0,
        })
    }

    /// Certify a torus under a true similarity transform.
    ///
    /// Self-intersecting spindle tori have no one global implicit orientation.
    #[must_use]
    pub fn torus(torus: &Processor<Torus, Matrix4>) -> Option<Self> {
        let transform = *torus.transform();
        let scale = similarity_scale(transform)?;
        let center = transform.transform_point(torus.entity().center());
        let axis = canonical_direction(transform.transform_vector(Vector3::unit_z()))?;
        let major_radius = torus.entity().large_radius() * scale;
        let minor_radius = torus.entity().small_radius() * scale;
        if !finite_point(center)
            || !major_radius.is_finite()
            || !minor_radius.is_finite()
            || minor_radius <= KERNEL_EPSILON
            || major_radius <= KERNEL_EPSILON
            || major_radius + KERNEL_EPSILON < minor_radius
        {
            return None;
        }
        Some(Self {
            support: AnalyticSupport::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            },
            orientation: processor_orientation(torus.orientation(), transform)? > 0.0,
        })
    }

    /// Certify a full linear surface of revolution as a cylinder or a cone.
    ///
    /// The generic revolution model may contain a line, spline or composite
    /// curve. Only a *structural line* permits this proof. `profile`, `origin`
    /// and `axis` are in the untransformed revolution coordinate system.
    #[must_use]
    pub fn linear_revolution(
        profile: Line<Point3>,
        origin: Point3,
        axis: Vector3,
        transform: Matrix4,
        chart_orientation: bool,
    ) -> Option<Self> {
        similarity_scale(transform)?;
        // The revolution chart's angular sense uses the ORIGINAL axis
        // direction, not the sign-canonicalized support axis.
        let axis_magnitude = axis.magnitude();
        if !axis_magnitude.is_finite() || axis_magnitude <= KERNEL_EPSILON {
            return None;
        }
        let local_axis = axis / axis_magnitude;
        let Line(local_front, local_back) = profile;
        let axial_step = (local_back - local_front).dot(local_axis);
        if !axial_step.is_finite() || axial_step.abs() <= KERNEL_EPSILON {
            return None;
        }
        let orientation =
            -axial_step.signum() * processor_orientation(chart_orientation, transform)?;
        let origin = transform.transform_point(origin);
        let axis = canonical_direction(transform.transform_vector(axis))?;
        let front = transform.transform_point(local_front);
        let back = transform.transform_point(local_back);
        if !finite_point(origin) || !finite_point(front) || !finite_point(back) {
            return None;
        }

        let direction = back - front;
        let radial0 = (front - origin) - axis * (front - origin).dot(axis);
        let radial_direction = direction - axis * direction.dot(axis);
        let scale = direction.magnitude().max(radial0.magnitude()).max(1.0);
        let support = if radial_direction.magnitude()
            <= DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity * scale
        {
            let radius = radial0.magnitude();
            if !radius.is_finite() || radius <= KERNEL_EPSILON {
                return None;
            }
            AnalyticSupport::Cylinder {
                origin,
                axis,
                radius,
            }
        } else {
            let radial_direction2 = radial_direction.magnitude2();
            if !radial_direction2.is_finite() || radial_direction2 <= KERNEL_EPSILON.powi(2) {
                return None;
            }
            let parameter = -radial0.dot(radial_direction) / radial_direction2;
            let apex = front + direction * parameter;
            let apex_radial = (apex - origin) - axis * (apex - origin).dot(axis);
            let scale = direction
                .magnitude()
                .max((apex - origin).magnitude())
                .max(1.0);
            if !finite_point(apex)
                || apex_radial.magnitude()
                    > DEFAULT_NUMERICAL_TOLERANCES.transform_similarity * scale
            {
                return None;
            }
            let axial = direction.dot(axis).abs();
            if !axial.is_finite() || axial <= KERNEL_EPSILON {
                return None;
            }
            let slope = radial_direction.magnitude() / axial;
            if !slope.is_finite() {
                return None;
            }
            let front_axial = (front - apex).dot(axis);
            let back_axial = (back - apex).dot(axis);
            let tolerance = DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity * scale;
            if !front_axial.is_finite()
                || !back_axial.is_finite()
                || (front_axial > tolerance && back_axial < -tolerance)
                || (front_axial < -tolerance && back_axial > tolerance)
                || (front_axial.abs() <= tolerance && back_axial.abs() <= tolerance)
            {
                // A profile crossing the apex reverses its parametric
                // orientation; the whole curve is not one regular nappe.
                return None;
            }
            let positive_nappe = if front_axial.abs() > tolerance {
                front_axial > 0.0
            } else {
                back_axial > 0.0
            };
            AnalyticSupport::Cone {
                apex,
                axis,
                slope,
                positive_nappe,
            }
        };
        Some(Self {
            support,
            orientation: orientation > 0.0,
        })
    }

    /// Compare two certified analytic supports and their *global* orientation.
    ///
    /// This is not a claim that trimmed boundaries or solids are equal.
    #[must_use]
    pub fn same_domain(
        self,
        other: Self,
        linear_tolerance: f64,
        angular_tolerance: f64,
    ) -> Option<CarrierDomainRelation> {
        if !same_support(
            self.support,
            other.support,
            linear_tolerance,
            angular_tolerance,
        ) {
            return None;
        }
        Some(if self.orientation == other.orientation {
            CarrierDomainRelation::SameOrientation
        } else {
            CarrierDomainRelation::ReversedOrientation
        })
    }
}

fn parallel(left: Vector3, right: Vector3, angular_tolerance: f64) -> bool {
    left.cross(right).magnitude() <= angular_tolerance.sin().abs().max(64.0 * f64::EPSILON)
}

fn axis_line_distance(origin: Point3, axis: Vector3, other: Point3) -> f64 {
    (other - origin).cross(axis).magnitude()
}

/// Compare certified unoriented analytic supports, independently of charts.
#[must_use]
pub fn same_support(
    left: AnalyticSupport,
    right: AnalyticSupport,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> bool {
    if !linear_tolerance.is_finite()
        || linear_tolerance < 0.0
        || !angular_tolerance.is_finite()
        || angular_tolerance < 0.0
    {
        return false;
    }
    match (left, right) {
        (
            AnalyticSupport::Plane {
                origin: o0,
                normal: n0,
            },
            AnalyticSupport::Plane {
                origin: o1,
                normal: n1,
            },
        ) => parallel(n0, n1, angular_tolerance) && (o1 - o0).dot(n0).abs() <= linear_tolerance,
        (
            AnalyticSupport::Sphere {
                center: c0,
                radius: r0,
            },
            AnalyticSupport::Sphere {
                center: c1,
                radius: r1,
            },
        ) => c0.distance(c1) <= linear_tolerance && (r0 - r1).abs() <= linear_tolerance,
        (
            AnalyticSupport::Cylinder {
                origin: o0,
                axis: a0,
                radius: r0,
            },
            AnalyticSupport::Cylinder {
                origin: o1,
                axis: a1,
                radius: r1,
            },
        ) => {
            parallel(a0, a1, angular_tolerance)
                && axis_line_distance(o0, a0, o1) <= linear_tolerance
                && (r0 - r1).abs() <= linear_tolerance
        }
        (
            AnalyticSupport::Cone {
                apex: p0,
                axis: a0,
                slope: s0,
                positive_nappe: n0,
            },
            AnalyticSupport::Cone {
                apex: p1,
                axis: a1,
                slope: s1,
                positive_nappe: n1,
            },
        ) => {
            n0 == n1
                && p0.distance(p1) <= linear_tolerance
                && parallel(a0, a1, angular_tolerance)
                && (s0.atan() - s1.atan()).abs() <= angular_tolerance
        }
        (
            AnalyticSupport::Torus {
                center: c0,
                axis: a0,
                major_radius: r0,
                minor_radius: m0,
            },
            AnalyticSupport::Torus {
                center: c1,
                axis: a1,
                major_radius: r1,
                minor_radius: m1,
            },
        ) => {
            c0.distance(c1) <= linear_tolerance
                && parallel(a0, a1, angular_tolerance)
                && (r0 - r1).abs() <= linear_tolerance
                && (m0 - m1).abs() <= linear_tolerance
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const fn point(x: f64, y: f64, z: f64) -> Point3 { Point3::new(x, y, z) }

    fn revolution(first: Point3, second: Point3, axis: Vector3) -> Option<AnalyticCarrier> {
        AnalyticCarrier::linear_revolution(
            Line(first, second),
            point(0.0, 0.0, 0.0),
            axis,
            Matrix4::from_scale(1.0),
            true,
        )
    }

    #[test]
    fn plane_chart_reversal_changes_orientation_not_support() -> anyhow::Result<()> {
        let forward =
            AnalyticCarrier::plane(&Plane::xy()).ok_or_else(|| anyhow::anyhow!("valid plane"))?;
        let reversed = AnalyticCarrier::plane(&Plane::new(
            point(0.0, 0.0, 0.0),
            point(0.0, 1.0, 0.0),
            point(1.0, 0.0, 0.0),
        ))
        .ok_or_else(|| anyhow::anyhow!("valid plane"))?;
        assert_eq!(
            forward.same_domain(reversed, 1.0e-9, 1.0e-9),
            Some(CarrierDomainRelation::ReversedOrientation),
        );
        assert_eq!(
            forward.same_domain(forward, 1.0e-9, 1.0e-9),
            Some(CarrierDomainRelation::SameOrientation)
        );
        Ok(())
    }

    #[test]
    fn reversed_revolution_axis_changes_global_chart_orientation() -> anyhow::Result<()> {
        let forward = revolution(
            point(2.0, 0.0, 0.0),
            point(2.0, 0.0, 1.0),
            Vector3::unit_z(),
        )
        .ok_or_else(|| anyhow::anyhow!("cylinder"))?;
        let reversed_axis = revolution(
            point(2.0, 0.0, 0.0),
            point(2.0, 0.0, 1.0),
            -Vector3::unit_z(),
        )
        .ok_or_else(|| anyhow::anyhow!("cylinder with reversed axis"))?;
        assert_eq!(
            forward.same_domain(reversed_axis, 1.0e-9, 1.0e-9),
            Some(CarrierDomainRelation::ReversedOrientation),
        );
        Ok(())
    }

    #[test]
    fn cone_nappes_are_distinct_even_with_the_same_implicit_equation() -> anyhow::Result<()> {
        let upper = revolution(
            point(1.0, 0.0, 1.0),
            point(2.0, 0.0, 2.0),
            Vector3::unit_z(),
        )
        .ok_or_else(|| anyhow::anyhow!("upper cone"))?;
        let lower = revolution(
            point(1.0, 0.0, -1.0),
            point(2.0, 0.0, -2.0),
            Vector3::unit_z(),
        )
        .ok_or_else(|| anyhow::anyhow!("lower cone"))?;
        assert!(matches!(
            upper.support(),
            AnalyticSupport::Cone {
                positive_nappe: true,
                ..
            }
        ));
        assert!(matches!(
            lower.support(),
            AnalyticSupport::Cone {
                positive_nappe: false,
                ..
            }
        ));
        assert_eq!(upper.same_domain(lower, 1.0e-9, 1.0e-9), None);

        // A line spanning the apex cannot have one globally regular chart.
        assert!(
            revolution(
                point(1.0, 0.0, 1.0),
                point(-1.0, 0.0, -1.0),
                Vector3::unit_z(),
            )
            .is_none()
        );
        Ok(())
    }

    #[test]
    fn sheared_and_perspective_spheres_are_refused() -> anyhow::Result<()> {
        let sphere = Sphere::new(point(0.0, 0.0, 0.0), 2.0);
        let shear = Matrix4::from_cols(
            Vector4::new(1.0, 0.0, 0.0, 0.0),
            Vector4::new(0.5, 1.0, 0.0, 0.0),
            Vector4::new(0.0, 0.0, 1.0, 0.0),
            Vector4::new(0.0, 0.0, 0.0, 1.0),
        );
        assert!(AnalyticCarrier::sphere(&Processor::with_transform(sphere, shear)).is_none());
        let mut perspective = Matrix4::from_scale(1.0);
        perspective.z.w = 0.1;
        assert!(AnalyticCarrier::sphere(&Processor::with_transform(sphere, perspective)).is_none());
        Ok(())
    }

    #[test]
    fn reflection_reverses_spherical_orientation_but_not_its_support() -> anyhow::Result<()> {
        let sphere = Sphere::new(point(0.0, 0.0, 0.0), 2.0);
        let direct = AnalyticCarrier::sphere(&Processor::new(sphere))
            .ok_or_else(|| anyhow::anyhow!("sphere"))?;
        let reflected = AnalyticCarrier::sphere(&Processor::with_transform(
            sphere,
            Matrix4::from_nonuniform_scale(-1.0, 1.0, 1.0),
        ))
        .ok_or_else(|| anyhow::anyhow!("reflected sphere"))?;
        assert_eq!(
            direct.same_domain(reflected, 1.0e-9, 1.0e-9),
            Some(CarrierDomainRelation::ReversedOrientation),
        );
        Ok(())
    }

    #[test]
    fn ring_torus_has_global_orientation_but_spindle_does_not() -> anyhow::Result<()> {
        let ring = Processor::new(Torus::new_unchecked(point(0.0, 0.0, 0.0), 3.0, 1.0));
        let spindle = Processor::new(Torus::new_unchecked(point(0.0, 0.0, 0.0), 1.0, 2.0));
        assert!(AnalyticCarrier::torus(&ring).is_some());
        assert!(AnalyticCarrier::torus(&spindle).is_none());
        Ok(())
    }

    #[test]
    fn nonfinite_or_negative_tolerances_cannot_justify_carrier_equality() -> anyhow::Result<()> {
        let plane = AnalyticCarrier::plane(&Plane::xy()).ok_or_else(|| anyhow::anyhow!("plane"))?;
        for (linear, angular) in [
            (f64::NAN, 0.0),
            (0.0, f64::NAN),
            (f64::INFINITY, 0.0),
            (-1.0, 0.0),
            (0.0, -1.0),
        ] {
            assert_eq!(plane.same_domain(plane, linear, angular), None);
        }
        Ok(())
    }
}
