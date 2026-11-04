use super::{ElementarySurface, Surface};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use monstertruck_geometry::prelude::*;
use monstertruck_topology::repair::SurfaceDomainRelation;

use monstertruck_geometry::analytic_carrier::{
    AnalyticCarrier, AnalyticSupport, CarrierDomainRelation, similarity_scale,
};

fn carrier(surface: &Surface) -> Option<AnalyticCarrier> {
    match surface {
        Surface::Elementary(ElementarySurface::Plane(plane)) => AnalyticCarrier::plane(plane),
        Surface::Elementary(ElementarySurface::Sphere(sphere)) => AnalyticCarrier::sphere_pose(
            sphere.entity().0.center(),
            sphere.entity().0.radius(),
            *sphere.transform(),
            sphere.orientation(),
        ),
        Surface::Elementary(ElementarySurface::CylindricalSurface(revolution)) => {
            let proof = AnalyticCarrier::linear_revolution(
                *revolution.entity().entity_curve(),
                revolution.entity().origin(),
                revolution.entity().axis(),
                *revolution.transform(),
                revolution.orientation(),
            )?;
            matches!(proof.support(), AnalyticSupport::Cylinder { .. }).then_some(proof)
        }
        Surface::Elementary(ElementarySurface::ConicalSurface(revolution)) => {
            let proof = AnalyticCarrier::linear_revolution(
                *revolution.entity().entity_curve(),
                revolution.entity().origin(),
                revolution.entity().axis(),
                *revolution.transform(),
                revolution.orientation(),
            )?;
            matches!(proof.support(), AnalyticSupport::Cone { .. }).then_some(proof)
        }
        Surface::Elementary(ElementarySurface::ToroidalSurface(torus)) => {
            AnalyticCarrier::torus(torus)
        }
        // A spline's net can prove coplanarity, but unless its UV Jacobian is
        // certified to keep one sign everywhere it cannot supply a global
        // oriented carrier. Keep that distinction at this boundary.
        Surface::Bspline(_) | Surface::Nurbs(_) | Surface::Sweep(_) => None,
    }
}

#[inline]
fn unwrap_parameter_component(value: f64, reference: f64, period: Option<f64>) -> f64 {
    period
        .filter(|period| period.is_finite() && *period > KERNEL_EPSILON)
        .map_or(value, |period| {
            ((reference - value) / period)
                .round()
                .mul_add(period, value)
        })
}

#[inline]
fn unwrap_parameter_near(surface: &Surface, value: Point2, reference: Point2) -> Point2 {
    Point2::new(
        unwrap_parameter_component(value.x, reference.x, surface.period_u()),
        unwrap_parameter_component(value.y, reference.y, surface.period_v()),
    )
}

#[inline]
const fn processor_parameter(orientation: bool, internal: Point2) -> Point2 {
    if orientation {
        internal
    } else {
        Point2::new(internal.y, internal.x)
    }
}

fn exact_plane_parameter(plane: &Plane, point: Point3) -> Option<Point2> {
    let parameter = plane.parameter(point)?;
    let scale = point
        .to_vec()
        .magnitude()
        .max(plane.origin().to_vec().magnitude())
        .max(1.0);
    (parameter.z.abs() <= KERNEL_EPSILON * scale).then_some(Point2::new(parameter.x, parameter.y))
}

fn exact_linear_revolution_parameter(
    surface: &Processor<RevolutionSurface<Line<Point3>>, Matrix4>,
    point: Point3,
) -> Option<Point2> {
    similarity_scale(*surface.transform())?;
    let inverse = surface.transform().invert()?;
    let local = inverse.transform_point(point);
    let revolution = surface.entity();
    let origin = revolution.origin();
    let axis = revolution.axis();
    let Line(profile0, profile1) = *revolution.entity_curve();
    let direction = profile1 - profile0;
    let axial_step = direction.dot(axis);
    let scale = direction
        .magnitude()
        .max((profile0 - origin).magnitude())
        .max((local - origin).magnitude())
        .max(1.0);
    if axial_step.abs() <= KERNEL_EPSILON * scale {
        return None;
    }

    let axial0 = (profile0 - origin).dot(axis);
    let axial = (local - origin).dot(axis);
    let u = (axial - axial0) / axial_step;
    let profile = profile0 + direction * u;
    let profile_delta = profile - origin;
    let point_delta = local - origin;
    let profile_radial = profile_delta - axis * profile_delta.dot(axis);
    let point_radial = point_delta - axis * point_delta.dot(axis);
    let profile_radius = profile_radial.magnitude();
    let point_radius = point_radial.magnitude();
    let radial_scale = profile_radius.max(point_radius).max(scale);
    if profile_radius <= KERNEL_EPSILON * radial_scale
        || point_radius <= KERNEL_EPSILON * radial_scale
        || (profile_radius - point_radius).abs() > KERNEL_EPSILON * radial_scale
    {
        return None;
    }
    let angle = axis
        .dot(profile_radial.cross(point_radial))
        .atan2(profile_radial.dot(point_radial))
        .rem_euclid(std::f64::consts::TAU);
    let external = processor_parameter(surface.orientation(), Point2::new(u, angle));
    surface
        .evaluate(external.x, external.y)
        .near(&point)
        .then_some(external)
}

fn exact_torus_parameter(surface: &Processor<Torus, Matrix4>, point: Point3) -> Option<Point2> {
    similarity_scale(*surface.transform())?;
    let inverse = surface.transform().invert()?;
    let local = inverse.transform_point(point);
    let torus = surface.entity();
    let relative = local - torus.center();
    let rho = relative.x.hypot(relative.y);
    let radius = torus.small_radius();
    let scale = torus
        .large_radius()
        .max(radius)
        .max(relative.magnitude())
        .max(1.0);
    let sine = relative.z / radius;
    let epsilon = 128.0 * f64::EPSILON * scale;
    if !sine.is_finite() || sine < -1.0 - epsilon || sine > 1.0 + epsilon {
        return None;
    }

    let alpha = sine.clamp(-1.0, 1.0).asin();
    let azimuth = relative.y.atan2(relative.x);
    [
        alpha.rem_euclid(std::f64::consts::TAU),
        (std::f64::consts::PI - alpha).rem_euclid(std::f64::consts::TAU),
    ]
    .into_iter()
    .find_map(|v| {
        let signed_radial = radius.mul_add(v.cos(), torus.large_radius());
        if (signed_radial.abs() - rho).abs() > KERNEL_EPSILON * scale {
            return None;
        }
        let u = if signed_radial >= 0.0 {
            azimuth
        } else {
            azimuth + std::f64::consts::PI
        }
        .rem_euclid(std::f64::consts::TAU);
        let external = processor_parameter(surface.orientation(), Point2::new(u, v));
        surface
            .evaluate(external.x, external.y)
            .near(&point)
            .then_some(external)
    })
}

fn exact_supported_parameter(surface: &Surface, point: Point3) -> Option<Point2> {
    match surface {
        Surface::Elementary(ElementarySurface::Plane(plane)) => exact_plane_parameter(plane, point),
        Surface::Elementary(
            ElementarySurface::CylindricalSurface(revolution)
            | ElementarySurface::ConicalSurface(revolution),
        ) => exact_linear_revolution_parameter(revolution, point),
        Surface::Elementary(ElementarySurface::ToroidalSurface(torus)) => {
            exact_torus_parameter(torus, point)
        }
        Surface::Elementary(ElementarySurface::Sphere(_))
        | Surface::Bspline(_)
        | Surface::Nurbs(_)
        | Surface::Sweep(_) => None,
    }
}

fn safe_revolution_basis(
    surface: &Processor<RevolutionSurface<Line<Point3>>, Matrix4>,
) -> Option<[Point2; 3]> {
    let revolution = surface.entity();
    let origin = revolution.origin();
    let axis = revolution.axis();
    let Line(profile0, profile1) = *revolution.entity_curve();
    let direction = profile1 - profile0;
    let radial_magnitude = |u: f64| {
        let point = profile0 + direction * u;
        let delta = point - origin;
        (delta - axis * delta.dot(axis)).magnitude()
    };
    let scale = direction
        .magnitude()
        .max((profile0 - origin).magnitude())
        .max(1.0);
    let mut safe = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0].into_iter().filter(|u| {
        radial_magnitude(*u)
            > DEFAULT_NUMERICAL_TOLERANCES.same_domain_direction_floor_factor
                * KERNEL_EPSILON
                * scale
    });
    let u0 = safe.next()?;
    let u1 = safe.find(|u| (*u - u0).abs() > 0.25)?;
    let angle = 0.37;
    Some([
        processor_parameter(surface.orientation(), Point2::new(u0, angle)),
        processor_parameter(surface.orientation(), Point2::new(u1, angle)),
        processor_parameter(surface.orientation(), Point2::new(u0, angle + 1.0)),
    ])
}

fn affine_map_from_correspondences(source: [Point2; 3], target: [Point2; 3]) -> Option<Matrix3> {
    let source = Matrix3::from_cols(
        Vector3::new(source[0].x, source[0].y, 1.0),
        Vector3::new(source[1].x, source[1].y, 1.0),
        Vector3::new(source[2].x, source[2].y, 1.0),
    );
    let target = Matrix3::from_cols(
        Vector3::new(target[0].x, target[0].y, 1.0),
        Vector3::new(target[1].x, target[1].y, 1.0),
        Vector3::new(target[2].x, target[2].y, 1.0),
    );
    Some(target * source.invert()?)
}

/// Returns the affine map from `source` UV coordinates into `target` UV
/// coordinates when both carriers are the same analytic support and their chart
/// relationship is affine.
///
/// The map is derived from analytic inverse parameterizations, never Newton or
/// sampled fitting. Planes and the STEP elementary cylinder/cone/torus charts
/// have affine relationships by construction. Spheres are deliberately
/// excluded: rotating a sphere's polar axis generally induces a nonlinear
/// spherical-coordinate map.
pub fn affine_same_domain_parameter_map(source: &Surface, target: &Surface) -> Option<Matrix3> {
    same_surface_domain(
        source,
        target,
        KERNEL_EPSILON,
        DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity,
    )?;

    let source_parameters = match (source, target) {
        (
            Surface::Elementary(ElementarySurface::Plane(_)),
            Surface::Elementary(ElementarySurface::Plane(_)),
        ) => [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ],
        (
            Surface::Elementary(ElementarySurface::CylindricalSurface(source)),
            Surface::Elementary(ElementarySurface::CylindricalSurface(_)),
        )
        | (
            Surface::Elementary(ElementarySurface::ConicalSurface(source)),
            Surface::Elementary(ElementarySurface::ConicalSurface(_)),
        ) => safe_revolution_basis(source)?,
        (
            Surface::Elementary(ElementarySurface::ToroidalSurface(source)),
            Surface::Elementary(ElementarySurface::ToroidalSurface(_)),
        ) => [
            processor_parameter(source.orientation(), Point2::new(0.37, 0.0)),
            processor_parameter(source.orientation(), Point2::new(1.37, 0.0)),
            processor_parameter(source.orientation(), Point2::new(0.37, 1.0)),
        ],
        _ => return None,
    };

    let target_base = exact_supported_parameter(
        target,
        source.evaluate(source_parameters[0].x, source_parameters[0].y),
    )?;
    let mut target_parameters = [target_base; 3];
    for index in 1..3 {
        let raw = exact_supported_parameter(
            target,
            source.evaluate(source_parameters[index].x, source_parameters[index].y),
        )?;
        target_parameters[index] = unwrap_parameter_near(target, raw, target_base);
    }
    let map = affine_map_from_correspondences(source_parameters, target_parameters)?;

    // Defense in depth only: the carrier families above already prove that the
    // chart relation is affine and the coefficients came from analytic inverse
    // maps. These extra points catch future parameterization changes and must
    // never be used to promote an unsupported family into the exact path.
    for uv in [
        Point2::new(0.23, 0.71),
        Point2::new(-0.41, 1.37),
        Point2::new(1.83, -0.29),
        Point2::new(2.31, 2.17),
    ] {
        let mapped = map.transform_point(uv);
        if !source
            .evaluate(uv.x, uv.y)
            .near(&target.evaluate(mapped.x, mapped.y))
        {
            return None;
        }
    }
    Some(map)
}

/// Prove that two STEP surface carriers occupy the same analytic support.
///
/// The comparison is support-based rather than parameter-based: placements,
/// angular origins, and parameter directions may differ.  Only exact analytic
/// support classes currently understood by Monstertruck are admitted. Planar
/// B-spline/NURBS carriers can prove geometric coplanarity, but are deliberately
/// refused here until their UV Jacobian is also proven to keep one orientation
/// globally. Unsupported or orientation-ambiguous classes return `None` rather
/// than being sampled into an optimistic answer.
///
/// The returned relation includes the orientation induced by each carrier's
/// parameterization, allowing topology to distinguish a harmless reversed
/// parameterization from an oppositely oriented physical face.
#[must_use]
pub fn same_surface_domain(
    left: &Surface,
    right: &Surface,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> Option<SurfaceDomainRelation> {
    if !linear_tolerance.is_finite()
        || linear_tolerance < 0.0
        || !angular_tolerance.is_finite()
        || angular_tolerance < 0.0
    {
        return None;
    }
    let relation =
        carrier(left)?.same_domain(carrier(right)?, linear_tolerance, angular_tolerance)?;
    Some(match relation {
        CarrierDomainRelation::SameOrientation => SurfaceDomainRelation::SameOrientation,
        CarrierDomainRelation::ReversedOrientation => SurfaceDomainRelation::ReversedOrientation,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn surface(elementary: ElementarySurface) -> Surface { Surface::Elementary(elementary) }

    #[test]
    fn sheared_elementary_surfaces_are_not_analytic_same_domain_supports() {
        let shear = Matrix4::from_cols(
            Vector4::new(1.0, 0.0, 0.0, 0.0),
            Vector4::new(0.5, 3.0_f64.sqrt() * 0.5, 0.0, 0.0),
            Vector4::new(0.0, 0.0, 1.0, 0.0),
            Vector4::new(0.0, 0.0, 0.0, 1.0),
        );
        let sphere = surface(ElementarySurface::Sphere(Processor::with_transform(
            crate::step::load::step_geometry::Sphere(monstertruck_geometry::prelude::Sphere::new(
                Point3::origin(),
                2.0,
            )),
            shear,
        )));
        assert!(carrier(&sphere).is_none());

        let profile = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0));
        let cylinder = surface(ElementarySurface::CylindricalSurface(
            Processor::with_transform(
                RevolutionSurface::by_revolution(profile, Point3::origin(), Vector3::unit_z()),
                shear,
            ),
        ));
        assert!(carrier(&cylinder).is_none());
    }

    #[test]
    fn differently_parameterized_coplanar_planes_are_same_domain() {
        let a = surface(ElementarySurface::Plane(Plane::new(
            Point3::new(0.0, 0.0, 2.0),
            Point3::new(1.0, 0.0, 2.0),
            Point3::new(0.0, 1.0, 2.0),
        )));
        let b = surface(ElementarySurface::Plane(Plane::new(
            Point3::new(4.0, -7.0, 2.0),
            Point3::new(4.0, -6.0, 2.0),
            Point3::new(3.0, -7.0, 2.0),
        )));
        assert_eq!(
            same_surface_domain(&a, &b, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::SameOrientation)
        );
    }

    #[test]
    fn opposite_plane_parameterization_reports_reversed_orientation() {
        let a = surface(ElementarySurface::Plane(Plane::xy()));
        let b = surface(ElementarySurface::Plane(Plane::new(
            Point3::origin(),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        )));
        assert_eq!(
            same_surface_domain(&a, &b, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::ReversedOrientation)
        );
    }

    fn cylinder(origin: Point3, axis: Vector3, radius: f64) -> Surface {
        let radial = if axis.cross(Vector3::unit_x()).magnitude() > 0.5 {
            axis.cross(Vector3::unit_x()).normalize() * radius
        } else {
            axis.cross(Vector3::unit_y()).normalize() * radius
        };
        let profile = Line(origin + radial, origin + radial + axis.normalize());
        surface(ElementarySurface::CylindricalSurface(Processor::new(
            RevolutionSurface::by_revolution(profile, origin, axis),
        )))
    }

    fn cone(origin: Point3, axis: Vector3, radius: f64, slope: f64) -> Surface {
        let axis = axis.normalize();
        let radial = if axis.cross(Vector3::unit_x()).magnitude() > 0.5 {
            axis.cross(Vector3::unit_x()).normalize()
        } else {
            axis.cross(Vector3::unit_y()).normalize()
        };
        let profile = Line(
            origin + radial * radius,
            origin + radial * (radius + slope) + axis,
        );
        surface(ElementarySurface::ConicalSurface(Processor::new(
            RevolutionSurface::by_revolution(profile, origin, axis),
        )))
    }

    #[test]
    fn cylinders_ignore_axis_origin_and_angular_zero() {
        let a = cylinder(Point3::origin(), Vector3::unit_z(), 2.0);
        let b = cylinder(Point3::new(0.0, 0.0, 17.0), -Vector3::unit_z(), 2.0);
        assert!(same_surface_domain(&a, &b, 1.0e-9, 1.0e-9).is_some());
        let wrong = cylinder(Point3::new(0.01, 0.0, 0.0), Vector3::unit_z(), 2.0);
        assert_eq!(same_surface_domain(&a, &wrong, 1.0e-6, 1.0e-9), None);
    }

    #[test]
    fn cylinder_profile_direction_controls_global_orientation_exactly() {
        let axis = Vector3::unit_z();
        let forward_profile = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0));
        let reverse_profile = forward_profile.inverse();
        let source = surface(ElementarySurface::CylindricalSurface(Processor::new(
            RevolutionSurface::by_revolution(forward_profile, Point3::origin(), axis),
        )));
        let target = surface(ElementarySurface::CylindricalSurface(Processor::new(
            RevolutionSurface::by_revolution(reverse_profile, Point3::origin(), axis),
        )));
        assert_eq!(
            same_surface_domain(&source, &target, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::ReversedOrientation)
        );

        let mut compensated = Processor::new(RevolutionSurface::by_revolution(
            reverse_profile,
            Point3::origin(),
            axis,
        ));
        compensated.invert();
        let compensated = surface(ElementarySurface::CylindricalSurface(compensated));
        assert_eq!(
            same_surface_domain(&source, &compensated, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::SameOrientation)
        );
    }

    #[test]
    fn sphere_support_ignores_parameter_orientation_but_preserves_normal_relation() {
        let raw = crate::step::load::step_geometry::Sphere(
            monstertruck_geometry::prelude::Sphere::new(Point3::new(1.0, 2.0, 3.0), 4.0),
        );
        let a = surface(ElementarySurface::Sphere(Processor::new(raw)));
        let mut reversed = Processor::new(raw);
        reversed.invert();
        let b = surface(ElementarySurface::Sphere(reversed));
        assert_eq!(
            same_surface_domain(&a, &b, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::ReversedOrientation)
        );
    }

    #[test]
    fn rotated_sphere_chart_preserves_support_orientation() {
        let raw = crate::step::load::step_geometry::Sphere(
            monstertruck_geometry::prelude::Sphere::new(Point3::origin(), 4.0),
        );
        let source = surface(ElementarySurface::Sphere(Processor::new(raw)));
        let target = surface(ElementarySurface::Sphere(Processor::with_transform(
            raw,
            Matrix4::from_axis_angle(Vector3::unit_z(), Rad(std::f64::consts::PI)),
        )));
        assert_eq!(
            same_surface_domain(&source, &target, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::SameOrientation)
        );
    }

    #[test]
    fn cone_support_is_apex_axis_and_opening_angle() {
        let a = cone(Point3::new(0.0, 0.0, 5.0), Vector3::unit_z(), 2.0, 0.5);
        let b = cone(Point3::new(0.0, 0.0, 9.0), Vector3::unit_z(), 4.0, 0.5);
        assert!(same_surface_domain(&a, &b, 1.0e-8, 1.0e-9).is_some());
        let wrong = cone(Point3::new(0.0, 0.0, 5.0), Vector3::unit_z(), 2.0, 0.6);
        assert_eq!(same_surface_domain(&a, &wrong, 1.0e-8, 1.0e-6), None);
    }

    #[test]
    fn torus_support_compares_axis_and_both_radii() {
        let a = surface(ElementarySurface::ToroidalSurface(Processor::new(
            Torus::new(Point3::new(1.0, 2.0, 3.0), 5.0, 1.0).expect("valid test torus"),
        )));
        let b = a.clone();
        assert_eq!(
            same_surface_domain(&a, &b, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::SameOrientation)
        );
        let wrong = surface(ElementarySurface::ToroidalSurface(Processor::new(
            Torus::new(Point3::new(1.0, 2.0, 3.0), 5.0, 1.01).expect("valid comparison torus"),
        )));
        assert_eq!(same_surface_domain(&a, &wrong, 1.0e-6, 1.0e-9), None);
    }

    #[test]
    fn affine_parameter_map_transports_rotated_plane_chart() {
        let source = surface(ElementarySurface::Plane(Plane::xy()));
        let target = surface(ElementarySurface::Plane(Plane::new(
            Point3::origin(),
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(-3.0, 0.0, 0.0),
        )));
        let map = affine_same_domain_parameter_map(&source, &target).unwrap();
        for uv in [Point2::new(0.2, 0.7), Point2::new(-1.3, 2.1)] {
            let mapped = map.transform_point(uv);
            assert!(
                source
                    .evaluate(uv.x, uv.y)
                    .near(&target.evaluate(mapped.x, mapped.y))
            );
        }
    }

    #[test]
    fn affine_parameter_map_transports_shifted_cylinder_chart() {
        let axis = Vector3::unit_z();
        let source_profile = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0));
        let target_profile = Line(Point3::new(0.0, 2.0, 0.0), Point3::new(0.0, 2.0, 1.0));
        let mut source_cylinder = Processor::new(RevolutionSurface::by_revolution(
            source_profile,
            Point3::origin(),
            axis,
        ));
        let mut target_cylinder = Processor::new(RevolutionSurface::by_revolution(
            target_profile,
            Point3::origin(),
            axis,
        ));
        source_cylinder.invert();
        target_cylinder.invert();
        let source = surface(ElementarySurface::CylindricalSurface(source_cylinder));
        let target = surface(ElementarySurface::CylindricalSurface(target_cylinder));
        let map = affine_same_domain_parameter_map(&source, &target).unwrap();
        for uv in [Point2::new(0.2, -3.0), Point2::new(5.9, 7.0)] {
            let mapped = map.transform_point(uv);
            assert!(
                source
                    .evaluate(uv.x, uv.y)
                    .near(&target.evaluate(mapped.x, mapped.y))
            );
        }
    }

    #[test]
    fn affine_parameter_map_handles_processor_uv_axis_swap() {
        let axis = Vector3::unit_z();
        let source_profile = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0));
        let target_profile = Line(Point3::new(0.0, 2.0, 0.0), Point3::new(0.0, 2.0, 1.0));
        let source = surface(ElementarySurface::CylindricalSurface(Processor::new(
            RevolutionSurface::by_revolution(source_profile, Point3::origin(), axis),
        )));
        let mut target_carrier = Processor::new(RevolutionSurface::by_revolution(
            target_profile,
            Point3::origin(),
            axis,
        ));
        target_carrier.invert();
        let target = surface(ElementarySurface::CylindricalSurface(target_carrier));
        let map = affine_same_domain_parameter_map(&source, &target).unwrap();
        for uv in [
            Point2::new(-3.0, 0.2),
            Point2::new(7.0, 5.9),
            Point2::new(0.5, -0.7),
        ] {
            let mapped = map.transform_point(uv);
            assert!(
                source
                    .evaluate(uv.x, uv.y)
                    .near(&target.evaluate(mapped.x, mapped.y)),
                "processor UV swap transport failed for {uv:?} -> {mapped:?}",
            );
        }
    }

    #[test]
    fn affine_parameter_map_transports_shifted_cone_chart() {
        let source = cone(Point3::new(0.0, 0.0, 5.0), Vector3::unit_z(), 2.0, 0.5);
        let target = cone(Point3::new(0.0, 0.0, 9.0), Vector3::unit_z(), 4.0, 0.5);
        let map = affine_same_domain_parameter_map(&source, &target).unwrap();
        for uv in [Point2::new(-1.3, 0.2), Point2::new(0.7, 2.4)] {
            let mapped = map.transform_point(uv);
            assert!(
                source
                    .evaluate(uv.x, uv.y)
                    .near(&target.evaluate(mapped.x, mapped.y))
            );
        }
    }

    #[test]
    fn affine_parameter_map_transports_rotated_torus_chart() {
        let torus = Torus::new(Point3::origin(), 5.0, 1.0).expect("valid ring torus");
        let source = surface(ElementarySurface::ToroidalSurface(Processor::new(torus)));
        let target = surface(ElementarySurface::ToroidalSurface(
            Processor::with_transform(
                torus,
                Matrix4::from_axis_angle(Vector3::unit_z(), Rad(0.73)),
            ),
        ));
        let map = affine_same_domain_parameter_map(&source, &target).unwrap();
        for uv in [Point2::new(0.2, 0.4), Point2::new(5.8, 4.9)] {
            let mapped = map.transform_point(uv);
            assert!(
                source
                    .evaluate(uv.x, uv.y)
                    .near(&target.evaluate(mapped.x, mapped.y))
            );
        }
    }

    #[test]
    fn spindle_torus_same_domain_orientation_is_refused() {
        let torus = Torus::new(Point3::origin(), 1.0, 2.0).expect("valid spindle torus");
        let source = surface(ElementarySurface::ToroidalSurface(Processor::new(torus)));
        let target = surface(ElementarySurface::ToroidalSurface(
            Processor::with_transform(
                torus,
                Matrix4::from_axis_angle(Vector3::unit_z(), Rad(0.61)),
            ),
        ));
        assert!(same_surface_domain(&source, &target, 1.0e-9, 1.0e-9).is_none());
        assert!(affine_same_domain_parameter_map(&source, &target).is_none());
    }

    #[test]
    fn nearby_parallel_plane_is_not_same_domain() {
        let a = surface(ElementarySurface::Plane(Plane::xy()));
        let b = surface(ElementarySurface::Plane(Plane::new(
            Point3::new(0.0, 0.0, 1.0e-4),
            Point3::new(1.0, 0.0, 1.0e-4),
            Point3::new(0.0, 1.0, 1.0e-4),
        )));
        assert_eq!(same_surface_domain(&a, &b, 1.0e-6, 1.0e-9), None);
    }

    #[test]
    fn reflected_sphere_reports_reversed_geometric_orientation() {
        let raw = crate::step::load::step_geometry::Sphere(
            monstertruck_geometry::prelude::Sphere::new(Point3::origin(), 4.0),
        );
        let source = surface(ElementarySurface::Sphere(Processor::new(raw)));
        let reflection = Matrix4::from_cols(
            Vector4::new(-1.0, 0.0, 0.0, 0.0),
            Vector4::new(0.0, 1.0, 0.0, 0.0),
            Vector4::new(0.0, 0.0, 1.0, 0.0),
            Vector4::new(0.0, 0.0, 0.0, 1.0),
        );
        let reflected = surface(ElementarySurface::Sphere(Processor::with_transform(
            raw, reflection,
        )));

        assert_eq!(
            same_surface_domain(&source, &reflected, 1.0e-9, 1.0e-9),
            Some(SurfaceDomainRelation::ReversedOrientation)
        );
    }

    #[test]
    fn invalid_same_domain_tolerances_are_refused() {
        let plane = surface(ElementarySurface::Plane(Plane::xy()));
        assert_eq!(same_surface_domain(&plane, &plane, -1.0, 1.0e-9), None);
        assert_eq!(same_surface_domain(&plane, &plane, f64::NAN, 1.0e-9), None);
        assert_eq!(same_surface_domain(&plane, &plane, 1.0e-9, f64::NAN), None);
    }
}
