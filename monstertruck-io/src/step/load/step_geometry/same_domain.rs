use super::{ElementarySurface, Surface};
use monstertruck_geometry::prelude::*;
use monstertruck_topology::repair::SurfaceDomainRelation;

#[derive(Clone, Copy, Debug)]
enum Support {
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
    },
    Torus {
        center: Point3,
        axis: Vector3,
        major_radius: f64,
        minor_radius: f64,
    },
}

#[inline]
fn canonical_direction(mut vector: Vector3) -> Option<Vector3> {
    let magnitude = vector.magnitude();
    if !magnitude.is_finite() || magnitude <= TOLERANCE {
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

#[inline]
fn similarity_scale(transform: Matrix4) -> Option<f64> {
    const SIMILARITY_TOLERANCE: f64 = 1.0e-10;
    if transform.x.w.abs() > SIMILARITY_TOLERANCE
        || transform.y.w.abs() > SIMILARITY_TOLERANCE
        || transform.z.w.abs() > SIMILARITY_TOLERANCE
        || (transform.w.w - 1.0).abs() > SIMILARITY_TOLERANCE
    {
        return None;
    }

    let columns = [
        transform.x.truncate(),
        transform.y.truncate(),
        transform.z.truncate(),
    ];
    let scale = columns[0].magnitude();
    if !scale.is_finite() || scale <= TOLERANCE {
        return None;
    }
    let length_tolerance = SIMILARITY_TOLERANCE * scale.max(1.0);
    let dot_tolerance = SIMILARITY_TOLERANCE * scale * scale;
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

fn revolution_support(
    surface: &Processor<RevolutionSurface<Line<Point3>>, Matrix4>,
    cylinder: bool,
) -> Option<Support> {
    let transform = *surface.transform();
    similarity_scale(transform)?;
    let revolution = surface.entity();
    let origin = transform.transform_point(revolution.origin());
    let axis = canonical_direction(transform.transform_vector(revolution.axis()))?;
    let Line(p0, p1) = *revolution.entity_curve();
    let p0 = transform.transform_point(p0);
    let p1 = transform.transform_point(p1);
    let direction = p1 - p0;
    let radial0 = (p0 - origin) - axis * (p0 - origin).dot(axis);
    let radial_direction = direction - axis * direction.dot(axis);

    if cylinder {
        let scale = direction.magnitude().max(radial0.magnitude()).max(1.0);
        if radial_direction.magnitude() > 1.0e-10 * scale {
            return None;
        }
        let radius = radial0.magnitude();
        (radius > TOLERANCE).then_some(Support::Cylinder {
            origin,
            axis,
            radius,
        })
    } else {
        let radial_direction2 = radial_direction.magnitude2();
        if radial_direction2 <= TOLERANCE.powi(2) {
            return None;
        }
        // Intersect the generating line with the revolution axis in the radial
        // quotient. For a true cone the residual is zero; reject skew profiles
        // instead of quietly manufacturing an apex.
        let parameter = -radial0.dot(radial_direction) / radial_direction2;
        let apex = p0 + direction * parameter;
        let apex_radial = (apex - origin) - axis * (apex - origin).dot(axis);
        let scale = (p1 - p0)
            .magnitude()
            .max((apex - origin).magnitude())
            .max(1.0);
        if apex_radial.magnitude() > 1.0e-9 * scale {
            return None;
        }
        let axial = direction.dot(axis).abs();
        if axial <= TOLERANCE {
            return None;
        }
        Some(Support::Cone {
            apex,
            axis,
            slope: radial_direction.magnitude() / axial,
        })
    }
}

fn elementary_support(surface: &ElementarySurface) -> Option<Support> {
    match surface {
        ElementarySurface::Plane(plane) => Some(Support::Plane {
            origin: plane.origin(),
            normal: canonical_direction(plane.normal())?,
        }),
        ElementarySurface::Sphere(sphere) => {
            let transform = *sphere.transform();
            let scale = similarity_scale(transform)?;
            Some(Support::Sphere {
                center: transform.transform_point(sphere.entity().0.center()),
                radius: sphere.entity().0.radius() * scale,
            })
        }
        ElementarySurface::CylindricalSurface(cylinder) => revolution_support(cylinder, true),
        ElementarySurface::ConicalSurface(cone) => revolution_support(cone, false),
        ElementarySurface::ToroidalSurface(torus) => {
            let transform = *torus.transform();
            let scale = similarity_scale(transform)?;
            Some(Support::Torus {
                center: transform.transform_point(torus.entity().center()),
                axis: canonical_direction(transform.transform_vector(Vector3::unit_z()))?,
                major_radius: torus.entity().large_radius() * scale,
                minor_radius: torus.entity().small_radius() * scale,
            })
        }
    }
}

fn support(surface: &Surface) -> Option<Support> {
    match surface {
        Surface::ElementarySurface(surface) => elementary_support(surface),
        Surface::BsplineSurface(surface) => match surface.try_into_analytic_surface_kind()? {
            AnalyticSurfaceKind::Plane(plane) => Some(Support::Plane {
                origin: plane.origin(),
                normal: canonical_direction(plane.normal())?,
            }),
            _ => None,
        },
        Surface::NurbsSurface(surface) => match surface.try_into_analytic_surface_kind()? {
            AnalyticSurfaceKind::Plane(plane) => Some(Support::Plane {
                origin: plane.origin(),
                normal: canonical_direction(plane.normal())?,
            }),
            _ => None,
        },
        Surface::SweepSurface(_) => None,
    }
}

#[inline]
fn parallel(left: Vector3, right: Vector3, angular_tolerance: f64) -> bool {
    left.cross(right).magnitude() <= angular_tolerance.sin().abs().max(64.0 * f64::EPSILON)
}

#[inline]
fn axis_line_distance(origin0: Point3, axis0: Vector3, origin1: Point3) -> f64 {
    (origin1 - origin0).cross(axis0).magnitude()
}

fn same_support(
    left: Support,
    right: Support,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> bool {
    match (left, right) {
        (
            Support::Plane {
                origin: o0,
                normal: n0,
            },
            Support::Plane {
                origin: o1,
                normal: n1,
            },
        ) => parallel(n0, n1, angular_tolerance) && (o1 - o0).dot(n0).abs() <= linear_tolerance,
        (
            Support::Sphere {
                center: c0,
                radius: r0,
            },
            Support::Sphere {
                center: c1,
                radius: r1,
            },
        ) => c0.distance(c1) <= linear_tolerance && (r0 - r1).abs() <= linear_tolerance,
        (
            Support::Cylinder {
                origin: o0,
                axis: a0,
                radius: r0,
            },
            Support::Cylinder {
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
            Support::Cone {
                apex: p0,
                axis: a0,
                slope: s0,
            },
            Support::Cone {
                apex: p1,
                axis: a1,
                slope: s1,
            },
        ) => {
            p0.distance(p1) <= linear_tolerance
                && parallel(a0, a1, angular_tolerance)
                && (s0.atan() - s1.atan()).abs() <= angular_tolerance
        }
        (
            Support::Torus {
                center: c0,
                axis: a0,
                major_radius: r0,
                minor_radius: m0,
            },
            Support::Torus {
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

#[inline]
fn processor_orientation_sign<S>(surface: &Processor<S, Matrix4>) -> Option<f64> {
    similarity_scale(*surface.transform())?;
    let determinant = surface.transform().determinant();
    if !determinant.is_finite() || determinant.abs() <= TOLERANCE {
        return None;
    }
    let chart = if surface.orientation() { 1.0 } else { -1.0 };
    Some(chart * determinant.signum())
}

/// Global orientation of a carrier relative to the orientation-free analytic
/// support used by [`same_support`].
///
/// This is deliberately structural. A few sampled normals are not evidence of
/// a global orientation: spheres change Cartesian octants, and a spindle torus
/// actually changes the relation between its parametric normal and the implicit
/// gradient across the self-intersection. Unsupported/global-ambiguous cases
/// therefore return `None`.
fn analytic_orientation_sign(surface: &Surface, support: Support) -> Option<bool> {
    let sign = match (surface, support) {
        (
            Surface::ElementarySurface(ElementarySurface::Plane(plane)),
            Support::Plane { normal, .. },
        ) => {
            let dot = plane.normal().dot(normal);
            if dot.abs() <= TOLERANCE {
                return None;
            }
            dot.signum()
        }
        (Surface::ElementarySurface(ElementarySurface::Sphere(sphere)), Support::Sphere { .. }) => {
            processor_orientation_sign(sphere)?
        }
        (
            Surface::ElementarySurface(ElementarySurface::CylindricalSurface(revolution)),
            Support::Cylinder { .. },
        )
        | (
            Surface::ElementarySurface(ElementarySurface::ConicalSurface(revolution)),
            Support::Cone { .. },
        ) => {
            let Line(front, back) = *revolution.entity().entity_curve();
            let axial_step = (back - front).dot(revolution.entity().axis());
            if !axial_step.is_finite() || axial_step.abs() <= TOLERANCE {
                return None;
            }
            -axial_step.signum() * processor_orientation_sign(revolution)?
        }
        (
            Surface::ElementarySurface(ElementarySurface::ToroidalSurface(torus)),
            Support::Torus {
                major_radius,
                minor_radius,
                ..
            },
        ) => {
            // For R < r, the signed radial factor R + r cos(v) changes sign.
            // The regular parametric normal does not, while the implicit
            // gradient does, so there is no one global orientation relation on
            // the self-intersecting spindle support. Refuse it rather than
            // assigning a sign from whichever branch happened to be sampled.
            if major_radius + TOLERANCE < minor_radius {
                return None;
            }
            processor_orientation_sign(torus)?
        }
        // Coplanarity proves a spline surface's image lies on a plane, but not
        // that its UV Jacobian keeps one sign globally. Until that stronger
        // proof exists, generic spline planes are support evidence only, not
        // orientation evidence for a topology mutation.
        (Surface::BsplineSurface(_), Support::Plane { .. })
        | (Surface::NurbsSurface(_), Support::Plane { .. }) => return None,
        _ => return None,
    };
    Some(sign >= 0.0)
}

#[inline]
fn unwrap_parameter_component(value: f64, reference: f64, period: Option<f64>) -> f64 {
    period
        .filter(|period| period.is_finite() && *period > TOLERANCE)
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
    (parameter.z.abs() <= TOLERANCE * scale).then_some(Point2::new(parameter.x, parameter.y))
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
    if axial_step.abs() <= TOLERANCE * scale {
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
    if profile_radius <= TOLERANCE * radial_scale
        || point_radius <= TOLERANCE * radial_scale
        || (profile_radius - point_radius).abs() > TOLERANCE * radial_scale
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
        if (signed_radial.abs() - rho).abs() > TOLERANCE * scale {
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
        Surface::ElementarySurface(ElementarySurface::Plane(plane)) => {
            exact_plane_parameter(plane, point)
        }
        Surface::ElementarySurface(
            ElementarySurface::CylindricalSurface(revolution)
            | ElementarySurface::ConicalSurface(revolution),
        ) => exact_linear_revolution_parameter(revolution, point),
        Surface::ElementarySurface(ElementarySurface::ToroidalSurface(torus)) => {
            exact_torus_parameter(torus, point)
        }
        Surface::ElementarySurface(ElementarySurface::Sphere(_))
        | Surface::BsplineSurface(_)
        | Surface::NurbsSurface(_)
        | Surface::SweepSurface(_) => None,
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
    let mut safe = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
        .into_iter()
        .filter(|u| radial_magnitude(*u) > 16.0 * TOLERANCE * scale);
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
    same_surface_domain(source, target, TOLERANCE, 1.0e-10)?;

    let source_parameters = match (source, target) {
        (
            Surface::ElementarySurface(ElementarySurface::Plane(_)),
            Surface::ElementarySurface(ElementarySurface::Plane(_)),
        ) => [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ],
        (
            Surface::ElementarySurface(ElementarySurface::CylindricalSurface(source)),
            Surface::ElementarySurface(ElementarySurface::CylindricalSurface(_)),
        )
        | (
            Surface::ElementarySurface(ElementarySurface::ConicalSurface(source)),
            Surface::ElementarySurface(ElementarySurface::ConicalSurface(_)),
        ) => safe_revolution_basis(source)?,
        (
            Surface::ElementarySurface(ElementarySurface::ToroidalSurface(source)),
            Surface::ElementarySurface(ElementarySurface::ToroidalSurface(_)),
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
    let left_support = support(left)?;
    let right_support = support(right)?;
    if !same_support(
        left_support,
        right_support,
        linear_tolerance,
        angular_tolerance,
    ) {
        return None;
    }
    let left_orientation = analytic_orientation_sign(left, left_support)?;
    let right_orientation = analytic_orientation_sign(right, right_support)?;
    Some(if left_orientation == right_orientation {
        SurfaceDomainRelation::SameOrientation
    } else {
        SurfaceDomainRelation::ReversedOrientation
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn surface(elementary: ElementarySurface) -> Surface { Surface::ElementarySurface(elementary) }

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
        assert!(support(&sphere).is_none());

        let profile = Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0));
        let cylinder = surface(ElementarySurface::CylindricalSurface(
            Processor::with_transform(
                RevolutionSurface::by_revolution(profile, Point3::origin(), Vector3::unit_z()),
                shear,
            ),
        ));
        assert!(support(&cylinder).is_none());
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
}
