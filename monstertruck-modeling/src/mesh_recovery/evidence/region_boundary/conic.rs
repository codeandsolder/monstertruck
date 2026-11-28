//! Exact rational planar sections of affine-radius revolutions.
//! Solving the plane equation in azimuth gives a rational quadratic
//! parametrization of ellipse, parabola and hyperbola branches alike.

use super::SurfaceCarrier;
use crate::{BsplineCurve, Curve, KnotVector, NurbsCurve, Point3, Vector3, Vector4};
use std::f64::consts::{FRAC_PI_2, TAU};

/// Independently constrained noncircular conic on a linear surface of
/// revolution and either an observed or a derived section plane.
///
/// Its observed finite branch is certified by
/// source sample agreement and monotone azimuth, with no asymptote or
/// negative radial sheet between adjacent samples.
#[derive(Clone, Copy, Debug)]
pub struct ContactConicEvidence {
    /// Plane origin.
    pub plane_origin: Point3,
    /// Unit plane normal.
    pub plane_normal: Vector3,
    /// Rotational axis origin.
    pub origin: Point3,
    /// Unit rotational axis.
    pub axis: Vector3,
    /// Profile radius at the axis origin.
    pub radius_at_origin: f64,
    /// Affine profile slope.
    pub slope: f64,
    /// Signed unwrapped azimuth sweep of the observed path.
    pub signed_sweep_radians: f64,
    /// Maximum absolute adjacent source azimuth step.
    pub max_step_radians: f64,
    /// Maximum distance from a source sample to its conic point at the same azimuth.
    pub max_vertex_residual: f64,
    // A derived section plane may also represent a second revolution.
    secondary: Option<SurfaceCarrier>,
}

impl ContactConicEvidence {
    const fn from_carriers(carriers: &[SurfaceCarrier]) -> Option<Self> {
        let (plane, linear) = match carriers {
            [
                plane @ SurfaceCarrier::Plane { .. },
                linear @ SurfaceCarrier::LinearRevolution { .. },
            ]
            | [
                linear @ SurfaceCarrier::LinearRevolution { .. },
                plane @ SurfaceCarrier::Plane { .. },
            ] => (*plane, *linear),
            _ => return None,
        };
        let SurfaceCarrier::Plane {
            origin: plane_origin,
            normal: plane_normal,
        } = plane
        else {
            return None;
        };
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } = linear
        else {
            return None;
        };
        Some(Self {
            plane_origin,
            plane_normal,
            origin,
            axis,
            radius_at_origin,
            slope,
            signed_sweep_radians: 0.0,
            max_step_radians: 0.0,
            max_vertex_residual: 0.0,
            secondary: None,
        })
    }

    // For two equal-radius intersecting cylinders, subtraction of their
    // squared implicit equations factors into the two axial bisector planes.
    // Propose those planes even for numerically fitted carriers, but validate
    // the resulting complete rational sections against the actual second
    // carrier rather than assuming exact symmetry.
    fn cylinder_sections(carriers: &[SurfaceCarrier], tolerance: f64) -> Vec<Self> {
        let [
            SurfaceCarrier::LinearRevolution {
                origin: first_origin,
                axis: first_axis,
                radius_at_origin: first_radius,
                slope: first_slope,
            },
            SurfaceCarrier::LinearRevolution {
                origin: second_origin,
                axis: second_axis,
                radius_at_origin: second_radius,
                slope: second_slope,
            },
        ] = carriers
        else {
            return Vec::new();
        };
        let cross = first_axis.cross(second_axis);
        let denominator = cross.norm_squared();
        if denominator <= 64.0 * f64::EPSILON {
            return Vec::new();
        }
        let delta = second_origin - first_origin;
        let first_axial = delta.cross(second_axis).dot(&cross) / denominator;
        let second_axial = delta.cross(first_axis).dot(&cross) / denominator;
        let left = first_origin + first_axis * first_axial;
        let right = second_origin + second_axis * second_axial;
        let radius = first_slope.mul_add(first_axial, *first_radius);
        let other_radius = second_slope.mul_add(second_axial, *second_radius);
        if !radius.is_finite()
            || !other_radius.is_finite()
            || radius.min(other_radius) <= 0.0
            || (left - right).norm() > tolerance
            || (radius - other_radius).abs() > tolerance
            || first_slope.abs().max(second_slope.abs()) * radius.max(other_radius) > tolerance
        {
            return Vec::new();
        }
        let origin = left + (right - left) * 0.5;
        [-1.0, 1.0]
            .into_iter()
            .filter_map(|sign| {
                let plane = SurfaceCarrier::Plane {
                    origin,
                    normal: (*first_axis + *second_axis * sign).normalize(),
                };
                let mut section = Self::from_carriers(&[plane, carriers[0]])?;
                section.secondary = Some(carriers[1]);
                Some(section)
            })
            .collect()
    }

    pub(super) fn auxiliary_plane(self) -> Option<SurfaceCarrier> {
        self.secondary.map(|_| SurfaceCarrier::Plane {
            origin: self.plane_origin,
            normal: self.plane_normal,
        })
    }

    fn direction(self, point: Point3) -> Option<Vector3> {
        let delta = point - self.origin;
        let radial = delta - self.axis * delta.dot(&self.axis);
        let norm = radial.norm();
        (norm.is_finite() && norm > 0.0).then(|| radial / norm)
    }

    pub(super) fn angular_step(self, a: Point3, b: Point3) -> Option<f64> {
        let (u, v) = (self.direction(a)?, self.direction(b)?);
        Some(self.axis.dot(&u.cross(&v)).atan2(u.dot(&v)))
    }

    // Let D/W be the rational circle direction (homogeneous quadratic).
    // The plane equation gives z = Z/H with
    // H = a W + slope n·D; Z = -c W - r0 n·D.
    // Then r = (r0 a - slope c) W/H. The intersection point is
    // [origin H + axis Z + (r0 a - slope c) D] / H.
    // Thus a homogeneous linear map of the circular controls gives the
    // exact conic; no sampled spline approximation is needed.
    fn homogeneous(self, direction: Vector3, weight: f64) -> Option<Vector4> {
        let axial_normal = self.plane_normal.dot(&self.axis);
        let plane_offset = self.plane_normal.dot(&(self.origin - self.plane_origin));
        let radial_normal = self.plane_normal.dot(&direction);
        let denominator = self.slope.mul_add(radial_normal, axial_normal * weight);
        let axial_numerator = -self
            .radius_at_origin
            .mul_add(radial_normal, plane_offset * weight);
        let radial_numerator = self
            .radius_at_origin
            .mul_add(axial_normal, -self.slope * plane_offset);
        // Strict positive radial sheet. A zero numerator corresponds to
        // an apex or generator degeneracy, not axial_normal regular conic branch.
        let h_error = 16.0
            * f64::EPSILON
            * ((axial_normal * weight).abs() + (self.slope * radial_normal).abs());
        let k_error = 16.0
            * f64::EPSILON
            * ((self.radius_at_origin * axial_normal).abs() + (self.slope * plane_offset).abs());
        if !denominator.is_finite()
            || !radial_numerator.is_finite()
            || denominator.abs() <= h_error
            || radial_numerator.abs() <= k_error
            || radial_numerator / denominator <= 0.0
        {
            return None;
        }
        let point = self.origin.to_homogeneous() * denominator
            + self.axis.push(0.0) * axial_numerator
            + direction.push(0.0) * radial_numerator;
        (point.x.is_finite() && point.y.is_finite() && point.z.is_finite()).then_some(point)
    }

    fn projected(self, point: Point3) -> Option<Point3> {
        let p = self.homogeneous(self.direction(point)?, 1.0)?;
        let q = Point3::new(p.x / p.w, p.y / p.w, p.z / p.w);
        (q.x.is_finite() && q.y.is_finite() && q.z.is_finite()).then_some(q)
    }

    fn controls(self, a: Point3, b: Point3) -> Option<[Vector4; 3]> {
        let u = self.direction(a)?;
        let v = self.direction(b)?;
        let step = self.angular_step(a, b)?;
        if step.abs() > f64::EPSILON.mul_add(16.0, FRAC_PI_2) || step == 0.0 {
            return None;
        }
        let midpoint = (u + v).normalize();
        let controls = [
            self.homogeneous(u, 1.0)?,
            self.homogeneous(midpoint, (step * 0.5).cos())?,
            self.homogeneous(v, 1.0)?,
        ];
        // Bernstein convex hull guarantees a denominator of constant
        // nonzero sign on the full interval. Also guarantees r>0 by
        // the sign test above; this excludes branch/asymptote crossings.
        let sign = controls[0].w.signum();
        if controls
            .iter()
            .any(|p| p.w.is_sign_positive() != sign.is_sign_positive())
        {
            return None;
        }
        let scale = controls.iter().map(|p| p.w.abs()).fold(0.0_f64, f64::max);
        Some(controls.map(|p| p * (sign / scale)))
    }

    // Restrict the secondary surface's squared implicit equation to the
    // rational quadratic. In homogeneous coordinates it is a quartic
    // Bernstein polynomial. Its coefficients bound the full interval,
    // including unsampled interior points. Positive weights and a positive
    // affine radius then convert that bound to an ordinary radial residual.
    fn certifies_secondary(self, controls: &[Vector4; 3], tolerance: f64) -> bool {
        let Some(SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        }) = self.secondary
        else {
            return self.secondary.is_none();
        };
        let origin = origin.coords;
        let mut local_scale = 0.0_f64;
        let mut coordinate_scale = 0.0_f64;
        let mut min_weight = f64::INFINITY;
        let mut min_radius = f64::INFINITY;
        let data = controls.map(|p| {
            let delta = p.fixed_rows::<3>(0).into_owned() - origin * p.w;
            let axial = delta.dot(&axis);
            let radial = delta - axis * axial;
            let radius = slope.mul_add(axial, radius_at_origin * p.w);
            coordinate_scale = coordinate_scale.max(
                origin
                    .norm()
                    .mul_add(p.w, p.fixed_rows::<3>(0).into_owned().norm()),
            );
            local_scale = local_scale.max(delta.norm() + radius.abs());
            min_weight = min_weight.min(p.w);
            min_radius = min_radius.min(radius / p.w);
            (radial, radius)
        });
        if min_weight <= 0.0 || min_radius <= 0.0 {
            return false;
        }
        let q = |i: usize, j: usize| (-data[i].1).mul_add(data[j].1, data[i].0.dot(&data[j].0));
        let bound = [
            q(0, 0),
            q(0, 1),
            2.0_f64.mul_add(q(1, 1), q(0, 2)) / 3.0,
            q(1, 2),
            q(2, 2),
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
        // Account for cancellation when changing homogeneous coordinate
        // origins and evaluating the quadratic form.
        let roundoff = 128.0 * f64::EPSILON * coordinate_scale.max(local_scale) * local_scale;
        let budget = tolerance * min_radius * min_weight * min_weight;
        bound.is_finite()
            && roundoff.is_finite()
            && budget.is_finite()
            && bound + roundoff <= budget
    }

    pub(super) fn native_curve(
        self,
        a: Point3,
        b: Point3,
        sweep: f64,
        tolerance: f64,
    ) -> Option<Curve> {
        let step = self.angular_step(a, b)?;
        // The reconstructed endpoints must select the same observed
        // branch. Endpoint consistency is also checked by native Edge.
        let rho = |p: Point3| {
            let delta = p - self.origin;
            (delta - self.axis * delta.dot(&self.axis)).norm()
        };
        if step * sweep <= 0.0
            || (step - sweep).abs() * rho(a).max(rho(b)) > tolerance * 4.0
            || (self.projected(a)? - a).norm() > tolerance
            || (self.projected(b)? - b).norm() > tolerance
        {
            return None;
        }
        let controls = self.controls(a, b)?;
        if !self.certifies_secondary(&controls, tolerance) {
            return None;
        }
        Some(Curve::NurbsCurve(NurbsCurve::new(
            BsplineCurve::new(KnotVector::bezier_knot(2), controls.to_vec()).ok()?,
        )))
    }
}

pub(super) fn supplied_conic(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    tolerance: f64,
) -> Option<ContactConicEvidence> {
    if let Some(conic) = ContactConicEvidence::from_carriers(carriers) {
        return measure_conic(conic, points, ids, tolerance);
    }
    let mut accepted = ContactConicEvidence::cylinder_sections(carriers, tolerance)
        .into_iter()
        .filter_map(|conic| measure_conic(conic, points, ids, tolerance));
    let first = accepted.next()?;
    // Both branches can meet at the same endpoint. With no independent
    // samples selecting one, do not invent a path between the branches.
    accepted.next().is_none().then_some(first)
}

fn measure_conic(
    mut conic: ContactConicEvidence,
    points: &[Point3],
    ids: &[usize],
    tolerance: f64,
) -> Option<ContactConicEvidence> {
    let mut sign = 0.0_f64;
    for &id in ids {
        let point = *points.get(id)?;
        let deviation = (conic.projected(point)? - point).norm();
        if deviation > tolerance {
            return None;
        }
        conic.max_vertex_residual = conic.max_vertex_residual.max(deviation);
    }
    for pair in ids.windows(2) {
        let (a, b) = (*points.get(pair[0])?, *points.get(pair[1])?);
        let step = conic.angular_step(a, b)?;
        if step == 0.0 || sign * step < 0.0 {
            return None;
        }
        let controls = conic.controls(a, b)?;
        if !conic.certifies_secondary(&controls, tolerance) {
            return None;
        }
        sign = step.signum();
        conic.signed_sweep_radians += step;
        conic.max_step_radians = conic.max_step_radians.max(step.abs());
    }
    (conic.signed_sweep_radians.abs()
        <= (f64::EPSILON * crate::usize_to_f64(ids.len())).mul_add(16.0, TAU)
        && sign != 0.0)
        .then_some(conic)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        IncludeCurve, Line, ParametricCurve, Plane, Processor, RevolutionSurface, Surface,
    };

    #[test]
    fn finite_plane_sections_are_exact_rational_quadratics()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let axis = Vector3::z();
        // Oblique cylinder ellipse, cone ellipse, parabola and hyperbola.
        for (slope, tilt) in [(0.0, 0.5), (0.5, 0.5), (1.0, 1.0), (1.0, 2.0)] {
            let plane_origin = Point3::new(0.0, 0.0, 1.0);
            let normal = Vector3::new(tilt, 0.0, 1.0).normalize();
            let carriers = [
                SurfaceCarrier::Plane {
                    origin: plane_origin,
                    normal,
                },
                SurfaceCarrier::LinearRevolution {
                    origin,
                    axis,
                    radius_at_origin: 2.0,
                    slope,
                },
            ];
            let points = [0.1_f64, 0.3, 0.5].map(|theta| {
                let (sin, cos) = theta.sin_cos();
                let z = (1.0 - 2.0 * tilt * cos) / (1.0 + slope * tilt * cos);
                let radius = 2.0 + slope * z;
                Point3::new(radius * cos, radius * sin, z)
            });
            let section =
                supplied_conic(&points, &[0, 1, 2], &carriers, 1.0e-10).ok_or("missing conic")?;
            let plane = Surface::Plane(Plane::new(
                plane_origin,
                plane_origin + Vector3::y(),
                plane_origin + normal.cross(&Vector3::y()),
            ));
            let revolution = Surface::Revolution(Processor::new(RevolutionSurface::by_revolution(
                Curve::Line(Line(
                    Point3::new(2.0 - slope * 1.9, 0.0, -1.9),
                    Point3::new(2.0 + slope * 2.0, 0.0, 2.0),
                )),
                origin,
                axis,
            )));
            for pair in points.windows(2) {
                let step = section
                    .angular_step(pair[0], pair[1])
                    .ok_or("angular step")?;
                let curve = section
                    .native_curve(pair[0], pair[1], step, 1.0e-10)
                    .ok_or("native curve")?;
                for t in [0.0, 0.125, 0.5, 0.875, 1.0] {
                    let p = curve.evaluate(t);
                    for carrier in carriers {
                        assert!(
                            carrier.signed_constraint(p).ok_or("constraint")?.1.abs() < 1.0e-10
                        );
                    }
                }
                assert!(plane.include(&curve));
                assert!(revolution.include(&curve));
            }
        }
        Ok(())
    }

    #[test]
    fn intersecting_cylinders_certify_both_bisectors_at_arbitrary_angles()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(10.0, -5.0, 7.0);
        let a = Vector3::z();
        for angle in [0.3_f64, 0.8, 1.57, 2.2] {
            let b = Vector3::new(angle.sin(), 0.0, angle.cos());
            let carriers = [
                SurfaceCarrier::LinearRevolution {
                    origin: center + a * 2.0,
                    axis: a,
                    radius_at_origin: 0.7,
                    slope: 0.0,
                },
                SurfaceCarrier::LinearRevolution {
                    origin: center - b * 3.0,
                    axis: b,
                    radius_at_origin: 0.7,
                    slope: 0.0,
                },
            ];
            for sign in [-1.0, 1.0] {
                let normal = (a + b * sign).normalize();
                let points = [0.2_f64, 0.4, 0.7].map(|theta| {
                    let x = 0.7 * theta.cos();
                    let y = 0.7 * theta.sin();
                    center + Vector3::new(x, y, -normal.x * x / normal.z)
                });
                for ordered in [carriers, [carriers[1], carriers[0]]] {
                    let section = supplied_conic(&points, &[0, 1, 2], &ordered, 1.0e-10)
                        .ok_or("missing cylinder bisector")?;
                    assert!(section.auxiliary_plane().is_some());
                    for pair in points.windows(2) {
                        let sweep = section.angular_step(pair[0], pair[1]).ok_or("step")?;
                        let curve = section
                            .native_curve(pair[0], pair[1], sweep, 1.0e-10)
                            .ok_or("native cylinder intersection")?;
                        for t in [0.0, 0.125, 0.5, 0.875, 1.0] {
                            for carrier in ordered {
                                assert!(
                                    carrier
                                        .signed_constraint(curve.evaluate(t))
                                        .ok_or("constraint")?
                                        .1
                                        .abs()
                                        < 1.0e-10
                                );
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn secondary_certificate_rejects_an_arc_with_valid_endpoints()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let cylinder = SurfaceCarrier::LinearRevolution {
            origin,
            axis: Vector3::z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let other = SurfaceCarrier::LinearRevolution {
            origin,
            axis: Vector3::x(),
            radius_at_origin: std::f64::consts::FRAC_1_SQRT_2,
            slope: 0.0,
        };
        let plane = SurfaceCarrier::Plane {
            origin,
            normal: Vector3::z(),
        };
        let mut section =
            ContactConicEvidence::from_carriers(&[cylinder, plane]).ok_or("section")?;
        section.secondary = Some(other);
        let r = std::f64::consts::FRAC_1_SQRT_2;
        let a = Point3::new(r, -r, 0.0);
        let b = Point3::new(r, r, 0.0);
        for p in [a, b] {
            assert!(other.signed_constraint(p).ok_or("constraint")?.1.abs() < 1.0e-12);
        }
        let controls = section.controls(a, b).ok_or("controls")?;
        assert!(!section.certifies_secondary(&controls, 1.0e-10));
        assert!(section.native_curve(a, b, FRAC_PI_2, 1.0e-10).is_none());
        assert!(ContactConicEvidence::cylinder_sections(&[cylinder, other], 1.0e-10).is_empty());
        Ok(())
    }

    #[test]
    fn rejects_asymptote_crossing_and_apex_degeneracy() -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let cone = SurfaceCarrier::LinearRevolution {
            origin,
            axis: Vector3::z(),
            radius_at_origin: 2.0,
            slope: 1.0,
        };
        let plane = SurfaceCarrier::Plane {
            origin: Point3::new(0.0, 0.0, 1.0),
            normal: Vector3::new(2.0, 0.0, 1.0).normalize(),
        };
        let section = ContactConicEvidence::from_carriers(&[cone, plane]).ok_or("section")?;
        // Directions straddle the pole of the rational parametrization.
        let p = |theta: f64| Point3::new(theta.cos(), theta.sin(), 0.0);
        assert!(section.controls(p(1.9), p(2.3)).is_none());
        let through_apex = SurfaceCarrier::Plane {
            origin: Point3::new(0.0, 0.0, -2.0),
            normal: Vector3::new(1.0, 0.0, 1.0).normalize(),
        };
        let degenerate = ContactConicEvidence::from_carriers(&[cone, through_apex])
            .ok_or("degenerate section")?;
        assert!(degenerate.projected(p(0.5)).is_none());
        Ok(())
    }
}
