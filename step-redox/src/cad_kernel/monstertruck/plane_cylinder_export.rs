//! Construct the exact conic section of a structurally certified cylinder
//! and plane. The directed SSI samples select the generator or angular
//! interval; the emitted rational curve is independently checked in full.
use super::{step_polylines_match, validated_step_parameter_curve, validated_step_surface_support};
use monstertruck_core::tolerance::{
    DEFAULT_NUMERICAL_TOLERANCES, DEFAULT_STEP_RECOVERY_TOLERANCES,
};
use monstertruck_geometry::{
    analytic_carrier::AnalyticSupport,
    analytic_curve_incidence::{CurveIncidence, certify_nurbs_curve_on_support},
    analytic_cylinder_intersection::certified_circular_cylinder,
    prelude::TryIntoHomogeneousBsplineSurface,
};
use monstertruck_modeling::{
    BoundedCurve, BsplineCurve, Curve, InnerProduct, KnotVector, NurbsCurve, Plane, Point3,
    Surface, SurfaceCurve, Vector3, Vector4, algo,
};
use std::f64::consts::{FRAC_PI_2, TAU};

struct Section {
    origin: Point3,
    axis: Vector3,
    radius: f64,
    plane_origin: Point3,
    normal: Vector3,
}

impl Section {
    fn new(plane: &Plane, surface: &Surface) -> Option<Self> {
        if !surface.homogeneous_net_covers_domain() {
            return None;
        }
        let mut net = surface.try_into_homogeneous_bspline_surface()?;
        if net.degrees() == (1, 2) {
            net.swap_axes();
        }
        let (origin, axis, radius, _) = certified_circular_cylinder(&net)?;
        Some(Self {
            origin,
            axis,
            radius,
            plane_origin: plane.origin(),
            normal: plane.normal().normalize(),
        })
    }

    fn radial(&self, p: Point3) -> Vector3 {
        let delta = p - self.origin;
        delta - self.axis * delta.dot(&self.axis)
    }

    fn plane_project(&self, p: Point3) -> Point3 {
        p + self.axis * ((self.plane_origin - p).dot(&self.normal) / self.axis.dot(&self.normal))
    }

    // A degree-one line lies on an axial plane/cylinder section if its
    // whole locus follows one of the two parallel generators. Keep the
    // source vertices, and independently prove both support incidences.
    fn generator(&self, points: &[Point3], tolerance: f64) -> Option<NurbsCurve<Vector4>> {
        let start = *points.first()?;
        let end = *points.last()?;
        let length = (end - start).dot(&self.axis);
        if length.abs() <= tolerance {
            return None;
        }
        let radial = self.radial(start);
        let sign = length.signum();
        let mut previous = 0.0;
        for &p in points {
            let station = sign * (p - start).dot(&self.axis);
            if (self.radial(p) - radial).norm() > tolerance
                || station + tolerance < previous
                || station > length.abs() + tolerance
            {
                return None;
            }
            previous = station;
        }
        let controls = vec![start.to_homogeneous(), end.to_homogeneous()];
        Some(NurbsCurve::new(
            BsplineCurve::new(KnotVector::bezier_knot(1), controls).ok()?,
        ))
    }

    fn sweep(&self, points: &[Point3], tolerance: f64) -> Option<(Vector3, f64)> {
        let x = self.radial(*points.first()?).normalize();
        let mut previous = x;
        let mut sweep = 0.0_f64;
        let angular = tolerance / self.radius;
        for &point in points {
            let radial = self.radial(point);
            if (radial.norm() - self.radius).abs() > tolerance {
                return None;
            }
            let ray = radial.normalize();
            let step = self
                .axis
                .dot(&previous.cross(&ray))
                .atan2(previous.dot(&ray));
            if !step.is_finite()
                || step.abs() > FRAC_PI_2
                || (sweep.abs() > angular && step * sweep.signum() < -angular)
            {
                return None;
            }
            sweep += step;
            previous = ray;
        }
        (sweep.abs() > angular && sweep.abs() <= TAU + angular).then_some((x, sweep))
    }

    fn ellipse(&self, points: &[Point3], tolerance: f64) -> Option<NurbsCurve<Vector4>> {
        let (x, sweep) = self.sweep(points, tolerance)?;
        let y = self.axis.cross(&x);
        let project_vector =
            |v: Vector3| v - self.axis * (v.dot(&self.normal) / self.axis.dot(&self.normal));
        let curve = super::conic::affine_arc(
            self.plane_project(self.origin),
            project_vector(x * self.radius),
            project_vector(y * self.radius),
            sweep.clamp(-TAU, TAU),
        )?;
        if (curve.front() - *points.first()?).norm() > tolerance
            || (curve.back() - *points.last()?).norm() > tolerance
        {
            return None;
        }
        Some(curve)
    }

    fn certify(&self, curve: &NurbsCurve<Vector4>, tolerance: f64) -> bool {
        [
            AnalyticSupport::Plane {
                origin: self.plane_origin,
                normal: self.normal,
            },
            AnalyticSupport::Cylinder {
                origin: self.origin,
                axis: self.axis,
                radius: self.radius,
            },
        ]
        .into_iter()
        .all(|support| {
            certify_nurbs_curve_on_support(curve, support, tolerance) == CurveIncidence::Proven
        })
    }
}

pub(super) fn exact_plane_cylinder_step_guide(
    curve: &Curve,
    points: &[Point3],
    tolerance: f64,
) -> Option<Curve> {
    let Curve::Intersection(intersection) = curve else {
        return None;
    };
    let ((Surface::Plane(plane), surface) | (surface, Surface::Plane(plane))) = (
        intersection.surface0().as_ref(),
        intersection.surface1().as_ref(),
    ) else {
        return None;
    };
    let section = Section::new(plane, surface)?;
    let proof = DEFAULT_STEP_RECOVERY_TOLERANCES
        .solid_geometry_mm
        .min(tolerance);
    let axial = section.normal.dot(&section.axis).abs()
        <= DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit;
    let exact = if axial {
        section.generator(points, proof)?
    } else {
        section.ellipse(points, proof)?
    };
    if !section.certify(&exact, proof) {
        return None;
    }
    let leader = Curve::NurbsCurve(exact);
    let (_, samples) = algo::curve::parameter_division(&leader, leader.range_tuple(), tolerance);
    if !step_polylines_match(&samples, points, tolerance)
        || ![intersection.surface0(), intersection.surface1()]
            .iter()
            .all(|carrier| validated_step_surface_support(carrier, &samples, tolerance))
        || ![intersection.boundary0(), intersection.boundary1()]
            .into_iter()
            .flatten()
            .all(|boundary| validated_step_parameter_curve(boundary, &samples, tolerance))
    {
        return None;
    }
    let mut exported = SurfaceCurve::with_boundaries(
        intersection.surface0().clone(),
        intersection.surface1().clone(),
        Box::new(leader),
        intersection.boundary0().cloned(),
        intersection.boundary1().cloned(),
    );
    let (front, back) = intersection.endpoint_anchors();
    exported.set_endpoint_anchors(front, back);
    Some(Curve::Intersection(exported))
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Context, Result};
    use monstertruck_modeling::ParametricCurve;

    #[test]
    fn rational_ellipse_preserves_both_directed_branches() -> Result<()> {
        let section = Section {
            origin: Point3::new(1.0, 2.0, 3.0),
            axis: Vector3::z(),
            radius: 0.5,
            plane_origin: Point3::new(1.0, 2.0, 4.0),
            normal: Vector3::new(0.2, 0.4, 1.0).normalize(),
        };
        for sweep in [0.8, -0.8, 4.0, -4.0, TAU] {
            let trace = (0..=64)
                .map(|i| {
                    let angle = sweep * f64::from(i) / 64.0;
                    let p = section.origin
                        + Vector3::new(angle.cos(), angle.sin(), 0.0) * section.radius;
                    section.plane_project(p)
                })
                .collect::<Vec<_>>();
            let curve = section
                .ellipse(&trace, 1e-9)
                .context("missing exact section")?;
            assert!(section.certify(&curve, 1e-9));
            assert!((curve.front() - trace[0]).norm() < 1e-9);
            assert!((curve.back() - trace[64]).norm() < 1e-9);
            assert_eq!(curve.degree(), 2);
            let midpoint = curve.evaluate(0.5);
            let projected = section.radial(midpoint).normalize();
            let expected = Vector3::new((sweep * 0.5).cos(), (sweep * 0.5).sin(), 0.0);
            assert!((projected - expected).norm() < 1e-9);
        }
        Ok(())
    }

    #[test]
    fn axial_section_rejects_a_branch_switch() -> Result<()> {
        let section = Section {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::z(),
            radius: 1.0,
            plane_origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::y(),
        };
        let start = Point3::new(1.0, 0.0, 1.0);
        let end = Point3::new(1.0, 0.0, 4.0);
        let line = section
            .generator(&[start, end], 1e-9)
            .context("missing generator")?;
        assert!(section.certify(&line, 1e-9));
        assert!(
            section
                .generator(&[start, Point3::new(-1.0, 0.0, 2.0), end], 1e-9)
                .is_none()
        );
        Ok(())
    }
}
