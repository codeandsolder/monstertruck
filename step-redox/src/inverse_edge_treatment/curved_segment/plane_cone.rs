//! Bounded inverse extensions of the analytic axis-parallel plane/cone branch.
//! The supplier spline is checked over its complete active domain before the
//! reconstructed sharp-edge endpoint is allowed to move.
use super::{
    BSplineSupport, FIT_MM, GEOM_TOL_MM, KernelPatchCurve, PlaneSupport, distance, dot, sub,
};
use crate::brep::ConeSupport;
use anyhow::{Context, Result, ensure};
use monstertruck_geometry::analytic_plane_cone_hyperbola::AxisParallelPlaneConeHyperbola;
use monstertruck_geometry::analytic_plane_deviation::bound_nurbs_plane_deviation;
use monstertruck_geometry::analytic_radial_deviation::bound_nurbs_radial_deviation;
use monstertruck_geometry::prelude::{Point3, Vector3};

#[derive(Clone, Copy)]
pub(super) struct PlaneConeHyperbola {
    analytic: AxisParallelPlaneConeHyperbola,
    plane: PlaneSupport,
}

impl PlaneConeHyperbola {
    pub(super) fn new(plane: PlaneSupport, cone: ConeSupport) -> Result<Self> {
        let slope = cone.semi_angle_rad.tan();
        let axis = Vector3::from(cone.axis) * slope.signum();
        let analytic = AxisParallelPlaneConeHyperbola::from_carriers(
            Point3::from(plane.origin_mm),
            Vector3::from(plane.normal),
            Point3::from(cone.reference_origin_mm),
            axis,
            cone.reference_radius_mm,
            slope.abs(),
        )
        .context("retained plane/cone intersection is not a supported axis-parallel hyperbola")?;
        Ok(Self { analytic, plane })
    }

    fn transverse(self) -> Vector3 {
        self.analytic
            .cone_axis()
            .cross(&self.analytic.plane_normal())
    }

    pub(super) fn fitted_to_source(
        self,
        spline: &BSplineSupport,
        ends: [[f64; 3]; 2],
    ) -> Result<()> {
        let analytic = self.analytic;
        let curve =
            crate::source_edge_locus::complete_bspline_support(spline.clone(), ends, FIT_MM)?;
        // Positive rational weights and a strictly monotone transverse control
        // projection certify the same directed branch over the complete domain.
        let transverse: [f64; 3] = self.transverse().into();
        let progress = dot(sub(ends[1], ends[0]), transverse);
        let direction = progress.signum();
        let ordered = |sign: f64| {
            spline
                .control_points_mm
                .windows(2)
                .all(|pair| sign * dot(sub(pair[1], pair[0]), transverse) > GEOM_TOL_MM)
        };
        ensure!(
            progress.abs() > GEOM_TOL_MM && (ordered(direction) || ordered(-direction)),
            "retained plane/cone spline has no monotone transverse trim"
        );
        let normal = analytic.plane_normal();
        let origin = analytic.reference_origin();
        let plane = bound_nurbs_plane_deviation(
            &curve,
            origin + normal * analytic.signed_plane_offset(),
            normal,
            GEOM_TOL_MM,
        )
        .context("complete source hyperbola plane bound unsupported")?;
        let cone = bound_nurbs_radial_deviation(
            &curve,
            origin,
            analytic.cone_axis(),
            analytic.reference_radius(),
            analytic.slope(),
            GEOM_TOL_MM,
        )
        .context("complete source hyperbola cone bound unsupported")?;
        // The same transverse coordinate selects one unique analytic point.
        // Moving the normal coordinate by <= plane.upper changes radius by
        // <= plane.upper. The remaining axial displacement is bounded by
        // the radial cone residual divided by the cone slope.
        let fit = plane
            .upper
            .hypot((plane.upper + cone.radial_upper) / analytic.slope());
        ensure!(
            fit <= FIT_MM,
            "retained spline exceeds the existing plane/cone feature fit: {fit:e} mm"
        );
        ensure!(
            analytic
                .oriented_trim(Point3::from(ends[0]), Point3::from(ends[1]), FIT_MM)
                .is_some(),
            "source endpoints do not belong to the same analytic hyperbola branch"
        );
        Ok(())
    }

    pub(super) fn extend_to(self, endpoint: [f64; 3], direction: [f64; 3]) -> Result<[f64; 3]> {
        let analytic = self.analytic;
        let origin: [f64; 3] = analytic.reference_origin().into();
        let axis: [f64; 3] = analytic.cone_axis().into();
        let normal: [f64; 3] = analytic.plane_normal().into();
        let transverse: [f64; 3] = self.transverse().into();
        let relative = sub(endpoint, origin);
        ensure!(
            (dot(relative, normal) - analytic.signed_plane_offset()).abs() <= GEOM_TOL_MM
                && dot(direction, normal).abs() <= super::DIR_TOL,
            "virtual edge is not coplanar with the hyperbola"
        );
        let radius = analytic
            .slope()
            .mul_add(dot(relative, axis), analytic.reference_radius());
        let radius_rate = analytic.slope() * dot(direction, axis);
        let lateral = dot(relative, transverse);
        let lateral_rate = dot(direction, transverse);
        let offset = analytic.signed_plane_offset();
        let quadratic = radius_rate.mul_add(radius_rate, -lateral_rate * lateral_rate);
        let linear = 2.0 * radius.mul_add(radius_rate, -lateral * lateral_rate);
        let constant = radius.mul_add(radius, -lateral.mul_add(lateral, offset * offset));
        let mut roots = Vec::with_capacity(2);
        if quadratic == 0.0 {
            ensure!(
                linear != 0.0,
                "virtual line has no isolated hyperbola intersection"
            );
            roots.push(-constant / linear);
        } else {
            let discriminant = linear.mul_add(linear, -4.0 * quadratic * constant);
            ensure!(discriminant >= 0.0, "virtual line misses the hyperbola");
            let stable = -0.5 * (linear + discriminant.sqrt().copysign(linear));
            if stable == 0.0 {
                roots.push(-linear / (2.0 * quadratic));
            } else {
                roots.extend([stable / quadratic, constant / stable]);
            }
        }
        let mut candidates = roots
            .into_iter()
            .filter(|root| root.is_finite())
            .map(|root| super::add(endpoint, super::mul(direction, root)))
            .filter(|point| {
                analytic
                    .parameter_at_vertex(Point3::from(*point), GEOM_TOL_MM)
                    .is_some()
            })
            .collect::<Vec<_>>();
        candidates.sort_by(|a, b| distance(*a, endpoint).total_cmp(&distance(*b, endpoint)));
        let first = *candidates
            .first()
            .context("virtual edge has no intersection on the positive-radius hyperbola")?;
        let extension = distance(first, endpoint);
        ensure!(
            extension <= FIT_MM,
            "nearest positive-radius hyperbola intersection is {extension:e} mm from the inferred endpoint (fit budget {FIT_MM:e} mm)"
        );
        candidates.retain(|point| distance(*point, endpoint) <= FIT_MM);
        ensure!(
            candidates
                .iter()
                .all(|point| distance(*point, first) <= GEOM_TOL_MM),
            "virtual edge has ambiguous nearby hyperbola intersections"
        );
        Ok(first)
    }

    pub(super) fn trimmed(self, first: [f64; 3], second: [f64; 3]) -> Result<KernelPatchCurve> {
        let analytic = self.analytic;
        let trim = analytic
            .oriented_trim(Point3::from(first), Point3::from(second), GEOM_TOL_MM)
            .context("retained hyperbola endpoints do not lie on the exact branch")?;
        let residual = analytic
            .source_plane_residual_upper(
                Point3::from(self.plane.origin_mm),
                Vector3::from(self.plane.normal),
                trim,
            )
            .context("cannot bound the complete reconstructed trim against its source plane")?;
        ensure!(
            residual <= GEOM_TOL_MM,
            "reconstructed hyperbola leaves its original plane"
        );
        let offset = analytic.signed_plane_offset();
        let axis = analytic.cone_axis();
        let center = analytic.reference_origin() + analytic.plane_normal() * offset
            - axis * (analytic.reference_radius() / analytic.slope());
        Ok(KernelPatchCurve::HyperbolaArc {
            center_mm: center.into(),
            cosh_axis_mm: (axis * (offset.abs() / analytic.slope())).into(),
            sinh_axis_mm: (self.transverse() * offset.abs()).into(),
            start_parameter: trim.start_parameter,
            end_parameter: trim.end_parameter,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn conic() -> Result<PlaneConeHyperbola> {
        PlaneConeHyperbola::new(
            PlaneSupport {
                origin_mm: [0.5, 0.0, 0.0],
                normal: [1.0, 0.0, 0.0],
                max_residual_mm: 0.0,
            },
            ConeSupport {
                reference_origin_mm: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                x_direction: [1.0, 0.0, 0.0],
                reference_radius_mm: 1.0,
                semi_angle_rad: std::f64::consts::FRAC_PI_4,
            },
        )
    }

    #[test]
    fn hyperbola_inverse_line_reaches_exact_positive_branch_and_retains_reverse_trim() -> Result<()>
    {
        let conic = conic()?;
        let original: [f64; 3] = conic
            .analytic
            .evaluate(0.45)
            .context("invalid test point")?
            .into();
        let shifted = [original[0], original[1], original[2] + FIT_MM / 4.0];
        let restored = conic.extend_to(shifted, [0.0, 0.0, 1.0])?;
        assert!(distance(restored, original) <= GEOM_TOL_MM);
        let last = conic
            .analytic
            .evaluate(-0.25)
            .context("invalid test end")?
            .into();
        let KernelPatchCurve::HyperbolaArc {
            start_parameter,
            end_parameter,
            ..
        } = conic.trimmed(restored, last)?
        else {
            anyhow::bail!("retained conic lost its native hyperbola type");
        };
        assert!((start_parameter - 0.45).abs() <= GEOM_TOL_MM);
        assert!((end_parameter + 0.25).abs() <= GEOM_TOL_MM);
        assert!(
            conic
                .extend_to([shifted[0] + 0.1, shifted[1], shifted[2]], [0.0, 0.0, 1.0])
                .is_err()
        );
        assert!(
            conic
                .extend_to([shifted[0], shifted[1], shifted[2] + 1.0], [0.0, 0.0, 1.0])
                .is_err()
        );
        Ok(())
    }
}
