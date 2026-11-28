//! A regular ring-torus chart certified from a separable rational tensor net.
//!
//! The two clamped boundary curves propose the major and meridian circles.
//! Their complete knot-span incidence and directed coverage are certified;
//! the entire homogeneous net must then agree with their tensor product.
//! This is independent of degree, knot subdivision and parameter-axis order.

use crate::{
    analytic_curve_incidence::{CircularArcFrame, certify_circular_arc_coverage, propose_circle},
    analytic_cylinder_intersection::project_positive_homogeneous_point,
    prelude::*,
};

/// Geometric identity and chart orientation of a regular toroidal patch.
#[derive(Clone, Copy, Debug)]
pub struct ToroidalTensorChart {
    /// Center of the major circle.
    pub center: Point3,
    /// Unit major-circle axis.
    pub axis: Vector3,
    /// Radial reference direction.
    pub radial: Vector3,
    /// Major circle radius.
    pub major_radius: f64,
    /// Meridian circle radius.
    pub minor_radius: f64,
    /// Whether the original parameter normal points outward.
    pub outward: bool,
}

struct TensorProfiles {
    major: BsplineCurve<Vector4>,
    minor: BsplineCurve<Vector4>,
}

impl TensorProfiles {
    fn new(net: &BsplineSurface<Vector4>) -> Option<Self> {
        let controls = net.control_points();
        let first = controls.first()?;
        let major = BsplineCurve::new(
            net.knot_vector_u().clone(),
            controls.iter().map(|row| row[0]).collect(),
        )
        .ok()?;
        let minor = BsplineCurve::new(net.knot_vector_v().clone(), first.clone()).ok()?;
        Some(Self { major, minor })
    }

    fn candidate(&self, tolerance: f64) -> Option<(ToroidalTensorChart, f64)> {
        let major = propose_circle(&self.major, tolerance)?;
        let minor = propose_circle(&self.minor, tolerance)?;
        let axis = major.normal;
        let offset = minor.center - major.center;
        let center = major.center + axis * axis.dot(offset);
        let radial_vector = minor.center - center;
        let major_radius = radial_vector.magnitude();
        let minor_radius = minor.radius;
        if !major_radius.is_finite() || major_radius <= minor_radius + tolerance {
            return None;
        }
        let radial = radial_vector / major_radius;
        // The major boundary runs at the first minor parameter. It must
        // meet the meridian at its positive radial ray, not on the opposite
        // side of the torus. All remaining constraints are full-net proofs.
        let first = project_positive_homogeneous_point(self.minor.control_points()[0])?;
        let rho = (first - center).dot(radial);
        if rho <= tolerance || (rho - major.radius).abs() > tolerance {
            return None;
        }
        let major_error = tolerance / (8.0 * (major_radius + minor_radius) / major.radius).max(8.0);
        let minor_error = tolerance / 8.0;
        let major_frame =
            CircularArcFrame::new(major.center, axis, major.radius, first, major_error)?;
        let minor_frame = CircularArcFrame::new(
            minor.center,
            radial.cross(axis),
            minor_radius,
            first,
            minor_error,
        )?;
        let major_coverage =
            certify_circular_arc_coverage(&NurbsCurve::new(self.major.clone()), major_frame)?;
        let minor_coverage =
            certify_circular_arc_coverage(&NurbsCurve::new(self.minor.clone()), minor_frame)?;
        // U follows the major circle and V the meridian from radial toward
        // axial. Reversing either direction reverses the oriented normal.
        let outward = major_coverage.forward == minor_coverage.forward;
        Some((
            ToroidalTensorChart {
                center,
                axis,
                radial,
                major_radius,
                minor_radius,
                outward,
            },
            major.radius,
        ))
    }
}

/// Compare the authored tensor with the separable circular product. The
/// ideal and actual weights may differ by a global homogeneous multiplier.
/// Relative weight discrepancies perturb normalized positive averages by at
/// most delta/(1-delta) times the control-hull diameter. Spatial pole error
/// has its own convex-hull bound. Both are measured in the same length budget.
fn prove_tensor(
    net: &BsplineSurface<Vector4>,
    profiles: &TensorProfiles,
    chart: ToroidalTensorChart,
    boundary_radius: f64,
    tolerance: f64,
) -> Option<()> {
    let major = profiles.major.control_points();
    let minor = profiles.minor.control_points();
    let reference_weight = major.first()?.w;
    let major_center = chart.center
        + chart.axis
            * (project_positive_homogeneous_point(*major.first()?)? - chart.center).dot(chart.axis);
    let mut weight_error = 0.0_f64;
    let mut pole_error = 0.0_f64;
    let mut hull_radius = 0.0_f64;
    for (row, &major_control) in net.control_points().iter().zip(major) {
        let projected = project_positive_homogeneous_point(major_control)? - major_center;
        let direction = (projected - chart.axis * projected.dot(chart.axis)) / boundary_radius;
        for (&control, &minor_control) in row.iter().zip(minor) {
            let actual = project_positive_homogeneous_point(control)?;
            let meridian = project_positive_homogeneous_point(minor_control)? - chart.center;
            let rho = meridian.dot(chart.radial);
            let z = meridian.dot(chart.axis);
            let ideal_offset = direction * rho + chart.axis * z;
            let ratio = (control.w / major_control.w) / (minor_control.w / reference_weight);
            let delta = (ratio - 1.0).abs();
            if !delta.is_finite() || delta >= 1.0 {
                return None;
            }
            weight_error = weight_error.max(delta);
            pole_error = pole_error.max((actual - chart.center - ideal_offset).magnitude());
            hull_radius = hull_radius.max(ideal_offset.magnitude());
        }
    }
    let bound = (weight_error / (1.0 - weight_error)).mul_add(2.0 * hull_radius, pole_error);
    (bound.is_finite() && bound <= tolerance / 2.0).then_some(())
}

fn certify_major_first(
    net: &BsplineSurface<Vector4>,
    tolerance: f64,
) -> Option<ToroidalTensorChart> {
    let profiles = TensorProfiles::new(net)?;
    let (chart, boundary_radius) = profiles.candidate(tolerance)?;
    prove_tensor(net, &profiles, chart, boundary_radius, tolerance)?;
    Some(chart)
}

/// Prove a finite, single-cover separable rational ring-torus chart.
/// No replacement of its original parameter-space trimming curves is implied.
#[must_use = "retain a toroidal chart only when its whole-domain certificate succeeds"]
pub fn certified_toroidal_tensor(
    net: &BsplineSurface<Vector4>,
    tolerance: f64,
) -> Option<ToroidalTensorChart> {
    if !tolerance.is_finite()
        || tolerance <= 0.0
        || !net.is_clamped()
        || net.udegree() < 2
        || net.vdegree() < 2
    {
        return None;
    }
    certify_major_first(net, tolerance).or_else(|| {
        let mut transposed = net.clone();
        transposed.swap_axes();
        let mut chart = certify_major_first(&transposed, tolerance)?;
        chart.outward = !chart.outward;
        Some(chart)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    const TOL: f64 = 1.0e-7;

    fn example() -> Option<BsplineSurface<Vector4>> {
        Torus::new(Point3::new(12.0, -7.0, 3.0), 2.3, 0.4)
            .ok()?
            .try_into_homogeneous_bspline_surface()
    }

    fn assert_identity(chart: ToroidalTensorChart) {
        assert!(chart.center.distance(Point3::new(12.0, -7.0, 3.0)) < TOL);
        assert!((chart.axis.dot(Vector3::unit_z()).abs() - 1.0).abs() < TOL);
        assert!((chart.major_radius - 2.3).abs() < TOL);
        assert!((chart.minor_radius - 0.4).abs() < TOL);
    }

    #[test]
    fn full_chart_degree_elevation_subdivision_and_transposition() -> anyhow::Result<()> {
        let mut net = example().ok_or_else(|| anyhow::anyhow!("torus construction failed"))?;
        let original = certified_toroidal_tensor(&net, TOL)
            .ok_or_else(|| anyhow::anyhow!("original tensor not certified"))?;
        assert_identity(original);
        assert!(original.outward);
        net.elevate_udegree().elevate_vdegree();
        net.add_knot_u(0.13).add_knot_v(0.37);
        let elevated = certified_toroidal_tensor(&net, TOL)
            .ok_or_else(|| anyhow::anyhow!("elevated tensor not certified"))?;
        assert_identity(elevated);
        assert!(elevated.outward);
        net.swap_axes();
        let transposed = certified_toroidal_tensor(&net, TOL)
            .ok_or_else(|| anyhow::anyhow!("transposed tensor not certified"))?;
        assert_identity(transposed);
        assert!(!transposed.outward);
        Ok(())
    }

    #[test]
    fn interior_pole_and_relative_weight_errors_are_rejected() -> anyhow::Result<()> {
        let net = example().ok_or_else(|| anyhow::anyhow!("torus construction failed"))?;
        let mut points = net.control_points().clone();
        points[1][1].z += points[1][1].w * 0.001;
        let displaced = BsplineSurface::new_unchecked(net.knot_vectors().clone(), points);
        assert!(certified_toroidal_tensor(&displaced, TOL).is_none());
        let mut scaled = net.clone();
        scaled.control_points_mut().for_each(|p| *p *= 1.0e-30);
        assert!(certified_toroidal_tensor(&scaled, TOL).is_some());
        let mut points = scaled.control_points().clone();
        points[1][1] *= 1.01;
        let weighted = BsplineSurface::new_unchecked(scaled.knot_vectors().clone(), points);
        assert!(certified_toroidal_tensor(&weighted, TOL).is_none());
        Ok(())
    }
}
