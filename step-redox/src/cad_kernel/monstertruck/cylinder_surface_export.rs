//! Rechart exact rational cylindrical carriers for STEP exchange only.
//!
//! The original evaluated shape retains its native parameterization. All
//! exchanged spatial edge guides and their directed topology are unchanged;
//! parameter-space curves are never reused on a different chart.

use monstertruck_geometry::{
    analytic_curve_incidence::CurveIncidence,
    analytic_cylinder_intersection::certified_circular_cylinder_chart,
    prelude::TryIntoHomogeneousBsplineSurface,
};
use monstertruck_modeling::{
    Curve, Invertible, Line, Processor, RevolutionSurface, Surface,
    certify_analytic_curve_on_surface,
};

fn native_cylinder(surface: &Surface) -> Option<Surface> {
    if !matches!(
        surface,
        Surface::BsplineSurface(_) | Surface::NurbsSurface(_)
    ) || !surface.homogeneous_net_covers_domain()
    {
        return None;
    }
    let chart = certified_circular_cylinder_chart(surface.try_into_homogeneous_bspline_surface()?)?;
    let start = chart.origin + chart.radial_start * chart.radius;
    let end = start + chart.axis * chart.axial_length;
    let mut revolution = Processor::new(RevolutionSurface::by_revolution(
        Curve::Line(Line(start, end)),
        chart.origin,
        chart.axis * chart.sweep.signum(),
    ));
    // A native revolution is (axial, angle), unlike the usual tensor
    // product (angle, axial) order. Preserve the original oriented normal,
    // including reversed angular traversal and swapped parameter axes.
    if chart.angular_first {
        revolution.invert();
    }
    Some(Surface::Revolution(revolution))
}

fn guide_on_cylinder(guide: &Curve, surface: &Surface, tolerance: f64) -> bool {
    certify_analytic_curve_on_surface(guide, surface, tolerance) == CurveIncidence::Proven
}

fn spatial_guide(curve: &Curve) -> Option<&Curve> {
    match curve {
        Curve::Intersection(intersection)
            if intersection.boundary0().is_none() && intersection.boundary1().is_none() =>
        {
            Some(intersection.leader())
        }
        Curve::Intersection(_) | Curve::ParameterCurve(_) => None,
        _ => Some(curve),
    }
}

/// Rechart the STEP-only associated surfaces of a 3-D intersection guide.
/// Keep any parameterized association in its original chart.
pub(super) fn rechart_intersection(mut curve: Curve, tolerance: f64) -> Curve {
    let Curve::Intersection(intersection) = &mut curve else {
        return curve;
    };
    let candidate = |surface: &Surface| {
        native_cylinder(surface)
            .filter(|native| guide_on_cylinder(intersection.leader(), native, tolerance))
    };
    let first = intersection
        .boundary0()
        .is_none()
        .then(|| candidate(intersection.surface0()))
        .flatten();
    let second = intersection
        .boundary1()
        .is_none()
        .then(|| candidate(intersection.surface1()))
        .flatten();
    if let Some(native) = first {
        **intersection.surface0_mut() = native;
    }
    if let Some(native) = second {
        **intersection.surface1_mut() = native;
    }
    curve
}

/// Rechart a complete face with certified spatial boundaries.
///
/// Every exchanged spatial boundary must lie on the new carrier and no
/// UV trim may be reused.
pub(super) fn rechart_face(
    surface: &Surface,
    curves: impl IntoIterator<Item = Curve>,
    tolerance: f64,
) -> Option<Surface> {
    let native = native_cylinder(surface)?;
    curves
        .into_iter()
        .all(|curve| {
            spatial_guide(&curve).is_some_and(|guide| guide_on_cylinder(guide, &native, tolerance))
        })
        .then_some(native)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_modeling::{
        BsplineSurface, InnerProduct, KnotVector, NurbsSurface, ParametricSurface,
        ParametricSurface3D, Point3, Vector4,
    };

    #[test]
    fn circular_tensor_chart_preserves_both_axis_orders_and_traversals() {
        for reversed in [false, true] {
            let weight = 0.5_f64.sqrt();
            let mut profile = [
                Vector4::new(2.0, 0.0, 0.0, 1.0),
                Vector4::new(2.0 * weight, 2.0 * weight, 0.0, weight),
                Vector4::new(0.0, 2.0, 0.0, 1.0),
            ];
            if reversed {
                profile.reverse();
            }
            let controls = profile
                .into_iter()
                .map(|p| vec![p, p + Vector4::new(0.0, 0.0, p.w, 0.0)])
                .collect();
            let net = BsplineSurface::new_unchecked(
                (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
                controls,
            );
            for swapped in [false, true] {
                let mut chart = net.clone();
                if swapped {
                    chart.swap_axes();
                }
                let original = Surface::NurbsSurface(NurbsSurface::new(chart));
                let native = native_cylinder(&original).expect("single-cover circle");
                let angle = std::f64::consts::FRAC_PI_4;
                let (u, v) = if swapped { (0.5, angle) } else { (angle, 0.5) };
                assert!((original.evaluate(0.5, 0.5) - native.evaluate(u, v)).norm() < 1e-12);
                assert!(original.normal(0.5, 0.5).dot(&native.normal(u, v)) > 1.0 - 1e-12);
                let wrong =
                    Curve::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0)));
                assert!(rechart_face(&original, [wrong], 1e-7).is_none());
            }
        }
    }
}
