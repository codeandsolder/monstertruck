//! Native modeling adapter for global analytic curve-on-carrier certificates.

use super::surface::analytic_carrier;
use super::{BsplineCurve, Curve, KnotVector, Line, Surface};
use monstertruck_geometry::analytic_curve_incidence::{
    CurveIncidence, certify_bspline_curve_on_support, certify_nurbs_curve_on_support,
};

/// Non-sampling global incidence proof for an entire native curve on a surface.
///
/// A match does not certify face trims, edge incidence or occupied solid
/// equivalence. Unsupported UV parameter compositions are never guessed.
#[must_use]
pub fn certify_analytic_curve_on_surface(
    curve: &Curve,
    surface: &Surface,
    tolerance: f64,
) -> CurveIncidence {
    let Some(carrier) = analytic_carrier(surface) else {
        return CurveIncidence::Unsupported;
    };
    match curve {
        Curve::Line(line) => {
            let Line(front, back) = *line;
            let spline = BsplineCurve::new_unchecked(KnotVector::bezier_knot(1), vec![front, back]);
            certify_bspline_curve_on_support(&spline, carrier.support(), tolerance)
        }
        Curve::BsplineCurve(spline) => {
            certify_bspline_curve_on_support(spline, carrier.support(), tolerance)
        }
        Curve::NurbsCurve(spline) => {
            certify_nurbs_curve_on_support(spline, carrier.support(), tolerance)
        }
        // IntersectionCurve::evaluate runs a Newton projection onto its two
        // source surfaces (with a leader fallback). The evaluated curve is NOT
        // structurally identical to its leader, even when the leader itself
        // is globally on the requested support. Do not promote a leader
        // certificate into a claim about the projected curve.
        Curve::Intersection(_) => CurveIncidence::Unsupported,
        Curve::ParameterCurve(_) => CurveIncidence::Unsupported,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Plane, Point3, Processor, RevolutionSurface, Vector3};
    use monstertruck_geometry::prelude::{KnotVector, NurbsCurve, Vector4};
    const TOL: f64 = 1.0e-7;

    #[test]
    fn native_surface_adapter_proves_planar_line_and_rejects_offset() {
        let plane = Surface::Plane(Plane::xy());
        let exact = Curve::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)));
        let above = Curve::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.1)));
        assert_eq!(
            certify_analytic_curve_on_surface(&exact, &plane, TOL),
            CurveIncidence::Proven
        );
        assert_eq!(
            certify_analytic_curve_on_surface(&above, &plane, TOL),
            CurveIncidence::NotProven
        );
    }

    #[test]
    fn projected_intersection_cannot_inherit_a_leader_only_certificate() {
        use crate::SurfaceCurve;
        let xy = Surface::Plane(Plane::xy());
        let yz = Surface::Plane(Plane::yz());
        let leader = Curve::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)));
        assert_eq!(
            certify_analytic_curve_on_surface(&leader, &xy, TOL),
            CurveIncidence::Proven,
        );
        let projected = Curve::Intersection(SurfaceCurve::with_boundaries(
            Box::new(xy.clone()),
            Box::new(yz),
            Box::new(leader),
            None,
            None,
        ));
        assert_eq!(
            certify_analytic_curve_on_surface(&projected, &xy, TOL),
            CurveIncidence::Unsupported,
        );
    }

    #[test]
    fn native_cylinder_accepts_exact_rational_quarter_arc_not_perturbation() {
        let revolution = RevolutionSurface::by_revolution(
            Curve::Line(Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 1.0))),
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
        );
        let surface = Surface::Revolution(Processor::new(revolution));
        let weight = 0.5_f64.sqrt();
        let rational = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(2.0, 0.0, 0.0, 1.0),
                Vector4::new(2.0 * weight, 2.0 * weight, 0.0, weight),
                Vector4::new(0.0, 2.0, 0.0, 1.0),
            ],
        )));
        assert_eq!(
            certify_analytic_curve_on_surface(&rational, &surface, TOL),
            CurveIncidence::Proven
        );
        let mut wrong = rational;
        if let Curve::NurbsCurve(curve) = &mut wrong
            && let Some(middle) = curve.control_points_mut().nth(1)
        {
            middle.x += 0.01;
        }
        assert_eq!(
            certify_analytic_curve_on_surface(&wrong, &surface, TOL),
            CurveIncidence::NotProven
        );
    }
}
