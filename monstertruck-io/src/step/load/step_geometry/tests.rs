//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;

#[test]
fn projected_step_intersection_refuses_inexact_marching_guide() {
    let plane0 = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let plane1 = Surface::ElementarySurface(ElementarySurface::Plane(Plane::new(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    )));
    let guide = Curve3D::BsplineCurve(BsplineCurve::new(
        KnotVector::bezier_knot(2),
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.01),
            Point3::new(1.0, 0.0, 0.0),
        ],
    ));
    let guide_midpoint = guide.evaluate(0.5);
    assert!(guide.try_into_homogeneous_bspline_curve().is_some());
    // Unlike an intersection, a STEP surface_curve evaluates its leader.
    let direct = Curve3D::SurfaceCurve(SurfaceCurve3D::new(
        SurfaceCurveKind::SurfaceCurve,
        Box::new(guide.clone()),
        Vec::new(),
        SurfaceCurveRepresentation::Curve3D,
    ));
    assert!(direct.try_into_homogeneous_bspline_curve().is_some());
    let projected = Curve3D::IntersectionCurve(IntersectionCurve::new(
        Box::new(plane0),
        Box::new(plane1),
        Box::new(guide),
    ));
    let midpoint = projected.evaluate(0.5);
    assert!(
        (midpoint.z - guide_midpoint.z).abs() > 0.004,
        "the upstream STEP test needs a distinct marching guide"
    );
    assert!(
        projected.try_into_homogeneous_bspline_curve().is_none(),
        "the STEP projected curve must not pass off its leader as exact NURBS"
    );
}

#[test]
fn surface_curve_rejects_invalid_pcurve_on_identical_surface() {
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let leader = Curve3D::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)));
    let invalid_trim = ParameterCurve::new(
        Box::new(Curve2D::Line(Line(
            Point2::new(0.0, 1.0),
            Point2::new(1.0, 1.0),
        ))),
        Box::new(surface.clone()),
    );
    let surface_curve = SurfaceCurve3D::new(
        SurfaceCurveKind::SurfaceCurve,
        Box::new(leader),
        vec![SurfaceCurveAssociatedGeometry::ParameterCurve(invalid_trim)],
        SurfaceCurveRepresentation::Curve3D,
    );

    assert!(
        surface_curve.parameter_curve_on(&surface).is_none(),
        "Invalid face-local pcurves must not be accepted only because their surface entity matches."
    );
}

#[test]
fn public_parameter_curve_api_uses_descriptive_names() {
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let trim: StepParameterCurve = ParameterCurve::new(
        Box::new(Curve2D::Line(Line(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
        ))),
        Box::new(surface),
    );
    let curve = Curve3D::ParameterCurve(trim);

    assert!(matches!(curve, Curve3D::ParameterCurve(_)));
}
