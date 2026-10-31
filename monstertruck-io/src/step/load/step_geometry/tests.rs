//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;

#[test]
fn surface_curve_rejects_invalid_pcurve_on_identical_surface() {
    let surface = Surface::Elementary(ElementarySurface::Plane(Plane::xy()));
    let leader = Curve3D::Line(Line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)));
    let invalid_trim = ParameterCurve::new(
        Box::new(Curve2D::Line(Line(
            Point2::new(0.0, 1.0),
            Point2::new(1.0, 1.0),
        ))),
        Box::new(surface.clone()),
    );
    let surface_curve = SurfaceCurve3D::new(
        SurfaceCurveKind::Surface,
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
fn sampled_same_surface_match_does_not_enter_exact_reparameterization() {
    let source = Surface::Bspline(BsplineSurface::from(Plane::xy()));
    let target = Surface::Elementary(ElementarySurface::Plane(Plane::xy()));
    assert!(SurfaceCurve3D::same_surface(&source, &target));
    assert_ne!(source, target);

    let leader = Curve3D::Line(Line(
        Point3::new(0.0, 0.25, 0.0),
        Point3::new(1.0, 0.25, 0.0),
    ));
    let source_trim = ParameterCurve::new(
        Box::new(Curve2D::Line(Line(
            Point2::new(0.0, 0.25),
            Point2::new(1.0, 0.25),
        ))),
        Box::new(source),
    );
    let surface_curve = SurfaceCurve3D::new(
        SurfaceCurveKind::Surface,
        Box::new(leader),
        vec![SurfaceCurveAssociatedGeometry::ParameterCurve(source_trim)],
        SurfaceCurveRepresentation::ParameterCurve0,
    );

    assert!(
        surface_curve.parameter_curve_on(&target).is_some(),
        "legacy lookup deliberately retains sampled same-surface compatibility",
    );
    assert!(
        surface_curve
            .parameter_curve_reparameterized_on(&target)
            .is_none(),
        "proof-bearing reparameterization must not inherit sampled same_surface() acceptance",
    );
}

#[test]
fn public_parameter_curve_api_uses_descriptive_names() {
    let surface = Surface::Elementary(ElementarySurface::Plane(Plane::xy()));
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

#[test]
fn surface_curve_reuses_pcurve_across_equivalent_cylinder_charts() {
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
    let source = Surface::Elementary(ElementarySurface::CylindricalSurface(source_cylinder));
    let target = Surface::Elementary(ElementarySurface::CylindricalSurface(target_cylinder));

    let leader = Curve3D::Line(Line(
        Point3::new(2.0, 0.0, -3.0),
        Point3::new(2.0, 0.0, 4.0),
    ));
    let source_trim = ParameterCurve::new(
        Box::new(Curve2D::Line(Line(
            Point2::new(0.0, -3.0),
            Point2::new(0.0, 4.0),
        ))),
        Box::new(source),
    );
    let surface_curve = SurfaceCurve3D::new(
        SurfaceCurveKind::Surface,
        Box::new(leader.clone()),
        vec![SurfaceCurveAssociatedGeometry::ParameterCurve(source_trim)],
        SurfaceCurveRepresentation::ParameterCurve0,
    );

    let remapped = surface_curve
        .parameter_curve_reparameterized_on(&target)
        .expect("equivalent cylinder charts must transport the stored pcurve");
    assert_eq!(remapped.surface().as_ref(), &target);
    let (t0, t1) = remapped.curve().range_tuple();
    for fraction in [0.0, 0.2, 0.5, 0.8, 1.0] {
        let t = t0 + (t1 - t0) * fraction;
        let uv = remapped.curve().evaluate(t);
        let point = target.evaluate(uv.x, uv.y);
        assert!(leader.search_parameter(point, None, 100).is_some());
    }
}
