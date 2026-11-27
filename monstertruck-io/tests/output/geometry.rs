use monstertruck_geometry::prelude::*;
use monstertruck_io::step::save::*;

fn step_test<T: StepLength>(x: T, ans: &str, length: usize)
where for<'a> StepDisplay<&'a T>: std::fmt::Display {
    let display = StepDisplay::new(&x, 1);
    assert_eq!(&display.to_string(), ans);
    assert_eq!(x.step_length(), length);
    let step = CompleteStepDisplay::new(display, Default::default()).to_string();
    step_p21::parser::parse(&step).unwrap();
}

#[test]
fn geometry() {
    step_test::<Point2>(
        Point2::new(0.0, 1.0),
        "#1 = CARTESIAN_POINT('', (0.0, 1.0));\n",
        1,
    );
    step_test::<Point3>(
        Point3::new(0.0, 1.0, 2.453),
        "#1 = CARTESIAN_POINT('', (0.0, 1.0, 2.453));\n",
        1,
    );
    step_test::<Vector2>(
        Vector2::new(3.0, 4.0),
        "#1 = VECTOR('', #2, 5.0);\n#2 = DIRECTION('', (0.6, 0.8));\n",
        2,
    );
    step_test::<Vector3>(
        Vector3::new(3.0, 4.0, 3.75),
        "#1 = VECTOR('', #2, 6.25);\n#2 = DIRECTION('', (0.48, 0.64, 0.6));\n",
        2,
    );
    step_test::<BsplineCurve<Point2>>(
        BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 1.0),
                Point2::new(2.0, 0.0),
            ],
        ),
            "\
#1 = B_SPLINE_CURVE_WITH_KNOTS('', 2, (#2, #3, #4), .UNSPECIFIED., .U., .U., (3, 3), (0.0, 1.0), .UNSPECIFIED.);
#2 = CARTESIAN_POINT('', (0.0, 0.0));
#3 = CARTESIAN_POINT('', (1.0, 1.0));
#4 = CARTESIAN_POINT('', (2.0, 0.0));\n",
    4,
    );
    step_test::<NurbsCurve<Vector3>>(
        NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector3::new(0.0, 0.0, 1.0),
                Vector3::new(1.0, 1.0, 2.0),
                Vector3::new(2.0, 0.0, 4.0),
            ],
        )),
        "\
#1 = (
    BOUNDED_CURVE()
    B_SPLINE_CURVE(2, (#2, #3, #4), .UNSPECIFIED., .U., .U.)
    B_SPLINE_CURVE_WITH_KNOTS((3, 3), (0.0, 1.0), .UNSPECIFIED.)
    CURVE()
    GEOMETRIC_REPRESENTATION_ITEM()
    RATIONAL_B_SPLINE_CURVE((1.0, 2.0, 4.0))
    REPRESENTATION_ITEM('')
);
#2 = CARTESIAN_POINT('', (0.0, 0.0));
#3 = CARTESIAN_POINT('', (0.5, 0.5));
#4 = CARTESIAN_POINT('', (0.5, 0.0));\n",
        4,
    );
    step_test::<Processor<TrimmedCurve<UnitCircle<Point2>>, Matrix3>>(
        Processor::new(TrimmedCurve::new(UnitCircle::new(), (0.0, 1.0))).transformed(
            Matrix3::from_cols(
                Vector3::new(0.0, 3.0, 0.0),
                Vector3::new(-3.0, 0.0, 0.0),
                Vector3::new(1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = CIRCLE('', #2, 3.0);
#2 = AXIS2_PLACEMENT_2D('', #3, #4);
#3 = CARTESIAN_POINT('', (1.0, 2.0));
#4 = DIRECTION('', (0.0, 1.0));\n",
        4,
    );
    step_test::<Processor<TrimmedCurve<UnitCircle<Point2>>, Matrix3>>(
        Processor::new(TrimmedCurve::new(UnitCircle::new(), (0.0, 1.0))).transformed(
            Matrix3::from_cols(
                Vector3::new(0.0, 3.0, 0.0),
                Vector3::new(-8.0, 0.0, 0.0),
                Vector3::new(1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = ELLIPSE('', #2, 3.0, 8.0);
#2 = AXIS2_PLACEMENT_2D('', #3, #4);
#3 = CARTESIAN_POINT('', (1.0, 2.0));
#4 = DIRECTION('', (0.0, 1.0));\n",
        4,
    );
    step_test::<Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>>(
        Processor::new(TrimmedCurve::new(UnitCircle::new(), (0.0, 1.0))).transformed(
            Matrix4::from_cols(
                Vector4::new(0.0, 3.0, 0.0, 0.0),
                Vector4::new(-3.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 3.0, 0.0),
                Vector4::new(3.0, 1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = CIRCLE('', #2, 3.0);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (3.0, 1.0, 2.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (0.0, 1.0, 0.0));\n",
        5,
    );
    step_test::<Processor<TrimmedCurve<UnitCircle<Point3>>, Matrix4>>(
        Processor::new(TrimmedCurve::new(UnitCircle::new(), (0.0, 1.0))).transformed(
            Matrix4::from_cols(
                Vector4::new(0.0, 3.0, 0.0, 0.0),
                Vector4::new(-8.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 1.0, 0.0),
                Vector4::new(3.0, 1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = ELLIPSE('', #2, 3.0, 8.0);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (3.0, 1.0, 2.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (0.0, 1.0, 0.0));\n",
        5,
    );
    step_test::<Processor<TrimmedCurve<UnitHyperbola<Point2>>, Matrix3>>(
        Processor::new(TrimmedCurve::new(UnitHyperbola::new(), (0.0, 1.0))).transformed(
            Matrix3::from_cols(
                Vector3::new(0.0, 3.0, 0.0),
                Vector3::new(-8.0, 0.0, 0.0),
                Vector3::new(1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = HYPERBOLA('', #2, 3.0, 8.0);
#2 = AXIS2_PLACEMENT_2D('', #3, #4);
#3 = CARTESIAN_POINT('', (1.0, 2.0));
#4 = DIRECTION('', (0.0, 1.0));\n",
        4,
    );
    step_test::<Processor<TrimmedCurve<UnitHyperbola<Point3>>, Matrix4>>(
        Processor::new(TrimmedCurve::new(UnitHyperbola::new(), (0.0, 1.0))).transformed(
            Matrix4::from_cols(
                Vector4::new(0.0, 3.0, 0.0, 0.0),
                Vector4::new(-8.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 1.0, 0.0),
                Vector4::new(3.0, 1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = HYPERBOLA('', #2, 3.0, 8.0);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (3.0, 1.0, 2.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (0.0, 1.0, 0.0));\n",
        5,
    );
    step_test::<Processor<TrimmedCurve<UnitParabola<Point2>>, Matrix3>>(
        Processor::new(TrimmedCurve::new(UnitParabola::new(), (0.0, 1.0))).transformed(
            Matrix3::from_cols(
                Vector3::new(0.0, 4.0, 0.0),
                Vector3::new(-2.0, 0.0, 0.0),
                Vector3::new(1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = PARABOLA('', #2, 1.0);
#2 = AXIS2_PLACEMENT_2D('', #3, #4);
#3 = CARTESIAN_POINT('', (1.0, 2.0));
#4 = DIRECTION('', (0.0, 1.0));\n",
        4,
    );
    step_test::<Processor<TrimmedCurve<UnitParabola<Point3>>, Matrix4>>(
        Processor::new(TrimmedCurve::new(UnitParabola::new(), (0.0, 1.0))).transformed(
            Matrix4::from_cols(
                Vector4::new(0.0, 4.0, 0.0, 0.0),
                Vector4::new(-2.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 1.0, 0.0),
                Vector4::new(3.0, 1.0, 2.0, 1.0),
            ),
        ),
        "\
#1 = PARABOLA('', #2, 1.0);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (3.0, 1.0, 2.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (0.0, 1.0, 0.0));\n",
        5,
    );
    step_test::<ParameterCurve<Line<Point2>, Plane>>(
        ParameterCurve::new(
            Line(Point2::new(0.0, 0.0), Point2::new(3.0, 4.0)),
            Plane::new(
                Point3::new(1.0, 2.0, 3.0),
                Point3::new(2.0, 2.0, 3.0),
                Point3::new(1.0, 5.0, 3.0),
            ),
        ),
        "\
#1 = PCURVE('', #2, #7);
#2 = PLANE('', #3);
#3 = AXIS2_PLACEMENT_3D('', #4, #5, #6);
#4 = CARTESIAN_POINT('', (1.0, 2.0, 3.0));
#5 = DIRECTION('', (0.0, 0.0, 1.0));
#6 = DIRECTION('', (1.0, 0.0, 0.0));
#7 = DEFINITIONAL_REPRESENTATION('', (#9), #8);
#8 = (
    GEOMETRIC_REPRESENTATION_CONTEXT(2)
    PARAMETRIC_REPRESENTATION_CONTEXT()
    REPRESENTATION_CONTEXT('2D SPACE', '')
);
#9 = LINE('', #10, #11);
#10 = CARTESIAN_POINT('', (0.0, 0.0));
#11 = VECTOR('', #12, 5.0);
#12 = DIRECTION('', (0.6, 0.8));\n",
        12,
    );
    step_test::<Plane>(
        Plane::new(
            Point3::new(1.0, 2.0, 3.0),
            Point3::new(1.0, 2.0, 4.0),
            Point3::new(2.0, 2.0, 3.0),
        ),
        "\
#1 = PLANE('', #2);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (1.0, 2.0, 3.0));
#4 = DIRECTION('', (0.0, 1.0, 0.0));
#5 = DIRECTION('', (0.0, 0.0, 1.0));\n",
        5,
    );
    step_test::<Processor<Sphere, Matrix4>>(
        Processor::new(Sphere::new(Point3::new(1.0, 2.0, 3.0), 5.0)).transformed(
            Matrix4::from_cols(
                Vector4::new(0.0, 3.0, 0.0, 0.0),
                Vector4::new(-3.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 3.0, 0.0),
                Vector4::new(2.0, 1.0, 3.0, 1.0),
            ),
        ),
        "\
#1 = SPHERICAL_SURFACE('', #2, 15.0);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (-4.0, 4.0, 12.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (0.0, 1.0, 0.0));\n",
        5,
    );
    step_test::<Processor<Torus, Matrix4>>(
        Processor::new(Torus::new_unchecked(Point3::new(1.0, 2.0, 3.0), 5.0, 3.0)).transformed(
            Matrix4::from_cols(
                Vector4::new(0.0, 3.0, 0.0, 0.0),
                Vector4::new(-3.0, 0.0, 0.0, 0.0),
                Vector4::new(0.0, 0.0, 3.0, 0.0),
                Vector4::new(2.0, 1.0, 3.0, 1.0),
            ),
        ),
        "\
#1 = TOROIDAL_SURFACE('', #2, 15.0, 9.0);
#2 = AXIS2_PLACEMENT_3D('', #3, #4, #5);
#3 = CARTESIAN_POINT('', (-4.0, 4.0, 12.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (0.0, 1.0, 0.0));\n",
        5,
    );
    step_test::<BsplineSurface<Point2>>(
        BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(2), KnotVector::uniform_knot(2, 2)),
            vec![
                vec![
                    Point2::new(0.0, 0.0),
                    Point2::new(0.0, 1.0),
                    Point2::new(0.0, 2.0),
                    Point2::new(0.0, 3.0),
                ],
                vec![
                    Point2::new(1.0, 0.0),
                    Point2::new(1.0, 1.0),
                    Point2::new(1.0, 2.0),
                    Point2::new(1.0, 3.0),
                ],
                vec![
                    Point2::new(2.0, 0.0),
                    Point2::new(2.0, 1.0),
                    Point2::new(2.0, 2.0),
                    Point2::new(2.0, 3.0),
                ],
            ],
        ),
            "\
#1 = B_SPLINE_SURFACE_WITH_KNOTS('', 2, 2, ((#2, #3, #4, #5), (#6, #7, #8, #9), (#10, #11, #12, #13)), .UNSPECIFIED., \
.U., .U., .U., (3, 3), (3, 1, 3), (0.0, 1.0), (0.0, 0.5, 1.0), .UNSPECIFIED.);
#2 = CARTESIAN_POINT('', (0.0, 0.0));
#3 = CARTESIAN_POINT('', (0.0, 1.0));
#4 = CARTESIAN_POINT('', (0.0, 2.0));
#5 = CARTESIAN_POINT('', (0.0, 3.0));
#6 = CARTESIAN_POINT('', (1.0, 0.0));
#7 = CARTESIAN_POINT('', (1.0, 1.0));
#8 = CARTESIAN_POINT('', (1.0, 2.0));
#9 = CARTESIAN_POINT('', (1.0, 3.0));
#10 = CARTESIAN_POINT('', (2.0, 0.0));
#11 = CARTESIAN_POINT('', (2.0, 1.0));
#12 = CARTESIAN_POINT('', (2.0, 2.0));
#13 = CARTESIAN_POINT('', (2.0, 3.0));\n",
    13
    );

    step_test::<NurbsSurface<Vector3>>(
        NurbsSurface::new(BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(2), KnotVector::uniform_knot(2, 2)),
            vec![
                vec![
                    Vector3::new(0.0, 0.0, 4.0),
                    Vector3::new(0.0, 1.0, 2.0),
                    Vector3::new(0.0, 2.0, 2.0),
                    Vector3::new(0.0, 3.0, 1.0),
                ],
                vec![
                    Vector3::new(1.0, 0.0, 4.0),
                    Vector3::new(1.0, 1.0, 4.0),
                    Vector3::new(1.0, 2.0, 2.0),
                    Vector3::new(1.0, 3.0, 1.0),
                ],
                vec![
                    Vector3::new(2.0, 0.0, 4.0),
                    Vector3::new(2.0, 1.0, 4.0),
                    Vector3::new(2.0, 2.0, 4.0),
                    Vector3::new(2.0, 3.0, 1.0),
                ],
            ],
        )),
            "\
#1 = (
    BOUNDED_SURFACE()
    B_SPLINE_SURFACE(2, 2, ((#2, #3, #4, #5), (#6, #7, #8, #9), (#10, #11, #12, #13)), .UNSPECIFIED., .U., .U., .U.)
    B_SPLINE_SURFACE_WITH_KNOTS((3, 3), (3, 1, 3), (0.0, 1.0), (0.0, 0.5, 1.0), .UNSPECIFIED.)
    GEOMETRIC_REPRESENTATION_ITEM()
    RATIONAL_B_SPLINE_SURFACE(((4.0, 2.0, 2.0, 1.0), (4.0, 4.0, 2.0, 1.0), (4.0, 4.0, 4.0, 1.0)))
    REPRESENTATION_ITEM('')
    SURFACE()
);
#2 = CARTESIAN_POINT('', (0.0, 0.0));
#3 = CARTESIAN_POINT('', (0.0, 0.5));
#4 = CARTESIAN_POINT('', (0.0, 1.0));
#5 = CARTESIAN_POINT('', (0.0, 3.0));
#6 = CARTESIAN_POINT('', (0.25, 0.0));
#7 = CARTESIAN_POINT('', (0.25, 0.25));
#8 = CARTESIAN_POINT('', (0.5, 1.0));
#9 = CARTESIAN_POINT('', (1.0, 3.0));
#10 = CARTESIAN_POINT('', (0.5, 0.0));
#11 = CARTESIAN_POINT('', (0.5, 0.25));
#12 = CARTESIAN_POINT('', (0.5, 0.5));
#13 = CARTESIAN_POINT('', (2.0, 3.0));\n",
    13
    );
}

#[test]
fn native_round_surface_parameter_curves_preserve_world_locus() -> anyhow::Result<()> {
    use monstertruck_io::step::load::{PcurveHolder, Table, step_geometry::StepParameterCurve};
    use monstertruck_modeling::{Curve2D as ModelCurve, Surface as ModelSurface};
    use step_p21::tables::EntityTable;

    let center = Point3::new(1.0, 2.0, 3.0);
    for sign in [1.0, -1.0] {
        let transform = Matrix4::from_cols(
            Vector4::new(0.0, 3.0, 0.0, 0.0),
            Vector4::new(-3.0 * sign, 0.0, 0.0, 0.0),
            Vector4::new(0.0, 0.0, 3.0, 0.0),
            Vector4::new(2.0, 1.0, 3.0, 1.0),
        );
        for inverted in [false, true] {
            let mut sphere = Processor::new(Sphere::new(center, 5.0)).transformed(transform);
            let mut torus = Processor::new(Torus::new(center, 5.0, 2.0)?).transformed(transform);
            if inverted {
                sphere.invert();
                torus.invert();
            }
            for (surface, entity) in [
                (ModelSurface::Sphere(sphere), "SPHERICAL_SURFACE"),
                (ModelSurface::Torus(torus), "TOROIDAL_SURFACE"),
            ] {
                let uv = Line(Point2::new(1.02, 0.6), Point2::new(1.3, 1.2));
                let native = ParameterCurve::new(ModelCurve::Line(uv), Box::new(surface));
                let text =
                    CompleteStepDisplay::new(StepDisplay::new(&native, 1), Default::default())
                        .to_string();
                assert!(text.contains(entity));
                let table = Table::from_step(&text)?;
                let owned = EntityTable::<PcurveHolder>::get_owned(&table, 1)?;
                let restored: StepParameterCurve = (&owned)
                    .try_into()
                    .map_err(|error| anyhow::anyhow!("{error}"))?;
                for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
                    assert!(
                        native.evaluate(t).near(&restored.evaluate(t)),
                        "{entity}, reflected={}, inverted={inverted}, t={t}: {:?} != {:?}",
                        sign < 0.0,
                        native.evaluate(t),
                        restored.evaluate(t)
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn cylinder_export_preserves_uv_and_oriented_normal() -> anyhow::Result<()> {
    use monstertruck_io::step::load::{PcurveHolder, Table, step_geometry::StepParameterCurve};
    use monstertruck_modeling::{
        Curve as ModelCurve, Curve2D as ModelCurve2d, Surface as ModelSurface,
    };
    use step_p21::tables::EntityTable;
    for axial_sign in [-1., 1.] {
        let generator = Line(
            Point3::new(3., 0., 2.),
            Point3::new(3., 0., 2. + axial_sign * 4.),
        );
        for reflection in [-1., 1.] {
            let transform = Matrix4::from_translation(Vector3::new(1., 5., 3.))
                * Matrix4::from_nonuniform_scale(2., 2. * reflection, 2.);
            for inverted in [false, true] {
                let mut native_surface = Processor::new(RevolutionSurface::by_revolution(
                    ModelCurve::Line(generator),
                    Point3::origin(),
                    Vector3::unit_z(),
                ))
                .transformed(transform);
                if inverted {
                    native_surface.invert();
                }
                let surface = ModelSurface::Revolution(native_surface);
                let uv = Line(Point2::new(0.2, 0.3), Point2::new(0.7, 0.8));
                let parameter_curve =
                    ParameterCurve::new(ModelCurve2d::Line(uv), Box::new(surface.clone()));
                let text = CompleteStepDisplay::new(
                    StepDisplay::new(&parameter_curve, 1),
                    Default::default(),
                )
                .to_string();
                assert!(text.contains("CYLINDRICAL_SURFACE"));
                let table = Table::from_step(&text)?;
                let owned = EntityTable::<PcurveHolder>::get_owned(&table, 1)?;
                let restored: StepParameterCurve = (&owned)
                    .try_into()
                    .map_err(|error| anyhow::anyhow!("{error}"))?;
                for t in [0., 0.25, 0.5, 0.75, 1.] {
                    let point = parameter_curve.evaluate(t);
                    assert!(point.near(&restored.evaluate(t)));
                    let radial = Vector3::new(point.x - 1., point.y - 5., 0.);
                    let uv = uv.evaluate(t);
                    assert_eq!(
                        surface.normal(uv.x, uv.y).dot(radial) > 0.,
                        surface.same_sense()
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn spline_step_output_preserves_close_distinct_knots() {
    let a = 0.25_f64;
    let b = a + 1.0e-12;
    let knots = KnotVector::from(vec![a, a, a, b, 0.5, 0.75, 1.0, 1.0, 1.0]);
    let points = (0..6)
        .map(|i| Point3::new(f64::from(i), 0.0, 0.0))
        .collect::<Vec<_>>();
    let curve = BsplineCurve::new(knots.clone(), points.clone()).unwrap();
    let rational: NurbsCurve<Vector4> = NurbsCurve::new(
        BsplineCurve::new(knots, points.iter().map(|p| p.to_homogeneous()).collect()).unwrap(),
    );
    let expected = format!(
        "(3, 1, 1, 1, 3), ({}, {}, 0.5, 0.75, 1.0)",
        FloatDisplay(a),
        FloatDisplay(b)
    );
    let polynomial_step = StepDisplay::new(&curve, 1).to_string();
    let rational_step = StepDisplay::new(&rational, 1).to_string();
    assert!(polynomial_step.contains(&expected), "{polynomial_step}");
    assert!(rational_step.contains(&expected), "{rational_step}");
}
