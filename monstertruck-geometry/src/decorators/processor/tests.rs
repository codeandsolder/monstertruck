//! Unit tests for the parent module.
//!
//! Split out of the module file so the source stays readable; the module
//! path is unchanged.

use super::*;
use crate::prelude::{ParametricCurve, Sphere, TrimmedCurve, UnitCircle};
use std::f64::consts::FRAC_PI_2;

fn assert_division_samples_match_parameters<C, P>(
    curve: &C,
    (parameters, points): (Vec<f64>, Vec<P>),
) where
    C: BoundedCurve + ParameterDivision1D<Point = P> + ParametricCurve<Point = P>,
    P: Tolerance,
{
    assert_eq!(parameters.len(), points.len());
    assert!(parameters.windows(2).all(|window| window[0] <= window[1]));
    assert!(
        parameters
            .iter()
            .zip(points)
            .all(|(parameter, point)| curve.evaluate(*parameter).near(&point))
    );
}

fn assert_division_points_match_parameters<C, P>(curve: &C)
where
    C: BoundedCurve + ParameterDivision1D<Point = P> + ParametricCurve<Point = P>,
    P: Tolerance, {
    let range = curve.range_tuple();
    assert_division_samples_match_parameters(curve, curve.parameter_division(range, 0.01));
    let Some(samples) = curve.try_parameter_division(range, 0.01) else {
        panic!("fallible division should accept positive tolerance");
    };
    assert_division_samples_match_parameters(curve, samples);
}

#[test]
fn inverted_matrix3_curve_division_points_match_parameters() {
    let mut curve = Processor::new(TrimmedCurve::new(
        UnitCircle::<Point2>::new(),
        (0.0, FRAC_PI_2),
    ));
    curve.invert();

    assert_division_points_match_parameters(&curve);
}

#[test]
fn inverted_matrix4_curve_division_points_match_parameters() {
    let mut curve = Processor::new(TrimmedCurve::new(
        UnitCircle::<Point3>::new(),
        (0.0, FRAC_PI_2),
    ));
    curve.invert();

    assert_division_points_match_parameters(&curve);
}

#[test]
fn inverted_surface_bulk_derivatives_match_named_derivatives() {
    let mut surface: Processor<Sphere, Matrix4> =
        Processor::new(Sphere::new(Point3::origin(), 5.0));
    surface.invert();
    let (u, v) = (1.2, -0.3);
    let derivatives = surface.derivatives(1, u, v);
    assert!(derivatives[1][0].near(&surface.derivative_u(u, v)));
    assert!(derivatives[0][1].near(&surface.derivative_v(u, v)));
}

#[test]
fn inverted_surface_search_hints_swap_parameter_axes() {
    assert_eq!(
        processor_surface_hint(SearchParameterHint2D::Parameter(2.0, 7.0), false),
        SearchParameterHint2D::Parameter(7.0, 2.0)
    );
    assert_eq!(
        processor_surface_hint(SearchParameterHint2D::Range((1.0, 2.0), (3.0, 4.0)), false,),
        SearchParameterHint2D::Range((3.0, 4.0), (1.0, 2.0))
    );
    assert_eq!(
        processor_surface_hint(SearchParameterHint2D::Parameter(2.0, 7.0), true),
        SearchParameterHint2D::Parameter(2.0, 7.0)
    );
}
