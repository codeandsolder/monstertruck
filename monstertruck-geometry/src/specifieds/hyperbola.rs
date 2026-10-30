#![allow(
    clippy::many_single_char_names,
    reason = "a through d are polynomial coefficients in the analytic hyperbola solve"
)]

use super::{
    CurveParameter, D1, EuclideanSpace, HashGen, MetricSpace, Origin, ParameterDivision1D,
    ParametricCurve, Point2, Point3, SearchNearestParameter, SearchParameter,
    SearchParameterHint1D, TOLERANCE, UnitHyperbola, Vector2, Vector3, algo, solver,
};

impl<P> UnitHyperbola<P> {
    /// constructor
    #[inline]
    #[must_use]
    pub const fn new() -> Self { Self(std::marker::PhantomData) }
}

impl ParametricCurve for UnitHyperbola<Point2> {
    type Point = Point2;
    type Vector = Vector2;
    #[inline]
    fn derivative_n(&self, n: usize, t: f64) -> Vector2 {
        match n % 2 {
            0 => Vector2::new(f64::cosh(t), f64::sinh(t)),
            _ => Vector2::new(f64::sinh(t), f64::cosh(t)),
        }
    }
    #[inline]
    fn evaluate(&self, t: f64) -> Self::Point { Point2::from_vec(self.derivative_n(0, t)) }
    #[inline]
    fn derivative(&self, t: f64) -> Self::Vector { self.derivative_n(1, t) }
    #[inline]
    fn derivative_2(&self, t: f64) -> Self::Vector { self.derivative_n(2, t) }
}

impl ParametricCurve for UnitHyperbola<Point3> {
    type Point = Point3;
    type Vector = Vector3;
    #[inline]
    fn derivative_n(&self, n: usize, t: f64) -> Vector3 {
        match n % 2 {
            0 => Vector3::new(f64::cosh(t), f64::sinh(t), 0.0),
            _ => Vector3::new(f64::sinh(t), f64::cosh(t), 0.0),
        }
    }
    #[inline]
    fn evaluate(&self, t: f64) -> Self::Point { Point3::from_vec(self.derivative_n(0, t)) }
    #[inline]
    fn derivative(&self, t: f64) -> Self::Vector { self.derivative_n(1, t) }
    #[inline]
    fn derivative_2(&self, t: f64) -> Self::Vector { self.derivative_n(2, t) }
}

// -- v2 scalar-generic impls ------------------------------------------------

use monstertruck_traits::v2;

macro_rules! impl_v2_hyperbola {
    ($point:ty, $vector:ty) => {
        impl v2::ParametricCurve for UnitHyperbola<$point> {
            type Scalar = f64;
            type Point = $point;
            type Vector = $vector;

            #[inline]
            fn evaluate(&self, t: f64) -> $point { ParametricCurve::evaluate(self, t) }
            #[inline]
            fn derivative(&self, t: f64) -> $vector { ParametricCurve::derivative(self, t) }
            #[inline]
            fn derivative_2(&self, t: f64) -> $vector { ParametricCurve::derivative_2(self, t) }
            #[inline]
            fn derivative_n(&self, n: usize, t: f64) -> $vector {
                ParametricCurve::derivative_n(self, n, t)
            }
            #[inline]
            fn period(&self) -> Option<f64> { ParametricCurve::period(self) }
            #[inline]
            fn try_range_tuple(&self) -> Option<(f64, f64)> {
                ParametricCurve::try_range_tuple(self)
            }
        }
    };
}

impl_v2_hyperbola!(Point2, Vector2);
impl_v2_hyperbola!(Point3, Vector3);

macro_rules! impl_v2_hyperbola_search {
    ($point:ty) => {
        impl v2::SearchNearestParameter<v2::CurveParameter<f64>> for UnitHyperbola<$point> {
            type Point = $point;
            #[inline]
            fn search_nearest_parameter<H: Into<v2::SearchParameterHint1D<f64>>>(
                &self,
                pt: $point,
                _: H,
                _: usize,
            ) -> Option<f64> {
                SearchNearestParameter::<D1>::search_nearest_parameter(self, pt, None, 0)
            }
        }

        impl v2::SearchParameter<v2::CurveParameter<f64>> for UnitHyperbola<$point> {
            type Point = $point;
            #[inline]
            fn search_parameter<H: Into<v2::SearchParameterHint1D<f64>>>(
                &self,
                pt: $point,
                _: H,
                _: usize,
            ) -> Option<f64> {
                SearchParameter::<D1>::search_parameter(self, pt, None, 0)
            }
        }
    };
}

impl_v2_hyperbola_search!(Point2);
impl_v2_hyperbola_search!(Point3);

impl<P> ParameterDivision1D for UnitHyperbola<P>
where
    Self: ParametricCurve<Point = P>,
    P: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64> + HashGen<f64>,
{
    type Point = P;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<P>) {
        algo::curve::parameter_division(self, range, tol)
    }
}

impl SearchNearestParameter<CurveParameter> for UnitHyperbola<Point2> {
    type Point = Point2;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        p: Point2,
        _: H,
        _: usize,
    ) -> Option<f64> {
        let a = -p.y;
        let b = p.x.mul_add(-p.x, p.y * p.y) / 4.0 + 1.0;
        let c = -p.y;
        let d = p.y * p.y / 4.0;
        let y = solver::solve_quartic(a, b, c, d)
            .into_iter()
            .filter_map(|z| if z.im.so_small() { Some(z.re) } else { None })
            .min_by(|s, t| {
                p.distance2(self.subs(*s))
                    .total_cmp(&p.distance2(self.subs(*t)))
            })?;
        Some(f64::asinh(y))
    }
}

impl SearchNearestParameter<CurveParameter> for UnitHyperbola<Point3> {
    type Point = Point3;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        p: Point3,
        _: H,
        trials: usize,
    ) -> Option<f64> {
        UnitHyperbola::<Point2>::new().search_nearest_parameter(Point2::new(p.x, p.y), None, trials)
    }
}

impl SearchParameter<CurveParameter> for UnitHyperbola<Point2> {
    type Point = Point2;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        p: Point2,
        _: H,
        _: usize,
    ) -> Option<f64> {
        // Verify that p lies on the unit hyperbola (cosh(t), sinh(t)).
        // The naive check `p.near(&self.subs(asinh(p.y)))` fails for large
        // |t| because the absolute difference between two independently
        // computed cosh values can vastly exceed the fixed tolerance, even
        // when the relative error is at machine-epsilon level.
        //
        // Instead we check: x must be positive, and x must equal
        // sqrt(1 + y^2) to within a relative tolerance.
        if p.x < -TOLERANCE {
            return None;
        }
        let expected_x = f64::hypot(1.0, p.y);
        let scale = f64::max(1.0, expected_x);
        if (p.x - expected_x).abs() > TOLERANCE * scale {
            return None;
        }
        Some(f64::asinh(p.y))
    }
}

impl SearchParameter<CurveParameter> for UnitHyperbola<Point3> {
    type Point = Point3;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        p: Point3,
        _: H,
        _: usize,
    ) -> Option<f64> {
        if !p.z.so_small() {
            return None;
        }
        UnitHyperbola::<Point2>::new().search_parameter(Point2::new(p.x, p.y), None, 0)
    }
}
