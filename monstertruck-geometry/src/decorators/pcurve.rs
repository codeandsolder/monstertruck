use super::{
    BoundedCurve, CurveDerivatives, CurveParameter, Cut, ElementWise, EuclideanSpace, HashGen,
    InnerSpace, Invertible, MAX_DER_ORDER, MetricSpace, PRESEARCH_DIVISION, ParameterCurve,
    ParameterDivision1D, ParameterRange, ParametricCurve2D, ParametricSurface, Point2,
    SearchNearestParameter, SearchParameter, SearchParameterHint1D, SearchParameterHint2D,
    SurfaceParameter, Tolerance, Transformed, Vector2, VectorSpace, algo,
};
use monstertruck_core::cgmath64::control_point::ControlPoint;
use monstertruck_traits::ParametricCurve as ParametricCurveTrait;

impl<C, S> ParameterCurve<C, S> {
    /// Creates composited
    #[inline]
    pub const fn new(curve: C, surface: S) -> Self { Self { curve, surface } }

    /// Returns the reference to the parameter curve
    #[inline]
    pub const fn curve(&self) -> &C { &self.curve }

    /// Returns the reference to the surface
    #[inline]
    pub const fn surface(&self) -> &S { &self.surface }

    /// Decompose the pcurve into its parameter curve and surface.
    #[inline]
    pub fn decompose(self) -> (C, S) { (self.curve, self.surface) }
}

impl<C, S> ParameterCurve<C, S>
where
    C: ParametricCurve2D,
    S: ParametricSurface,
    S::Point: ControlPoint<f64, Diff = S::Vector>,
    S::Vector: VectorSpace<Scalar = f64> + ElementWise,
{
    fn der3(&self, t: f64) -> S::Vector {
        let cders = self.curve.derivatives(3, t);
        let Vector2 { x: u, y: v } = cders[0];
        let der = cders[1];
        let der2 = cders[2];
        let der3 = cders[3];
        let surface = self.surface();
        surface.derivative_mn(3, 0, u, v) * (der[0] * der[0] * der[0])
            + surface.derivative_mn(2, 1, u, v) * (der[0] * der[0] * der[1] * 3.0)
            + surface.derivative_mn(1, 2, u, v) * (der[0] * der[1] * der[1] * 3.0)
            + surface.derivative_mn(0, 3, u, v) * (der[1] * der[1] * der[1])
            + surface.derivative_uu(u, v) * (der2[0] * der[0] * 3.0)
            + surface.derivative_uv(u, v) * (der2[1].mul_add(der[0], der2[0] * der[1]) * 3.0)
            + surface.derivative_vv(u, v) * (der2[1] * der[1] * 3.0)
            + surface.derivative_u(u, v) * der3[0]
            + surface.derivative_v(u, v) * der3[1]
    }
}

impl<C, S> ParametricCurveTrait for ParameterCurve<C, S>
where
    C: ParametricCurve2D,
    S: ParametricSurface,
    S::Point: ControlPoint<f64, Diff = S::Vector>,
    S::Vector: VectorSpace<Scalar = f64> + ElementWise,
{
    type Point = S::Point;
    type Vector = S::Vector;
    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector {
        if n > MAX_DER_ORDER {
            let pt = self.curve.evaluate(t);
            return self.surface.derivative_u(pt.x, pt.y) * f64::NAN;
        }
        match n {
            0 => return self.evaluate(t).to_vec(),
            1 => return self.derivative(t),
            2 => return self.derivative_2(t),
            3 => return self.der3(t),
            _ => {}
        }
        self.derivatives(n, t)[n]
    }
    fn derivatives(&self, n: usize, t: f64) -> CurveDerivatives<Self::Vector> {
        let n = n.min(MAX_DER_ORDER);
        let cders = self.curve.derivatives(n, t);
        let Vector2 { x: u, y: v } = cders[0];
        let sders = self.surface.derivatives(n, u, v);
        sders.composite_derivatives(&cders)
    }
    #[inline]
    fn evaluate(&self, t: f64) -> Self::Point {
        let pt = self.curve.evaluate(t);
        self.surface.evaluate(pt[0], pt[1])
    }
    #[inline]
    fn derivative(&self, t: f64) -> Self::Vector {
        let cders = self.curve.derivatives(1, t);
        let pt = cders[0];
        let der = cders[1];
        self.surface.derivative_u(pt[0], pt[1]) * der[0]
            + self.surface.derivative_v(pt[0], pt[1]) * der[1]
    }
    #[inline]
    fn derivative_2(&self, t: f64) -> Self::Vector {
        let cders = self.curve.derivatives(2, t);
        let pt = cders[0];
        let der = cders[1];
        let der2 = cders[2];
        self.surface.derivative_uu(pt[0], pt[1]) * (der[0] * der[0])
            + self.surface.derivative_uv(pt[0], pt[1]) * (der[0] * der[1] * 2.0)
            + self.surface.derivative_vv(pt[0], pt[1]) * (der[1] * der[1])
            + self.surface.derivative_u(pt[0], pt[1]) * der2[0]
            + self.surface.derivative_v(pt[0], pt[1]) * der2[1]
    }
    #[inline]
    fn parameter_range(&self) -> ParameterRange { self.curve.parameter_range() }
    #[inline]
    fn period(&self) -> Option<f64> { self.curve.period() }
}

impl<C, S> BoundedCurve for ParameterCurve<C, S>
where
    C: BoundedCurve,
    Self: ParametricCurveTrait,
{
    fn range_tuple(&self) -> (f64, f64) { self.curve.range_tuple() }
}

impl<C, S> Cut for ParameterCurve<C, S>
where
    C: Cut,
    S: Clone,
    Self: ParametricCurveTrait,
{
    fn cut(&mut self, t: f64) -> Self {
        let curve = self.curve.cut(t);
        Self {
            curve,
            surface: self.surface.clone(),
        }
    }
}

impl<C, S> SearchParameter<CurveParameter> for ParameterCurve<C, S>
where
    C: ParametricCurve2D + SearchParameter<CurveParameter, Point = Point2>,
    S: SearchParameter<SurfaceParameter>,
{
    type Point = <S as SearchParameter<SurfaceParameter>>::Point;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        let hint = hint.into();
        let shint = match hint {
            SearchParameterHint1D::Parameter(hint) => {
                let p = self.curve.evaluate(hint);
                SearchParameterHint2D::Parameter(p.x, p.y)
            }
            SearchParameterHint1D::Range(x, y) => {
                let p = self.curve.evaluate(y);
                let ranges = (0..PRESEARCH_DIVISION).fold(
                    ((p.x, p.x), (p.y, p.y)),
                    |((x0, x1), (y0, y1)), i| {
                        let t = x
                            + (y - x) * crate::usize_to_f64(i)
                                / crate::usize_to_f64(PRESEARCH_DIVISION);
                        let p = self.curve.evaluate(t);
                        (
                            (f64::min(x0, p.x), f64::max(x1, p.x)),
                            (f64::min(y0, p.y), f64::max(y1, p.y)),
                        )
                    },
                );
                SearchParameterHint2D::Range(ranges.0, ranges.1)
            }
            SearchParameterHint1D::None => SearchParameterHint2D::None,
        };
        let (x, y) = self.surface.search_parameter(point, shint, trials)?;
        self.curve.search_parameter(Point2::new(x, y), hint, trials)
    }
}

impl<C, S> SearchNearestParameter<CurveParameter> for ParameterCurve<C, S>
where
    Self: BoundedCurve,
    <Self as ParametricCurveTrait>::Point: EuclideanSpace<Scalar = f64, Diff = <Self as ParametricCurveTrait>::Vector>
        + MetricSpace<Metric = f64>,
    <Self as ParametricCurveTrait>::Vector: InnerSpace<Scalar = f64> + Tolerance,
{
    type Point = <Self as ParametricCurveTrait>::Point;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        let hint = match hint.into() {
            SearchParameterHint1D::Parameter(hint) => hint,
            SearchParameterHint1D::Range(x, y) => {
                algo::curve::presearch(self, point, (x, y), PRESEARCH_DIVISION)
            }
            SearchParameterHint1D::None => {
                algo::curve::presearch(self, point, self.range_tuple(), PRESEARCH_DIVISION)
            }
        };
        algo::curve::search_nearest_parameter(self, point, hint, trials)
    }
}

impl<C, S> ParameterDivision1D for ParameterCurve<C, S>
where
    C: ParametricCurve2D,
    S: ParametricSurface,
    S::Point: EuclideanSpace<Scalar = f64, Diff = S::Vector>
        + MetricSpace<Metric = f64>
        + HashGen<f64>
        + ControlPoint<f64, Diff = S::Vector>,
    S::Vector: VectorSpace<Scalar = f64> + ElementWise,
{
    type Point = S::Point;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<S::Point>) {
        algo::curve::parameter_division(self, range, tol)
    }
}

impl<C, S> Invertible for ParameterCurve<C, S>
where
    C: Invertible,
    S: Clone,
{
    fn invert(&mut self) { self.curve.invert(); }
}

impl<C: Clone, S, T> Transformed<T> for ParameterCurve<C, S>
where S: Clone + Transformed<T>
{
    #[inline]
    fn transform_by(&mut self, trans: T) { self.surface.transform_by(trans); }
    fn transformed(&self, trans: T) -> Self {
        Self {
            curve: self.curve.clone(),
            surface: self.surface.transformed(trans),
        }
    }
}
