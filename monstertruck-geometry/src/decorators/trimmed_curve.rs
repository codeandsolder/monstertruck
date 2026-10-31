use super::{
    Bound, BoundedCurve, CurveParameter, Cut, ParameterDivision1D, ParameterRange,
    SearchNearestParameter, SearchParameter, SearchParameterHint1D, TrimmedCurve,
};
use monstertruck_traits::ParametricCurve as ParametricCurveTrait;

impl<C> TrimmedCurve<C> {
    /// constructor
    #[inline]
    pub const fn new(curve: C, range: (f64, f64)) -> Self { Self { curve, range } }
    /// Returns the reference of non-trimmed curve
    #[inline]
    pub const fn curve(&self) -> &C { &self.curve }
    /// Returns the mutable reference of non-trimmed curve
    #[inline]
    pub const fn curve_mut(&mut self) -> &mut C { &mut self.curve }
    /// Returns the trim range as `(start, end)`.
    #[inline]
    pub const fn range(&self) -> (f64, f64) { self.range }
}

impl<C: ParametricCurveTrait> ParametricCurveTrait for TrimmedCurve<C> {
    type Point = C::Point;
    type Vector = C::Vector;
    #[inline]
    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector { self.curve.derivative_n(n, t) }
    #[inline]
    fn evaluate(&self, t: f64) -> Self::Point { self.curve.evaluate(t) }
    #[inline]
    fn derivative(&self, t: f64) -> Self::Vector { self.curve.derivative(t) }
    #[inline]
    fn derivative_2(&self, t: f64) -> Self::Vector { self.curve.derivative_2(t) }
    #[inline]
    fn period(&self) -> Option<f64> { self.curve.period() }
    #[inline]
    fn parameter_range(&self) -> ParameterRange {
        (Bound::Included(self.range.0), Bound::Included(self.range.1))
    }
}

impl<C: ParametricCurveTrait> BoundedCurve for TrimmedCurve<C> {
    fn range_tuple(&self) -> (f64, f64) { self.range }
}

impl<C: ParametricCurveTrait> Cut for TrimmedCurve<C> {
    fn cut(&mut self, t: f64) -> Self {
        let (t0, t1) = self.range;
        self.range = (t0, t);
        Self::new(self.curve.clone(), (t, t1))
    }
}

impl<C: SearchNearestParameter<CurveParameter>> SearchNearestParameter<CurveParameter>
    for TrimmedCurve<C>
{
    type Point = C::Point;
    #[inline]
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        pt: C::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        self.curve.search_nearest_parameter(pt, hint, trials)
    }
}

impl<C: SearchParameter<CurveParameter>> SearchParameter<CurveParameter> for TrimmedCurve<C> {
    type Point = C::Point;
    #[inline]
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        pt: C::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        self.curve.search_parameter(pt, hint, trials)
    }
}

impl<C: ParameterDivision1D> ParameterDivision1D for TrimmedCurve<C> {
    type Point = C::Point;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        self.curve.parameter_division(range, tol)
    }
}
