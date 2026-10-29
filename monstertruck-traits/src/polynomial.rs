#![allow(
    clippy::similar_names,
    reason = "polynomial surface composition uses canonical paired u/v component names"
)]

use crate::{
    BoundedCurve, BoundedSurface, CurveParameter, ParameterDivision1D, ParameterDivision2D,
    ParameterRange, ParametricCurve, ParametricSurface, ParametricSurface3D,
    SearchNearestParameter, SearchParameter, SearchParameterHint1D, SurfaceParameterRange, algo,
};
use monstertruck_core::{
    cgmath64::{
        ElementWise, EuclideanSpace, InnerSpace, MetricSpace, Point1, Point2, Point3, Vector1,
        Vector3, Zero,
    },
    hash::HashGen,
    tolerance::Tolerance,
};
use std::ops::Bound;

const PRESEARCH_DIVISION: usize = 50;

fn scalar_from_usize(mut value: usize) -> f64 {
    let mut result = 0.0;
    let mut bit = 1.0;
    while value != 0 {
        if value & 1 == 1 {
            result += bit;
        }
        value >>= 1;
        bit += bit;
    }
    result
}

fn falling_factorial(value: usize, order: usize) -> f64 {
    (0..order).fold(1.0, |factor, offset| {
        factor * scalar_from_usize(value - offset)
    })
}

/// polynomial curve
#[derive(Clone, Debug)]
pub struct PolynomialCurve<P: EuclideanSpace<Scalar = f64>>(pub Vec<P::Diff>);

impl<P: EuclideanSpace<Scalar = f64>> ParametricCurve for PolynomialCurve<P> {
    type Point = P;
    type Vector = P::Diff;
    fn derivative_n(&self, n: usize, t: f64) -> P::Diff {
        let iter = self.0.iter().enumerate().skip(n);
        let (_, res) = iter.fold((1.0, P::Diff::zero()), |(s, res), (deg, a)| {
            let factor = falling_factorial(deg, n);
            (s * t, res + *a * s * factor)
        });
        res
    }
    #[inline]
    fn evaluate(&self, t: f64) -> P { P::from_vec(self.derivative_n(0, t)) }
    #[inline]
    fn derivative(&self, t: f64) -> P::Diff { self.derivative_n(1, t) }
    #[inline]
    fn derivative_2(&self, t: f64) -> P::Diff { self.derivative_n(2, t) }
    #[inline]
    fn parameter_range(&self) -> ParameterRange {
        (Bound::Included(-100.0), Bound::Included(100.0))
    }
}

impl<P: EuclideanSpace<Scalar = f64>> BoundedCurve for PolynomialCurve<P> {
    fn range_tuple(&self) -> (f64, f64) { (-100.0, 100.0) }
}

impl<P> ParameterDivision1D for PolynomialCurve<P>
where P: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64> + HashGen<f64>
{
    type Point = P;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        algo::curve::parameter_division(self, range, tol)
    }
}

impl<P> SearchNearestParameter<CurveParameter> for PolynomialCurve<P>
where
    P: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64> + Tolerance,
    <P as EuclideanSpace>::Diff: InnerSpace<Scalar = f64> + Tolerance,
{
    type Point = P;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: P,
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

impl<P> SearchParameter<CurveParameter> for PolynomialCurve<P>
where
    P: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64> + Tolerance,
    <P as EuclideanSpace>::Diff: InnerSpace<Scalar = f64> + Tolerance,
{
    type Point = P;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: P,
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
        algo::curve::search_parameter(self, point, hint, trials)
    }
}

impl<P> std::ops::Add for &PolynomialCurve<P>
where
    P: EuclideanSpace<Scalar = f64>,
    P::Diff: ElementWise,
{
    type Output = PolynomialCurve<P>;
    fn add(self, rhs: Self) -> PolynomialCurve<P> {
        if self.0.len() < rhs.0.len() {
            let mut vec = rhs.0.clone();
            vec.iter_mut().zip(&self.0).for_each(|(x, &y)| *x = *x + y);
            PolynomialCurve(vec)
        } else {
            let mut vec = self.0.clone();
            vec.iter_mut().zip(&rhs.0).for_each(|(x, &y)| *x = *x + y);
            PolynomialCurve(vec)
        }
    }
}

impl<P> std::ops::Mul for &PolynomialCurve<P>
where
    P: EuclideanSpace<Scalar = f64>,
    P::Diff: ElementWise,
{
    type Output = PolynomialCurve<P>;
    fn mul(self, rhs: Self) -> PolynomialCurve<P> {
        let n = self.0.len() + rhs.0.len() - 1;
        let iter = (0..n).map(|n| {
            (0..=n).fold(P::Diff::zero(), |sum, i| {
                if i < self.0.len() && n - i < rhs.0.len() {
                    sum + self.0[i].mul_element_wise(rhs.0[n - i])
                } else {
                    sum
                }
            })
        });
        PolynomialCurve(iter.collect())
    }
}

/// polynomial surface
#[derive(Clone, Debug)]
pub struct PolynomialSurface<P: EuclideanSpace<Scalar = f64>>(pub Vec<Vec<P::Diff>>);

impl<P: EuclideanSpace<Scalar = f64>> ParametricSurface for PolynomialSurface<P> {
    type Point = P;
    type Vector = P::Diff;
    #[inline]
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        let iter = self.0.iter().enumerate().skip(m);
        let (_, res) = iter.fold((1.0, P::Diff::zero()), |(us, res), (udeg, vec)| {
            let up = falling_factorial(udeg, m);
            let iter = vec.iter().enumerate().skip(n);
            let (_, sum) = iter.fold((1.0, P::Diff::zero()), |(vs, res), (vdeg, &a)| {
                let vp = falling_factorial(vdeg, n);
                (vs * v, res + a * (vs * vp))
            });
            (us * u, res + sum * (us * up))
        });
        res
    }
    #[inline]
    fn evaluate(&self, u: f64, v: f64) -> P { P::from_vec(self.derivative_mn(0, 0, u, v)) }
    #[inline]
    fn derivative_u(&self, u: f64, v: f64) -> P::Diff { self.derivative_mn(1, 0, u, v) }
    #[inline]
    fn derivative_v(&self, u: f64, v: f64) -> P::Diff { self.derivative_mn(0, 1, u, v) }
    #[inline]
    fn derivative_uu(&self, u: f64, v: f64) -> P::Diff { self.derivative_mn(2, 0, u, v) }
    #[inline]
    fn derivative_uv(&self, u: f64, v: f64) -> P::Diff { self.derivative_mn(1, 1, u, v) }
    #[inline]
    fn derivative_vv(&self, u: f64, v: f64) -> P::Diff { self.derivative_mn(0, 2, u, v) }
    #[inline]
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        (
            (Bound::Included(-100.0), Bound::Included(100.0)),
            (Bound::Included(-50.0), Bound::Included(50.0)),
        )
    }
}

impl ParametricSurface3D for PolynomialSurface<Point3> {
    #[inline]
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        self.derivative_u(u, v)
            .cross(self.derivative_v(u, v))
            .normalize()
    }
}

impl<P: EuclideanSpace<Scalar = f64>> BoundedSurface for PolynomialSurface<P> {
    fn range_tuple(&self) -> SurfaceParameterRange { ((-100.0, 100.0), (-50.0, 50.0)) }
}

impl<P> ParameterDivision2D for PolynomialSurface<P>
where P: EuclideanSpace<Scalar = f64> + MetricSpace<Metric = f64> + HashGen<f64>
{
    fn parameter_division(&self, range: SurfaceParameterRange, tol: f64) -> (Vec<f64>, Vec<f64>) {
        algo::surface::parameter_division(self, range, tol)
    }
}

impl<P> PolynomialSurface<P>
where
    P: EuclideanSpace<Scalar = f64>,
    P::Diff: ElementWise,
{
    /// create polysurface by tensor product e.g. (u^2 + 2u + 2) * (v^2 + 3v - 2)
    #[must_use]
    pub fn by_tensor(curve0: PolynomialCurve<P>, curve1: PolynomialCurve<P>) -> Self {
        let PolynomialCurve(curve0) = curve0;
        let PolynomialCurve(curve1) = curve1;
        let vec = curve0
            .into_iter()
            .map(|p| curve1.iter().map(|&q| p.mul_element_wise(q)).collect())
            .collect();
        Self(vec)
    }

    /// composite `curve` and `self`.
    #[must_use]
    pub fn composite(&self, curve: &PolynomialCurve<Point2>) -> PolynomialCurve<P> {
        let ucurve = PolynomialCurve::<Point1>(curve.0.iter().map(|p| Vector1::new(p.x)).collect());
        let vcurve = PolynomialCurve::<Point1>(curve.0.iter().map(|p| Vector1::new(p.y)).collect());

        let mut s = PolynomialCurve::<Point1>(vec![Vector1::new(1.0)]);
        let ucurve_pow = (0..self.0.len())
            .map(|_| {
                let reserve = s.clone();
                s = &s * &ucurve;
                reserve
            })
            .collect::<Vec<_>>();
        let mut t = PolynomialCurve::<Point1>(vec![Vector1::new(1.0)]);
        let v_max = self.0.iter().map(Vec::len).max().unwrap_or(1);
        let vcurve_pow = (0..v_max)
            .map(|_| {
                let reserve = t.clone();
                t = &t * &vcurve;
                reserve
            })
            .collect::<Vec<_>>();

        let mut sum = PolynomialCurve::<P>(Vec::new());
        for (i, vec) in self.0.iter().enumerate() {
            for (j, &p) in vec.iter().enumerate() {
                let mult = &ucurve_pow[i] * &vcurve_pow[j];
                let vec = mult.0.iter().map(|a| p * a.x).collect();
                sum = &sum + &PolynomialCurve::<P>(vec);
            }
        }
        sum
    }
}
