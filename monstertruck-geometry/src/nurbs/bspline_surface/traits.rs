use super::{
    BoundedCurve, BoundedSurface, BsplineCurve, BsplineSurface, ControlPoint, D2, Deserialize,
    EuclideanSpace, HashGen, INCLUDE_CURVE_TRIALS, IncludeCurve, Invertible, KnotVector, Matrix2,
    Matrix3, Matrix4, MetricSpace, NurbsCurve, PRESEARCH_DIVISION, ParameterDivision2D,
    ParametricCurve, ParametricSurface, ParametricSurface3D, Point2, Point3,
    SearchNearestParameter, SearchParameter, SearchParameterHint2D, SurfaceParameter,
    SurfaceParameterRange, TOLERANCE, Tolerance, Transform, Transformed, Vector2, Vector3, Vector4,
    algo,
};
use algo::surface::{SearchNearestParameterVector, SearchParameterVector};
use std::ops::Mul;

impl<P> ParameterDivision2D for BsplineSurface<P>
where P: ControlPoint<f64>
        + EuclideanSpace<Scalar = f64, Diff = <P as ControlPoint<f64>>::Diff>
        + MetricSpace<Metric = f64>
        + HashGen<f64>
{
    #[inline]
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        algo::surface::parameter_division(self, range, tol)
    }
}

impl ParametricSurface3D for BsplineSurface<Point3> {}

impl<V> BoundedSurface for BsplineSurface<V>
where Self: ParametricSurface
{
    fn range_tuple(&self) -> SurfaceParameterRange {
        let u_degree = self.knot_vecs.0.len() - self.control_points.len() - 1;
        let v_degree = self.knot_vecs.1.len() - self.control_points[0].len() - 1;
        (
            (
                self.knot_vecs.0[u_degree],
                self.knot_vecs.0[self.control_points.len()],
            ),
            (
                self.knot_vecs.1[v_degree],
                self.knot_vecs.1[self.control_points[0].len()],
            ),
        )
    }
}

// -- v2 scalar-generic impls ------------------------------------------------

use monstertruck_core::scalar::HasScalar;
use monstertruck_traits::v2;

impl<P> v2::ParametricSurface for BsplineSurface<P>
where P: HasScalar<Scalar = f64> + ControlPoint<f64>
{
    type Scalar = f64;
    type Point = P;
    type Vector = P::Diff;

    #[inline]
    fn evaluate(&self, u: Self::Scalar, v: Self::Scalar) -> P {
        ParametricSurface::evaluate(self, u, v)
    }
    #[inline]
    fn derivative_u(&self, u: Self::Scalar, v: Self::Scalar) -> P::Diff {
        ParametricSurface::derivative_u(self, u, v)
    }
    #[inline]
    fn derivative_v(&self, u: Self::Scalar, v: Self::Scalar) -> P::Diff {
        ParametricSurface::derivative_v(self, u, v)
    }
    #[inline]
    fn derivative_uu(&self, u: Self::Scalar, v: Self::Scalar) -> P::Diff {
        ParametricSurface::derivative_uu(self, u, v)
    }
    #[inline]
    fn derivative_uv(&self, u: Self::Scalar, v: Self::Scalar) -> P::Diff {
        ParametricSurface::derivative_uv(self, u, v)
    }
    #[inline]
    fn derivative_vv(&self, u: Self::Scalar, v: Self::Scalar) -> P::Diff {
        ParametricSurface::derivative_vv(self, u, v)
    }
    #[inline]
    fn period_u(&self) -> Option<Self::Scalar> { ParametricSurface::period_u(self) }
    #[inline]
    fn period_v(&self) -> Option<Self::Scalar> { ParametricSurface::period_v(self) }
}

impl<P> v2::BoundedSurface for BsplineSurface<P>
where
    P: HasScalar<Scalar = f64> + ControlPoint<f64>,
    Self: ParametricSurface + BoundedSurface,
{
    #[inline]
    fn range_tuple(&self) -> ((Self::Scalar, Self::Scalar), (Self::Scalar, Self::Scalar)) {
        BoundedSurface::range_tuple(self)
    }
}

impl v2::ParametricSurface3D for BsplineSurface<Point3> {}

impl<P> v2::SearchNearestParameter<v2::SurfaceParameter<f64>> for BsplineSurface<P>
where
    P: HasScalar<Scalar = f64>
        + ControlPoint<f64>
        + EuclideanSpace<Scalar = f64, Diff = <P as ControlPoint<f64>>::Diff>
        + MetricSpace<Metric = f64>,
    <P as ControlPoint<f64>>::Diff: SearchNearestParameterVector<Point = P>,
{
    type Point = P;
    #[inline]
    fn search_nearest_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: P,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchNearestParameter::<D2>::search_nearest_parameter(self, pt, None, trials)
    }
}

impl<P> v2::SearchParameter<v2::SurfaceParameter<f64>> for BsplineSurface<P>
where
    P: HasScalar<Scalar = f64>
        + ControlPoint<f64>
        + EuclideanSpace<Scalar = f64, Diff = <P as ControlPoint<f64>>::Diff>
        + MetricSpace<Metric = f64>
        + Tolerance,
    <P as ControlPoint<f64>>::Diff: SearchParameterVector<Point = P>,
{
    type Point = P;
    #[inline]
    fn search_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: P,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchParameter::<D2>::search_parameter(self, pt, None, trials)
    }
}

impl<V: Clone> Invertible for BsplineSurface<V> {
    #[inline]
    fn invert(&mut self) { self.swap_axes(); }
}

impl<P, V> SearchParameter<SurfaceParameter> for BsplineSurface<P>
where
    P: ControlPoint<f64, Diff = V>
        + EuclideanSpace<Scalar = f64, Diff = V>
        + MetricSpace<Metric = f64>
        + Tolerance,
    V: SearchParameterVector<Point = P>,
{
    type Point = P;
    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: P,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(x, y) => (x, y),
            SearchParameterHint2D::Range(range0, range1) => {
                algo::surface::presearch(self, point, (range0, range1), PRESEARCH_DIVISION)
            }
            SearchParameterHint2D::None => {
                algo::surface::presearch(self, point, self.range_tuple(), PRESEARCH_DIVISION)
            }
        };
        algo::surface::search_parameter(self, point, hint, trials)
    }
}

impl<P> SearchNearestParameter<SurfaceParameter> for BsplineSurface<P>
where
    P: ControlPoint<f64>
        + EuclideanSpace<Scalar = f64, Diff = <P as ControlPoint<f64>>::Diff>
        + MetricSpace<Metric = f64>,
    <P as ControlPoint<f64>>::Diff: SearchNearestParameterVector<Point = P>,
{
    type Point = P;
    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: P,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let hint = match hint.into() {
            SearchParameterHint2D::Parameter(x, y) => (x, y),
            SearchParameterHint2D::Range(range0, range1) => {
                algo::surface::presearch(self, point, (range0, range1), PRESEARCH_DIVISION)
            }
            SearchParameterHint2D::None => {
                algo::surface::presearch(self, point, self.range_tuple(), PRESEARCH_DIVISION)
            }
        };
        algo::surface::search_nearest_parameter(self, point, hint, trials)
    }
}

impl IncludeCurve<BsplineCurve<Point2>> for BsplineSurface<Point2> {
    fn include(&self, curve: &BsplineCurve<Point2>) -> bool {
        let pt = curve.front();
        let mut hint = algo::surface::presearch(self, pt, self.range_tuple(), PRESEARCH_DIVISION);
        hint = match algo::surface::search_parameter(self, pt, hint, INCLUDE_CURVE_TRIALS) {
            Some(got) => got,
            None => return false,
        };
        let knot_vector_u = self.knot_vector_u();
        let knot_vector_v = self.knot_vector_v();
        let degree = curve.degree() * 6;
        let (knots, _) = curve.knot_vector().to_single_multi();
        for i in 1..knots.len() {
            for j in 1..=degree {
                let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
                let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
                let pt = ParametricCurve::subs(curve, t);
                hint = match algo::surface::search_parameter(self, pt, hint, INCLUDE_CURVE_TRIALS) {
                    Some(got) => got,
                    None => return false,
                };
                if !ParametricSurface::subs(self, hint.0, hint.1).near(&pt)
                    || hint.0 < knot_vector_u[0] - TOLERANCE
                    || hint.0 - knot_vector_u[0] > knot_vector_u.range_length() + TOLERANCE
                    || hint.1 < knot_vector_v[0] - TOLERANCE
                    || hint.1 - knot_vector_v[0] > knot_vector_v.range_length() + TOLERANCE
                {
                    return false;
                }
            }
        }
        true
    }
}

impl IncludeCurve<BsplineCurve<Point3>> for BsplineSurface<Point3> {
    fn include(&self, curve: &BsplineCurve<Point3>) -> bool {
        let pt = curve.front();
        let mut hint = algo::surface::presearch(self, pt, self.range_tuple(), PRESEARCH_DIVISION);
        hint = match algo::surface::search_parameter(self, pt, hint, INCLUDE_CURVE_TRIALS) {
            Some(got) => got,
            None => return false,
        };
        let knot_vector_u = self.knot_vector_u();
        let knot_vector_v = self.knot_vector_v();
        let degree = curve.degree() * 6;
        let (knots, _) = curve.knot_vector().to_single_multi();
        for i in 1..knots.len() {
            for j in 1..=degree {
                let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
                let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
                let pt = ParametricCurve::subs(curve, t);
                hint = match algo::surface::search_parameter(self, pt, hint, INCLUDE_CURVE_TRIALS) {
                    Some(got) => got,
                    None => return false,
                };
                if !ParametricSurface::subs(self, hint.0, hint.1).near(&pt)
                    || hint.0 < knot_vector_u[0] - TOLERANCE
                    || hint.0 - knot_vector_u[0] > knot_vector_u.range_length() + TOLERANCE
                    || hint.1 < knot_vector_v[0] - TOLERANCE
                    || hint.1 - knot_vector_v[0] > knot_vector_v.range_length() + TOLERANCE
                {
                    return false;
                }
            }
        }
        true
    }
}

impl IncludeCurve<NurbsCurve<Vector4>> for BsplineSurface<Point3> {
    fn include(&self, curve: &NurbsCurve<Vector4>) -> bool {
        let pt = curve.subs(curve.knot_vector()[0]);
        let mut hint = algo::surface::presearch(self, pt, self.range_tuple(), PRESEARCH_DIVISION);
        hint = match algo::surface::search_parameter(self, pt, hint, INCLUDE_CURVE_TRIALS) {
            Some(got) => got,
            None => return false,
        };
        let knot_vector_u = self.knot_vector_u();
        let knot_vector_v = self.knot_vector_v();
        let degree = curve.degree() * 6;
        let (knots, _) = curve.knot_vector().to_single_multi();
        for i in 1..knots.len() {
            for j in 1..=degree {
                let p = crate::usize_to_f64(j) / crate::usize_to_f64(degree);
                let t = f64::mul_add(knots[i], p, knots[i - 1] * (1.0 - p));
                let pt = curve.subs(t);
                hint = match algo::surface::search_parameter(self, pt, hint, INCLUDE_CURVE_TRIALS) {
                    Some(got) => got,
                    None => return false,
                };
                if !ParametricSurface::subs(self, hint.0, hint.1).near(&pt)
                    || hint.0 < knot_vector_u[0] - TOLERANCE
                    || hint.0 - knot_vector_u[0] > knot_vector_u.range_length() + TOLERANCE
                    || hint.1 < knot_vector_v[0] - TOLERANCE
                    || hint.1 - knot_vector_v[0] > knot_vector_v.range_length() + TOLERANCE
                {
                    return false;
                }
            }
        }
        true
    }
}

macro_rules! impl_mat_multi {
    ($vector: ty, $matrix: ty) => {
        impl Mul<BsplineSurface<$vector>> for $matrix {
            type Output = BsplineSurface<$vector>;
            fn mul(self, mut spline: BsplineSurface<$vector>) -> Self::Output {
                spline
                    .control_points
                    .iter_mut()
                    .flat_map(|vec| vec.iter_mut())
                    .for_each(|vec| *vec = self * *vec);
                spline
            }
        }
        impl Mul<&BsplineSurface<$vector>> for $matrix {
            type Output = BsplineSurface<$vector>;
            fn mul(self, spline: &BsplineSurface<$vector>) -> Self::Output { self * spline.clone() }
        }
    };
}

macro_rules! impl_scalar_multi {
    ($vector: ty, $scalar: ty) => {
        impl_mat_multi!($vector, $scalar);
        impl Mul<$scalar> for &BsplineSurface<$vector> {
            type Output = BsplineSurface<$vector>;
            fn mul(self, scalar: $scalar) -> Self::Output { scalar * self }
        }
        impl Mul<$scalar> for BsplineSurface<$vector> {
            type Output = BsplineSurface<$vector>;
            fn mul(self, scalar: $scalar) -> Self::Output { scalar * self }
        }
    };
}

impl_mat_multi!(Vector2, Matrix2);
impl_scalar_multi!(Vector2, f64);
impl_mat_multi!(Vector3, Matrix3);
impl_scalar_multi!(Vector3, f64);
impl_mat_multi!(Vector4, Matrix4);
impl_scalar_multi!(Vector4, f64);

impl<M, P: EuclideanSpace<Scalar = f64>> Transformed<M> for BsplineSurface<P>
where M: Transform<P>
{
    #[inline]
    fn transform_by(&mut self, trans: M) {
        self.control_points
            .iter_mut()
            .flatten()
            .for_each(|p| *p = trans.transform_point(*p));
    }
}

impl<'de, P> Deserialize<'de> for BsplineSurface<P>
where P: Deserialize<'de>
{
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        #[derive(Deserialize)]
        struct BsplineSurface_<P> {
            knot_vecs: (KnotVector, KnotVector),
            control_points: Vec<Vec<P>>,
        }
        let BsplineSurface_ {
            knot_vecs,
            control_points,
        } = BsplineSurface_::<P>::deserialize(deserializer)?;
        Self::new(knot_vecs, control_points).map_err(serde::de::Error::custom)
    }
}
