#![allow(
    clippy::similar_names,
    reason = "paired U/V ranges and divisions are standard surface-processing notation"
)]

use super::{
    BoundedCurve, BoundedSurface, CurveParameter, Cut, D1, D2, Decomposition, Deref, DerefMut,
    EuclideanSpace, IncludeCurve, InnerSpace, Invertible, Matrix3, Matrix4, MetricSpace, Mul, One,
    ParameterDivision1D, ParameterDivision2D, ParameterRange, ParametricSurface,
    ParametricSurface3D, Point2, Point3, Processor, SearchNearestParameter, SearchParameter,
    SearchParameterHint1D, SearchParameterHint2D, SquareMatrix, SurfaceParameter,
    SurfaceParameterRange, ToSameGeometry, Tolerance, Transform, Transformed, Vector3, VectorSpace,
    algo,
};
use algo::surface::SearchNearestParameterVector;
use monstertruck_traits::ParametricCurve as ParametricCurveTrait;

#[inline]
const fn undefined_vector3() -> Vector3 { Vector3::new(f64::NAN, f64::NAN, f64::NAN) }

fn matrix3_scale_bound(a: Matrix3) -> f64 {
    if let Some((_, k, _)) = a.iwasawa_decomposition() {
        return f64::abs(k[0][0])
            .max(f64::abs(k[1][1]))
            .max(f64::abs(k[2][2]))
            .max(f64::EPSILON);
    }
    let linear_norm = a[1][1]
        .mul_add(
            a[1][1],
            a[1][0].mul_add(a[1][0], a[0][0].mul_add(a[0][0], a[0][1].powi(2))),
        )
        .sqrt();
    if linear_norm.is_finite() {
        linear_norm.max(f64::EPSILON)
    } else {
        1.0
    }
}

fn matrix4_scale_bound(a: Matrix4) -> f64 {
    if let Some((_, k, _)) = a.iwasawa_decomposition() {
        let homogeneous = f64::abs(k[3][3]);
        if homogeneous > f64::EPSILON {
            return f64::abs(k[0][0])
                .max(f64::abs(k[1][1]))
                .max(f64::abs(k[2][2]))
                / homogeneous;
        }
    }
    let mut sum = 0.0;
    for column in 0..3 {
        for row in 0..3 {
            sum = a[column][row].mul_add(a[column][row], sum);
        }
    }
    let homogeneous = f64::abs(a[3][3]);
    let linear_norm = sum.sqrt();
    if linear_norm.is_finite() && homogeneous > f64::EPSILON {
        (linear_norm / homogeneous).max(f64::EPSILON)
    } else {
        1.0
    }
}

impl<E, T: One> Processor<E, T> {
    /// Creates new processor
    #[inline]
    pub fn new(entity: E) -> Self {
        Self {
            entity,
            transform: T::one(),
            orientation: true,
        }
    }

    /// Creates new transformed processor
    #[inline]
    pub const fn with_transform(entity: E, transform: T) -> Self {
        Self {
            entity,
            transform,
            orientation: true,
        }
    }

    /// Returns the reference of entity
    #[inline]
    pub const fn entity(&self) -> &E { &self.entity }

    /// Returns the reference of transform
    #[inline]
    pub const fn transform(&self) -> &T { &self.transform }

    /// Returns the orientation of surface
    #[inline]
    pub const fn orientation(&self) -> bool { self.orientation }

    #[inline]
    const fn sign(&self) -> f64 { if self.orientation { 1.0 } else { -1.0 } }

    /// apply the function to the entity geometry
    #[inline]
    pub fn map<G, F: FnOnce(E) -> G>(self, f: F) -> Processor<G, T> {
        Processor {
            entity: f(self.entity),
            transform: self.transform,
            orientation: self.orientation,
        }
    }

    /// apply the function to the entity geometry
    #[inline]
    pub fn map_ref<G, F: FnOnce(&E) -> G>(&self, f: F) -> Processor<G, T>
    where T: Copy {
        Processor {
            entity: f(&self.entity),
            transform: self.transform,
            orientation: self.orientation,
        }
    }

    /// apply the transform and inverse
    pub fn contract(self) -> E
    where E: Transformed<T> + Invertible {
        let mut res = self.entity;
        res.transform_by(self.transform);
        if !self.orientation {
            res.invert();
        }
        res
    }
}

impl<E: Clone, T: Clone> Invertible for Processor<E, T> {
    #[inline]
    fn invert(&mut self) { self.orientation = !self.orientation; }
    #[inline]
    fn inverse(&self) -> Self {
        Self {
            entity: self.entity.clone(),
            transform: self.transform.clone(),
            orientation: !self.orientation,
        }
    }
}

impl<C: BoundedCurve, T> Processor<C, T> {
    #[inline]
    fn curve_parameter(&self, t: f64) -> f64 {
        let (t0, t1) = self.range_tuple();
        if self.orientation { t } else { t0 + t1 - t }
    }
}

impl<C, T> ParametricCurveTrait for Processor<C, T>
where
    C: BoundedCurve,
    C::Point: EuclideanSpace<Diff = C::Vector>,
    C::Vector: VectorSpace<Scalar = f64>,
    T: Transform<C::Point> + Clone,
{
    type Point = C::Point;
    type Vector = C::Vector;
    #[inline]
    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector {
        if n == 0 {
            self.evaluate(t).to_vec()
        } else {
            let t = self.curve_parameter(t);
            self.transform
                .transform_vector(self.entity.derivative_n(n, t))
        }
    }
    #[inline]
    fn evaluate(&self, t: f64) -> C::Point {
        let t = self.curve_parameter(t);
        self.transform.transform_point(self.entity.evaluate(t))
    }
    #[inline]
    fn derivative(&self, t: f64) -> Self::Vector {
        let t = self.curve_parameter(t);
        self.transform.transform_vector(self.entity.derivative(t)) * self.sign()
    }
    #[inline]
    fn derivative_2(&self, t: f64) -> Self::Vector {
        let t = self.curve_parameter(t);
        self.transform.transform_vector(self.entity.derivative_2(t))
    }
    #[inline]
    fn parameter_range(&self) -> ParameterRange { self.entity.parameter_range() }
    #[inline]
    fn period(&self) -> Option<f64> { self.entity.period() }
}

impl<C, T> BoundedCurve for Processor<C, T>
where
    C: BoundedCurve,
    C::Point: EuclideanSpace<Diff = C::Vector>,
    C::Vector: VectorSpace<Scalar = f64>,
    T: Transform<C::Point> + Clone,
{
    fn range_tuple(&self) -> (f64, f64) { self.entity.range_tuple() }
}

impl<C, T> Cut for Processor<C, T>
where
    C: BoundedCurve + Cut,
    C::Point: EuclideanSpace<Diff = C::Vector>,
    C::Vector: VectorSpace<Scalar = f64>,
    T: Transform<C::Point> + Clone,
{
    fn cut(&mut self, t: f64) -> Self {
        let t = self.curve_parameter(t);
        let mut entity = self.entity.cut(t);
        if !self.orientation {
            std::mem::swap(&mut entity, &mut self.entity);
        }
        Self {
            entity,
            transform: self.transform.clone(),
            orientation: self.orientation,
        }
    }
}

impl<S, T> ParametricSurface for Processor<S, T>
where
    S: ParametricSurface,
    S::Point: EuclideanSpace<Scalar = f64, Diff = S::Vector>,
    T: Transform<S::Point> + SquareMatrix<Scalar = f64> + Clone,
{
    type Point = S::Point;
    type Vector = S::Vector;
    #[inline]
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        if (m, n) == (0, 0) {
            self.evaluate(u, v).to_vec()
        } else if self.orientation {
            self.transform
                .transform_vector(self.entity.derivative_mn(m, n, u, v))
        } else {
            self.transform
                .transform_vector(self.entity.derivative_mn(n, m, v, u))
        }
    }
    #[inline]
    fn evaluate(&self, u: f64, v: f64) -> Self::Point {
        if self.orientation {
            self.transform.transform_point(self.entity.evaluate(u, v))
        } else {
            self.transform.transform_point(self.entity.evaluate(v, u))
        }
    }
    #[inline]
    fn derivative_u(&self, u: f64, v: f64) -> Self::Vector {
        if self.orientation {
            self.transform
                .transform_vector(self.entity.derivative_u(u, v))
        } else {
            self.transform
                .transform_vector(self.entity.derivative_v(v, u))
        }
    }
    #[inline]
    fn derivative_v(&self, u: f64, v: f64) -> Self::Vector {
        if self.orientation {
            self.transform
                .transform_vector(self.entity.derivative_v(u, v))
        } else {
            self.transform
                .transform_vector(self.entity.derivative_u(v, u))
        }
    }
    #[inline]
    fn derivative_uu(&self, u: f64, v: f64) -> Self::Vector {
        if self.orientation {
            self.transform
                .transform_vector(self.entity.derivative_uu(u, v))
        } else {
            self.transform
                .transform_vector(self.entity.derivative_vv(v, u))
        }
    }
    #[inline]
    fn derivative_uv(&self, u: f64, v: f64) -> Self::Vector {
        if self.orientation {
            self.transform
                .transform_vector(self.entity.derivative_uv(u, v))
        } else {
            self.transform
                .transform_vector(self.entity.derivative_uv(v, u))
        }
    }
    #[inline]
    fn derivative_vv(&self, u: f64, v: f64) -> Self::Vector {
        if self.orientation {
            self.transform
                .transform_vector(self.entity.derivative_vv(u, v))
        } else {
            self.transform
                .transform_vector(self.entity.derivative_uu(v, u))
        }
    }
    #[inline]
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        let (urange, vrange) = self.entity.parameter_range();
        if self.orientation {
            (urange, vrange)
        } else {
            (vrange, urange)
        }
    }
    #[inline]
    fn period_u(&self) -> Option<f64> {
        if self.orientation {
            self.entity.period_u()
        } else {
            self.entity.period_v()
        }
    }
    #[inline]
    fn period_v(&self) -> Option<f64> {
        if self.orientation {
            self.entity.period_v()
        } else {
            self.entity.period_u()
        }
    }
}

impl<S, T> ParametricSurface3D for Processor<S, T>
where
    S: ParametricSurface3D,
    T: Transform<Point3> + SquareMatrix<Scalar = f64> + Clone,
{
    #[inline]
    fn normal(&self, u: f64, v: f64) -> Self::Vector {
        let transtrans = self.transform.transpose();
        let n = if self.orientation {
            self.entity.normal(u, v)
        } else {
            -self.entity.normal(v, u)
        };
        let Some(n) = transtrans.inverse_transform_vector(n) else {
            return undefined_vector3();
        };
        let determinant = self.transform.determinant();
        if determinant.abs() <= f64::EPSILON {
            return undefined_vector3();
        }
        (n / determinant).normalize()
    }
}

impl<S, T> BoundedSurface for Processor<S, T>
where
    S: BoundedSurface<Point = Point3, Vector = Vector3>,
    T: Transform<S::Point> + SquareMatrix<Scalar = f64> + Clone,
{
    fn range_tuple(&self) -> SurfaceParameterRange {
        let (u_range, v_range) = self.entity.range_tuple();
        if self.orientation {
            (u_range, v_range)
        } else {
            (v_range, u_range)
        }
    }
}

// -- v2 scalar-generic impls ------------------------------------------------

use monstertruck_traits::v2;

impl<S, T> v2::ParametricSurface for Processor<S, T>
where
    S: ParametricSurface<Point = Point3, Vector = Vector3>,
    T: Transform<Point3> + SquareMatrix<Scalar = f64> + Clone,
{
    type Scalar = f64;
    type Point = Point3;
    type Vector = Vector3;

    #[inline]
    fn evaluate(&self, u: f64, v: f64) -> Point3 { ParametricSurface::evaluate(self, u, v) }
    #[inline]
    fn derivative_u(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_u(self, u, v)
    }
    #[inline]
    fn derivative_v(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_v(self, u, v)
    }
    #[inline]
    fn derivative_uu(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_uu(self, u, v)
    }
    #[inline]
    fn derivative_uv(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_uv(self, u, v)
    }
    #[inline]
    fn derivative_vv(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_vv(self, u, v)
    }
    #[inline]
    fn period_u(&self) -> Option<f64> { ParametricSurface::period_u(self) }
    #[inline]
    fn period_v(&self) -> Option<f64> { ParametricSurface::period_v(self) }
}

impl<S, T> v2::BoundedSurface for Processor<S, T>
where
    S: BoundedSurface<Point = Point3, Vector = Vector3>,
    T: Transform<Point3> + SquareMatrix<Scalar = f64> + Clone,
{
    #[inline]
    fn range_tuple(&self) -> ((f64, f64), (f64, f64)) { BoundedSurface::range_tuple(self) }
}

impl<S, T> v2::ParametricSurface3D for Processor<S, T>
where
    S: ParametricSurface3D,
    T: Transform<Point3> + SquareMatrix<Scalar = f64> + Clone,
{
    #[inline]
    fn normal(&self, u: f64, v: f64) -> Vector3 { ParametricSurface3D::normal(self, u, v) }
}

impl<E, T> v2::SearchParameter<v2::CurveParameter<f64>> for Processor<E, T>
where
    E: SearchParameter<CurveParameter> + BoundedCurve,
    <E as SearchParameter<CurveParameter>>::Point: EuclideanSpace,
    T: Transform<<E as SearchParameter<CurveParameter>>::Point>,
{
    type Point = <E as SearchParameter<CurveParameter>>::Point;
    #[inline]
    fn search_parameter<H: Into<v2::SearchParameterHint1D<f64>>>(
        &self,
        pt: Self::Point,
        _: H,
        trials: usize,
    ) -> Option<f64> {
        SearchParameter::<D1>::search_parameter(self, pt, None, trials)
    }
}

impl<P, E, T> v2::SearchNearestParameter<v2::CurveParameter<f64>> for Processor<E, T>
where
    E: BoundedCurve<Point = P> + SearchNearestParameter<CurveParameter, Point = P>,
    P: EuclideanSpace<Scalar = f64, Diff = E::Vector>,
    E::Vector: InnerSpace<Scalar = f64> + Tolerance,
    T: Transform<P> + Clone,
{
    type Point = P;
    #[inline]
    fn search_nearest_parameter<H: Into<v2::SearchParameterHint1D<f64>>>(
        &self,
        pt: P,
        _: H,
        trials: usize,
    ) -> Option<f64> {
        SearchNearestParameter::<D1>::search_nearest_parameter(self, pt, None, trials)
    }
}

impl<E, T> v2::SearchParameter<v2::SurfaceParameter<f64>> for Processor<E, T>
where
    E: SearchParameter<SurfaceParameter>,
    E::Point: EuclideanSpace,
    T: Transform<E::Point>,
{
    type Point = E::Point;
    #[inline]
    fn search_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: E::Point,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchParameter::<D2>::search_parameter(self, pt, None, trials)
    }
}

impl<P, E, T> v2::SearchNearestParameter<v2::SurfaceParameter<f64>> for Processor<E, T>
where
    E: ParametricSurface<Point = P> + SearchNearestParameter<SurfaceParameter, Point = P>,
    P: EuclideanSpace<Scalar = f64, Diff = E::Vector> + MetricSpace<Metric = f64> + Tolerance,
    E::Vector: SearchNearestParameterVector<Point = P>,
    T: Transform<P> + SquareMatrix<Scalar = f64> + Clone,
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

impl<E, T> Deref for Processor<E, T> {
    type Target = E;
    #[inline]
    fn deref(&self) -> &E { &self.entity }
}

impl<E, T> DerefMut for Processor<E, T> {
    #[inline]
    fn deref_mut(&mut self) -> &mut E { &mut self.entity }
}

impl<E, T> Transformed<T> for Processor<E, T>
where
    T: Mul<T, Output = T> + Copy,
    E: Clone,
{
    #[inline]
    fn transform_by(&mut self, trans: T) { self.transform = trans * self.transform; }
    #[inline]
    fn transformed(&self, trans: T) -> Self {
        Self {
            entity: self.entity.clone(),
            transform: trans * self.transform,
            orientation: self.orientation,
        }
    }
}

impl<E, T, C> IncludeCurve<C> for Processor<E, T>
where
    C: ParametricCurveTrait + Transformed<T> + Clone,
    C::Point: EuclideanSpace,
    E: IncludeCurve<C>,
    T: Transform<C::Point>,
{
    fn include(&self, curve: &C) -> bool {
        let Some(inv) = self.transform.inverse_transform() else {
            return false;
        };
        let curve = curve.clone().transformed(inv);
        self.entity.include(&curve)
    }
}

impl<C> ParameterDivision1D for Processor<C, Matrix3>
where C: ParameterDivision1D<Point = Point2> + BoundedCurve<Point = Point2>
{
    type Point = Point2;

    fn try_parameter_division(
        &self,
        range: (f64, f64),
        tol: f64,
    ) -> Option<(Vec<f64>, Vec<Self::Point>)> {
        let a = self.transform;
        let range = if self.orientation {
            range
        } else {
            (self.curve_parameter(range.1), self.curve_parameter(range.0))
        };
        let (_, k, _) = a.iwasawa_decomposition()?;
        let n = f64::abs(k[0][0])
            .max(f64::abs(k[1][1]))
            .max(f64::abs(k[2][2]));
        let (mut params, mut points) = self.entity.try_parameter_division(range, tol / n)?;
        for pt in &mut points {
            *pt = a.transform_point(*pt);
        }
        if !self.orientation {
            for t in &mut params {
                *t = self.curve_parameter(*t);
            }
            params.reverse();
            points.reverse();
        }
        Some((params, points))
    }

    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        let a = self.transform;
        let range = if self.orientation {
            range
        } else {
            (self.curve_parameter(range.1), self.curve_parameter(range.0))
        };
        let n = matrix3_scale_bound(a);
        let (mut params, mut points) = self.entity.parameter_division(range, tol / n);
        for pt in &mut points {
            *pt = a.transform_point(*pt);
        }
        if !self.orientation {
            for t in &mut params {
                *t = self.curve_parameter(*t);
            }
            params.reverse();
            points.reverse();
        }
        (params, points)
    }
}

impl<C> ParameterDivision1D for Processor<C, Matrix4>
where C: ParameterDivision1D<Point = Point3> + BoundedCurve<Point = Point3>
{
    type Point = Point3;

    fn try_parameter_division(
        &self,
        range: (f64, f64),
        tol: f64,
    ) -> Option<(Vec<f64>, Vec<Self::Point>)> {
        let a = self.transform;
        let range = if self.orientation {
            range
        } else {
            (self.curve_parameter(range.1), self.curve_parameter(range.0))
        };
        let (_, k, _) = a.iwasawa_decomposition()?;
        let n = f64::abs(k[0][0])
            .max(f64::abs(k[1][1]))
            .max(f64::abs(k[2][2]))
            / f64::abs(k[3][3]);
        let (mut params, mut points) = self.entity.try_parameter_division(range, tol / n)?;
        for pt in &mut points {
            *pt = a.transform_point(*pt);
        }
        if !self.orientation {
            for t in &mut params {
                *t = self.curve_parameter(*t);
            }
            params.reverse();
            points.reverse();
        }
        Some((params, points))
    }

    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        let a = self.transform;
        let range = if self.orientation {
            range
        } else {
            (self.curve_parameter(range.1), self.curve_parameter(range.0))
        };
        let n = matrix4_scale_bound(a);
        let (mut params, mut points) = self.entity.parameter_division(range, tol / n);
        for pt in &mut points {
            *pt = a.transform_point(*pt);
        }
        if !self.orientation {
            for t in &mut params {
                *t = self.curve_parameter(*t);
            }
            params.reverse();
            points.reverse();
        }
        (params, points)
    }
}

impl<S: ParameterDivision2D> ParameterDivision2D for Processor<S, Matrix3> {
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let a = self.transform;
        let range = if self.orientation {
            range
        } else {
            (range.1, range.0)
        };
        let n = matrix3_scale_bound(a);
        let (udiv, vdiv) = self.entity.parameter_division(range, tol / n);
        if self.orientation {
            (udiv, vdiv)
        } else {
            (vdiv, udiv)
        }
    }
}

impl<S: ParameterDivision2D> ParameterDivision2D for Processor<S, Matrix4> {
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let a = self.transform;
        let range = if self.orientation {
            range
        } else {
            (range.1, range.0)
        };
        let n = matrix4_scale_bound(a);
        let (udiv, vdiv) = self.entity.parameter_division(range, tol / n);
        if self.orientation {
            (udiv, vdiv)
        } else {
            (vdiv, udiv)
        }
    }
}

impl<E, T> SearchParameter<CurveParameter> for Processor<E, T>
where
    E: SearchParameter<CurveParameter> + BoundedCurve,
    <E as SearchParameter<CurveParameter>>::Point: EuclideanSpace,
    T: Transform<<E as SearchParameter<CurveParameter>>::Point>,
{
    type Point = <E as SearchParameter<CurveParameter>>::Point;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: <E as SearchParameter<CurveParameter>>::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        let inv = self.transform.inverse_transform()?;
        let t = self
            .entity
            .search_parameter(inv.transform_point(point), hint, trials)?;
        Some(self.curve_parameter(t))
    }
}

impl<E, T> SearchParameter<SurfaceParameter> for Processor<E, T>
where
    E: SearchParameter<SurfaceParameter>,
    E::Point: EuclideanSpace,
    T: Transform<E::Point>,
{
    type Point = E::Point;
    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: E::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let inv = self.transform.inverse_transform()?;
        let (u, v) = self
            .entity
            .search_parameter(inv.transform_point(point), hint, trials)?;
        if self.orientation {
            Some((u, v))
        } else {
            Some((v, u))
        }
    }
}

impl<P, E, T> SearchNearestParameter<CurveParameter> for Processor<E, T>
where
    E: BoundedCurve<Point = P> + SearchNearestParameter<CurveParameter, Point = P>,
    P: EuclideanSpace<Scalar = f64, Diff = E::Vector>,
    E::Vector: InnerSpace<Scalar = f64> + Tolerance,
    T: Transform<P> + Clone,
{
    type Point = P;
    fn search_nearest_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        let inv = self.transform.inverse_transform()?;
        let hint =
            self.entity
                .search_nearest_parameter(inv.transform_point(point), hint, trials)?;
        let hint = self.curve_parameter(hint);
        algo::curve::search_nearest_parameter(self, point, hint, trials)
    }
}

impl<P, E, T> SearchNearestParameter<SurfaceParameter> for Processor<E, T>
where
    E: ParametricSurface<Point = P> + SearchNearestParameter<SurfaceParameter, Point = P>,
    P: EuclideanSpace<Scalar = f64, Diff = E::Vector> + MetricSpace<Metric = f64> + Tolerance,
    E::Vector: SearchNearestParameterVector<Point = P>,
    T: Transform<P> + SquareMatrix<Scalar = f64> + Clone,
{
    type Point = P;
    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let inv = self.transform.inverse_transform()?;
        let hint =
            self.entity
                .search_nearest_parameter(inv.transform_point(point), hint, trials)?;
        let hint = if self.orientation {
            hint
        } else {
            (hint.1, hint.0)
        };
        algo::surface::search_nearest_parameter(self, point, hint, trials)
    }
}

impl<E, T, U> ToSameGeometry<U> for Processor<E, T>
where
    E: ToSameGeometry<U>,
    T: Copy,
    U: Transformed<T> + Invertible,
{
    fn to_same_geometry(&self) -> U {
        let Self {
            entity,
            transform,
            orientation,
        } = self;
        let mut u = entity.to_same_geometry();
        u.transform_by(*transform);
        if !orientation {
            u.invert();
        }
        u
    }
}

#[cfg(test)]
mod tests;
