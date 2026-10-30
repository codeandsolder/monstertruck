#![allow(
    clippy::similar_names,
    reason = "paired U/V parameter-division names are standard surface notation"
)]

use super::{
    Angle, BoundedCurve, BoundedSurface, BsplineCurve, CurveDerivatives, CurveParameter, Debug,
    EuclideanSpace, HashGen, InnerSpace, KnotVector, MAX_DER_ORDER, Matrix, Matrix2, Matrix3,
    Matrix4, MetricSpace, Mul, NurbsCurve, ParameterDivision1D, ParameterDivision2D,
    ParameterRange, ParametricCurve, ParametricCurve3D, ParametricSurface, ParametricSurface3D,
    Point1, Point2, Point3, Rad, RollingBallFilletContactCurve, RollingBallFilletSurface,
    SearchNearestParameter, SearchParameter, SearchParameterHint1D, SearchParameterHint2D,
    SquareMatrix, SurfaceDerivatives, SurfaceParameter, SurfaceParameterRange, TOLERANCE,
    ToSameGeometry, Tolerance, Vector2, Vector3, Vector4, Zero, nonpositive_tolerance,
};
use monstertruck_traits::ParametricCurve as ParametricCurveTrait;
use std::f64::consts::PI;

impl<C, S0, S1, R> RollingBallFilletSurface<C, S0, S1, R> {
    /// constructor
    #[inline]
    pub const fn new(edge_curve: C, surface0: S0, surface1: S1, radius: R) -> Self {
        Self {
            edge_curve,
            surface0,
            surface1,
            radius,
        }
    }

    /// returns edge curve
    #[inline]
    pub const fn edge_curve(&self) -> &C { &self.edge_curve }
    /// returns first surface
    #[inline]
    pub const fn surface0(&self) -> &S0 { &self.surface0 }
    /// returns second surface
    #[inline]
    pub const fn surface1(&self) -> &S1 { &self.surface1 }
    /// returns radius function
    #[inline]
    pub const fn radius(&self) -> &R { &self.radius }

    /// returns the orbit curve of contact point with `surface0`.
    #[inline]
    pub fn contact_curve0(&self) -> RollingBallFilletContactCurve<C, S0, S1, R>
    where Self: Clone {
        RollingBallFilletContactCurve {
            surface: self.clone(),
            index: 0,
        }
    }
    /// returns the orbit curve of contact point with `surface1`.
    #[inline]
    pub fn contact_curve1(&self) -> RollingBallFilletContactCurve<C, S0, S1, R>
    where Self: Clone {
        RollingBallFilletContactCurve {
            surface: self.clone(),
            index: 1,
        }
    }
}

/// Trait for radius functions.
pub trait RadiusFunction: Clone {
    /// Returns the `n`-th derivative at parameter `t`.
    fn derivative_n(&self, n: usize, t: f64) -> f64;
    /// Evaluates the radius at parameter `t`.
    #[inline]
    fn evaluate(&self, t: f64) -> f64 { self.derivative_n(0, t) }
    /// Returns the first derivative at parameter `t`.
    #[inline]
    fn derivative(&self, t: f64) -> f64 { self.derivative_n(1, t) }
    /// Returns the second derivative at parameter `t`.
    #[inline]
    fn derivative_2(&self, t: f64) -> f64 { self.derivative_n(2, t) }
    /// Returns all derivatives at parameter `t` with order `0..=max_order`.
    #[inline]
    fn derivatives(&self, max_order: usize, t: f64) -> CurveDerivatives<f64> {
        let mut derivatives = CurveDerivatives::new(max_order);
        for (order, derivative) in derivatives.iter_mut().enumerate() {
            *derivative = self.derivative_n(order, t);
        }
        derivatives
    }
}

impl RadiusFunction for f64 {
    #[inline]
    fn derivative_n(&self, n: usize, _: f64) -> f64 {
        match n {
            0 => *self,
            _ => 0.0,
        }
    }
}

macro_rules! impl_radius_1dim {
    ($ty: ty) => {
        impl RadiusFunction for $ty {
            #[inline]
            fn derivative_n(&self, n: usize, t: f64) -> f64 {
                ParametricCurveTrait::derivative_n(self, n, t).x
            }
        }
    };
}
impl_radius_1dim!(BsplineCurve<Point1>);
impl_radius_1dim!(NurbsCurve<Vector2>);

/// Contact point of the rolling ball and surface
#[derive(Clone, Copy, Debug)]
pub struct ContactPoint {
    /// the 3d-coordinate of contact point
    pub point: Point3,
    /// the parameter on the surface
    pub uv: Point2,
}

impl From<(Point3, Point2)> for ContactPoint {
    #[inline]
    fn from((point, uv): (Point3, Point2)) -> Self { Self { point, uv } }
}

impl From<ContactPoint> for (Point3, (f64, f64)) {
    #[inline]
    fn from(cp: ContactPoint) -> Self { (cp.point, (cp.uv.x, cp.uv.y)) }
}

/// Contact circle for rolling ball fillet.
#[derive(Clone, Copy, Debug)]
pub struct ContactCircle {
    center: Point3,
    axis: Vector3,
    angle: Rad<f64>,
    t: f64,
    contact_point0: ContactPoint,
    contact_point1: ContactPoint,
}

mod algo;
mod contact_circle;

#[inline]
const fn undefined_vector() -> Vector3 { Vector3::new(f64::NAN, f64::NAN, f64::NAN) }

#[inline]
const fn undefined_point() -> Point3 { Point3::new(f64::NAN, f64::NAN, f64::NAN) }

// `ParametricCurve`/`ParametricSurface` are legacy infallible traits. The
// rolling-ball solve is genuinely fallible at singular contact configurations,
// so only those trait adapters translate `None` into NaN geometry. New direct
// callers should use the `try_*` methods below.
fn undefined_curve_derivatives(max_order: usize) -> CurveDerivatives<Vector3> {
    let mut derivatives = CurveDerivatives::new(max_order);
    for order in 0..=derivatives.max_order() {
        derivatives[order] = undefined_vector();
    }
    derivatives
}

fn undefined_surface_derivatives(max_order: usize) -> SurfaceDerivatives<Vector3> {
    let mut derivatives = SurfaceDerivatives::new(max_order);
    for row in derivatives.slice_iter_mut() {
        for derivative in row {
            *derivative = undefined_vector();
        }
    }
    derivatives
}

impl<C, S0, S1, R> RollingBallFilletSurface<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    fn sub_der_mn(&self, m: usize, n: usize, u: f64, cc: ContactCircle) -> Option<Vector3> {
        if n == 0 {
            return Some(cc.derivative_n(m, u));
        }
        Some(self.vder_info(cc, n)?.derivative_mn(m, n, u))
    }

    /// Evaluates the fillet surface when a valid rolling contact circle exists.
    /// Returns `None` at singular or otherwise unsolved contact configurations.
    #[inline]
    pub fn try_evaluate(&self, u: f64, v: f64) -> Option<Point3> {
        self.contact_circle(v).map(|cc| cc.evaluate(u))
    }

    /// Evaluates one mixed derivative when the rolling-contact solve is regular.
    /// Returns `None` if the contact or derivative solve is singular.
    #[inline]
    pub fn try_derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Option<Vector3> {
        self.contact_circle(v)
            .and_then(|cc| self.sub_der_mn(m, n, u, cc))
    }

    /// Evaluates all representable derivatives up to `max_order`.
    ///
    /// The fixed derivative container caps the order at [`MAX_DER_ORDER`].
    /// Returns `None` if any requested derivative crosses a singular rolling-contact state.
    pub fn try_derivatives(
        &self,
        max_order: usize,
        u: f64,
        v: f64,
    ) -> Option<SurfaceDerivatives<Vector3>> {
        let cc = self.contact_circle(v)?;
        let mut out = SurfaceDerivatives::new(max_order);
        for (i, row) in out.slice_iter_mut().enumerate() {
            for (j, derivative) in row.iter_mut().enumerate() {
                *derivative = self.sub_der_mn(i, j, u, cc)?;
            }
        }
        Some(out)
    }
}

impl<C, S0, S1, R> ParametricSurface for RollingBallFilletSurface<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    type Point = Point3;
    type Vector = Vector3;
    fn derivatives(&self, max_order: usize, u: f64, v: f64) -> SurfaceDerivatives<Vector3> {
        self.try_derivatives(max_order, u, v)
            .unwrap_or_else(|| undefined_surface_derivatives(max_order))
    }
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Vector3 {
        self.try_derivative_mn(m, n, u, v)
            .unwrap_or_else(undefined_vector)
    }
    fn evaluate(&self, u: f64, v: f64) -> Point3 {
        self.try_evaluate(u, v).unwrap_or_else(undefined_point)
    }
    fn derivative_u(&self, u: f64, v: f64) -> Vector3 {
        self.contact_circle(v)
            .map_or_else(undefined_vector, |cc| cc.derivative(u))
    }
    fn derivative_v(&self, u: f64, v: f64) -> Vector3 {
        self.contact_circle(v)
            .and_then(|cc| self.vder_info(cc, 1))
            .map_or_else(undefined_vector, |info| info.derivative_v(u))
    }
    fn derivative_uu(&self, u: f64, v: f64) -> Self::Vector {
        self.contact_circle(v)
            .map_or_else(undefined_vector, |cc| cc.derivative_2(u))
    }
    fn derivative_uv(&self, u: f64, v: f64) -> Self::Vector {
        self.contact_circle(v)
            .and_then(|cc| self.vder_info(cc, 1))
            .map_or_else(undefined_vector, |info| info.derivative_uv(u))
    }
    fn derivative_vv(&self, u: f64, v: f64) -> Self::Vector {
        self.contact_circle(v)
            .and_then(|cc| self.vder_info(cc, 2))
            .map_or_else(undefined_vector, |info| info.derivative_vv(u))
    }
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        use std::ops::Bound::Included;
        (
            (Included(0.0), Included(1.0)),
            self.edge_curve.parameter_range(),
        )
    }
    fn period_v(&self) -> Option<f64> { self.edge_curve.period() }
}

impl<C, S0, S1, R> ParametricSurface3D for RollingBallFilletSurface<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
}

impl<C, S0, S1, R> BoundedSurface for RollingBallFilletSurface<C, S0, S1, R>
where
    C: ParametricCurve3D + BoundedCurve,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    fn range_tuple(&self) -> SurfaceParameterRange { ((0.0, 1.0), self.edge_curve.range_tuple()) }
}

impl<C, S0, S1, R> ParameterDivision2D for RollingBallFilletSurface<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        nonpositive_tolerance!(tol);
        let udiv = self
            .parameter_division_u(range, tol)
            .unwrap_or_else(|| vec![range.0.0, range.0.1]);
        let mut vdiv = vec![range.1.0, range.1.1];
        algo::parameter_division_v_for_fillet(self, &udiv, &mut vdiv, tol);
        (udiv, vdiv)
    }
}

impl<C, S0, S1, R> SearchParameter<SurfaceParameter> for RollingBallFilletSurface<C, S0, S1, R>
where
    C: ParametricCurve3D + SearchNearestParameter<CurveParameter, Point = Point3>,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    type Point = Point3;
    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        let curve_hint = match hint.into() {
            SearchParameterHint2D::Parameter(_, v) => SearchParameterHint1D::Parameter(v),
            SearchParameterHint2D::Range(_, (v0, v1)) => SearchParameterHint1D::Range(v0, v1),
            SearchParameterHint2D::None => SearchParameterHint1D::None,
        };
        let edge_curve = &self.edge_curve;
        let v = edge_curve.search_nearest_parameter(point, curve_hint, trials)?;
        let cc = self.contact_circle(v)?;

        let cp0 = cc.contact_point0.point - cc.center;
        let cp = point - cc.center;
        let u = cp.angle(cp0).0 / cc.angle.0;
        if cp.magnitude2().near(&cp0.magnitude2()) {
            Some((u, v))
        } else {
            None
        }
    }
}

impl<C, S0, S1, R> RollingBallFilletContactCurve<C, S0, S1, R> {
    /// original fillet surface
    #[inline]
    pub const fn fillet_surface(&self) -> &RollingBallFilletSurface<C, S0, S1, R> { &self.surface }
    /// curve index: curve on `surface0` => 0, curve on `surface1` => 1.
    #[inline]
    pub const fn index(&self) -> usize { self.index }
}

impl<C, S0, S1, R> RollingBallFilletContactCurve<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    /// Evaluates the contact point when the rolling-contact solve is regular.
    /// Returns `None` at singular or otherwise unsolved contact configurations.
    pub fn try_evaluate(&self, t: f64) -> Option<Point3> {
        let cc = self.surface.contact_circle(t)?;
        Some(match self.index {
            0 => cc.contact_point0.point,
            _ => cc.contact_point1.point,
        })
    }

    /// Evaluates all representable derivatives up to `max_order`.
    ///
    /// The fixed derivative container caps the order at [`MAX_DER_ORDER`].
    /// Returns `None` if the rolling-contact derivative solve is singular.
    pub fn try_derivatives(&self, max_order: usize, t: f64) -> Option<CurveDerivatives<Vector3>> {
        let max_order = max_order.min(MAX_DER_ORDER);
        if max_order == 0 {
            let mut derivatives = CurveDerivatives::new(0);
            derivatives[0] = self.try_evaluate(t)?.to_vec();
            return Some(derivatives);
        }
        let cc = self.surface.contact_circle(t)?;
        let rders = self.surface.radius.derivatives(max_order, t);
        let cc_ders = self
            .surface
            .sub_center_contacts_ders(cc, &rders, max_order)?;
        Some(match self.index {
            0 => cc_ders.contact0,
            _ => cc_ders.contact1,
        })
    }

    /// Evaluates one derivative when its order is representable and the solve is regular.
    /// Returns `None` for orders above [`MAX_DER_ORDER`] or a singular contact solve.
    pub fn try_derivative_n(&self, n: usize, t: f64) -> Option<Vector3> {
        if n > MAX_DER_ORDER {
            return None;
        }
        Some(self.try_derivatives(n, t)?[n])
    }
}

impl<C, S0, S1, R> ParametricCurveTrait for RollingBallFilletContactCurve<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    type Point = Point3;
    type Vector = Vector3;
    fn derivatives(&self, n: usize, t: f64) -> CurveDerivatives<Vector3> {
        self.try_derivatives(n, t)
            .unwrap_or_else(|| undefined_curve_derivatives(n))
    }
    fn derivative_n(&self, n: usize, t: f64) -> Self::Vector {
        self.try_derivative_n(n, t).unwrap_or_else(undefined_vector)
    }
    fn evaluate(&self, t: f64) -> Self::Point {
        self.try_evaluate(t).unwrap_or_else(undefined_point)
    }
    fn derivative(&self, t: f64) -> Self::Vector { self.derivative_n(1, t) }
    fn derivative_2(&self, t: f64) -> Self::Vector { self.derivative_n(2, t) }
    #[inline]
    fn parameter_range(&self) -> ParameterRange { self.surface.edge_curve.parameter_range() }
    #[inline]
    fn period(&self) -> Option<f64> { self.surface.edge_curve.period() }
}

impl<C, S0, S1, R> BoundedCurve for RollingBallFilletContactCurve<C, S0, S1, R>
where
    C: ParametricCurve3D + BoundedCurve,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    fn range_tuple(&self) -> (f64, f64) { self.surface.edge_curve.range_tuple() }
}

impl<C, S0, S1, R> ParameterDivision1D for RollingBallFilletContactCurve<C, S0, S1, R>
where
    C: ParametricCurve3D,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    type Point = Point3;
    fn parameter_division(&self, range: (f64, f64), tol: f64) -> (Vec<f64>, Vec<Self::Point>) {
        monstertruck_traits::algo::curve::parameter_division(self, range, tol)
    }
}

impl<C, S0, S1, R> SearchParameter<CurveParameter> for RollingBallFilletContactCurve<C, S0, S1, R>
where
    C: ParametricCurve3D + SearchNearestParameter<CurveParameter, Point = Point3>,
    S0: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    S1: ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>,
    R: RadiusFunction,
{
    type Point = Point3;
    fn search_parameter<H: Into<SearchParameterHint1D>>(
        &self,
        point: Self::Point,
        hint: H,
        trials: usize,
    ) -> Option<f64> {
        let edge_curve = &self.surface.edge_curve;
        let t = edge_curve.search_nearest_parameter(point, hint, trials)?;
        let cc = self.surface.contact_circle(t)?;
        let q = match self.index {
            0 => cc.contact_point0.point,
            _ => cc.contact_point1.point,
        };
        if point.near(&q) { Some(t) } else { None }
    }
}
