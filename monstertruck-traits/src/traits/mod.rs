use monstertruck_core::cgmath64::{
    CurveDerivatives, InnerSpace, Matrix3, Matrix4, Point2, Point3, SurfaceDerivatives, Transform,
    Vector2, Vector3, Zero,
};
use std::ops::Bound;

mod curve;
pub use curve::*;
mod surface;
pub use surface::*;
mod search_parameter;
pub use search_parameter::*;

/// parameter range
pub type ParameterRange = (Bound<f64>, Bound<f64>);
fn bound2opt<T>(x: Bound<T>) -> Option<T> {
    match x {
        Bound::Included(x) | Bound::Excluded(x) => Some(x),
        Bound::Unbounded => None,
    }
}

/// Oriented and reversible
pub trait Invertible: Clone {
    /// Inverts `self`
    fn invert(&mut self);
    /// Returns the inverse.
    #[inline]
    #[must_use]
    fn inverse(&self) -> Self {
        let mut res = self.clone();
        res.invert();
        res
    }
}

/// Implementation for the test of topological methods.
impl Invertible for (usize, usize) {
    fn invert(&mut self) { *self = (self.1, self.0); }
    fn inverse(&self) -> Self { (self.1, self.0) }
}

impl<P: Clone> Invertible for Vec<P> {
    #[inline]
    fn invert(&mut self) { self.reverse(); }
    #[inline]
    fn inverse(&self) -> Self { self.iter().rev().cloned().collect() }
}

impl<T: Invertible> Invertible for Box<T> {
    #[inline]
    fn invert(&mut self) { (**self).invert(); }
    #[inline]
    fn inverse(&self) -> Self { Self::new((**self).inverse()) }
}

/// Transform geometry
pub trait Transformed<T>: Clone {
    /// transform by `trans`.
    fn transform_by(&mut self, trans: T);
    /// transformed geometry by `trans`.
    #[inline]
    #[must_use]
    fn transformed(&self, trans: T) -> Self {
        let mut res = self.clone();
        res.transform_by(trans);
        res
    }
}

impl<S: Transformed<T>, T> Transformed<T> for Box<S> {
    #[inline]
    fn transform_by(&mut self, trans: T) { (**self).transform_by(trans); }
    #[inline]
    fn transformed(&self, trans: T) -> Self { Self::new((**self).transformed(trans)) }
}

macro_rules! impl_transformed {
    ($point: ty, $matrix: ty) => {
        impl Transformed<$matrix> for $point {
            #[inline]
            fn transform_by(&mut self, trans: $matrix) { *self = trans.transform_point(*self) }
            #[inline]
            fn transformed(&self, trans: $matrix) -> Self { trans.transform_point(*self) }
        }
    };
}
impl_transformed!(Point2, Matrix3);
impl_transformed!(Point3, Matrix3);
impl_transformed!(Point3, Matrix4);

/// Obtain a curve or surface that gives the same image as a given curve or surface.
pub trait ToSameGeometry<T> {
    /// Obtain a curve or surface that gives the same image as a given curve or surface.
    fn to_same_geometry(&self) -> T;
}
