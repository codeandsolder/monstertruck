//! Native nalgebra types and the algebraic contracts needed by generic CAD algorithms.

pub use crate::derivatives::*;
pub use crate::differential::*;
pub use crate::numerical::{affine_qr3, affine_qr4, monic_cubic_roots, monic_quartic_roots};
pub use approx::{AbsDiffEq, RelativeEq, UlpsEq};
pub use nalgebra;
use nalgebra::{ComplexField, Point, RealField, SMatrix, SVector};
pub use nalgebra::{Rotation2, Rotation3, Unit};
pub use num_traits::{One, Zero};
use std::{
    fmt::Debug,
    iter::Sum,
    ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign},
};

macro_rules! f64_type {
    ($($name:ident),* $(,)?) => { $(
        /// Native double-precision geometry type.
        pub type $name = nalgebra::$name<f64>;
    )* };
}
f64_type!(
    Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4, Matrix2, Matrix3, Matrix4
);

/// Unsigned angle between two 3D directions in radians.
/// Using atan2 retains angular resolution near parallel and antiparallel directions.
#[must_use]
pub fn angle_between(a: Vector3, b: Vector3) -> f64 { a.cross(&b).norm().atan2(a.dot(&b)) }

/// Real scalar supported by the fixed-size geometric kernel.
pub trait KernelScalar: RealField + Copy {
    /// Positive infinity for an empty bounding extent.
    fn infinity() -> Self;
    /// Negative infinity for an empty bounding extent.
    fn neg_infinity() -> Self;
}
impl KernelScalar for f32 {
    fn infinity() -> Self { f32::INFINITY }
    fn neg_infinity() -> Self { f32::NEG_INFINITY }
}
impl KernelScalar for f64 {
    fn infinity() -> Self { f64::INFINITY }
    fn neg_infinity() -> Self { f64::NEG_INFINITY }
}

/// Finite-dimensional linear vector used in differentiation and interpolation.
pub trait LinearVector:
    Copy
    + Debug
    + Zero
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + Mul<Self::Scalar, Output = Self>
    + Div<Self::Scalar, Output = Self>
    + MulAssign<Self::Scalar>
    + DivAssign<Self::Scalar>
    + Sum<Self> {
    /// Coordinate scalar.
    type Scalar: KernelScalar;
}
impl<S: KernelScalar, const D: usize> LinearVector for SVector<S, D> {
    type Scalar = S;
}

/// Metric used by proximity searches and geometric subdivision.
pub trait MetricSpace: Sized {
    /// Distance scalar.
    type Metric: KernelScalar;
    /// Squared Euclidean distance.
    fn distance2(self, other: Self) -> Self::Metric;
    /// Euclidean distance.
    fn distance(self, other: Self) -> Self::Metric { self.distance2(other).sqrt() }
}
impl<S: KernelScalar, const D: usize> MetricSpace for SVector<S, D> {
    type Metric = S;
    fn distance2(self, other: Self) -> S { (self - other).norm_squared() }
}
impl<S: KernelScalar, const D: usize> MetricSpace for Point<S, D> {
    type Metric = S;
    fn distance2(self, other: Self) -> S { (self - other).norm_squared() }
}

/// Vector with a Euclidean inner product, used in generic differential geometry.
pub trait InnerProduct: LinearVector + MetricSpace<Metric = Self::Scalar> {
    /// Dot product.
    fn dot(&self, other: &Self) -> Self::Scalar;
    /// Squared Euclidean norm.
    fn norm_squared(&self) -> Self::Scalar { self.dot(self) }
    /// Euclidean norm.
    fn norm(&self) -> Self::Scalar { self.norm_squared().sqrt() }
    /// Unit vector in the same direction; undefined for zero length.
    fn normalize(&self) -> Self { *self / self.norm() }
}
impl<S: KernelScalar, const D: usize> InnerProduct for SVector<S, D> {
    fn dot(&self, other: &Self) -> S { SVector::dot(self, other) }
}

/// Point of an affine space. Subtracting points yields a displacement, not a point.
pub trait AffinePoint:
    Copy
    + Debug
    + Add<Self::Diff, Output = Self>
    + Sub<Self::Diff, Output = Self>
    + Sub<Self, Output = Self::Diff>
    + Mul<Self::Scalar, Output = Self>
    + Div<Self::Scalar, Output = Self>
    + Index<usize, Output = Self::Scalar>
    + IndexMut<usize> {
    /// Coordinate scalar.
    type Scalar: KernelScalar;
    /// Associated displacement vector.
    type Diff: LinearVector<Scalar = Self::Scalar>;
    /// Origin of the coordinate system.
    fn origin() -> Self;
    /// Construct a point from coordinates relative to the origin.
    fn from_coordinates(vector: Self::Diff) -> Self;
    /// Coordinates relative to the origin.
    fn position_vector(self) -> Self::Diff;
    /// Midpoint along the affine segment.
    fn midpoint(self, other: Self) -> Self {
        self + (other - self) / (Self::Scalar::one() + Self::Scalar::one())
    }
}
impl<S: KernelScalar, const D: usize> AffinePoint for Point<S, D> {
    type Scalar = S;
    type Diff = SVector<S, D>;
    fn origin() -> Self { Point::origin() }
    fn from_coordinates(vector: Self::Diff) -> Self { Point::from(vector) }
    fn position_vector(self) -> Self::Diff { self.coords }
}

/// Componentwise arithmetic required by polynomial and attribute interpolation.
pub trait ComponentWise<Rhs = Self>: Sized {
    /// Product of matching coordinate components.
    fn component_product(self, rhs: Rhs) -> Self;
    /// Sum of matching coordinate components.
    fn component_sum(self, rhs: Rhs) -> Self;
}
impl<S: KernelScalar, const D: usize> ComponentWise for SVector<S, D> {
    fn component_product(self, rhs: Self) -> Self { self.component_mul(&rhs) }
    fn component_sum(self, rhs: Self) -> Self { self + rhs }
}
impl<S: KernelScalar, const D: usize> ComponentWise for Point<S, D> {
    fn component_product(self, rhs: Self) -> Self {
        Point::from(self.coords.component_mul(&rhs.coords))
    }
    fn component_sum(self, rhs: Self) -> Self { Point::from(self.coords + rhs.coords) }
}

/// Square linear mapping contract for generic transformed geometry.
pub trait SquareLinearMap: Copy + One {
    /// Coordinate scalar.
    type Scalar: KernelScalar;
    /// Determinant of the linear mapping.
    fn determinant(&self) -> Self::Scalar;
    /// Transposed mapping.
    fn transpose(&self) -> Self;
}
macro_rules! square_map {
    ($($n:literal),*) => { $(impl<S:KernelScalar> SquareLinearMap for SMatrix<S,$n,$n> {
        type Scalar=S;
        fn determinant(&self)->S { SMatrix::determinant(self) }
        fn transpose(&self)->Self { SMatrix::transpose(self) }
    })* };
}
square_map!(2, 3, 4);

/// Apply a CAD coordinate mapping to points and tangent/displacement vectors.
/// Unlike a point, a tangent must not receive the translational component.
pub trait GeometryTransform<P: AffinePoint>: Copy + One {
    /// Map an affine point.
    fn apply_point(&self, point: P) -> P;
    /// Map a displacement or derivative.
    fn apply_vector(&self, vector: P::Diff) -> P::Diff;
    /// Inverse coordinate mapping, if nonsingular.
    fn inverse_map(&self) -> Option<Self>;
    /// Apply inverse to a displacement, if nonsingular.
    fn apply_inverse_vector(&self, vector: P::Diff) -> Option<P::Diff> {
        Some(self.inverse_map()?.apply_vector(vector))
    }
}
impl<S: KernelScalar> GeometryTransform<Point<S, 2>> for SMatrix<S, 3, 3> {
    fn apply_point(&self, p: Point<S, 2>) -> Point<S, 2> {
        Point::from((self * p.to_homogeneous()).fixed_rows::<2>(0).into_owned())
    }
    fn apply_vector(&self, v: SVector<S, 2>) -> SVector<S, 2> {
        (self * v.push(S::zero())).fixed_rows::<2>(0).into_owned()
    }
    fn inverse_map(&self) -> Option<Self> { self.try_inverse() }
}
impl<S: KernelScalar> GeometryTransform<Point<S, 3>> for SMatrix<S, 3, 3> {
    fn apply_point(&self, p: Point<S, 3>) -> Point<S, 3> { Point::from(self * p.coords) }
    fn apply_vector(&self, v: SVector<S, 3>) -> SVector<S, 3> { self * v }
    fn inverse_map(&self) -> Option<Self> { self.try_inverse() }
}
impl<S: KernelScalar> GeometryTransform<Point<S, 3>> for SMatrix<S, 4, 4> {
    fn apply_point(&self, p: Point<S, 3>) -> Point<S, 3> {
        let h = self * p.to_homogeneous();
        Point::from(h.fixed_rows::<3>(0).into_owned() / h[3])
    }
    fn apply_vector(&self, v: SVector<S, 3>) -> SVector<S, 3> {
        (self * v.push(S::zero())).fixed_rows::<3>(0).into_owned()
    }
    fn inverse_map(&self) -> Option<Self> { self.try_inverse() }
}

#[cfg(test)]
mod native_contract_tests {
    use super::*;

    #[test]
    fn translation_affects_points_but_not_derivatives() {
        let transform = Matrix4::new_translation(&Vector3::new(1.0, -2.0, 3.0));
        let point = Point3::new(3.0, 4.0, 5.0);
        let vector = Vector3::new(3.0, 4.0, 5.0);
        assert_eq!(transform.apply_point(point), Point3::new(4.0, 2.0, 8.0));
        assert_eq!(transform.apply_vector(vector), vector);
        assert_eq!(
            transform
                .inverse_map()
                .map(|m| m.apply_point(transform.apply_point(point))),
            Some(point)
        );
    }

    #[test]
    fn projective_weights_and_native_column_order() {
        let transform = Matrix4::new(
            2.0, 0.0, 0.0, 8.0, 0.0, 3.0, 0.0, 6.0, 0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 0.0, 2.0,
        );
        assert_eq!(
            transform.apply_point(Point3::new(1.0, 2.0, 3.0)),
            Point3::new(5.0, 6.0, 8.0)
        );
        assert_eq!(
            transform.apply_vector(Vector3::new(1.0, 2.0, 3.0)),
            Vector3::new(2.0, 6.0, 12.0)
        );
        assert_eq!(
            transform.column(3).into_owned(),
            Vector4::new(8.0, 6.0, 4.0, 2.0)
        );
    }

    #[test]
    fn native_storage_preserves_generic_single_precision() {
        let transform =
            nalgebra::Matrix3::<f32>::new_translation(&nalgebra::Vector2::new(2.0, -3.0));
        let point = nalgebra::Point2::new(1.0, 1.0);
        let vector = nalgebra::Vector2::new(1.0, 1.0);
        assert_eq!(
            transform.apply_point(point),
            nalgebra::Point2::new(3.0, -2.0)
        );
        assert_eq!(transform.apply_vector(vector), vector);
        assert_eq!(
            std::mem::size_of::<nalgebra::Vector3<f32>>(),
            3 * std::mem::size_of::<f32>()
        );
        assert_eq!(
            std::mem::size_of::<Vector3>(),
            3 * std::mem::size_of::<f64>()
        );
    }
}
