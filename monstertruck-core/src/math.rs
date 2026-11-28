//! Native nalgebra types and the algebraic contracts needed by generic CAD algorithms.

pub use crate::derivatives::*;
pub use crate::differential::*;
pub use crate::numerical::{affine_qr, monic_cubic_roots, monic_quartic_roots};
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
#[must_use = "the returned angle is used for geometric orientation tests"]
pub fn angle_between(a: Vector3, b: Vector3) -> f64 { a.cross(&b).norm().atan2(a.dot(&b)) }

/// Derivative of the unit normal along a specified parameter direction.
/// Both tangents and their directional derivatives must be expressed in the same frame.
#[must_use = "the returned derivative is used in differential geometry"]
pub fn surface_normal_derivative([du, dv]: [Vector3; 2], [ddu, ddv]: [Vector3; 2]) -> Vector3 {
    let cross = du.cross(&dv);
    let derivative = ddu.cross(&dv) + du.cross(&ddv);
    let magnitude = cross.norm();
    let magnitude_derivative = cross.dot(&derivative) / magnitude;
    (derivative * magnitude - cross * magnitude_derivative) / magnitude.powi(2)
}

/// Circumcenter of three noncollinear points, obtained from their 2D Gram system.
/// The result is undefined for a collinear or repeated triple.
#[must_use = "the returned point is the center of the circumcircle"]
pub fn circumcenter(first: Point3, second: Point3, third: Point3) -> Point3 {
    let (a, b) = (second - first, third - first);
    let (aa, ab, bb) = (a.dot(&a), a.dot(&b), b.dot(&b));
    let determinant = ab.mul_add(-ab, aa * bb);
    let u = (aa - ab) * bb;
    let v = (bb - ab) * aa;
    first + u / (2.0 * determinant) * a + v / (2.0 * determinant) * b
}

/// Tangent sphere center and contact points for two oriented planes, constrained
/// to a third plane through the supplied point.
///
/// Surface normals must be unit vectors.
/// A negative signed radius selects the opposite contact side.
#[must_use = "the returned center and contacts must be used to construct the fillet"]
pub fn tangent_sphere_contacts(
    (point, section_normal): (Point3, Vector3),
    [(point0, normal0), (point1, normal1)]: [(Point3, Vector3); 2],
    signed_radius: f64,
) -> Option<(Point3, Point3, Point3)> {
    let matrix = Matrix3::from_columns(&[section_normal, normal0, normal1]).transpose();
    let rhs = Vector3::new(
        section_normal.dot(&point.coords),
        normal0.dot(&point0.coords) - signed_radius,
        normal1.dot(&point1.coords) - signed_radius,
    );
    let center = Point3::from(matrix.try_inverse()? * rhs);
    Some((
        center,
        center + signed_radius * normal0,
        center + signed_radius * normal1,
    ))
}

/// Arbitrary angular derivative of Rodrigues' rotation matrix about a fixed unit axis.
#[must_use = "the returned matrix is the angular derivative of the rotation"]
pub fn axis_angle_derivative(n: usize, axis: Vector3, angle: f64) -> Matrix3 {
    let (s, c) = f64::sin_cos(angle);
    let (s, c) = match n % 4 {
        0 => (s, c),
        1 => (c, -s),
        2 => (-s, -c),
        _ => (-c, s),
    };
    let one_minus_cosine_derivative = match n {
        0 => 1.0 - c,
        _ => -c,
    };

    Matrix3::new(
        (one_minus_cosine_derivative * axis.x).mul_add(axis.x, c),
        f64::mul_add(s, -axis.z, one_minus_cosine_derivative * axis.x * axis.y),
        f64::mul_add(s, axis.y, one_minus_cosine_derivative * axis.x * axis.z),
        f64::mul_add(s, axis.z, one_minus_cosine_derivative * axis.x * axis.y),
        (one_minus_cosine_derivative * axis.y).mul_add(axis.y, c),
        f64::mul_add(s, -axis.x, one_minus_cosine_derivative * axis.y * axis.z),
        f64::mul_add(s, -axis.y, one_minus_cosine_derivative * axis.x * axis.z),
        f64::mul_add(s, axis.x, one_minus_cosine_derivative * axis.y * axis.z),
        (one_minus_cosine_derivative * axis.z).mul_add(axis.z, c),
    )
}

/// Real scalar supported by the fixed-size geometric kernel.
pub trait KernelScalar: RealField + Copy {
    /// Positive infinity for an empty bounding extent.
    fn infinity() -> Self;
    /// Negative infinity for an empty bounding extent.
    fn neg_infinity() -> Self;
}
impl KernelScalar for f32 {
    fn infinity() -> Self { Self::INFINITY }
    fn neg_infinity() -> Self { Self::NEG_INFINITY }
}
impl KernelScalar for f64 {
    fn infinity() -> Self { Self::INFINITY }
    fn neg_infinity() -> Self { Self::NEG_INFINITY }
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
    #[must_use = "the returned vector is normalized without modifying the input"]
    fn normalize(&self) -> Self { *self / self.norm() }
}
impl<S: KernelScalar, const D: usize> InnerProduct for SVector<S, D> {
    fn dot(&self, other: &Self) -> S { Self::dot(self, other) }
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
    #[must_use = "the returned point is the segment midpoint"]
    fn midpoint(self, other: Self) -> Self {
        let one = Self::Scalar::one();
        let two = one + one;
        self + (other - self) / two
    }
}
impl<S: KernelScalar, const D: usize> AffinePoint for Point<S, D> {
    type Scalar = S;
    type Diff = SVector<S, D>;
    fn origin() -> Self { Self::origin() }
    fn from_coordinates(vector: Self::Diff) -> Self { Self::from(vector) }
    fn position_vector(self) -> Self::Diff { self.coords }
}

/// Componentwise arithmetic required by polynomial and attribute interpolation.
pub trait ComponentWise<Rhs = Self>: Sized {
    /// Product of matching coordinate components.
    #[must_use = "the result contains the products of corresponding coordinates"]
    fn component_product(self, rhs: Rhs) -> Self;
    /// Sum of matching coordinate components.
    #[must_use = "the result contains the sums of corresponding coordinates"]
    fn component_sum(self, rhs: Rhs) -> Self;
}
impl<S: KernelScalar, const D: usize> ComponentWise for SVector<S, D> {
    fn component_product(self, rhs: Self) -> Self { self.component_mul(&rhs) }
    fn component_sum(self, rhs: Self) -> Self { self + rhs }
}
impl<S: KernelScalar, const D: usize> ComponentWise for Point<S, D> {
    fn component_product(self, rhs: Self) -> Self {
        Self::from(self.coords.component_mul(&rhs.coords))
    }
    fn component_sum(self, rhs: Self) -> Self { Self::from(self.coords + rhs.coords) }
}

/// Square linear mapping contract for generic transformed geometry.
pub trait SquareLinearMap: Copy + One {
    /// Coordinate scalar.
    type Scalar: KernelScalar;
    /// Determinant of the linear mapping.
    fn determinant(&self) -> Self::Scalar;
    /// Transposed mapping.
    #[must_use = "the original mapping remains unchanged"]
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
        assert_eq!(
            <nalgebra::Matrix3<f32> as GeometryTransform<nalgebra::Point2<f32>>>::apply_vector(
                &transform, vector
            ),
            vector
        );
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
