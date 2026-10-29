use cgmath::{
    BaseFloat, BaseNum, EuclideanSpace, InnerSpace, Point1, Point2, Point3, Vector1, Vector2,
    Vector3, Vector4, VectorSpace, Zero,
};

/// Shape/capacity error for derivative-array conversion helpers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DerivativeBufferError {
    /// The outer output buffer cannot hold every requested derivative row/item.
    OutputTooShort { required: usize, actual: usize },
    /// An input surface-derivative row has a different width from the first row.
    RaggedInput {
        row: usize,
        expected: usize,
        actual: usize,
    },
    /// One output surface-derivative row cannot hold the full input row width.
    OutputRowTooShort {
        row: usize,
        required: usize,
        actual: usize,
    },
}

impl std::fmt::Display for DerivativeBufferError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Self::OutputTooShort { required, actual } => {
                write!(
                    f,
                    "derivative output has length {actual}, requires at least {required}"
                )
            }
            Self::RaggedInput {
                row,
                expected,
                actual,
            } => write!(
                f,
                "derivative input row {row} has width {actual}, expected {expected}"
            ),
            Self::OutputRowTooShort {
                row,
                required,
                actual,
            } => write!(
                f,
                "derivative output row {row} has width {actual}, requires at least {required}"
            ),
        }
    }
}

impl std::error::Error for DerivativeBufferError {}

pub(crate) fn scalar_from_u128<S: BaseNum>(mut value: u128) -> S {
    let mut result = S::zero();
    let mut bit = S::one();
    while value != 0 {
        if value & 1 == 1 {
            result += bit;
        }
        value >>= 1;
        if value != 0 {
            bit += bit;
        }
    }
    result
}

/// declare control point
pub mod control_point {
    use super::{
        BaseFloat, EuclideanSpace, Point1, Point2, Point3, Vector1, Vector2, Vector3, Vector4, Zero,
    };
    use std::fmt::Debug;
    use std::ops::{
        Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign,
    };
    /// trait for abstract control points of polylines and B-splines
    pub trait ControlPoint<S>:
        Add<Self::Diff, Output = Self>
        + Sub<Self::Diff, Output = Self>
        + Sub<Self, Output = Self::Diff>
        + Mul<S, Output = Self>
        + Div<S, Output = Self>
        + AddAssign<Self::Diff>
        + SubAssign<Self::Diff>
        + MulAssign<S>
        + DivAssign<S>
        + Copy
        + Clone
        + Debug
        + Index<usize, Output = S>
        + IndexMut<usize, Output = S> {
        /// differential vector
        type Diff: Add<Self::Diff, Output = Self::Diff>
            + Sub<Self::Diff, Output = Self::Diff>
            + Mul<S, Output = Self::Diff>
            + Div<S, Output = Self::Diff>
            + AddAssign<Self::Diff>
            + SubAssign<Self::Diff>
            + MulAssign<S>
            + DivAssign<S>
            + Zero
            + Copy
            + Clone
            + Debug
            + Index<usize, Output = S>
            + IndexMut<usize, Output = S>;
        /// dimension
        const DIM: usize;
        /// origin
        fn origin() -> Self;
        /// into the vector
        fn to_vec(self) -> Self::Diff;
        /// from the vector
        fn from_vec(vec: Self::Diff) -> Self;
    }

    impl<S: BaseFloat> ControlPoint<S> for Point1<S> {
        type Diff = Vector1<S>;
        const DIM: usize = 1;
        fn origin() -> Self { EuclideanSpace::origin() }
        fn to_vec(self) -> Self::Diff { EuclideanSpace::to_vec(self) }
        fn from_vec(vec: Self::Diff) -> Self { EuclideanSpace::from_vec(vec) }
    }
    impl<S: BaseFloat> ControlPoint<S> for Point2<S> {
        type Diff = Vector2<S>;
        const DIM: usize = 2;
        fn origin() -> Self { EuclideanSpace::origin() }
        fn to_vec(self) -> Self::Diff { EuclideanSpace::to_vec(self) }
        fn from_vec(vec: Self::Diff) -> Self { EuclideanSpace::from_vec(vec) }
    }
    impl<S: BaseFloat> ControlPoint<S> for Point3<S> {
        type Diff = Vector3<S>;
        const DIM: usize = 3;
        fn origin() -> Self { EuclideanSpace::origin() }
        fn to_vec(self) -> Self::Diff { EuclideanSpace::to_vec(self) }
        fn from_vec(vec: Self::Diff) -> Self { EuclideanSpace::from_vec(vec) }
    }
    impl<S: BaseFloat> ControlPoint<S> for Vector1<S> {
        type Diff = Self;
        const DIM: usize = 1;
        fn origin() -> Self { Zero::zero() }
        fn to_vec(self) -> Self { self }
        fn from_vec(vec: Self::Diff) -> Self { vec }
    }
    impl<S: BaseFloat> ControlPoint<S> for Vector2<S> {
        type Diff = Self;
        const DIM: usize = 2;
        fn origin() -> Self { Zero::zero() }
        fn to_vec(self) -> Self { self }
        fn from_vec(vec: Self::Diff) -> Self { vec }
    }
    impl<S: BaseFloat> ControlPoint<S> for Vector3<S> {
        type Diff = Self;
        const DIM: usize = 3;
        fn origin() -> Self { Zero::zero() }
        fn to_vec(self) -> Self { self }
        fn from_vec(vec: Self::Diff) -> Self { vec }
    }
    impl<S: BaseFloat> ControlPoint<S> for Vector4<S> {
        type Diff = Self;
        const DIM: usize = 4;
        fn origin() -> Self { Zero::zero() }
        fn to_vec(self) -> Self { self }
        fn from_vec(vec: Self::Diff) -> Self { vec }
    }
}

/// Homogeneous coordinate of an Euclidean space and a vector space.
/// # Examples
/// ```
/// use monstertruck_core::cgmath64::*;
/// assert_eq!(Vector4::new(8.0, 6.0, 4.0, 2.0).truncate(), Vector3::new(8.0, 6.0, 4.0));
/// assert_eq!(Vector4::new(8.0, 6.0, 4.0, 2.0).weight(), 2.0);
/// assert_eq!(Vector4::new(8.0, 6.0, 4.0, 2.0).to_point(), Point3::new(4.0, 3.0, 2.0));
/// assert_eq!(Vector4::from_point(Point3::new(4.0, 3.0, 2.0)), Vector4::new(4.0, 3.0, 2.0, 1.0));
/// ```
pub trait Homogeneous: VectorSpace {
    /// The point expressed by homogeneous coordinate
    type Point: EuclideanSpace<Scalar = Self::Scalar>;
    /// Returns the first dim - 1 components.
    fn truncate(self) -> <Self::Point as EuclideanSpace>::Diff;
    /// Returns the last component.
    fn weight(self) -> Self::Scalar;
    /// Returns homogeneous coordinate.
    fn from_point(point: Self::Point) -> Self;
    /// Returns homogeneous coordinate from point and weight.
    #[inline]
    fn from_point_weight(point: Self::Point, weight: Self::Scalar) -> Self {
        Self::from_point(point) * weight
    }
    /// Returns the projection to the plane whose the last component is `1.0`.
    #[inline]
    fn to_point(self) -> Self::Point { Self::Point::from_vec(self.truncate() / self.weight()) }
}

/// Returns the higher order derivation of the rational curve.
/// # Examples
/// ```
/// use monstertruck_core::cgmath64::*;
/// // calculate the derivation at t = 1.5
/// let t = 1.5;
///
/// // the curve: c(t) = (t^2, t^3, t^4, t)
/// let pt = Vector4::new(t * t, t * t * t, t * t * t * t, t);
/// // the derivation: c'(t) = (2t, 3t^2, 4t^3, 1)
/// let der = Vector4::new(2.0 * t, 3.0 * t * t, 4.0 * t * t * t, 1.0);
/// // the 2nd ord. deri.: c''(t) = (2, 6t, 12t^2, 0)
/// let der2 = Vector4::new(2.0, 6.0 * t, 12.0 * t * t, 0.0);
/// // the 3rd ord. deri.: c''(t) = (0, 6, 24t, 0)
/// let der3 = Vector4::new(0.0, 6.0, 24.0 * t, 0.0);
///
/// // the projected curve: \bar{c}(t) = (t, t^2, t^3)
/// // the 1st ord. deri. of the proj'ed curve: \bar{c}''(t) = (1, 2t, 3t^2)
/// let ans_der1 = Vector3::new(1.0, 2.0 * t, 3.0 * t * t);
/// assert_eq!(rational_derivative(&[pt, der]), ans_der1);
/// // the 2nd ord. deri. of the proj'ed curve: \bar{c}''(t) = (0, 2, 6t)
/// let ans_der2 = Vector3::new(0.0, 2.0, 6.0 * t);
/// assert_eq!(rational_derivative(&[pt, der, der2]), ans_der2);
/// // the 3rd ord. deri. of the proj'ed curve: \bar{c}''(t) = (0, 0, 6)
/// let ans_der3 = Vector3::new(0.0, 0.0, 6.0);
/// assert_eq!(rational_derivative(&[pt, der, der2, der3]), ans_der3);
/// ```
pub fn rational_derivative<V: Homogeneous>(ders: &[V]) -> <V::Point as EuclideanSpace>::Diff {
    let zero = <V::Point as EuclideanSpace>::Diff::zero();
    let len = ders.len();
    if len == 0 {
        zero
    } else if len == 1 {
        ders[0].to_point().to_vec()
    } else if len == 2 {
        let (s, sw) = (ders[0].truncate(), ders[0].weight());
        let (d, dw) = (ders[1].truncate(), ders[1].weight());
        let sw2 = sw * sw;
        (d * sw - s * dw) / sw2
    } else if len == 3 {
        let (s, sw) = (ders[0].truncate(), ders[0].weight());
        let (d, dw) = (ders[1].truncate(), ders[1].weight());
        let (d2, d2w) = (ders[2].truncate(), ders[2].weight());
        // SAFETY: The integer 2 is representable in any numeric type.
        let two = scalar_from_u128(2);
        let sw2 = sw * sw;
        d2 / sw - d * (dw / sw2 * two) + s * (dw * dw * two / (sw2 * sw) - d2w / sw2)
    } else if len < 32 {
        let mut evals = [zero; 32];
        rational_derivatives_into(ders, &mut evals);
        evals[ders.len() - 1]
    } else {
        let mut evals = vec![zero; ders.len()];
        rational_derivatives_into(ders, &mut evals);
        evals[ders.len() - 1]
    }
}

/// Store the multi-orders derivations of the rational curve in `evals`.
/// # Examples
/// ```
/// use monstertruck_core::cgmath64::*;
/// // calculate the derivation at t = 1.5
/// let t = 1.5;
///
/// // the curve: c(t) = (t^2, t^3, t^4, t)
/// let ders = [
///     Vector4::new(t * t, t * t * t, t * t * t * t, t), // 0th-order
///     Vector4::new(2.0 * t, 3.0 * t * t, 4.0 * t * t * t, 1.0), // 1st-order
///     Vector4::new(2.0, 6.0 * t, 12.0 * t * t, 0.0), // 2nd-order
///     Vector4::new(0.0, 6.0, 24.0 * t, 0.0), // 3rd-order
/// ];
/// let mut evals = [Vector3::zero(); 4];
/// assert!(rational_derivatives(&ders, &mut evals).is_ok());
///
/// // the projected curve: \bar{c}(t) = (t, t^2, t^3)
/// let ans = [
///     Vector3::new(t, t * t, t * t * t), // 0th-order
///     Vector3::new(1.0, 2.0 * t, 3.0 * t * t), // 1st-order
///     Vector3::new(0.0, 2.0, 6.0 * t), // 2nd-order
///     Vector3::new(0.0, 0.0, 6.0), // 3rd-order
/// ];
/// assert_eq!(evals, ans);
/// ```
fn rational_derivatives_into<V: Homogeneous>(
    ders: &[V],
    evals: &mut [<V::Point as EuclideanSpace>::Diff],
) {
    for i in 0..ders.len() {
        let mut c = 1;
        let sum = (1..i).fold(evals[0] * ders[i].weight(), |sum, j| {
            c = c * (i - j + 1) / j;
            sum + evals[j] * (ders[i - j].weight() * scalar_from_u128(c as u128))
        });
        evals[i] = (ders[i].truncate() - sum) / ders[0].weight();
    }
}

/// Stores the projected derivatives of a rational curve in `evals`.
///
/// # Errors
/// Returns [`DerivativeBufferError::OutputTooShort`] if `evals` cannot hold every input order.
pub fn rational_derivatives<V: Homogeneous>(
    ders: &[V],
    evals: &mut [<V::Point as EuclideanSpace>::Diff],
) -> Result<(), DerivativeBufferError> {
    if evals.len() < ders.len() {
        return Err(DerivativeBufferError::OutputTooShort {
            required: ders.len(),
            actual: evals.len(),
        });
    }
    rational_derivatives_into(ders, evals);
    Ok(())
}

/// Returns the multi-orders derivations of the rational surface.
/// # Examples
/// ```
/// use monstertruck_core::cgmath64::*;
/// // calculate the derivation at (u, v) = (1.0, 2.0).
/// let (u, v) = (1.0, 2.0);
/// // the curve: s(u, v) = (u^3 v^2, u^2 v^3, u v, u)
/// let ders: [[Vector4; 3]; 3] = [
///     [
///         // u-rank = 0, v-rank = 0
///         (u * u * u * v * v, u * u * v * v * v, u * v, u).into(),
///         // u-rank = 0, v-rank = 1
///         (2.0 * u * u * u * v, 3.0 * u * u * v * v, u, 0.0).into(),
///         // u-rank = 0, v-rank = 2
///         (2.0 * u * u * u, 6.0 * u * u * v, 0.0, 0.0).into(),
///     ],
///     [
///         // u-rank = 1, v-rank = 0
///         (3.0 * u * u * v * v, 2.0 * u * v * v * v, v, 1.0).into(),
///         // u-rank = 1, v-rank = 1
///         (6.0 * u * u * v, 6.0 * u * v * v, 1.0, 0.0).into(),
///         // u-rank = 1, v-rank = 2
///         (6.0 * u * u, 12.0 * u * v, 0.0, 0.0).into(),
///     ],
///     [
///         // u-rank = 2, v-rank = 0
///         (6.0 * u * v * v, 2.0 * v * v * v, 0.0, 0.0).into(),
///         // u-rank = 2, v-rank = 1
///         (12.0 * u * v, 6.0 * v * v, 0.0, 0.0).into(),
///         // u-rank = 2, v-rank = 2
///         (12.0 * u, 12.0 * v, 0.0, 0.0).into(),
///     ],
/// ];
///
/// // the projected surface: \bar{s}(u, v) = (u^2 v^2, u v^3, v)
/// assert_eq!(multi_rational_derivative(&ders), Ok(Vector3::new(4.0, 0.0, 0.0)));
/// ```
///
/// # Errors
/// Returns [`DerivativeBufferError::RaggedInput`] when the input rows do not all have the same width.
pub fn multi_rational_derivative<V, A>(
    ders: &[A],
) -> Result<<V::Point as EuclideanSpace>::Diff, DerivativeBufferError>
where
    V: Homogeneous,
    A: AsRef<[V]>, {
    let zero = <V::Point as EuclideanSpace>::Diff::zero();
    if ders.is_empty() {
        return Ok(zero);
    }
    let m = ders.len();
    let n = ders[0].as_ref().len();
    for (row, values) in ders.iter().enumerate().skip(1) {
        let actual = values.as_ref().len();
        if actual != n {
            return Err(DerivativeBufferError::RaggedInput {
                row,
                expected: n,
                actual,
            });
        }
    }
    let value = if n == 0 {
        zero
    } else if (m, n) == (1, 1) {
        ders[0].as_ref()[0].to_point().to_vec()
    } else if m == 1 {
        rational_derivative(ders[0].as_ref())
    } else if (m, n) == (2, 1) {
        rational_derivative(&[ders[0].as_ref()[0], ders[1].as_ref()[0]])
    } else if n == 1 && m < 32 {
        let mut vders = [V::zero(); 32];
        for (vder, array) in vders.iter_mut().zip(ders) {
            *vder = array.as_ref()[0];
        }
        rational_derivative(&vders[..m])
    } else if n == 1 {
        let mut vders = vec![V::zero(); m];
        for (vder, array) in vders.iter_mut().zip(ders) {
            *vder = array.as_ref()[0];
        }
        rational_derivative(&vders)
    } else if (m, n) == (2, 2) {
        let two = scalar_from_u128(2);
        let (value_row, u_row) = (ders[0].as_ref(), ders[1].as_ref());
        let surface_value = value_row[0];
        let u_derivative = u_row[0];
        let v_derivative = value_row[1];
        let cross_derivative = u_row[1];
        let (surface_point, surface_weight) = (surface_value.truncate(), surface_value.weight());
        let (u_vector, u_weight) = (u_derivative.truncate(), u_derivative.weight());
        let (v_vector, v_weight) = (v_derivative.truncate(), v_derivative.weight());
        let (cross_vector, cross_weight) = (cross_derivative.truncate(), cross_derivative.weight());
        let weight_squared = surface_weight * surface_weight;
        cross_vector / surface_weight
            - u_vector * (v_weight / weight_squared)
            - v_vector * (u_weight / weight_squared)
            + surface_point
                * (u_weight * v_weight * two / (weight_squared * surface_weight)
                    - cross_weight / weight_squared)
    } else if m < 8 && n < 8 {
        let mut evals = [[zero; 8]; 8];
        multi_rational_derivatives_into(ders, &mut evals, m, n);
        evals[m - 1][n - 1]
    } else {
        let mut evals = vec![vec![zero; n]; m];
        multi_rational_derivatives_into(ders, &mut evals, m, n);
        evals[m - 1][n - 1]
    };
    Ok(value)
}

/// Store the multi-orders derivations of the rational surface in `evals`.
/// # Remarks
/// - `evals` must be initialized by zero.
/// - Jagged or short inputs are evaluated over their common rectangular prefix.
/// - Extra rows or columns in either input are ignored.
/// # Examples
/// ```
/// use monstertruck_core::cgmath64::*;
/// // calculate the derivation at (u, v) = (1.0, 2.0).
/// let (u, v) = (1.0, 2.0);
/// // the curve: s(u, v) = (u^3 v^2, u^2 v^3, u v, u)
/// let ders: [[Vector4; 3]; 3] = [
///     [
///         // u-rank = 0, v-rank = 0
///         (u * u * u * v * v, u * u * v * v * v, u * v, u).into(),
///         // u-rank = 0, v-rank = 1
///         (2.0 * u * u * u * v, 3.0 * u * u * v * v, u, 0.0).into(),
///         // u-rank = 0, v-rank = 2
///         (2.0 * u * u * u, 6.0 * u * u * v, 0.0, 0.0).into(),
///     ],
///     [
///         // u-rank = 1, v-rank = 0
///         (3.0 * u * u * v * v, 2.0 * u * v * v * v, v, 1.0).into(),
///         // u-rank = 1, v-rank = 1
///         (6.0 * u * u * v, 6.0 * u * v * v, 1.0, 0.0).into(),
///         // u-rank = 1, v-rank = 2
///         (6.0 * u * u, 12.0 * u * v, 0.0, 0.0).into(),
///     ],
///     [
///         // u-rank = 2, v-rank = 0
///         (6.0 * u * v * v, 2.0 * v * v * v, 0.0, 0.0).into(),
///         // u-rank = 2, v-rank = 1
///         (12.0 * u * v, 6.0 * v * v, 0.0, 0.0).into(),
///         // u-rank = 2, v-rank = 2
///         (12.0 * u, 12.0 * v, 0.0, 0.0).into(),
///     ],
/// ];
/// // evals must be initialized by zero.
/// let mut evals = [[Vector3::zero(); 3]; 3];
/// multi_assert!(rational_derivatives(&ders, &mut evals).is_ok());
///
/// // the projected surface: \bar{s}(u, v) = (u^2 v^2, u v^3, v)
/// let ans: [[Vector3; 3]; 3] = [
///     [
///         // u-rank = 0, v-rank = 0
///         (u * u * v * v, u * v * v * v, v).into(),
///         // u-rank = 0, v-rank = 1
///         (2.0 * u * u * v, 3.0 * u * v * v, 1.0).into(),
///         // u-rank = 0, v-rank = 2
///         (2.0 * u * u, 6.0 * u * v, 0.0).into(),
///     ],
///     [
///         // u-rank = 1, v-rank = 0
///         (2.0 * u * v * v, v * v * v, 0.0).into(),
///         // u-rank = 1, v-rank = 1
///         (4.0 * u * v, 3.0 * v * v, 0.0).into(),
///         // u-rank = 1, v-rank = 2
///         (4.0 * u, 6.0 * v, 0.0).into(),
///     ],
///     [
///         // u-rank = 2, v-rank = 0
///         (2.0 * v * v, 0.0, 0.0).into(),
///         // u-rank = 2, v-rank = 1
///         (4.0 * v, 0.0, 0.0).into(),
///         // u-rank = 2, v-rank = 2
///         (4.0, 0.0, 0.0).into(),
///     ],
/// ];
/// assert_eq!(evals, ans);
/// ```
fn multi_rational_derivatives_into<V, A0, A1>(
    ders: &[A0],
    evals: &mut [A1],
    m_max: usize,
    n_max: usize,
) where
    V: Homogeneous,
    A0: AsRef<[V]>,
    A1: AsMut<[<V::Point as EuclideanSpace>::Diff]>,
{
    for m in 0..m_max {
        for n in 0..n_max {
            let mut sum = <V::Point as EuclideanSpace>::Diff::zero();
            let mut c0 = 1;
            for i in 0..=m {
                let mut c1 = 1;
                let (evals, ders) = (evals[i].as_mut(), ders[m - i].as_ref());
                for j in 0..=n {
                    let (c0_s, c1_s) = (scalar_from_u128(c0 as u128), scalar_from_u128(c1 as u128));
                    sum = sum + evals[j] * (ders[n - j].weight() * c0_s * c1_s);
                    c1 = c1 * (n - j) / (j + 1);
                }
                c0 = c0 * (m - i) / (i + 1);
            }
            let (eval_mn, der_mn) = (&mut evals[m].as_mut()[n], ders[m].as_ref()[n]);
            *eval_mn = (der_mn.truncate() - sum) / ders[0].as_ref()[0].weight();
        }
    }
}

/// Stores projected derivatives of a rational surface in `evals`.
///
/// # Errors
/// Returns [`DerivativeBufferError`] if `ders` is ragged or `evals` cannot hold the full
/// rectangular derivative grid.
pub fn multi_rational_derivatives<V, A0, A1>(
    ders: &[A0],
    evals: &mut [A1],
) -> Result<(), DerivativeBufferError>
where
    V: Homogeneous,
    A0: AsRef<[V]>,
    A1: AsMut<[<V::Point as EuclideanSpace>::Diff]>,
{
    if evals.len() < ders.len() {
        return Err(DerivativeBufferError::OutputTooShort {
            required: ders.len(),
            actual: evals.len(),
        });
    }
    if ders.is_empty() {
        return Ok(());
    }
    let n_max = ders[0].as_ref().len();
    for (row, values) in ders.iter().enumerate().skip(1) {
        let actual = values.as_ref().len();
        if actual != n_max {
            return Err(DerivativeBufferError::RaggedInput {
                row,
                expected: n_max,
                actual,
            });
        }
    }
    for (row, values) in evals.iter_mut().take(ders.len()).enumerate() {
        let actual = values.as_mut().len();
        if actual < n_max {
            return Err(DerivativeBufferError::OutputRowTooShort {
                row,
                required: n_max,
                actual,
            });
        }
    }
    multi_rational_derivatives_into(ders, evals, ders.len(), n_max);
    Ok(())
}

impl<S: BaseFloat> Homogeneous for Vector2<S> {
    type Point = Point1<S>;
    #[inline]
    fn truncate(self) -> Vector1<S> { Vector1::new(self[0]) }
    #[inline]
    fn weight(self) -> S { self[1] }
    #[inline]
    fn from_point(point: Self::Point) -> Self { Self::new(point[0], S::one()) }
}

impl<S: BaseFloat> Homogeneous for Vector3<S> {
    type Point = Point2<S>;
    #[inline]
    fn truncate(self) -> Vector2<S> { self.truncate() }
    #[inline]
    fn weight(self) -> S { self[2] }
    #[inline]
    fn from_point(point: Self::Point) -> Self { Self::new(point[0], point[1], S::one()) }
}

impl<S: BaseFloat> Homogeneous for Vector4<S> {
    type Point = Point3<S>;
    #[inline]
    fn truncate(self) -> Vector3<S> { self.truncate() }
    #[inline]
    fn weight(self) -> S { self[3] }
    #[inline]
    fn from_point(point: Self::Point) -> Self { point.to_homogeneous() }
}

/// Store the multi-orders derivations of the magnitude of the curve of vector.
/// # Examples
/// ```
/// use monstertruck_core::{assert_near, cgmath64::*};
///
/// let t = 0.3;
///
/// // c(t) = (2t, 1 - t^2) -> |c(t)| = 1 + t^2
/// let ders = [
///     Vector2::new(2.0 * t, 1.0 - t * t),
///     Vector2::new(2.0, -2.0 * t),
///     Vector2::new(0.0, -2.0),
///     Vector2::new(0.0, 0.0),
/// ];
/// let mut evals = [0.0; 4];
///
/// assert!(absolute_derivatives(&ders, &mut evals).is_ok());
///
/// assert_near!(evals[0], 1.0 + t * t);
/// assert_near!(evals[1], 2.0 * t);
/// assert_near!(evals[2], 2.0);
/// assert_near!(evals[3], 0.0);
/// ```
///
/// # Errors
/// Returns [`DerivativeBufferError::OutputTooShort`] when `evals` cannot hold every input derivative order.
pub fn absolute_derivatives<V>(ders: &[V], evals: &mut [f64]) -> Result<(), DerivativeBufferError>
where V: InnerSpace<Scalar = f64> {
    if evals.len() < ders.len() {
        return Err(DerivativeBufferError::OutputTooShort {
            required: ders.len(),
            actual: evals.len(),
        });
    }
    let n = ders.len();
    for output in evals.iter_mut() {
        *output = 0.0;
    }

    if n == 0 {
        return Ok(());
    }
    evals[0] = ders[0].magnitude();
    for m in 1..n {
        let mut c = 1;
        let sum = (0..m).fold(0.0, |mut sum, i| {
            let x = ders[i + 1].dot(ders[m - 1 - i]);
            let y = evals[i + 1] * evals[m - 1 - i];
            let coefficient = scalar_from_u128(c as u128);
            sum = (x - y).mul_add(coefficient, sum);
            c = c * (m - 1 - i) / (i + 1);
            sum
        });
        evals[m] = sum / evals[0];
    }
    Ok(())
}

#[test]
fn multi_rational_derivatives_reject_short_jagged_buffers() {
    let ders = [
        vec![Vector3::new(2.0, 4.0, 2.0), Vector3::new(1.0, 0.0, 0.0)],
        vec![Vector3::new(0.0, 2.0, 0.0)],
    ];
    let mut evals = [
        vec![Vector2::zero()],
        vec![Vector2::zero(), Vector2::zero()],
    ];
    assert!(matches!(
        multi_rational_derivatives(&ders, &mut evals),
        Err(DerivativeBufferError::RaggedInput { row: 1, .. })
    ));
}
