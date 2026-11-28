pub use crate::cgmath_extend_traits::*;
pub use crate::derivatives::*;
pub use crate::numerical::{affine_qr3, affine_qr4, monic_cubic_roots, monic_quartic_roots};
pub use cgmath;
pub use cgmath::prelude::*;
pub use cgmath::{Deg, Rad, frustum, ortho, perspective};
macro_rules! f64_type {
        ($typename: ident) => {
            /// redefinition, scalar = f64
            pub type $typename = cgmath::$typename<f64>;
        };
        ($a: ident, $($b: ident), *) => { f64_type!($a); f64_type!($($b),*); }
    }
f64_type!(
    Vector1, Vector2, Vector3, Vector4, Matrix2, Matrix3, Matrix4, Point1, Point2, Point3
);
