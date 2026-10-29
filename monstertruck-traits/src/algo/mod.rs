#![allow(
    clippy::many_single_char_names,
    reason = "numerical geometry algorithms use paired u/v notation where similar names are semantically canonical"
)]
#![allow(
    clippy::similar_names,
    reason = "numerical geometry algorithms use paired u/v notation where similar names are semantically canonical"
)]

use crate::traits::{
    ParameterRange, ParametricCurve, ParametricCurve3D, ParametricSurface, ParametricSurface3D,
    SurfaceParameterRange,
};
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

use monstertruck_core::{
    cgmath64::{
        EuclideanSpace, InnerSpace, Matrix2, Matrix3, MetricSpace, Point2, Point3, Vector2,
        Vector3, VectorSpace, Zero,
    },
    hash::HashGen,
    newton::{self, CalcOutput},
    tolerance::{TOLERANCE, Tolerance},
};

/// curve algorithms
pub mod curve;
/// surface algorithms
pub mod surface;
