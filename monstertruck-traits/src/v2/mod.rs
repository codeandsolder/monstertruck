//! Scalar-generic v2 trait family.
//!
//! These traits are the scalar-generic future of the monstertruck geometry
//! kernel. They mirror the existing `f64`-hardcoded traits in [`crate::traits`]
//! but parameterize the scalar type via associated types.
//!
//! The established kernel remains primarily `f64`, while selected geometry
//! types also implement these traits for alternate scalar types. Implementations
//! are native rather than bridged through blanket compatibility adapters.

mod bounded_curve;
mod bounded_surface;
mod curve;
mod curve3d;
mod cut;
mod division;
mod parameter_boundary;
mod search_parameter;
mod surface;
mod surface3d;

/// Scalar-generic algorithm helpers.
pub mod algo;

pub use bounded_curve::BoundedCurve;
pub use bounded_surface::BoundedSurface;
pub use curve::ParametricCurve;
pub use curve3d::ParametricCurve3D;
pub use cut::Cut;
pub use division::{ParameterDivision1D, ParameterDivision2D, SurfaceParameterRange};
pub use parameter_boundary::ParameterBoundary2D;
pub use search_parameter::{
    CurveParameter, ParameterRange, SearchNearestParameter, SearchParameter,
    SearchParameterDimension, SearchParameterHint1D, SearchParameterHint2D, SurfaceParameter,
};
pub use surface::ParametricSurface;
pub use surface3d::ParametricSurface3D;
