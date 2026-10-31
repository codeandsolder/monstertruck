//! STEP geometry conversions, split by the kind of carrier being
//! converted. The private inherent helpers on
//! [`SurfaceCurveAssociatedGeometry`] live here because two of the parts
//! below need them.

use super::{
    BoundaryCurveFromSamples, BoundedCurve, BsplineCurve, BsplineSurface, Conic2D, Conic3D,
    Curve2D, Curve3D, CurveParameter, CurveTrimRef, Cut, ElementarySurface, Ellipse,
    EuclideanSpace, ExactParameterBoundary2D, ExtrusionSurface, Hyperbola, IncludeCurve,
    InnerSpace, IntersectionCurve, Invertible, KERNEL_EPSILON, KnotVector, Line, Matrix3, Matrix4,
    MetricSpace, NurbsCurve, NurbsSurface, Origin, Parabola, ParameterBoundary2D, ParameterCurve,
    ParameterDivision1D, ParameterRange, ParametricCurve, ParametricCurve3D, ParametricSurface,
    Plane, Point2, Point3, Processor, RevolutionSurface, SearchNearestParameter, SearchParameter,
    SearchParameterHint1D, SearchParameterHint2D, Sphere, SphericalSurface, SquareMatrix,
    StepConvertingError, StepParameterCurve, StepRevolutionSurface, Surface, SurfaceCurve,
    SurfaceCurve3D, SurfaceCurveAssociatedGeometry, SurfaceCurveTrimRef, SweepSurface,
    ToSameGeometry, Tolerance, ToroidalSurface, Torus, Transform, Transformed, TrimmedCurve,
    UnitCircle, UnitHyperbola, UnitParabola, Vector3, Vector4, uv_clamp,
};
use monstertruck_geometry::prelude::{
    HomogeneousSurfaceConversion, SupportsExactPatchDomains, SurfaceParameterRectangle,
    TryIntoBsplineSurface, TryIntoHomogeneousBsplineCurve, TryIntoHomogeneousBsplineSurface,
};
use monstertruck_modeling::{
    Conic2D as ModelingConic2D, Curve as ModelingCurve, Curve2D as ModelingCurve2D,
    Surface as ModelingSurface,
};
use monstertruck_traits::SnapCurveEndpoints;
use std::{cmp::Ordering, env};
// `std::time::Instant::now()` panics on `wasm32-unknown-unknown`. The
// `web_time` crate is std-compatible on native and falls back to
// `performance.now()` in the browser, so all the `Instant::now()` /
// `.elapsed()` call sites below stay unchanged.
use web_time::Instant;

impl SurfaceCurveAssociatedGeometry {
    fn surface(&self) -> &Surface {
        match self {
            Self::ParameterCurve(curve) => curve.surface().as_ref(),
            Self::Surface(surface) => surface,
        }
    }
}

impl SurfaceCurveAssociatedGeometry {
    fn split_at(&mut self, t: f64) -> Self {
        match self {
            Self::ParameterCurve(curve) => Self::ParameterCurve(curve.cut(t)),
            Self::Surface(surface) => Self::Surface(surface.clone()),
        }
    }
}

/// The B-spline / rational-net flattening of STEP geometry.
mod bspline_nets;
/// Curve behaviour on the STEP curve enums.
mod curve_traits;
/// Parameter-space trims of a STEP curve on a STEP surface.
mod parameter_trims;
/// Lifting concrete geometry back into the STEP enums.
mod same_geometry;
/// Conversion into `monstertruck-modeling` geometry.
mod to_modeling;

pub use to_modeling::ROUTE_ANALYTIC_SPHERE;
