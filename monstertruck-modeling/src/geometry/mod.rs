use super::{
    BoundaryCurveFromSamples, BoundedCurve, BoundedSurface, BoundingBox, Cut, EuclideanSpace,
    ExactParameterBoundary2D, Homogeneous, IncludeCurve, InnerSpace, Invertible, KERNEL_EPSILON,
    Matrix3, Matrix4, MetricSpace, ParameterBoundary2D, ParameterDivision1D, ParameterDivision2D,
    ParametricCurve, ParametricCurve3D, ParametricSurface, ParametricSurface3D, Point2, Point3,
    Rad, SearchNearestParameter, SearchNearestParameterD1, SearchParameter, SearchParameterD1,
    SearchParameterD2, SearchParameterHint2D, Shell, Solid, SurfaceParameter, ToSameGeometry,
    Tolerance, Transform, Transformed, Vector3, Vector4, bounding,
};
use derive_more::From;
use monstertruck_core::{ContentHasher, DeterministicContentHash};
use monstertruck_geometry::prelude::{
    AnalyticSurfaceKind, BoundaryCurve2D, HomogeneousSurfaceConversion, SupportsExactPatchDomains,
    SurfaceParameterRectangle, TryIntoAnalyticSurfaceKind, TryIntoBsplineSurface,
    TryIntoHomogeneousBsplineCurve, TryIntoHomogeneousBsplineSurface,
};
#[doc(hidden)]
pub use monstertruck_geometry::prelude::{algo, inv_or_zero};
pub use monstertruck_geometry::{decorators::*, nurbs::*, specifieds::*, t_spline::*};
pub use monstertruck_mesh::PolylineCurve;
use monstertruck_topology::compress::{CompressedTrimmedShell, CompressedTrimmedSolid};
use monstertruck_topology::trimmed::{TrimmedShell, TrimmedSolid};
// Only the rayon-parallel (native) `to_trimmed_with_parameter_curves` builds
// faces directly; the wasm32 arm goes through `to_trimmed_with_face_trims`.
#[cfg(not(target_arch = "wasm32"))]
use monstertruck_topology::trimmed::TrimmedFace;
use monstertruck_traits::SnapCurveEndpoints;
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;
use rustc_hash::FxHashMap as HashMap;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::hash::Hasher;
use std::{env, iter};
// `web_time::Instant` is `std::time::Instant` on native and falls back to the
// browser performance clock on wasm32, where `std::time::Instant::now()` panics.
use web_time::Instant;

mod content_hash;
mod curve;
mod domain_projection;
mod parameter_curves;
mod surface;

pub use curve::{Conic2D, Curve, Curve2D};
pub use parameter_curves::{ToCompressedTrimmedParameterCurves, ToTrimmedParameterCurves};
pub use surface::Surface;
