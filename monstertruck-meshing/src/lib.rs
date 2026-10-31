//! Tessellation and meshing algorithms for B-rep shapes.
//!
//! # Examples
//!
//! ```
//! use monstertruck_meshing::prelude::*;
//! use monstertruck_modeling::*;
//!
//! // Build a cube and tessellate it
//! let v = builder::vertex(Point3::origin());
//! let cube: Solid = builder::extrude(
//!     &builder::extrude(&builder::extrude(&v, Vector3::unit_x()), Vector3::unit_y()),
//!     Vector3::unit_z(),
//! );
//!
//! let mesh = cube.triangulation(0.01).to_polygon();
//! ```

#![cfg_attr(not(debug_assertions), deny(warnings))]
#![deny(clippy::all, rust_2018_idioms)]

use common::{Adjacency, FaceNormal};

/// re-export polymesh
pub mod rexport_polymesh {
    pub use monstertruck_mesh::*;
}

#[expect(
    clippy::cast_precision_loss,
    reason = "mesh collection/sampling counts are bounded by actually iterable/allocatable data and are converted only for normalized f64 geometry formulas"
)]
pub(crate) const fn count_as_f64(value: usize) -> f64 { value as f64 }
#[cfg(all(feature = "tessellation", not(feature = "analyzers")))]
use monstertruck_mesh::KERNEL_EPSILON;
#[cfg(any(feature = "analyzers", feature = "tessellation"))]
use monstertruck_mesh::Origin;
#[cfg(feature = "analyzers")]
use monstertruck_mesh::PolylineCurve;
#[cfg(feature = "tessellation")]
use monstertruck_mesh::{
    BoundedCurve, CurveParameter, ExactParameterBoundary2D, HashGen, Invertible,
    ParameterBoundary2D, ParameterDivision1D, ParameterDivision2D, ParametricCurve2D,
    ParametricCurve3D, ParametricSurface3D, Point2, SearchNearestParameter, SearchParameter,
    SurfaceParameter, algo,
};
#[cfg(any(feature = "analyzers", feature = "filters"))]
use monstertruck_mesh::{BoundingBox, Matrix3, MetricSpace, SquareMatrix};
#[cfg(feature = "filters")]
use monstertruck_mesh::{ElementWise, Tolerance, Vector2, polygon_mesh};
use monstertruck_mesh::{
    EuclideanSpace, Faces, InnerSpace, Point3, PolygonMesh, StandardVertex as Vertex, Vector3, Zero,
};
#[cfg(feature = "analyzers")]
use monstertruck_mesh::{KERNEL_EPSILON, KERNEL_EPSILON_SQUARED, Vector4, hash};
#[cfg(any(feature = "analyzers", feature = "filters", feature = "tessellation"))]
use monstertruck_mesh::{StandardAttributes, StandardVertex, nonpositive_tolerance};

/// polygon mesh analyzers, including
///
/// - determines topological properties: connectivity, boundary extraction, or shell conditions (closed or oriented)
/// - detects collisions between two meshes and extracts interference lines
/// - investigates positional relations between mesh and point clouds.
#[cfg(feature = "analyzers")]
pub mod analyzers;
mod common;
/// Edits meshes. Add normals, optimizing data, and so on.
#[cfg(feature = "filters")]
pub mod filters;
/// Tessellates shapes.
#[cfg(feature = "tessellation")]
pub mod tessellation;

/// VTK Output
#[cfg(feature = "vtk")]
#[cfg(not(target_arch = "wasm32"))]
pub mod vtk;

/// This module contains all traits and re-exports `monstertruck_mesh`.
pub mod prelude {
    #[cfg(feature = "analyzers")]
    pub use crate::analyzers::*;
    #[cfg(feature = "filters")]
    pub use crate::filters::*;
    pub use crate::rexport_polymesh::*;
    #[cfg(feature = "tessellation")]
    pub use crate::tessellation::*;
    #[cfg(feature = "vtk")]
    #[cfg(not(target_arch = "wasm32"))]
    pub use crate::vtk::*;
}
