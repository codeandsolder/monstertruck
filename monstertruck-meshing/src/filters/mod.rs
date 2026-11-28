use crate::{
    Adjacency, AffinePoint, BoundingBox, ComponentWise, FaceNormal, Faces, Matrix3, MetricSpace,
    Point3, PolygonMesh, StandardAttributes, StandardVertex, Tolerance, Vector2, Vector3, Vertex,
    Zero, nonpositive_tolerance, polygon_mesh,
};
use polygon_mesh::PolygonMeshEditor;

mod normal_filters;
mod optimizing;
mod structuring;
mod subdivision;

pub use normal_filters::NormalFilters;
pub use optimizing::OptimizingFilter;
pub use structuring::StructuringFilter;
pub use subdivision::Subdivision;
