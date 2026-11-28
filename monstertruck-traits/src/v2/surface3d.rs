use monstertruck_core::{Point3, Vector3};

/// Scalar-generic 3D parametric surface.
///
/// Adds normal-vector computation to [`super::ParametricSurface`] for surfaces
/// embedded in 3D Euclidean space.
///
/// # Bridge-phase note
///
/// This trait is **parameter-scalar generic** but fixes the ambient geometry to
/// `Point3` / `Vector3` (i.e. `f64`).  Making the ambient space scalar-generic
/// as well is deferred to a later phase.
pub trait ParametricSurface3D: super::ParametricSurface<Point = Point3, Vector = Vector3> {
    /// Returns the unit normal vector at `(u, v)`.
    fn normal(&self, u: Self::Scalar, v: Self::Scalar) -> Vector3 {
        self.derivative_u(u, v)
            .cross(&self.derivative_v(u, v))
            .normalize()
    }
    /// Returns the derivative of the normal vector w.r.t. `u` at `(u, v)`.
    fn normal_uder(&self, u: Self::Scalar, v: Self::Scalar) -> Vector3 {
        monstertruck_core::surface_normal_derivative(
            [self.derivative_u(u, v), self.derivative_v(u, v)],
            [self.derivative_uu(u, v), self.derivative_uv(u, v)],
        )
    }
    /// Returns the derivative of the normal vector w.r.t. `v` at `(u, v)`.
    fn normal_vder(&self, u: Self::Scalar, v: Self::Scalar) -> Vector3 {
        monstertruck_core::surface_normal_derivative(
            [self.derivative_u(u, v), self.derivative_v(u, v)],
            [self.derivative_uv(u, v), self.derivative_vv(u, v)],
        )
    }
}
