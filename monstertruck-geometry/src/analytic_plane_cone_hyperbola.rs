//! The analytic hyperbola at an axis-parallel, axis-offset plane/cone intersection.
//!
//! This is a deliberately named conic regime. Oblique planes and planes
//! containing the cone axis are not misclassified as hyperbolas.

use crate::prelude::{
    InnerSpace, Matrix4, Point3, Processor, TrimmedCurve, UnitHyperbola, Vector3, Vector4,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

/// Positive-radius cone/plane intersection, with a unique real branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisParallelPlaneConeHyperbola {
    origin: Point3,
    axis: Vector3,
    plane_normal: Vector3,
    transverse: Vector3,
    offset: f64,
    radius: f64,
    slope: f64,
    /// Source normal/axis orthogonality residual (numerical, not silently zero).
    pub normal_axis_misalignment: f64,
}

/// Directed trim of a newly constructed analytic conic, not the original cubic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedHyperbolaTrim {
    pub start_parameter: f64,
    pub end_parameter: f64,
    pub start_residual: f64,
    pub end_residual: f64,
}

impl OrientedHyperbolaTrim {
    #[must_use]
    pub fn forward(self) -> bool { self.start_parameter < self.end_parameter }
}

fn valid_point(point: Point3) -> bool {
    [point.x, point.y, point.z].into_iter().all(f64::is_finite)
}

fn unit(vector: Vector3) -> bool {
    vector.magnitude2().is_finite()
        && (vector.magnitude2() - 1.0).abs()
            <= DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
}

impl AxisParallelPlaneConeHyperbola {
    /// Derive a single-nappe hyperbola from independently authored carriers.
    ///
    /// The source plane is parallel to the cone axis within a reported
    /// floating-point misalignment, and must not contain that axis. Other
    /// intersection regimes are rejected; there is no silent healing.
    #[must_use]
    pub fn from_carriers(
        plane_origin: Point3,
        plane_normal: Vector3,
        cone_origin: Point3,
        cone_axis: Vector3,
        reference_radius: f64,
        slope: f64,
    ) -> Option<Self> {
        if !valid_point(plane_origin)
            || !valid_point(cone_origin)
            || !unit(plane_normal)
            || !unit(cone_axis)
            || !reference_radius.is_finite()
            || reference_radius <= 0.0
            || !slope.is_finite()
            || slope <= 0.0
        {
            return None;
        }
        let misalignment = plane_normal.dot(cone_axis);
        if !misalignment.is_finite()
            || misalignment.abs() > DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
        {
            return None;
        }
        // Project away small source-axis direction noise; retain its amount
        // so that users can account for non-exact input plane parallelism.
        let axis = cone_axis.normalize();
        let normal = (plane_normal - axis * plane_normal.dot(axis)).normalize();
        let transverse = axis.cross(normal).normalize();
        let offset = (plane_origin - cone_origin).dot(normal);
        let scale = [
            plane_origin.x.abs(),
            plane_origin.y.abs(),
            plane_origin.z.abs(),
            cone_origin.x.abs(),
            cone_origin.y.abs(),
            cone_origin.z.abs(),
            reference_radius,
            1.0,
        ]
        .into_iter()
        .fold(0.0_f64, f64::max);
        let floor = 256.0 * f64::EPSILON * scale;
        // A plane through the cone axis produces intersecting generators,
        // not a hyperbola; a near-axis plane is not numerically resolvable.
        if !offset.is_finite() || offset.abs() <= floor || !floor.is_finite() {
            return None;
        }
        let semireal = offset.abs() / slope;
        let vertex_axial = (offset.abs() - reference_radius) / slope;
        if !semireal.is_finite() || semireal <= 0.0 || !vertex_axial.is_finite() {
            return None;
        }
        Some(Self {
            origin: cone_origin,
            axis,
            plane_normal: normal,
            transverse,
            offset,
            radius: reference_radius,
            slope,
            normal_axis_misalignment: misalignment,
        })
    }

    /// Evaluate the analytic hyperbola, including its smooth vertex at t=0.
    #[must_use]
    pub fn evaluate(self, t: f64) -> Option<Point3> {
        if !t.is_finite() {
            return None;
        }
        let transverse = self.offset.abs() * t.sinh();
        let radius = self.offset.abs() * t.cosh();
        let axial = (radius - self.radius) / self.slope;
        let point = self.origin
            + self.plane_normal * self.offset
            + self.transverse * transverse
            + self.axis * axial;
        valid_point(point).then_some(point)
    }

    /// Invert an original source topological vertex; reject any off-branch
    /// point rather than projecting it onto an unrelated nearby conic.
    #[must_use]
    pub fn parameter_at_vertex(self, point: Point3, tolerance: f64) -> Option<(f64, f64)> {
        if !valid_point(point) || !tolerance.is_finite() || tolerance < 0.0 {
            return None;
        }
        let projected = (point - self.origin).dot(self.transverse);
        let parameter = (projected / self.offset.abs()).asinh();
        let actual = self.evaluate(parameter)?;
        let residual = (actual - point).magnitude();
        (residual.is_finite() && residual <= tolerance).then_some((parameter, residual))
    }

    /// Maintain source edge start/end direction, including decreasing conic
    /// parameters across the vertex; no periodic wrapping is permitted.
    #[must_use]
    pub fn oriented_trim(
        self,
        start: Point3,
        end: Point3,
        tolerance: f64,
    ) -> Option<OrientedHyperbolaTrim> {
        let (start_parameter, start_residual) = self.parameter_at_vertex(start, tolerance)?;
        let (end_parameter, end_residual) = self.parameter_at_vertex(end, tolerance)?;
        if !(start - end).magnitude().is_finite()
            || (start - end).magnitude() <= tolerance
            || start_parameter.total_cmp(&end_parameter).is_eq()
        {
            return None;
        }
        Some(OrientedHyperbolaTrim {
            start_parameter,
            end_parameter,
            start_residual,
            end_residual,
        })
    }

    /// Exact kernel-native hyperbola carrier, exportable as STEP HYPERBOLA.
    ///
    /// The source cubic is not reused and is not claimed equivalent. The
    /// returned trim retains the original directed endpoint order.
    #[must_use]
    pub fn modeling_curve(
        self,
        trim: OrientedHyperbolaTrim,
    ) -> Processor<TrimmedCurve<UnitHyperbola<Point3>>, Matrix4> {
        let real = self.offset.abs() / self.slope;
        let imaginary = self.offset.abs();
        let center =
            self.origin + self.plane_normal * self.offset - self.axis * (self.radius / self.slope);
        let normal = self.axis.cross(self.transverse);
        let transform = Matrix4::from_cols(
            (self.axis * real).extend(0.0),
            (self.transverse * imaginary).extend(0.0),
            (normal * real.min(imaginary)).extend(0.0),
            Vector4::new(center.x, center.y, center.z, 1.0),
        );
        Processor::with_transform(
            TrimmedCurve::new(
                UnitHyperbola::<Point3>::new(),
                (trim.start_parameter, trim.end_parameter),
            ),
            transform,
        )
    }

    /// Upper bound on the maximum distance of the *entire finite trim*
    /// from its original authored source plane, not only the endpoints.
    ///
    /// For the constructed analytic branch the signed plane residual is
    /// `C + A*cosh(t) + B*sinh(t)`. Bounding *both* hyperbolic functions
    /// over the complete closed parameter interval also covers floating-
    /// point transverse-frame nonorthogonality and small source plane tilt.
    /// A finite-f64 arithmetic allowance is added; this is not formally
    /// outward-rounded IEEE interval arithmetic.
    #[must_use]
    pub fn source_plane_residual_upper(
        self,
        original_plane_origin: Point3,
        original_plane_normal: Vector3,
        trim: OrientedHyperbolaTrim,
    ) -> Option<f64> {
        if !valid_point(original_plane_origin)
            || !unit(original_plane_normal)
            || !trim.start_parameter.is_finite()
            || !trim.end_parameter.is_finite()
        {
            return None;
        }
        let largest_parameter = trim.start_parameter.abs().max(trim.end_parameter.abs());
        let largest_cosh = largest_parameter.cosh();
        let largest_sinh = largest_parameter.sinh();
        let real = self.offset.abs() / self.slope;
        let imaginary = self.offset.abs();
        let center =
            self.origin + self.plane_normal * self.offset - self.axis * (self.radius / self.slope);
        let constant = (center - original_plane_origin).dot(original_plane_normal);
        let cosh_coefficient = real * self.axis.dot(original_plane_normal);
        let sinh_coefficient = imaginary * self.transverse.dot(original_plane_normal);
        let algebraic_upper = sinh_coefficient.abs().mul_add(
            largest_sinh,
            cosh_coefficient.abs().mul_add(largest_cosh, constant.abs()),
        );
        let coordinate_scale = [
            original_plane_origin.x.abs(),
            original_plane_origin.y.abs(),
            original_plane_origin.z.abs(),
            self.origin.x.abs(),
            self.origin.y.abs(),
            self.origin.z.abs(),
            center.x.abs(),
            center.y.abs(),
            center.z.abs(),
            real,
            imaginary,
            1.0,
        ]
        .into_iter()
        .fold(0.0_f64, f64::max);
        let margin = 512.0 * f64::EPSILON * coordinate_scale * (1.0 + largest_cosh + largest_sinh);
        let upper = algebraic_upper + margin;
        upper.is_finite().then_some(upper)
    }

    #[must_use]
    pub const fn reference_origin(self) -> Point3 { self.origin }
    #[must_use]
    pub const fn signed_plane_offset(self) -> f64 { self.offset }
    #[must_use]
    pub const fn reference_radius(self) -> f64 { self.radius }
    #[must_use]
    pub const fn slope(self) -> f64 { self.slope }
    #[must_use]
    pub const fn plane_normal(self) -> Vector3 { self.plane_normal }
    #[must_use]
    pub const fn cone_axis(self) -> Vector3 { self.axis }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::ParametricCurve;

    fn example() -> AxisParallelPlaneConeHyperbola {
        AxisParallelPlaneConeHyperbola::from_carriers(
            Point3::new(0.0, -0.3, 2.0),
            Vector3::new(0.0, -1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            0.5,
            1.5,
        )
        .unwrap()
    }

    #[test]
    fn full_hyperbola_lies_on_both_exact_carriers() {
        let curve = example();
        for i in -32..=32 {
            let t = f64::from(i) / 16.0;
            let p = curve.evaluate(t).unwrap();
            assert!((p.x.hypot(p.y) - (0.5 + 1.5 * p.z)).abs() < 2.0e-14);
            assert!((p.y + 0.3).abs() < 1.0e-14);
            let (inverted, residual) = curve.parameter_at_vertex(p, 1.0e-12).unwrap();
            assert!((inverted - t).abs() < 1.0e-12);
            assert!(residual < 1.0e-12);
        }
    }

    #[test]
    fn modeling_curve_preserves_reversed_nonperiodic_parameters() {
        let hyperbola = example();
        let start = hyperbola.evaluate(1.2).unwrap();
        let end = hyperbola.evaluate(-0.4).unwrap();
        let trim = hyperbola.oriented_trim(start, end, 1.0e-12).unwrap();
        assert!(!trim.forward());
        let output = hyperbola.modeling_curve(trim);
        assert!((output.evaluate(trim.start_parameter) - start).magnitude() < 1.0e-12);
        assert!((output.evaluate(trim.end_parameter) - end).magnitude() < 1.0e-12);
        assert!(
            hyperbola
                .parameter_at_vertex(end + Vector3::new(0.01, 0.0, 0.0), 1.0e-8)
                .is_none()
        );
    }

    #[test]
    fn full_trim_source_plane_bound_includes_transverse_projection_residual() {
        let curve = example();
        let source_origin = Point3::new(0.0, -0.3, 2.0);
        let source_normal = Vector3::new(1.0e-9, -1.0, 0.0);
        let trim = curve
            .oriented_trim(
                curve.evaluate(-1.3).unwrap(),
                curve.evaluate(1.8).unwrap(),
                1.0e-12,
            )
            .unwrap();
        let upper = curve
            .source_plane_residual_upper(source_origin, source_normal, trim)
            .unwrap();
        assert!(upper.is_finite() && upper > 0.0);
        for sample in 0..=4096 {
            let t = -1.3 + 3.1 * f64::from(sample) / 4096.0;
            let point = curve.evaluate(t).unwrap();
            let deviation = (point - source_origin).dot(source_normal).abs();
            assert!(
                deviation <= upper,
                "t={t}, deviation={deviation}, upper={upper}"
            );
        }
    }

    #[test]
    fn other_plane_cone_intersection_regimes_are_rejected() {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let axis = Vector3::unit_z();
        assert!(
            AxisParallelPlaneConeHyperbola::from_carriers(
                origin,
                Vector3::unit_y(),
                origin,
                axis,
                0.5,
                1.0,
            )
            .is_none()
        );
        assert!(
            AxisParallelPlaneConeHyperbola::from_carriers(
                Point3::new(0.0, 1.0, 0.0),
                axis,
                origin,
                axis,
                0.5,
                1.0,
            )
            .is_none()
        );
        assert!(
            AxisParallelPlaneConeHyperbola::from_carriers(
                Point3::new(0.0, 1.0, 0.0),
                Vector3::unit_y(),
                origin,
                axis,
                0.5,
                0.0,
            )
            .is_none()
        );
    }
}
