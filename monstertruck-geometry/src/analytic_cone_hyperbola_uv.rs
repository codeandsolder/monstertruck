//! Continuous exact analytic UV of a plane–cone hyperbola.
//!
//! A spatial hyperbola is not generally rational in a cone's angle/axial
//! chart. This module provides an analytic parameter curve, not a falsely
//! exact STEP PCURVE or a fitted/sample-based substitute. Azimuth remains
//! unwrapped across the periodic seam.

use crate::analytic_plane_cone_hyperbola::{AxisParallelPlaneConeHyperbola, OrientedHyperbolaTrim};
use crate::prelude::{InnerSpace, Point2, Point3, Vector2, Vector3};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

/// The source conical surface's axial U / azimuth V chart for one hyperbola.
#[derive(Clone, Copy, Debug)]
pub struct ConeHyperbolaUv {
    hyperbola: AxisParallelPlaneConeHyperbola,
    x: Vector3,
    y: Vector3,
    azimuth_at_vertex: f64,
    /// Source reference-x and cone-axis dot before orthogonal projection.
    pub source_x_axis_misalignment: f64,
}

/// Source-directed UV range, never modulo reduced across the angle seam.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConeUvTrim {
    /// Axial millimetres, unwrapped radians.
    pub start: [f64; 2],
    /// Axial millimetres, unwrapped radians.
    pub end: [f64; 2],
    pub azimuth_span_rad: f64,
    pub increasing_azimuth: bool,
}

fn unit(vector: Vector3) -> bool {
    vector.magnitude2().is_finite()
        && (vector.magnitude2() - 1.0).abs()
            <= DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
}

fn finite(point: Point3) -> bool { [point.x, point.y, point.z].into_iter().all(f64::is_finite) }

impl ConeHyperbolaUv {
    /// Use the literal STEP cone x direction, not a guessed world axis.
    /// Reject reference directions outside the source axis orthogonality
    /// allowance and preserve the misalignment for error provenance.
    #[must_use = "The newly constructed value is returned rather than retained by this constructor."]
    pub fn new(hyperbola: AxisParallelPlaneConeHyperbola, source_x: Vector3) -> Option<Self> {
        let axis = hyperbola.cone_axis();
        if !unit(source_x) {
            return None;
        }
        let misalignment = source_x.dot(axis);
        if !misalignment.is_finite()
            || misalignment.abs() > DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
        {
            return None;
        }
        let x = (source_x - axis * misalignment).normalize();
        let y = axis.cross(x).normalize();
        let radial_at_vertex = hyperbola.plane_normal() * hyperbola.signed_plane_offset().signum();
        let azimuth_at_vertex = radial_at_vertex.dot(y).atan2(radial_at_vertex.dot(x));
        azimuth_at_vertex.is_finite().then_some(Self {
            hyperbola,
            x,
            y,
            azimuth_at_vertex,
            source_x_axis_misalignment: misalignment,
        })
    }

    /// The source-aligned orthonormal reference vectors used by the native
    /// cone-revolution surface (with the original STEP axis preserved).
    #[must_use = "The source-aligned cone reference axes are returned; dropping them loses chart orientation."]
    pub const fn reference_axes(self) -> (Vector3, Vector3) { (self.x, self.y) }

    /// Analytic UV at the existing 3D hyperbola parameter, with no wrap.
    ///
    /// Axial coordinate = (|d| cosh(t) - r0) / slope.
    /// Angle = `vertex_angle + sign(d) * atan(sinh(t))`.
    #[must_use = "The optional evaluate outcome must be inspected to distinguish a value from absence."]
    pub fn evaluate(self, t: f64) -> Option<Point2> {
        if !t.is_finite() {
            return None;
        }
        let d = self.hyperbola.signed_plane_offset();
        let u = d
            .abs()
            .mul_add(t.cosh(), -self.hyperbola.reference_radius())
            / self.hyperbola.slope();
        let v = d.signum().mul_add(t.sinh().atan(), self.azimuth_at_vertex);
        (u.is_finite() && v.is_finite()).then_some(Point2::new(u, v))
    }

    /// Analytic parameter tangent; azimuth is strictly monotone.
    #[must_use = "The optional derivative outcome must be inspected to distinguish a value from absence."]
    pub fn derivative(self, t: f64) -> Option<Vector2> {
        if !t.is_finite() {
            return None;
        }
        let d = self.hyperbola.signed_plane_offset();
        let du = d.abs() * t.sinh() / self.hyperbola.slope();
        let dv = d.signum() / t.cosh();
        (du.is_finite() && dv.is_finite()).then_some(Vector2::new(du, dv))
    }

    /// Evaluate the source cone's independent generator-line revolution
    /// chart, not a projection onto the proposed healed hyperbola.
    #[must_use = "The source cone point is returned, if defined, without changing the cone."]
    pub fn cone_evaluate(self, uv: Point2) -> Option<Point3> {
        if !uv.x.is_finite() || !uv.y.is_finite() {
            return None;
        }
        let radius = self
            .hyperbola
            .slope()
            .mul_add(uv.x, self.hyperbola.reference_radius());
        if !radius.is_finite() || radius <= 0.0 {
            return None;
        }
        let radial = self.x * uv.y.cos() + self.y * uv.y.sin();
        let point =
            self.hyperbola.reference_origin() + self.hyperbola.cone_axis() * uv.x + radial * radius;
        finite(point).then_some(point)
    }

    /// The finite branch's angle is monotone and spans strictly less than
    /// pi radians, including if its source-cone chart seam lies inside it.
    #[must_use = "The optional trim outcome must be inspected to distinguish a value from absence."]
    pub fn trim(self, trim: OrientedHyperbolaTrim) -> Option<ConeUvTrim> {
        let first = self.evaluate(trim.start_parameter)?;
        let last = self.evaluate(trim.end_parameter)?;
        let span = (last.y - first.y).abs();
        if !span.is_finite() || span == 0.0 || span >= std::f64::consts::PI {
            return None;
        }
        Some(ConeUvTrim {
            start: [first.x, first.y],
            end: [last.x, last.y],
            azimuth_span_rad: span,
            increasing_azimuth: last.y > first.y,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn chart_reproduces_independent_cone_geometry() {
        let h = example();
        let chart = ConeHyperbolaUv::new(h, Vector3::unit_x()).unwrap();
        for i in -512..=512 {
            let t = f64::from(i) / 256.0;
            let uv = chart.evaluate(t).unwrap();
            let generated = chart.cone_evaluate(uv).unwrap();
            assert!((generated - h.evaluate(t).unwrap()).magnitude() < 1.0e-12);
            assert!(chart.derivative(t).unwrap().y > 0.0);
        }
    }

    #[test]
    fn seam_crossing_has_continuous_unwrapped_azimuth() {
        let h = example();
        let source_x = Vector3::new(0.01, 0.999_949_998_749_937_5, 0.0);
        let chart = ConeHyperbolaUv::new(h, source_x).unwrap();
        let trim = h
            .oriented_trim(h.evaluate(-1.3).unwrap(), h.evaluate(1.3).unwrap(), 1.0e-12)
            .unwrap();
        let range = chart.trim(trim).unwrap();
        assert!(range.azimuth_span_rad < std::f64::consts::PI);
        assert!(range.increasing_azimuth);
        assert!(range.start[1] < -std::f64::consts::PI);
        assert!(range.end[1] > -std::f64::consts::PI);
        for i in 0..=256 {
            let t = f64::from(i).mul_add(2.6 / 256.0, -1.3);
            let p = chart.cone_evaluate(chart.evaluate(t).unwrap()).unwrap();
            assert!((p - h.evaluate(t).unwrap()).magnitude() < 1.0e-12);
        }
    }

    #[test]
    fn source_reverse_trim_reverses_azimuth_without_wrapping() {
        let h = example();
        let chart = ConeHyperbolaUv::new(h, Vector3::unit_x()).unwrap();
        let trim = h
            .oriented_trim(h.evaluate(1.1).unwrap(), h.evaluate(-0.7).unwrap(), 1.0e-12)
            .unwrap();
        let range = chart.trim(trim).unwrap();
        assert!(!range.increasing_azimuth);
    }

    #[test]
    fn opposite_signed_plane_offset_reverses_angular_chart_direction() {
        let h = AxisParallelPlaneConeHyperbola::from_carriers(
            Point3::new(0.0, 0.3, 2.0),
            Vector3::new(0.0, -1.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            0.5,
            1.5,
        )
        .unwrap();
        let source_x = Vector3::new(1.0, 0.0, 1.0e-10).normalize();
        let chart = ConeHyperbolaUv::new(h, source_x).unwrap();
        assert!(chart.source_x_axis_misalignment > 0.0);
        let trim = h
            .oriented_trim(h.evaluate(-1.1).unwrap(), h.evaluate(0.9).unwrap(), 1.0e-12)
            .unwrap();
        let range = chart.trim(trim).unwrap();
        assert!(!range.increasing_azimuth);
        for sample in 0..=100 {
            let t = f64::from(sample).mul_add(2.0 / 100.0, -1.1);
            let uv = chart.evaluate(t).unwrap();
            assert!(chart.derivative(t).unwrap().y < 0.0);
            assert!(
                (chart.cone_evaluate(uv).unwrap() - h.evaluate(t).unwrap()).magnitude() < 1.0e-12
            );
        }
    }

    #[test]
    fn invalid_source_frame_and_nonfinite_parameters_fail() {
        let h = example();
        assert!(ConeHyperbolaUv::new(h, Vector3::unit_z()).is_none());
        assert!(ConeHyperbolaUv::new(h, Vector3::new(f64::NAN, 0.0, 0.0),).is_none());
        let chart = ConeHyperbolaUv::new(h, Vector3::unit_x()).unwrap();
        assert!(chart.evaluate(f64::INFINITY).is_none());
        assert!(chart.evaluate(1024.0).is_none());
        assert!(chart.cone_evaluate(Point2::new(-1.0, 0.0)).is_none());
    }
}
