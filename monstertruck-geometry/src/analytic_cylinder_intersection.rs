//! Analytic branches of two nonparallel circular cylinder intersections.
//!
//! The first cylinder's angular chart parametrizes its axial generators.
//! Intersecting each generator with the second cylinder is a quadratic;
//! a strictly positive discriminant on the entire requested interval
//! separates its two roots. This constructs new geometry, not a claim
//! that any approximate source trim is equivalent to it.

use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use std::f64::consts::{PI, TAU};
mod carrier;
mod leader;
mod source;
pub use carrier::{certified_circular_cylinder, project_positive_homogeneous_point};
pub use leader::IntersectionLeader;
pub use source::{BoundedSourceIntersection, SourceIntersectionProposal};

/// A circular cylindrical carrier in physical coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylinderCarrier {
    origin: Point3,
    axis: Vector3,
    radius: f64,
}

fn finite_point(point: Point3) -> bool { [point.x, point.y, point.z].iter().all(|v| v.is_finite()) }

impl CylinderCarrier {
    /// Validate a finite, positive-radius cylinder and its axial direction.
    /// Only a numerically unit direction is accepted.
    #[must_use = "inspect whether the original carrier has a valid frame"]
    pub fn new(origin: Point3, axis: Vector3, radius: f64) -> Option<Self> {
        let length = axis.magnitude();
        if !finite_point(origin)
            || !length.is_finite()
            || (length - 1.0).abs() > DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
            || !radius.is_finite()
            || radius <= 0.0
        {
            return None;
        }
        Some(Self {
            origin,
            axis: axis / length,
            radius,
        })
    }

    /// Certify a rational quadratic circular extrusion, accepting either
    /// parameter-axis order. The input must describe the complete surface
    /// domain, not merely an approximate or local control net.
    #[must_use = "inspect the structurally recognized cylindrical support"]
    pub fn from_homogeneous_surface(mut surface: BsplineSurface<Vector4>) -> Option<Self> {
        if surface.degrees() == (1, 2) {
            surface.swap_axes();
        }
        let (origin, axis, radius, _) = certified_circular_cylinder(&surface)?;
        Self::new(origin, axis, radius)
    }

    /// Euclidean distance to the unbounded cylinder.
    #[must_use = "the measured carrier deviation is returned"]
    pub fn distance(self, point: Point3) -> f64 {
        let delta = point - self.origin;
        (delta.cross(self.axis).magnitude() - self.radius).abs()
    }
}

/// The sign of the secondary cylinder's transverse root, relative to the
/// projection of the primary cylinder's axial direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntersectionBranch {
    Negative,
    Positive,
}

impl IntersectionBranch {
    const fn sign(self) -> f64 {
        match self {
            Self::Negative => -1.0,
            Self::Positive => 1.0,
        }
    }
}

/// Finite directed arc of the explicit, separated cylinder intersection.
///
/// Construction certifies the entire closed angular domain rather than
/// accepting a series of successfully solved sample points. Tangent roots,
/// parallel axes and non-finite or badly conditioned frames are rejected.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CylinderIntersectionArc {
    primary: CylinderCarrier,
    secondary: CylinderCarrier,
    radial_x: Vector3,
    radial_y: Vector3,
    normal: Vector3,
    axial_projection: Vector3,
    axis_sine: f64,
    sweep: f64,
    branch: IntersectionBranch,
    min_root: f64,
    arithmetic_margin: f64,
}

fn max_abs_trig(constant: f64, cosine: f64, sine: f64, sweep: f64) -> f64 {
    let value = |angle: f64| sine.mul_add(angle.sin(), cosine.mul_add(angle.cos(), constant));
    let mut result = value(0.0).abs().max(value(sweep).abs());
    let phase = sine.atan2(cosine);
    let lo = sweep.min(0.0);
    let hi = sweep.max(0.0);
    let padding = 64.0 * f64::EPSILON * TAU;
    // The interval is anchored at zero and never exceeds one revolution.
    // The stationary points are phase+k*pi, both positive and negative.
    for k in -3..=3 {
        let angle = f64::from(k).mul_add(PI, phase);
        if angle >= lo - padding && angle <= hi + padding {
            result = result.max(value(angle).abs());
        }
    }
    result
}

impl CylinderIntersectionArc {
    /// Construct a branch over a directed angular interval beginning at
    /// the radial direction of `start` on the primary cylinder.
    ///
    /// The original start vertex must lie on the selected root within the
    /// explicit endpoint tolerance. The branch sign does not permit a root
    /// switch in the interior. No STEP or source geometry is changed.
    #[must_use = "inspect the validated analytic intersection or unsupported result"]
    pub fn new(
        primary: CylinderCarrier,
        secondary: CylinderCarrier,
        start: Point3,
        sweep: f64,
        branch: IntersectionBranch,
        endpoint_tolerance: f64,
    ) -> Option<Self> {
        if !finite_point(start)
            || !sweep.is_finite()
            || sweep.abs() <= 0.0
            || sweep.abs() > TAU
            || !endpoint_tolerance.is_finite()
            || endpoint_tolerance <= 0.0
            || primary.distance(start) > endpoint_tolerance
            || secondary.distance(start) > endpoint_tolerance
        {
            return None;
        }
        let cross = primary.axis.cross(secondary.axis);
        let axis_sine = cross.magnitude();
        if !axis_sine.is_finite() || axis_sine <= f64::EPSILON.sqrt() {
            return None;
        }
        let normal = cross / axis_sine;
        let axial_projection = secondary.axis.cross(normal).normalize();
        let offset = start - primary.origin;
        let radial = offset - primary.axis * offset.dot(primary.axis);
        let length = radial.magnitude();
        if !length.is_finite() || length <= 0.0 {
            return None;
        }
        let radial_x = radial / length;
        let radial_y = primary.axis.cross(radial_x).normalize();
        // Moving the primary origin along its own axis does not change the
        // carrier and avoids cancellation when the source chart is remote.
        let primary = CylinderCarrier {
            origin: primary.origin + primary.axis * offset.dot(primary.axis),
            ..primary
        };
        let relative = primary.origin - secondary.origin;
        let scale = [
            primary.origin.x.abs(),
            primary.origin.y.abs(),
            primary.origin.z.abs(),
            secondary.origin.x.abs(),
            secondary.origin.y.abs(),
            secondary.origin.z.abs(),
            primary.radius,
            secondary.radius,
            relative.magnitude(),
            1.0,
        ]
        .into_iter()
        .fold(1.0_f64, f64::max);
        let radial_margin = 512.0 * f64::EPSILON * scale;
        let max_distance = max_abs_trig(
            relative.dot(normal),
            primary.radius * radial_x.dot(normal),
            primary.radius * radial_y.dot(normal),
            sweep,
        ) + radial_margin;
        let gap = secondary.radius - max_distance;
        if !gap.is_finite() || gap <= 0.0 {
            return None;
        }
        let min_root = (gap * (secondary.radius + max_distance)).sqrt();
        if !min_root.is_finite() || min_root <= radial_margin {
            return None;
        }
        // Finite-f64 conditioning allowance, not outward-rounded interval
        // arithmetic. Account for the inverse transverse root derivative
        // and the nearly-parallel axial projection.
        let arithmetic_margin = radial_margin * (1.0 + secondary.radius / min_root) / axis_sine;
        if !arithmetic_margin.is_finite() || arithmetic_margin > endpoint_tolerance {
            return None;
        }
        let result = Self {
            primary,
            secondary,
            radial_x,
            radial_y,
            normal,
            axial_projection,
            axis_sine,
            sweep,
            branch,
            min_root,
            arithmetic_margin,
        };
        ((result.evaluate(0.0)? - start).magnitude() <= endpoint_tolerance).then_some(result)
    }

    /// Directed angular sweep of the bounded arc, in radians.
    #[must_use = "the angular interval is part of the arc identity"]
    pub const fn sweep(self) -> f64 { self.sweep }

    /// Selected separated transverse branch in this ordered carrier frame.
    #[must_use = "the branch is part of the intersection identity"]
    pub const fn branch(self) -> IntersectionBranch { self.branch }

    /// Lower bound on the magnitude of the separated transverse root.
    #[must_use = "the branch margin quantifies its conditioning"]
    pub const fn minimum_root(self) -> f64 { self.min_root }

    /// Numerical conditioning allowance in world-coordinate length units.
    #[must_use = "the reported conditioning margin is part of the evidence"]
    pub const fn arithmetic_margin(self) -> f64 { self.arithmetic_margin }

    /// Evaluate the explicit analytic intersection and its angular first
    /// and second derivatives at a parameter in the directed arc.
    #[must_use = "the evaluated point and derivatives are returned"]
    pub fn jet(self, angle: f64) -> Option<(Point3, Vector3, Vector3)> {
        if !angle.is_finite() || angle < self.sweep.min(0.0) || angle > self.sweep.max(0.0) {
            return None;
        }
        let (sin, cos) = angle.sin_cos();
        let radial = (self.radial_x * cos + self.radial_y * sin) * self.primary.radius;
        let first = (self.radial_y * cos - self.radial_x * sin) * self.primary.radius;
        let second = -radial;
        let generator = self.primary.origin + radial;
        let offset = generator - self.secondary.origin;
        let distance = offset.dot(self.normal);
        let distance_first = first.dot(self.normal);
        let distance_second = second.dot(self.normal);
        let radius = self.secondary.radius;
        let root = ((radius - distance.abs()) * (radius + distance.abs())).sqrt();
        if !root.is_finite() || root <= 0.0 {
            return None;
        }
        let root_first = -distance * distance_first / root;
        let root_second = -distance.mul_add(distance_second, distance_first.powi(2)) / root
            - distance * distance * distance_first * distance_first / root.powi(3);
        let axial = |root: f64, projection: f64| {
            self.branch.sign().mul_add(root, -projection) / self.axis_sine
        };
        let point = generator + self.primary.axis * axial(root, offset.dot(self.axial_projection));
        let derivative =
            first + self.primary.axis * axial(root_first, first.dot(self.axial_projection));
        let second_derivative =
            second + self.primary.axis * axial(root_second, second.dot(self.axial_projection));
        let finite = |v: Vector3| [v.x, v.y, v.z].iter().all(|x| x.is_finite());
        (finite_point(point) && finite(derivative) && finite(second_derivative)).then_some((
            point,
            derivative,
            second_derivative,
        ))
    }

    /// Evaluate the selected root without a numerical surface intersection.
    #[must_use = "the evaluated point is returned"]
    pub fn evaluate(self, angle: f64) -> Option<Point3> { self.jet(angle).map(|jet| jet.0) }

    /// Verify an independently authored terminal vertex without changing it.
    #[must_use = "inspect the directed endpoint residual"]
    pub fn end_residual(self, end: Point3, tolerance: f64) -> Option<f64> {
        let residual = (self.evaluate(self.sweep)? - end).magnitude();
        (residual.is_finite() && tolerance.is_finite() && tolerance >= 0.0 && residual <= tolerance)
            .then_some(residual)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn carrier(origin: Point3, axis: Vector3, radius: f64) -> CylinderCarrier {
        CylinderCarrier::new(origin, axis, radius).expect("valid cylinder")
    }

    #[test]
    fn skew_cylinder_intersection_satisfies_both_carriers_and_derivatives() {
        let primary = carrier(Point3::new(2., -1., 3.), Vector3::unit_z(), 0.4);
        let secondary = carrier(
            Point3::new(2., -1., 3.),
            Vector3::new(1., 0., 1.).normalize(),
            1.2,
        );
        // At angle zero, the primary generator is at x=2.4, y=-1.
        // Its positive secondary transverse root is z=3.4+sqrt(2)*1.2.
        let start = Point3::new(2.4, -1., 3.4 + 2.0_f64.sqrt() * 1.2);
        let arc = CylinderIntersectionArc::new(
            primary,
            secondary,
            start,
            TAU,
            IntersectionBranch::Positive,
            1e-9,
        )
        .expect("two separated branches over a full turn");
        for i in 1..200 {
            let angle = TAU * f64::from(i) / 200.;
            let (point, first, second) = arc.jet(angle).expect("valid interior jet");
            assert!(primary.distance(point) < 1e-12 && secondary.distance(point) < 1e-12);
            let h = 1e-5;
            let (left, left_first, _) = arc.jet(angle - h).expect("left");
            let (right, right_first, _) = arc.jet(angle + h).expect("right");
            assert!(((right - left) / (2. * h) - first).magnitude() < 1e-8);
            assert!(((right_first - left_first) / (2. * h) - second).magnitude() < 1e-8);
        }
        assert!(arc.end_residual(start, 1e-9).is_some());
    }

    #[test]
    fn interior_discriminant_extremum_cannot_be_hidden_by_valid_endpoints() {
        let primary = carrier(Point3::origin(), Vector3::unit_z(), 1.);
        let secondary = carrier(Point3::origin(), Vector3::unit_x(), 0.5);
        let start = Point3::new(1., 0., 0.5);
        assert!(
            CylinderIntersectionArc::new(
                primary,
                secondary,
                start,
                PI,
                IntersectionBranch::Positive,
                1e-9
            )
            .is_none()
        );
        assert!(
            CylinderIntersectionArc::new(
                primary,
                secondary,
                start,
                -0.25,
                IntersectionBranch::Positive,
                1e-9
            )
            .is_some()
        );
        assert!(
            CylinderIntersectionArc::new(
                primary,
                secondary,
                start,
                -0.25,
                IntersectionBranch::Negative,
                1e-9
            )
            .is_none()
        );
    }

    #[test]
    fn parallel_axes_and_tangent_roots_are_refused() {
        let primary = carrier(Point3::origin(), Vector3::unit_z(), 1.);
        let start = Point3::new(1., 0., 0.);
        let parallel = carrier(Point3::origin(), Vector3::unit_z(), 1.);
        let tangent = carrier(Point3::origin(), Vector3::unit_x(), 1.);
        assert!(
            CylinderIntersectionArc::new(
                primary,
                parallel,
                start,
                PI,
                IntersectionBranch::Positive,
                1e-9
            )
            .is_none()
        );
        assert!(
            CylinderIntersectionArc::new(
                primary,
                tangent,
                Point3::new(1., 0., 1.),
                PI,
                IntersectionBranch::Positive,
                1e-9
            )
            .is_none()
        );
    }
}
