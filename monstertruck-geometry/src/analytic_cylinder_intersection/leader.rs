//! Certified cubic Hermite leader for the analytical intersection.
//!
//! STEP uses a finite spline as the preferred 3D representation of a
//! general intersection. Its approximation error is bounded independently
//! of subsequent surface projection or the source's approximation error.

use super::CylinderIntersectionArc;
use crate::prelude::*;

/// A piecewise cubic approximation with a complete-domain error bound.
#[derive(Clone, Debug)]
pub struct IntersectionLeader {
    pub curve: BsplineCurve<Point3>,
    pub segments: u32,
    pub maximum_error: f64,
    pub arithmetic_margin: f64,
}

impl CylinderIntersectionArc {
    /// Bound the Euclidean norm of the fourth angular derivative.
    ///
    /// With q² = R² - d², successive differentiation bounds q derivatives
    /// recursively using the previously certified positive minimum q.
    /// All derivatives of the sinusoidal part d are bounded by the
    /// primary cylinder radius. The zeroth derivative has its own bound.
    fn fourth_derivative_upper(self) -> Option<f64> {
        let radius = self.primary.radius;
        let distance = self.secondary.radius;
        let root = self.min_root;
        let denominator = 2.0 * root;
        let first = distance * radius / root;
        let second = (distance.mul_add(radius, radius.powi(2)) + first.powi(2)) / root;
        let third = (6.0_f64.mul_add(
            first * second,
            (2.0 * distance).mul_add(radius, 6.0 * radius.powi(2)),
        )) / denominator;
        // Leibniz coefficients for (d²)'''' and (q²)''''.
        let fourth = (6.0_f64.mul_add(
            second.powi(2),
            8.0_f64.mul_add(
                first * third,
                (2.0 * distance).mul_add(radius, 14.0 * radius.powi(2)),
            ),
        )) / denominator;
        let upper = radius + (fourth + radius) / self.axis_sine;
        (upper.is_finite() && upper > 0.0).then_some(upper)
    }

    /// Construct a finite, directed cubic leader on the normalized
    /// parameter interval [0,1].
    ///
    /// The Hermite interpolation remainder is bounded by M*h^4/384,
    /// where M bounds the fourth derivative over the entire arc.
    /// This includes a finite-f64 arithmetic allowance, not formal
    /// outward-rounded interval arithmetic. Failure to meet the explicit
    /// error or segment budget returns None.
    #[must_use = "inspect the bounded intersection leader before exporting it"]
    pub fn cubic_leader(self, tolerance: f64, max_segments: u32) -> Option<IntersectionLeader> {
        if !tolerance.is_finite() || tolerance <= 0.0 || max_segments == 0 {
            return None;
        }
        let fourth = self.fourth_derivative_upper()?;
        // The analytical evaluation allowance also covers cancellation in
        // the Hermite handles, with a separate allowance for derivative
        // conditioning as the transverse root approaches zero.
        let speed = self.primary.radius
            * (1.0 + (1.0 + self.secondary.radius / self.min_root) / self.axis_sine);
        let conditioning = (1.0 + self.secondary.radius / self.min_root).powi(2);
        let base_margin = 4.0 * self.arithmetic_margin;
        if !base_margin.is_finite() || base_margin >= tolerance {
            return None;
        }
        let mut segments = 1_u32;
        let (maximum_error, margin) = loop {
            let width = self.sweep / f64::from(segments);
            let margin = width.abs().mul_add(
                self.arithmetic_margin
                    .mul_add(conditioning, 256.0 * f64::EPSILON * speed),
                base_margin,
            );
            let error = (fourth / 384.0).mul_add(width.powi(4), margin);
            if error.is_finite() && error <= tolerance {
                break (error, margin);
            }
            segments = segments.checked_mul(2)?;
            if segments > max_segments {
                return None;
            }
        };
        let mut controls = Vec::new();
        let mut knots = vec![0.0; 4];
        let width = self.sweep / f64::from(segments);
        let (mut start, mut tangent, _) = self.jet(0.0)?;
        controls.push(start);
        for segment in 1..=segments {
            let parameter = f64::from(segment) / f64::from(segments);
            let (end, derivative, _) = self.jet(self.sweep * parameter)?;
            controls.extend([
                start + tangent * (width / 3.0),
                end - derivative * (width / 3.0),
                end,
            ]);
            if segment < segments {
                knots.extend([parameter; 3]);
            }
            start = end;
            tangent = derivative;
        }
        knots.extend([1.0; 4]);
        let curve = BsplineCurve::new(KnotVector::try_from(knots).ok()?, controls).ok()?;
        Some(IntersectionLeader {
            curve,
            segments,
            maximum_error,
            arithmetic_margin: margin,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analytic_cylinder_intersection::{CylinderCarrier, IntersectionBranch};

    #[test]
    fn bounded_cubic_leader_matches_directed_analytical_arc() {
        let primary =
            CylinderCarrier::new(Point3::origin(), Vector3::unit_z(), 0.4).expect("primary");
        let secondary =
            CylinderCarrier::new(Point3::origin(), Vector3::unit_x(), 1.2).expect("secondary");
        for sweep in [-3.0, 3.0] {
            let arc = CylinderIntersectionArc::new(
                primary,
                secondary,
                Point3::new(0.4, 0.0, 1.2),
                sweep,
                IntersectionBranch::Positive,
                1e-9,
            )
            .expect("arc");
            let leader = arc.cubic_leader(1e-9, 8192).expect("bounded leader");
            assert!(leader.maximum_error <= 1e-9);
            assert!(leader.segments < 8192);
            for i in 0..=4096 {
                let t = f64::from(i) / 4096.0;
                let expected = arc.evaluate(t * sweep).expect("point");
                assert!((leader.curve.evaluate(t) - expected).magnitude() <= leader.maximum_error);
            }
            assert!(arc.cubic_leader(1e-20, 8192).is_none());
            assert!(arc.cubic_leader(1e-9, 1).is_none());
        }
    }
}
