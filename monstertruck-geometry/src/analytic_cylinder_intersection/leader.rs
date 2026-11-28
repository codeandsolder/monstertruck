//! Certified compact Hermite leader for the analytical intersection.
//!
//! STEP uses a finite spline as the preferred 3D representation of a
//! general intersection. Its approximation error is bounded independently
//! of subsequent surface projection or the source's approximation error.

use super::CylinderIntersectionArc;
use crate::prelude::*;

/// A piecewise polynomial approximation with a complete-domain error bound.
#[derive(Clone, Debug)]
pub struct IntersectionLeader {
    pub curve: BsplineCurve<Point3>,
    pub segments: u32,
    pub maximum_error: f64,
    pub arithmetic_margin: f64,
}

#[derive(Clone, Copy)]
enum HermiteOrder {
    Cubic,
    Quintic,
}

impl HermiteOrder {
    const fn degree(self) -> usize {
        match self {
            Self::Cubic => 3,
            Self::Quintic => 5,
        }
    }

    const fn remainder_divisor(self) -> f64 {
        // (2m)! * 4^m for m-fold Hermite conditions at each endpoint.
        match self {
            Self::Cubic => 384.0,
            Self::Quintic => 46080.0,
        }
    }
}

#[derive(Clone, Copy)]
struct LeaderPlan {
    order: HermiteOrder,
    segments: u32,
    maximum_error: f64,
    margin: f64,
}

impl LeaderPlan {
    fn control_count(self) -> u64 {
        u64::from(self.segments)
            * (if matches!(self.order, HermiteOrder::Cubic) {
                3
            } else {
                5
            })
            + 1
    }
}

/// Leibniz coefficients needed for up to six angular derivatives.
const BINOMIAL: [[f64; 7]; 7] = [
    [1., 0., 0., 0., 0., 0., 0.],
    [1., 1., 0., 0., 0., 0., 0.],
    [1., 2., 1., 0., 0., 0., 0.],
    [1., 3., 3., 1., 0., 0., 0.],
    [1., 4., 6., 4., 1., 0., 0.],
    [1., 5., 10., 10., 5., 1., 0.],
    [1., 6., 15., 20., 15., 6., 1.],
];

impl CylinderIntersectionArc {
    /// Bound the norms of angular derivatives through order six.
    ///
    /// With q² = R² - d², differentiation gives
    /// 2q q^(n) = -(d²)^(n) - Σ binom(n,j) q^(j) q^(n-j).
    /// The complete-domain positive minimum q is already certified.
    /// The sinusoid d has |d| <= R and |d^(j)| <= r for j > 0.
    fn derivative_upper(self, order: usize) -> Option<f64> {
        let radius = self.primary.radius;
        let distance = self.secondary.radius;
        let root = self.min_root;
        let mut q = [0.0; 7];
        for n in 1..=order {
            let mut numerator = 2.0 * distance * radius;
            for j in 1..n {
                numerator =
                    BINOMIAL[n][j].mul_add(radius.mul_add(radius, q[j] * q[n - j]), numerator);
            }
            // Positive arithmetic only. Include finite-f64 evaluation
            // slack in each bound before using it in the next recurrence.
            q[n] = numerator / (2.0 * root) * 64.0_f64.mul_add(f64::EPSILON, 1.0);
        }
        let upper = radius + (q[order] + radius) / self.axis_sine;
        (upper.is_finite() && upper > 0.0).then_some(upper)
    }

    fn leader_plan(
        self,
        order: HermiteOrder,
        tolerance: f64,
        max_segments: u32,
    ) -> Option<LeaderPlan> {
        let derivative_order = match order {
            HermiteOrder::Cubic => 4,
            HermiteOrder::Quintic => 6,
        };
        let derivative = self.derivative_upper(derivative_order)?;
        let speed = self.derivative_upper(1)?;
        let ratio = self.secondary.radius / self.min_root;
        let conditioning = (1.0 + ratio).powi(2);
        let base_margin = 4.0 * self.arithmetic_margin;
        let first_margin = self
            .arithmetic_margin
            .mul_add(conditioning, 256.0 * f64::EPSILON * speed);
        // q'' = -(d'^2+d d'')/q - d² d'^2/q³. Its partial
        // derivatives in d, d' and d'' are bounded by a polynomial
        // in R/q and r/q; the following positive envelope also
        // accounts for projection onto the primary axial direction.
        let second_margin = if matches!(order, HermiteOrder::Quintic) {
            let envelope =
                8.0 * (1.0 + ratio).powi(3) * (1.0 + self.primary.radius / self.min_root).powi(2);
            self.arithmetic_margin
                .mul_add(envelope, 512.0 * f64::EPSILON * self.derivative_upper(2)?)
        } else {
            0.0
        };
        let measure = |segments| {
            let width = (self.sweep / f64::from(segments)).abs();
            let margin = (width * width / 20.0)
                .mul_add(second_margin, width.mul_add(first_margin, base_margin));
            let remainder = match order {
                HermiteOrder::Cubic => width.powi(4),
                HermiteOrder::Quintic => width.powi(6),
            };
            let maximum_error = (derivative / order.remainder_divisor()).mul_add(remainder, margin);
            LeaderPlan {
                order,
                segments,
                maximum_error,
                margin,
            }
        };
        let fits =
            |plan: LeaderPlan| plan.maximum_error.is_finite() && plan.maximum_error <= tolerance;
        if !fits(measure(max_segments)) {
            return None;
        }
        // The bound is monotone in the uniform interval width. Find the
        // minimum feasible integral segment count, not just a power of two.
        let (mut low, mut high) = (1, max_segments);
        while low < high {
            let middle = low + (high - low) / 2;
            if fits(measure(middle)) {
                high = middle;
            } else {
                low = middle + 1;
            }
        }
        Some(measure(low))
    }

    fn build_leader(self, plan: LeaderPlan) -> Option<IntersectionLeader> {
        let degree = plan.order.degree();
        let mut controls = Vec::new();
        let mut knots = vec![0.0; degree + 1];
        let width = self.sweep / f64::from(plan.segments);
        let (mut start, mut tangent, mut acceleration) = self.jet(0.0)?;
        controls.push(start);
        for segment in 1..=plan.segments {
            let parameter = f64::from(segment) / f64::from(plan.segments);
            let (end, derivative, second) = self.jet(self.sweep * parameter)?;
            match plan.order {
                HermiteOrder::Cubic => controls.extend([
                    start + tangent * (width / 3.0),
                    end - derivative * (width / 3.0),
                    end,
                ]),
                HermiteOrder::Quintic => controls.extend([
                    start + tangent * (width / 5.0),
                    start + tangent * (2.0 * width / 5.0) + acceleration * (width * width / 20.0),
                    end - derivative * (2.0 * width / 5.0) + second * (width * width / 20.0),
                    end - derivative * (width / 5.0),
                    end,
                ]),
            }
            if segment < plan.segments {
                knots.extend(std::iter::repeat_n(parameter, degree));
            }
            start = end;
            tangent = derivative;
            acceleration = second;
        }
        knots.extend(std::iter::repeat_n(1.0, degree + 1));
        let curve = BsplineCurve::new(KnotVector::try_from(knots).ok()?, controls).ok()?;
        Some(IntersectionLeader {
            curve,
            segments: plan.segments,
            maximum_error: plan.maximum_error,
            arithmetic_margin: plan.margin,
        })
    }

    /// Construct the smallest certified cubic or quintic Hermite leader
    /// on the normalized parameter interval [0,1].
    ///
    /// The uniform per-span interpolation remainder is bounded by
    /// M*h^4/384 or M*h^6/46080, respectively. The latter matches both
    /// derivatives available from the analytic jet. For a vector curve,
    /// apply the scalar remainder in any unit projection direction.
    /// The bound includes finite-f64 arithmetic allowances, not formal
    /// outward-rounded interval arithmetic. No feasible plan means None.
    #[must_use = "inspect the bounded intersection leader before exporting it"]
    pub fn bounded_leader(self, tolerance: f64, max_segments: u32) -> Option<IntersectionLeader> {
        if !tolerance.is_finite() || tolerance <= 0.0 || max_segments == 0 {
            return None;
        }
        let plan = [HermiteOrder::Cubic, HermiteOrder::Quintic]
            .into_iter()
            .filter_map(|order| self.leader_plan(order, tolerance, max_segments))
            .min_by_key(|plan| plan.control_count())?;
        self.build_leader(plan)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analytic_cylinder_intersection::{CylinderCarrier, IntersectionBranch};
    use std::f64::consts::PI;

    #[test]
    fn compact_plan_uses_the_smallest_feasible_integral_segment_count() {
        let primary = CylinderCarrier::new(Point3::origin(), Vector3::z(), 0.4).expect("primary");
        let secondary =
            CylinderCarrier::new(Point3::origin(), Vector3::x(), 1.2).expect("secondary");
        let arc = CylinderIntersectionArc::new(
            primary,
            secondary,
            Point3::new(0.4, 0.0, 1.2),
            3.0,
            IntersectionBranch::Positive,
            1e-9,
        )
        .expect("arc");
        let tolerance = 1e-9;
        let cubic = arc
            .leader_plan(HermiteOrder::Cubic, tolerance, 8192)
            .expect("cubic");
        let quintic = arc
            .leader_plan(HermiteOrder::Quintic, tolerance, 8192)
            .expect("quintic");
        assert!(quintic.control_count() < cubic.control_count());
        for plan in [cubic, quintic] {
            assert!(plan.segments > 1);
            assert!(
                arc.leader_plan(plan.order, tolerance, plan.segments - 1)
                    .is_none()
            );
            let spline = arc.build_leader(plan).expect("spline");
            assert_eq!(spline.curve.degree(), plan.order.degree());
            for i in 0..=8192 {
                let t = f64::from(i) / 8192.0;
                let expected = arc.evaluate(t * arc.sweep()).expect("analytical point");
                assert!(spline.curve.evaluate(t).distance(expected) <= plan.maximum_error);
            }
        }
        let chosen = arc
            .bounded_leader(tolerance, 8192)
            .expect("selected leader");
        assert_eq!(chosen.curve.degree(), 5);
        assert_eq!(chosen.segments, quintic.segments);
    }

    #[test]
    fn compact_leaders_preserve_skew_separated_branches_near_tangency() {
        let origin = Point3::new(2., -1., 3.);
        let primary = CylinderCarrier::new(origin, Vector3::z(), 0.4).expect("primary");
        let secondary = CylinderCarrier::new(origin, Vector3::new(1., 0., 1.).normalize(), 0.41)
            .expect("secondary");
        for (branch, sign) in [
            (IntersectionBranch::Positive, 1.0),
            (IntersectionBranch::Negative, -1.0),
        ] {
            let start =
                origin + Vector3::new(0.4, 0.0, 0.41_f64.mul_add(sign * 2.0_f64.sqrt(), 0.4));
            for sweep in [-PI, PI] {
                let arc =
                    CylinderIntersectionArc::new(primary, secondary, start, sweep, branch, 1e-8)
                        .expect("skew branch");
                let leader = arc.bounded_leader(1e-8, 8192).expect("compact skew leader");
                for i in 0..=4096 {
                    let t = f64::from(i) / 4096.0;
                    let point = arc.evaluate(t * sweep).expect("point");
                    assert!(leader.curve.evaluate(t).distance(point) <= leader.maximum_error);
                }
            }
        }
    }

    #[test]
    fn bounded_leader_matches_directed_analytical_arc() {
        let primary = CylinderCarrier::new(Point3::origin(), Vector3::z(), 0.4).expect("primary");
        let secondary =
            CylinderCarrier::new(Point3::origin(), Vector3::x(), 1.2).expect("secondary");
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
            let leader = arc.bounded_leader(1e-9, 8192).expect("bounded leader");
            assert!(leader.maximum_error <= 1e-9);
            assert!(leader.segments < 8192);
            for i in 0..=4096 {
                let t = f64::from(i) / 4096.0;
                let expected = arc.evaluate(t * sweep).expect("point");
                assert!((leader.curve.evaluate(t) - expected).norm() <= leader.maximum_error);
            }
            assert!(arc.bounded_leader(1e-20, 8192).is_none());
            assert!(arc.bounded_leader(1e-9, 1).is_none());
        }
    }
}
