//! Regular intersections of two affine-radius revolution carriers. Fix an
//! azimuth on either carrier; the other carrier's implicit equation becomes
//! a quadratic in axial height along the first carrier's generator.
//! A source path selects a root of that quadratic, not a guessed conic.
//! Interval bounds separate roots continuously; native cubic approximants
//! retain the same root and are independently bounded on both carriers.

use super::super::super::sample_fit::normal_basis;
use super::SurfaceCarrier;
use crate::{BsplineCurve, Curve, InnerProduct, KnotVector, NurbsCurve, Point3, Vector3};
use monstertruck_geometry::{
    analytic_carrier::AnalyticSupport,
    analytic_curve_incidence::{CurveIncidence, certify_bspline_curve_on_support},
    analytic_radial_deviation::bound_nurbs_radial_deviation,
};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[derive(Clone, Copy, Debug)]
struct Linear {
    origin: Point3,
    axis: Vector3,
    radius: f64,
    slope: f64,
}

impl Linear {
    fn new(carrier: SurfaceCarrier) -> Option<Self> {
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: radius,
            slope,
        } = carrier
        else {
            return None;
        };
        (super::containment::unit(axis) && radius.is_finite() && slope.is_finite()).then_some(
            Self {
                origin,
                axis,
                radius,
                slope,
            },
        )
    }

    fn radial(self, delta: Vector3) -> Vector3 { delta - self.axis * delta.dot(&self.axis) }

    fn contains(self, curve: &BsplineCurve<Point3>, tolerance: f64) -> bool {
        if self.slope == 0.0 {
            return certify_bspline_curve_on_support(
                curve,
                AnalyticSupport::Cylinder {
                    origin: self.origin,
                    axis: self.axis,
                    radius: self.radius,
                },
                tolerance,
            ) == CurveIncidence::Proven;
        }
        // Recenter at a positive-radius observed point rather than form a
        // distant apex. Reverse the axis for a decreasing radius.
        let Some(first) = curve.control_points().first() else {
            return false;
        };
        let height = (*first - self.origin).dot(&self.axis);
        let radius = self.slope.mul_add(height, self.radius);
        bound_nurbs_radial_deviation(
            &NurbsCurve::from(curve.clone()),
            self.origin + self.axis * height,
            self.axis * self.slope.signum(),
            radius,
            self.slope.abs(),
            tolerance * 0.25,
        )
        .is_some_and(|bounds| bounds.radial_upper <= tolerance)
    }
}

#[derive(Clone, Copy, Debug)]
struct Frame {
    first: Linear,
    second: Linear,
    u: Vector3,
    v: Vector3,
}

#[derive(Clone, Copy)]
struct Quadratic {
    a: f64,
    b: f64,
    c: f64,
}

#[derive(Clone, Copy)]
struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    fn new(lo: f64, hi: f64) -> Self {
        let padding = 16.0 * f64::EPSILON * lo.abs().max(hi.abs());
        Self {
            lo: lo - padding,
            hi: hi + padding,
        }
    }
    fn around(x: f64, error: f64) -> Self { Self::new(x - error, x + error) }
    fn add(self, other: Self) -> Self { Self::new(self.lo + other.lo, self.hi + other.hi) }
    fn neg(self) -> Self { Self::new(-self.hi, -self.lo) }
    fn mul(self, other: Self) -> Self {
        let values = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        Self::new(
            values.into_iter().fold(f64::INFINITY, f64::min),
            values.into_iter().fold(f64::NEG_INFINITY, f64::max),
        )
    }
    fn square(self) -> Self {
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            self.lo.abs().min(self.hi.abs()).powi(2)
        };
        Self::new(lo, self.lo.abs().max(self.hi.abs()).powi(2))
    }
    fn reciprocal(self) -> Option<Self> {
        (self.lo > 0.0 || self.hi < 0.0).then(|| Self::new(self.hi.recip(), self.lo.recip()))
    }
    fn finite(self) -> bool { self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi }
}

struct Generator {
    base: Vector3,
    direction: Vector3,
    transverse: Vector3,
    transverse_direction: Vector3,
    radius: f64,
    radius_derivative: f64,
}

impl Generator {
    fn coefficients(&self) -> Quadratic {
        Quadratic {
            a: self.radius_derivative.mul_add(
                -self.radius_derivative,
                self.transverse_direction.norm_squared(),
            ),
            b: self.radius.mul_add(
                -self.radius_derivative,
                self.transverse.dot(&self.transverse_direction),
            ),
            c: self
                .radius
                .mul_add(-self.radius, self.transverse.norm_squared()),
        }
    }
}

impl Frame {
    fn new(first: Linear, second: Linear) -> Self {
        let (u, v) = normal_basis(first.axis);
        Self {
            first,
            second,
            u,
            v,
        }
    }

    fn angle(self, point: Point3) -> Option<f64> {
        let delta = point - self.first.origin;
        let x = delta.dot(&self.u);
        let y = delta.dot(&self.v);
        (x.is_finite() && y.is_finite() && x.hypot(y) > 0.0).then(|| y.atan2(x))
    }

    fn direction(self, angle: f64) -> (Vector3, Vector3) {
        let (sin, cos) = angle.sin_cos();
        (self.u * cos + self.v * sin, self.v * cos - self.u * sin)
    }

    fn generator(self, radial: Vector3) -> Generator {
        let base = (self.first.origin - self.second.origin) + radial * self.first.radius;
        let direction = self.first.axis + radial * self.first.slope;
        Generator {
            base,
            direction,
            transverse: self.second.radial(base),
            transverse_direction: self.second.radial(direction),
            radius: self
                .second
                .slope
                .mul_add(base.dot(&self.second.axis), self.second.radius),
            radius_derivative: self.second.slope * direction.dot(&self.second.axis),
        }
    }

    fn evaluate(self, positive: bool, angle: f64) -> Option<(Point3, Vector3)> {
        let (radial, tangent) = self.direction(angle);
        let g = self.generator(radial);
        let Quadratic { a, b, c } = g.coefficients();
        let disc = b.mul_add(b, -a * c);
        if !disc.is_finite() || disc <= 0.0 || a == 0.0 {
            return None;
        }
        let root = disc.sqrt() * if positive { 1.0 } else { -1.0 };
        // Stable quadratic root, avoiding subtraction of nearly equal
        // values when the actual axial coordinate is near zero.
        let z = if b * root > 0.0 {
            c / (-b - root)
        } else {
            (-b + root) / a
        };
        let radius = self.first.slope.mul_add(z, self.first.radius);
        let second_radius = g.radius_derivative.mul_add(z, g.radius);
        if !z.is_finite() || radius <= 0.0 || second_radius <= 0.0 {
            return None;
        }
        let delta = g.base + g.direction * z;
        let normal =
            self.second.radial(delta) - self.second.axis * (self.second.slope * second_radius);
        let angular = tangent * radius;
        let denominator = normal.dot(&g.direction);
        if !denominator.is_finite() || denominator == 0.0 {
            return None;
        }
        let derivative = angular - g.direction * (normal.dot(&angular) / denominator);
        let point = self.second.origin + delta;
        Some((point, derivative))
    }

    // Coefficient variation is bounded by the maximum movement of the
    // generator over the angular interval. This handles arbitrary axis
    // skew and signed affine slopes without a far-away cone apex.
    fn roots(self, start: f64, end: f64, positive: bool, tolerance: f64) -> Option<[Interval; 2]> {
        let (radial, _) = self.direction(start.midpoint(end));
        let generator = self.generator(radial);
        let Quadratic { a, b, c } = generator.coefficients();
        let angular = (end - start).abs() * 0.5;
        let move_bound = 2.0 * (angular * 0.5).sin().abs();
        let point_motion = move_bound * self.first.radius.abs();
        let direction_motion = move_bound * self.first.slope.abs();
        let radius_motion = point_motion * self.second.slope.abs();
        let slope_motion = direction_motion * self.second.slope.abs();
        let point_norm = generator.transverse.norm();
        let direction_norm = generator.transverse_direction.norm();
        let scale = generator
            .base
            .norm()
            .max(generator.direction.norm())
            .max(generator.radius.abs())
            .max(generator.radius_derivative.abs())
            .max(1.0);
        let roundoff = 256.0 * f64::EPSILON * scale * scale;
        let square_change =
            |magnitude: f64, delta: f64| delta.mul_add(delta, 2.0 * magnitude * delta);
        let ai = Interval::around(
            a,
            square_change(direction_norm, direction_motion)
                + square_change(generator.radius_derivative.abs(), slope_motion)
                + roundoff,
        );
        let product_change =
            |left: f64, right: f64, dl: f64, dr: f64| dl.mul_add(dr, left.mul_add(dr, right * dl));
        let bi = Interval::around(
            b,
            product_change(point_norm, direction_norm, point_motion, direction_motion)
                + product_change(
                    generator.radius.abs(),
                    generator.radius_derivative.abs(),
                    radius_motion,
                    slope_motion,
                )
                + roundoff,
        );
        let ci = Interval::around(
            c,
            square_change(point_norm, point_motion)
                + square_change(generator.radius.abs(), radius_motion)
                + roundoff,
        );
        let disc = bi.square().add(ai.mul(ci).neg());
        if disc.lo <= roundoff || !disc.finite() {
            return None;
        }
        let sqrt = Interval::new(disc.lo.sqrt(), disc.hi.sqrt());
        let reciprocal = ai.reciprocal()?;
        let root = |positive| {
            bi.neg()
                .add(if positive { sqrt } else { sqrt.neg() })
                .mul(reciprocal)
        };
        let selected = root(positive);
        let other = root(!positive);
        let first_radius = Interval::around(self.first.radius, roundoff)
            .add(selected.mul(Interval::around(self.first.slope, roundoff)));
        let second_radius =
            Interval::around(generator.radius, radius_motion + roundoff).add(selected.mul(
                Interval::around(generator.radius_derivative, slope_motion + roundoff),
            ));
        (selected.finite()
            && other.finite()
            && first_radius.lo > tolerance
            && second_radius.lo > tolerance
            && (2.0f64.mul_add(tolerance, selected.hi) < other.lo
                || 2.0f64.mul_add(tolerance, other.hi) < selected.lo))
            .then_some([selected, other])
    }

    fn regular(self, start: f64, end: f64, positive: bool, tolerance: f64, depth: usize) -> bool {
        if self.roots(start, end, positive, tolerance).is_some() {
            return true;
        }
        if depth == 0 {
            return false;
        }
        let mid = start.midpoint(end);
        self.regular(start, mid, positive, tolerance, depth - 1)
            && self.regular(mid, end, positive, tolerance, depth - 1)
    }

    fn cubic_on_branch(
        self,
        controls: [Point3; 4],
        start: f64,
        end: f64,
        positive: bool,
        tolerance: f64,
    ) -> bool {
        let Some([selected, other]) = self.roots(start, end, positive, tolerance) else {
            return false;
        };
        let local = controls.map(|p| p - self.first.origin);
        let heights = local.map(|p| p.dot(&self.first.axis));
        // The positive-weight polynomial control hull must not approach the
        // competing quadratic root, even if it lies on both carriers.
        let lower = heights.into_iter().fold(f64::INFINITY, f64::min);
        let upper = heights.into_iter().fold(f64::NEG_INFINITY, f64::max);
        if lower > selected.hi + tolerance
            || upper < selected.lo - tolerance
            || !(upper < other.lo - tolerance || lower > other.hi + tolerance)
        {
            return false;
        }
        let radial = local.map(|p| self.first.radial(p));
        let sign = (end - start).signum();
        // A strict common angular direction for all control vectors and
        // derivative controls proves azimuth monotonicity on the full cubic.
        let derivatives = std::array::from_fn::<_, 3, _>(|i| radial[i + 1] - radial[i]);
        let scale = radial.iter().map(|p| p.norm()).fold(1.0_f64, f64::max);
        let margin = 128.0 * f64::EPSILON * scale * scale;
        let (middle, _) = self.direction(start.midpoint(end));
        radial.iter().all(|r| {
            r.dot(&middle) > tolerance && {
                derivatives
                    .iter()
                    .all(|d| sign * self.first.axis.dot(&r.cross(&*d)) > margin)
            }
        })
    }
}

/// Source-backed regular branch of two affine-radius revolution carriers.
#[derive(Clone, Copy, Debug)]
pub struct ContactLinearSectionEvidence {
    frame: Frame,
    positive: bool,
    /// Signed unwrapped azimuth of the observed source path.
    pub signed_sweep_radians: f64,
    /// Largest source azimuth step.
    pub max_step_radians: f64,
    /// Maximum source vertex displacement onto the selected branch.
    pub max_vertex_residual: f64,
}

impl ContactLinearSectionEvidence {
    pub(super) fn angular_step(self, a: Point3, b: Point3) -> Option<f64> {
        Some((self.frame.angle(b)? - self.frame.angle(a)? + PI).rem_euclid(TAU) - PI)
    }
    pub(super) fn projected(self, point: Point3) -> Option<Point3> {
        self.frame
            .evaluate(self.positive, self.frame.angle(point)?)
            .map(|(p, _)| p)
    }
    fn measure(mut self, points: &[Point3], ids: &[usize], tolerance: f64) -> Option<Self> {
        for &id in ids {
            let p = *points.get(id)?;
            let residual = (self.projected(p)? - p).norm();
            if !residual.is_finite() || residual > tolerance {
                return None;
            }
            self.max_vertex_residual = self.max_vertex_residual.max(residual);
        }
        let mut sign = 0.0;
        for pair in ids.windows(2) {
            let a = points[pair[0]];
            let b = points[pair[1]];
            let start = self.frame.angle(a)?;
            let step = self.angular_step(a, b)?;
            if step == 0.0
                || step.abs() > FRAC_PI_2
                || step * sign < 0.0
                || !self
                    .frame
                    .regular(start, start + step, self.positive, tolerance, 8)
            {
                return None;
            }
            let (_, da) = self.frame.evaluate(self.positive, start)?;
            let (_, db) = self.frame.evaluate(self.positive, start + step)?;
            if da.dot(&(b - a)) * step <= 0.0 || db.dot(&(b - a)) * step <= 0.0 {
                return None;
            }
            sign = step.signum();
            self.signed_sweep_radians += step;
            self.max_step_radians = self.max_step_radians.max(step.abs());
        }
        (sign != 0.0
            && self.signed_sweep_radians.abs()
                <= (64.0 * f64::EPSILON).mul_add(crate::usize_to_f64(ids.len()), TAU))
        .then_some(self)
    }
    fn same_branch(self, other: Self, points: &[Point3], ids: &[usize], tolerance: f64) -> bool {
        ids.windows(2).all(|pair| {
            let a = points[pair[0]];
            let b = points[pair[1]];
            let Some(start) = self.frame.angle(a) else {
                return false;
            };
            let Some(step) = self.angular_step(a, b) else {
                return false;
            };
            [0.0, 0.5, 1.0].into_iter().all(|fraction| {
                let Some((p, dp)) = self
                    .frame
                    .evaluate(self.positive, step.mul_add(fraction, start))
                else {
                    return false;
                };
                let Some(angle) = other.frame.angle(p) else {
                    return false;
                };
                let Some((q, dq)) = other.frame.evaluate(other.positive, angle) else {
                    return false;
                };
                (p - q).norm() <= tolerance
                    && dp.dot(&dq) * self.signed_sweep_radians * other.signed_sweep_radians > 0.0
            })
        })
    }

    fn cubic(self, start: f64, end: f64, tolerance: f64) -> Option<[Point3; 4]> {
        let (a, da) = self.frame.evaluate(self.positive, start)?;
        let (b, db) = self.frame.evaluate(self.positive, end)?;
        let step = end - start;
        let controls = [a, a + da * (step / 3.0), b - db * (step / 3.0), b];
        if !self
            .frame
            .cubic_on_branch(controls, start, end, self.positive, tolerance)
        {
            return None;
        }
        let curve = BsplineCurve::new(KnotVector::bezier_knot(3), controls.to_vec()).ok()?;
        (self.frame.first.contains(&curve, tolerance)
            && self.frame.second.contains(&curve, tolerance))
        .then_some(controls)
    }

    pub(super) fn native_curve(
        self,
        a: Point3,
        b: Point3,
        sweep: f64,
        tolerance: f64,
    ) -> Option<Curve> {
        let start = self.frame.angle(a)?;
        let step = self.angular_step(a, b)?;
        let end = start + step;
        if step * sweep <= 0.0
            || step.abs() > FRAC_PI_2
            || (self.projected(a)? - a).norm() > tolerance
            || (self.projected(b)? - b).norm() > tolerance
            || !self.frame.regular(start, end, self.positive, tolerance, 8)
        {
            return None;
        }
        // Certify the actual assembled B-spline, not just the individual
        // input cubics. Reserve error budget for the knot-span extraction
        // and for the global maximum proof's conditioning allowance.
        let curve = super::section_curve::section_curve(start, end, |a, b| {
            self.cubic(a, b, tolerance * 0.5)
        })?;
        let Curve::BsplineCurve(spline) = &curve else {
            return None;
        };
        (self.frame.first.contains(spline, tolerance)
            && self.frame.second.contains(spline, tolerance))
        .then_some(curve)
    }
}

pub(super) fn supplied_section(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    tolerance: f64,
) -> Option<ContactLinearSectionEvidence> {
    let [a, b] = carriers else {
        return None;
    };
    let (a, b) = (Linear::new(*a)?, Linear::new(*b)?);
    if ids.len() < 2 || !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let mut selected: Option<ContactLinearSectionEvidence> = None;
    for frame in [Frame::new(a, b), Frame::new(b, a)] {
        for positive in [false, true] {
            let candidate = ContactLinearSectionEvidence {
                frame,
                positive,
                signed_sweep_radians: 0.0,
                max_step_radians: 0.0,
                max_vertex_residual: 0.0,
            };
            let Some(candidate) = candidate.measure(points, ids, tolerance) else {
                continue;
            };
            if let Some(previous) = selected {
                if !candidate.same_branch(previous, points, ids, tolerance)
                    || !previous.same_branch(candidate, points, ids, tolerance)
                {
                    return None;
                }
                if !candidate
                    .max_step_radians
                    .total_cmp(&previous.max_step_radians)
                    .then(
                        candidate
                            .max_vertex_residual
                            .total_cmp(&previous.max_vertex_residual),
                    )
                    .is_lt()
                {
                    continue;
                }
            }
            selected = Some(candidate);
        }
    }
    selected
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ParametricCurve;

    fn carrier(origin: Point3, axis: Vector3, radius: f64, slope: f64) -> SurfaceCarrier {
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin: radius,
            slope,
        }
    }

    #[test]
    fn unequal_cylinders_have_a_certified_nonconic_branch() -> Result<(), Box<dyn std::error::Error>>
    {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let carriers = [
            carrier(origin, Vector3::z(), 1.0, 0.0),
            carrier(origin, Vector3::x(), 2.0, 0.0),
        ];
        // Independent formula: x²+y²=1, y²+z²=4.
        let points = [0.1_f64, 0.22, 0.38].map(|angle| {
            let (y, x) = angle.sin_cos();
            Point3::new(x, y, (4.0 - y * y).sqrt())
        });
        for ids in [[0, 1, 2], [2, 1, 0]] {
            let section =
                supplied_section(&points, &ids, &carriers, 1.0e-8).ok_or("missing branch")?;
            for pair in ids.windows(2) {
                let (a, b) = (points[pair[0]], points[pair[1]]);
                let step = section.angular_step(a, b).ok_or("step")?;
                let curve = section.native_curve(a, b, step, 1.0e-8).ok_or("native")?;
                let Curve::BsplineCurve(curve) = curve else {
                    return Err("expected native spline".into());
                };
                for carrier in carriers {
                    assert!(
                        Linear::new(carrier)
                            .ok_or("carrier")?
                            .contains(&curve, 1.0e-8)
                    );
                }
                for i in 0..=20 {
                    let p = curve.evaluate(f64::from(i) / 20.0);
                    assert!(p.z > 0.0);
                    assert!((p.x.hypot(p.y) - 1.0).abs() <= 1.0e-8);
                    assert!((p.y.hypot(p.z) - 2.0).abs() <= 1.0e-8);
                }
            }
        }
        // Same azimuth, opposing roots: there is no uniquely witnessed
        // short regular arc connecting these two points.
        let wrong = [
            points[0],
            Point3::new(points[1].x, points[1].y, -points[1].z),
        ];
        assert!(supplied_section(&wrong, &[0, 1], &carriers, 1.0e-8).is_none());
        Ok(())
    }

    #[test]
    fn sparse_catalog_contact_preserves_small_signed_fitted_slopes()
    -> Result<(), Box<dyn std::error::Error>> {
        // Source endpoints and independently fitted carriers from the
        // 2091967a development model; no original CAD labels.
        let carriers = [
            carrier(
                Point3::new(3.4274859842246683, -2.35, 0.1625140157753312),
                Vector3::y(),
                0.17251401577533082,
                -1.2578156450763003e-32,
            ),
            carrier(
                Point3::new(
                    3.4301855673676367,
                    -0.15000000000000002,
                    -6.713379852030243e-16,
                ),
                Vector3::new(2.35716712806807e-16, 2.3255621459172727e-16, 1.0),
                0.2,
                4.158069842179508e-16,
            ),
        ];
        let points = [
            Point3::new(3.4301855673676367, -0.35, -0.009978876512395629),
            Point3::new(
                3.4274859842246683,
                -0.349981779797196,
                -0.009999999999999593,
            ),
        ];
        let section =
            supplied_section(&points, &[0, 1], &carriers, 1.0e-7).ok_or("catalog branch")?;
        let step = section.angular_step(points[0], points[1]).ok_or("step")?;
        let Curve::BsplineCurve(curve) = section
            .native_curve(points[0], points[1], step, 1.0e-7)
            .ok_or("catalog native edge")?
        else {
            return Err("expected spline".into());
        };
        for carrier in carriers {
            assert!(
                Linear::new(carrier)
                    .ok_or("carrier")?
                    .contains(&curve, 1.0e-7)
            );
        }
        Ok(())
    }

    #[test]
    fn skew_affine_revolutions_keep_both_signed_slopes() -> Result<(), Box<dyn std::error::Error>> {
        let carriers = [
            carrier(Point3::new(0.0, 0.0, 0.0), Vector3::z(), 1.0, 0.08),
            carrier(
                Point3::new(0.3, -0.2, 0.1),
                Vector3::new(1.0, 0.2, 0.1).normalize(),
                1.8,
                -0.05,
            ),
        ];
        let frame = Frame::new(
            Linear::new(carriers[0]).ok_or("first")?,
            Linear::new(carriers[1]).ok_or("second")?,
        );
        let points = [0.1, 0.12, 0.14].map(|angle| frame.evaluate(true, angle).map(|(p, _)| p));
        let [Some(a), Some(b), Some(c)] = points else {
            return Err("fixture".into());
        };
        let points = [a, b, c];
        let section = supplied_section(&points, &[0, 1, 2], &carriers, 1.0e-8).ok_or("branch")?;
        let step = section.angular_step(a, c).ok_or("step")?;
        let Curve::BsplineCurve(curve) =
            section.native_curve(a, c, step, 1.0e-8).ok_or("native")?
        else {
            return Err("expected spline".into());
        };
        for carrier in carriers {
            assert!(
                Linear::new(carrier)
                    .ok_or("carrier")?
                    .contains(&curve, 1.0e-8)
            );
        }
        Ok(())
    }
}
