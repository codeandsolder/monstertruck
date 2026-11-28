//! Regular plane/positive-radial toroidal intersection branches. Both angular charts
//! are available; a branch is admitted only on intervals with a resolved
//! discriminant. Polynomial native approximants carry a continuous
//! two-carrier Bernstein certificate, not merely a sampled fit.

use super::super::super::sample_fit::normal_basis;
use super::SurfaceCarrier;
use crate::{BsplineCurve, Curve, KnotVector, Point3, Vector3};
use monstertruck_geometry::{
    analytic_carrier::AnalyticSupport,
    analytic_curve_incidence::{CurveIncidence, certify_bspline_curve_on_support},
};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[derive(Clone, Copy, Debug)]
enum Chart {
    Meridian,
    Azimuth,
}

/// Source-backed regular branch of a plane/toroidal profile intersection.
///
/// This is generally not a conic. Native polynomial edges approximate the
/// selected analytic branch within an independently verified tolerance.
#[derive(Clone, Copy, Debug)]
pub struct ContactToroidalSectionEvidence {
    frame: Frame,
    chart: Chart,
    positive: bool,
    /// Signed unwrapped sweep in the selected angular chart.
    pub signed_sweep_radians: f64,
    /// Largest adjacent source step in that chart.
    pub max_step_radians: f64,
    /// Maximum source sample displacement onto the selected branch.
    pub max_vertex_residual: f64,
}

#[derive(Clone, Copy, Debug)]
struct Frame {
    plane_origin: Point3,
    normal: Vector3,
    center: Point3,
    axis: Vector3,
    major: f64,
    minor: f64,
    u: Vector3,
    v: Vector3,
    axial_normal: f64,
    transverse_normal: f64,
    offset: f64,
}

#[derive(Clone, Copy)]
struct Range {
    low: f64,
    high: f64,
}

impl Range {
    const fn point(x: f64) -> Self { Self { low: x, high: x } }

    fn from_values(values: impl IntoIterator<Item = f64>) -> Self {
        values.into_iter().fold(
            Self {
                low: f64::INFINITY,
                high: f64::NEG_INFINITY,
            },
            |range, x| Self {
                low: range.low.min(x),
                high: range.high.max(x),
            },
        )
    }

    fn mul(self, other: Self) -> Self {
        Self::from_values([
            self.low * other.low,
            self.low * other.high,
            self.high * other.low,
            self.high * other.high,
        ])
    }

    const fn add(self, other: Self) -> Self {
        Self {
            low: self.low + other.low,
            high: self.high + other.high,
        }
    }

    const fn neg(self) -> Self {
        Self {
            low: -self.high,
            high: -self.low,
        }
    }

    fn scale(self, scale: f64) -> Self { self.mul(Self::point(scale)) }

    const fn abs_max(self) -> f64 { self.low.abs().max(self.high.abs()) }

    fn strict_sign(self, positive: bool) -> bool {
        let aligned = if positive { self } else { self.neg() };
        aligned.low.is_finite()
            && aligned.high.is_finite()
            && aligned.low > 64.0 * f64::EPSILON * aligned.abs_max()
    }
}

// A range over a short unwrapped angular interval, including any interior
// extrema. Expanding by a rounding allowance keeps the branch test away
// from a numerically unresolved tangency.
fn sine_range(a: f64, b: f64) -> Range {
    let (lo, hi) = (a.min(b), a.max(b));
    let mut range = Range::from_values([lo.sin(), hi.sin()]);
    for (phase, value) in [(FRAC_PI_2, 1.0), (-FRAC_PI_2, -1.0)] {
        if (phase - lo).rem_euclid(TAU) <= 16.0_f64.mul_add(f64::EPSILON, hi - lo) {
            range.low = range.low.min(value);
            range.high = range.high.max(value);
        }
    }
    Range {
        low: 16.0f64.mul_add(-f64::EPSILON, range.low).max(-1.0),
        high: 16.0f64.mul_add(f64::EPSILON, range.high).min(1.0),
    }
}

impl Frame {
    fn new(carriers: &[SurfaceCarrier]) -> Option<Self> {
        let (plane, torus) = match carriers {
            [
                plane @ SurfaceCarrier::Plane { .. },
                torus @ SurfaceCarrier::CircularRevolution { .. },
            ]
            | [
                torus @ SurfaceCarrier::CircularRevolution { .. },
                plane @ SurfaceCarrier::Plane { .. },
            ] => (*plane, *torus),
            _ => return None,
        };
        let SurfaceCarrier::Plane {
            origin: plane_origin,
            normal,
        } = plane
        else {
            return None;
        };
        let SurfaceCarrier::CircularRevolution {
            axial_center: center,
            axis,
            radial_center: major,
            profile_radius: minor,
        } = torus
        else {
            return None;
        };
        if !super::containment::unit(axis)
            || !super::containment::unit(normal)
            || !major.is_finite()
            || !minor.is_finite()
            || major <= 0.0
            || minor <= 0.0
        {
            return None;
        }
        let axial_normal = normal.dot(&axis);
        let transverse = normal - axis * axial_normal;
        let transverse_normal = transverse.norm();
        let (u, v) = if transverse_normal > 64.0 * f64::EPSILON {
            let u = transverse / transverse_normal;
            (u, axis.cross(&u).normalize())
        } else {
            normal_basis(axis)
        };
        let offset = (center - plane_origin).dot(&normal);
        Some(Self {
            plane_origin,
            normal,
            center,
            axis,
            major,
            minor,
            u,
            v,
            axial_normal,
            transverse_normal,
            offset,
        })
    }

    const fn analytic_supports(self) -> [AnalyticSupport; 2] {
        [
            AnalyticSupport::Plane {
                origin: self.plane_origin,
                normal: self.normal,
            },
            AnalyticSupport::Torus {
                center: self.center,
                axis: self.axis,
                major_radius: self.major,
                minor_radius: self.minor,
            },
        ]
    }

    fn radial(self, point: Point3) -> (f64, Vector3) {
        let delta = point - self.center;
        let z = delta.dot(&self.axis);
        (z, delta - self.axis * z)
    }

    fn angle(self, chart: Chart, point: Point3) -> Option<f64> {
        let (z, radial) = self.radial(point);
        let (x, y) = match chart {
            Chart::Meridian => (radial.norm() - self.major, z),
            Chart::Azimuth => (radial.dot(&self.u), radial.dot(&self.v)),
        };
        (x.is_finite() && y.is_finite() && x.hypot(y) > 0.0).then(|| y.atan2(x))
    }

    fn evaluate(self, chart: Chart, positive: bool, t: f64) -> Option<(Point3, Vector3)> {
        match chart {
            Chart::Meridian => self.meridian(positive, t),
            Chart::Azimuth => self.azimuth(positive, t),
        }
    }

    // Fix the meridian angle. The plane fixes one transverse coordinate;
    // intersect that line with the positive-radius parallel circle.
    fn meridian(self, positive: bool, angle: f64) -> Option<(Point3, Vector3)> {
        if self.transverse_normal <= 64.0 * f64::EPSILON {
            return None;
        }
        let (sin, cos) = angle.sin_cos();
        let axial = self.minor * sin;
        let profile_offset = self.minor * cos;
        let rho = self.major + profile_offset;
        let plane_x = -self.axial_normal.mul_add(axial, self.offset) / self.transverse_normal;
        let squared = (rho - plane_x) * (rho + plane_x);
        if rho <= 0.0 || squared <= 0.0 || !squared.is_finite() {
            return None;
        }
        let plane_y = squared.sqrt() * if positive { 1.0 } else { -1.0 };
        let x_derivative = -self.axial_normal * profile_offset / self.transverse_normal;
        let y_derivative = rho.mul_add(-axial, -plane_x * x_derivative) / plane_y;
        Some((
            self.center + self.axis * axial + self.u * plane_x + self.v * plane_y,
            self.axis * profile_offset + self.u * x_derivative + self.v * y_derivative,
        ))
    }

    // Fix azimuth. The plane restricts the meridian circle by a line.
    // The sign of a*profile_offset-b*axial selects one of its two intersections.
    fn azimuth(self, positive: bool, angle: f64) -> Option<(Point3, Vector3)> {
        let (sin, cos) = angle.sin_cos();
        let axial_coeff = self.axial_normal;
        let radial_coeff = self.transverse_normal * cos;
        let normal_squared = axial_coeff.mul_add(axial_coeff, radial_coeff * radial_coeff);
        let constant = -radial_coeff.mul_add(self.major, self.offset);
        let discriminant = self
            .minor
            .mul_add(self.minor * normal_squared, -constant * constant);
        if normal_squared <= 0.0 || discriminant <= 0.0 || !discriminant.is_finite() {
            return None;
        }
        let determinant = discriminant.sqrt() * if positive { 1.0 } else { -1.0 };
        let axial = axial_coeff.mul_add(constant, -radial_coeff * determinant) / normal_squared;
        let profile_offset =
            radial_coeff.mul_add(constant, axial_coeff * determinant) / normal_squared;
        let rho = self.major + profile_offset;
        if rho <= 0.0 {
            return None;
        }
        let radial_coeff_derivative = -self.transverse_normal * sin;
        let axial_derivative = -rho * radial_coeff_derivative * profile_offset / determinant;
        let radial_derivative = rho * radial_coeff_derivative * axial / determinant;
        let direction = self.u * cos + self.v * sin;
        let tangent = self.v * cos - self.u * sin;
        Some((
            self.center + self.axis * axial + direction * rho,
            self.axis * axial_derivative + direction * radial_derivative + tangent * rho,
        ))
    }

    // The azimuth discriminant is a quadratic in b. Ring tori give a
    // concave quadratic; spindle tori give a convex one. Both extrema must
    // be included, since the selected branch radius also depends on sqrt(D).
    fn azimuth_discriminant(self, b: Range) -> Range {
        let evaluate = |value: f64| {
            let normal_squared = self.axial_normal.mul_add(self.axial_normal, value * value);
            let offset = value.mul_add(self.major, self.offset);
            self.minor
                .mul_add(self.minor * normal_squared, -offset * offset)
        };
        let mut result = Range::from_values([evaluate(b.low), evaluate(b.high)]);
        let coefficient = (self.minor - self.major) * (self.minor + self.major);
        if coefficient != 0.0 {
            let critical = self.major * self.offset / coefficient;
            if (b.low..=b.high).contains(&critical) {
                let d = evaluate(critical);
                result.low = result.low.min(d);
                result.high = result.high.max(d);
            }
        }
        result
    }

    fn regular(
        self,
        chart: Chart,
        positive: bool,
        start: f64,
        end: f64,
        tolerance: f64,
        depth: usize,
    ) -> bool {
        let sine = sine_range(start, end);
        let cosine = sine_range(start + FRAC_PI_2, end + FRAC_PI_2);
        let scale = self.major.max(self.minor).max(self.offset.abs());
        let roundoff = 256.0 * f64::EPSILON * scale * scale;
        let margin = tolerance.mul_add(tolerance, roundoff);
        let admitted = match chart {
            Chart::Meridian => {
                if self.transverse_normal <= 64.0 * f64::EPSILON {
                    return false;
                }
                let rho = self.minor.mul_add(cosine.low, self.major);
                let x = sine
                    .scale(self.minor * self.axial_normal)
                    .add(Range::point(self.offset))
                    .abs_max()
                    / self.transverse_normal;
                rho > tolerance && (rho - x) * (rho + x) > margin
            }
            Chart::Azimuth => {
                let b = cosine.scale(self.transverse_normal);
                let d = self.azimuth_discriminant(b);
                if d.low <= margin || !d.high.is_finite() {
                    false
                } else if self.major - self.minor > tolerance {
                    true
                } else {
                    // rho = (R*a^2 - b*offset + a*det)/(a^2+b^2).
                    // Use an interval for the complete selected root, not
                    // just its positive endpoint radii. Neither a crossing
                    // of the axis nor the reflected spindle sheet is valid.
                    let a = self.axial_normal;
                    let determinant = Range {
                        low: (d.low - roundoff).sqrt(),
                        high: (d.high + roundoff).sqrt(),
                    }
                    .scale(if positive { 1.0 } else { -1.0 });
                    let numerator = Range::point(self.major * a * a)
                        .add(b.scale(-self.offset))
                        .add(determinant.scale(a));
                    let denominator = a.mul_add(a, b.abs_max().powi(2));
                    let rho_margin = 256.0 * f64::EPSILON * scale;
                    numerator.low > tolerance.mul_add(denominator, rho_margin)
                }
            }
        };
        if admitted {
            return true;
        }
        if depth == 0 {
            return false;
        }
        let mid = start.midpoint(end);
        self.regular(chart, positive, start, mid, tolerance, depth - 1)
            && self.regular(chart, positive, mid, end, tolerance, depth - 1)
    }

    // Full control-hull proof that the polynomial stays on the selected
    // branch and progresses in the observed angular direction.
    fn chart_monotone(
        self,
        chart: Chart,
        positive: bool,
        forward: bool,
        controls: [Point3; 4],
    ) -> bool {
        let local = controls.map(|p| self.radial(p));
        let derivative = std::array::from_fn::<_, 3, _>(|i| {
            let (z0, r0) = local[i];
            let (z1, r1) = local[i + 1];
            (3.0 * (z1 - z0), (r1 - r0) * 3.0)
        });
        let Some((_, radial)) = local
            .iter()
            .copied()
            .max_by(|a, b| a.1.norm_squared().total_cmp(&b.1.norm_squared()))
        else {
            return false;
        };
        if radial.norm_squared() <= 0.0 {
            return false;
        }
        let direction = radial.normalize();
        let rho = Range {
            low: local
                .iter()
                .map(|(_, r)| r.dot(&direction))
                .fold(f64::INFINITY, f64::min),
            high: local.iter().map(|(_, r)| r.norm()).fold(0.0_f64, f64::max),
        };
        if rho.low <= 0.0 {
            return false;
        }
        let q = rho.add(Range::point(-self.major));
        let z = Range::from_values(local.iter().map(|(z, _)| *z));
        match chart {
            Chart::Meridian => {
                let branch = Range::from_values(local.iter().map(|(_, r)| r.dot(&self.v)));
                let dz = Range::from_values(derivative.iter().map(|(z, _)| *z));
                let radial_dot = Range::from_values(
                    local
                        .iter()
                        .flat_map(|(_, r)| derivative.iter().map(move |(_, dr)| r.dot(dr))),
                );
                let progress = rho.mul(q).mul(dz).add(z.mul(radial_dot).neg());
                branch.strict_sign(positive) && progress.strict_sign(forward)
            }
            Chart::Azimuth => {
                let projected = Range::from_values(local.iter().map(|(_, r)| r.dot(&self.normal)));
                let branch = rho
                    .mul(q)
                    .scale(self.axial_normal)
                    .add(z.mul(projected).neg());
                let progress = Range::from_values(local.iter().flat_map(|(_, r)| {
                    derivative
                        .iter()
                        .map(move |(_, dr)| self.axis.dot(&r.cross(dr)))
                }));
                branch.strict_sign(positive) && progress.strict_sign(forward)
            }
        }
    }
}

impl ContactToroidalSectionEvidence {
    pub(super) fn angular_step(self, a: Point3, b: Point3) -> Option<f64> {
        let difference = self.frame.angle(self.chart, b)? - self.frame.angle(self.chart, a)?;
        Some((difference + PI).rem_euclid(TAU) - PI)
    }

    fn evaluated(self, t: f64) -> Option<(Point3, Vector3)> {
        self.frame.evaluate(self.chart, self.positive, t)
    }

    fn project(self, p: Point3) -> Option<(Point3, Vector3)> {
        self.evaluated(self.frame.angle(self.chart, p)?)
    }

    // The regular analytic branch supplies an independent candidate for
    // a shared native junction. The caller must still check all incident
    // constraints and the original source displacement budget.
    pub(super) fn projected(self, point: Point3) -> Option<Point3> {
        self.project(point).map(|(p, _)| p)
    }

    fn measure(mut self, points: &[Point3], ids: &[usize], tolerance: f64) -> Option<Self> {
        for &id in ids {
            let point = *points.get(id)?;
            let (projected, _) = self.project(point)?;
            let residual = (projected - point).norm();
            if !residual.is_finite() || residual > tolerance {
                return None;
            }
            self.max_vertex_residual = self.max_vertex_residual.max(residual);
        }
        let mut sign = 0.0;
        for pair in ids.windows(2) {
            let a = points[pair[0]];
            let b = points[pair[1]];
            let start = self.frame.angle(self.chart, a)?;
            let step = self.angular_step(a, b)?;
            if step == 0.0
                || step.abs() > FRAC_PI_2
                || step * sign < 0.0
                || !self
                    .frame
                    .regular(self.chart, self.positive, start, start + step, tolerance, 8)
            {
                return None;
            }
            let (_, first_tangent) = self.evaluated(start)?;
            let (_, last_tangent) = self.evaluated(start + step)?;
            let chord = b - a;
            if first_tangent.dot(&chord) * step <= 0.0 || last_tangent.dot(&chord) * step <= 0.0 {
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
        ids.iter().all(|&id| {
            let (Some((a, da)), Some((b, db))) =
                (self.project(points[id]), other.project(points[id]))
            else {
                return false;
            };
            (a - b).norm() <= tolerance
                && da.dot(&db) * self.signed_sweep_radians * other.signed_sweep_radians > 0.0
        })
    }

    fn cubic(self, start: f64, end: f64, tolerance: f64) -> Option<[Point3; 4]> {
        let (a, da) = self.evaluated(start)?;
        let (b, db) = self.evaluated(end)?;
        let step = end - start;
        let controls = [a, a + da * (step / 3.0), b - db * (step / 3.0), b];
        if !self
            .frame
            .chart_monotone(self.chart, self.positive, step > 0.0, controls)
        {
            return None;
        }
        let curve = BsplineCurve::new(KnotVector::bezier_knot(3), controls.to_vec()).ok()?;
        self.frame
            .analytic_supports()
            .into_iter()
            .all(|support| {
                certify_bspline_curve_on_support(&curve, support, tolerance)
                    == CurveIncidence::Proven
            })
            .then_some(controls)
    }

    pub(super) fn native_curve(
        self,
        a: Point3,
        b: Point3,
        sweep: f64,
        tolerance: f64,
    ) -> Option<Curve> {
        let start = self.frame.angle(self.chart, a)?;
        let step = self.angular_step(a, b)?;
        let end = start + step;
        if step * sweep <= 0.0
            || step.abs() > FRAC_PI_2
            || (self.evaluated(start)?.0 - a).norm() > tolerance
            || (self.evaluated(end)?.0 - b).norm() > tolerance
            || !self
                .frame
                .regular(self.chart, self.positive, start, end, tolerance, 8)
        {
            return None;
        }
        let curve = super::section_curve::section_curve(start, end, |a, b| {
            self.cubic(a, b, tolerance * 0.5)
        })?;
        let Curve::BsplineCurve(spline) = &curve else {
            return None;
        };
        self.frame
            .analytic_supports()
            .into_iter()
            .all(|support| {
                certify_bspline_curve_on_support(spline, support, tolerance)
                    == CurveIncidence::Proven
            })
            .then_some(curve)
    }
}

pub(super) fn supplied_section(
    points: &[Point3],
    ids: &[usize],
    carriers: &[SurfaceCarrier],
    tolerance: f64,
) -> Option<ContactToroidalSectionEvidence> {
    let frame = Frame::new(carriers)?;
    if ids.len() < 2 || !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let mut selected: Option<ContactToroidalSectionEvidence> = None;
    for chart in [Chart::Meridian, Chart::Azimuth] {
        for positive in [false, true] {
            let candidate = ContactToroidalSectionEvidence {
                frame,
                chart,
                positive,
                signed_sweep_radians: 0.0,
                max_step_radians: 0.0,
                max_vertex_residual: 0.0,
            };
            let Some(candidate) = candidate.measure(points, ids, tolerance) else {
                continue;
            };
            if let Some(previous) = selected {
                if !candidate.same_branch(previous, points, ids, tolerance) {
                    return None;
                }
                if !candidate
                    .max_vertex_residual
                    .total_cmp(&previous.max_vertex_residual)
                    .then(
                        candidate
                            .max_step_radians
                            .total_cmp(&previous.max_step_radians),
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

    fn carriers(normal: Vector3, origin: Point3) -> [SurfaceCarrier; 2] {
        [
            SurfaceCarrier::Plane {
                origin,
                normal: normal.normalize(),
            },
            SurfaceCarrier::CircularRevolution {
                axial_center: Point3::new(0.0, 0.0, 0.0),
                axis: Vector3::z(),
                radial_center: 2.0,
                profile_radius: 0.5,
            },
        ]
    }

    fn check_section(
        carriers: [SurfaceCarrier; 2],
        points: &[Point3],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let ids = (0..points.len()).collect::<Vec<_>>();
        let section = supplied_section(points, &ids, &carriers, 1.0e-8).ok_or("missing section")?;
        for pair in points.windows(2) {
            let step = section.angular_step(pair[0], pair[1]).ok_or("step")?;
            let curve = section
                .native_curve(pair[0], pair[1], step, 1.0e-8)
                .ok_or("native curve")?;
            let Curve::BsplineCurve(curve) = curve else {
                return Err("expected spline".into());
            };
            for support in section.frame.analytic_supports() {
                assert_eq!(
                    certify_bspline_curve_on_support(&curve, support, 1.0e-8),
                    CurveIncidence::Proven
                );
            }
            for i in 0..=20 {
                let p = curve.evaluate(f64::from(i) / 20.0);
                for carrier in carriers {
                    assert!(carrier.signed_constraint(p).ok_or("constraint")?.1.abs() < 1.0e-8);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn vertical_plane_section_crosses_an_azimuth_extremum() -> Result<(), Box<dyn std::error::Error>>
    {
        let carriers = carriers(Vector3::x(), Point3::new(0.7, 0.0, 0.0));
        let points = [-0.6_f64, -0.3, 0.0, 0.3, 0.6].map(|t| {
            let rho = 0.5_f64.mul_add(t.cos(), 2.0);
            Point3::new(0.7, rho.mul_add(rho, -0.49).sqrt(), 0.5 * t.sin())
        });
        check_section(carriers, &points)?;
        check_section(
            [carriers[1], carriers[0]],
            &points.into_iter().rev().collect::<Vec<_>>(),
        )?;
        Ok(())
    }

    #[test]
    fn spindle_and_horn_sections_retain_the_positive_radial_chart()
    -> Result<(), Box<dyn std::error::Error>> {
        for major in [0.7, 1.0] {
            let surfaces = [
                SurfaceCarrier::Plane {
                    origin: Point3::new(0.4, 0.0, 0.0),
                    normal: Vector3::x(),
                },
                SurfaceCarrier::CircularRevolution {
                    axial_center: Point3::new(0.0, 0.0, 0.0),
                    axis: Vector3::z(),
                    radial_center: major,
                    profile_radius: 1.0,
                },
            ];
            let points = [-0.4_f64, -0.2, 0.0, 0.2, 0.4].map(|angle| {
                let rho = major + angle.cos();
                Point3::new(0.4, (rho * rho - 0.16).sqrt(), angle.sin())
            });
            check_section(surfaces, &points)?;
            // The regular local section has no continuation through the
            // axis into the reflected spindle sheet.
            let frame = Frame::new(&surfaces).ok_or("frame")?;
            assert!(!frame.regular(Chart::Meridian, true, PI - 0.05, PI + 0.05, 1.0e-8, 8));
        }
        Ok(())
    }

    #[test]
    fn oblique_plane_both_azimuth_branches() -> Result<(), Box<dyn std::error::Error>> {
        let carriers = carriers(Vector3::new(-0.1, 0.0, 1.0), Point3::new(0.0, 0.0, 0.05));
        let frame = Frame::new(&carriers).ok_or("frame")?;
        for positive in [false, true] {
            let points = [0.2, 0.4, 0.6, 0.8]
                .into_iter()
                .map(|angle| frame.azimuth(positive, angle).map(|p| p.0).ok_or("azimuth"))
                .collect::<Result<Vec<_>, _>>()?;
            check_section(carriers, &points)?;
            let mut wrong = points.clone();
            wrong[1].z += 0.01;
            assert!(supplied_section(&wrong, &[0, 1, 2, 3], &carriers, 1.0e-8).is_none());
        }
        Ok(())
    }
}
