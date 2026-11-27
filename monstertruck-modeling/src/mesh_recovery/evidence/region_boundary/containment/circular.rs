//! Continuous restrictions of a circular meridian to a line or circle.
//! Squaring the meridian equation produces a quartic whose circle restriction
//! is only degree two in the angular Fourier basis. The unsquared branch is
//! checked separately; a small quartic residual alone is insufficient.

use super::{CarrierResidual, residual, unit};
use crate::{InnerSpace, Point3, Vector3};
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy)]
pub(super) struct CircularProfile {
    center: Point3,
    axis: Vector3,
    radial_center: f64,
    profile_radius: f64,
}

impl CircularProfile {
    pub(super) fn new(
        center: Point3,
        axis: Vector3,
        radial_center: f64,
        profile_radius: f64,
    ) -> Option<Self> {
        (unit(axis)
            && [center.x, center.y, center.z, radial_center, profile_radius]
                .into_iter()
                .all(f64::is_finite)
            && profile_radius > 0.0)
            .then_some(Self {
                center,
                axis,
                radial_center,
                profile_radius,
            })
    }

    fn transverse(self, delta: Vector3) -> Vector3 { delta - self.axis * delta.dot(self.axis) }

    fn error(self, delta: Vector3) -> f64 {
        let axial = delta.dot(self.axis);
        let radial = self.transverse(delta).magnitude();
        ((radial - self.radial_center).hypot(axial) - self.profile_radius).abs()
    }

    fn offset(self) -> f64 {
        self.radial_center
            .mul_add(self.radial_center, -self.profile_radius.powi(2))
    }

    fn bound(self, a_range: [f64; 2], squared_bound: f64, witness: f64) -> Option<CarrierResidual> {
        if !a_range
            .into_iter()
            .chain([squared_bound, witness])
            .all(f64::is_finite)
        {
            return None;
        }
        // The desired residual before the last square root is
        // A - 2*R*rho, where A = |delta|² + R² - r².
        // The polynomial Q = A² - 4*R²*rho² contains both radial branches.
        // Sign-align A with R to select the observed meridian branch.
        let meridian_squared_error = if self.radial_center == 0.0 {
            a_range[0].abs().max(a_range[1].abs())
        } else {
            let aligned_min = if self.radial_center > 0.0 {
                a_range[0]
            } else {
                -a_range[1]
            };
            let root = squared_bound.sqrt();
            if aligned_min > 0.0 {
                root.min(squared_bound / aligned_min)
            } else {
                2.0f64.mul_add(-aligned_min, root)
            }
        };
        // |distance-r| <= min(sqrt(|distance²-r²|), |distance²-r²|/r).
        let length_error = meridian_squared_error
            .sqrt()
            .min(meridian_squared_error / self.profile_radius);
        residual(length_error, witness, None)
    }

    pub(super) fn line(self, start: Point3, end: Point3) -> Option<CarrierResidual> {
        let origin = start - self.center;
        let direction = end - start;
        let a = [
            origin.magnitude2() + self.offset(),
            2.0 * origin.dot(direction),
            direction.magnitude2(),
        ];
        let transverse_origin = self.transverse(origin);
        let transverse_direction = self.transverse(direction);
        let radial_squared = [
            transverse_origin.magnitude2(),
            2.0 * transverse_origin.dot(transverse_direction),
            transverse_direction.magnitude2(),
        ];
        let mut squared = [0.0_f64; 5];
        for (i, &left) in a.iter().enumerate() {
            for (j, &right) in a.iter().enumerate() {
                squared[i + j] = left.mul_add(right, squared[i + j]);
            }
        }
        let factor = -4.0 * self.radial_center.powi(2);
        for (coefficient, radial) in squared.iter_mut().zip(radial_squared) {
            *coefficient = factor.mul_add(radial, *coefficient);
        }
        let witness = (0..=4)
            .map(|i| self.error(origin + direction * (f64::from(i) * 0.25)))
            .fold(0.0_f64, f64::max);
        let mut measured = self.bound(
            quadratic_range(a)?,
            bernstein_quartic_bound(squared)?,
            witness,
        )?;
        // Cylindrical-coordinate projection (axial, radial distance) is
        // nonexpansive, and distance to the meridian circle is 1-Lipschitz.
        // Each point on this finite segment lies within length/8 of one
        // of the five evaluated witnesses. Intersect the two independent
        // upper bounds, especially when the squared branch bound is loose.
        let lipschitz = direction.magnitude().mul_add(0.125, witness);
        measured.upper_bound = measured.upper_bound.min(lipschitz).max(witness);
        Some(measured)
    }

    // Bound the meridian distance over the convex hull of the source
    // samples, not just the perimeter. For x = sum(w_i x_i), rho(x)
    // lies between a fixed transverse projection and sum(w_i rho_i).
    // Both surrogates are affine in the weights. Convexity of the
    // meridian norm gives an upper bound from the two sets of vertex
    // values; one supporting meridian direction gives a lower bound.
    fn hull_bound(self, points: &[Vector3]) -> Option<CarrierResidual> {
        let first = *points.first()?;
        let center = first
            + points
                .iter()
                .map(|&p| p - first)
                .fold(Vector3::new(0.0, 0.0, 0.0), |sum, p| sum + p)
                / crate::usize_to_f64(points.len());
        let radial = self.transverse(center);
        let radial_direction = if radial.magnitude2() > 0.0 {
            radial.normalize()
        } else {
            let farthest = points
                .iter()
                .map(|&p| self.transverse(p))
                .max_by(|a, b| a.magnitude2().total_cmp(&b.magnitude2()))?;
            if farthest.magnitude2() > 0.0 {
                farthest.normalize()
            } else {
                super::normal_basis(self.axis).0
            }
        };
        let axial = center.dot(self.axis);
        let offset = radial.magnitude() - self.radial_center;
        let distance = axial.hypot(offset);
        let (direction_z, direction_r) = if distance > 0.0 {
            (axial / distance, offset / distance)
        } else {
            (0.0, 1.0)
        };
        let mut upper = 0.0_f64;
        let mut lower = f64::INFINITY;
        let mut witnessed = self.error(center);
        let mut scale = self.radial_center.abs().max(self.profile_radius);
        for &point in points {
            let z = point.dot(self.axis);
            let radial = self.transverse(point);
            let low = radial.dot(radial_direction) - self.radial_center;
            let high = radial.magnitude() - self.radial_center;
            upper = upper.max(z.hypot(low)).max(z.hypot(high));
            lower = lower.min(
                direction_z.mul_add(z, direction_r * if direction_r >= 0.0 { low } else { high }),
            );
            witnessed = witnessed.max(self.error(point));
            scale = scale.max(point.magnitude());
        }
        // The supporting affine expression may be negative, but the
        // actual meridian distance is always nonnegative.
        let roundoff = 64.0 * f64::EPSILON * scale;
        let bound = (upper - self.profile_radius)
            .max(self.profile_radius - lower.max(0.0))
            .max(witnessed)
            + roundoff;
        residual(bound, witnessed, None)
    }

    fn triangle_to_budget(
        self,
        vertices: [Vector3; 3],
        budget: f64,
        depth: usize,
    ) -> Option<CarrierResidual> {
        let measured = self.hull_bound(&vertices)?;
        if depth == 0 || measured.upper_bound <= budget || measured.witnessed_max > budget {
            return Some(measured);
        }
        // Split the longest edge, so even a skinny source triangle
        // converges without wasting subdivision on its short edges.
        let (a, b, c) =
            [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
                .into_iter()
                .max_by(|&(a, b, _), &(i, j, _)| {
                    (vertices[b] - vertices[a])
                        .magnitude2()
                        .total_cmp(&(vertices[j] - vertices[i]).magnitude2())
                })?;
        let middle = vertices[a] + (vertices[b] - vertices[a]) * 0.5;
        let left =
            self.triangle_to_budget([vertices[a], middle, vertices[c]], budget, depth - 1)?;
        if left.witnessed_max > budget {
            return residual(
                measured.upper_bound,
                measured.witnessed_max.max(left.witnessed_max),
                None,
            );
        }
        let right =
            self.triangle_to_budget([middle, vertices[b], vertices[c]], budget, depth - 1)?;
        residual(
            left.upper_bound.max(right.upper_bound),
            left.witnessed_max.max(right.witnessed_max),
            None,
        )
    }

    // Source polygon interiors use the same vertex-zero fan semantics as
    // the mesh triangle iterator. First certify the entire convex hull
    // when possible; only subdivide triangles for inconclusive bounds.
    pub(super) fn facet(self, points: &[Point3], budget: f64) -> Option<CarrierResidual> {
        if points.len() < 3 || budget.is_nan() || budget < 0.0 {
            return None;
        }
        let local = points.iter().map(|&p| p - self.center).collect::<Vec<_>>();
        let measured = self.hull_bound(&local)?;
        if measured.upper_bound <= budget {
            let witness = local.iter().map(|&p| self.error(p)).fold(0.0_f64, f64::max);
            return residual(measured.upper_bound, witness, None);
        }
        // The hull centroid need not lie in a concave source polygon.
        // Do not use its witness to reject the original filled fan.
        let mut upper = 0.0_f64;
        let mut witnessed = 0.0_f64;
        for pair in local[1..].windows(2) {
            let measured = self.triangle_to_budget([local[0], pair[0], pair[1]], budget, 8)?;
            upper = upper.max(measured.upper_bound);
            witnessed = witnessed.max(measured.witnessed_max);
        }
        residual(upper, witnessed, None)
    }

    pub(super) fn circle(
        self,
        center: Point3,
        basis_u: Vector3,
        basis_v: Vector3,
        radius: f64,
    ) -> Option<CarrierResidual> {
        self.circle_over(center, basis_u, basis_v, radius, std::f64::consts::TAU)
    }

    // A spindle meridian's supporting circle crosses both signed
    // radial sheets. Only the observed finite arc needs to lie on
    // the selected sheet. Restrict the unsquared branch test and
    // both harmonic ranges to the actual signed angular interval.
    pub(super) fn circle_over(
        self,
        center: Point3,
        basis_u: Vector3,
        basis_v: Vector3,
        radius: f64,
        sweep: f64,
    ) -> Option<CarrierResidual> {
        if !sweep.is_finite() || sweep == 0.0 {
            return None;
        }
        let origin = center - self.center;
        // For orthogonal equal-length circle basis vectors, |delta|² has
        // only constant and first harmonics, even on a tilted circle.
        let constant = radius.mul_add(radius, origin.magnitude2()) + self.offset();
        let cosine = 2.0 * origin.dot(basis_u);
        let sine = 2.0 * origin.dot(basis_v);
        let a_squared = Harmonics::square_first(constant, cosine, sine);
        let radial_squared = Harmonics::norm_squared(
            self.transverse(origin),
            self.transverse(basis_u),
            self.transverse(basis_v),
        );
        let polynomial = a_squared.add_scaled(radial_squared, -4.0 * self.radial_center.powi(2));
        let witness = (0..=8)
            .map(|i| {
                let (sin, cos) = (f64::from(i) * 0.125 * sweep).sin_cos();
                self.error(origin + basis_u * cos + basis_v * sin)
            })
            .fold(0.0_f64, f64::max);
        let [lo, hi] = harmonic_range(cosine, sine, sweep);
        self.bound(
            [constant + lo, constant + hi],
            polynomial.arc_bound(sweep),
            witness,
        )
    }
}

/// Exact quadratic range over the finite parameter interval [0,1].
fn quadratic_range([constant, linear, quadratic]: [f64; 3]) -> Option<[f64; 2]> {
    if ![constant, linear, quadratic]
        .into_iter()
        .all(f64::is_finite)
    {
        return None;
    }
    let end = constant + linear + quadratic;
    let mut range = [constant.min(end), constant.max(end)];
    if quadratic != 0.0 {
        let stationary = -linear / (2.0 * quadratic);
        if (0.0..=1.0).contains(&stationary) {
            let value = quadratic
                .mul_add(stationary, linear)
                .mul_add(stationary, constant);
            range = [range[0].min(value), range[1].max(value)];
        }
    }
    Some(range)
}

/// Convex-hull property of the Bernstein representation over [0,1].
/// Rows map ascending power coefficients to degree-four Bernstein weights.
fn bernstein_quartic_bound(power: [f64; 5]) -> Option<f64> {
    const TRANSFORM: [[f64; 5]; 5] = [
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 0.25, 0.0, 0.0, 0.0],
        [1.0, 0.5, 1.0 / 6.0, 0.0, 0.0],
        [1.0, 0.75, 0.5, 0.25, 0.0],
        [1.0, 1.0, 1.0, 1.0, 1.0],
    ];
    if !power.into_iter().all(f64::is_finite) {
        return None;
    }
    Some(
        TRANSFORM
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .zip(power)
                    .fold(0.0_f64, |sum, (weight, coefficient)| {
                        weight.mul_add(coefficient, sum)
                    })
                    .abs()
            })
            .fold(0.0_f64, f64::max),
    )
}

// Exact extrema of c*cos(t)+s*sin(t) on the unwrapped interval
// between zero and sweep, including stationary points in its interior.
fn harmonic_range(cosine: f64, sine: f64, sweep: f64) -> [f64; 2] {
    let amplitude = cosine.hypot(sine);
    if sweep.abs() >= TAU {
        return [-amplitude, amplitude];
    }
    let (lo, hi) = (sweep.min(0.0), sweep.max(0.0));
    let phase = sine.atan2(cosine);
    let (sin, cos) = sweep.sin_cos();
    let end = cosine.mul_add(cos, sine * sin);
    let mut lower = cosine.min(end);
    let mut upper = cosine.max(end);
    let allowance = 32.0 * f64::EPSILON;
    if (phase - lo).rem_euclid(TAU) <= hi - lo + allowance {
        upper = amplitude;
    }
    if (phase + PI - lo).rem_euclid(TAU) <= hi - lo + allowance {
        lower = -amplitude;
    }
    let roundoff = allowance * amplitude;
    [lower - roundoff, upper + roundoff]
}

#[derive(Clone, Copy)]
struct Harmonics([f64; 5]);

impl Harmonics {
    fn square_first(constant: f64, cosine: f64, sine: f64) -> Self {
        Self([
            constant.mul_add(constant, cosine.powi(2).midpoint(sine.powi(2))),
            2.0 * constant * cosine,
            2.0 * constant * sine,
            0.5 * cosine.mul_add(cosine, -sine.powi(2)),
            cosine * sine,
        ])
    }

    fn norm_squared(origin: Vector3, cosine: Vector3, sine: Vector3) -> Self {
        let cos_squared = cosine.magnitude2();
        let sin_squared = sine.magnitude2();
        Self([
            origin.magnitude2() + cos_squared.midpoint(sin_squared),
            2.0 * origin.dot(cosine),
            2.0 * origin.dot(sine),
            0.5 * (cos_squared - sin_squared),
            cosine.dot(sine),
        ])
    }

    fn add_scaled(self, other: Self, scale: f64) -> Self {
        Self(std::array::from_fn(|i| {
            scale.mul_add(other.0[i], self.0[i])
        }))
    }

    fn arc_bound(self, sweep: f64) -> f64 {
        let [constant, cosine, sine, double_cosine, double_sine] = self.0;
        let [lo, hi] = harmonic_range(cosine, sine, sweep);
        let [lo2, hi2] = harmonic_range(double_cosine, double_sine, 2.0 * sweep);
        (constant + lo + lo2).abs().max((constant + hi + hi2).abs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::sample_fit::normal_basis;
    use std::f64::consts::TAU;

    #[test]
    fn toroidal_inner_ring_perimeter_does_not_certify_filled_disk()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(0.0, 0.0, 0.0);
        let profile = CircularProfile::new(center, Vector3::unit_z(), 2.0, 1.0).ok_or("profile")?;
        let points = (0..64)
            .map(|i| {
                let (sin, cos) = (TAU * f64::from(i) / 64.0).sin_cos();
                Point3::new(cos, sin, 0.0)
            })
            .collect::<Vec<_>>();
        for (&a, &b) in points.iter().zip(points.iter().cycle().skip(1)) {
            assert!(profile.line(a, b).ok_or("edge")?.upper_bound < 0.01);
        }
        let filled = profile.facet(&points, 0.01).ok_or("filled disk")?;
        assert!(filled.witnessed_max > 0.5);
        assert!(filled.upper_bound >= 1.0 - 1.0e-12);
        assert!(filled.upper_bound >= filled.witnessed_max);
        Ok(())
    }

    #[test]
    fn filled_facet_bounds_enclose_scattered_toroid_samples()
    -> Result<(), Box<dyn std::error::Error>> {
        // All meridian quadrants, nonuniform azimuth/meridian sampling,
        // and a large world-coordinate translation.
        let center = Point3::new(1000.0, -2000.0, 3000.0);
        let axis = Vector3::new(1.0, 2.0, 3.0).normalize();
        let (u, v) = normal_basis(axis);
        let profile = CircularProfile::new(center, axis, 4.0, 1.5).ok_or("profile")?;
        let point = |a: f64, b: f64| {
            let (sa, ca) = a.sin_cos();
            let (sb, cb) = b.sin_cos();
            center + (u * ca + v * sa) * 1.5_f64.mul_add(cb, 4.0) + axis * (1.5 * sb)
        };
        for i in 0..16 {
            let a = 0.33 * f64::from(i);
            let b = TAU * f64::from(i) / 16.0;
            let points = [
                point(a, b),
                point(a + 0.2, b + 0.04),
                point(a + 0.06, b + 0.3),
            ];
            let bound = profile.facet(&points, 0.1).ok_or("triangle bound")?;
            assert!(bound.upper_bound <= 0.1, "{bound:?}");
            for i in 0..=20 {
                for j in 0..=20 - i {
                    let sample = points[0]
                        + (points[1] - points[0]) * (f64::from(i) / 20.0)
                        + (points[2] - points[0]) * (f64::from(j) / 20.0);
                    assert!(profile.error(sample - center) <= bound.upper_bound + 1e-12);
                }
            }
        }
        Ok(())
    }

    #[test]
    fn toroidal_parallel_and_meridian_circles_are_contained()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(100.0, 200.0, -300.0);
        let axis = Vector3::unit_z();
        let profile = CircularProfile::new(center, axis, 4.0, 1.5).ok_or("profile")?;
        let outer = profile
            .circle(
                center,
                Vector3::unit_x() * 5.5,
                Vector3::unit_y() * 5.5,
                5.5,
            )
            .ok_or("outer")?;
        assert!(outer.upper_bound < 1.0e-10);
        let cross_section = profile
            .circle(
                center + Vector3::unit_x() * 4.0,
                Vector3::unit_x() * 1.5,
                axis * 1.5,
                1.5,
            )
            .ok_or("meridian")?;
        assert!(cross_section.upper_bound < 1.0e-10);
        let chord = profile
            .line(
                center + Vector3::unit_x() * 5.5,
                center - Vector3::unit_x() * 5.5,
            )
            .ok_or("chord")?;
        assert!(chord.witnessed_max > 2.0);
        assert!(chord.upper_bound >= chord.witnessed_max);
        Ok(())
    }

    #[test]
    fn squared_profile_does_not_silently_accept_wrong_radial_branch()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(0.0, 0.0, 0.0);
        let axis = Vector3::unit_z();
        let profile = CircularProfile::new(center, axis, 1.0, 2.0).ok_or("profile")?;
        let inner_circle = profile
            .circle(center, Vector3::unit_x(), Vector3::unit_y(), 1.0)
            .ok_or("wrong branch")?;
        assert!(inner_circle.witnessed_max > 1.99);
        assert!(inner_circle.upper_bound >= inner_circle.witnessed_max);
        let negative = CircularProfile::new(center, axis, -1.0, 2.0).ok_or("negative profile")?;
        let fit = negative
            .circle(center, Vector3::unit_x(), Vector3::unit_y(), 1.0)
            .ok_or("right branch")?;
        assert!(fit.upper_bound < 1.0e-12);
        Ok(())
    }

    #[test]
    fn arbitrary_circle_and_segment_are_continuously_bounded()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(1.0, 2.0, 3.0);
        let profile =
            CircularProfile::new(center, Vector3::new(0.2, 0.3, 1.0).normalize(), 2.0, 1.0)
                .ok_or("profile")?;
        let curve_center = Point3::new(2.0, 1.0, 5.0);
        let (u, v) = normal_basis(Vector3::new(0.3, -0.4, 0.5).normalize());
        let radius = 1.2;
        let u = u * radius;
        let v = v * radius;
        let circle = profile.circle(curve_center, u, v, radius).ok_or("circle")?;
        let start = curve_center + u;
        let end = curve_center + v;
        let line = profile.line(start, end).ok_or("line")?;
        for i in 0..=4096 {
            let t = f64::from(i) / 4096.0;
            let (sine, cosine) = (t * TAU).sin_cos();
            assert!(
                profile.error(curve_center - center + u * cosine + v * sine)
                    <= circle.upper_bound + 1.0e-12
            );
            assert!(
                profile.error(start - center + (end - start) * t) <= line.upper_bound + 1.0e-12
            );
        }
        Ok(())
    }
}
