//! Continuous restrictions of a circular meridian to a line or circle.
//! Squaring the meridian equation produces a quartic whose circle restriction
//! is only degree two in the angular Fourier basis. The unsquared branch is
//! checked separately; a small quartic residual alone is insufficient.

use super::{CarrierResidual, residual, unit};
use crate::{InnerSpace, Point3, Vector3};
use std::f64::consts::FRAC_PI_4;

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

    pub(super) fn circle(
        self,
        center: Point3,
        basis_u: Vector3,
        basis_v: Vector3,
        radius: f64,
    ) -> Option<CarrierResidual> {
        let origin = center - self.center;
        // For orthogonal equal-length circle basis vectors, |delta|² has
        // only constant and first harmonics, even on a tilted circle.
        let constant = radius.mul_add(radius, origin.magnitude2()) + self.offset();
        let cosine = 2.0 * origin.dot(basis_u);
        let sine = 2.0 * origin.dot(basis_v);
        let amplitude = cosine.hypot(sine);
        let a_squared = Harmonics::square_first(constant, cosine, sine);
        let radial_squared = Harmonics::norm_squared(
            self.transverse(origin),
            self.transverse(basis_u),
            self.transverse(basis_v),
        );
        let polynomial = a_squared.add_scaled(radial_squared, -4.0 * self.radial_center.powi(2));
        let witness = (0..8)
            .map(|i| {
                let (sin, cos) = (f64::from(i) * FRAC_PI_4).sin_cos();
                self.error(origin + basis_u * cos + basis_v * sin)
            })
            .fold(0.0_f64, f64::max);
        self.bound(
            [constant - amplitude, constant + amplitude],
            polynomial.abs_bound(),
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

    fn abs_bound(self) -> f64 {
        let [constant, cosine, sine, double_cosine, double_sine] = self.0;
        constant.abs() + cosine.hypot(sine) + double_cosine.hypot(double_sine)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::sample_fit::normal_basis;
    use std::f64::consts::TAU;

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
