#![expect(
    clippy::similar_names,
    reason = "paired U/V ranges, tolerances, and divisions are standard torus parameterization notation"
)]

use super::{
    Bound, BoundedSurface, ElementWise, Error, EuclideanSpace, InnerSpace, ParameterDivision1D,
    ParameterDivision2D, ParameterRange, ParametricSurface, ParametricSurface3D, Point2, Point3,
    Result, SearchNearestParameter, SearchParameter, SearchParameterHint2D, SurfaceParameter,
    SurfaceParameterRange, TOLERANCE, Torus, UnitCircle, Vector2, Vector3, Zero,
};
use std::f64::consts::PI;

impl Torus {
    /// constructor
    /// # Errors
    /// Returns [`Error::InvalidTorusRadii`] unless both radii are finite and strictly positive.
    ///
    #[inline]
    pub fn new(center: Point3, large_radius: f64, small_radius: f64) -> Result<Self> {
        if !large_radius.is_finite()
            || !small_radius.is_finite()
            || large_radius <= 0.0
            || small_radius <= 0.0
        {
            return Err(Error::InvalidTorusRadii(large_radius, small_radius));
        }
        Ok(Self {
            center,
            large_radius,
            small_radius,
        })
    }

    /// Constructs a torus without validating its radii.
    ///
    /// Use this only when the caller has already established that both radii
    /// are finite and strictly positive (for example, fixed test geometry).
    #[inline]
    #[must_use]
    pub const fn new_unchecked(center: Point3, large_radius: f64, small_radius: f64) -> Self {
        Self {
            center,
            large_radius,
            small_radius,
        }
    }

    /// get center
    #[inline]
    #[must_use]
    pub const fn center(&self) -> Point3 { self.center }

    /// get large radius
    #[inline]
    #[must_use]
    pub const fn large_radius(&self) -> f64 { self.large_radius }

    /// get small radius
    #[inline]
    #[must_use]
    pub const fn small_radius(&self) -> f64 { self.small_radius }
}

impl ParametricSurface for Torus {
    type Point = Point3;
    type Vector = Vector3;
    #[inline]
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        let ((su, cu), (sv, cv)) = (u.sin_cos(), v.sin_cos());
        let center = match (m, n) {
            (0, 0) => self.center.to_vec(),
            _ => Vector3::zero(),
        };
        let u_z = if m == 0 { 1.0 } else { 0.0 };
        let u_part = match m % 4 {
            0 => Vector3::new(cu, su, u_z),
            1 => Vector3::new(-su, cu, 0.0),
            2 => Vector3::new(-cu, -su, 0.0),
            _ => Vector3::new(su, -cu, 0.0),
        };
        let r0 = if n == 0 { self.large_radius } else { 0.0 };
        let r1 = self.small_radius;
        let v_part_d2 = match n % 4 {
            0 => Vector2::new(r1.mul_add(cv, r0), r1 * sv),
            1 => Vector2::new(-r1 * sv, r1 * cv),
            2 => Vector2::new(-r1 * cv, -r1 * sv),
            _ => Vector2::new(r1 * sv, -r1 * cv),
        };
        let v_part = Vector3::new(v_part_d2.x, v_part_d2.x, v_part_d2.y);
        center + u_part.mul_element_wise(v_part)
    }
    #[inline]
    fn evaluate(&self, u: f64, v: f64) -> Point3 {
        let sr = self.small_radius() * Vector2::new(f64::cos(v), f64::sin(v));
        let lr = (self.large_radius() + sr.x) * Vector2::new(f64::cos(u), f64::sin(u));
        self.center() + Vector3::new(lr.x, lr.y, sr.y)
    }
    #[inline]
    fn derivative_u(&self, u: f64, v: f64) -> Vector3 {
        let sr = self.small_radius() * f64::cos(v);
        let lr = (self.large_radius() + sr) * Vector2::new(f64::cos(u), f64::sin(u));
        Vector3::new(-lr.y, lr.x, 0.0)
    }
    #[inline]
    fn derivative_v(&self, u: f64, v: f64) -> Vector3 {
        let sv = self.small_radius() * Vector2::new(-f64::sin(v), f64::cos(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y)
    }
    #[inline]
    fn derivative_uu(&self, u: f64, v: f64) -> Vector3 {
        let sr = self.small_radius() * f64::cos(v);
        let lr = (self.large_radius() + sr) * Vector2::new(f64::cos(u), f64::sin(u));
        Vector3::new(-lr.x, -lr.y, 0.0)
    }
    #[inline]
    fn derivative_uv(&self, u: f64, v: f64) -> Vector3 {
        let sr = -self.small_radius() * f64::sin(v);
        let lr = sr * Vector2::new(f64::cos(u), f64::sin(u));
        Vector3::new(-lr.y, lr.x, 0.0)
    }
    #[inline]
    fn derivative_vv(&self, u: f64, v: f64) -> Vector3 {
        let sv = -self.small_radius() * Vector2::new(f64::cos(v), f64::sin(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y)
    }
    #[inline]
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        const RANGE: (Bound<f64>, Bound<f64>) = (Bound::Included(0.0), Bound::Excluded(2.0 * PI));
        (RANGE, RANGE)
    }
    #[inline]
    fn period_u(&self) -> Option<f64> { Some(2.0 * PI) }
    #[inline]
    fn period_v(&self) -> Option<f64> { Some(2.0 * PI) }
}

impl ParametricSurface3D for Torus {
    #[inline]
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        let sv = Vector2::new(f64::cos(v), f64::sin(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y)
    }
    #[inline]
    fn normal_uder(&self, u: f64, v: f64) -> Vector3 {
        let sv = Vector2::new(f64::cos(v), f64::sin(v));
        Vector3::new(-sv.x * f64::sin(u), sv.x * f64::cos(u), sv.y)
    }
    #[inline]
    fn normal_vder(&self, u: f64, v: f64) -> Vector3 {
        let sv = Vector2::new(-f64::sin(v), f64::cos(v));
        Vector3::new(sv.x * f64::cos(u), sv.x * f64::sin(u), sv.y)
    }
}

impl BoundedSurface for Torus {
    fn range_tuple(&self) -> SurfaceParameterRange {
        let range = (0.0, 2.0 * PI);
        (range, range)
    }
}

// -- v2 scalar-generic impls ------------------------------------------------

use monstertruck_traits::v2;

impl v2::ParametricSurface for Torus {
    type Scalar = f64;
    type Point = Point3;
    type Vector = Vector3;

    #[inline]
    fn evaluate(&self, u: f64, v: f64) -> Point3 { ParametricSurface::evaluate(self, u, v) }
    #[inline]
    fn derivative_u(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_u(self, u, v)
    }
    #[inline]
    fn derivative_v(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_v(self, u, v)
    }
    #[inline]
    fn derivative_uu(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_uu(self, u, v)
    }
    #[inline]
    fn derivative_uv(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_uv(self, u, v)
    }
    #[inline]
    fn derivative_vv(&self, u: f64, v: f64) -> Vector3 {
        ParametricSurface::derivative_vv(self, u, v)
    }
    #[inline]
    fn period_u(&self) -> Option<f64> { ParametricSurface::period_u(self) }
    #[inline]
    fn period_v(&self) -> Option<f64> { ParametricSurface::period_v(self) }
}

impl v2::BoundedSurface for Torus {
    #[inline]
    fn range_tuple(&self) -> ((f64, f64), (f64, f64)) { BoundedSurface::range_tuple(self) }
}

impl v2::ParametricSurface3D for Torus {}

impl v2::SearchNearestParameter<v2::SurfaceParameter<f64>> for Torus {
    type Point = Point3;
    #[inline]
    fn search_nearest_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: Point3,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchNearestParameter::<SurfaceParameter>::search_nearest_parameter(self, pt, None, trials)
    }
}

impl v2::SearchParameter<v2::SurfaceParameter<f64>> for Torus {
    type Point = Point3;
    #[inline]
    fn search_parameter<H: Into<v2::SearchParameterHint2D<f64>>>(
        &self,
        pt: Point3,
        _: H,
        trials: usize,
    ) -> Option<(f64, f64)> {
        SearchParameter::<SurfaceParameter>::search_parameter(self, pt, None, trials)
    }
}

/// The `(u, v)` a [`SearchParameterHint2D`] carries, or `None`.
///
/// `Range` answers `None`: it names a rectangle, not a point, and the seam
/// disambiguation below wants an anchor. Every caller that needs the
/// disambiguation (the boundary projection chain, the SSI's parameter refiner)
/// passes `Parameter`.
#[inline]
const fn hint_parameter(hint: SearchParameterHint2D) -> Option<(f64, f64)> {
    match hint {
        SearchParameterHint2D::Parameter(u, v) => Some((u, v)),
        SearchParameterHint2D::Range(..) | SearchParameterHint2D::None => None,
    }
}

/// The co-periodic representative of `angle` nearest `hint` -- `angle + 2 k PI`
/// for the integer `k` that minimises the distance -- or `angle` unchanged when
/// no hint is available.
///
/// # The defect this closes, measured (spec 012 W1)
///
/// A torus is exactly `2 * PI`-periodic in both parameters, so every angle has
/// infinitely many correct spellings and this routine's branch arithmetic picks
/// one by a rule that reads only the CURRENT point: `2 * PI - acos(x)` is
/// exactly `2 * PI` when the point sits on the `u = 0` seam with a negative-zero
/// `y`, and exactly `0` when the same point carries a positive-zero `y`.
/// Ordinary for a seam vertex of a real STEP face, and the caller was never
/// given a way to say which spelling it wanted -- the `hint` argument was
/// literally discarded (`_: H`).
///
/// That is invisible until a torus reaches a caller that PROJECTS a boundary.
/// `sampled_parameter_boundary` (`monstertruck-modeling/src/geometry/`) walks
/// a face's boundary polyline through this routine, feeding each answer forward
/// as the next hint, so one seam vertex spelled a whole period away from its
/// neighbours leaves the face's parameter LOOP with a spike across the entire
/// domain. Measured on ap224 with spec 012's analytic-torus routing on -- two
/// `R = 0.9107, r = 0.0312` fillet tori, seam-touching from OPPOSITE sides:
///
/// | face | true `u` extent (from the rational net) | analytic answer, before |
/// |---|---|---|
/// | 15 | `(0, PI)` | `(0.0645, 6.2832)` -- the whole ring |
/// | 19 | `(PI, 2 PI)` | `(0, 6.2187)` -- the whole ring |
///
/// The SSI then traced the RIGHT curves against the box cutter's planes and the
/// trim filter discarded every one of them (`side0 = 0` on all 8 segments),
/// which surfaced as `CreateLoopsStoreFailed{IntersectionCurvesFailed{(15,4),
/// SsiFailed}}`. **So this is ledger class C4, not a missing analytic SSI path:
/// the kernel intersected the analytic torus perfectly well and was handed a
/// face whose own trim loop excluded its interior.**
///
/// **A fixed rule cannot work here, and that is a measurement, not a guess.**
/// The first attempt at this reduced an exact `2 * PI` to `0` -- correct for
/// face 15, and it moved the refusal straight onto face 19, whose seam vertex
/// needs the `2 * PI` spelling. The two faces are the same shape on opposite
/// sides of one seam, so only the CALLER's context distinguishes them, and the
/// hint is that context.
///
/// **Unhinted callers are byte-identical.** `SearchParameterHint2D::None`
/// returns `angle` untouched, so every call site that passes no hint -- which is
/// every call site that existed before spec 012 routed tori here -- sees exactly
/// the arithmetic it saw before.
#[inline]
fn nearest_periodic_angle(angle: f64, hint: Option<f64>) -> f64 {
    hint.map_or(angle, |hint| {
        (2.0 * PI).mul_add(((hint - angle) / (2.0 * PI)).round(), angle)
    })
}

#[inline]
fn periodic_parameter_distance2((u, v): (f64, f64), hint: (f64, f64)) -> f64 {
    fn delta(a: f64, b: f64) -> f64 {
        let mut d = (a - b).rem_euclid(2.0 * PI);
        if d > PI {
            d -= 2.0 * PI;
        }
        d
    }
    delta(u, hint.0).powi(2) + delta(v, hint.1).powi(2)
}

fn choose_parameter_candidate(
    mut candidates: Vec<(f64, f64, f64)>,
    hint: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    if candidates.is_empty() {
        return None;
    }
    candidates.sort_by(|a, b| {
        let distance_order = a.2.total_cmp(&b.2);
        if distance_order != std::cmp::Ordering::Equal {
            return distance_order;
        }
        match hint {
            Some(hint) => periodic_parameter_distance2((a.0, a.1), hint)
                .total_cmp(&periodic_parameter_distance2((b.0, b.1), hint)),
            None => a.1.total_cmp(&b.1).then_with(|| a.0.total_cmp(&b.0)),
        }
    });
    Some((candidates[0].0, candidates[0].1))
}

impl Torus {
    /// Exact inverse candidates for the torus parameterisation.
    ///
    /// A spindle torus is not one-to-one in the meridian projection: when
    /// `R + r cos(v)` is negative, the same world-space azimuth is represented
    /// by `u + pi`. Solving only from `rho = hypot(x, y)` therefore loses the
    /// inner branch. Enumerating the two sine-compatible `v` values and using
    /// the signed radial factor recovers both ring and spindle regimes without
    /// weakening the final point-on-surface check.
    fn exact_parameter_candidates(
        &self,
        point: Point3,
        hint: Option<(f64, f64)>,
    ) -> Vec<(f64, f64, f64)> {
        let relative = point - self.center();
        let rho = Vector2::new(relative.x, relative.y).magnitude();
        let scale = (self.large_radius() + self.small_radius() + relative.magnitude()).max(1.0);
        let eps = 128.0 * f64::EPSILON * scale;
        let sine = relative.z / self.small_radius();
        if !sine.is_finite() || sine < -1.0 - eps || sine > 1.0 + eps {
            return Vec::new();
        }
        let alpha = f64::asin(sine.clamp(-1.0, 1.0));
        let mut v_candidates = vec![
            alpha.rem_euclid(2.0 * PI),
            (PI - alpha).rem_euclid(2.0 * PI),
        ];
        if (v_candidates[0] - v_candidates[1]).abs() <= eps {
            v_candidates.truncate(1);
        }
        let azimuth = if rho > eps {
            relative.y.atan2(relative.x).rem_euclid(2.0 * PI)
        } else {
            hint.map(|h| h.0).unwrap_or(0.0)
        };
        let mut candidates = Vec::with_capacity(v_candidates.len());
        for v0 in v_candidates {
            let signed_radius = self.large_radius() + self.small_radius() * v0.cos();
            if (rho - signed_radius.abs()).abs() > eps.max(TOLERANCE * scale) {
                continue;
            }
            let u0 = if signed_radius < 0.0 {
                azimuth + PI
            } else {
                azimuth
            };
            let u = nearest_periodic_angle(u0.rem_euclid(2.0 * PI), hint.map(|h| h.0));
            let v = nearest_periodic_angle(v0, hint.map(|h| h.1));
            let distance2 = self.evaluate(u, v).distance2(point);
            if distance2 <= (TOLERANCE * scale).powi(2) {
                candidates.push((u, v, distance2));
            }
        }
        candidates
    }

    /// Global nearest-point candidates reduced to the generating circle.
    ///
    /// For a world point with cylindrical radius `rho`, minimising over `u`
    /// turns the torus distance into the minimum of two 2D circle problems:
    /// target signed radii `+rho` and `-rho`. Each branch is valid only where
    /// the generating circle has the matching sign; if its unconstrained
    /// optimum crosses the axis, the branch boundary `R + r cos(v) = 0` is also
    /// tested. This remains exact in the spindle regime and degenerates to the
    /// familiar single positive-radial branch for a ring torus.
    fn nearest_parameter_candidates(
        &self,
        point: Point3,
        hint: Option<(f64, f64)>,
    ) -> Vec<(f64, f64, f64)> {
        let relative = point - self.center();
        let rho = Vector2::new(relative.x, relative.y).magnitude();
        let scale = (self.large_radius() + self.small_radius() + relative.magnitude()).max(1.0);
        let eps = 128.0 * f64::EPSILON * scale;
        let azimuth = if rho > eps {
            relative.y.atan2(relative.x).rem_euclid(2.0 * PI)
        } else {
            hint.map(|h| h.0).unwrap_or(0.0)
        };
        let mut raw = Vec::<(f64, f64)>::new();
        for sign in [1.0_f64, -1.0] {
            let target_q = sign * rho;
            let dx = target_q - self.large_radius();
            let dz = relative.z;
            let norm = dx.hypot(dz);
            if norm > eps {
                let v = dz.atan2(dx).rem_euclid(2.0 * PI);
                let q = self.large_radius() + self.small_radius() * v.cos();
                if sign * q >= -eps {
                    let u = if sign > 0.0 { azimuth } else { azimuth + PI };
                    raw.push((u.rem_euclid(2.0 * PI), v));
                }
            } else {
                let v = hint.map(|h| h.1.rem_euclid(2.0 * PI)).unwrap_or(0.0);
                let q = self.large_radius() + self.small_radius() * v.cos();
                if sign * q >= -eps {
                    let u = if sign > 0.0 { azimuth } else { azimuth + PI };
                    raw.push((u.rem_euclid(2.0 * PI), v));
                }
            }
        }
        if self.large_radius() <= self.small_radius() + eps {
            let c = (-self.large_radius() / self.small_radius()).clamp(-1.0, 1.0);
            let axis_v = c.acos();
            raw.push((azimuth, axis_v));
            raw.push((azimuth, (2.0 * PI - axis_v).rem_euclid(2.0 * PI)));
        }

        let mut candidates = Vec::with_capacity(raw.len());
        for (u0, v0) in raw {
            let u = nearest_periodic_angle(u0, hint.map(|h| h.0));
            let v = nearest_periodic_angle(v0, hint.map(|h| h.1));
            let distance2 = self.evaluate(u, v).distance2(point);
            if distance2.is_finite()
                && !candidates.iter().any(|candidate: &(f64, f64, f64)| {
                    periodic_parameter_distance2((candidate.0, candidate.1), (u, v)) <= 1.0e-24
                })
            {
                candidates.push((u, v, distance2));
            }
        }
        candidates
    }
}

impl SearchParameter<SurfaceParameter> for Torus {
    type Point = Point3;
    fn search_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Point3,
        hint: H,
        _: usize,
    ) -> Option<(f64, f64)> {
        let hint = hint_parameter(hint.into());
        choose_parameter_candidate(self.exact_parameter_candidates(point, hint), hint)
    }
}

impl SearchNearestParameter<SurfaceParameter> for Torus {
    type Point = Point3;
    fn search_nearest_parameter<H: Into<SearchParameterHint2D>>(
        &self,
        point: Point3,
        hint: H,
        _: usize,
    ) -> Option<(f64, f64)> {
        let hint = hint_parameter(hint.into());
        choose_parameter_candidate(self.nearest_parameter_candidates(point, hint), hint)
    }
}

impl ParameterDivision2D for Torus {
    fn parameter_division(
        &self,
        (urange, vrange): ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let circle = UnitCircle::<Point2>::new();
        let utol = tol / (self.small_radius() + self.large_radius());
        let (udiv, _) = circle.parameter_division(urange, utol);
        let vtol = tol / self.small_radius();
        let (vdiv, _) = circle.parameter_division(vrange, vtol);
        (udiv, vdiv)
    }
}
