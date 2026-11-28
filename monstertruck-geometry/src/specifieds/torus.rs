#![expect(
    clippy::similar_names,
    reason = "paired U/V ranges, tolerances, and divisions are standard torus parameterization notation"
)]

use super::{
    AffinePoint, Bound, BoundedSurface, ComponentWise, Error, KERNEL_EPSILON, MetricSpace,
    ParameterDivision2D, ParameterRange, ParametricSurface, ParametricSurface3D, Point3, Result,
    SearchNearestParameter, SearchParameter, SearchParameterHint2D, SurfaceParameter,
    SurfaceParameterRange, Torus, Vector2, Vector3, Zero,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
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
    #[must_use = "The newly constructed value is returned rather than retained by this constructor."]
    pub const fn new_unchecked(center: Point3, large_radius: f64, small_radius: f64) -> Self {
        Self {
            center,
            large_radius,
            small_radius,
        }
    }

    /// get center
    #[inline]
    #[must_use = "The geometric center is returned; ignoring it discards the placement."]
    pub const fn center(&self) -> Point3 { self.center }

    /// get large radius
    #[inline]
    #[must_use = "The computed large radius quantity is reported through this return value."]
    pub const fn large_radius(&self) -> f64 { self.large_radius }

    /// get small radius
    #[inline]
    #[must_use = "The computed small radius quantity is reported through this return value."]
    pub const fn small_radius(&self) -> f64 { self.small_radius }
}

impl ParametricSurface for Torus {
    type Point = Point3;
    type Vector = Vector3;
    #[inline]
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Self::Vector {
        let ((su, cu), (sv, cv)) = (u.sin_cos(), v.sin_cos());
        let center = match (m, n) {
            (0, 0) => self.center.position_vector(),
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
        center + u_part.component_product(v_part)
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
            d = 2.0f64.mul_add(-PI, d);
        }
        d
    }
    let du = delta(u, hint.0);
    let dv = delta(v, hint.1);
    du.mul_add(du, dv * dv)
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
        hint.map_or_else(
            || a.1.total_cmp(&b.1).then_with(|| a.0.total_cmp(&b.0)),
            |hint| {
                periodic_parameter_distance2((a.0, a.1), hint)
                    .total_cmp(&periodic_parameter_distance2((b.0, b.1), hint))
            },
        )
    });
    Some((candidates[0].0, candidates[0].1))
}

impl Torus {
    /// Exact inversion uses the same signed-meridian projection as the
    /// nearest-point solver, with the existing geometric residual gate.
    /// Both meridian coordinates determine the angle through `atan2`;
    /// using only `asin(z / r)` is ill-conditioned at the axial extrema
    /// and can reject roundoff-sized deviations before the residual test.
    fn exact_parameter_candidates(
        &self,
        point: Point3,
        hint: Option<(f64, f64)>,
    ) -> Vec<(f64, f64, f64)> {
        let relative = point - self.center();
        let scale = (self.large_radius() + self.small_radius() + relative.norm()).max(1.0);
        let tolerance2 = (KERNEL_EPSILON * scale).powi(2);
        if !tolerance2.is_finite() {
            return Vec::new();
        }
        let mut candidates = self.nearest_parameter_candidates(point, hint);
        candidates.retain(|candidate| candidate.2 <= tolerance2);
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
        let rho = Vector2::new(relative.x, relative.y).norm();
        let scale = (self.large_radius() + self.small_radius() + relative.norm()).max(1.0);
        let eps = 128.0 * f64::EPSILON * scale;
        let azimuth = if rho > eps {
            relative.y.atan2(relative.x).rem_euclid(2.0 * PI)
        } else {
            hint.map_or(0.0, |h| h.0)
        };
        let mut raw = Vec::<(f64, f64)>::new();
        for sign in [1.0_f64, -1.0] {
            let target_q = sign * rho;
            let dx = target_q - self.large_radius();
            let dz = relative.z;
            let norm = dx.hypot(dz);
            if norm > eps {
                let v = dz.atan2(dx).rem_euclid(2.0 * PI);
                let q = self.small_radius().mul_add(v.cos(), self.large_radius());
                if sign * q >= -eps {
                    let u = if sign > 0.0 { azimuth } else { azimuth + PI };
                    raw.push((u.rem_euclid(2.0 * PI), v));
                }
            } else {
                let v = hint.map_or(0.0, |h| h.1.rem_euclid(2.0 * PI));
                let q = self.small_radius().mul_add(v.cos(), self.large_radius());
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
            raw.push((azimuth, 2.0f64.mul_add(PI, -axis_v).rem_euclid(2.0 * PI)));
        }

        let mut candidates = Vec::with_capacity(raw.len());
        for (u0, v0) in raw {
            let u = nearest_periodic_angle(u0, hint.map(|h| h.0));
            let v = nearest_periodic_angle(v0, hint.map(|h| h.1));
            let distance2 = self.evaluate(u, v).distance2(point);
            if distance2.is_finite()
                && !candidates.iter().any(|candidate: &(f64, f64, f64)| {
                    periodic_parameter_distance2((candidate.0, candidate.1), (u, v))
                        <= DEFAULT_NUMERICAL_TOLERANCES.structural.powi(2)
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
        // Estimate the complete grid before allocating either axis. An
        // unwrapped or otherwise huge chart must not allocate an unbounded
        // circular polyline. Use the common budgeted adaptive divider if
        // the closed-form grid does not fit.
        let count = |span: f64, radius: f64| {
            let ratio = tol / radius;
            if !ratio.is_finite() || ratio <= 0.0 {
                return None;
            }
            // Stable even when the chord/radius ratio is too small for
            // subtracting it from one: acos(1-x) = 2 asin(sqrt(x/2)).
            let step = 4.0 * (ratio.min(0.8) * 0.5).sqrt().asin();
            crate::strict_segment_count(span, step)
        };
        let counts = count(
            urange.1 - urange.0,
            self.small_radius() + self.large_radius(),
        )
        .zip(count(vrange.1 - vrange.0, self.small_radius()));
        let max_cells =
            usize::try_from(monstertruck_traits::algo::surface::MAX_PARAMETER_DIVISION_CELLS)
                .unwrap_or(usize::MAX);
        let Some((nu, nv)) =
            counts.filter(|&(nu, nv)| nu.checked_mul(nv).is_some_and(|n| n <= max_cells))
        else {
            return monstertruck_traits::algo::surface::parameter_division(
                self,
                (urange, vrange),
                tol,
            );
        };
        let division = |range: (f64, f64), count: usize| {
            (0..=count)
                .map(|i| {
                    let t = crate::usize_to_f64(i) / crate::usize_to_f64(count);
                    range.1.mul_add(t, range.0 * (1.0 - t))
                })
                .collect()
        };
        (division(urange, nu), division(vrange, nv))
    }
}
