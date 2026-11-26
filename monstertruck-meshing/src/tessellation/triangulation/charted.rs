//! Exact C0 chart tiling for closed, intrinsically bounded spline surfaces.
//! Never use polynomial extrapolation as a substitute for periodic evaluation.

use super::PreMeshableSurface;
use monstertruck_geometry::prelude::{
    ParameterDivision2D, ParameterRange, ParametricSurface, ParametricSurface3D, Point2, Point3,
    Vector3, algo,
};

#[derive(Clone, Copy)]
pub(super) struct TrimChart<'a, S> {
    original: &'a S,
    /// Only newly certified finite chart periods, not analytic periods.
    wrap: [Option<(f64, f64)>; 2],
}

impl<'a, S: PreMeshableSurface> TrimChart<'a, S> {
    pub(super) fn new(original: &'a S) -> Self {
        let (u, v) = original.try_range_tuple();
        let (period_u, period_v) = original.closed_chart_periods();
        let (real_u, real_v) = original.parameter_bounds_are_intrinsic();
        let wrap_axis =
            |range: Option<(f64, f64)>, period: Option<f64>, analytic: Option<f64>, real: bool| {
                let (a, b) = range?;
                let period = period?;
                (real
                    && analytic.is_none()
                    && period.is_finite()
                    && period > 0.0
                    && (b - a).to_bits() == period.to_bits())
                .then_some((a, b))
            };
        Self {
            original,
            wrap: [
                wrap_axis(u, period_u, original.period_u(), real_u),
                wrap_axis(v, period_v, original.period_v(), real_v),
            ],
        }
    }

    pub(super) const fn original(&self) -> &'a S { self.original }

    pub(super) fn is_wrapped(&self) -> bool { self.wrap.iter().any(Option::is_some) }

    fn canonical_axis(value: f64, range: Option<(f64, f64)>) -> f64 {
        if let Some((a, b)) = range
            && !(a..=b).contains(&value)
        {
            a + (value - a).rem_euclid(b - a)
        } else {
            value
        }
    }

    pub(super) fn canonical_uv(&self, (u, v): (f64, f64)) -> (f64, f64) {
        (
            Self::canonical_axis(u, self.wrap[0]),
            Self::canonical_axis(v, self.wrap[1]),
        )
    }

    fn seeds(range: (f64, f64), wrap: Option<(f64, f64)>) -> Option<Vec<f64>> {
        let mut result = vec![range.0];
        if let Some((lo, hi)) = wrap {
            let period = hi - lo;
            let mut lap = ((range.0 - lo) / period).floor() + 1.0;
            let limit = f64::from(u32::try_from(algo::surface::MAX_PARAMETER_DIVISION_CELLS).ok()?);
            if !lap.is_finite() || (range.1 - range.0) / period > limit {
                return None;
            }
            loop {
                let seam = lap.mul_add(period, lo);
                if seam >= range.1 {
                    break;
                }
                if seam > range.0 {
                    result.push(seam);
                }
                let next = lap + 1.0;
                if next <= lap
                    || result.len()
                        >= usize::try_from(algo::surface::MAX_PARAMETER_DIVISION_CELLS).ok()?
                {
                    return None;
                }
                lap = next;
            }
        }
        result.push(range.1);
        Some(result)
    }
}

impl<S: PreMeshableSurface> ParametricSurface for TrimChart<'_, S> {
    type Point = Point3;
    type Vector = Vector3;

    fn evaluate(&self, u: f64, v: f64) -> Point3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.evaluate(u, v)
    }
    fn derivative_u(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.derivative_u(u, v)
    }
    fn derivative_v(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.derivative_v(u, v)
    }
    fn derivative_uu(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.derivative_uu(u, v)
    }
    fn derivative_uv(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.derivative_uv(u, v)
    }
    fn derivative_vv(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.derivative_vv(u, v)
    }
    fn derivative_mn(&self, m: usize, n: usize, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.derivative_mn(m, n, u, v)
    }
    fn evaluate_with_first_derivatives(&self, u: f64, v: f64) -> (Point3, Vector3, Vector3) {
        let (u, v) = self.canonical_uv((u, v));
        self.original.evaluate_with_first_derivatives(u, v)
    }
    fn evaluate_with_second_derivatives(
        &self,
        u: f64,
        v: f64,
    ) -> (Point3, Vector3, Vector3, Vector3, Vector3, Vector3) {
        let (u, v) = self.canonical_uv((u, v));
        self.original.evaluate_with_second_derivatives(u, v)
    }
    fn parameter_range(&self) -> (ParameterRange, ParameterRange) {
        self.original.parameter_range()
    }
    fn parameter_bounds_are_intrinsic(&self) -> (bool, bool) {
        self.original.parameter_bounds_are_intrinsic()
    }
    fn period_u(&self) -> Option<f64> {
        self.wrap[0]
            .map(|(a, b)| b - a)
            .or_else(|| self.original.period_u())
    }
    fn period_v(&self) -> Option<f64> {
        self.wrap[1]
            .map(|(a, b)| b - a)
            .or_else(|| self.original.period_v())
    }
}

impl<S: PreMeshableSurface> ParametricSurface3D for TrimChart<'_, S> {
    fn normal(&self, u: f64, v: f64) -> Vector3 {
        let (u, v) = self.canonical_uv((u, v));
        self.original.normal(u, v)
    }

    fn collapsed_boundary_parameter(
        &self,
        point: Point3,
        hint: Point2,
        tolerance: f64,
    ) -> Option<Point2> {
        let (u, v) = self.canonical_uv((hint.x, hint.y));
        let mut uv =
            self.original
                .collapsed_boundary_parameter(point, Point2::new(u, v), tolerance)?;
        for axis in 0..2 {
            if let Some((a, b)) = self.wrap[axis] {
                let period = b - a;
                uv[axis] = ((hint[axis] - uv[axis]) / period)
                    .round()
                    .mul_add(period, uv[axis]);
            }
        }
        Some(uv)
    }
}

impl<S: PreMeshableSurface> ParameterDivision2D for TrimChart<'_, S> {
    fn parameter_division(
        &self,
        range: ((f64, f64), (f64, f64)),
        tol: f64,
    ) -> (Vec<f64>, Vec<f64>) {
        let outside = |requested: (f64, f64), wrap: Option<(f64, f64)>| {
            wrap.is_some_and(|(a, b)| requested.0 < a || requested.1 > b)
        };
        if !outside(range.0, self.wrap[0]) && !outside(range.1, self.wrap[1]) {
            return self.original.parameter_division(range, tol);
        }
        let budget = algo::surface::MAX_PARAMETER_DIVISION_CELLS;
        let grid = Self::seeds(range.0, self.wrap[0])
            .zip(Self::seeds(range.1, self.wrap[1]))
            .filter(|(u, v)| ((u.len() - 1) as u64).saturating_mul((v.len() - 1) as u64) <= budget);
        grid.map_or_else(
            || {
                // Report exhaustion through the existing division accounting
                // without attempting to allocate an unbounded seed grid.
                algo::surface::parameter_division_with_budget(self, range, tol, 0).division
            },
            |grid| {
                algo::surface::parameter_division_from_grid_with_budget(self, grid, tol, budget)
                    .division
            },
        )
    }
}
