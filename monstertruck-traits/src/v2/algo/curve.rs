use crate::v2;
use monstertruck_core::MetricSpace;
use num_traits::Float;

/// Divides the domain into equal parts and returns `t` such that
/// `curve.evaluate(t)` is closest to `point`.
///
/// Scalar-generic version of [`crate::algo::curve::presearch`].
pub fn presearch<S, C>(curve: &C, point: C::Point, range: (S, S), division: usize) -> S
where
    S: Float,
    C: v2::ParametricCurve<Scalar = S>,
    C::Point: MetricSpace<Metric = S> + Copy, {
    let (t0, t1) = range;
    if division == 0 {
        return t0;
    }
    let mut res = t0;
    let mut min = S::infinity();
    let n = S::from(division).unwrap();
    for i in 0..=division {
        let p = S::from(i).unwrap() / n;
        let t = t0 * (S::one() - p) + t1 * p;
        let dist = curve.evaluate(t).distance2(point);
        if dist < min {
            min = dist;
            res = t;
        }
    }
    res
}


#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_core::cgmath64::{Point2, Vector2};

    #[derive(Clone)]
    struct PanicCurve;

    impl v2::ParametricCurve for PanicCurve {
        type Scalar = f64;
        type Point = Point2;
        type Vector = Vector2;

        fn evaluate(&self, _: f64) -> Point2 { panic!("zero-division presearch must not evaluate") }
        fn derivative(&self, _: f64) -> Vector2 { unreachable!() }
        fn derivative_2(&self, _: f64) -> Vector2 { unreachable!() }
        fn derivative_n(&self, _: usize, _: f64) -> Vector2 { unreachable!() }
        fn period(&self) -> Option<f64> { None }
        fn try_range_tuple(&self) -> Option<(f64, f64)> { None }
    }

    #[test]
    fn zero_division_presearch_returns_range_start() {
        let point = Point2::new(0.0, 0.0);
        assert_eq!(presearch(&PanicCurve, point, (2.5, 9.0), 0), 2.5);
    }
}
