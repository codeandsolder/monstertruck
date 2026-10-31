use super::scalar_from_usize;
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
        let _ = curve.evaluate(t0);
        return t0;
    }
    let mut res = t0;
    let mut min = S::infinity();
    let n: S = scalar_from_usize(division);
    for i in 0..=division {
        let p = scalar_from_usize::<S>(i) / n;
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
    use std::{cell::Cell, rc::Rc};

    #[derive(Clone)]
    struct RecordingCurve(Rc<Cell<Option<f64>>>);

    impl v2::ParametricCurve for RecordingCurve {
        type Scalar = f64;
        type Point = Point2;
        type Vector = Vector2;

        fn evaluate(&self, t: f64) -> Point2 {
            self.0.set(Some(t));
            Point2::new(t, 0.0)
        }
        fn derivative(&self, _: f64) -> Vector2 { Vector2::new(0.0, 0.0) }
        fn derivative_2(&self, _: f64) -> Vector2 { Vector2::new(0.0, 0.0) }
        fn derivative_n(&self, _: usize, _: f64) -> Vector2 { Vector2::new(0.0, 0.0) }
        fn period(&self) -> Option<f64> { None }
        fn try_range_tuple(&self) -> Option<(f64, f64)> { None }
    }

    #[test]
    fn zero_division_presearch_samples_range_start_once() {
        let seen = Rc::new(Cell::new(None));
        let curve = RecordingCurve(seen.clone());
        let point = Point2::new(2.5, 0.0);

        assert_eq!(presearch(&curve, point, (2.5, 9.0), 0), 2.5);
        assert_eq!(seen.get(), Some(2.5));
    }
}
