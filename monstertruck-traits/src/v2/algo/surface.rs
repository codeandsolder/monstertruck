use super::scalar_from_usize;
use crate::v2;
use monstertruck_core::MetricSpace;
use num_traits::Float;
/// Divides the domain into equal parts and returns `(u, v)` such that
/// `surface.evaluate(u, v)` is closest to `point`.
///
/// Scalar-generic version of [`crate::algo::surface::presearch`].
pub fn presearch<S, Surf>(
    surface: &Surf,
    point: Surf::Point,
    (urange, vrange): ((S, S), (S, S)),
    division: usize,
) -> (S, S)
where
    S: Float,
    Surf: v2::ParametricSurface<Scalar = S>,
    Surf::Point: MetricSpace<Metric = S> + Copy,
{
    let ((u0, u1), (v0, v1)) = (urange, vrange);
    if division == 0 {
        let _ = surface.evaluate(u0, v0);
        return (u0, v0);
    }
    let mut res = (u0, v0);
    let mut min = S::infinity();
    let n: S = scalar_from_usize(division);
    for i in 0..=division {
        for j in 0..=division {
            let p = scalar_from_usize::<S>(i) / n;
            let q = scalar_from_usize::<S>(j) / n;
            let u = u0 * (S::one() - p) + u1 * p;
            let v = v0 * (S::one() - q) + v1 * q;
            let dist = surface.evaluate(u, v).distance2(point);
            if dist < min {
                min = dist;
                res = (u, v);
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_core::cgmath64::{Point3, Vector3};
    use std::{cell::Cell, rc::Rc};

    #[derive(Clone)]
    struct RecordingSurface(Rc<Cell<Option<(f64, f64)>>>);

    impl v2::ParametricSurface for RecordingSurface {
        type Scalar = f64;
        type Point = Point3;
        type Vector = Vector3;

        fn evaluate(&self, u: f64, v: f64) -> Point3 {
            self.0.set(Some((u, v)));
            Point3::new(u, v, 0.0)
        }
        fn derivative_u(&self, _: f64, _: f64) -> Vector3 { Vector3::new(0.0, 0.0, 0.0) }
        fn derivative_v(&self, _: f64, _: f64) -> Vector3 { Vector3::new(0.0, 0.0, 0.0) }
        fn derivative_uu(&self, _: f64, _: f64) -> Vector3 { Vector3::new(0.0, 0.0, 0.0) }
        fn derivative_uv(&self, _: f64, _: f64) -> Vector3 { Vector3::new(0.0, 0.0, 0.0) }
        fn derivative_vv(&self, _: f64, _: f64) -> Vector3 { Vector3::new(0.0, 0.0, 0.0) }
        fn period_u(&self) -> Option<f64> { None }
        fn period_v(&self) -> Option<f64> { None }
    }

    #[test]
    fn zero_division_presearch_samples_range_start_once() {
        let seen = Rc::new(Cell::new(None));
        let surface = RecordingSurface(seen.clone());
        let point = Point3::new(f64::NAN, 0.0, 0.0);

        assert_eq!(
            presearch(&surface, point, ((2.0, 3.0), (4.0, 5.0)), 0),
            (2.0, 4.0)
        );
        assert_eq!(seen.get(), Some((2.0, 4.0)));
    }
}
