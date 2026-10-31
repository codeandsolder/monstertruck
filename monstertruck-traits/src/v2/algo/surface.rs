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
        return (u0, v0);
    }
    let mut res = (S::zero(), S::zero());
    let mut min = S::infinity();
    let n = S::from(division).unwrap();
    for i in 0..=division {
        for j in 0..=division {
            let p = S::from(i).unwrap() / n;
            let q = S::from(j).unwrap() / n;
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

    #[derive(Clone)]
    struct PanicSurface;

    impl v2::ParametricSurface for PanicSurface {
        type Scalar = f64;
        type Point = Point3;
        type Vector = Vector3;

        fn evaluate(&self, _: f64, _: f64) -> Point3 { panic!("zero-division presearch must not evaluate") }
        fn derivative_u(&self, _: f64, _: f64) -> Vector3 { unreachable!() }
        fn derivative_v(&self, _: f64, _: f64) -> Vector3 { unreachable!() }
        fn derivative_uu(&self, _: f64, _: f64) -> Vector3 { unreachable!() }
        fn derivative_uv(&self, _: f64, _: f64) -> Vector3 { unreachable!() }
        fn derivative_vv(&self, _: f64, _: f64) -> Vector3 { unreachable!() }
        fn period_u(&self) -> Option<f64> { None }
        fn period_v(&self) -> Option<f64> { None }
    }

    #[test]
    fn zero_division_presearch_returns_range_start() {
        let point = Point3::new(0.0, 0.0, 0.0);
        assert_eq!(presearch(&PanicSurface, point, ((2.0, 3.0), (4.0, 5.0)), 0), (2.0, 4.0));
    }
}
