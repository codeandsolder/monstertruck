use crate::v2;
use monstertruck_core::MetricSpace;
use num_traits::Float;

fn scalar_from_usize<S: Float>(mut value: usize) -> S {
    let mut result = S::zero();
    let mut bit = S::one();
    while value != 0 {
        if value & 1 == 1 {
            result = result + bit;
        }
        value >>= 1;
        if value != 0 {
            bit = bit + bit;
        }
    }
    result
}

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
