//! Intersection rejection against a convex footprint prism. Inflate by the
//! kernel incidence allowance: uncertain/touching unrelated geometry is not
//! silently removed. The original source faces are never modified here.
use super::planar;
use crate::KERNEL_EPSILON;

fn clip(poly: &[[f64; 3]], signed: impl Fn([f64; 3]) -> f64) -> Option<Vec<[f64; 3]>> {
    let Some(&last) = poly.last() else {
        return Some(Vec::new());
    };
    let mut previous = last;
    let mut prev = signed(last) + KERNEL_EPSILON;
    let mut result = Vec::new();
    for &point in poly {
        let current = signed(point) + KERNEL_EPSILON;
        if !prev.is_finite() || !current.is_finite() {
            return None;
        }
        if (prev >= 0.0) != (current >= 0.0) {
            let t = prev / (prev - current);
            if !t.is_finite() {
                return None;
            }
            result.push(std::array::from_fn(|i| {
                previous[i] + t * (point[i] - previous[i])
            }));
        }
        if current >= 0.0 {
            result.push(point);
        }
        previous = point;
        prev = current;
    }
    Some(result)
}

fn triangle(tri: [[f64; 3]; 3], points: &[[f64; 2]], ring: &[usize], depth: f64) -> Option<bool> {
    let mut poly = clip(&tri, |p| p[2])?;
    poly = clip(&poly, |p| depth - p[2])?;
    for (a, b) in planar::edges(ring) {
        if poly.is_empty() {
            return Some(false);
        }
        let [a, b] = [points[a], points[b]];
        let length = (b[0] - a[0]).hypot(b[1] - a[1]);
        if length <= 0.0 || !length.is_finite() {
            return None;
        }
        poly = clip(&poly, |p| planar::turn(a, b, [p[0], p[1]]) / length)?;
    }
    Some(!poly.is_empty())
}

pub(super) fn intersects(
    coords: &[[f64; 3]],
    polygon: &[usize],
    points: &[[f64; 2]],
    ring: &[usize],
    depth: f64,
) -> bool {
    let Some((&first, rest)) = polygon.split_first() else {
        return true;
    };
    if polygon.iter().all(|&id| coords[id][2] < -KERNEL_EPSILON)
        || polygon
            .iter()
            .all(|&id| coords[id][2] > depth + KERNEL_EPSILON)
    {
        return false;
    }
    if (0..2).any(|axis| {
        let lo = ring
            .iter()
            .map(|&id| points[id][axis])
            .fold(f64::INFINITY, f64::min);
        let hi = ring
            .iter()
            .map(|&id| points[id][axis])
            .fold(f64::NEG_INFINITY, f64::max);
        polygon
            .iter()
            .all(|&id| coords[id][axis] < lo - KERNEL_EPSILON)
            || polygon
                .iter()
                .all(|&id| coords[id][axis] > hi + KERNEL_EPSILON)
    }) {
        return false;
    }
    rest.windows(2).any(|p| {
        triangle(
            [coords[first], coords[p[0]], coords[p[1]]],
            points,
            ring,
            depth,
        ) != Some(false)
    })
}
