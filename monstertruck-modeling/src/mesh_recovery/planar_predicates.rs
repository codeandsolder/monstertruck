//! Shared exact-orientation predicates for polygonal planar source evidence.
//! Use robust predicates for the actual projected coordinates; no geometric
//! tolerance is hidden in winding, crossing or orientation decisions.
use robust::{Coord, orient2d};

type Point = Coord<f64>;

pub(super) fn in_bounds(a: Point, b: Point, p: Point) -> bool {
    p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
}

pub(super) fn on_segment(a: Point, b: Point, p: Point) -> bool {
    orient2d(a, b, p) == 0.0 && in_bounds(a, b, p)
}

fn opposite(a: f64, b: f64) -> bool { (a > 0.0 && b < 0.0) || (a < 0.0 && b > 0.0) }

fn turns(a: Point, b: Point, c: Point, d: Point) -> [f64; 4] {
    [
        orient2d(a, b, c),
        orient2d(a, b, d),
        orient2d(c, d, a),
        orient2d(c, d, b),
    ]
}

pub(super) fn proper_crossing(a: Point, b: Point, c: Point, d: Point) -> bool {
    let t = turns(a, b, c, d);
    t.iter().any(|v| !v.is_finite()) || (opposite(t[0], t[1]) && opposite(t[2], t[3]))
}

pub(super) fn intersects(a: Point, b: Point, c: Point, d: Point) -> bool {
    let t = turns(a, b, c, d);
    t.iter().any(|v| !v.is_finite())
        || (opposite(t[0], t[1]) && opposite(t[2], t[3]))
        || (t[0] == 0.0 && in_bounds(a, b, c))
        || (t[1] == 0.0 && in_bounds(a, b, d))
        || (t[2] == 0.0 && in_bounds(c, d, a))
        || (t[3] == 0.0 && in_bounds(c, d, b))
}

/// Exact signed orientation at a lexicographic extreme vertex of a simple
/// polygon. Straight-through edges survive; retracing and intersections do
/// not. The caller supplies the coordinate projection, avoiding a copy of
/// indexed geometry just to inspect its planar boundary.
pub(super) fn polygon_sign<T: Copy>(points: &[T], project: impl Fn(T) -> Point) -> Option<bool> {
    let count = points.len();
    if count < 3 {
        return None;
    }
    for (i, &id) in points.iter().enumerate() {
        let a = project(id);
        let b = project(points[(i + 1) % count]);
        let c = project(points[(i + 2) % count]);
        let t = orient2d(a, b, c);
        if ![a.x, a.y, b.x, b.y, c.x, c.y, t]
            .into_iter()
            .all(f64::is_finite)
            || ((a.x - b.x).abs() == 0.0 && (a.y - b.y).abs() == 0.0)
            || (t == 0.0 && (in_bounds(a, b, c) || in_bounds(b, c, a)))
        {
            return None;
        }
        for (j, &other) in points.iter().enumerate().skip(i + 2) {
            if i == 0 && j + 1 == count {
                continue;
            }
            if intersects(a, b, project(other), project(points[(j + 1) % count])) {
                return None;
            }
        }
    }
    let extreme = points
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let a = project(**a);
            let b = project(**b);
            a.x.partial_cmp(&b.x)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
        })?
        .0;
    let t = orient2d(
        project(points[(extreme + count - 1) % count]),
        project(points[extreme]),
        project(points[(extreme + 1) % count]),
    );
    (t.is_finite() && t != 0.0).then_some(t > 0.0)
}
