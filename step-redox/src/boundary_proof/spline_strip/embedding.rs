//! Injectivity of a positive-weight Bernstein profile transverse to a sweep.
//! A span lies in its complete control hull. Separation certificates are
//! therefore whole-domain proofs, not polyline or sampled-point tests.
use super::{SplineCurveKey, POINT_Q_MM, decode};
use crate::math3::{cross, dot, mul, normalize, sub};

type Point = [f64; 3];

fn transverse_direction(direction: Point, axis: Point) -> Option<Point> {
    normalize(sub(direction, mul(axis, dot(axis, direction))), POINT_Q_MM)
}

/// Candidate separating axes from the two projected control hulls.
/// Every convex-polygon edge is a pair of control points. Include chord
/// directions as well for degenerate line hulls.
fn directions(points: &[Point], axis: Point) -> impl Iterator<Item = Point> + '_ {
    points.iter().enumerate().flat_map(move |(i, &a)| {
        points[i + 1..].iter().flat_map(move |&b| {
            let delta = sub(b, a);
            [delta, cross(axis, delta)]
                .into_iter()
                .filter_map(move |d| transverse_direction(d, axis))
        })
    })
}

fn monotone_span(points: &[Point], axis: Point) -> bool {
    let (Some(first), Some(last)) = (points.first(), points.last()) else {
        return false;
    };
    let chord = sub(*last, *first);
    let monotone_on = |d: Point| {
        let d = if dot(chord, d) < 0. { mul(d, -1.) } else { d };
        dot(chord, d) > POINT_Q_MM
            && points
                .windows(2)
                .all(|pair| dot(sub(pair[1], pair[0]), d) >= 0.)
    };
    transverse_direction(chord, axis).is_some_and(monotone_on)
        || directions(points, axis).any(monotone_on)
}

fn projected_range(points: &[Point], origin: Point, d: Point) -> (f64, f64) {
    points
        .iter()
        .map(|&point| dot(sub(point, origin), d))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), value| {
            (lo.min(value), hi.max(value))
        })
}

fn separated(left: &[Point], right: &[Point], adjacent: bool, axis: Point) -> bool {
    let (Some(&origin), Some(&next)) = (left.last(), right.first()) else {
        return false;
    };
    if adjacent && origin.map(f64::to_bits) != next.map(f64::to_bits) {
        return false;
    }
    let separates = |d| {
        let (a_lo, a_hi) = projected_range(left, origin, d);
        let (b_lo, b_hi) = projected_range(right, origin, d);
        if adjacent {
            // The common endpoint is exactly on the plane. At least one
            // strictly sided control on each span guarantees every interior
            // parameter is strictly sided by positive Bernstein weights.
            (a_hi <= 0. && b_lo >= 0. && a_lo < -POINT_Q_MM && b_hi > POINT_Q_MM)
                || (b_hi <= 0. && a_lo >= 0. && b_lo < -POINT_Q_MM && a_hi > POINT_Q_MM)
        } else {
            a_hi + POINT_Q_MM < b_lo || b_hi + POINT_Q_MM < a_lo
        }
    };
    directions(left, axis)
        .chain(directions(right, axis))
        .chain(transverse_direction(sub(next, origin), axis))
        .any(separates)
}

/// Prove that projection along the extrusion axis cannot identify distinct
/// curve parameters. Each span is strictly monotone on some transverse
/// coordinate; adjacent hulls meet only on a plane through their common
/// endpoint, and all nonadjacent hulls are strictly separated.
pub(super) fn is_embedded(key: &SplineCurveKey, axis: Point) -> bool {
    let Some(spans) = key
        .spans
        .iter()
        .map(|span| {
            span.points_mm
                .iter()
                .copied()
                .map(decode)
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()
    else {
        return false;
    };
    !spans.is_empty()
        && spans.iter().all(|span| monotone_span(span, axis))
        && spans.iter().enumerate().all(|(i, left)| {
            spans
                .iter()
                .enumerate()
                .skip(i + 1)
                .all(|(j, right)| separated(left, right, j == i + 1, axis))
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn turning_spans_are_separated_without_a_global_monotone_direction() {
        let left = [[0., 0., 0.], [1., -0.001, 0.], [2., 1., 0.], [2., 2., 0.]];
        let right = [[2., 2., 0.], [2., 3., 0.], [1., 4.001, 0.], [0., 4., 0.]];
        let axis = [0., 0., 1.];
        assert!(monotone_span(&left, axis));
        assert!(monotone_span(&right, axis));
        assert!(separated(&left, &right, true, axis));
        // A return across an earlier hull cannot establish injectivity.
        let crossing = [[2., 2., 0.], [1., 1., 0.], [0., 0., 0.]];
        assert!(!separated(&left, &crossing, true, axis));
        assert!(!separated(&left, &crossing, false, axis));
        // Coincidence after projection must also fail, despite 3D separation.
        let lifted = left.map(|[x, y, z]| [x, y, z + 1.]);
        assert!(!separated(&left, &lifted, false, axis));
    }
}
