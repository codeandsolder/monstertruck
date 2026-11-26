//! Common piecewise cubic native representation of a certified analytic
//! intersection branch. The caller proves every entire cubic on both
//! carriers and selects the continuous branch; no interpolant is admitted
//! merely because it matches its endpoints.

use crate::{BsplineCurve, Curve, KnotVector, Point3};

pub(super) fn section_curve(
    start: f64,
    end: f64,
    cubic: impl Fn(f64, f64) -> Option<[Point3; 4]>,
) -> Option<Curve> {
    let mut pieces = Vec::new();
    append_cubics(start, end, &cubic, 8, &mut pieces)?;
    let mut controls = Vec::with_capacity(3 * pieces.len() + 1);
    let mut knots = vec![0.0; 4];
    controls.push(pieces.first()?.1[0]);
    let last = pieces.len() - 1;
    for (i, (angle, piece)) in pieces.into_iter().enumerate() {
        controls.extend_from_slice(&piece[1..]);
        if i == last {
            knots.extend([1.0; 4]);
        } else {
            knots.extend([(angle - start) / (end - start); 3]);
        }
    }
    Some(Curve::BsplineCurve(
        BsplineCurve::new(KnotVector::from(knots), controls).ok()?,
    ))
}

fn append_cubics(
    start: f64,
    end: f64,
    cubic: &impl Fn(f64, f64) -> Option<[Point3; 4]>,
    depth: usize,
    pieces: &mut Vec<(f64, [Point3; 4])>,
) -> Option<()> {
    if let Some(controls) = cubic(start, end) {
        pieces.push((end, controls));
        return Some(());
    }
    if depth == 0 {
        return None;
    }
    let mid = start.midpoint(end);
    append_cubics(start, mid, cubic, depth - 1, pieces)?;
    append_cubics(mid, end, cubic, depth - 1, pieces)
}
