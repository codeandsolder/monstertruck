//! Common piecewise cubic native representation of a certified analytic
//! intersection branch. The caller proves every entire cubic on both
//! carriers and selects the continuous branch; no interpolant is admitted
//! merely because it matches its endpoints.

use crate::{BsplineCurve, Curve, KnotVector, Point3, Vector3};

// Keep the reconciled topology vertex as the exact spline endpoint. The
// adjacent control uses the analytic branch tangent, so pinning does not
// change the derivative at that endpoint. The caller must certify the
// modified entire cubic against both carriers and its selected branch.
pub(super) fn anchored_hermite(
    a: (Point3, Vector3),
    b: (Point3, Vector3),
    step: f64,
    anchors: [Option<Point3>; 2],
) -> [Point3; 4] {
    let front = anchors[0].unwrap_or(a.0);
    let back = anchors[1].unwrap_or(b.0);
    [
        front,
        front + a.1 * (step / 3.0),
        back - b.1 * (step / 3.0),
        back,
    ]
}

pub(super) fn section_curve(
    start: f64,
    end: f64,
    endpoints: [Point3; 2],
    cubic: impl Fn(f64, f64, [Option<Point3>; 2]) -> Option<[Point3; 4]>,
) -> Option<Curve> {
    let mut pieces = Vec::new();
    append_cubics(start, end, endpoints.map(Some), &cubic, 8, &mut pieces)?;
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
    anchors: [Option<Point3>; 2],
    cubic: &impl Fn(f64, f64, [Option<Point3>; 2]) -> Option<[Point3; 4]>,
    depth: usize,
    pieces: &mut Vec<(f64, [Point3; 4])>,
) -> Option<()> {
    if let Some(controls) = cubic(start, end, anchors) {
        pieces.push((end, controls));
        return Some(());
    }
    if depth == 0 {
        return None;
    }
    let mid = start.midpoint(end);
    append_cubics(start, mid, [anchors[0], None], cubic, depth - 1, pieces)?;
    append_cubics(mid, end, [None, anchors[1]], cubic, depth - 1, pieces)
}
