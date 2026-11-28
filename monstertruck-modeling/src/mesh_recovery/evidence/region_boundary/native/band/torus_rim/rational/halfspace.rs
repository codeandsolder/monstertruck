//! Strict open-halfspace certificates for complete rational curve spans.
use super::{Control, Controls, Interval, MAX_DEPTH, MAX_PIECES, controls};
use crate::{Edge, Point3, Vector3};
use monstertruck_geometry::prelude::TryIntoHomogeneousBsplineCurve;

fn signed(control: Control, origin: Point3, normal: Vector3, margin: f64) -> Interval {
    let w = control[3];
    (0..3)
        .map(|i| control[i].sub(w.scale(origin[i])).scale(normal[i]))
        .fold(Interval::point(0.0), Interval::add)
        .sub(w.scale(margin))
}

pub(in crate::mesh_recovery) fn point_in_halfspace(
    point: Point3,
    origin: Point3,
    normal: Vector3,
    margin: f64,
) -> bool {
    let bound = signed(
        [point.x, point.y, point.z, 1.0].map(Interval::point),
        origin,
        normal,
        margin,
    );
    bound.finite() && bound.lo > 0.0
}

fn span_in_halfspace(
    span: &Controls,
    origin: Point3,
    normal: Vector3,
    margin: f64,
    depth: usize,
    remaining: &mut usize,
) -> bool {
    if *remaining == 0 {
        return false;
    }
    *remaining -= 1;
    let bounds = span
        .values
        .iter()
        .copied()
        .map(|control| signed(control, origin, normal, margin))
        .collect::<Vec<_>>();
    if bounds.iter().all(|bound| bound.finite() && bound.lo > 0.0) {
        return true;
    }
    if depth == 0
        || bounds.first().is_none_or(|b| !b.finite() || b.lo <= 0.0)
        || bounds.last().is_none_or(|b| !b.finite() || b.lo <= 0.0)
    {
        return false;
    }
    span.split()
        .iter()
        .all(|part| span_in_halfspace(part, origin, normal, margin, depth - 1, remaining))
}

pub(in crate::mesh_recovery) fn curve_in_halfspace(
    edge: &Edge,
    origin: Point3,
    normal: Vector3,
    margin: f64,
) -> bool {
    if ![edge.front().point(), edge.back().point()]
        .iter()
        .all(|&point| point_in_halfspace(point, origin, normal, margin))
    {
        return false;
    }
    let Some(spline) = edge.oriented_curve().try_into_homogeneous_bspline_curve() else {
        return false;
    };
    let Some(spans) = controls::spans(&spline) else {
        return false;
    };
    let mut remaining = MAX_PIECES;
    !spans.is_empty()
        && spans
            .iter()
            .all(|span| span_in_halfspace(span, origin, normal, margin, MAX_DEPTH, &mut remaining))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subdivision_proves_positive_curve_but_rejects_interior_exit() {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let normal = Vector3::new(1.0, 0.0, 0.0);
        let make = |middle: f64| Controls {
            values: [1.0, middle, 1.0]
                .into_iter()
                .map(|x| [x, 0.0, 0.0, 1.0].map(Interval::point))
                .collect(),
            parameter: [0.0, 1.0],
        };
        let mut remaining = MAX_PIECES;
        assert!(span_in_halfspace(
            &make(-0.2),
            origin,
            normal,
            0.1,
            MAX_DEPTH,
            &mut remaining
        ));
        let mut remaining = MAX_PIECES;
        assert!(!span_in_halfspace(
            &make(-2.0),
            origin,
            normal,
            0.1,
            MAX_DEPTH,
            &mut remaining
        ));
    }
}
