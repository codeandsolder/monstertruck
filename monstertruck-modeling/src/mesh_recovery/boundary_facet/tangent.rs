//! Sample exact tangent directions of positive-weight rational Bezier spans.
//! Used only to propose fan parameters; continuous certification follows.

use super::{Point3, Vector3, Vector4, rational_point, subdivide_homogeneous_bezier};

pub(super) fn visit_span_tangents(
    controls: &[[f64; 4]],
    depth: usize,
    remaining: &mut usize,
    visit: &mut impl FnMut(Point3, Vector3) -> bool,
) -> Option<bool> {
    *remaining = remaining.checked_sub(1)?;
    if depth > 0 {
        let (left, right) = subdivide_homogeneous_bezier(controls)?;
        if !visit_span_tangents(&left, depth - 1, remaining, visit)? {
            return Some(false);
        }
        return visit_span_tangents(&right, depth - 1, remaining, visit);
    }
    let point = |[x, y, z, w]: [f64; 4]| rational_point(Vector4::new(x, y, z, w));
    let [a, b, ..] = controls else {
        return None;
    };
    let [.., c, d] = controls else {
        return None;
    };
    let (a, b, c, d) = (point(*a), point(*b), point(*c), point(*d));
    // Positive endpoint weights preserve the direction of the adjacent
    // Euclidean control segment in the exact rational derivative.
    Some(visit(a, b - a) && visit(d, d - c))
}
