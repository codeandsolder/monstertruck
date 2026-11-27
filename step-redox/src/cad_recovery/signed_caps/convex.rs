//! Geometrically equivalent convex prism descriptions and exterior openings.
//! Moving a cutter plane outside a supporting halfspace of the host changes
//! no material inside that host. The complete boundary proof still applies.

use crate::{
    math2, math3, profile_curves::RecoveredProfileCurve, solid_extrusions::RecoveredSolidExtrusion,
};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;

const TOL: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIR: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

type Plane2 = ([f64; 2], f64);

fn polygon(prism: &RecoveredSolidExtrusion) -> Option<(Vec<[f64; 2]>, Vec<Plane2>)> {
    if !prism.inner_profile_loops.is_empty() || prism.profile_curves.len() < 3 {
        return None;
    }
    let points = prism
        .profile_curves
        .iter()
        .map(|curve| match curve {
            RecoveredProfileCurve::Line { start_mm, .. } => Some(*start_mm),
            _ => None,
        })
        .collect::<Option<Vec<_>>>()?;
    let planes = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(&a, &b)| {
            let [x, y] = math2::normalize(math2::sub(b, a), TOL)?;
            let normal = [y, -x];
            Some((normal, math2::dot(normal, a)))
        })
        .collect::<Option<Vec<_>>>()?;
    // Canonical outer loops are CCW. Refuse nonconvex, degenerate and
    // self-intersecting profiles rather than treating their vertices as a hull.
    if planes
        .iter()
        .any(|&(normal, offset)| points.iter().any(|&p| math2::dot(normal, p) > offset + TOL))
        || planes
            .iter()
            .zip(planes.iter().cycle().skip(1))
            .any(|(a, b)| math2::cross(a.0, b.0) <= DIR)
    {
        return None;
    }
    Some((points, planes))
}

fn world(prism: &RecoveredSolidExtrusion, p: [f64; 2], station: f64) -> [f64; 3] {
    use math3::{add, mul};
    add(
        add(
            add(prism.origin_mm, mul(prism.x_axis, p[0])),
            mul(prism.y_axis, p[1]),
        ),
        mul(prism.z_axis, station),
    )
}

fn vertices(prism: &RecoveredSolidExtrusion) -> Option<Vec<[f64; 3]>> {
    let (points, _) = polygon(prism)?;
    Some(
        points
            .into_iter()
            .flat_map(|p| [world(prism, p, 0.0), world(prism, p, prism.height_mm)])
            .collect(),
    )
}

/// Convex hulls with a bijective vertex correspondence describe the same
/// material. Nonconvex or curved shapes are not merged by this shortcut.
pub(super) fn equivalent(a: &RecoveredSolidExtrusion, b: &RecoveredSolidExtrusion) -> bool {
    let (Some(a), Some(mut b)) = (vertices(a), vertices(b)) else {
        return false;
    };
    if a.len() != b.len() {
        return false;
    }
    for p in a {
        let mut matching = b
            .iter()
            .enumerate()
            .filter(|(_, q)| math3::distance(p, **q) <= TOL);
        let Some((i, _)) = matching.next() else {
            return false;
        };
        if matching.next().is_some() {
            return false;
        }
        b.swap_remove(i);
    }
    true
}

fn supporting(host: &[[f64; 3]], normal: [f64; 3], point: [f64; 3]) -> bool {
    let maximum = host
        .iter()
        .map(|&v| math3::dot(normal, math3::sub(v, point)))
        .fold(f64::NEG_INFINITY, f64::max);
    maximum.is_finite() && maximum.abs() <= TOL
}

/// Extend only planes coinciding with an exterior supporting plane.
/// All other cutter planes remain fixed. The positive clearance scales
/// with the cutter rather than with a modeling tolerance.
pub(super) fn through_openings(
    host: &RecoveredSolidExtrusion,
    cutter: &RecoveredSolidExtrusion,
) -> Option<RecoveredSolidExtrusion> {
    use math3::{add, mul};
    let host = vertices(host)?;
    let (points, mut planes) = polygon(cutter)?;
    let clearance = points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .map(|(&a, &b)| math2::distance(a, b))
        .fold(cutter.height_mm, f64::min);
    if !clearance.is_finite() || clearance <= TOL {
        return None;
    }
    let mut changed = false;
    for (&point, (normal, offset)) in points.iter().zip(&mut planes) {
        let outward = add(mul(cutter.x_axis, normal[0]), mul(cutter.y_axis, normal[1]));
        if supporting(&host, outward, world(cutter, point, 0.0)) {
            *offset += clearance;
            changed = true;
        }
    }
    let mut result = cutter.clone();
    if supporting(&host, mul(cutter.z_axis, -1.0), cutter.origin_mm) {
        result.origin_mm = add(cutter.origin_mm, mul(cutter.z_axis, -clearance));
        result.height_mm += clearance;
        changed = true;
    }
    if supporting(
        &host,
        cutter.z_axis,
        world(cutter, [0.0, 0.0], cutter.height_mm),
    ) {
        result.height_mm += clearance;
        changed = true;
    }
    if !changed {
        return None;
    }
    let intersection = |(a, u): Plane2, (b, v): Plane2| {
        let det = math2::cross(a, b);
        let p = [
            u.mul_add(b[1], -a[1] * v) / det,
            a[0].mul_add(v, -u * b[0]) / det,
        ];
        p.into_iter().all(f64::is_finite).then_some(p)
    };
    let expanded = planes
        .iter()
        .cycle()
        .skip(planes.len() - 1)
        .zip(&planes)
        .map(|(&a, &b)| intersection(a, b))
        .collect::<Option<Vec<_>>>()?;
    // Reject if offsetting planes made any original non-opening constraint
    // inactive or changed the convex polygon's combinatorics.
    if expanded.iter().any(|&p| {
        planes
            .iter()
            .any(|&(normal, offset)| math2::dot(normal, p) > offset + TOL)
    }) {
        return None;
    }
    for ((curve, &start), &end) in result
        .profile_curves
        .iter_mut()
        .zip(&expanded)
        .zip(expanded.iter().cycle().skip(1))
    {
        let RecoveredProfileCurve::Line {
            start_mm, end_mm, ..
        } = curve
        else {
            return None;
        };
        *start_mm = start;
        *end_mm = end;
    }
    polygon(&result)?;
    Some(result)
}
