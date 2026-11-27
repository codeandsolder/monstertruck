//! Source-cap separation proposals. An overlapping cover is not unique:
//! a supported end plane of one prism can cut redundant material from
//! another. Only the complete source-boundary proof admits such proposals.
use crate::{
    cad_ir::Profile2d,
    math2::{cross as cross2, distance as distance2, dot as dot2, norm as norm2, sub as sub2},
    math3::{add, cross, dot, mul, normalize, sub},
    profile_curves::RecoveredProfileCurve,
    solid_extrusions::{RecoveredSourceCapExtrusion, angle_on_arc, signed_area_curves},
};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;

const LENGTH: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_geometry_mm;
const DIRECTION: f64 = DEFAULT_STEP_RECOVERY_TOLERANCES.solid_direction;

#[derive(Clone, Copy)]
struct Plane {
    point: [f64; 3],
    normal: [f64; 3],
}

fn world(cap: &RecoveredSourceCapExtrusion, point: [f64; 2]) -> [f64; 3] {
    add(
        cap.origin_mm,
        add(mul(cap.x_axis, point[0]), mul(cap.y_axis, point[1])),
    )
}

fn curve_maximum(curve: &RecoveredProfileCurve, direction: [f64; 2]) -> Option<f64> {
    let endpoints = dot2(curve.start_point()?, direction).max(dot2(curve.end_point()?, direction));
    match curve {
        RecoveredProfileCurve::Line { .. } => Some(endpoints),
        RecoveredProfileCurve::CircleArc {
            center_mm,
            radius_mm,
            start_angle_rad,
            end_angle_rad,
            ..
        } => {
            let angle = direction[1].atan2(direction[0]);
            Some(
                if angle_on_arc(angle, *start_angle_rad, *end_angle_rad, DIRECTION) {
                    endpoints.max(radius_mm.mul_add(norm2(direction), dot2(*center_mm, direction)))
                } else {
                    endpoints
                },
            )
        }
        // Positive-weight control hulls bound even noncircular profiles.
        RecoveredProfileCurve::Bezier {
            control_points_mm, ..
        } => Some(
            control_points_mm
                .iter()
                .map(|&p| dot2(p, direction))
                .fold(endpoints, f64::max),
        ),
        RecoveredProfileCurve::BSpline {
            control_points_mm,
            weights,
            ..
        } if weights
            .as_ref()
            .is_none_or(|w| w.iter().all(|x| x.is_finite() && *x > 0.0)) =>
        {
            Some(
                control_points_mm
                    .iter()
                    .map(|&p| dot2(p, direction))
                    .fold(endpoints, f64::max),
            )
        }
        RecoveredProfileCurve::BSpline { .. } => None,
    }
}

fn support_height(cap: &RecoveredSourceCapExtrusion, plane: Plane) -> Option<f64> {
    let direction = [dot(cap.x_axis, plane.normal), dot(cap.y_axis, plane.normal)];
    let base = dot(sub(cap.origin_mm, plane.point), plane.normal);
    let axial = (dot(cap.z_axis, plane.normal) * cap.height_mm).max(0.0);
    let maximum = cap
        .profile_curves
        .iter()
        .chain(cap.inner_profile_loops.iter().flatten())
        .try_fold(f64::NEG_INFINITY, |maximum, curve| {
            Some(maximum.max(curve_maximum(curve, direction)?))
        })?;
    let height = base + maximum + axial;
    height.is_finite().then_some(height)
}

fn supporting_planes(cap: &RecoveredSourceCapExtrusion) -> Vec<Plane> {
    let mut planes = vec![
        Plane {
            point: cap.origin_mm,
            normal: mul(cap.z_axis, -1.0),
        },
        Plane {
            point: add(cap.origin_mm, mul(cap.z_axis, cap.height_mm)),
            normal: cap.z_axis,
        },
    ];
    let sign = signed_area_curves(&cap.profile_curves).signum();
    for curve in &cap.profile_curves {
        let RecoveredProfileCurve::Line {
            start_mm, end_mm, ..
        } = curve
        else {
            continue;
        };
        let point = world(cap, *start_mm);
        let direction = sub(world(cap, *end_mm), point);
        if let Some(normal) = normalize(mul(cross(direction, cap.z_axis), sign), LENGTH) {
            planes.push(Plane { point, normal });
        }
    }
    planes.retain(|&plane| support_height(cap, plane).is_some_and(|h| h <= LENGTH));
    planes
}

// A convex polygon remains one connected profile after a half-plane
// intersection. Check convexity against every edge rather than merely
// requiring consistent local turns, which would admit star polygons.
fn convex_profile(cap: &RecoveredSourceCapExtrusion) -> Option<Vec<[f64; 2]>> {
    if !cap.inner_profile_loops.is_empty() {
        return None;
    }
    let points = cap
        .profile_curves
        .iter()
        .map(|curve| {
            let RecoveredProfileCurve::Line { start_mm, .. } = curve else {
                return None;
            };
            Some(*start_mm)
        })
        .collect::<Option<Vec<_>>>()?;
    let sign = signed_area_curves(&cap.profile_curves).signum();
    if points.len() < 3 || !sign.is_finite() || sign.abs() < 0.5 {
        return None;
    }
    for (i, &a) in points.iter().enumerate() {
        let delta = sub2(points[(i + 1) % points.len()], a);
        let length = norm2(delta);
        if length <= LENGTH
            || points
                .iter()
                .any(|&p| sign * cross2(delta, sub2(p, a)) < -LENGTH * length)
        {
            return None;
        }
    }
    Some(points)
}

fn clip_polygon(points: &[[f64; 2]], normal: [f64; 2], offset: f64) -> Option<Vec<[f64; 2]>> {
    let heights = points
        .iter()
        .map(|&p| dot2(p, normal) + offset)
        .collect::<Vec<_>>();
    if !heights.iter().any(|&h| h > LENGTH) || !heights.iter().any(|&h| h < -LENGTH) {
        return None;
    }
    let mut clipped = Vec::new();
    let mut push = |point| {
        if clipped.last().is_none_or(|&p| distance2(p, point) > LENGTH) {
            clipped.push(point);
        }
    };
    for (i, &a) in points.iter().enumerate() {
        let j = (i + 1) % points.len();
        let b = points[j];
        let [ha, hb] = [heights[i], heights[j]];
        if ha >= 0.0 {
            push(a);
        }
        if (ha >= 0.0) != (hb >= 0.0) {
            let t = ha / (ha - hb);
            push([t.mul_add(b[0] - a[0], a[0]), t.mul_add(b[1] - a[1], a[1])]);
        }
    }
    if clipped.len() > 1 && distance2(clipped[0], *clipped.last()?) <= LENGTH {
        clipped.pop();
    }
    (clipped.len() >= 3).then_some(clipped)
}

fn separated_profile(cap: &RecoveredSourceCapExtrusion, plane: Plane) -> Option<Profile2d> {
    // Cutting a prism by a plane parallel to its generators is a profile
    // intersection, not a 3-D Boolean. The resulting profile is a proposal
    // in the original local frame; the kernel must still certify contact.
    if dot(cap.z_axis, plane.normal).abs() > DIRECTION {
        return None;
    }
    let normal = [dot(cap.x_axis, plane.normal), dot(cap.y_axis, plane.normal)];
    let offset = dot(sub(cap.origin_mm, plane.point), plane.normal);
    let points = convex_profile(cap)?;
    Profile2d::polygon(clip_polygon(&points, normal, offset)?).ok()
}

pub(super) fn proposals(caps: &[RecoveredSourceCapExtrusion]) -> Vec<(usize, Profile2d)> {
    let mut proposals = Vec::new();
    for (i, host) in caps.iter().enumerate() {
        for plane in supporting_planes(host) {
            for (j, target) in caps.iter().enumerate() {
                if i != j
                    && let Some(profile) = separated_profile(target, plane)
                    && !proposals
                        .iter()
                        .any(|(index, existing)| *index == j && *existing == profile)
                {
                    proposals.push((j, profile));
                }
            }
        }
    }
    proposals
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn clips_oblique_convex_profile_without_moving_retained_vertices() {
        let points = [[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]];
        let clipped = clip_polygon(&points, [1.0, 1.0], -2.0);
        assert_eq!(clipped, Some(vec![[2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]));
        assert!(clip_polygon(&points, [1.0, 0.0], 1.0).is_none());
    }
}
