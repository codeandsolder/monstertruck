//! Structural identity of rational ruled-extrusion support surfaces.
//!
//! This is an algebraic control-net certificate, NOT a sampled-surface
//! proximity test. Identical degree/knots/normalized profile weights and
//! control points, together with parallel constant linear sweep rows, prove
//! the same infinite ruled carrier under an affine sweep-coordinate change.
//! A positive common sweep interval is independently required below.

use monstertruck_core::{
    cgmath64::{InnerSpace, Point3},
    tolerance::DEFAULT_NUMERICAL_TOLERANCES,
};
use monstertruck_geometry::prelude::{BsplineSurface, Vector4};

fn project(point: Vector4) -> Option<Point3> {
    (point.w.is_finite() && point.w > 0.0)
        .then(|| Point3::new(point.x / point.w, point.y / point.w, point.z / point.w))
}

/// Certify a positive-area shared carrier between two complete rational
/// ruled-extrusion patches (same profile, independently trimmed sweep spans).
/// This does not prove a boolean face relation, wire trim membership or
/// orientation; downstream topology partition and classification still apply.
pub(super) fn common_ruled_extrusion_strip(
    a: &BsplineSurface<Vector4>,
    b: &BsplineSurface<Vector4>,
) -> bool {
    let epsilon = DEFAULT_NUMERICAL_TOLERANCES.same_domain_similarity;
    let (u_degree_a, v_degree_a) = a.degrees();
    let (u_degree_b, v_degree_b) = b.degrees();
    if u_degree_a != u_degree_b || v_degree_a != 1 || v_degree_b != 1 {
        return false;
    }
    let net_a = a.control_points();
    let net_b = b.control_points();
    if net_a.len() < 2
        || net_a.len() != net_b.len()
        || net_a.iter().any(|row| row.len() != 2)
        || net_b.iter().any(|row| row.len() != 2)
    {
        return false;
    }
    let (u_a, v_a) = a.knot_vectors();
    let (u_b, v_b) = b.knot_vectors();
    if u_a.len() != u_b.len()
        || v_a.len() != 4
        || v_b.len() != 4
        || u_a
            .iter()
            .zip(u_b.iter())
            .any(|(a, b)| (a - b).abs() > epsilon)
        || (v_a[0] - v_a[1]).abs() > epsilon
        || (v_a[2] - v_a[3]).abs() > epsilon
        || (v_b[0] - v_b[1]).abs() > epsilon
        || (v_b[2] - v_b[3]).abs() > epsilon
    {
        return false;
    }

    let Some(origin_a) = project(net_a[0][0]) else {
        return false;
    };
    let Some(end_a) = project(net_a[0][1]) else {
        return false;
    };
    let Some(origin_b) = project(net_b[0][0]) else {
        return false;
    };
    let Some(end_b) = project(net_b[0][1]) else {
        return false;
    };
    let sweep_a = end_a - origin_a;
    let sweep_b = end_b - origin_b;
    if sweep_a.magnitude() <= epsilon || sweep_b.magnitude() <= epsilon {
        return false;
    }
    let axis = sweep_a.normalize();
    if sweep_a.cross(sweep_b).magnitude() > epsilon
        || (origin_b - origin_a).cross(axis).magnitude() > epsilon
    {
        return false;
    }
    let offset = origin_b - origin_a;
    let mut profile_is_nondegenerate = false;
    let weight_a0 = net_a[0][0].w;
    let weight_b0 = net_b[0][0].w;
    for (row_a, row_b) in net_a.iter().zip(net_b) {
        let Some(aa) = project(row_a[0]) else {
            return false;
        };
        let Some(ab) = project(row_a[1]) else {
            return false;
        };
        let Some(ba) = project(row_b[0]) else {
            return false;
        };
        let Some(bb) = project(row_b[1]) else {
            return false;
        };
        if (ab - aa - sweep_a).magnitude() > epsilon
            || (bb - ba - sweep_b).magnitude() > epsilon
            || (ba - aa - offset).magnitude() > epsilon
        {
            return false;
        }
        for (p, q) in row_a.iter().zip(row_b.iter()) {
            if (p.w / weight_a0 - q.w / weight_b0).abs() > epsilon {
                return false;
            }
        }
        if (row_a[0].w - row_a[1].w).abs() > epsilon * weight_a0
            || (row_b[0].w - row_b[1].w).abs() > epsilon * weight_b0
        {
            return false;
        }
        profile_is_nondegenerate |= (aa - origin_a).cross(axis).magnitude() > epsilon;
    }
    if !profile_is_nondegenerate {
        return false;
    }
    // The complete two-control-point sweeps each describe an exact finite
    // axial interval; unlike a broad-phase bbox this is a support certificate.
    let offset_on_axis = offset.dot(axis);
    let sweep_b_on_axis = sweep_b.dot(axis);
    let (b0, b1) = (
        offset_on_axis.min(offset_on_axis + sweep_b_on_axis),
        offset_on_axis.max(offset_on_axis + sweep_b_on_axis),
    );
    0.0f64.min(sweep_a.magnitude()).max(b0) + epsilon < sweep_a.magnitude().max(0.0).min(b1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::KnotVector;

    fn extruded(
        xs: (f64, f64),
        ydelta: f64,
        changed_weight: bool,
        skew: bool,
    ) -> BsplineSurface<Vector4> {
        let profiles = [(0.0, 0.0, 1.0), (0.35, 0.6, 0.8), (1.0, 1.0, 1.0)];
        let net = profiles
            .into_iter()
            .enumerate()
            .map(|(i, (y, z, w))| {
                let w0 = if changed_weight && i == 1 {
                    w + 0.025
                } else {
                    w
                };
                let drift = if skew && i == 1 { 0.01 } else { 0.0 };
                vec![
                    Vector4::new(xs.0 * w0, (y + ydelta) * w0, z * w0, w0),
                    Vector4::new(xs.1 * w0, (y + ydelta + drift) * w0, z * w0, w0),
                ]
            })
            .collect();
        BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
            net,
        )
    }

    #[test]
    fn nested_rational_ruled_profile_with_different_sweep_span_is_same_carrier() {
        let outer = extruded((0.0, 2.0), 0.0, false, false);
        let inner = extruded((0.5, 1.0), 0.0, false, false);
        assert!(common_ruled_extrusion_strip(&outer, &inner));
        assert!(common_ruled_extrusion_strip(&inner, &outer));
        assert!(common_ruled_extrusion_strip(
            &outer,
            &extruded((1.0, 0.5), 0.0, false, false)
        ));
    }
    #[test]
    fn ruled_carrier_refuses_different_profiles_twist_or_disjoint_spans() {
        let outer = extruded((0.0, 2.0), 0.0, false, false);
        assert!(!common_ruled_extrusion_strip(
            &outer,
            &extruded((0.5, 1.0), 0.01, false, false)
        ));
        assert!(!common_ruled_extrusion_strip(
            &outer,
            &extruded((0.5, 1.0), 0.0, true, false)
        ));
        assert!(!common_ruled_extrusion_strip(
            &outer,
            &extruded((0.5, 1.0), 0.0, false, true)
        ));
        assert!(!common_ruled_extrusion_strip(
            &outer,
            &extruded((2.01, 3.0), 0.0, false, false)
        ));
    }
}
