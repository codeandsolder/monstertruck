//! Structural identity of rational ruled-extrusion support surfaces.
//!
//! This is an algebraic control-net certificate, NOT a sampled-surface
//! proximity test. Identical degree/knots/normalized profile weights and
//! control points, together with parallel constant linear sweep rows, prove
//! the same infinite ruled carrier under an affine sweep-coordinate change.
//! A positive common sweep interval is independently required below.

use monstertruck_core::{
    cgmath64::{EuclideanSpace, InnerSpace, Point3, Vector3},
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
pub(in crate::transversal::classic) fn common_ruled_extrusion_strip(
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

/// Algebraically certify a circular-cylinder carrier represented by a
/// quadratic rational B-spline profile with a linear extruded second axis.
///
/// Each nondegenerate knot span is converted to a quadratic Bezier span.
/// The quadratic circle equation in homogeneous coordinates then has degree
/// four: checking all five *polynomial coefficients* proves the entire span
/// lies on the circle, not merely finitely many sampled profile points.
pub(in crate::transversal::classic) fn certified_circular_cylinder(
    surface: &BsplineSurface<Vector4>,
) -> Option<(Point3, Vector3, f64, (f64, f64))> {
    use monstertruck_geometry::prelude::{BsplineCurve, Cut, ParametricCurve};
    let policy = DEFAULT_NUMERICAL_TOLERANCES;
    let spatial_eps = policy.boolean_cylinder_structural_length;
    let weight_eps = policy.boolean_cylinder_homogeneous_weight;
    let coefficient_eps = policy.boolean_cylinder_polynomial_relative;
    let (degree_u, degree_v) = surface.degrees();
    let net = surface.control_points();
    if degree_u != 2 || degree_v != 1 || net.len() < 3 || net.iter().any(|row| row.len() != 2) {
        return None;
    }
    let origin = project(net[0][0])?;
    let end = project(net[0][1])?;
    let sweep = end - origin;
    if sweep.magnitude() <= spatial_eps {
        return None;
    }
    let axis = sweep.normalize();
    for row in net {
        let start = project(row[0])?;
        let finish = project(row[1])?;
        if (finish - start - sweep).magnitude() > spatial_eps
            || (row[0].w - row[1].w).abs() > weight_eps
        {
            return None;
        }
    }

    let reference = if axis.x.abs() < policy.boolean_cylinder_reference_axis_component_limit {
        Vector3::unit_x()
    } else {
        Vector3::unit_y()
    };
    let basis_u = axis.cross(reference).normalize();
    let basis_v = axis.cross(basis_u).normalize();
    let profile = BsplineCurve::new_unchecked(
        surface.knot_vector_u().clone(),
        net.iter().map(|row| row[0]).collect(),
    );
    let knots = profile.knot_vector().clone();
    let lo = *knots.first()?;
    let hi = *knots.last()?;
    let xy = |h: Vector4| -> Option<(f64, f64)> {
        let p = project(h)?;
        let delta = p - origin;
        Some((delta.dot(basis_u), delta.dot(basis_v)))
    };
    let p0 = xy(profile.evaluate(lo))?;
    let p1 = xy(profile.evaluate(f64::midpoint(lo, hi)))?;
    let p2 = xy(profile.evaluate(hi))?;
    let twice_det = 2.0 * ((p1.0 - p0.0) * (p2.1 - p0.1) - (p2.0 - p0.0) * (p1.1 - p0.1));
    if twice_det.abs() <= spatial_eps * spatial_eps {
        return None;
    }
    let d1 = p1.0.mul_add(p1.0, p1.1 * p1.1) - p0.0.mul_add(p0.0, p0.1 * p0.1);
    let d2 = p2.0.mul_add(p2.0, p2.1 * p2.1) - p0.0.mul_add(p0.0, p0.1 * p0.1);
    let cx = (d1 * (p2.1 - p0.1) - d2 * (p1.1 - p0.1)) / twice_det;
    let cy = ((p1.0 - p0.0) * d2 - (p2.0 - p0.0) * d1) / twice_det;
    let radius2 = (p0.0 - cx).mul_add(p0.0 - cx, (p0.1 - cy).powi(2));
    if !radius2.is_finite() || radius2 <= spatial_eps * spatial_eps {
        return None;
    }

    // Every degree-two homogeneous Bezier span must satisfy
    // (X - cx W)^2 + (Y - cy W)^2 - r² W² == 0 identically.
    // Convert this product to power-basis coefficients, not to a collection of
    // point samples. Nonzero knot intervals are split by exact knot insertion.
    let certify = |segment: &BsplineCurve<Vector4>| -> Option<()> {
        let cp = segment.control_points();
        if cp.len() != 3 {
            return None;
        }
        let mut x = [0.0; 3];
        let mut y = [0.0; 3];
        let mut w = [0.0; 3];
        for (i, h) in cp.iter().copied().enumerate() {
            if !h.w.is_finite() || h.w <= 0.0 {
                return None;
            }
            let relative = Vector3::new(h.x, h.y, h.z) - origin.to_vec() * h.w;
            x[i] = relative.dot(basis_u) - cx * h.w;
            y[i] = relative.dot(basis_v) - cy * h.w;
            w[i] = h.w;
        }
        let power = |c: [f64; 3]| [c[0], 2.0 * (c[1] - c[0]), c[0] - 2.0 * c[1] + c[2]];
        let x = power(x);
        let y = power(y);
        let w = power(w);
        let mut residual = [0.0f64; 5];
        let mut coefficient_scale = [0.0f64; 5];
        for i in 0..3 {
            for j in 0..3 {
                let xx = x[i] * x[j];
                let yy = y[i] * y[j];
                let rr = radius2 * w[i] * w[j];
                residual[i + j] += xx + yy - rr;
                coefficient_scale[i + j] += xx.abs() + yy.abs() + rr.abs();
            }
        }
        // All five coefficients come from the same polynomial identity.
        // Zero or near-zero coefficients still accumulate roundoff from
        // the entire rational circle fit, not just their own cancellation.
        // Normalize by the maximum *actual squared-length* coefficient
        // contribution, keeping one dimensionally consistent scale.
        let scale = coefficient_scale.into_iter().fold(0.0_f64, f64::max);
        residual
            .iter()
            .all(|value| value.abs() <= coefficient_eps * scale)
            .then_some(())
    };

    let mut remaining = profile;
    let knot_slack = policy.boolean_cylinder_knot_relative * (hi - lo).abs()
        + policy.boolean_iso_parameter_roundoff_ulps * f64::EPSILON * lo.abs().max(hi.abs());
    let mut interior = knots
        .iter()
        .copied()
        .filter(|&k| k > lo + knot_slack && k < hi - knot_slack)
        .collect::<Vec<_>>();
    interior.dedup_by(|a, b| (*a - *b).abs() <= knot_slack);
    for knot in interior {
        let next = remaining.cut(knot);
        certify(&remaining)?;
        remaining = next;
    }
    certify(&remaining)?;
    let center = origin + basis_u * cx + basis_v * cy;
    let axial = (origin.to_vec().dot(axis), end.to_vec().dot(axis));
    Some((center, axis, radius2.sqrt(), axial))
}

/// Two distinct rational quadratic spline arcs can encode overlapping pieces
/// of the SAME exact circular cylinder even when their knot/control nets
/// differ. Refuse mere proximity: each curve must independently satisfy its
/// fitted circle polynomial on every knot span.
pub(in crate::transversal::classic) fn common_circular_extrusion_cylinder(
    a: &BsplineSurface<Vector4>,
    b: &BsplineSurface<Vector4>,
) -> bool {
    let (Some((ca, axis_a, ra, span_a)), Some((cb, axis_b, rb, span_b))) = (
        certified_circular_cylinder(a),
        certified_circular_cylinder(b),
    ) else {
        return false;
    };
    let spatial_eps = DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_structural_length;
    let axis_eps = DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_axis_direction;
    if axis_a.cross(axis_b).magnitude() > axis_eps
        || (cb - ca).cross(axis_a).magnitude() > spatial_eps
        || (ra - rb).abs() > spatial_eps
    {
        return false;
    }
    let a = (span_a.0.min(span_a.1), span_a.0.max(span_a.1));
    // Axial coordinates must share A's frame, even when B's sweep runs
    // in the opposite direction.
    let sign = if axis_a.dot(axis_b) >= 0.0 { 1.0 } else { -1.0 };
    let bx0 = sign * span_b.0;
    let bx1 = sign * span_b.1;
    let b = (bx0.min(bx1), bx0.max(bx1));
    a.0.max(b.0) + spatial_eps < a.1.min(b.1)
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

    fn rational_circular_arc(
        sweep: (f64, f64),
        angles: (f64, f64),
        radius: f64,
        center_y: f64,
        distortion: f64,
    ) -> BsplineSurface<Vector4> {
        let mid = f64::midpoint(angles.0, angles.1);
        let half = (angles.1 - angles.0) / 2.0;
        let w = half.cos();
        let points = [
            (angles.0.cos() * radius, angles.0.sin() * radius, 1.0),
            (
                mid.cos() * radius / w,
                mid.sin() * radius / w + distortion,
                w,
            ),
            (angles.1.cos() * radius, angles.1.sin() * radius, 1.0),
        ];
        let control = points
            .into_iter()
            .map(|(y, z, w)| {
                vec![
                    Vector4::new(sweep.0 * w, (center_y + y) * w, z * w, w),
                    Vector4::new(sweep.1 * w, (center_y + y) * w, z * w, w),
                ]
            })
            .collect();
        BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
            control,
        )
    }

    #[test]
    fn rational_circular_carrier_is_invariant_to_native_knot_reparameterization() {
        let original = rational_circular_arc(
            (0.0, 2.0),
            (0.0, std::f64::consts::FRAC_PI_2),
            0.6731,
            5.364572,
            0.0,
        );
        let mut u_knots = original.knot_vector_u().clone();
        u_knots.transform(1.0e-6, 1.0e-3);
        let reparameterized = BsplineSurface::new_unchecked(
            (u_knots, original.knot_vector_v().clone()),
            original.control_points().to_vec(),
        );
        assert!(certified_circular_cylinder(&original).is_some());
        assert!(certified_circular_cylinder(&reparameterized).is_some());
        assert!(common_circular_extrusion_cylinder(
            &original,
            &reparameterized
        ));
    }

    #[test]
    fn distinct_quadratic_rational_arcs_prove_same_exact_cylinder() {
        let outer = rational_circular_arc(
            (0.0, 2.0),
            (0.0, std::f64::consts::FRAC_PI_2),
            0.6731,
            5.364572,
            0.0,
        );
        let inner = rational_circular_arc(
            (0.5, 1.0),
            (0.0, std::f64::consts::FRAC_PI_4),
            0.6731,
            5.364572,
            0.0,
        );
        let other = rational_circular_arc(
            (0.5, 1.0),
            (std::f64::consts::FRAC_PI_4, std::f64::consts::FRAC_PI_2),
            0.6731,
            5.364572,
            0.0,
        );
        assert!(common_circular_extrusion_cylinder(&outer, &inner));
        assert!(common_circular_extrusion_cylinder(&inner, &outer));
        assert!(common_circular_extrusion_cylinder(&outer, &other));
        assert!(common_circular_extrusion_cylinder(
            &outer,
            &rational_circular_arc(
                (1.0, 0.5),
                (0.0, std::f64::consts::FRAC_PI_4),
                0.6731,
                5.364572,
                0.0
            ),
        ));
    }

    #[test]
    fn circular_carrier_rejects_nearby_but_different_analytic_cylinders() {
        let outer = rational_circular_arc(
            (0.0, 2.0),
            (0.0, std::f64::consts::FRAC_PI_2),
            0.6731,
            5.364572,
            0.0,
        );
        assert!(!common_circular_extrusion_cylinder(
            &outer,
            &rational_circular_arc(
                (0.5, 1.0),
                (0.0, std::f64::consts::FRAC_PI_4),
                0.67311,
                5.364572,
                0.0
            )
        ));
        assert!(!common_circular_extrusion_cylinder(
            &outer,
            &rational_circular_arc(
                (0.5, 1.0),
                (0.0, std::f64::consts::FRAC_PI_4),
                0.6731,
                5.364582,
                0.0
            )
        ));
        assert!(!common_circular_extrusion_cylinder(
            &outer,
            &rational_circular_arc(
                (0.5, 1.0),
                (0.0, std::f64::consts::FRAC_PI_4),
                0.6731,
                5.364572,
                0.001
            )
        ));
        assert!(!common_circular_extrusion_cylinder(
            &outer,
            &rational_circular_arc(
                (2.01, 3.0),
                (0.0, std::f64::consts::FRAC_PI_4),
                0.6731,
                5.364572,
                0.0
            )
        ));
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
