//! Structural identity of rational ruled-extrusion support surfaces.
//!
//! This is an algebraic control-net certificate, NOT a sampled-surface
//! proximity test. Identical degree/knots/normalized profile weights and
//! control points, together with parallel constant linear sweep rows, prove
//! the same infinite ruled carrier under an affine sweep-coordinate change.
//! A positive common sweep interval is independently required below.

use monstertruck_core::{math::InnerProduct, tolerance::DEFAULT_NUMERICAL_TOLERANCES};
pub(in crate::transversal::classic) use monstertruck_geometry::analytic_cylinder_intersection::certified_circular_cylinder;
use monstertruck_geometry::{
    analytic_cylinder_intersection::project_positive_homogeneous_point as project,
    prelude::{BsplineSurface, KnotVector, Vector4},
};

/// Compare the native knot positions in their respective affine parameter
/// domains. Equal normalized profiles describe the same spline even when a
/// STEP producer uses different scale and origin for its U/V parameters.
fn normalized_ruled_knots_match(a: &KnotVector, b: &KnotVector, relative: f64) -> bool {
    let (Some(&a0), Some(&a1), Some(&b0), Some(&b1)) = (a.first(), a.last(), b.first(), b.last())
    else {
        return false;
    };
    let span_a = a1 - a0;
    let span_b = b1 - b0;
    if a.len() != b.len()
        || !span_a.is_finite()
        || !span_b.is_finite()
        || span_a <= 0.0
        || span_b <= 0.0
    {
        return false;
    }
    a.iter().zip(b.iter()).all(|(&ka, &kb)| {
        let na = (ka - a0) / span_a;
        let nb = (kb - b0) / span_b;
        na.is_finite() && nb.is_finite() && (na - nb).abs() <= relative
    })
}

/// Prove that the complete homogeneous control nets represent congruent
/// ruled carriers with a positive common axial interval.
fn congruent_ruled_sweeps(
    net_a: &[Vec<Vector4>],
    net_b: &[Vec<Vector4>],
    spatial_tol: f64,
    direction_tol: f64,
    weight_tol: f64,
) -> bool {
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
    let length_a = sweep_a.norm();
    let length_b = sweep_b.norm();
    if !length_a.is_finite()
        || !length_b.is_finite()
        || length_a <= spatial_tol
        || length_b <= spatial_tol
    {
        return false;
    }
    let axis = sweep_a.normalize();
    // Both the angle and its model-space displacement along the complete
    // sweep must agree. A tiny angle over a long sweep can still separate
    // supporting carriers by much more than the permitted spatial residual.
    if axis.cross(&(sweep_b / length_b)).norm() > direction_tol
        || axis.cross(&sweep_b).norm() > spatial_tol
        || (origin_b - origin_a).cross(&axis).norm() > spatial_tol
    {
        return false;
    }
    let offset = origin_b - origin_a;
    let mut profile_is_nondegenerate = false;
    let reference_weight = net_a[0][0].w;
    let candidate_weight = net_b[0][0].w;
    if !reference_weight.is_finite()
        || !candidate_weight.is_finite()
        || reference_weight <= 0.0
        || candidate_weight <= 0.0
    {
        return false;
    }
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
        if (ab - aa - sweep_a).norm() > spatial_tol
            || (bb - ba - sweep_b).norm() > spatial_tol
            || (ba - aa - offset).norm() > spatial_tol
        {
            return false;
        }
        for (p, q) in row_a.iter().zip(row_b.iter()) {
            if (p.w / reference_weight - q.w / candidate_weight).abs() > weight_tol {
                return false;
            }
        }
        if (row_a[0].w - row_a[1].w).abs() > weight_tol * row_a[0].w
            || (row_b[0].w - row_b[1].w).abs() > weight_tol * row_b[0].w
        {
            return false;
        }
        profile_is_nondegenerate |= (aa - origin_a).cross(&axis).norm() > spatial_tol;
    }
    if !profile_is_nondegenerate {
        return false;
    }
    // The complete two-control-point sweeps each describe an exact finite
    // axial interval; unlike a broad-phase bbox this is a support certificate.
    let offset_on_axis = offset.dot(&axis);
    let sweep_b_on_axis = sweep_b.dot(&axis);
    let (b0, b1) = (
        offset_on_axis.min(offset_on_axis + sweep_b_on_axis),
        offset_on_axis.max(offset_on_axis + sweep_b_on_axis),
    );
    0.0f64.max(b0) + spatial_tol < length_a.min(b1)
}

/// Certify a positive-area shared carrier between two complete rational
/// ruled-extrusion patches (same profile, independently trimmed sweep spans).
/// This does not prove a boolean face relation, wire trim membership or
/// orientation; downstream topology partition and classification still apply.
pub(in crate::transversal::classic) fn common_ruled_extrusion_strip(
    a: &BsplineSurface<Vector4>,
    b: &BsplineSurface<Vector4>,
) -> bool {
    let policy = DEFAULT_NUMERICAL_TOLERANCES;
    let spatial_tol = policy.boolean_ruled_structural_length;
    let direction_tol = policy.boolean_ruled_direction;
    let weight_tol = policy.boolean_ruled_weight_relative;
    let knot_tol = policy.boolean_ruled_knot_relative;
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
    if !normalized_ruled_knots_match(u_a, u_b, knot_tol)
        || v_a.len() != 4
        || v_b.len() != 4
        || !normalized_ruled_knots_match(v_a, v_b, knot_tol)
        || !((v_a[0] - v_a[1]).abs() / v_a.range_length() <= knot_tol
            && (v_a[2] - v_a[3]).abs() / v_a.range_length() <= knot_tol
            && (v_b[0] - v_b[1]).abs() / v_b.range_length() <= knot_tol
            && (v_b[2] - v_b[3]).abs() / v_b.range_length() <= knot_tol)
    {
        return false;
    }

    congruent_ruled_sweeps(net_a, net_b, spatial_tol, direction_tol, weight_tol)
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
    if axis_a.cross(&axis_b).norm() > axis_eps
        || (cb - ca).cross(&axis_a).norm() > spatial_eps
        || (ra - rb).abs() > spatial_eps
    {
        return false;
    }
    let a = (span_a.0.min(span_a.1), span_a.0.max(span_a.1));
    // Axial coordinates must share A's frame, even when B's sweep runs
    // in the opposite direction.
    let sign = if axis_a.dot(&axis_b) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let bx0 = sign * span_b.0;
    let bx1 = sign * span_b.1;
    let b = (bx0.min(bx1), bx0.max(bx1));
    a.0.max(b.0) + spatial_eps < a.1.min(b.1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::{KnotVector, Point3, Vector3};

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
    fn closed_rational_full_circle_is_certified_without_three_distinct_endpoints() {
        use monstertruck_geometry::prelude::{
            Line, RevolutionSurface, TryIntoHomogeneousBsplineSurface,
        };
        let radial = 0.4826;
        let origin = Point3::new(59.384_373_707_781_98, -5.257_812_234_042_557, 0.0);
        let start = origin + Vector3::x() * radial;
        let cylinder = RevolutionSurface::by_revolution(
            Line(start, start + Vector3::z()),
            origin,
            Vector3::z(),
        );
        let mut net = cylinder.try_into_homogeneous_bspline_surface().unwrap();
        net.swap_axes();
        let (center, axis, radius, _) = certified_circular_cylinder(&net)
            .expect("a closed 360-degree rational circle is structurally a cylinder");
        assert!((radius - radial).abs() < 1e-10);
        assert!((center - origin).cross(&axis).norm() < 1e-10);
        // Break a control point well beyond the coefficient proof tolerance.
        net.control_point_mut(1, 0).x += 0.002;
        assert!(
            certified_circular_cylinder(&net).is_none(),
            "a structurally perturbed rational circle is not a certified circle",
        );
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
    fn ruled_carrier_preserves_identity_across_affine_knot_domains() {
        let original = extruded((0.0, 2.0), 0.0, false, false);
        let mut u = original.knot_vector_u().clone();
        let mut v = original.knot_vector_v().clone();
        u.transform(1.0e-7, 1.0e-4);
        v.transform(300.0, -10.0);
        let reparameterized =
            BsplineSurface::new_unchecked((u, v), original.control_points().to_vec());
        assert!(common_ruled_extrusion_strip(&original, &reparameterized));
        assert!(common_ruled_extrusion_strip(&reparameterized, &original));
        let shared = extruded((0.5, 1.0), 0.0, false, false);
        assert!(common_ruled_extrusion_strip(&shared, &reparameterized));
    }

    #[test]
    fn ruled_carrier_rejects_small_angle_if_long_sweep_separates_in_3d() {
        let source = extruded((0.0, 10_000.0), 0.0, false, false);
        let tilted = |offset: f64| {
            let mut net = source.control_points().to_vec();
            for row in &mut net {
                row[1].y += offset * row[1].w;
            }
            BsplineSurface::new_unchecked(
                (
                    source.knot_vector_u().clone(),
                    source.knot_vector_v().clone(),
                ),
                net,
            )
        };
        // Large sweep times small angle must still obey world-space agreement.
        assert!(common_ruled_extrusion_strip(&source, &tilted(5.0e-11)));
        assert!(!common_ruled_extrusion_strip(&source, &tilted(5.0e-8)));
    }

    #[test]
    fn ruled_profile_knot_comparison_rejects_non_affine_distortion() {
        let a = KnotVector::from(vec![0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0]);
        let b = KnotVector::from(vec![4.0, 4.0, 4.0, 4.48, 5.2, 5.2, 5.2]);
        assert!(normalized_ruled_knots_match(&a, &b, 1.0e-10));
        let changed = KnotVector::from(vec![4.0, 4.0, 4.0, 4.481, 5.2, 5.2, 5.2]);
        assert!(!normalized_ruled_knots_match(&a, &changed, 1.0e-10));
        assert!(!normalized_ruled_knots_match(
            &KnotVector::from(vec![0.0, 0.0]),
            &KnotVector::from(vec![1.0, 1.0]),
            1.0e-10,
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
