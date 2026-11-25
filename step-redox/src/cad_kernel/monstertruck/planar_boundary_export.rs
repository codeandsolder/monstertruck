//! Preserve authored planar boundary isocurves during STEP export.
//! The carrier certificate and directed trim validation are independent of
//! the marching leader used by the CAD intersection solver.

use super::{step_polylines_match, validated_step_parameter_curve};
use monstertruck_core::tolerance::{
    DEFAULT_NUMERICAL_TOLERANCES, DEFAULT_STEP_RECOVERY_TOLERANCES,
};
use monstertruck_geometry::{
    analytic_curve_incidence::CurveIncidence, prelude::TryIntoHomogeneousBsplineSurface,
};
use monstertruck_modeling::{
    BoundedCurve, Curve, Cut, InnerSpace, Invertible, NurbsCurve, NurbsSurface, ParametricCurve,
    Plane, Point3, SearchNearestParameter, SearchParameter, Surface, SurfaceCurve, Vector4, algo,
    certify_analytic_curve_on_surface,
};

/// Preserve an exact authored boundary isocurve of a positive-weight
/// tensor-product spline. Clamped transverse knots make the outer control
/// row a true surface boundary. A planar certificate covers the entire
/// row; all other controls must lie strictly on one side of the plane.
/// For exact coplanarity this excludes any interior intersection branch.
/// For bounded coplanarity the caller must also match the evaluated trim.
fn certified_planar_boundary_isocurve(
    planar: &Plane,
    surface: &Surface,
    tolerance: f64,
) -> Option<NurbsCurve<Vector4>> {
    if !surface.homogeneous_net_covers_domain() || !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let spline = NurbsSurface::new(surface.try_into_homogeneous_bspline_surface()?);
    let net = spline.control_points();
    let first = net.first()?;
    if first.len() < 2 || net.len() < 2 || net.iter().any(|row| row.len() != first.len()) {
        return None;
    }
    let normal = planar.normal();
    let distances = net
        .iter()
        .map(|row| {
            row.iter()
                .map(|point| {
                    if point.w <= 0.0 || !point.w.is_finite() {
                        return None;
                    }
                    let p = Point3::new(point.x / point.w, point.y / point.w, point.z / point.w);
                    let distance = normal.dot(p - planar.origin());
                    distance.is_finite().then_some(distance)
                })
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()?;
    let mut selected = None;
    for transverse_u in [true, false] {
        let (count, degree, knots) = if transverse_u {
            (net.len(), spline.udegree(), spline.knot_vector_u())
        } else {
            (first.len(), spline.vdegree(), spline.knot_vector_v())
        };
        if degree == 0 || !knots.is_clamped(degree) {
            continue;
        }
        for boundary in [0, count - 1] {
            let mut on_plane = true;
            let mut positive = true;
            let mut negative = true;
            for (u, row) in distances.iter().enumerate() {
                for (v, &distance) in row.iter().enumerate() {
                    if (if transverse_u { u } else { v }) == boundary {
                        on_plane &= distance.abs() <= tolerance;
                    } else {
                        positive &= distance > tolerance;
                        negative &= distance < -tolerance;
                    }
                }
            }
            if !on_plane || !(positive || negative) {
                continue;
            }
            let curve = if transverse_u {
                spline.curve_v(boundary)
            } else {
                spline.curve_u(boundary)
            };
            if certify_analytic_curve_on_surface(
                &Curve::NurbsCurve(curve.clone()),
                &Surface::Plane(*planar),
                tolerance,
            ) != CurveIncidence::Proven
            {
                continue;
            }
            if selected.replace(curve).is_some() {
                return None;
            }
        }
    }
    selected
}

/// Replace only the export leader of a plane/spline intersection. The
/// source trace chooses a connected, directed subinterval of the
/// independently certified boundary row; both traces and any existing
/// face-local pcurves must still agree at the export accuracy.
pub(super) fn certified_planar_boundary_step_guide(
    original: &Curve,
    points: &[Point3],
    tolerance: f64,
) -> Option<Curve> {
    let Curve::Intersection(intersection) = original else {
        return None;
    };
    let ((Surface::Plane(plane), surface) | (surface, Surface::Plane(plane))) = (
        intersection.surface0().as_ref(),
        intersection.surface1().as_ref(),
    ) else {
        return None;
    };
    let proof_tolerance = DEFAULT_STEP_RECOVERY_TOLERANCES
        .solid_geometry_mm
        .min(tolerance);
    let mut exact = certified_planar_boundary_isocurve(plane, surface, proof_tolerance)?;
    if points.len() < 2 {
        return None;
    }
    let mut positions = Vec::with_capacity(points.len());
    for point in points {
        let parameter = exact
            .search_parameter(*point, None, 160)
            .or_else(|| exact.search_nearest_parameter(*point, None, 160))?;
        if (exact.evaluate(parameter) - *point).magnitude() > proof_tolerance {
            return None;
        }
        positions.push(parameter);
    }
    let start = *positions.first()?;
    let end = *positions.last()?;
    let (lo, hi) = exact.range_tuple();
    if ![start, end, lo, hi].iter().all(|v| v.is_finite())
        || (start - end).abs() <= DEFAULT_NUMERICAL_TOLERANCES.structural * (hi - lo).abs()
    {
        return None;
    }
    let increasing = end > start;
    let policy = DEFAULT_NUMERICAL_TOLERANCES;
    let param_tol = (policy.trim_parameter_roundoff_ulps * f64::EPSILON).mul_add(
        lo.abs().max(hi.abs()),
        policy.parameter_projection * (hi - lo).abs(),
    );
    if !positions.windows(2).all(|pair| {
        let delta = if increasing {
            pair[1] - pair[0]
        } else {
            pair[0] - pair[1]
        };
        delta >= -param_tol
    }) {
        return None;
    }
    let lower = start.min(end);
    let upper = start.max(end);
    if lower < lo - param_tol || upper > hi + param_tol {
        return None;
    }
    if lower > lo + param_tol {
        exact = exact.cut(lower);
    }
    if upper < hi - param_tol {
        let _discarded_tail = exact.cut(upper);
    }
    if !increasing {
        exact.invert();
    }
    let [Some(&front), Some(&back)] = [points.first(), points.last()] else {
        return None;
    };
    let (s, e) = exact.range_tuple();
    if (exact.evaluate(s) - front).magnitude() > proof_tolerance
        || (exact.evaluate(e) - back).magnitude() > proof_tolerance
    {
        return None;
    }
    let guide = Curve::NurbsCurve(exact);
    let (_, re_sampled) = algo::curve::parameter_division(&guide, guide.range_tuple(), tolerance);
    if !step_polylines_match(&re_sampled, points, tolerance)
        || ![intersection.boundary0(), intersection.boundary1()]
            .into_iter()
            .flatten()
            .all(|boundary| validated_step_parameter_curve(boundary, &re_sampled, tolerance))
    {
        return None;
    }
    let mut exported = SurfaceCurve::with_boundaries(
        intersection.surface0().clone(),
        intersection.surface1().clone(),
        Box::new(guide),
        intersection.boundary0().cloned(),
        intersection.boundary1().cloned(),
    );
    let (front, back) = intersection.endpoint_anchors();
    exported.set_endpoint_anchors(front, back);
    Some(Curve::Intersection(exported))
}

#[cfg(test)]
mod tests {
    use super::super::corrected_step_curve;
    use super::*;
    use anyhow::{Context, Result, bail};
    use monstertruck_modeling::{BsplineSurface, KnotVector};

    fn ruled_quarter_circle_surface(x_offsets: [f64; 2]) -> NurbsSurface<Vector4> {
        let w = std::f64::consts::FRAC_1_SQRT_2;
        let curve_points = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
        let weights = [1.0, w, 1.0];
        let control_points = (0..3)
            .map(|i| {
                (0..2)
                    .map(|j| {
                        let x = x_offsets[j];
                        let y = curve_points[i][0];
                        let other = curve_points[i][1];
                        Vector4::new(
                            x * weights[i],
                            y * weights[i],
                            other * weights[i],
                            weights[i],
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        NurbsSurface::new(BsplineSurface::new_unchecked(
            (
                KnotVector::from(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]),
                KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]),
            ),
            control_points,
        ))
    }

    #[test]
    fn step_export_prefers_exact_ruled_plane_isocurve_over_sampled_ssi() -> Result<()> {
        let ruled = ruled_quarter_circle_surface([0.0, 0.2]);
        let certified = certified_planar_boundary_isocurve(
            &Plane::yz(),
            &Surface::NurbsSurface(ruled.clone()),
            1.0e-9,
        )
        .context("failed to find exact authored rational row")?;
        let source = Curve::Intersection(SurfaceCurve::with_boundaries(
            Box::new(Surface::Plane(Plane::yz())),
            Box::new(Surface::NurbsSurface(ruled)),
            Box::new(Curve::NurbsCurve(certified)),
            None,
            None,
        ));
        let corrected = corrected_step_curve(&source, 1.0e-5)?;
        let Curve::Intersection(rewritten) = corrected else {
            bail!("lost original two-surface intersection");
        };
        let Curve::NurbsCurve(exact) = rewritten.leader().as_ref() else {
            bail!("exact source rational row was inappropriately sampled");
        };
        assert_eq!(exact.degree(), 2);
        assert!(exact.control_points().iter().any(|p| p.w < 1.0));
        for i in 0..=32 {
            let t = f64::from(i) / 32.0;
            let point = exact.evaluate(t);
            assert!(point.x.abs() <= 1.0e-12);
            assert!((point.y.hypot(point.z) - 1.0).abs() <= 1.0e-12);
            assert!((point - source.evaluate(t)).magnitude() <= 1.0e-10);
        }
        Ok(())
    }

    #[test]
    fn complete_ruled_plane_intersection_yields_exact_rational_iso_arc() -> Result<()> {
        let surface = ruled_quarter_circle_surface([0.0, 0.2]);
        let curve = certified_planar_boundary_isocurve(
            &Plane::yz(),
            &Surface::NurbsSurface(surface),
            1.0e-9,
        )
        .context("unique exact boundary row must be proved")?;
        assert_eq!(curve.degree(), 2);
        let p0 = curve.evaluate(0.0);
        let p1 = curve.evaluate(1.0);
        let middle = curve.evaluate(0.5);
        assert!((p0 - Point3::new(0.0, 1.0, 0.0)).magnitude() < 1.0e-12);
        assert!((p1 - Point3::new(0.0, 0.0, 1.0)).magnitude() < 1.0e-12);
        assert!(
            (middle
                - Point3::new(
                    0.0,
                    std::f64::consts::FRAC_1_SQRT_2,
                    std::f64::consts::FRAC_1_SQRT_2
                ))
            .magnitude()
                < 1.0e-12
        );
        // No amount of endpoint fitting can justify a second row crossing
        // or touching the plane: otherwise the planar intersection is not
        // a uniquely certified boundary row.
        assert!(
            certified_planar_boundary_isocurve(
                &Plane::yz(),
                &Surface::NurbsSurface(ruled_quarter_circle_surface([0.0, 0.0])),
                1e-9
            )
            .is_none()
        );
        Ok(())
    }

    #[test]
    fn planar_boundary_certifies_multispan_higher_degree_in_either_axis() -> Result<()> {
        let net = [0.3, 0.2, 0.1, 0.0]
            .into_iter()
            .map(|x| {
                vec![
                    Point3::new(x, 1.0, 0.0),
                    Point3::new(x, 1.0, 1.0),
                    Point3::new(x, 0.0, 1.0),
                ]
            })
            .collect();
        let mut surface = BsplineSurface::new_unchecked(
            (
                KnotVector::from(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]),
                KnotVector::bezier_knot(2),
            ),
            net,
        );
        for _ in 0..2 {
            let candidate = certified_planar_boundary_isocurve(
                &Plane::yz(),
                &Surface::BsplineSurface(surface.clone()),
                1.0e-9,
            )
            .context("missing exact multi-span polynomial surface boundary")?;
            assert_eq!(candidate.degree(), 2);
            for i in 0..=16 {
                let t = f64::from(i) / 16.0;
                let actual = candidate.evaluate(t);
                let expected = Point3::new(0.0, t.mul_add(-t, 1.0), t * (2.0 - t));
                assert!((actual - expected).magnitude() <= 1.0e-12);
            }
            surface.swap_axes();
        }
        Ok(())
    }

    #[test]
    fn step_export_exact_boundary_keeps_reversed_partial_trim() -> Result<()> {
        let ruled = Surface::NurbsSurface(ruled_quarter_circle_surface([0.0, 0.2]));
        let mut leader = certified_planar_boundary_isocurve(&Plane::yz(), &ruled, 1.0e-9)
            .context("missing original rational row")?;
        leader = leader.cut(0.2);
        let _tail = leader.cut(0.8);
        leader.invert();
        let source = Curve::Intersection(SurfaceCurve::with_boundaries(
            Box::new(ruled),
            Box::new(Surface::Plane(Plane::yz())),
            Box::new(Curve::NurbsCurve(leader.clone())),
            None,
            None,
        ));
        let corrected = corrected_step_curve(&source, 1.0e-5)?;
        let Curve::Intersection(rewritten) = corrected else {
            bail!("lost original intersection");
        };
        let Curve::NurbsCurve(exact) = rewritten.leader().as_ref() else {
            bail!("partial exact row was sampled");
        };
        assert_eq!(exact.degree(), 2);
        let (start, end) = exact.range_tuple();
        let (leader_start, leader_end) = leader.range_tuple();
        assert!((exact.evaluate(start) - leader.evaluate(leader_start)).magnitude() < 1.0e-10);
        assert!((exact.evaluate(end) - leader.evaluate(leader_end)).magnitude() < 1.0e-10);
        for i in 0..=16 {
            let t = f64::from(i) / 16.0;
            let point = exact.evaluate((end - start).mul_add(t, start));
            assert!(point.x.abs() < 1.0e-12);
            assert!((point.y.hypot(point.z) - 1.0).abs() < 1.0e-12);
        }
        Ok(())
    }
}
