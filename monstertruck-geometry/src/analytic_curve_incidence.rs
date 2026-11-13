//! Global curve-on-carrier proof using homogeneous Bernstein coefficient bounds.
//!
//! Every rational spline knot interval is Bezierized; no sampled points can
//! promote an off-carrier curve into geometric equality.

use crate::analytic_bezier_partition::proof_bezier_spans;
use crate::analytic_carrier::AnalyticSupport;
use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use num_traits::ToPrimitive;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurveIncidence {
    /// All points on the entire supplied curve are within the stated tolerance.
    Proven,
    /// The coefficient bound does not establish incidence.
    NotProven,
    /// Input is malformed, underconstrained or cannot be handled safely.
    Unsupported,
}

// Explicit work budget: even a torus quartic at degree 64 has 257 coefficients.
const MAX_PROOF_DEGREE: usize = 64;

fn choose(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    (1..=k).fold(1.0, |acc, i| {
        // All degrees at this call site are bounded to 4 * MAX_PROOF_DEGREE.
        let numerator = (n - k + i).to_f64().unwrap_or(f64::INFINITY);
        let denominator = i.to_f64().unwrap_or(f64::INFINITY);
        acc * numerator / denominator
    })
}

fn multiply(first: &[f64], second: &[f64]) -> Option<Vec<f64>> {
    let n = first.len().checked_sub(1)?;
    let m = second.len().checked_sub(1)?;
    let degree = n.checked_add(m)?;
    let a = (0..=n).map(|i| choose(n, i)).collect::<Vec<_>>();
    let b = (0..=m).map(|i| choose(m, i)).collect::<Vec<_>>();
    let denominator = (0..=degree).map(|i| choose(degree, i)).collect::<Vec<_>>();
    let mut out = vec![0.0; degree + 1];
    for (i, &v) in first.iter().enumerate() {
        for (j, &w) in second.iter().enumerate() {
            out[i + j] = (v * w).mul_add(a[i] * b[j] / denominator[i + j], out[i + j]);
        }
    }
    out.iter().all(|v| v.is_finite()).then_some(out)
}

fn add_scaled(first: &[f64], second: &[f64], scale: f64) -> Option<Vec<f64>> {
    if first.len() != second.len() || !scale.is_finite() {
        return None;
    }
    first
        .iter()
        .zip(second)
        .map(|(&a, &b)| {
            let value = scale.mul_add(b, a);
            value.is_finite().then_some(value)
        })
        .collect()
}
fn squared(poly: &[f64]) -> Option<Vec<f64>> { multiply(poly, poly) }

fn centered(points: &[[f64; 4]], origin: Point3, direction: Vector3) -> Option<Vec<f64>> {
    points
        .iter()
        .map(|p| {
            let value = direction.x.mul_add(
                origin.x.mul_add(-p[3], p[0]),
                direction.y.mul_add(
                    origin.y.mul_add(-p[3], p[1]),
                    direction.z * (origin.z.mul_add(-p[3], p[2])),
                ),
            );
            value.is_finite().then_some(value)
        })
        .collect()
}

fn coordinates(points: &[[f64; 4]], origin: Point3) -> Option<[Vec<f64>; 4]> {
    Some([
        centered(points, origin, Vector3::unit_x())?,
        centered(points, origin, Vector3::unit_y())?,
        centered(points, origin, Vector3::unit_z())?,
        points.iter().map(|p| p[3]).collect(),
    ])
}

fn dot3(a: &[Vec<f64>; 3]) -> Option<Vec<f64>> {
    let first = squared(&a[0])?;
    let with_y = add_scaled(&first, &squared(&a[1])?, 1.0)?;
    add_scaled(&with_y, &squared(&a[2])?, 1.0)
}

fn quadric(
    points: &[[f64; 4]],
    origin: Point3,
    axis: Option<Vector3>,
    radius_squared: f64,
) -> Option<Vec<f64>> {
    let [x, y, z, w] = coordinates(points, origin)?;
    let mut norm = dot3(&[x, y, z])?;
    if let Some(axis) = axis {
        norm = add_scaled(&norm, &squared(&centered(points, origin, axis)?)?, -1.0)?;
    }
    add_scaled(&norm, &squared(&w)?, -radius_squared)
}

fn cone(points: &[[f64; 4]], apex: Point3, axis: Vector3, slope: f64) -> Option<Vec<f64>> {
    let [x, y, z, _] = coordinates(points, apex)?;
    let axial2 = squared(&centered(points, apex, axis)?)?;
    add_scaled(&dot3(&[x, y, z])?, &axial2, slope.mul_add(-slope, -1.0))
}

fn torus(
    points: &[[f64; 4]],
    center: Point3,
    axis: Vector3,
    major: f64,
    minor: f64,
) -> Option<Vec<f64>> {
    let [x, y, z, w] = coordinates(points, center)?;
    let norm = dot3(&[x, y, z])?;
    let w2 = squared(&w)?;
    let radial = add_scaled(&norm, &squared(&centered(points, center, axis)?)?, -1.0)?;
    let lhs = add_scaled(&norm, &w2, minor.mul_add(-minor, major * major))?;
    add_scaled(
        &squared(&lhs)?,
        &multiply(&radial, &w2)?,
        -4.0 * major * major,
    )
}

pub(crate) fn valid(curve: &BsplineCurve<Vector4>) -> bool {
    let cp = curve.control_points();
    let knots = curve.knot_vector();
    let n = cp.len();
    if n < 2 || knots.len() <= n + 1 {
        return false;
    }
    let degree = knots.len() - n - 1;
    degree > 0
        && degree < n
        && degree <= MAX_PROOF_DEGREE
        && knots.iter().all(|x| x.is_finite())
        && knots.iter().zip(knots.iter().skip(1)).all(|(a, b)| a <= b)
        && knots[degree] < knots[n]
        && cp
            .iter()
            .all(|p| [p.x, p.y, p.z, p.w].iter().all(|x| x.is_finite()) && p.w > 0.0)
}

fn valid_axis(axis: Vector3) -> bool {
    axis.magnitude2().is_finite()
        && (axis.magnitude2() - 1.0).abs() <= DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
}
fn finite_point(p: Point3) -> bool { [p.x, p.y, p.z].iter().all(|x| x.is_finite()) }

// Refuse cases in which floating-point conditioning could swallow a genuine
// geometric departure smaller than the arithmetic error of the implicit
// polynomial. The Bernstein hull bounds are algebraic; this check supplies
// a separate arithmetic safety margin before relying on those bounds.
fn well_conditioned(points: &[[f64; 4]], support: AnalyticSupport, tolerance: f64) -> bool {
    let min_weight = points.iter().map(|p| p[3]).fold(f64::INFINITY, f64::min);
    let max_weight = points.iter().map(|p| p[3]).fold(0.0_f64, f64::max);
    if min_weight <= 0.0 || !min_weight.is_finite() {
        return false;
    }
    let weight_ratio = max_weight / min_weight;
    let location_scale = points
        .iter()
        .flat_map(|p| p[..3].iter().map(move |v| v.abs() / p[3]))
        .fold(1.0_f64, f64::max);
    let reference_scale = match support {
        AnalyticSupport::Plane { origin, .. } => [origin.x, origin.y, origin.z, 0.0, 0.0],
        AnalyticSupport::Sphere { center, radius } => [center.x, center.y, center.z, radius, 0.0],
        AnalyticSupport::Cylinder { origin, radius, .. } => {
            [origin.x, origin.y, origin.z, radius, 0.0]
        }
        AnalyticSupport::Cone { apex, slope, .. } => [apex.x, apex.y, apex.z, slope, 0.0],
        AnalyticSupport::Torus {
            center,
            major_radius,
            minor_radius,
            ..
        } => [center.x, center.y, center.z, major_radius, minor_radius],
    }
    .into_iter()
    .map(f64::abs)
    .fold(1.0_f64, f64::max);
    let scale = location_scale.max(reference_scale);
    let degree = points
        .len()
        .saturating_sub(1)
        .to_f64()
        .unwrap_or(f64::INFINITY);
    let operation_budget = 256.0 * (degree + 1.0).powi(2);
    let resolution_floor = (operation_budget * f64::EPSILON) * scale * weight_ratio;
    resolution_floor.is_finite() && tolerance > resolution_floor
}

fn cone_bound(
    points: &[[f64; 4]],
    apex: Point3,
    axis: Vector3,
    slope: f64,
    positive_nappe: bool,
    tolerance: f64,
    min_w: f64,
) -> Option<(Vec<f64>, f64)> {
    if !finite_point(apex) || !valid_axis(axis) || !slope.is_finite() || slope <= 0.0 {
        return None;
    }
    let axial = centered(points, apex, axis)?;
    let min_axial = axial
        .iter()
        .zip(points)
        .map(|(a, p)| {
            if positive_nappe {
                *a / p[3]
            } else {
                -*a / p[3]
            }
        })
        .fold(f64::INFINITY, f64::min);
    let min_r = slope * min_axial;
    if !min_r.is_finite() || min_r <= tolerance {
        return None;
    }
    cone(points, apex, axis, slope).map(|p| {
        (
            p,
            tolerance * 2.0_f64.mul_add(min_r, -tolerance) * min_w * min_w,
        )
    })
}

fn torus_bound(
    points: &[[f64; 4]],
    center: Point3,
    axis: Vector3,
    major_radius: f64,
    minor_radius: f64,
    tolerance: f64,
    min_w: f64,
) -> Option<(Vec<f64>, f64)> {
    if !finite_point(center)
        || !valid_axis(axis)
        || !major_radius.is_finite()
        || !minor_radius.is_finite()
        || major_radius <= minor_radius
        || minor_radius <= tolerance
    {
        return None;
    }
    // F = (tube_distance^2-r^2)*((rho+R)^2+z^2-r^2).
    // The second factor is >= R^2-r^2 for a genuine ring torus.
    let limit = tolerance
        * 2.0_f64.mul_add(minor_radius, -tolerance)
        * minor_radius.mul_add(-minor_radius, major_radius * major_radius)
        * min_w.powi(4);
    torus(points, center, axis, major_radius, minor_radius).map(|p| (p, limit))
}

fn prove_span(points: &[[f64; 4]], support: AnalyticSupport, tolerance: f64) -> CurveIncidence {
    let min_w = points.iter().map(|p| p[3]).fold(f64::INFINITY, f64::min);
    if !min_w.is_finite() || min_w <= 0.0 || !well_conditioned(points, support, tolerance) {
        return CurveIncidence::Unsupported;
    }
    let result = match support {
        AnalyticSupport::Plane { origin, normal } => {
            if !finite_point(origin) || !valid_axis(normal) {
                return CurveIncidence::Unsupported;
            }
            centered(points, origin, normal).map(|p| (p, tolerance * min_w))
        }
        AnalyticSupport::Cylinder {
            origin,
            axis,
            radius,
        } => {
            if !finite_point(origin)
                || !valid_axis(axis)
                || !radius.is_finite()
                || radius <= tolerance
            {
                return CurveIncidence::Unsupported;
            }
            quadric(points, origin, Some(axis), radius * radius).map(|p| {
                (
                    p,
                    tolerance * (2.0_f64.mul_add(radius, -tolerance)) * min_w * min_w,
                )
            })
        }
        AnalyticSupport::Sphere { center, radius } => {
            if !finite_point(center) || !radius.is_finite() || radius <= tolerance {
                return CurveIncidence::Unsupported;
            }
            quadric(points, center, None, radius * radius).map(|p| {
                (
                    p,
                    tolerance * (2.0_f64.mul_add(radius, -tolerance)) * min_w * min_w,
                )
            })
        }
        AnalyticSupport::Cone {
            apex,
            axis,
            slope,
            positive_nappe,
        } => cone_bound(points, apex, axis, slope, positive_nappe, tolerance, min_w),
        AnalyticSupport::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => torus_bound(
            points,
            center,
            axis,
            major_radius,
            minor_radius,
            tolerance,
            min_w,
        ),
    };
    let Some((poly, bound)) = result else {
        return CurveIncidence::Unsupported;
    };
    if !bound.is_finite() || bound <= 0.0 {
        return CurveIncidence::Unsupported;
    }
    if poly.iter().all(|x| x.is_finite() && x.abs() <= bound) {
        CurveIncidence::Proven
    } else {
        CurveIncidence::NotProven
    }
}

/// Proves a global bounded-deviation claim for every knot span of a NURBS curve.
///
/// The curve and analytic support must be independently validated. A global
/// proof is sufficient for a trimmed edge, but not necessary: trimming a
/// partially off-carrier parent must instead use a separate range restriction.
#[must_use]
pub fn certify_nurbs_curve_on_support(
    curve: &NurbsCurve<Vector4>,
    support: AnalyticSupport,
    tolerance: f64,
) -> CurveIncidence {
    if !tolerance.is_finite() || tolerance <= 0.0 || !valid(curve.non_rationalized()) {
        return CurveIncidence::Unsupported;
    }
    let Some(pieces) = proof_bezier_spans(curve.non_rationalized()) else {
        return CurveIncidence::Unsupported;
    };
    for piece in &pieces {
        let result = prove_span(&piece.controls, support, tolerance);
        if result != CurveIncidence::Proven {
            return result;
        }
    }
    CurveIncidence::Proven
}

/// Weight-one embedding of a polynomial B-spline; same complete-domain proof.
#[must_use]
pub fn certify_bspline_curve_on_support(
    curve: &BsplineCurve<Point3>,
    support: AnalyticSupport,
    tolerance: f64,
) -> CurveIncidence {
    certify_nurbs_curve_on_support(&NurbsCurve::from(curve.clone()), support, tolerance)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(x: f64, y: f64, z: f64) -> Point3 { Point3::new(x, y, z) }
    fn point_support() -> AnalyticSupport {
        AnalyticSupport::Plane {
            origin: p(0.0, 0.0, 0.0),
            normal: Vector3::unit_z(),
        }
    }
    fn rational_circle(radius: f64) -> NurbsCurve<Vector4> {
        let q = 0.5_f64.sqrt();
        NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(radius, 0.0, 0.0, 1.0),
                Vector4::new(radius * q, radius * q, 0.0, q),
                Vector4::new(0.0, radius, 0.0, 1.0),
            ],
        ))
    }
    #[test]
    fn cubic_polynomial_on_plane_proves_whole_curve() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![
                p(0.0, 0.0, 0.0),
                p(1.0, 5.0, 0.0),
                p(2.0, -5.0, 0.0),
                p(3.0, 0.0, 0.0),
            ],
        );
        assert_eq!(
            certify_bspline_curve_on_support(&curve, point_support(), 1.0e-7),
            CurveIncidence::Proven
        );
        let mut bumped = curve;
        if let Some(second) = bumped.control_points_mut().nth(1) {
            second.z = 0.01;
        }
        assert_eq!(
            certify_bspline_curve_on_support(&bumped, point_support(), 1.0e-7),
            CurveIncidence::NotProven
        );
    }

    #[test]
    fn circular_rational_quadratic_proves_plane_cylinder_and_sphere() {
        let curve = rational_circle(2.0);
        let support = [
            point_support(),
            AnalyticSupport::Cylinder {
                origin: p(0.0, 0.0, 0.0),
                axis: Vector3::unit_z(),
                radius: 2.0,
            },
            AnalyticSupport::Sphere {
                center: p(0.0, 0.0, 0.0),
                radius: 2.0,
            },
        ];
        for s in support {
            assert_eq!(
                certify_nurbs_curve_on_support(&curve, s, 1.0e-7),
                CurveIncidence::Proven
            );
        }
        let wrong = AnalyticSupport::Cylinder {
            origin: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radius: 2.01,
        };
        assert_eq!(
            certify_nurbs_curve_on_support(&curve, wrong, 1.0e-7),
            CurveIncidence::NotProven
        );
    }

    #[test]
    fn circle_between_sample_deviation_cannot_pass_bernstein_bound() {
        let mut curve = rational_circle(2.0);
        curve.control_points_mut().nth(1).expect("middle control").x += 0.01;
        let cylinder = AnalyticSupport::Cylinder {
            origin: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radius: 2.0,
        };
        assert_eq!(
            certify_nurbs_curve_on_support(&curve, cylinder, 1.0e-7),
            CurveIncidence::NotProven
        );
    }

    #[test]
    fn cone_generatrix_is_proven_only_on_correct_nappe() {
        let curve = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![p(1.0, 0.0, 1.0), p(2.0, 0.0, 2.0)],
        );
        let cone = AnalyticSupport::Cone {
            apex: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            slope: 1.0,
            positive_nappe: true,
        };
        assert_eq!(
            certify_bspline_curve_on_support(&curve, cone, 1.0e-7),
            CurveIncidence::Proven
        );
        let opposite = AnalyticSupport::Cone {
            apex: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            slope: 1.0,
            positive_nappe: false,
        };
        assert_eq!(
            certify_bspline_curve_on_support(&curve, opposite, 1.0e-7),
            CurveIncidence::Unsupported
        );
    }

    #[test]
    fn ring_torus_outer_equatorial_circle_is_proven() {
        let circle = rational_circle(4.0);
        let ring = AnalyticSupport::Torus {
            center: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            major_radius: 3.0,
            minor_radius: 1.0,
        };
        assert_eq!(
            certify_nurbs_curve_on_support(&circle, ring, 1.0e-7),
            CurveIncidence::Proven
        );
        let wrong = AnalyticSupport::Torus {
            center: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            major_radius: 3.2,
            minor_radius: 1.0,
        };
        assert_eq!(
            certify_nurbs_curve_on_support(&circle, wrong, 1.0e-7),
            CurveIncidence::NotProven
        );
    }

    #[test]
    fn multiple_rational_quadratic_spans_are_all_checked() {
        let q = 0.5_f64.sqrt();
        let spline = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]),
            vec![
                Vector4::new(1.0, 0.0, 0.0, 1.0),
                Vector4::new(q, q, 0.0, q),
                Vector4::new(0.0, 1.0, 0.0, 1.0),
                Vector4::new(-q, q, 0.0, q),
                Vector4::new(-1.0, 0.0, 0.0, 1.0),
            ],
        ));
        let circle = AnalyticSupport::Cylinder {
            origin: p(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radius: 1.0,
        };
        assert_eq!(
            certify_nurbs_curve_on_support(&spline, circle, 1.0e-7),
            CurveIncidence::Proven
        );
        let mut changed = spline;
        changed
            .control_points_mut()
            .nth(3)
            .expect("fourth control")
            .x -= 0.1;
        assert_eq!(
            certify_nurbs_curve_on_support(&changed, circle, 1.0e-7),
            CurveIncidence::NotProven
        );
    }

    #[test]
    fn malformed_and_nonfinite_splines_refuse() {
        let wrong = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 1.0, 1.0]),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 0.0),
                Vector4::new(1.0, 0.0, 0.0, 1.0),
            ],
        ));
        assert_eq!(
            certify_nurbs_curve_on_support(&wrong, point_support(), 1.0e-7),
            CurveIncidence::Unsupported
        );
        assert_eq!(
            certify_nurbs_curve_on_support(&rational_circle(1.0), point_support(), f64::NAN),
            CurveIncidence::Unsupported
        );
        assert_eq!(
            certify_nurbs_curve_on_support(&rational_circle(1.0), point_support(), -1.0),
            CurveIncidence::Unsupported
        );
    }

    #[test]
    fn near_distinct_knots_prove_and_overmultiplicity_refuses() {
        let controls = vec![
            p(0.0, 0.0, 0.0),
            p(0.25, 0.0, 0.0),
            p(0.75, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
        ];
        let near_distinct = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5 + 1.0e-9, 1.0, 1.0]),
            controls.clone(),
        );
        assert_eq!(
            certify_bspline_curve_on_support(&near_distinct, point_support(), 1.0e-7),
            CurveIncidence::Proven,
        );
        let too_many = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0),
                p(0.2, 0.0, 0.0),
                p(0.3, 0.0, 0.0),
                p(0.7, 0.0, 0.0),
                p(1.0, 0.0, 0.0),
            ],
        );
        assert_eq!(
            certify_bspline_curve_on_support(&too_many, point_support(), 1.0e-7),
            CurveIncidence::Unsupported,
        );
    }

    #[test]
    fn unresolvable_large_coordinates_do_not_masquerade_as_proofs() {
        let huge = 1.0e12;
        let curve = rational_circle(huge);
        let sphere = AnalyticSupport::Sphere {
            center: p(0.0, 0.0, 0.0),
            radius: huge,
        };
        assert_eq!(
            certify_nurbs_curve_on_support(&curve, sphere, 1.0e-9),
            CurveIncidence::Unsupported,
        );
        // A range wide enough to resolve at this scale is still allowed.
        assert_eq!(
            certify_nurbs_curve_on_support(&curve, sphere, 1.0),
            CurveIncidence::Proven,
        );
    }
}
