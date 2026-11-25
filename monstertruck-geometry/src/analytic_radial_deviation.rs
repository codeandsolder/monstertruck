//! Global bracket on the radial deviation of a complete rational spline from
//! a cylindrical or single-nappe conical surface.
//!
//! Every genuinely positive knot span is decomposed into homogeneous Bezier
//! pieces. The degree-d radial-square residual is bounded by a degree-2d
//! Bernstein polynomial; positive weights and expected radii give a global
//! bound on the *radial* (not shortest Euclidean) deviation. Worst-first
//! subdivision tightens the bracket without assuming parameter correspondence.

use crate::analytic_bezier_partition::proof_bezier_spans;
use crate::analytic_bezier_refinement::refine_maximum;
use crate::analytic_curve_incidence::{reference_cone_polynomial, valid};
use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

/// Bracket on maximum shortest Euclidean distance to a cylinder or single-nappe cone.
///
/// At every supported curve point the expected radius is positive. The closest
/// point on the physical nappe is therefore at a positive generatrix height,
/// and distance = `|rho - (reference_radius + slope*axial)| / hypot(1,slope)`.
/// Radial bounds are exposed separately for physical source diagnostics.
///
/// The interval includes a conservative f64 conditioning allowance but does
/// not claim formally outward-rounded IEEE interval arithmetic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadialDeviationBounds {
    /// Maximum 3-D shortest-distance lower witness across the complete curve.
    pub lower: f64,
    /// Complete-curve 3-D shortest-distance upper bound.
    pub upper: f64,
    /// Corresponding radial-direction maximum lower witness.
    pub radial_lower: f64,
    /// Corresponding radial-direction complete-curve maximum upper bound.
    pub radial_upper: f64,
    /// Conservatively padded 3-D shortest-distance floating-point allowance.
    pub arithmetic_margin: f64,
    pub subdivisions: usize,
}

impl RadialDeviationBounds {
    #[must_use = "The boolean result of `within` must be observed to use the predicate."]
    pub fn within(self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.upper <= tolerance
    }

    #[must_use = "The boolean result of `exceeds` must be observed to use the predicate."]
    pub fn exceeds(self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance >= 0.0 && self.lower > tolerance
    }
}

fn expected_radius(
    point: &[f64; 4],
    origin: Point3,
    axis: Vector3,
    reference_radius: f64,
    slope: f64,
) -> Option<(Vector3, f64)> {
    let [x, y, z, w] = *point;
    if ![x, y, z, w].iter().all(|value| value.is_finite()) || w <= 0.0 {
        return None;
    }
    let delta = Vector3::new(x / w - origin.x, y / w - origin.y, z / w - origin.z);
    let axial = delta.dot(axis);
    let radius = slope.mul_add(axial, reference_radius);
    (radius.is_finite() && radius > 0.0 && axial.is_finite()).then_some((delta, radius))
}

fn residual(
    point: &[f64; 4],
    origin: Point3,
    axis: Vector3,
    reference_radius: f64,
    slope: f64,
) -> Option<f64> {
    let (delta, radius) = expected_radius(point, origin, axis, reference_radius, slope)?;
    let radial = (delta - axis * delta.dot(axis)).magnitude();
    let result = (radial - radius).abs();
    result.is_finite().then_some(result)
}

fn span_upper(
    controls: &[[f64; 4]],
    origin: Point3,
    axis: Vector3,
    reference_radius: f64,
    slope: f64,
) -> Option<f64> {
    let min_weight = controls
        .iter()
        .map(|point| point[3])
        .fold(f64::INFINITY, f64::min);
    let min_radius = controls.iter().try_fold(f64::INFINITY, |minimum, point| {
        Some(minimum.min(expected_radius(point, origin, axis, reference_radius, slope)?.1))
    })?;
    // Project the entire radial control polygon onto the first control's
    // transverse direction. The positive-weight rational convex hull proves
    // rho >= min(projected radius) on the entire subspan. This bound tightens
    // under subdivision, unlike rho >= 0, which leaves a factor-of-two gap.
    let first_delta = expected_radius(controls.first()?, origin, axis, reference_radius, slope)?.0;
    let radial = first_delta - axis * first_delta.dot(axis);
    let radial_length = radial.magnitude();
    let min_radial = if radial_length > 0.0 && radial_length.is_finite() {
        let direction = radial / radial_length;
        controls
            .iter()
            .try_fold(f64::INFINITY, |minimum, control| {
                let delta = expected_radius(control, origin, axis, reference_radius, slope)?.0;
                let projection = (delta - axis * delta.dot(axis)).dot(direction);
                projection.is_finite().then_some(minimum.min(projection))
            })?
            .max(0.0)
    } else {
        0.0
    };
    // The residual polynomial is (rho² - R²)*w². Positive-weight convexity
    // gives rho >= min_radial, R >= min_radius and w >= min_weight everywhere,
    // hence |rho - R| <= max(|F_i|)/(min_w²*(min_radial+min_R)).
    let polynomial = reference_cone_polynomial(controls, origin, axis, reference_radius, slope)?;
    let maximum = polynomial.iter().try_fold(0.0_f64, |max, coefficient| {
        coefficient
            .is_finite()
            .then_some(max.max(coefficient.abs()))
    })?;
    let denominator = min_weight.powi(2) * (min_radius + min_radial);
    let upper = maximum / denominator;
    upper.is_finite().then_some(upper)
}

fn arithmetic_margin(
    controls: &[[f64; 4]],
    origin: Point3,
    axis: Vector3,
    reference_radius: f64,
    slope: f64,
) -> Option<f64> {
    let mut min_weight = f64::INFINITY;
    let mut max_weight = 0.0_f64;
    let mut min_radius = f64::INFINITY;
    let mut scale = reference_radius.max(1.0);
    let origin_scale = [origin.x.abs(), origin.y.abs(), origin.z.abs()]
        .into_iter()
        .fold(1.0_f64, f64::max);
    for point in controls {
        min_weight = min_weight.min(point[3]);
        max_weight = max_weight.max(point[3]);
        let (delta, radius) = expected_radius(point, origin, axis, reference_radius, slope)?;
        min_radius = min_radius.min(radius);
        scale = scale.max(delta.magnitude()).max(radius);
    }
    let degree = crate::usize_to_f64(controls.len());
    let ratio = max_weight / min_weight;
    // Account separately for subtracting large global coordinates from the
    // local cone frame, and the homogeneous Bernstein quadratic convolution.
    // Like the existing plane bracket, this is a finite-f64 conditioning
    // allowance, not an outward-rounded IEEE interval certificate.
    let roundoff = 512.0
        * f64::EPSILON
        * degree
        * degree
        * ratio
        * ratio
        * (origin_scale + (1.0 + slope.abs()).powi(2) * scale * scale / min_radius);
    roundoff.is_finite().then_some(roundoff)
}

/// Bound maximum shortest Euclidean distance of a complete rational B-spline
/// to the radial carrier `rho = reference_radius + slope * axial`.
///
/// The axis must be unit length, slope nonnegative (zero is a cylinder), and the expected
/// radius positive at all dehomogenized Bezier controls, keeping the entire
/// curve on a single nappe. The reference plane need not pass through the
/// apex; this avoids cancellation from constructing a distant cone apex.
/// The upper bound is global across every positive original knot interval;
/// the lower bound is attained at a subspan endpoint. When the requested
/// gap cannot be reached, the valid wider bracket is returned.
///
/// # Returns
/// None for malformed geometry, invalid axis/radius/slope, an apex-crossing
/// control hull, or non-finite/unsupported conditioning.
///
/// The tolerance itself is NOT an input: callers apply independently sourced
/// geometric uncertainty to these bounds rather than altering tolerances.
#[must_use = "The optional bound nurbs radial deviation outcome must be inspected to distinguish a value from absence."]
pub fn bound_nurbs_radial_deviation(
    curve: &NurbsCurve<Vector4>,
    reference_origin: Point3,
    axis: Vector3,
    reference_radius: f64,
    slope: f64,
    target_gap: f64,
) -> Option<RadialDeviationBounds> {
    if ![reference_origin.x, reference_origin.y, reference_origin.z]
        .iter()
        .all(|v| v.is_finite())
        || !axis.magnitude2().is_finite()
        || (axis.magnitude2() - 1.0).abs() > DEFAULT_NUMERICAL_TOLERANCES.analytic_curve_axis_unit
        || !reference_radius.is_finite()
        || reference_radius <= 0.0
        || !slope.is_finite()
        || slope < 0.0
        || !valid(curve.non_rationalized())
    {
        return None;
    }
    let cone_normal = 1.0_f64.hypot(slope);
    let radial_target_gap = target_gap * cone_normal;
    if !radial_target_gap.is_finite() || radial_target_gap <= 0.0 {
        return None;
    }
    let spans = proof_bezier_spans(curve.non_rationalized())?;
    let max_margin = spans.iter().try_fold(0.0_f64, |margin, span| {
        Some(margin.max(arithmetic_margin(
            &span.controls,
            reference_origin,
            axis,
            reference_radius,
            slope,
        )?))
    })?;
    let bracket = refine_maximum(spans, radial_target_gap, |controls| {
        let first = residual(
            controls.first()?,
            reference_origin,
            axis,
            reference_radius,
            slope,
        )?;
        let last = residual(
            controls.last()?,
            reference_origin,
            axis,
            reference_radius,
            slope,
        )?;
        let upper = span_upper(controls, reference_origin, axis, reference_radius, slope)?;
        Some((first.max(last), upper))
    })?;
    let radial_margin = max_margin * crate::usize_to_f64(bracket.max_depth + 1);
    let radial_upper = bracket.upper + radial_margin;
    let radial_lower = (bracket.lower - radial_margin).max(0.0);
    let lower = radial_lower / cone_normal;
    let upper = radial_upper / cone_normal;
    if !upper.is_finite() || !lower.is_finite() || lower > upper {
        return None;
    }
    Some(RadialDeviationBounds {
        lower,
        upper,
        radial_lower,
        radial_upper,
        arithmetic_margin: radial_margin / cone_normal,
        subdivisions: bracket.subdivisions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn linear_cone_generator() -> NurbsCurve<Vector4> {
        NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![
                Vector4::new(1.0, 0.0, 1.0, 1.0),
                Vector4::new(4.0, 0.0, 4.0, 2.0),
            ],
        ))
    }

    #[test]
    fn cylindrical_distance_bounds_include_interior_departure() {
        let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(1.0, 0.0, 0.0, 1.0),
                Vector4::new(1.25, 0.0, 1.0, 1.0),
                Vector4::new(1.0, 0.0, 2.0, 1.0),
            ],
        ));
        let bound = bound_nurbs_radial_deviation(
            &curve,
            Point3::origin(),
            Vector3::unit_z(),
            1.0,
            0.0,
            1e-10,
        )
        .expect("valid cylinder");
        assert!(bound.lower <= 0.125 && bound.upper >= 0.125, "{bound:?}");
        assert!(bound.upper - bound.lower < 1e-8, "{bound:?}");
        assert_eq!(bound.radial_lower, bound.lower);
        assert_eq!(bound.radial_upper, bound.upper);
        assert!(
            bound_nurbs_radial_deviation(
                &curve,
                Point3::origin(),
                Vector3::unit_z(),
                1.0,
                -0.1,
                1e-10
            )
            .is_none()
        );
    }

    #[test]
    fn entire_rational_generator_is_on_cone_and_wrong_axis_is_rejected() {
        let cone = bound_nurbs_radial_deviation(
            &linear_cone_generator(),
            Point3::new(0.0, 0.0, 1.0),
            Vector3::unit_z(),
            1.0,
            1.0,
            1.0e-9,
        )
        .expect("rational linear generator is bounded");
        assert!(cone.within(1.0e-8), "{cone:?}");
        assert!(!cone.exceeds(1.0e-8));
        assert_eq!(cone.subdivisions, 0);
        assert!(
            bound_nurbs_radial_deviation(
                &linear_cone_generator(),
                Point3::new(0.0, 0.0, 1.0),
                Vector3::new(0.0, 0.0, 2.0),
                1.0,
                1.0,
                1.0e-9,
            )
            .is_none()
        );
    }

    #[test]
    fn mid_span_departure_is_not_hidden_by_on_cone_endpoints() {
        let bumped = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![
                Vector4::new(1.0, 0.0, 1.0, 1.0),
                Vector4::new(1.5, 0.0, 1.5, 1.0),
                Vector4::new(1.9, 0.0, 1.8, 1.0),
                Vector4::new(2.0, 0.0, 2.0, 1.0),
            ],
        ));
        let cone = bound_nurbs_radial_deviation(
            &bumped,
            Point3::new(0.0, 0.0, 1.0),
            Vector3::unit_z(),
            1.0,
            1.0,
            1.0e-10,
        )
        .expect("cubic departure is bounded");
        assert!(cone.exceeds(0.01), "{cone:?}");
        assert!(cone.upper - cone.lower < 1.0e-7, "{cone:?}");
        assert!((cone.radial_lower / cone.lower - 2.0_f64.sqrt()).abs() < 1.0e-8);
    }

    #[test]
    fn rational_multispan_global_bound_contains_independent_curve_evaluations() {
        let euclidean = [
            (1.0, 1.0, 1.0),
            (1.4, 1.4, 1.5),
            (2.02, 2.0, 0.7),
            (2.6, 2.65, 1.7),
            (3.0, 3.0, 1.0),
        ];
        let controls = euclidean
            .into_iter()
            .map(|(x, z, weight)| Vector4::new(x * weight, 0.0, z * weight, weight))
            .collect();
        let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0]),
            controls,
        ));
        let bounds = bound_nurbs_radial_deviation(
            &curve,
            Point3::new(0.0, 0.0, 1.0),
            Vector3::unit_z(),
            1.0,
            1.0,
            1.0e-9,
        )
        .expect("entire multi-span rational cone residual is bounded");
        assert!(bounds.lower > 0.0, "{bounds:?}");
        assert!(bounds.radial_lower >= bounds.lower);
        assert!(bounds.radial_upper >= bounds.upper);
        for index in 0..=1000 {
            let parameter = f64::from(index) / 1000.0;
            let point = curve.evaluate(parameter);
            let actual_radius = point.x.hypot(point.y);
            let expected = point.z;
            let actual_shortest_distance = (actual_radius - expected).abs() / 2.0_f64.sqrt();
            assert!(
                actual_shortest_distance <= bounds.upper + 1.0e-10,
                "t={parameter}: {actual_shortest_distance} > {bounds:?}"
            );
        }
    }

    #[test]
    fn do_not_certify_apex_crossing_curve_by_squared_cone_equation() {
        let crossing = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![
                Vector4::new(1.0, 0.0, 1.0, 1.0),
                Vector4::new(-1.0, 0.0, -1.0, 1.0),
            ],
        ));
        assert!(
            bound_nurbs_radial_deviation(
                &crossing,
                Point3::new(0.0, 0.0, 0.0),
                Vector3::unit_z(),
                0.1,
                1.0,
                1.0e-9,
            )
            .is_none()
        );
    }
}
