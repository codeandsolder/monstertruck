//! Source-guided selection of a separated analytic branch.
//!
//! Positive-weight clamped spline control rays select the angular interval,
//! independently of radial incidence. Ordered control rays in each proper
//! local wedge prove single directed angular traversal. The original spline
//! remains a separate geometric object, not an exact intersection certificate.

use super::{CylinderCarrier, CylinderIntersectionArc, IntersectionBranch};
use crate::{
    analytic_curve_incidence::valid,
    analytic_radial_deviation::{RadialDeviationBounds, bound_nurbs_radial_deviation},
    prelude::*,
};
use std::f64::consts::{PI, TAU};

/// A newly proposed analytic intersection, selected using the entire source
/// control polygon. The source may be measurably off either carrier.
#[derive(Clone, Copy, Debug)]
pub struct SourceIntersectionProposal {
    pub arc: CylinderIntersectionArc,
    /// Positive minimum signed transverse branch coordinate over the full
    /// source control hull, including a finite-f64 conditioning allowance.
    pub source_branch_margin: f64,
    pub source_start_residual: f64,
    pub source_end_residual: f64,
    pub angular_arithmetic_margin: f64,
}

/// Whole-locus deviation from the original spline to a new analytical arc.
///
/// Bounds include finite-f64 conditioning allowances, not directed-rounding
/// interval arithmetic. Neither bound certifies the adjacent trimmed faces.
#[derive(Clone, Copy, Debug)]
pub struct BoundedSourceIntersection {
    pub proposal: SourceIntersectionProposal,
    pub primary_deviation: RadialDeviationBounds,
    pub secondary_deviation: RadialDeviationBounds,
    /// Lower witness to the Hausdorff distance.
    pub separation_lower: f64,
    /// Upper bound in both directions, using continuous angular coverage.
    pub separation_upper: f64,
}

impl BoundedSourceIntersection {
    /// Select and bound the same source curve in one operation, so incidence
    /// bounds cannot accidentally be attached to an unrelated branch.
    #[must_use = "inspect the whole-locus separation bounds"]
    pub fn from_spline(
        primary: CylinderCarrier,
        secondary: CylinderCarrier,
        source: &NurbsCurve<Vector4>,
        endpoint_tolerance: f64,
        target_gap: f64,
    ) -> Option<Self> {
        let proposal = SourceIntersectionProposal::from_spline(
            primary,
            secondary,
            source,
            endpoint_tolerance,
        )?;
        let bound = |carrier: CylinderCarrier| {
            bound_nurbs_radial_deviation(
                source,
                carrier.origin,
                carrier.axis,
                carrier.radius,
                0.0,
                target_gap,
            )
        };
        let primary_deviation = bound(primary)?;
        let secondary_deviation = bound(secondary)?;
        let arc = proposal.arc;
        let a = primary_deviation.upper;
        let e = a + secondary_deviation.upper;
        // Radially project a source point onto the primary cylinder while
        // preserving angle and axial position. Its displacement is <= a,
        // and its distance to the secondary cylinder is <= e. Its signed
        // transverse root coordinate remains positive if the hull margin
        // exceeds a. Factoring s²-s_exact² then bounds the axial correction.
        let source_root = proposal.source_branch_margin - a - arc.arithmetic_margin;
        if source_root <= 0.0 || !source_root.is_finite() {
            return None;
        }
        let axial = e * (2.0_f64.mul_add(secondary.radius, e))
            / ((source_root + arc.min_root) * arc.axis_sine);
        // Angular projection covers the entire directed interval, so the
        // same pointwise correspondence supplies the reverse bound too.
        let speed_bound =
            primary.radius * (1.0 + (1.0 + secondary.radius / arc.min_root) / arc.axis_sine);
        let rounding =
            speed_bound.mul_add(proposal.angular_arithmetic_margin, arc.arithmetic_margin);
        let separation_upper = a + axial + rounding;
        let separation_lower = (primary_deviation.lower.max(secondary_deviation.lower)
            - arc.arithmetic_margin)
            .max(0.0);
        if !separation_upper.is_finite() || separation_lower > separation_upper {
            return None;
        }
        Some(Self {
            proposal,
            primary_deviation,
            secondary_deviation,
            separation_lower,
            separation_upper,
        })
    }
}

fn clamped_continuous(curve: &BsplineCurve<Vector4>) -> bool {
    if !valid(curve) {
        return false;
    }
    let degree = curve.degree();
    let knots = curve.knot_vector();
    let count = curve.control_points().len();
    let lo = knots[degree];
    let hi = knots[count];
    let same = |a: f64, b: f64| {
        a.total_cmp(&b).is_eq() || a.abs().to_bits() == 0 && b.abs().to_bits() == 0
    };
    knots.iter().take(degree + 1).all(|&k| same(k, lo))
        && knots.iter().skip(count).all(|&k| same(k, hi))
        && knots[degree + 1] > lo
        && knots[count - 1] < hi
        && !(degree + 1..count).any(|i| same(knots[i - degree], knots[i]))
}

fn euclidean_controls(curve: &BsplineCurve<Vector4>) -> Option<Vec<Point3>> {
    curve
        .control_points()
        .iter()
        .map(|p| {
            let point = Point3::new(p.x / p.w, p.y / p.w, p.z / p.w);
            super::finite_point(point).then_some(point)
        })
        .collect()
}

fn ray_angle(first: Vector3, second: Vector3, axis: Vector3) -> f64 {
    axis.dot(first.cross(second)).atan2(first.dot(second))
}

/// Every active control set is contained in a strict angular wedge of width
/// less than π. Within that chart, increasing polar angle is the same as
/// increasing slope relative to a positive half-plane direction. Positive
/// rational B-spline bases preserve ordered control slopes; consequently
/// the projected curve traverses one directed interval, with no axis crossing,
/// extra turn, or hidden backward span. The exact endpoint ray difference,
/// lifted onto this interval, is used rather than an accumulated angle sum.
fn ordered_sweep(
    carrier: CylinderCarrier,
    curve: &BsplineCurve<Vector4>,
    points: &[Point3],
) -> Option<(f64, f64)> {
    let origin = carrier.origin;
    let axis = carrier.axis;
    let mut scale = [origin.x.abs(), origin.y.abs(), origin.z.abs(), 1.0_f64]
        .into_iter()
        .fold(1.0_f64, f64::max);
    let mut min_radius = f64::INFINITY;
    let mut rays = Vec::with_capacity(points.len());
    for &point in points {
        scale = scale
            .max(point.x.abs())
            .max(point.y.abs())
            .max(point.z.abs());
        let delta = point - origin;
        let radial = delta - axis * delta.dot(axis);
        let length = radial.magnitude();
        if !length.is_finite() || length <= 0.0 {
            return None;
        }
        min_radius = min_radius.min(length);
        rays.push(radial / length);
    }
    let margin = 1024.0 * f64::EPSILON * scale / min_radius;
    if !margin.is_finite() || margin >= 0.01 {
        return None;
    }
    let mut angles = Vec::with_capacity(points.len());
    angles.push(0.0);
    let mut sign = 0.0_f64;
    for pair in rays.windows(2) {
        let delta = ray_angle(pair[0], pair[1], axis);
        if !delta.is_finite() || delta.abs() <= margin || delta.abs() >= PI - margin {
            return None;
        }
        if sign.abs() > 0.0 && sign * delta <= 0.0 {
            return None;
        }
        sign = delta.signum();
        angles.push(angles.last()? + delta);
    }
    let degree = curve.degree();
    let knots = curve.knot_vector();
    for i in degree..points.len() {
        if knots[i] < knots[i + 1] && (angles[i] - angles[i - degree]).abs() >= PI - margin {
            return None;
        }
    }
    let total = *angles.last()?;
    let principal = ray_angle(*rays.first()?, *rays.last()?, axis);
    let sweep = ((total - principal) / TAU).round().mul_add(TAU, principal);
    // The accumulated angle is only used to identify an integral turn. A
    // margin proportional to the count accounts for finite-f64 unwrapping.
    let accumulated_margin = margin * crate::usize_to_f64(points.len());
    if !sweep.is_finite()
        || (sweep - total).abs() > accumulated_margin
        || sweep.abs() <= margin
        || sweep.abs() > TAU
    {
        return None;
    }
    Some((sweep, accumulated_margin))
}

impl SourceIntersectionProposal {
    /// Select the regular intersection branch and its complete angular
    /// interval from a positive-weight clamped spline and two independent
    /// cylindrical carriers. Both source endpoints must match the newly
    /// constructed analytic arc. Intermediate source points need not lie
    /// on either cylinder; incidence must be established separately.
    #[must_use = "inspect the independently constructed source-guided arc"]
    pub fn from_spline(
        primary: CylinderCarrier,
        secondary: CylinderCarrier,
        source: &NurbsCurve<Vector4>,
        endpoint_tolerance: f64,
    ) -> Option<Self> {
        let curve = source.non_rationalized();
        if !clamped_continuous(curve) {
            return None;
        }
        let points = euclidean_controls(curve)?;
        let (sweep, angular_arithmetic_margin) = ordered_sweep(primary, curve, &points)?;
        let start = *points.first()?;
        let end = *points.last()?;
        let axis_cross = primary.axis.cross(secondary.axis);
        let length = axis_cross.magnitude();
        if !length.is_finite() || length <= f64::EPSILON.sqrt() {
            return None;
        }
        let axial_projection = secondary.axis.cross(axis_cross / length).normalize();
        let branch = if (start - secondary.origin).dot(axial_projection) > 0.0 {
            IntersectionBranch::Positive
        } else {
            IntersectionBranch::Negative
        };
        let arc = CylinderIntersectionArc::new(
            primary,
            secondary,
            start,
            sweep,
            branch,
            endpoint_tolerance,
        )?;
        let source_start_residual = (arc.evaluate(0.0)? - start).magnitude();
        let source_end_residual = arc.end_residual(end, endpoint_tolerance)?;
        let scale = points
            .iter()
            .map(|&p| (p - secondary.origin).magnitude())
            .fold(1.0_f64, f64::max)
            .max(secondary.origin.to_vec().magnitude());
        let margin = 512.0 * f64::EPSILON * scale;
        let source_branch_margin = points
            .iter()
            .map(|&point| branch.sign() * (point - secondary.origin).dot(axial_projection))
            .fold(f64::INFINITY, f64::min)
            - margin;
        if !source_branch_margin.is_finite() || source_branch_margin <= 0.0 {
            return None;
        }
        Some(Self {
            arc,
            source_branch_margin,
            source_start_residual,
            source_end_residual,
            angular_arithmetic_margin,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_control_rays_choose_the_half_turn_and_reject_branch_switch() {
        let primary =
            CylinderCarrier::new(Point3::origin(), Vector3::unit_z(), 0.2).expect("primary");
        let secondary =
            CylinderCarrier::new(Point3::origin(), Vector3::unit_x(), 1.0).expect("secondary");
        // An ordinary polynomial cubic approximates the semicircle. Its
        // endpoints are exact; its interior is deliberately off both carriers.
        let controls = vec![
            Vector4::new(0.2, 0.0, 1.0, 1.0),
            Vector4::new(0.2, 0.2, 1.1, 1.0),
            Vector4::new(-0.2, 0.2, 1.1, 1.0),
            Vector4::new(-0.2, 0.0, 1.0, 1.0),
        ];
        // The active hull spans a full π, so no strict half-plane proof
        // exists for this single cubic. Two regular quarter spans do.
        let one = NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            controls,
        ));
        assert!(SourceIntersectionProposal::from_spline(primary, secondary, &one, 1e-9).is_none());
        let weight = std::f64::consts::FRAC_1_SQRT_2;
        let mid = (1.0_f64 - 0.2 * 0.2).sqrt();
        let mut points = vec![
            Vector4::new(0.2, 0.0, 1.0, 1.0),
            Vector4::new(0.2 * weight, 0.2 * weight, 1.0 * weight, weight),
            Vector4::new(0.0, 0.2, mid, 1.0),
            Vector4::new(-0.2 * weight, 0.2 * weight, 1.0 * weight, weight),
            Vector4::new(-0.2, 0.0, 1.0, 1.0),
        ];
        let knots = KnotVector::from(vec![0., 0., 0., 0.5, 0.5, 1., 1., 1.]);
        let source = NurbsCurve::new(BsplineCurve::new_unchecked(knots.clone(), points.clone()));
        let proposal = SourceIntersectionProposal::from_spline(primary, secondary, &source, 1e-9)
            .expect("two regular source quarter spans");
        assert!((proposal.arc.sweep() - PI).abs() < 1e-12);
        assert!(proposal.source_branch_margin > 0.9);
        let bounded =
            BoundedSourceIntersection::from_spline(primary, secondary, &source, 1e-9, 1e-9)
                .expect("whole-locus bounds");
        assert!(bounded.separation_lower > 0.0);
        assert!(bounded.separation_upper < 0.1);
        for i in 0..=200 {
            let p = source.evaluate(f64::from(i) / 200.0);
            let angle = p.y.atan2(p.x).clamp(0.0, PI);
            let q = bounded.proposal.arc.evaluate(angle).expect("arc point");
            assert!((p - q).magnitude() <= bounded.separation_upper);
        }
        points[2].z = -1.0;
        let bad = NurbsCurve::new(BsplineCurve::new_unchecked(knots, points));
        assert!(SourceIntersectionProposal::from_spline(primary, secondary, &bad, 1e-9).is_none());
    }
}
