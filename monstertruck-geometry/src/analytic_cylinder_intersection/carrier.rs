//! Structural circle and straight-sweep certificates for rational cylinder nets.
use crate::prelude::*;
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;

/// Dehomogenize a finite positive-weight spatial control point.
#[must_use = "invalid homogeneous weights must be handled"]
pub fn project_positive_homogeneous_point(point: Vector4) -> Option<Point3> {
    if !point.w.is_finite() || point.w <= 0.0 {
        return None;
    }
    let projected = Point3::new(point.x / point.w, point.y / point.w, point.z / point.w);
    (projected.x.is_finite() && projected.y.is_finite() && projected.z.is_finite())
        .then_some(projected)
}

/// Certify a full quadratic rational Bezier span algebraically, not from
/// sampled points. Every coefficient of the degree-four circle polynomial
/// must vanish within its accumulated homogeneous roundoff budget.
/// The single homogeneous circle identity shared by all Bezier knot spans.
#[derive(Clone, Copy)]
struct CircleEvidence {
    origin: Point3,
    basis_u: Vector3,
    basis_v: Vector3,
    cx: f64,
    cy: f64,
    radius2: f64,
    coefficient_eps: f64,
}

impl CircleEvidence {
    fn certify_span(&self, segment: &BsplineCurve<Vector4>) -> Option<()> {
        let Self {
            origin,
            basis_u,
            basis_v,
            cx,
            cy,
            radius2,
            coefficient_eps,
        } = *self;
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
            x[i] = cx.mul_add(-h.w, relative.dot(basis_u));
            y[i] = cy.mul_add(-h.w, relative.dot(basis_v));
            w[i] = h.w;
        }
        let power = |c: [f64; 3]| {
            [
                c[0],
                2.0 * (c[1] - c[0]),
                (-2.0_f64).mul_add(c[1], c[0] + c[2]),
            ]
        };
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
        (scale.is_finite()
            && scale > 0.0
            && residual
                .iter()
                .all(|value| value.abs() <= coefficient_eps * scale))
        .then_some(())
    }
}

fn column_weight_scale(net: &[Vec<Vector4>], column: usize) -> Option<f64> {
    let scale = net.iter().try_fold(0.0_f64, |maximum, row| {
        let weight = row.get(column)?.w;
        (weight.is_finite() && weight > 0.0).then_some(maximum.max(weight))
    })?;
    (scale > 0.0).then_some(scale)
}

fn normalized_profile(surface: &BsplineSurface<Vector4>) -> Option<BsplineCurve<Vector4>> {
    let net = surface.control_points();
    let scale = column_weight_scale(net, 0)?;
    let controls = net
        .iter()
        .map(|row| {
            let normalized = row[0] / scale;
            project_positive_homogeneous_point(normalized)?;
            Some(normalized)
        })
        .collect::<Option<Vec<_>>>()?;
    BsplineCurve::new(surface.knot_vector_u().clone(), controls).ok()
}

/// The second axial column must be a translated copy with proportional
/// positive weights. Homogeneous scale is arbitrary; comparing absolute
/// weight differences would accept unrelated tiny-weight rational profiles.
/// Bound the effect of the residual relative weight perturbation on the
/// full tensor-product control hull, not just the two endpoint columns.
fn certify_axial_tensor(
    net: &[Vec<Vector4>],
    origin: Point3,
    sweep: Vector3,
    spatial_eps: f64,
    weight_eps: f64,
) -> Option<()> {
    let scales = [column_weight_scale(net, 0)?, column_weight_scale(net, 1)?];
    let mut weight_error = 0.0_f64;
    let mut displacement_error = 0.0_f64;
    let mut hull_radius = 0.0_f64;
    for row in net {
        let start = project_positive_homogeneous_point(row[0])?;
        let finish = project_positive_homogeneous_point(row[1])?;
        let first = row[0].w / scales[0];
        let second = row[1].w / scales[1];
        if first <= 0.0 || second <= 0.0 {
            return None;
        }
        let relative = (second - first).abs() / first;
        if !relative.is_finite() || relative > weight_eps {
            return None;
        }
        weight_error = weight_error.max(relative);
        displacement_error = displacement_error.max((finish - start - sweep).magnitude());
        hull_radius = hull_radius
            .max((start - origin).magnitude())
            .max((start + sweep - origin).magnitude());
    }
    // Relative weight errors change any normalized positive control average
    // by at most δ/(1-δ) times the ideal control-hull diameter.
    let bound =
        (weight_error / (1.0 - weight_error)).mul_add(2.0 * hull_radius, displacement_error);
    (bound.is_finite() && bound <= spatial_eps).then_some(())
}

/// Algebraically certify a circular-cylinder carrier represented by a
/// quadratic rational B-spline profile with a linear extruded second axis.
///
/// Each nondegenerate knot span is converted to a quadratic Bezier span.
/// The quadratic circle equation in homogeneous coordinates then has degree
/// four: checking all five *polynomial coefficients* proves the entire span
/// lies on the circle, not merely finitely many sampled profile points.
#[must_use = "inspect the structurally certified cylinder and its finite axial interval"]
pub fn certified_circular_cylinder(
    surface: &BsplineSurface<Vector4>,
) -> Option<(Point3, Vector3, f64, (f64, f64))> {
    let policy = DEFAULT_NUMERICAL_TOLERANCES;
    let spatial_eps = policy.boolean_cylinder_structural_length;
    let weight_eps = policy.boolean_cylinder_homogeneous_weight;
    let coefficient_eps = policy.boolean_cylinder_polynomial_relative;
    let (degree_u, degree_v) = surface.degrees();
    let net = surface.control_points();
    if degree_u != 2 || degree_v != 1 || net.len() < 3 || net.iter().any(|row| row.len() != 2) {
        return None;
    }
    let origin = project_positive_homogeneous_point(net[0][0])?;
    let end = project_positive_homogeneous_point(net[0][1])?;
    let sweep = end - origin;
    if sweep.magnitude() <= spatial_eps {
        return None;
    }
    let axis = sweep.normalize();
    certify_axial_tensor(net, origin, sweep, spatial_eps, weight_eps)?;

    let reference = if axis.x.abs() < policy.boolean_cylinder_reference_axis_component_limit {
        Vector3::unit_x()
    } else {
        Vector3::unit_y()
    };
    let basis_u = axis.cross(reference).normalize();
    let basis_v = axis.cross(basis_u).normalize();
    let profile = normalized_profile(surface)?;
    let knots = profile.knot_vector().clone();
    let lo = *knots.first()?;
    let hi = *knots.last()?;
    let xy = |h: Vector4| -> Option<(f64, f64)> {
        let p = project_positive_homogeneous_point(h)?;
        let delta = p - origin;
        Some((delta.dot(basis_u), delta.dot(basis_v)))
    };
    // Endpoints of a *closed* circular B-spline coincide. Fitting a center
    // from (front, middle, back) is then singular, even though every rational
    // quadratic knot span can be algebraically certified. Prefer the widest
    // triple for open arcs, then use separated interior samples for closed
    // profiles. This is only a center proposal; every homogeneous polynomial
    // coefficient on every span is still required to prove a circle below.
    let span = hi - lo;
    let (cx, cy, radius2) = [
        [0.0, 0.5, 1.0],
        [0.0, 0.25, 0.5],
        [0.25, 0.5, 0.75],
        [0.0, 0.375, 0.75],
    ]
    .into_iter()
    .find_map(|fractions| {
        let [f0, f1, f2] = fractions;
        let p0 = xy(profile.evaluate(f64::mul_add(span, f0, lo)))?;
        let p1 = xy(profile.evaluate(f64::mul_add(span, f1, lo)))?;
        let p2 = xy(profile.evaluate(f64::mul_add(span, f2, lo)))?;
        let twice_det = 2.0 * (p2.0 - p0.0).mul_add(-(p1.1 - p0.1), (p1.0 - p0.0) * (p2.1 - p0.1));
        if !twice_det.is_finite() || twice_det.abs() <= spatial_eps * spatial_eps {
            return None;
        }
        let d1 = p1.0.mul_add(p1.0, p1.1 * p1.1) - p0.0.mul_add(p0.0, p0.1 * p0.1);
        let d2 = p2.0.mul_add(p2.0, p2.1 * p2.1) - p0.0.mul_add(p0.0, p0.1 * p0.1);
        let cx = d2.mul_add(-(p1.1 - p0.1), d1 * (p2.1 - p0.1)) / twice_det;
        let cy = (p2.0 - p0.0).mul_add(-d1, (p1.0 - p0.0) * d2) / twice_det;
        let radius2 = (p0.0 - cx).mul_add(p0.0 - cx, (p0.1 - cy).powi(2));
        (radius2.is_finite() && radius2 > spatial_eps * spatial_eps).then_some((cx, cy, radius2))
    })?;

    // Every degree-two homogeneous Bezier span must satisfy
    // (X - cx W)^2 + (Y - cy W)^2 - r² W² == 0 identically.
    // Convert this product to power-basis coefficients, not to a collection of
    // point samples. Nonzero knot intervals are split by exact knot insertion.

    let evidence = CircleEvidence {
        origin,
        basis_u,
        basis_v,
        cx,
        cy,
        radius2,
        coefficient_eps,
    };
    let mut remaining = profile;
    let knot_slack = (policy.boolean_iso_parameter_roundoff_ulps * f64::EPSILON).mul_add(
        lo.abs().max(hi.abs()),
        policy.boolean_cylinder_knot_relative * (hi - lo).abs(),
    );
    let mut interior = knots
        .iter()
        .copied()
        .filter(|&k| k > lo + knot_slack && k < hi - knot_slack)
        .collect::<Vec<_>>();
    interior.dedup_by(|a, b| (*a - *b).abs() <= knot_slack);
    for knot in interior {
        let tail_segment = remaining.cut(knot);
        evidence.certify_span(&remaining)?;
        remaining = tail_segment;
    }
    evidence.certify_span(&remaining)?;
    let center = origin + basis_u * cx + basis_v * cy;
    let axial = (origin.to_vec().dot(axis), end.to_vec().dot(axis));
    Some((center, axis, radius2.sqrt(), axial))
}

/// Orientation and finite domain of an independently certified circular sweep.
#[derive(Clone, Copy, Debug)]
pub struct CircularCylinderChart {
    pub origin: Point3,
    pub axis: Vector3,
    pub radius: f64,
    pub axial_length: f64,
    pub radial_start: Vector3,
    pub sweep: f64,
    pub angular_first: bool,
}

/// Certify a single-cover directed circular profile and its linear axial
/// extrusion. A point-set cylinder certificate alone is insufficient for
/// recharting an oriented trimmed surface: it could cover a circle twice
/// or reverse the angular direction in an interior knot interval.
#[must_use = "inspect the complete directed cylinder chart certificate"]
pub fn certified_circular_cylinder_chart(
    mut net: BsplineSurface<Vector4>,
) -> Option<CircularCylinderChart> {
    use super::{
        CylinderCarrier,
        source::{clamped_continuous, euclidean_controls, ordered_sweep},
    };
    use crate::analytic_curve_incidence::{CircularArcFrame, certify_circular_arc_coverage};
    let angular_first = match net.degrees() {
        (2, 1) => true,
        (1, 2) => {
            net.swap_axes();
            false
        }
        _ => return None,
    };
    if !net.is_clamped() {
        return None;
    }
    let (origin, axis, radius, _) = certified_circular_cylinder(&net)?;
    let profile = normalized_profile(&net)?;
    if !clamped_continuous(&profile) {
        return None;
    }
    let controls = euclidean_controls(&profile)?;
    let carrier = CylinderCarrier::new(origin, axis, radius)?;
    let (sweep, _) = ordered_sweep(carrier, &profile, &controls)?;
    let start = *controls.first()?;
    let frame = CircularArcFrame::new(
        origin,
        axis,
        radius,
        start,
        DEFAULT_NUMERICAL_TOLERANCES.boolean_cylinder_structural_length,
    )?;
    // Also prove the profile plane; a circular axial projection of a
    // nonplanar spline is not the circular sweep being reparameterized.
    certify_circular_arc_coverage(&NurbsCurve::new(profile), frame)?;
    let offset = start - origin;
    let radial = offset - axis * offset.dot(axis);
    let front = project_positive_homogeneous_point(net.control_points()[0][0])?;
    let back = project_positive_homogeneous_point(net.control_points()[0][1])?;
    let axial_length = (back - front).magnitude();
    (axial_length.is_finite() && axial_length > 0.0).then(|| CircularCylinderChart {
        origin,
        axis,
        radius,
        axial_length,
        radial_start: radial.normalize(),
        sweep,
        angular_first,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quarter_net(scales: [f64; 2], wrong_middle_weight: bool) -> BsplineSurface<Vector4> {
        let w = 0.5_f64.sqrt();
        let profile = [
            Vector4::new(2.0, 0.0, 0.0, 1.0),
            Vector4::new(2.0 * w, 2.0 * w, 0.0, w),
            Vector4::new(0.0, 2.0, 0.0, 1.0),
        ];
        let controls = profile
            .into_iter()
            .enumerate()
            .map(|(i, p)| {
                let end_scale = if wrong_middle_weight && i == 1 {
                    scales[1] * 2.0
                } else {
                    scales[1]
                };
                vec![
                    p * scales[0],
                    (p + Vector4::new(0.0, 0.0, p.w, 0.0)) * end_scale,
                ]
            })
            .collect();
        BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
            controls,
        )
    }

    #[test]
    fn axial_tensor_certificate_is_projectively_scale_invariant() {
        for scales in [[1.0, 1.0], [1.0, 3.0], [1e-100, 1e-80], [1e100, 1e80]] {
            let net = quarter_net(scales, false);
            let chart = certified_circular_cylinder_chart(net).expect("proportional columns");
            assert!((chart.radius - 2.0).abs() < 1e-12);
            assert!((chart.axial_length - 1.0).abs() < 1e-12);
            assert!((chart.sweep - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
            let wrong = quarter_net(scales, true);
            assert!(certified_circular_cylinder(&wrong).is_none());
        }
    }

    #[test]
    fn full_revolution_has_one_directed_cylinder_chart() {
        let source = RevolutionSurface::by_revolution(
            Line(Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 0.0, 3.0)),
            Point3::origin(),
            Vector3::unit_z(),
        );
        let net = source
            .try_into_homogeneous_bspline_surface()
            .expect("cylinder net");
        let chart = certified_circular_cylinder_chart(net).expect("full turn");
        assert!((chart.sweep.abs() - std::f64::consts::TAU).abs() < 1e-12);
        assert!((chart.axial_length - 3.0).abs() < 1e-12);
    }
}
