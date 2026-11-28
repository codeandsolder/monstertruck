//! Whole-domain angular coverage of positive-weight rational circular curves.
use super::{CurveIncidence, certify_nurbs_curve_on_support};
use crate::{
    analytic_bezier_partition::proof_bezier_spans, analytic_carrier::AnalyticSupport, prelude::*,
};

const TURN: f64 = std::f64::consts::TAU;

/// Validated local frame for measuring circular-arc coverage.
#[derive(Clone, Copy)]
pub struct CircularArcFrame {
    center: Point3,
    normal: Vector3,
    x: Vector3,
    y: Vector3,
    radius: f64,
    tolerance: f64,
}

impl CircularArcFrame {
    /// Validate a circle and its angular origin. A starting point may differ
    /// from the exact circle by at most the requested geometric tolerance.
    #[must_use = "The newly constructed value is returned rather than retained by this constructor."]
    pub fn new(
        center: Point3,
        normal: Vector3,
        radius: f64,
        start: Point3,
        tolerance: f64,
    ) -> Option<Self> {
        let length = normal.magnitude();
        if !length.is_finite()
            || length <= 0.0
            || !radius.is_finite()
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || radius <= tolerance
        {
            return None;
        }
        let normal = normal / length;
        let offset = start - center;
        let axial = offset.dot(normal);
        let radial = offset - normal * axial;
        let length = radial.magnitude();
        if !length.is_finite()
            || length <= 0.0
            || axial.abs() > tolerance
            || (length - radius).abs() > tolerance
        {
            return None;
        }
        let x = radial / length;
        Some(Self {
            center,
            normal,
            x,
            y: normal.cross(x),
            radius,
            tolerance,
        })
    }

    /// Geometric tolerance expressed in radians on this circle.
    #[must_use = "The computed angle tolerance quantity is reported through this return value."]
    pub fn angle_tolerance(self) -> f64 { self.tolerance / self.radius }

    /// Evaluate a point in the positive angular direction around the axis.
    #[must_use = "The evaluated point on the circular carrier is returned, not stored."]
    pub fn point(self, angle: f64) -> Point3 {
        self.center + (self.x * angle.cos() + self.y * angle.sin()) * self.radius
    }

    fn ray(self, point: [f64; 4]) -> Option<[f64; 2]> {
        let offset = Vector3::new(
            (-self.center.x).mul_add(point[3], point[0]),
            (-self.center.y).mul_add(point[3], point[1]),
            (-self.center.z).mul_add(point[3], point[2]),
        );
        let ray = [offset.dot(self.x), offset.dot(self.y)];
        let length = ray[0].hypot(ray[1]);
        (length.is_finite() && length > 0.0).then(|| ray.map(|v| v / length))
    }
}

fn angle_between(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
        .atan2(a[0].mul_add(b[0], a[1] * b[1]))
}

/// The ordered Bernstein control rays must lie in a single proper angular
/// wedge. Positive weights keep the whole rational span within this wedge.
/// Pairwise ray ordering checks angular traversal within tolerance (the numerator of
/// the angular derivative is the sum of ordered Bernstein cross products).
/// This also excludes a hidden extra revolution between endpoint samples.
fn span_sweep(rays: &[[f64; 2]], tolerance: f64) -> Option<f64> {
    let first = *rays.first()?;
    let last = *rays.last()?;
    let sweep = rays
        .windows(2)
        .map(|pair| angle_between(pair[0], pair[1]))
        .sum::<f64>();
    if !sweep.is_finite() || sweep == 0.0 || sweep.abs() > std::f64::consts::PI + tolerance {
        return None;
    }
    // A common half-plane is required; do not accept angular backtracking
    // or a Bernstein polygon with a turn hidden by angle normalization.
    let direction = sweep.signum();
    for (i, &a) in rays.iter().enumerate() {
        for &b in &rays[i + 1..] {
            if direction * angle_between(a, b) < -tolerance {
                return None;
            }
        }
    }
    let end_angle = angle_between(first, last);
    let delta = (sweep - end_angle).abs();
    if delta.min((TURN - delta).abs()) > tolerance {
        return None;
    }
    Some(sweep)
}

/// Unwrapped angular interval on a certified circular support.
#[derive(Clone, Copy)]
pub struct CircularArcCoverage {
    /// Increasing angular endpoint, measured from the supplied arc origin.
    pub start_rad: f64,
    /// Positive angular length, at most one revolution.
    pub sweep_rad: f64,
    /// Whether the authored parameter increases the frame angle.
    pub forward: bool,
}

/// Prove every rational span lies on the supplied circle and traverses one
/// unambiguous angular interval. Supports either curve orientation.
#[must_use = "The optional certify circular arc coverage outcome must be inspected to distinguish a value from absence."]
pub fn certify_circular_arc_coverage(
    rational: &NurbsCurve<Vector4>,
    frame: CircularArcFrame,
) -> Option<CircularArcCoverage> {
    let plane = AnalyticSupport::Plane {
        origin: frame.center,
        normal: frame.normal,
    };
    let sphere = AnalyticSupport::Sphere {
        center: frame.center,
        radius: frame.radius,
    };
    if [plane, sphere].into_iter().any(|support| {
        certify_nurbs_curve_on_support(rational, support, frame.tolerance) != CurveIncidence::Proven
    }) {
        return None;
    }
    angular_coverage(rational, frame)
}

/// Prove continuous, single-turn angular coverage on a cylinder.
///
/// The frame
/// defines the transverse circle and the cylinder axis through its center.
/// Unlike the circular certificate this makes no claim about axial height:
/// an oblique planar ellipse additionally requires whole-curve plane incidence.
#[must_use = "The optional certify cylindrical arc coverage outcome must be inspected to distinguish a value from absence."]
pub fn certify_cylindrical_arc_coverage(
    rational: &NurbsCurve<Vector4>,
    frame: CircularArcFrame,
) -> Option<CircularArcCoverage> {
    let cylinder = AnalyticSupport::Cylinder {
        origin: frame.center,
        axis: frame.normal,
        radius: frame.radius,
    };
    if certify_nurbs_curve_on_support(rational, cylinder, frame.tolerance) != CurveIncidence::Proven
    {
        return None;
    }
    angular_coverage(rational, frame)
}

fn angular_coverage(
    rational: &NurbsCurve<Vector4>,
    frame: CircularArcFrame,
) -> Option<CircularArcCoverage> {
    let spans = proof_bezier_spans(rational.non_rationalized())?;
    let tolerance = frame.angle_tolerance();
    let mut total: f64 = 0.0;
    let mut first = None;
    let mut last = None;
    for span in spans {
        let rays = span
            .controls
            .into_iter()
            .map(|p| frame.ray(p))
            .collect::<Option<Vec<_>>>()?;
        let a = *rays.first()?;
        let b = *rays.last()?;
        let sweep = span_sweep(&rays, tolerance)?;
        if last.is_some_and(|previous| angle_between(previous, a).abs() > tolerance)
            || (total != 0.0 && total.signum() != sweep.signum())
        {
            return None;
        }
        first.get_or_insert(a);
        last = Some(b);
        total += sweep;
    }
    if total.abs() > TURN + tolerance {
        return None;
    }
    let start = if total > 0.0 { first? } else { last? };
    let mut from = start[1].atan2(start[0]).rem_euclid(TURN);
    if TURN - from <= tolerance {
        from = 0.0;
    }
    Some(CircularArcCoverage {
        start_rad: from,
        sweep_rad: total.abs(),
        forward: total > 0.0,
    })
}

/// Candidate circle from well-separated points of an authored spline.
/// This is NOT a geometric certificate. Callers must certify full-domain
/// incidence and, when needed, angular coverage before accepting it.
#[derive(Clone, Copy, Debug)]
pub struct CircleProposal {
    /// Candidate center.
    pub center: Point3,
    /// Candidate unit normal.
    pub normal: Vector3,
    /// Candidate radius.
    pub radius: f64,
}

/// Propose a circle without assuming the endpoints are distinct or the
/// parameterization is uniform. Every positive span supplies witnesses.
#[must_use = "the proposal requires a subsequent full-domain certificate"]
pub fn propose_circle(curve: &BsplineCurve<Vector4>, tolerance: f64) -> Option<CircleProposal> {
    use crate::analytic_cylinder_intersection::project_positive_homogeneous_point as project;
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let spans = proof_bezier_spans(curve)?;
    let point = |p: [f64; 4]| project(Vector4::new(p[0], p[1], p[2], p[3]));
    let mut witnesses = Vec::new();
    for span in spans {
        witnesses.push(point(*span.controls.first()?)?);
        witnesses.push(project(
            curve.evaluate(f64::midpoint(span.start, span.end)),
        )?);
        witnesses.push(point(*span.controls.last()?)?);
    }
    let first = *witnesses.first()?;
    let u = witnesses
        .iter()
        .map(|p| *p - first)
        .max_by(|a, b| a.magnitude2().total_cmp(&b.magnitude2()))?;
    let v = witnesses.iter().map(|p| *p - first).max_by(|a, b| {
        u.cross(*a)
            .magnitude2()
            .total_cmp(&u.cross(*b).magnitude2())
    })?;
    let normal = u.cross(v);
    let denominator = 2.0 * normal.magnitude2();
    if !denominator.is_finite() || denominator <= 0.0 {
        return None;
    }
    let displacement =
        (v.cross(normal) * u.magnitude2() + normal.cross(u) * v.magnitude2()) / denominator;
    let center = first + displacement;
    let radius = displacement.magnitude();
    (radius.is_finite() && radius > tolerance).then_some(CircleProposal {
        center,
        normal: normal.normalize(),
        radius,
    })
}
