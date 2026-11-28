//! Representation-independent bounded circle trims.
//! A rational edge becomes a circular edge only after complete incidence,
//! endpoint and ordered angular coverage have all been certified.

use crate::brep::{CircleSupport, CurveSupport, OrientedEdgeUse};
use crate::math3::{distance, normalize, sub};
use crate::solid_extrusions::circle_edge_sweep;
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES as TOL;
use monstertruck_geometry::{
    analytic_curve_incidence::{CircularArcFrame, certify_circular_arc_coverage, propose_circle},
    prelude::{BoundedCurve, NurbsCurve, Point3, Vector3, Vector4},
};

const GEOM_TOL_MM: f64 = TOL.solid_geometry_mm;

/// Propose a circle; the returned support is not yet a curve certificate.
pub(crate) fn propose_spline_circle(
    curve: &NurbsCurve<Vector4>,
    tolerance: f64,
) -> Option<CircleSupport> {
    let candidate = propose_circle(curve.non_rationalized(), tolerance)?;
    let center = candidate.center;
    let normal = candidate.normal;
    let first = curve.front();
    let x = normalize(
        [first.x - center.x, first.y - center.y, first.z - center.z],
        0.0,
    )?;
    Some(CircleSupport {
        center_mm: [center.x, center.y, center.z],
        radius_mm: candidate.radius,
        normal: [normal.x, normal.y, normal.z],
        x_direction: x,
    })
}

/// Normalize a bounded edge on a circular carrier. An independent known
/// circle is preferred when available; otherwise propose and prove one from
/// the actual rational curve. Parameter orientation is recovered from the
/// global angular proof, not derivative samples.
pub(crate) fn normalize_circular_edge(
    edge: &OrientedEdgeUse,
    expected: Option<CircleSupport>,
) -> Option<OrientedEdgeUse> {
    let CurveSupport::BSpline(spline) = &edge.support else {
        return matches!(edge.support, CurveSupport::Circle(_)).then(|| edge.clone());
    };
    let curve = spline.rational_curve()?;
    let circle = expected.or_else(|| propose_spline_circle(&curve, GEOM_TOL_MM))?;
    let front = curve.front();
    let back = curve.back();
    let point = |p: Point3| [p.x, p.y, p.z];
    let (start, end) = if edge.parameter_forward {
        (front, back)
    } else {
        (back, front)
    };
    if distance(point(start), edge.start_mm) > GEOM_TOL_MM
        || distance(point(end), edge.end_mm) > GEOM_TOL_MM
    {
        return None;
    }
    let frame = CircularArcFrame::new(
        Point3::from(circle.center_mm),
        Vector3::from(circle.normal),
        circle.radius_mm,
        front,
        GEOM_TOL_MM,
    )?;
    let coverage = certify_circular_arc_coverage(&curve, frame)?;
    let circle = CircleSupport {
        x_direction: normalize(sub(point(front), circle.center_mm), 0.0)?,
        ..circle
    };
    let mut normalized = edge.clone();
    normalized.support = CurveSupport::Circle(circle);
    normalized.parameter_forward = coverage.forward == edge.parameter_forward;
    let sweep = circle_edge_sweep(&normalized, circle)?;
    let angle_tolerance = TOL.source_sweep_angle_rad.max(frame.angle_tolerance());
    ((sweep.abs() - coverage.sweep_rad).abs() <= angle_tolerance).then_some(normalized)
}
