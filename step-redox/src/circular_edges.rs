//! Representation-independent bounded circle trims.
//! A rational edge becomes a circular edge only after complete incidence,
//! endpoint and ordered angular coverage have all been certified.

use crate::brep::{CircleSupport, CurveSupport, OrientedEdgeUse};
use crate::math3::{cross, distance, dot, norm, normalize, sub};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES as TOL;
use monstertruck_geometry::{
    analytic_curve_incidence::{CircularArcFrame, certify_circular_arc_coverage, propose_circle},
    prelude::{BoundedCurve, NurbsCurve, Point3, Vector3, Vector4},
};
use std::f64::consts::TAU;

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

pub(crate) fn circle_edge_sweep(edge: &OrientedEdgeUse, circle: CircleSupport) -> Option<f64> {
    let start = circle_angle_3d(edge.start_mm, circle)?;
    let end = circle_angle_3d(edge.end_mm, circle)?;
    let closed = distance(edge.start_mm, edge.end_mm) <= GEOM_TOL_MM;
    Some(directional_sweep(
        start,
        end,
        edge.parameter_forward,
        closed,
    ))
}

pub(crate) fn directional_sweep(start: f64, end: f64, increasing: bool, closed: bool) -> f64 {
    if closed {
        return if increasing { TAU } else { -TAU };
    }
    if increasing {
        (end - start).rem_euclid(TAU)
    } else {
        -((start - end).rem_euclid(TAU))
    }
}

pub(crate) fn circle_angle_3d(point: [f64; 3], circle: CircleSupport) -> Option<f64> {
    let y_direction = normalize(
        cross(circle.normal, circle.x_direction),
        TOL.solid_direction,
    )?;
    let relative = sub(point, circle.center_mm);
    if dot(relative, circle.normal).abs() > GEOM_TOL_MM {
        return None;
    }
    let radial = norm(relative);
    if (radial - circle.radius_mm).abs() > GEOM_TOL_MM {
        return None;
    }
    Some(dot(relative, y_direction).atan2(dot(relative, circle.x_direction)))
}
