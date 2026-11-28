//! Geometric sample fitting shared across mesh-recovery hypothesis generators.

use super::evidence::canonical_unit;
use crate::{Point3, Vector3};

/// Least-squares circle in a known plane. The relative frame avoids
/// cancellation for large world offsets.
pub(super) fn fit_circle(points: &[Point3], ids: &[usize], axis: Vector3) -> Option<(Point3, f64)> {
    let start = points[*ids.first()?];
    let delta = points[*ids.get(1)?] - start;
    let projected = delta - axis * delta.dot(&axis);
    let norm = projected.norm();
    if !norm.is_finite() || norm <= 0.0 {
        return None;
    }
    let xaxis = projected / norm;
    let yaxis = axis.cross(&xaxis);
    let samples = ids
        .iter()
        .map(|&id| {
            let delta = points[id] - start;
            (delta.dot(&xaxis), delta.dot(&yaxis))
        })
        .collect::<Vec<_>>();
    let ([cx, cy], radius) = fit_circle_2d(&samples)?;
    Some((start + xaxis * cx + yaxis * cy, radius))
}

/// Algebraic least-squares circle in a source-local two-dimensional frame.
pub(super) fn fit_circle_2d(samples: &[(f64, f64)]) -> Option<([f64; 2], f64)> {
    if samples.len() < 3 {
        return None;
    }
    let n = crate::usize_to_f64(samples.len());
    let mx = samples.iter().map(|&(x, _)| x).sum::<f64>() / n;
    let my = samples.iter().map(|&(_, y)| y).sum::<f64>() / n;
    let mqq = samples
        .iter()
        .map(|&(x, y)| x.mul_add(x, y * y))
        .sum::<f64>()
        / n;
    let (mut xx, mut xy, mut yy, mut qx, mut qy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y) in samples {
        let dx = x - mx;
        let dy = y - my;
        let dq = x.mul_add(x, y * y) - mqq;
        xx = dx.mul_add(dx, xx);
        xy = dx.mul_add(dy, xy);
        yy = dy.mul_add(dy, yy);
        qx = (dx * 0.5).mul_add(dq, qx);
        qy = (dy * 0.5).mul_add(dq, qy);
    }
    let det = xx.mul_add(yy, -(xy * xy));
    if !det.is_finite() || det <= f64::EPSILON * (xx + yy).powi(2) {
        return None;
    }
    let cx = (qx.mul_add(yy, -(qy * xy))) / det;
    let cy = (qy.mul_add(xx, -(qx * xy))) / det;
    let radius = samples
        .iter()
        .map(|&(x, y)| (x - cx).hypot(y - cy))
        .sum::<f64>()
        / n;
    (radius > 0.0 && radius.is_finite()).then_some(([cx, cy], radius))
}

pub(super) fn normal_basis(axis: Vector3) -> (Vector3, Vector3) {
    // Select a transverse reference away from the axis for conditioning.
    let reference = if axis.z.abs() <= 0.75 {
        Vector3::z()
    } else {
        Vector3::x()
    };
    let u = axis.cross(&reference).normalize();
    let v = axis.cross(&u);
    (u, v)
}

/// Placed circle fitted to original samples without requiring a path order.
#[derive(Clone, Copy)]
pub(super) struct CircleCarrierFit {
    pub center: Point3,
    pub normal: Vector3,
    pub radius: f64,
    pub plane_span_strength: f64,
    pub source_span: f64,
    pub max_plane_residual: f64,
    pub max_radial_residual: f64,
}

fn source_plane(points: &[Point3], ids: &[usize]) -> Option<(Vector3, f64, usize, f64)> {
    let origin = points[*ids.first()?];
    let (second, delta) = ids
        .iter()
        .enumerate()
        .map(|(i, &id)| (i, points[id] - origin))
        .max_by(|a, b| a.1.norm_squared().total_cmp(&b.1.norm_squared()))?;
    let span_squared = delta.norm_squared();
    if !span_squared.is_finite() || span_squared <= 0.0 {
        return None;
    }
    let cross = ids
        .iter()
        .map(|&id| delta.cross(&(points[id] - origin)))
        .max_by(|a, b| a.norm_squared().total_cmp(&b.norm_squared()))?;
    let strength = cross.norm() / span_squared;
    if !strength.is_finite() || strength <= 64.0 * f64::EPSILON {
        return None;
    }
    Some((canonical_unit(cross), strength, second, span_squared.sqrt()))
}

pub(super) fn fit_circle_carrier(points: &[Point3], ids: &[usize]) -> Option<CircleCarrierFit> {
    if ids.len() < 3 {
        return None;
    }
    let (normal, plane_span_strength, second, source_span) = source_plane(points, ids)?;
    let mut fit_ids = ids.to_vec();
    fit_ids.swap(1, second);
    let (center, radius) = fit_circle(points, &fit_ids, normal)?;
    let (max_plane_residual, max_radial_residual) =
        ids.iter().fold((0.0_f64, 0.0_f64), |(plane, radial), &id| {
            let relative = points[id] - center;
            let distance = relative.dot(&normal);
            let projected = relative - normal * distance;
            (
                plane.max(distance.abs()),
                radial.max((projected.norm() - radius).abs()),
            )
        });
    Some(CircleCarrierFit {
        center,
        normal,
        radius,
        plane_span_strength,
        source_span,
        max_plane_residual,
        max_radial_residual,
    })
}
