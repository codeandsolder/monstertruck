//! Geometric sample fitting shared across mesh-recovery hypothesis generators.

use crate::{InnerSpace, Point3, Vector3};

/// Least-squares circle in a known plane. The relative frame avoids
/// cancellation for large world offsets.
pub(super) fn fit_circle(points: &[Point3], ids: &[usize], axis: Vector3) -> Option<(Point3, f64)> {
    let start = points[*ids.first()?];
    let delta = points[*ids.get(1)?] - start;
    let projected = delta - axis * delta.dot(axis);
    let norm = projected.magnitude();
    if !norm.is_finite() || norm <= 0.0 {
        return None;
    }
    let xaxis = projected / norm;
    let yaxis = axis.cross(xaxis);
    let samples = ids
        .iter()
        .map(|&id| {
            let delta = points[id] - start;
            (delta.dot(xaxis), delta.dot(yaxis))
        })
        .collect::<Vec<_>>();
    let n = crate::usize_to_f64(samples.len());
    let mx = samples.iter().map(|&(x, _)| x).sum::<f64>() / n;
    let my = samples.iter().map(|&(_, y)| y).sum::<f64>() / n;
    let mqq = samples
        .iter()
        .map(|&(x, y)| x.mul_add(x, y * y))
        .sum::<f64>()
        / n;
    let (mut xx, mut xy, mut yy, mut qx, mut qy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y) in &samples {
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
    let center = start + xaxis * cx + yaxis * cy;
    let radius = samples
        .iter()
        .map(|&(x, y)| (x - cx).hypot(y - cy))
        .sum::<f64>()
        / n;
    (radius > 0.0 && radius.is_finite()).then_some((center, radius))
}

pub(super) fn normal_basis(axis: Vector3) -> (Vector3, Vector3) {
    // Select a transverse reference away from the axis for conditioning.
    let reference = if axis.z.abs() <= 0.75 {
        Vector3::unit_z()
    } else {
        Vector3::unit_x()
    };
    let u = axis.cross(reference).normalize();
    let v = axis.cross(u);
    (u, v)
}
