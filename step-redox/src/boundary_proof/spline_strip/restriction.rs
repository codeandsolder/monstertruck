//! Whole-span restriction of certified rational Bernstein cap profiles.
//! Partition vertices are parameter proposals; complete homogeneous
//! coefficient comparison and ordered coverage are the geometric evidence.
use super::{
    BernsteinSpanKey, CurveEvidence, POINT_Q_MM, PointKey, SplineCurveKey, TOL, TrimEndpoint,
    decode, point_key,
};
use crate::math3::{distance, dot, mul, normalize, sub};
use monstertruck_geometry::prelude::Vector4;
use num_traits::ToPrimitive;
use std::collections::BTreeMap;

/// A position on the canonical parent, with its actual topological point.
/// A knot belongs to the following span at fraction zero; this makes a
/// parent knot and an edge subdivision join identically.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(in crate::boundary_proof) struct ProfilePosition {
    span: usize,
    fraction: u64,
    pub(super) point: PointKey,
}

impl ProfilePosition {
    const fn new(span: usize, fraction: f64, point: PointKey) -> Self {
        if fraction.to_bits() == 1.0_f64.to_bits() {
            Self {
                span: span + 1,
                fraction: 0,
                point,
            }
        } else {
            Self {
                span,
                fraction: fraction.to_bits(),
                point,
            }
        }
    }
}

fn blossom(controls: &[Vector4], lo: f64, hi: f64, right: usize) -> Vector4 {
    let degree = controls.len() - 1;
    let mut work = controls.to_vec();
    for round in 0..degree {
        let t = if round < degree - right { lo } else { hi };
        for i in 0..degree - round {
            work[i] = work[i] * (1. - t) + work[i + 1] * t;
        }
    }
    work[0]
}

/// Every coefficient of a restricted rational Bernstein span comes from
/// blossom evaluation with left/right-repeated interval parameters.
pub(super) fn restrict_span(controls: &[Vector4], lo: f64, hi: f64) -> Vec<Vector4> {
    (0..controls.len())
        .map(|right| blossom(controls, lo, hi, right))
        .collect()
}

fn span_controls(span: &BernsteinSpanKey) -> Option<Vec<Vector4>> {
    span.points_mm
        .iter()
        .zip(&span.weights)
        .map(|(&point, &weight)| {
            let [x, y, z] = decode(point)?;
            let w = f64::from_bits(weight);
            Some(Vector4::new(x * w, y * w, z * w, w))
        })
        .collect()
}

fn projected(p: Vector4) -> [f64; 3] { [p.x / p.w, p.y / p.w, p.z / p.w] }

/// Locate endpoints on a monotone Bernstein span. The coordinate is only a
/// parameter proposal; full restricted coefficient equality below is the
/// certificate. Identical endpoint keys always obtain the same parameter.
fn span_parameter(controls: &[Vector4], direction: [f64; 3], target: PointKey) -> Option<f64> {
    let start = dot(projected(*controls.first()?), direction);
    let end = dot(projected(*controls.last()?), direction);
    let value = dot(decode(target)?, direction);
    if value < start - TOL || value > end + TOL {
        return None;
    }
    if target == point_key(projected(*controls.first()?))? {
        return Some(0.);
    }
    if target == point_key(projected(*controls.last()?))? {
        return Some(1.);
    }
    if value <= start || value >= end {
        return None;
    }
    let (mut lo, mut hi) = (0., 1.);
    for _ in 0..48 {
        let t = f64::midpoint(lo, hi);
        let p = blossom(controls, t, t, 0);
        if dot(projected(p), direction) < value {
            lo = t;
        } else {
            hi = t;
        }
    }
    Some(f64::midpoint(lo, hi))
}

/// Convex-hull error between equal-degree rational Bernstein nets under
/// a common affine chart. Both nets have positive weights. Remove their
/// independent homogeneous scales, then bound the change in the normalized
/// basis using the minimum relative weight and the control-hull radius.
/// The actual (unrounded) child net is compared with the canonical parent.
/// Reserve its extraction and subdivision arithmetic error budgets.
fn same_restricted_span(a: &[Vector4], b: &[Vector4]) -> bool {
    if a.len() != b.len() || a.is_empty() {
        return false;
    }
    let scale = |c: &[Vector4]| c.iter().map(|p| p.w).fold(0., f64::max);
    let (sa, sb) = (scale(a), scale(b));
    let mut weight_min = 1.0_f64;
    let mut weight_error = 0.0_f64;
    let mut point_error = 0.0_f64;
    let mut radius = 0.0_f64;
    let origin = projected(b[0]);
    for (&p, &q) in a.iter().zip(b) {
        let (w, v) = (p.w / sa, q.w / sb);
        weight_min = weight_min.min(w).min(v);
        weight_error = weight_error.max((w - v).abs());
        point_error = point_error.max(distance(projected(p), projected(q)));
        radius = radius.max(distance(projected(q), origin));
    }
    let Some(order) = a.len().to_f64() else {
        return false;
    };
    let coordinate_scale = a
        .iter()
        .chain(b)
        .map(|&p| crate::math3::norm(projected(p)))
        .fold(1., f64::max);
    let subdivision_error =
        32. * order * order * f64::EPSILON * coordinate_scale / weight_min.powi(2);
    let input_error = POINT_Q_MM / 4.;
    let total =
        point_error + 2. * radius * weight_error / weight_min + subdivision_error + input_error;
    weight_min > 0. && total.is_finite() && total <= TOL / 2.
}

/// Pick one deterministic monotone coordinate for the parent span, never
/// a different parameter estimator for each incident piece. The largest
/// chord projection avoids poorly conditioned nearly-stationary axes.
fn span_direction(controls: &[Vector4]) -> Option<[f64; 3]> {
    let chord = sub(projected(*controls.last()?), projected(*controls.first()?));
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], chord]
        .into_iter()
        .filter_map(|d| normalize(d, TOL))
        .map(|axis| {
            if dot(chord, axis) < 0. {
                mul(axis, -1.)
            } else {
                axis
            }
        })
        .filter(|&d| {
            dot(chord, d) > TOL
                && controls
                    .windows(2)
                    .all(|p| dot(sub(projected(p[1]), projected(p[0])), d) >= 0.)
        })
        .max_by(|&a, &b| dot(chord, a).total_cmp(&dot(chord, b)))
}

fn restricted_span(
    parent: &BernsteinSpanKey,
    piece: &BernsteinSpanKey,
    actual: &[Vector4],
) -> Option<[f64; 2]> {
    if parent == piece {
        return Some([0., 1.]);
    }
    if parent.points_mm.len() != actual.len() {
        return None;
    }
    for axis in 0..3 {
        let minimum = parent.points_mm.iter().map(|p| p[axis]).min()?;
        let maximum = parent.points_mm.iter().map(|p| p[axis]).max()?;
        if [piece.points_mm.first()?, piece.points_mm.last()?]
            .into_iter()
            .any(|p| {
                let delta = p[axis].abs_diff(p[axis].clamp(minimum, maximum));
                delta
                    .to_f64()
                    .is_none_or(|distance| distance * POINT_Q_MM > TOL)
            })
        {
            return None;
        }
    }
    let controls = span_controls(parent)?;
    if controls.len() != actual.len() {
        return None;
    }
    let direction = span_direction(&controls)?;
    let lo = span_parameter(&controls, direction, *piece.points_mm.first()?)?;
    let hi = span_parameter(&controls, direction, *piece.points_mm.last()?)?;
    (lo < hi && same_restricted_span(&restrict_span(&controls, lo, hi), actual)).then_some([lo, hi])
}

/// A restriction may subdivide any individual original knot span. Its
/// own span boundaries must join at identical parent coordinates, not
/// merely happen to have nearby endpoints in world space.
fn restriction_in_order(
    parent: &SplineCurveKey,
    piece: &CurveEvidence,
    ends: [TrimEndpoint; 2],
) -> Option<[ProfilePosition; 2]> {
    let mut interval: Option<[ProfilePosition; 2]> = None;
    for (j, (span, actual)) in piece.key.spans.iter().zip(&piece.controls).enumerate() {
        let points = [
            if j == 0 {
                ends[0]
            } else {
                TrimEndpoint::new(projected(*actual.first()?))?
            },
            if j + 1 == piece.key.spans.len() {
                ends[1]
            } else {
                TrimEndpoint::new(projected(*actual.last()?))?
            },
        ];
        let mut matches = parent
            .spans
            .iter()
            .enumerate()
            // Subsequent spans must begin exactly at the preceding
            // normalized parent position, so only one span can qualify.
            .filter(|(i, _)| interval.is_none_or(|prior| prior[1].span == *i))
            .filter_map(|(i, candidate)| {
                let [lo, hi] = restricted_span(candidate, span, actual)?;
                // Topological vertices may differ slightly from the spline
                // trim's evaluated endpoints, but must still be on the
                // parent at their certified parameter position.
                let expected = [*span.points_mm.first()?, *span.points_mm.last()?];
                if points.map(|point| point.key) != expected {
                    let controls = span_controls(candidate)?;
                    for (t, point) in [lo, hi].into_iter().zip(points) {
                        let error = distance(projected(blossom(&controls, t, t, 0)), point.mm);
                        if !error.is_finite() || error > TOL {
                            return None;
                        }
                    }
                }
                Some([
                    ProfilePosition::new(i, lo, points[0].key),
                    ProfilePosition::new(i, hi, points[1].key),
                ])
            });
        let range = matches.next()?;
        if matches.next().is_some() {
            return None;
        }
        if let Some(accumulated) = &mut interval {
            if accumulated[1] != range[0] {
                return None;
            }
            accumulated[1] = range[1];
        } else {
            interval = Some(range);
        }
    }
    interval
}

fn full_profile_interval(parent: &SplineCurveKey) -> Option<[ProfilePosition; 2]> {
    let [start, end] = parent.ends()?;
    Some([
        ProfilePosition::new(0, 0., start),
        ProfilePosition::new(parent.spans.len(), 0., end),
    ])
}

pub(super) fn profile_piece(
    parent: &SplineCurveKey,
    piece: &CurveEvidence,
    mut ends: [TrimEndpoint; 2],
) -> Option<[ProfilePosition; 2]> {
    if *parent == piece.key && parent.ends()? == ends.map(|point| point.key) {
        return full_profile_interval(parent);
    }
    if let Some(interval) = restriction_in_order(parent, piece, ends) {
        return Some(interval);
    }
    let mut reversed = piece.clone();
    reversed.key.spans.reverse();
    reversed.controls.reverse();
    for (span, controls) in reversed.key.spans.iter_mut().zip(&mut reversed.controls) {
        span.points_mm.reverse();
        span.weights.reverse();
        controls.reverse();
    }
    ends.reverse();
    restriction_in_order(parent, &reversed, ends)
}

/// Check actual carrier incidence without charging a fictitious rounding
/// error for vertices whose original coordinates remain available.
pub(super) fn raw_parent_incidence(
    parent: &CurveEvidence,
    interval: [ProfilePosition; 2],
    mut ends: [TrimEndpoint; 2],
) -> Option<()> {
    if interval[0].point != ends[0].key {
        ends.reverse();
    }
    for (position, endpoint) in interval.into_iter().zip(ends) {
        if position.point != endpoint.key {
            return None;
        }
        let (span, t) = if position.span == parent.controls.len() {
            (parent.controls.last()?, 1.)
        } else {
            (
                parent.controls.get(position.span)?,
                f64::from_bits(position.fraction),
            )
        };
        let weight_max = span.iter().map(|p| p.w).fold(0., f64::max);
        let weight_min = span.iter().map(|p| p.w / weight_max).fold(1., f64::min);
        let order = span.len().to_f64()?;
        let coordinate_scale = span
            .iter()
            .map(|&p| crate::math3::norm(projected(p)))
            .fold(1., f64::max);
        let arithmetic = 32. * order * order * f64::EPSILON * coordinate_scale / weight_min.powi(2);
        let error = distance(projected(blossom(span, t, t, 0)), endpoint.mm)
            + arithmetic
            + parent.extraction_error_mm;
        if !error.is_finite() || error > TOL {
            return None;
        }
    }
    Some(())
}

/// Exact, ordered coverage of the complete canonical parent. Each
/// adjacent pair must agree on both parameter and topological endpoint.
/// Overlaps, skipped knot intervals and duplicate traversals refuse.
pub(in crate::boundary_proof) fn complete_profile_cover(
    parent: &SplineCurveKey,
    pieces: &[[ProfilePosition; 2]],
) -> Option<()> {
    if pieces.iter().any(|[a, b]| a >= b) {
        return None;
    }
    complete_directed_cover(full_profile_interval(parent)?, pieces)
}

/// Each directed piece has an already certified geometric locus. Establish
/// connected, exactly-once topological coverage of the two parent endpoints.
pub(super) fn complete_directed_cover<T: Copy + Ord>(
    ends: [T; 2],
    pieces: &[[T; 2]],
) -> Option<()> {
    let [mut point, end] = ends;
    let mut outgoing = BTreeMap::new();
    for &[a, b] in pieces {
        if a == b || outgoing.insert(a, b).is_some() {
            return None;
        }
    }
    for _ in 0..pieces.len() {
        if point == end {
            return None;
        }
        point = outgoing.remove(&point)?;
    }
    (point == end && outgoing.is_empty()).then_some(())
}

#[cfg(test)]
mod tests {
    use super::super::curve_evidence;
    use super::*;
    use crate::brep::BSplineSupport;
    use anyhow::{Context, Result};

    fn restricted(controls: &[Vector4], lo: f64, hi: f64) -> BSplineSupport {
        let controls = restrict_span(controls, lo, hi);
        let order = controls.len();
        BSplineSupport {
            degree: order - 1,
            control_points_mm: controls.iter().copied().map(projected).collect(),
            knots: std::iter::repeat_n(-2., order)
                .chain(std::iter::repeat_n(7., order))
                .collect(),
            weights: Some(controls.iter().map(|p| p.w).collect()),
        }
    }

    #[test]
    fn raw_trim_incidence_keeps_the_actual_geometric_budget() -> Result<()> {
        let original = BSplineSupport {
            degree: 3,
            control_points_mm: vec![
                [0., 0., 0.],
                [0.25, 0.25, 0.],
                [0.75, 0.25, 0.],
                [1., 0., 0.],
            ],
            knots: vec![0., 0., 0., 0., 1., 1., 1., 1.],
            weights: None,
        };
        let parent = curve_evidence(&original).context("parent")?;
        let child = curve_evidence(&restricted(&parent.controls[0], 0., 0.5)).context("child")?;
        let start = TrimEndpoint::new(projected(child.controls[0][0])).context("start")?;
        let mut split = projected(*child.controls[0].last().context("split")?);
        // The authored point remains within the tolerance, but the old
        // reconstructed-coordinate + half-cell-diagonal bound rejects it.
        split[2] = 0.9 * TOL;
        let end = TrimEndpoint::new(split).context("end")?;
        let reconstructed = decode(end.key).context("rounded point")?;
        assert!(
            distance(
                reconstructed,
                projected(*child.controls[0].last().context("split")?)
            ) + POINT_Q_MM * 3.0_f64.sqrt() / 2.
                > TOL
        );
        let interval = profile_piece(&parent.key, &child, [start, end]).context("raw incidence")?;
        assert!(raw_parent_incidence(&parent, interval, [start, end]).is_some());
        split[2] = 1.05 * TOL;
        let outside = TrimEndpoint::new(split).context("outside")?;
        assert!(profile_piece(&parent.key, &child, [start, outside]).is_none());
        assert!(raw_parent_incidence(&parent, interval, [start, outside]).is_none());
        Ok(())
    }

    #[test]
    fn rational_multi_span_partitions_cover_once_and_refuse_wrong_interior() -> Result<()> {
        let original = BSplineSupport {
            degree: 2,
            control_points_mm: vec![[0., 0., 0.], [0.25, 0.5, 0.], [0.75, 0.5, 0.], [1., 0., 0.]],
            knots: vec![0., 0., 0., 0.5, 1., 1., 1.],
            weights: Some(vec![1., 2., 2., 1.]),
        };
        let parent = curve_evidence(&original).context("rational parent")?;
        let children = parent
            .controls
            .iter()
            .flat_map(|span| [(0., 0.375), (0.375, 1.)].map(|(lo, hi)| restricted(span, lo, hi)))
            .collect::<Vec<_>>();
        let interval = |curve: &BSplineSupport| {
            let piece = curve_evidence(curve)?;
            profile_piece(
                &parent.key,
                &piece,
                [
                    TrimEndpoint::new(projected(*piece.controls.first()?.first()?))?,
                    TrimEndpoint::new(projected(*piece.controls.last()?.last()?))?,
                ],
            )
        };
        let pieces = children
            .iter()
            .map(interval)
            .collect::<Option<Vec<_>>>()
            .context("all restrictions")?;
        assert_eq!(pieces.len(), 4);
        assert!(complete_profile_cover(&parent.key, &pieces).is_some());
        assert!(
            complete_profile_cover(&parent.key, &[pieces[3], pieces[1], pieces[0], pieces[2]])
                .is_some()
        );
        assert!(complete_profile_cover(&parent.key, &pieces[..3]).is_none());
        assert!(
            complete_profile_cover(
                &parent.key,
                &[pieces[0], pieces[1], pieces[1], pieces[2], pieces[3]]
            )
            .is_none()
        );
        let mut gap = pieces.clone();
        gap[1][0].fraction += 1;
        assert!(complete_profile_cover(&parent.key, &gap).is_none());
        let mut wrong = children[0].clone();
        wrong.control_points_mm[1][1] += 1e-4;
        assert!(interval(&wrong).is_none());
        let mut wrong = children[0].clone();
        wrong.weights.as_mut().context("rational weights")?[1] *= 1.1;
        assert!(interval(&wrong).is_none());
        let mut reversed = children[2].clone();
        reversed.control_points_mm.reverse();
        reversed
            .weights
            .as_mut()
            .context("rational weights")?
            .reverse();
        assert_eq!(interval(&reversed), Some(pieces[2]));
        Ok(())
    }
}
