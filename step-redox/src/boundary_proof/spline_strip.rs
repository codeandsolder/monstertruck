//! Full-domain spline-strip certificates. Neither fitting sampled points nor
//! matching endpoint pairs establishes a curved edge. Retain complete positive-
//! weight rational Bernstein control nets, ordered spans, and actual trims.
use super::{
    ANALYTIC_BOUNDARY_PROOF_TOL_MM as TOL, BSplineSupport, CurveSupport, DIRECTION_Q, EdgeGeometry,
    OrientedSupportKey, POINT_Q_MM, PointKey, SupportLocusKey, bounded, dequantize_scalar,
    dequantize_vector, point_key,
};
use crate::{
    brep,
    math3::{distance, dot, mul, normalize, sub},
};
use monstertruck_geometry::{
    analytic_line_segment::{SegmentLocus, certify_nurbs_segment_locus},
    prelude::Point3,
    proof_bezier_spans,
};
use num_traits::ToPrimitive;
use serde::Serialize;
use std::collections::HashMap;
use step_p21::ast::EntityInstance;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct BernsteinSpanKey {
    points_mm: Vec<PointKey>,
    // Positive weights with a common homogeneous scale removed. Do not
    // quantize weights or silently collapse unequal rational basis ratios.
    weights: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub struct SplineCurveKey {
    spans: Vec<BernsteinSpanKey>,
}

impl SplineCurveKey {
    fn ends(&self) -> Option<[PointKey; 2]> {
        Some([
            *self.spans.first()?.points_mm.first()?,
            *self.spans.last()?.points_mm.last()?,
        ])
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub struct SplineStripKey {
    profiles: [SplineCurveKey; 2],
}

fn decode(p: PointKey) -> Option<[f64; 3]> { dequantize_vector(p, POINT_Q_MM) }

fn weight_ratios(weights: &[f64]) -> Option<Vec<u64>> {
    if weights.iter().any(|w| !w.is_finite() || *w <= 0.) {
        return None;
    }
    let scale = weights.iter().copied().reduce(f64::max)?;
    weights
        .iter()
        .map(|w| {
            let ratio = w / scale;
            (ratio.is_finite() && ratio > 0.).then_some(ratio.to_bits())
        })
        .collect()
}

/// Bound numerical conditioning before using rounded Bernstein control keys.
/// Each blossom round makes at most sixteen elementary operations including
/// knot normalization and homogeneous interpolation. Positive convex
/// combinations do not amplify earlier absolute coefficient errors. The
/// reciprocal minimum normalized knot spacing bounds perturbations of an
/// interpolation fraction; two reciprocal minimum relative weights bound
/// dehomogenization and rational evaluation. Reject rather than claim a
/// sub-quantum geometric certificate for very ill-conditioned input.
fn well_conditioned(spline: &BSplineSupport) -> Option<()> {
    let (lo, hi) = super::bspline_domain(spline)?;
    let n = spline.control_points_mm.len();
    let order = spline.degree.checked_add(1)?;
    if n < order
        || spline.knots[..order]
            .iter()
            .any(|k| k.total_cmp(&lo).is_ne())
        || spline.knots[n..].iter().any(|k| k.total_cmp(&hi).is_ne())
    {
        return None;
    }
    let min_span = spline
        .knots
        .windows(2)
        .filter_map(|w| (w[1] > w[0]).then_some((w[1] - w[0]) / (hi - lo)))
        .reduce(f64::min)?;
    let coordinate_scale = spline
        .control_points_mm
        .iter()
        .flatten()
        .map(|x| x.abs())
        .reduce(f64::max)?
        .max(1.);
    let weight_ratio = if let Some(weights) = &spline.weights {
        weights.iter().copied().reduce(f64::min)? / weights.iter().copied().reduce(f64::max)?
    } else {
        1.
    };
    let operation_error = 16. * order.to_f64()? * f64::EPSILON / min_span;
    let error_mm = 2. * coordinate_scale * operation_error / weight_ratio.powi(2);
    (error_mm.is_finite() && error_mm <= POINT_Q_MM / 4.).then_some(())
}

/// Canonicalize the whole curve, allowing affine reversal and different
/// knot interval lengths. Each positive span retains every homogeneous
/// Bernstein coefficient, with its own common weight scale removed.
fn curve_key(spline: &BSplineSupport) -> Option<(SplineCurveKey, bool)> {
    well_conditioned(spline)?;
    let rational = bounded::rational(spline)?;
    let mut spans = proof_bezier_spans(rational.non_rationalized())?;
    let start = point_key(*spline.control_points_mm.first()?)?;
    let end = point_key(*spline.control_points_mm.last()?)?;
    if start == end {
        return None;
    }
    let forward = start < end;
    if !forward {
        spans.reverse();
        for span in &mut spans {
            span.controls.reverse();
        }
    }
    let polynomial = spline
        .weights
        .as_ref()
        .is_none_or(|weights| weights.iter().all(|w| w.to_bits() == weights[0].to_bits()));
    let spans = spans
        .into_iter()
        .map(|span| {
            let points_mm = span
                .controls
                .iter()
                .map(|p| point_key([p[0] / p[3], p[1] / p[3], p[2] / p[3]]))
                .collect::<Option<_>>()?;
            let weights = if polynomial {
                vec![1.0_f64.to_bits(); span.controls.len()]
            } else {
                weight_ratios(&span.controls.iter().map(|p| p[3]).collect::<Vec<_>>())?
            };
            Some(BernsteinSpanKey { points_mm, weights })
        })
        .collect::<Option<Vec<_>>>()?;
    let key = SplineCurveKey { spans };
    if key.ends()? != if forward { [start, end] } else { [end, start] } {
        return None;
    }
    // Discontinuous pieces may not be admitted as one complete edge.
    if key
        .spans
        .windows(2)
        .any(|pair| pair[0].points_mm.last() != pair[1].points_mm.first())
    {
        return None;
    }
    Some((key, forward))
}

/// Establish monotonicity of the canonical Bernstein locus, not of an arbitrary
/// source control polygon. The same rounded control keys establish geometric
/// equivalence; exact positive weight ratios retain the rational basis.
/// This permits equal coordinates in the proof grid without introducing a
/// second tolerance for almost-monotone control polygons.
fn monotone_profile(key: &SplineCurveKey, direction: [f64; 3]) -> bool {
    key.spans.iter().all(|span| {
        let Some(points) = span
            .points_mm
            .iter()
            .map(|&p| decode(p))
            .collect::<Option<Vec<_>>>()
        else {
            return false;
        };
        let increments = points
            .windows(2)
            .map(|pair| dot(sub(pair[1], pair[0]), direction));
        let (minimum, maximum) = increments
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
                (lo.min(d), hi.max(d))
            });
        minimum >= 0. && maximum > TOL
    })
}

type DirectedProfile = (SplineCurveKey, bool);

fn translated_profiles(profiles: &[BSplineSupport; 2]) -> Option<[DirectedProfile; 2]> {
    let [first, second] = profiles;
    if first.degree != second.degree
        || first.knots != second.knots
        || first.control_points_mm.len() != second.control_points_mm.len()
    {
        return None;
    }
    let weights = |s: &BSplineSupport| {
        weight_ratios(
            &s.weights
                .clone()
                .unwrap_or_else(|| vec![1.; s.control_points_mm.len()]),
        )
    };
    if weights(first)? != weights(second)? {
        return None;
    }
    let delta = sub(
        *second.control_points_mm.first()?,
        *first.control_points_mm.first()?,
    );
    let axis = normalize(delta, TOL)?;
    if first
        .control_points_mm
        .iter()
        .zip(&second.control_points_mm)
        .any(|(&a, &b)| distance(sub(b, a), delta) > TOL)
    {
        return None;
    }
    let chord = sub(
        *first.control_points_mm.last()?,
        *first.control_points_mm.first()?,
    );
    let candidates = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], chord];
    let keys = [curve_key(first)?, curve_key(second)?];
    candidates
        .into_iter()
        .filter_map(|e| normalize(sub(e, mul(axis, dot(axis, e))), TOL))
        .flat_map(|d| [d, mul(d, -1.)])
        .any(|d| keys.iter().all(|(key, _)| monotone_profile(key, d)))
        .then_some(keys)
}

fn curve_matches_edge(spline: &BSplineSupport, edge: &EdgeGeometry) -> Option<SplineCurveKey> {
    let (key, forward) = curve_key(spline)?;
    let [start, end] = key.ends()?;
    let actual = if edge.parameter_forward == forward {
        [edge.start, edge.end]
    } else {
        [edge.end, edge.start]
    };
    (actual == [start, end]).then_some(key)
}

fn segment_matches_edge(curve: &CurveSupport, ends: [PointKey; 2]) -> Option<()> {
    let [a, b] = ends.map(decode);
    let (a, b) = (a?, b?);
    match curve {
        CurveSupport::Line(line) => {
            (bounded::on_line(a, *line) && bounded::on_line(b, *line)).then_some(())
        }
        CurveSupport::BSpline(spline) => {
            let rational = bounded::rational(spline)?;
            (certify_nurbs_segment_locus(
                &rational,
                Point3::new(a[0], a[1], a[2]),
                Point3::new(b[0], b[1], b[2]),
                TOL,
            ) == SegmentLocus::WithinTolerance)
                .then_some(())
        }
        _ => None,
    }
}

/// Verify one complete rectangle in parameter space. The physical boundary
/// must traverse both complete cap curves and both end generators exactly
/// once in the sense of its oriented surface; holes, partial trims, and
/// unrelated curves sharing only endpoints cannot establish this certificate.
fn full_trim(
    profiles: &[SplineCurveKey; 2],
    forward: bool,
    loops: &[brep::FaceLoop],
) -> Option<()> {
    let [boundary] = loops else {
        return None;
    };
    if boundary.edges.len() != 4 {
        return None;
    }
    let [a, b] = profiles[0].ends()?;
    let [c, d] = profiles[1].ends()?;
    let mut expected = [[a, b], [b, d], [d, c], [c, a]];
    if forward != boundary.orientation {
        for ends in &mut expected {
            ends.reverse();
        }
    }
    let mut seen = [false; 4];
    for edge in &boundary.edges {
        let ends = [point_key(edge.start_mm)?, point_key(edge.end_mm)?];
        let i = expected.iter().position(|e| *e == ends)?;
        if std::mem::replace(&mut seen[i], true) {
            return None;
        }
        if i % 2 == 1 {
            segment_matches_edge(&edge.support, ends)?;
        } else {
            let CurveSupport::BSpline(spline) = &edge.support else {
                return None;
            };
            let geometry = EdgeGeometry {
                start: ends[0],
                end: ends[1],
                start_vertex: edge.start_vertex,
                end_vertex: edge.end_vertex,
                parameter_forward: edge.parameter_forward,
                curve_id: edge.curve_id,
                support: edge.support.clone(),
            };
            if curve_matches_edge(spline, &geometry)? != profiles[i / 2] {
                return None;
            }
        }
    }
    seen.into_iter().all(|s| s).then_some(())
}

pub(super) fn oriented_key(
    face_id: u64,
    surface_id: u64,
    forward: bool,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<OrientedSupportKey> {
    let (profiles, swapped) =
        crate::surface_recovery::boundary_strip_profiles(surface_id, entities, index)?;
    let [(a, u_forward), (b, other_forward)] = translated_profiles(&profiles)?;
    if u_forward != other_forward || a == b {
        return None;
    }
    let mut profiles = [a, b];
    let v_forward = profiles[0] < profiles[1];
    if !v_forward {
        profiles.reverse();
    }
    let outward_positive = forward ^ swapped ^ !u_forward ^ !v_forward;
    full_trim(
        &profiles,
        outward_positive,
        &brep::face_loops(face_id, entities, index)?,
    )?;
    Some(OrientedSupportKey::SplineStrip {
        strip: Box::new(SplineStripKey { profiles }),
        outward_positive,
    })
}

pub(super) fn line_on_strip(line: brep::LineSupport, strip: &SplineStripKey) -> bool {
    let Some([a, b]) = strip.profiles[0].ends() else {
        return false;
    };
    let Some([c, d]) = strip.profiles[1].ends() else {
        return false;
    };
    [[a, c], [b, d]].into_iter().any(|ends| {
        ends.into_iter()
            .all(|p| decode(p).is_some_and(|p| bounded::on_line(p, line)))
    })
}

pub(super) fn certified_edge(
    edge: &EdgeGeometry,
    spline: &BSplineSupport,
    loci: &[SupportLocusKey; 2],
) -> Option<SplineCurveKey> {
    if !loci
        .iter()
        .any(|locus| matches!(locus, SupportLocusKey::SplineStrip { .. }))
    {
        return None;
    }
    let key = curve_matches_edge(spline, edge)?;
    for locus in loci {
        match locus {
            SupportLocusKey::SplineStrip { strip } if strip.profiles.contains(&key) => {}
            SupportLocusKey::Plane { normal, offset_mm } => {
                let normal = dequantize_vector(*normal, DIRECTION_Q)?;
                let offset = dequantize_scalar(*offset_mm)?;
                // The whole positive-weight control hull lies in the plane.
                if spline
                    .control_points_mm
                    .iter()
                    .any(|&p| (dot(p, normal) - offset).abs() > TOL)
                {
                    return None;
                }
            }
            _ => return None,
        }
    }
    Some(key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Context, Result};

    #[test]
    fn original_relief_box_has_certified_spline_strips() -> Result<()> {
        let input = include_bytes!("../../validation/fixtures/native_relief_box_0b5fc9e7.step");
        let exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let diagnostic =
            super::super::analytic_boundary_diagnostic(1202, &exchange.data[0].entities);
        let signature = diagnostic
            .signature
            .with_context(|| format!("source refused: {:?}", diagnostic.reason))?;
        assert_eq!(
            signature
                .supports
                .iter()
                .filter(|key| matches!(key, OrientedSupportKey::SplineStrip { .. }))
                .count(),
            98
        );
        assert_eq!(
            signature
                .chains
                .iter()
                .filter(|key| matches!(key, super::super::BoundaryChainKey::Spline { .. }))
                .count(),
            196
        );
        Ok(())
    }

    fn profile(z: f64) -> BSplineSupport {
        BSplineSupport {
            degree: 2,
            control_points_mm: vec![[0., 0., z], [0.25, 0.5, z], [0.75, 0.5, z], [1., 0., z]],
            knots: vec![0., 0., 0., 0.5, 1., 1., 1.],
            weights: None,
        }
    }

    #[test]
    fn full_curve_key_retains_all_spans_and_is_reversal_invariant() -> Result<()> {
        let original = profile(0.);
        let (key, forward) = curve_key(&original).context("original")?;
        assert!(forward);
        assert_eq!(key.spans.len(), 2);
        let mut reversed = original.clone();
        reversed.control_points_mm.reverse();
        reversed.knots.reverse();
        for k in &mut reversed.knots {
            *k = 16. - 4. * *k;
        }
        let (other, forward) = curve_key(&reversed).context("reversed and recharted")?;
        assert!(!forward);
        assert_eq!(key, other);
        let mut changed = original;
        changed.control_points_mm[1][1] += 0.001;
        let (changed, _) = curve_key(&changed).context("changed curvature")?;
        assert_ne!(changed, key);
        assert_eq!(changed.ends(), key.ends());
        Ok(())
    }

    #[test]
    fn strip_requires_translated_positive_weight_monotone_profiles() {
        let mut profiles = [profile(0.), profile(0.01)];
        assert!(translated_profiles(&profiles).is_some());
        profiles[0].weights = Some(vec![1., 2., 2., 1.]);
        profiles[1].weights = Some(vec![2., 4., 4., 2.]);
        assert!(translated_profiles(&profiles).is_some());
        profiles[1].weights = Some(vec![2., 4., 3., 2.]);
        assert!(translated_profiles(&profiles).is_none());
        profiles[1].weights = Some(vec![2., 4., 4., 2.]);
        profiles[1].control_points_mm[1][2] += 0.001;
        assert!(translated_profiles(&profiles).is_none());
        for p in &mut profiles {
            p.control_points_mm[1] = [2., 0., p.control_points_mm[0][2]];
            p.control_points_mm[2] = [-1., 0., p.control_points_mm[0][2]];
        }
        assert!(translated_profiles(&profiles).is_none());
    }

    #[test]
    fn canonical_monotonicity_accepts_stationary_roundoff_but_not_real_backtracking() {
        let profile = |z| BSplineSupport {
            degree: 3,
            control_points_mm: vec![
                [0., 0., z],
                [0.3, -1e-15, z],
                [0.6, 0.25, z],
                [0.6, 0.5, z],
                [0.6, 0.75, z],
                [0.3, 1. + 1e-15, z],
                [0., 1., z],
            ],
            knots: vec![0., 0., 0., 0., 0.5, 0.5, 0.5, 1., 1., 1., 1.],
            weights: None,
        };
        let mut profiles = [profile(0.), profile(0.01)];
        assert!(translated_profiles(&profiles).is_some());
        for p in &mut profiles {
            p.control_points_mm[1][1] = -0.0001;
            p.control_points_mm[5][1] = 1.0001;
        }
        assert!(translated_profiles(&profiles).is_none());
    }

    fn edge(
        start: [f64; 3],
        end: [f64; 3],
        support: CurveSupport,
        forward: bool,
        id: u64,
    ) -> brep::OrientedEdgeUse {
        brep::OrientedEdgeUse {
            oriented_edge_id: id,
            edge_id: id,
            curve_id: id,
            curve_same_sense: forward,
            parameter_forward: forward,
            start_vertex: id,
            end_vertex: id + 1,
            start_mm: start,
            end_mm: end,
            support,
        }
    }

    #[test]
    fn full_trim_rejects_wrong_curve_interior_wrong_sense_and_extra_loop() -> Result<()> {
        let low = profile(0.);
        let high = profile(0.01);
        let keys = [
            curve_key(&low).context("low")?.0,
            curve_key(&high).context("high")?.0,
        ];
        let a = low.control_points_mm[0];
        let b = *low.control_points_mm.last().context("end")?;
        let c = high.control_points_mm[0];
        let d = *high.control_points_mm.last().context("end")?;
        let line = |origin_mm| {
            CurveSupport::Line(brep::LineSupport {
                origin_mm,
                direction: [0., 0., 1.],
            })
        };
        let mut boundary = brep::FaceLoop {
            bound_id: 10,
            loop_id: 11,
            outer: true,
            orientation: true,
            edges: vec![
                edge(a, b, CurveSupport::BSpline(low.clone()), true, 1),
                edge(b, d, line(b), true, 2),
                edge(d, c, CurveSupport::BSpline(high), false, 3),
                edge(c, a, line(a), false, 4),
            ],
        };
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary)).is_some());
        assert!(full_trim(&keys, false, std::slice::from_ref(&boundary)).is_none());
        assert!(full_trim(&keys, true, &[boundary.clone(), boundary.clone()]).is_none());
        boundary.edges[0].parameter_forward = false;
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary)).is_none());
        boundary.edges[0].parameter_forward = true;
        let mut changed = low;
        changed.control_points_mm[1][1] += 0.01;
        boundary.edges[0].support = CurveSupport::BSpline(changed);
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary)).is_none());
        Ok(())
    }
}
