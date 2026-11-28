//! Full-domain spline-strip certificates. Neither fitting sampled points nor
//! matching endpoint pairs establishes a curved edge. Retain complete positive-
//! weight rational Bernstein control nets, ordered spans, and actual trims.
use super::{
    BSplineSupport, BoundaryProofContext, CurveSupport,
    DIRECTION_Q, EdgeGeometry, OrientedSupportKey, POINT_Q_MM, PointKey, SupportLocusKey, bounded,
    dequantize_scalar, dequantize_vector, point_key,
};
use crate::{
    brep,
    math3::{distance, dot, mul, normalize, sub},
};
use monstertruck_geometry::{
    analytic_line_segment::{SegmentLocus, certify_nurbs_segment_locus},
    prelude::{Point3, Vector4},
    proof_bezier_spans,
};
use num_traits::ToPrimitive;
use serde::Serialize;
use std::collections::HashMap;
use step_p21::ast::EntityInstance;

#[cfg(feature = "cad-kernel-monstertruck")]
mod cap_trim;
#[cfg(any(test, feature = "cad-kernel-monstertruck"))]
pub(super) mod coincidence;
mod embedding;
mod restriction;
pub(super) use restriction::{ProfilePosition, complete_profile_cover};
use restriction::{complete_directed_cover, profile_piece};

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
    pub(super) fn map_points(&mut self, mut map: impl FnMut(&mut PointKey)) {
        for point in self.spans.iter_mut().flat_map(|span| &mut span.points_mm) {
            map(point);
        }
    }

    pub(super) fn ends(&self) -> Option<[PointKey; 2]> {
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

impl SplineStripKey {
    pub(super) fn map_points(&mut self, mut map: impl FnMut(&mut PointKey)) {
        for profile in &mut self.profiles {
            profile.map_points(&mut map);
        }
    }
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
    certified_bspline_extraction_error_mm(spline).map(|_| ())
}

/// The actual (input-dependent) conservative arithmetic margin for the
/// polynomial/rational B-spline to homogeneous Bernstein conversion.
/// The strict key proof only needs to know whether this fits a quarter of
/// its coordinate quantum; a geometric proximity proof also needs the value.
fn certified_bspline_extraction_error_mm(spline: &BSplineSupport) -> Option<f64> {
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
    (error_mm.is_finite() && error_mm <= POINT_Q_MM / 4.).then_some(error_mm)
}

/// A finite union of planes and elliptic cylinders cannot contain an
/// open nonplanar polynomial surface patch. Substituting its polynomial
/// parameterization into a carrier's implicit equation must give the
/// zero polynomial on an open set. For an elliptic cylinder, the sum of
/// squares of two transverse coordinate polynomials would be constant,
/// forcing both to be constant and the patch to have rank at most one.
/// For a plane, every Bernstein control point must be coplanar.
fn nonplanar_polynomial_profiles(profiles: &[BSplineSupport; 2], tolerance_mm: f64) -> bool {
    if profiles.iter().any(|profile| {
        profile.weights.as_ref().is_some_and(|weights| {
            weights.first().is_none_or(|first| {
                !first.is_finite()
                    || *first <= 0.
                    || weights.iter().any(|w| w.to_bits() != first.to_bits())
            })
        })
    }) {
        return false;
    }
    // Constant weights on each row also permit a common rational linear
    // reparameterization in the sweep direction; the geometric patch
    // remains the polynomial homotopy of the two profile rows.
    let Some(a) = curve_evidence(&profiles[0]) else {
        return false;
    };
    let Some(b) = curve_evidence(&profiles[1]) else {
        return false;
    };
    if a.forward != b.forward || a.controls.len() != b.controls.len() {
        return false;
    }
    a.controls.iter().zip(&b.controls).any(|(first, second)| {
        let points = first
            .iter()
            .chain(second)
            .map(|p| [p.x / p.w, p.y / p.w, p.z / p.w])
            .collect::<Vec<_>>();
        nonplanar_control_net(&points, tolerance_mm)
    })
}

fn nonplanar_control_net(points: &[[f64; 3]], tolerance_mm: f64) -> bool {
    let Some(&origin) = points.first() else {
        return false;
    };
    let vectors = points.iter().map(|&p| sub(p, origin)).collect::<Vec<_>>();
    let norm = crate::math3::norm;
    let cross = crate::math3::cross;
    let Some(&a) = vectors
        .iter()
        .max_by(|&&a, &&b| norm(a).total_cmp(&norm(b)))
    else {
        return false;
    };
    let Some(&b) = vectors
        .iter()
        .max_by(|&&b, &&c| norm(cross(a, b)).total_cmp(&norm(cross(a, c))))
    else {
        return false;
    };
    let normal = cross(a, b);
    // Reserve the existing source proof length for each control point.
    // Differences can have twice that error. Bound all terms of the
    // perturbed scalar triple product, including arithmetic roundoff.
    let e = 2. * tolerance_mm;
    let (na, nb) = (norm(a), norm(b));
    vectors.iter().any(|&c| {
        let nc = norm(c);
        let pairwise = nb.mul_add(nc, na.mul_add(nc, na * nb));
        let bound = (e * e).mul_add(na + nb + nc, e.mul_add(pairwise, e.powi(3)));
        let bound = (32. * f64::EPSILON * na * nb).mul_add(nc, bound);
        let determinant = dot(normal, c).abs();
        determinant.is_finite() && determinant > bound
    })
}

/// Require a complete certified strip, not just a nonplanar untrimmed
/// carrier which could have been trimmed to one of its planar spans.
pub(super) fn has_nonplanar_polynomial_patch(
    face_id: u64,
    surface_id: u64,
    forward: bool,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>, tolerance_mm: f64,
) -> bool {
    crate::surface_recovery::boundary_strip_profiles(surface_id, entities, index).is_some_and(
        |(profiles, _)| {
            nonplanar_polynomial_profiles(&profiles, tolerance_mm)
                && oriented_key(
                    face_id,
                    surface_id,
                    forward,
                    entities,
                    index,
                    &mut BoundaryProofContext::new(tolerance_mm),
                )
                .is_some()
        },
    )
}

/// Canonicalize the whole curve, allowing affine reversal and different
/// knot interval lengths. Each positive span retains every homogeneous
/// Bernstein coefficient, with its own common weight scale removed.
#[derive(Clone)]
struct CurveEvidence {
    key: SplineCurveKey,
    forward: bool,
    controls: Vec<Vec<Vector4>>,
}

fn curve_evidence(spline: &BSplineSupport) -> Option<CurveEvidence> {
    well_conditioned(spline)?;
    let rational = spline.rational_curve()?;
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
    let controls = spans
        .iter()
        .map(|span| {
            span.controls
                .iter()
                .map(|p| Vector4::new(p[0], p[1], p[2], p[3]))
                .collect()
        })
        .collect();
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
    Some(CurveEvidence {
        key,
        forward,
        controls,
    })
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
        minimum >= 0. && maximum > POINT_Q_MM
    })
}

type DirectedProfile = (SplineCurveKey, bool);

fn translated_geometry(
    profiles: &[BSplineSupport; 2],
    cells: &mut BoundaryProofContext,
) -> Option<([DirectedProfile; 2], [f64; 3], [f64; 3])> {
    let tolerance_mm = cells.tolerance_mm();
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
    let axis = normalize(delta, POINT_Q_MM)?;
    if first
        .control_points_mm
        .iter()
        .zip(&second.control_points_mm)
        .any(|(&a, &b)| distance(sub(b, a), delta) > tolerance_mm)
    {
        return None;
    }
    let chord = sub(
        *first.control_points_mm.last()?,
        *first.control_points_mm.first()?,
    );
    let evidence = [curve_evidence(first)?, curve_evidence(second)?];
    if cells.enabled() {
        for profile in &evidence {
            for (span, controls) in profile.key.spans.iter().zip(&profile.controls) {
                for (&key, p) in span.points_mm.iter().zip(controls) {
                    cells.record_point(key, [p.x / p.w, p.y / p.w, p.z / p.w]);
                }
            }
        }
    }
    let keys = evidence.map(|evidence| (evidence.key, evidence.forward));
    Some((keys, axis, chord))
}

fn translated_profiles(
    profiles: &[BSplineSupport; 2],
    cells: &mut BoundaryProofContext,
) -> Option<[DirectedProfile; 2]> {
    let (keys, axis, chord) = translated_geometry(profiles, cells)?;
    [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], chord]
        .into_iter()
        .filter_map(|e| normalize(sub(e, mul(axis, dot(axis, e))), 0.0))
        .flat_map(|d| [d, mul(d, -1.)])
        .any(|d| keys.iter().all(|(key, _)| monotone_profile(key, d)))
        .then_some(keys)
}

/// The global monotonicity check is a cheap certificate for most profiles.
/// For turning profiles prove each span injective and every projected
/// span pair separated, rather than weakening the monotonicity inequalities.
fn embedded_translated_profiles(
    profiles: &[BSplineSupport; 2],
    cells: &mut BoundaryProofContext,
) -> Option<[DirectedProfile; 2]> {
    translated_profiles(profiles, cells).or_else(|| {
        let (keys, axis, _) = translated_geometry(profiles, cells)?;
        keys.iter()
            .all(|(key, _)| embedding::is_embedded(key, axis))
            .then_some(keys)
    })
}

/// The topological trim vertices are shared by the incident faces and
/// adjacent pieces. The spline's endpoint evaluations may straddle a
/// coordinate-key cell boundary despite being within the proof budget.
fn curve_matches_edge(
    spline: &BSplineSupport,
    edge: &EdgeGeometry, tolerance_mm: f64,
) -> Option<(CurveEvidence, [PointKey; 2])> {
    let evidence = curve_evidence(spline)?;
    let actual = if edge.parameter_forward == evidence.forward {
        [edge.start, edge.end]
    } else {
        [edge.end, edge.start]
    };
    let mut expected = [
        *spline.control_points_mm.first()?,
        *spline.control_points_mm.last()?,
    ];
    if !evidence.forward {
        expected.reverse();
    }
    let expected_keys = evidence.key.ends()?;
    let bound = |point: PointKey, end: [f64; 3]| {
        decode(point)
            .is_some_and(|point| distance(point, end) + POINT_Q_MM * 3.0_f64.sqrt() / 2. <= tolerance_mm)
    };
    if actual[0] == actual[1] || !bound(actual[0], expected[0]) || !bound(actual[1], expected[1]) {
        return None;
    }
    // Preserve the authored curve sense. If quantized endpoints do not
    // identify the direction exactly and both matches are plausible, refuse.
    if actual != expected_keys && bound(actual[0], expected[1]) && bound(actual[1], expected[0]) {
        return None;
    }
    Some((evidence, actual))
}

/// Coverage of an independently certified restriction of a canonical cap.
#[derive(Debug, Clone)]
pub(super) struct ProfilePiece {
    pub profile: Box<SplineCurveKey>,
    pub interval: [ProfilePosition; 2],
}

fn generator_piece(parent: [PointKey; 2], edge: &EdgeGeometry, tolerance_mm: f64) -> Option<()> {
    let [a, b] = parent.map(decode);
    let (a, b) = (a?, b?);
    let delta = sub(b, a);
    let direction = normalize(delta, POINT_Q_MM)?;
    let [start, end] = [decode(edge.start)?, decode(edge.end)?];
    let line = brep::LineSupport {
        origin_mm: a,
        direction,
    };
    if !bounded::on_line(start, line, tolerance_mm) || !bounded::on_line(end, line, tolerance_mm) {
        return None;
    }
    let lo = dot(sub(start, a), direction);
    let hi = dot(sub(end, a), direction);
    if lo < -tolerance_mm || hi > crate::math3::norm(delta) + tolerance_mm || lo >= hi {
        return None;
    }
    match &edge.support {
        CurveSupport::Line(line) => {
            let directed = dot(sub(end, start), line.direction);
            (bounded::on_line(start, *line, tolerance_mm)
                && bounded::on_line(end, *line, tolerance_mm)
                && (directed > 0.) == edge.parameter_forward)
                .then_some(())
        }
        CurveSupport::BSpline(spline) => {
            let rational = spline.rational_curve()?;
            if certify_nurbs_segment_locus(
                &rational,
                Point3::new(start[0], start[1], start[2]),
                Point3::new(end[0], end[1], end[2]),
                tolerance_mm,
            ) != SegmentLocus::WithinTolerance
            {
                return None;
            }
            let (evidence, ends) = curve_matches_edge(spline, edge, tolerance_mm)?;
            let direction = if ends[0] == edge.start {
                direction
            } else {
                mul(direction, -1.)
            };
            monotone_profile(&evidence.key, direction).then_some(())
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
    loops: &[brep::FaceLoop], tolerance_mm: f64,
) -> Option<()> {
    let [boundary] = loops else {
        return None;
    };
    let [first_start, first_end] = profiles[0].ends()?;
    let [second_start, second_end] = profiles[1].ends()?;
    let mut generators = [[first_end, second_end], [second_start, first_start]];
    let cap_forward = [
        forward == boundary.orientation,
        forward != boundary.orientation,
    ];
    if forward != boundary.orientation {
        for ends in &mut generators {
            ends.reverse();
        }
    }
    let mut generator_covers = [Vec::new(), Vec::new()];
    let mut covers = [Vec::new(), Vec::new()];
    for edge in &boundary.edges {
        let ends = [point_key(edge.start_mm)?, point_key(edge.end_mm)?];
        let geometry = EdgeGeometry {
            start: ends[0],
            end: ends[1],
            start_vertex: edge.start_vertex,
            end_vertex: edge.end_vertex,
            parameter_forward: edge.parameter_forward,
            curve_id: edge.curve_id,
            support: edge.support.clone(),
        };
        let matching_generators = generators
            .iter()
            .enumerate()
            .filter_map(|(i, &generator)| generator_piece(generator, &geometry, tolerance_mm).map(|()| i))
            .collect::<Vec<_>>();
        match matching_generators.as_slice() {
            [i] => {
                generator_covers[*i].push(ends);
                continue;
            }
            [] => {}
            _ => return None,
        }
        let CurveSupport::BSpline(spline) = &geometry.support else {
            return None;
        };
        let (key, actual) = curve_matches_edge(spline, &geometry, tolerance_mm)?;
        let matches = profiles
            .iter()
            .enumerate()
            .filter_map(|(i, parent)| {
                let interval = profile_piece(parent, &key, actual, tolerance_mm)?;
                let mut directed = interval.map(|position| position.point);
                if !cap_forward[i] {
                    directed.reverse();
                }
                (directed == ends).then_some((i, interval))
            })
            .collect::<Vec<_>>();
        let [(i, interval)] = matches.as_slice() else {
            return None;
        };
        covers[*i].push(*interval);
    }
    for (generator, cover) in generators.into_iter().zip(&generator_covers) {
        complete_directed_cover(generator, cover)?;
    }
    for (parent, cover) in profiles.iter().zip(&covers) {
        complete_profile_cover(parent, cover)?;
    }
    Some(())
}

fn canonical_strip(
    surface_id: u64,
    forward: bool,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    cells: &mut BoundaryProofContext,
) -> Option<([SplineCurveKey; 2], bool)> {
    let (profiles, swapped) =
        crate::surface_recovery::boundary_strip_profiles(surface_id, entities, index)?;
    let [(a, u_forward), (b, other_forward)] = embedded_translated_profiles(&profiles, cells)?;
    if u_forward != other_forward || a == b {
        return None;
    }
    let mut profiles = [a, b];
    let v_forward = profiles[0] < profiles[1];
    if !v_forward {
        profiles.reverse();
    }
    Some((profiles, forward ^ swapped ^ !u_forward ^ !v_forward))
}

pub(super) fn oriented_key(
    face_id: u64,
    surface_id: u64,
    forward: bool,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    cells: &mut BoundaryProofContext,
) -> Option<OrientedSupportKey> {
    let tolerance_mm = cells.tolerance_mm();
    let (profiles, outward_positive) =
        canonical_strip(surface_id, forward, entities, index, cells)?;
    full_trim(
        &profiles,
        outward_positive,
        &brep::face_loops(face_id, entities, index)?, tolerance_mm,
    )?;
    Some(OrientedSupportKey::SplineStrip {
        strip: Box::new(SplineStripKey { profiles }),
        outward_positive,
    })
}

/// Check that an already-authored complete trim traverses an embedded strip
/// in exactly the opposite sense to its declared `FACE_SURFACE.same_sense`.
/// This does not by itself establish an oriented closed shell; callers must
/// separately validate the entire reconciled boundary, including coedges.
pub(super) fn reversed_face_sense_is_proven(
    face_id: u64,
    surface_id: u64,
    forward: bool,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>, tolerance_mm: f64,
) -> bool {
    let Some((profiles, outward_positive)) = canonical_strip(
        surface_id,
        forward,
        entities,
        index,
        &mut BoundaryProofContext::new(tolerance_mm),
    ) else {
        return false;
    };
    let Some(loops) = brep::face_loops(face_id, entities, index) else {
        return false;
    };
    full_trim(&profiles, outward_positive, &loops, tolerance_mm).is_none()
        && full_trim(&profiles, !outward_positive, &loops, tolerance_mm).is_some()
}

pub(super) fn line_on_strip(line: brep::LineSupport, strip: &SplineStripKey, tolerance_mm: f64) -> bool {
    let Some([a, b]) = strip.profiles[0].ends() else {
        return false;
    };
    let Some([c, d]) = strip.profiles[1].ends() else {
        return false;
    };
    [[a, c], [b, d]].into_iter().any(|ends| {
        ends.into_iter()
            .all(|p| decode(p).is_some_and(|p| bounded::on_line(p, line, tolerance_mm)))
    })
}

pub(super) fn certified_edge(
    edge: &EdgeGeometry,
    spline: &BSplineSupport,
    loci: &[SupportLocusKey; 2], tolerance_mm: f64,
) -> Option<ProfilePiece> {
    if !loci
        .iter()
        .any(|locus| matches!(locus, SupportLocusKey::SplineStrip { .. }))
    {
        return None;
    }
    let (key, actual) = curve_matches_edge(spline, edge, tolerance_mm)?;
    let mut selected: Option<ProfilePiece> = None;
    for locus in loci {
        match locus {
            SupportLocusKey::SplineStrip { strip } => {
                let matches = strip
                    .profiles
                    .iter()
                    .filter_map(|parent| {
                        Some(ProfilePiece {
                            profile: Box::new(parent.clone()),
                            interval: profile_piece(parent, &key, actual, tolerance_mm)?,
                        })
                    })
                    .collect::<Vec<_>>();
                let [piece] = matches.as_slice() else {
                    return None;
                };
                if selected.as_ref().is_some_and(|other| {
                    other.profile != piece.profile || other.interval != piece.interval
                }) {
                    return None;
                }
                selected = Some(piece.clone());
            }
            SupportLocusKey::Plane { normal, offset_mm } => {
                let normal = dequantize_vector(*normal, DIRECTION_Q)?;
                let offset = dequantize_scalar(*offset_mm)?;
                if spline
                    .control_points_mm
                    .iter()
                    .any(|&p| (dot(p, normal) - offset).abs() > tolerance_mm)
                {
                    return None;
                }
            }
            _ => return None,
        }
    }
    selected
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Context, Result};

    #[test]
    fn original_relief_box_has_certified_spline_strips() -> Result<()> {
        let input = include_bytes!("../../validation/fixtures/native_relief_box_0b5fc9e7.step");
        let exchange = monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange;
        let diagnostic = super::super::analytic_boundary_diagnostic(
            1202,
            &exchange.data[0].entities,
            &crate::step_graph::build_index(&exchange.data[0].entities),
        );
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
    fn polynomial_exclusion_requires_a_single_nonplanar_full_span() {
        let mut profiles = [profile(0.), profile(0.01)];
        assert!(nonplanar_polynomial_profiles(&profiles, crate::boundary_proof::TEST_BOUNDARY_TOL_MM));
        // Different common weights in the two sweep rows only rechart v.
        profiles[0].weights = Some(vec![1.; 4]);
        profiles[1].weights = Some(vec![2.; 4]);
        assert!(nonplanar_polynomial_profiles(&profiles, crate::boundary_proof::TEST_BOUNDARY_TOL_MM));
        // A rational arc is not ruled out by polynomial algebra.
        profiles[0].weights = Some(vec![1., 0.8, 0.8, 1.]);
        profiles[1].weights = Some(vec![2., 1.6, 1.6, 2.]);
        assert!(!nonplanar_polynomial_profiles(&profiles, crate::boundary_proof::TEST_BOUNDARY_TOL_MM));

        // A bent piecewise-linear profile has a nonplanar complete
        // control net but each individual swept span is a plane.
        let bent = |z| BSplineSupport {
            degree: 1,
            control_points_mm: vec![[0., 0., z], [1., 1., z], [2., 0., z]],
            knots: vec![0., 0., 0.5, 1., 1.],
            weights: None,
        };
        assert!(!nonplanar_polynomial_profiles(&[bent(0.), bent(1.)], crate::boundary_proof::TEST_BOUNDARY_TOL_MM));

        // A degree-elevated line remains planar.
        let mut profiles = [profile(0.), profile(0.01)];
        for p in &mut profiles {
            for point in &mut p.control_points_mm {
                point[1] = 0.;
            }
        }
        assert!(!nonplanar_polynomial_profiles(&profiles, crate::boundary_proof::TEST_BOUNDARY_TOL_MM));
    }

    #[test]
    fn full_curve_key_retains_all_spans_and_is_reversal_invariant() -> Result<()> {
        let original = profile(0.);
        let CurveEvidence { key, forward, .. } = curve_evidence(&original).context("original")?;
        assert!(forward);
        assert_eq!(key.spans.len(), 2);
        let mut reversed = original.clone();
        reversed.control_points_mm.reverse();
        reversed.knots.reverse();
        for k in &mut reversed.knots {
            *k = 16. - 4. * *k;
        }
        let CurveEvidence {
            key: other,
            forward,
            ..
        } = curve_evidence(&reversed).context("reversed and recharted")?;
        assert!(!forward);
        assert_eq!(key, other);
        let mut changed = original;
        changed.control_points_mm[1][1] += 0.001;
        let changed = curve_evidence(&changed).context("changed curvature")?.key;
        assert_ne!(changed, key);
        assert_eq!(changed.ends(), key.ends());
        Ok(())
    }

    #[test]
    fn strip_requires_translated_positive_weight_monotone_profiles() {
        let mut profiles = [profile(0.), profile(0.01)];
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_some());
        profiles[0].weights = Some(vec![1., 2., 2., 1.]);
        profiles[1].weights = Some(vec![2., 4., 4., 2.]);
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_some());
        profiles[1].weights = Some(vec![2., 4., 3., 2.]);
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_none());
        profiles[1].weights = Some(vec![2., 4., 4., 2.]);
        profiles[1].control_points_mm[1][2] += 0.001;
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_none());
        for p in &mut profiles {
            p.control_points_mm[1] = [2., 0., p.control_points_mm[0][2]];
            p.control_points_mm[2] = [-1., 0., p.control_points_mm[0][2]];
        }
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_none());
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
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_some());
        for p in &mut profiles {
            p.control_points_mm[1][1] = -0.0001;
            p.control_points_mm[5][1] = 1.0001;
        }
        assert!(translated_profiles(&profiles, &mut BoundaryProofContext::default()).is_none());
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
            curve_evidence(&low).context("low")?.key,
            curve_evidence(&high).context("high")?.key,
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
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary), crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_some());
        let unpartitioned = boundary.clone();
        let midpoint = [b[0], b[1], f64::midpoint(b[2], d[2])];
        boundary.edges.splice(
            1..2,
            [
                edge(b, midpoint, line(b), true, 2),
                edge(midpoint, d, line(b), true, 5),
            ],
        );
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary), crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_some());
        boundary.edges[2].start_mm = b;
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary), crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_none());
        boundary = unpartitioned;
        assert!(full_trim(&keys, false, std::slice::from_ref(&boundary), crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_none());
        assert!(full_trim(&keys, true, &[boundary.clone(), boundary.clone()], crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_none());
        boundary.edges[0].parameter_forward = false;
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary), crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_none());
        boundary.edges[0].parameter_forward = true;
        let mut changed = low;
        changed.control_points_mm[1][1] += 0.01;
        boundary.edges[0].support = CurveSupport::BSpline(changed);
        assert!(full_trim(&keys, true, std::slice::from_ref(&boundary), crate::boundary_proof::TEST_BOUNDARY_TOL_MM).is_none());
        Ok(())
    }
}
