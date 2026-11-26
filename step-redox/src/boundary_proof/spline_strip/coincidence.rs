//! Precision-qualified geometric containment of two *complete* rational curves.
//!
//! Equal STEP parameter values are not equal physical points if authored knot
//! charts differ. Instead, a common endpoint-spanned projected coordinate
//! pairs geometric positions by continuity without sampling a parameter grid. Bernstein
//! convex hulls then enclose the entire residual at that coordinate.
//! Subdivision tightens the bound and has a finite fail-closed work budget.
//! This does not establish exact source equivalence or repair any STEP entity.

use super::{BSplineSupport, certified_bspline_extraction_error_mm};
use crate::math3::{add, distance, dot, mul, norm, normalize, sub};
use monstertruck_geometry::{proof_bezier_spans, subdivide_homogeneous_bezier};

const MAX_DEPTH: usize = 32;
const MAX_COMPARISONS: usize = 262_144;

#[derive(Clone)]
struct Span {
    controls: Vec<[f64; 4]>,
    lo: f64,
    hi: f64,
    floor: f64,
    extraction_error_mm: f64,
    depth: usize,
}

fn cartesian(control: &[f64; 4]) -> Option<[f64; 3]> {
    let [x, y, z, w] = *control;
    let result = [x / w, y / w, z / w];
    (w > 0.0 && result.iter().all(|v| v.is_finite())).then_some(result)
}

fn roundoff_floor(controls: &[[f64; 4]], depth: usize, extraction_error_mm: f64) -> Option<f64> {
    let mut smallest = f64::INFINITY;
    let mut largest = 0.0_f64;
    let mut scale = 1.0_f64;
    for control in controls {
        let p = cartesian(control)?;
        smallest = smallest.min(control[3]);
        largest = largest.max(control[3]);
        scale = scale.max(norm(p));
    }
    let n = f64::from(u32::try_from(controls.len()).ok()?);
    let depth = f64::from(u32::try_from(depth + 1).ok()?);
    let ratio = largest / smallest;
    // Retain the input-specific, fully certified NURBS extraction bound
    // through every subdivision rather than reserving an entire quantum.
    let arithmetic = 2048.0 * f64::EPSILON * scale * ratio.powi(2) * n.powi(3);
    let floor = arithmetic.mul_add(depth, extraction_error_mm);
    floor.is_finite().then_some(floor)
}

impl Span {
    fn new(
        controls: Vec<[f64; 4]>,
        axis: [f64; 3],
        depth: usize,
        extraction_error_mm: f64,
    ) -> Option<Self> {
        if controls.len() < 2 || depth > MAX_DEPTH {
            return None;
        }
        let floor = roundoff_floor(&controls, depth, extraction_error_mm)?;
        // A rational curve lies in the convex hull of its positive-weight
        // Cartesian controls, regardless of knot parameterization, control
        // projection order or local backtracking. Every possible matching
        // projected position is included in this *entire* hull interval.
        let (lo, hi) =
            controls
                .iter()
                .try_fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                    let z = dot(cartesian(p)?, axis);
                    z.is_finite().then_some((lo.min(z), hi.max(z)))
                })?;
        if !lo.is_finite() || !hi.is_finite() {
            return None;
        }
        Some(Self {
            controls,
            lo,
            hi,
            floor,
            extraction_error_mm,
            depth,
        })
    }

    fn split(&self, axis: [f64; 3]) -> Option<[Self; 2]> {
        let (a, b) = subdivide_homogeneous_bezier(&self.controls)?;
        Some([
            Self::new(a, axis, self.depth + 1, self.extraction_error_mm)?,
            Self::new(b, axis, self.depth + 1, self.extraction_error_mm)?,
        ])
    }

    fn end(&self, first: bool) -> Option<[f64; 3]> {
        cartesian(if first {
            self.controls.first()?
        } else {
            self.controls.last()?
        })
    }
}

fn spans(spline: &BSplineSupport, axis: [f64; 3]) -> Option<Vec<Span>> {
    let extraction_error_mm = certified_bspline_extraction_error_mm(spline)?;
    let rational = spline.rational_curve()?;
    let pieces = proof_bezier_spans(rational.non_rationalized())?;
    let spans = pieces
        .into_iter()
        .map(|piece| Span::new(piece.controls, axis, 0, extraction_error_mm))
        .collect::<Option<Vec<_>>>()?;
    for pair in spans.windows(2) {
        if distance(pair[0].end(false)?, pair[1].end(true)?) > pair[0].floor + pair[1].floor {
            return None;
        }
    }
    Some(spans)
}

fn reverse(mut spans: Vec<Span>, axis: [f64; 3]) -> Option<Vec<Span>> {
    spans.reverse();
    spans
        .into_iter()
        .map(|span| {
            let mut controls = span.controls;
            controls.reverse();
            Span::new(controls, axis, span.depth, span.extraction_error_mm)
        })
        .collect()
}

fn extension(
    endpoint: [f64; 3],
    step: f64,
    axis: [f64; 3],
    extraction_error_mm: f64,
) -> Option<Span> {
    let other = add(endpoint, mul(axis, step));
    let controls = [endpoint, other].map(|p| [p[0], p[1], p[2], 1.0]);
    Span::new(controls.to_vec(), axis, 0, extraction_error_mm)
}

/// All rational controls have positive weights. Under any affine map, their
/// dehomogenized curve stays in the Euclidean convex hull of the mapped
/// controls. Subtract a common affine function of the projected coordinate;
/// the map differences therefore bound every pair at equal projection.
fn pair_bound(a: &Span, b: &Span, axis: [f64; 3]) -> Option<f64> {
    let (anchor, slope) = [a, b]
        .into_iter()
        .find_map(|span| {
            let anchor = span.end(true)?;
            let delta = sub(span.end(false)?, anchor);
            let dz = dot(delta, axis);
            (dz.abs() > span.floor).then_some((anchor, mul(delta, 1.0 / dz)))
        })
        .unwrap_or((a.end(true)?, [0.0; 3]));
    let residual = |point| {
        let relative = sub(point, anchor);
        sub(relative, mul(slope, dot(relative, axis)))
    };
    let aa = a
        .controls
        .iter()
        .map(cartesian)
        .map(|p| p.map(residual))
        .collect::<Option<Vec<_>>>()?;
    let bb = b
        .controls
        .iter()
        .map(cartesian)
        .map(|p| p.map(residual))
        .collect::<Option<Vec<_>>>()?;
    let radius = aa
        .iter()
        .flat_map(|&p| bb.iter().map(move |&q| norm(sub(p, q))))
        .fold(0.0_f64, f64::max);
    // The affine residual map has operator norm at most 1+|slope|;
    // propagate both dehomogenization/blossom errors through that map.
    let conditioning = 16.0 * (a.floor + b.floor) * (1.0 + norm(slope));
    (radius.is_finite() && conditioning.is_finite()).then_some(radius + conditioning)
}

fn certify_pair(
    a: &Span,
    b: &Span,
    axis: [f64; 3],
    limit: f64,
    work: &mut usize,
    proven: &mut f64,
) -> Option<()> {
    *work = work.checked_sub(1)?;
    let begin = a.lo.max(b.lo);
    let end = a.hi.min(b.hi);
    // No common projected coordinate exists inside disjoint intervals.
    if end + a.floor + b.floor < begin {
        return Some(());
    }
    let upper = pair_bound(a, b, axis)?;
    if upper <= limit {
        *proven = proven.max(upper);
        return Some(());
    }
    if (a.hi - a.lo) >= (b.hi - b.lo) {
        for child in a.split(axis)? {
            certify_pair(&child, b, axis, limit, work, proven)?;
        }
    } else {
        for child in b.split(axis)? {
            certify_pair(a, &child, axis, limit, work, proven)?;
        }
    }
    Some(())
}

/// Complete authored NURBS with the same control net, knots and rational
/// weights represent the same mathematical locus (independent of evaluation
/// rounding). Check the exact coefficients before expensive interval
/// subdivision, including the only equally exact alternative: reversing
/// every coefficient and reflecting the full authored knot vector.
fn exactly_same_bspline(a: &BSplineSupport, b: &BSplineSupport) -> bool {
    if a == b {
        return true;
    }
    if a.degree != b.degree
        || a.control_points_mm.len() != b.control_points_mm.len()
        || a.knots.len() != b.knots.len()
        || !a
            .control_points_mm
            .iter()
            .eq(b.control_points_mm.iter().rev())
    {
        return false;
    }
    let weights_reversed = match (&a.weights, &b.weights) {
        (None, None) => true,
        (Some(a), Some(b)) => a.iter().eq(b.iter().rev()),
        _ => false,
    };
    if !weights_reversed {
        return false;
    }
    let (Some(&start), Some(&end)) = (a.knots.first(), a.knots.last()) else {
        return false;
    };
    a.knots
        .iter()
        .zip(b.knots.iter().rev())
        .all(|(&x, &y)| x.to_bits() == ((start + end) - y).to_bits())
}

/// A two-sided full-domain Hausdorff *upper bound*, at a user-supplied source
/// precision. Returns `None` for unproved geometry, badly conditioned inputs,
/// unresolved geometry or exhausted proof work. It never edits the source.
///
/// This certificate proves proximity, **not exact identity** or equivalent
/// occupied material. Keep it separate from strict boundary signatures.
#[must_use = "The optional certified coincident spline loci outcome must be inspected to distinguish a value from absence."]
pub fn certified_coincident_spline_loci(
    a: &BSplineSupport,
    b: &BSplineSupport,
    precision_mm: f64,
) -> Option<f64> {
    if !precision_mm.is_finite() || precision_mm <= 0.0 {
        return None;
    }
    let start_a = *a.control_points_mm.first()?;
    let end_a = *a.control_points_mm.last()?;
    let axis = normalize(sub(end_a, start_a), precision_mm)?;
    let mut aa = spans(a, axis)?;
    let mut bb = spans(b, axis)?;
    if exactly_same_bspline(a, b) {
        return Some(0.0);
    }
    let a_start = aa.first()?.end(true)?;
    let a_end = aa.last()?.end(false)?;
    let b_start = bb.first()?.end(true)?;
    let b_end = bb.last()?.end(false)?;
    if distance(a_start, b_end) + distance(a_end, b_start)
        < distance(a_start, b_start) + distance(a_end, b_end)
    {
        // The source basis may enumerate the same entire curve backward.
        bb = reverse(bb, axis)?;
    }
    let za = [dot(a_start, axis), dot(a_end, axis)];
    let zb = [
        dot(bb.first()?.end(true)?, axis),
        dot(bb.last()?.end(false)?, axis),
    ];
    // Extend each curve's endpoint with a virtual axial segment until
    // *both complete control hulls* are in the common projected interval.
    // A point on such a virtual extension differs from its actual endpoint
    // by at most `extension_error`. Thus pairing against the extended curve
    // costs that much in the final two-sided Hausdorff bound.
    let zlo = aa
        .iter()
        .chain(&bb)
        .map(|span| span.lo)
        .fold(za[0].min(zb[0]), f64::min);
    let zhi = aa
        .iter()
        .chain(&bb)
        .map(|span| span.hi)
        .fold(za[1].max(zb[1]), f64::max);
    let margin = [za, zb]
        .into_iter()
        .flat_map(|ends| [ends[0] - zlo, zhi - ends[1]])
        .fold(0.0_f64, f64::max);
    if !margin.is_finite() || margin >= precision_mm {
        return None;
    }
    if za[0] > zlo {
        aa.insert(
            0,
            extension(a_start, zlo - za[0], axis, aa.first()?.extraction_error_mm)?,
        );
    }
    if zb[0] > zlo {
        bb.insert(
            0,
            extension(
                bb.first()?.end(true)?,
                zlo - zb[0],
                axis,
                bb.first()?.extraction_error_mm,
            )?,
        );
    }
    if za[1] < zhi {
        aa.push(extension(
            a_end,
            zhi - za[1],
            axis,
            aa.last()?.extraction_error_mm,
        )?);
    }
    if zb[1] < zhi {
        bb.push(extension(
            bb.last()?.end(false)?,
            zhi - zb[1],
            axis,
            bb.last()?.extraction_error_mm,
        )?);
    }
    let mut work = MAX_COMPARISONS;
    let mut bound = 0.0_f64;
    for a in &aa {
        for b in &bb {
            if b.lo > a.hi + a.floor + b.floor || a.lo > b.hi + a.floor + b.floor {
                continue;
            }
            certify_pair(a, b, axis, precision_mm - margin, &mut work, &mut bound)?;
        }
    }
    Some(bound + margin)
}

/// Independent finite-precision evidence for one *authored* STEP trim. This
/// does not modify original entities or authorize strict source recovery.
#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct SplineCapPrecisionWitness {
    pub face_id: u64,
    pub surface_id: u64,
    pub source_edge_id: u64,
    pub surface_row: usize,
    pub bound_mm: f64,
}

/// Why a particular original source spline face has no two-cap certificate.
/// No failure is converted into a guessed cap or an implicit tolerance increase.
#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SplineCapPrecisionGapReason {
    FaceLoopsUnavailable,
    CapNotProven { surface_row: usize },
    CapAmbiguous { surface_row: usize },
    SameCapEdgeUsedTwice,
    IncompleteFaceTrim,
}

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct SplineCapPrecisionGap {
    pub face_id: u64,
    pub surface_id: u64,
    pub reason: SplineCapPrecisionGapReason,
}

/// Every recorded row is complete and unique within its source face; omitted
/// faces have no independently established dual-cap coincidence certificate.
#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct SplineCapPrecisionAudit {
    pub solid_id: u64,
    pub authored_uncertainty_mm: f64,
    pub target_bound_mm: f64,
    pub strip_faces: usize,
    pub fully_certified_strip_faces: usize,
    pub cap_witnesses: Vec<SplineCapPrecisionWitness>,
    pub unresolved_strip_faces: Vec<SplineCapPrecisionGap>,
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn verify_source_strip_cap_rows(
    face_id: u64,
    surface_id: u64,
    rows: &[BSplineSupport; 2],
    entities: &[step_p21::ast::EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
    coordinate_tolerance: f64,
    mm_per_unit: f64,
) -> Result<[SplineCapPrecisionWitness; 2], SplineCapPrecisionGapReason> {
    use crate::brep;
    use SplineCapPrecisionGapReason as Gap;

    let loops = brep::face_loops(face_id, entities, index).ok_or(Gap::FaceLoopsUnavailable)?;
    let mut witnesses = Vec::with_capacity(2);
    for (surface_row, profile) in rows.iter().enumerate() {
        let mut matches = loops
            .iter()
            .flat_map(|wire| &wire.edges)
            .filter_map(|edge| {
                let brep::CurveSupport::BSpline(curve) = &edge.support else {
                    return None;
                };
                let bound = certified_coincident_spline_loci(profile, curve, coordinate_tolerance)?;
                Some(SplineCapPrecisionWitness {
                    face_id,
                    surface_id,
                    source_edge_id: edge.edge_id,
                    surface_row,
                    bound_mm: bound * mm_per_unit,
                })
            });
        let witness = matches.next().ok_or(Gap::CapNotProven { surface_row })?;
        if matches.next().is_some() {
            return Err(Gap::CapAmbiguous { surface_row });
        }
        witnesses.push(witness);
    }
    let [first, second] = witnesses
        .try_into()
        .map_err(|_| Gap::FaceLoopsUnavailable)?;
    if first.source_edge_id == second.source_edge_id {
        return Err(Gap::SameCapEdgeUsedTwice);
    }
    if !super::cap_trim::complete_ruled_face_trim(
        &loops,
        rows,
        [first.source_edge_id, second.source_edge_id],
        coordinate_tolerance,
    ) {
        return Err(Gap::IncompleteFaceTrim);
    }
    Ok([first, second])
}

/// Read-only STEP/source-context precision audit for complete B-spline strips.
///
/// Cap curves are independently bounded against the two source surface rows.
/// Source uncertainty must be explicit, unit-correct, and subject to the
/// existing revolution uncertainty ceiling; no default precision is assumed.
/// This is not the strict complete-source boundary proof.
///
/// # Errors
/// Invalid STEP structure, missing source solid or uncertainty context, or
/// unresolvable/conflicting source length unit.
#[cfg(feature = "cad-kernel-monstertruck")]
pub fn audit_source_spline_cap_precision(
    source_step: &[u8],
    data_section: usize,
    solid_id: u64,
    requested_bound_mm: Option<f64>,
) -> anyhow::Result<SplineCapPrecisionAudit> {
    use crate::{brep, step_graph::build_index};
    use anyhow::{Context, ensure};
    use monstertruck_io::step::load::ParsedStepSource;
    let exchange = ParsedStepSource::parse(source_step)?.exchange;
    crate::step_io::require_supported_sections(&exchange)?;
    let entities = &exchange
        .data
        .get(data_section)
        .context("missing STEP DATA section")?
        .entities;
    let index = build_index(entities);
    let faces = super::super::strict_closed_shell_faces(solid_id, entities, &index)
        .map_err(anyhow::Error::msg)?;
    let declared = crate::solid_revolutions::declared_source_uncertainty_by_representation_item(
        entities, &index,
    );
    let authored_uncertainty_mm = *declared
        .get(&solid_id)
        .context("missing trusted source-declared precision")?;
    let target_bound_mm = requested_bound_mm.unwrap_or(authored_uncertainty_mm);
    ensure!(
        target_bound_mm.is_finite()
            && target_bound_mm > 0.0
            && target_bound_mm <= authored_uncertainty_mm,
        "requested bound must be positive and no larger than author-declared STEP uncertainty",
    );
    let mm_per_unit = crate::source_edge_locus::source_solid_length_scale_mm(solid_id, entities)?;
    let coordinate_tolerance = target_bound_mm / mm_per_unit;
    ensure!(
        coordinate_tolerance.is_finite() && coordinate_tolerance > 0.0,
        "invalid source units/precision"
    );
    let mut audit = SplineCapPrecisionAudit {
        solid_id,
        authored_uncertainty_mm,
        target_bound_mm,
        strip_faces: 0,
        fully_certified_strip_faces: 0,
        cap_witnesses: Vec::new(),
        unresolved_strip_faces: Vec::new(),
    };
    for face_id in faces {
        let Some(surface_id) = brep::face_surface(face_id, entities, &index) else {
            continue;
        };
        let Some((rows, _)) =
            crate::surface_recovery::boundary_strip_profiles(surface_id, entities, &index)
        else {
            continue;
        };
        audit.strip_faces += 1;
        match verify_source_strip_cap_rows(
            face_id,
            surface_id,
            &rows,
            entities,
            &index,
            coordinate_tolerance,
            mm_per_unit,
        ) {
            Ok([first, second]) => {
                audit.fully_certified_strip_faces += 1;
                audit.cap_witnesses.extend([first, second]);
            }
            Err(reason) => audit.unresolved_strip_faces.push(SplineCapPrecisionGap {
                face_id,
                surface_id,
                reason,
            }),
        }
    }
    audit
        .cap_witnesses
        .sort_unstable_by_key(|w| (w.face_id, w.surface_row));
    audit
        .unresolved_strip_faces
        .sort_unstable_by_key(|gap| gap.face_id);
    Ok(audit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep;
    use monstertruck_io::step::load::ParsedStepSource;

    fn curve(points: &[[f64; 3]], degree: usize, weights: Option<Vec<f64>>) -> BSplineSupport {
        BSplineSupport {
            degree,
            control_points_mm: points.to_vec(),
            knots: std::iter::repeat_n(0.0, degree + 1)
                .chain(std::iter::repeat_n(1.0, degree + 1))
                .collect(),
            weights,
        }
    }

    #[test]
    fn independent_degree_and_rational_speed_preserve_straight_locus() {
        let a = curve(
            &[[0., 0., 0.], [0.2, 0., 0.], [0.7, 0., 0.], [1., 0., 0.]],
            3,
            None,
        );
        let mut b = curve(
            &[[0., 0., 0.], [0.7, 0., 0.], [1., 0., 0.]],
            2,
            Some(vec![1., 0.2, 1.]),
        );
        let proven = certified_coincident_spline_loci(&a, &b, 1.0e-5);
        assert!(proven.is_some_and(|bound| bound <= 1.0e-5), "{proven:?}");
        b.control_points_mm.reverse();
        b.weights = b.weights.map(|mut weights| {
            weights.reverse();
            weights
        });
        assert!(certified_coincident_spline_loci(&a, &b, 1.0e-5).is_some());
        b.control_points_mm[1][1] = 0.001;
        assert!(certified_coincident_spline_loci(&a, &b, 1.0e-5).is_none());
        b.control_points_mm[1][1] = 0.0;
        b.control_points_mm[1][0] = -0.5;
        assert!(certified_coincident_spline_loci(&a, &b, 1.0e-5).is_none());
    }

    #[test]
    fn exact_authored_spline_with_reversed_knots_is_identical_without_subdivision() {
        let a = curve(
            &[[0., 0., 0.], [0.25, 0.5, 0.], [0.75, 0.5, 0.], [1., 0., 0.]],
            3,
            Some(vec![1., 2., 2., 1.]),
        );
        let mut b = a.clone();
        assert_eq!(certified_coincident_spline_loci(&a, &b, 1.0e-9), Some(0.0));
        b.control_points_mm.reverse();
        b.weights = b.weights.map(|mut weights| {
            weights.reverse();
            weights
        });
        assert_eq!(certified_coincident_spline_loci(&a, &b, 1.0e-9), Some(0.0));
        b.control_points_mm[1][1] += 0.01;
        assert!(certified_coincident_spline_loci(&a, &b, 1.0e-6).is_none());
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn original_round_network_precision_qualified_spline_caps() -> anyhow::Result<()> {
        let original = include_bytes!(
            "../../../validation/fixtures/native_drafted_round_network_0b5a3544.step"
        );
        let source = ParsedStepSource::parse(original)?.exchange;
        let [section] = source.data.as_slice() else {
            anyhow::bail!("expected one STEP data section");
        };
        let solids = brep::manifold_solid_ids(&section.entities);
        let [solid] = solids.as_slice() else {
            anyhow::bail!("expected a single source solid");
        };
        let report = audit_source_spline_cap_precision(original, 0, *solid, None)?;
        assert_eq!(report.strip_faces, 13);
        assert_eq!(report.fully_certified_strip_faces, 13);
        assert_eq!(report.cap_witnesses.len(), 26);
        assert_eq!(report.unresolved_strip_faces.len(), 0);
        assert_eq!(
            report.authored_uncertainty_mm.to_bits(),
            1.0e-5_f64.to_bits()
        );
        let maximum = report
            .cap_witnesses
            .iter()
            .map(|e| e.bound_mm)
            .reduce(f64::max)
            .ok_or_else(|| anyhow::anyhow!("missing precision evidence"))?;
        assert!(maximum <= report.authored_uncertainty_mm);

        let strict =
            audit_source_spline_cap_precision(original, 0, *solid, Some(super::super::TOL))?;
        eprintln!(
            "strict boundary quantum={} mm, certified={}/{} caps={}, worst={} mm",
            super::super::TOL,
            strict.fully_certified_strip_faces,
            strict.strip_faces,
            strict.cap_witnesses.len(),
            strict
                .cap_witnesses
                .iter()
                .map(|w| w.bound_mm)
                .fold(0.0_f64, f64::max),
        );
        assert_eq!(strict.fully_certified_strip_faces, 13);
        assert_eq!(strict.cap_witnesses.len(), 26);
        assert!(
            strict
                .cap_witnesses
                .iter()
                .all(|cap| cap.bound_mm <= super::super::TOL)
        );
        // The bound must be meaningful, not a tolerance-independent
        // acceptance of every author-declared STEP curve mismatch.
        let sub_angstrom = audit_source_spline_cap_precision(original, 0, *solid, Some(1.0e-10))?;
        assert_eq!(sub_angstrom.strip_faces, 13);
        assert!(sub_angstrom.fully_certified_strip_faces < strict.fully_certified_strip_faces);
        assert!(audit_source_spline_cap_precision(original, 0, *solid, Some(0.0)).is_err());
        assert!(
            audit_source_spline_cap_precision(
                original,
                0,
                *solid,
                Some(2.0 * report.authored_uncertainty_mm)
            )
            .is_err()
        );

        let original_utf8 = std::str::from_utf8(original)?;
        let micron_precision =
            original_utf8.replacen("LENGTH_MEASURE(0.00001)", "LENGTH_MEASURE(0.000000001)", 1);
        assert_ne!(micron_precision, original_utf8);
        let narrow =
            audit_source_spline_cap_precision(micron_precision.as_bytes(), 0, *solid, None)?;
        assert_eq!(narrow.authored_uncertainty_mm, 1.0e-9);
        assert!(narrow.fully_certified_strip_faces < report.fully_certified_strip_faces);
        assert_eq!(
            narrow.unresolved_strip_faces.len(),
            narrow.strip_faces - narrow.fully_certified_strip_faces
        );
        let absent_context = original_utf8.replacen(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
            "NOT_AN_UNCERTAINTY_CONTEXT",
            1,
        );
        assert!(
            audit_source_spline_cap_precision(absent_context.as_bytes(), 0, *solid, None).is_err()
        );
        let housing =
            include_bytes!("../../../validation/fixtures/native_relief_housing_0914dda4.step");
        let independent = audit_source_spline_cap_precision(housing, 0, 4138, None)?;
        assert_eq!(independent.strip_faces, 112);
        assert_eq!(independent.fully_certified_strip_faces, 112);
        assert!(independent.unresolved_strip_faces.is_empty());
        let oblique =
            include_bytes!("../../../validation/fixtures/native_relief_oblique_8ddfec04.step");
        let independent_oblique = audit_source_spline_cap_precision(oblique, 0, 5442, None)?;
        assert_eq!(independent_oblique.strip_faces, 98);
        assert_eq!(independent_oblique.fully_certified_strip_faces, 98);
        assert!(independent_oblique.unresolved_strip_faces.is_empty());

        eprintln!(
            "manufacturer: {} source strip faces, {} independently bounded caps, max certified {} mm / source {} mm",
            report.fully_certified_strip_faces,
            report.cap_witnesses.len(),
            maximum,
            report.authored_uncertainty_mm
        );
        Ok(())
    }
}
