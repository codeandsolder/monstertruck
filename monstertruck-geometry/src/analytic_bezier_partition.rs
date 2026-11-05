//! Proof-specific rational B-spline knot-span extraction.
//!
//! The generic modeling decomposition groups nearby knots with a tolerance.
//! Algebraic geometry certificates must instead visit EVERY positive input
//! knot interval, even one narrower than the kernel modeling tolerance.
//! Bernstein coefficients are obtained by the symmetric B-spline blossom,
//! without tolerant knot grouping, snapping, or a generic Cut operation.

use crate::prelude::{BsplineCurve, Vector4};

const MAX_PROOF_DEGREE: usize = 64;
const MAX_PROOF_SPANS: usize = 8192;
const MAX_BLOSSOM_STEPS: usize = 32_000_000;

/// Original positive-length knot span and its homogeneous Bernstein controls.
pub(super) struct ProofBezierSpan {
    pub(super) start: f64,
    pub(super) end: f64,
    pub(super) controls: Vec<[f64; 4]>,
}

fn finite_homogeneous(point: Vector4) -> bool {
    [point.x, point.y, point.z, point.w]
        .iter()
        .all(|value| value.is_finite())
        && point.w > 0.0
}

/// Extract the full active domain into exactly partitioned Bernstein pieces.
///
/// None indicates invalid knot/control layout, impossible multiplicity,
/// arithmetic overflow/underflow, or exceeded explicit proof-work budget.
/// This numerical step requires a separate arithmetic error bound when used
/// to certify geometric tolerances.
pub(super) fn proof_bezier_spans(curve: &BsplineCurve<Vector4>) -> Option<Vec<ProofBezierSpan>> {
    let points = curve.control_points();
    let knots = curve.knot_vector();
    let count = points.len();
    if count < 2 || knots.len() <= count + 1 {
        return None;
    }
    let degree = knots.len().checked_sub(count + 1)?;
    if degree == 0 || degree >= count || degree > MAX_PROOF_DEGREE {
        return None;
    }
    if !points.iter().copied().all(finite_homogeneous)
        || !knots.iter().all(|knot| knot.is_finite())
        || !knots.iter().zip(knots.iter().skip(1)).all(|(a, b)| a <= b)
        || knots[degree] >= knots[count]
    {
        return None;
    }
    // Only genuinely equal knots contribute to multiplicity.
    let mut multiplicity = 1_usize;
    for (left, right) in knots.iter().zip(knots.iter().skip(1)) {
        if left.total_cmp(right).is_eq() {
            multiplicity += 1;
            if multiplicity > degree + 1 {
                return None;
            }
        } else {
            multiplicity = 1;
        }
    }

    let nonempty = (degree..count)
        .filter(|&index| knots[index] < knots[index + 1])
        .count();
    if nonempty == 0 || nonempty > MAX_PROOF_SPANS {
        return None;
    }
    let steps_per_span = degree
        .checked_mul(degree + 1)?
        .checked_mul(degree + 1)?
        .checked_div(2)?;
    if nonempty.checked_mul(steps_per_span)? > MAX_BLOSSOM_STEPS {
        return None;
    }

    let mut spans = Vec::with_capacity(nonempty);
    for index in degree..count {
        let start = knots[index];
        let end = knots[index + 1];
        if start >= end {
            continue;
        }
        let source = &points[index - degree..=index];
        let mut controls = Vec::with_capacity(degree + 1);
        for right_repeats in 0..=degree {
            let mut values = source.to_vec();
            // Blossom: degree-right_repeats arguments at start and the
            // remaining arguments at end. Every knot denominator brackets
            // this positive interval and therefore is strictly positive.
            for round in 1..=degree {
                let parameter = if round <= degree - right_repeats {
                    start
                } else {
                    end
                };
                for local in (round..=degree).rev() {
                    let left_knot = knots[index - degree + local];
                    let right_knot = knots[index + 1 + local - round];
                    let width = right_knot - left_knot;
                    let numerator = parameter - left_knot;
                    if !width.is_finite() || width <= 0.0 || !numerator.is_finite() {
                        return None;
                    }
                    let alpha = numerator / width;
                    if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
                        return None;
                    }
                    values[local] = values[local - 1] * (1.0 - alpha) + values[local] * alpha;
                    if !finite_homogeneous(values[local]) {
                        return None;
                    }
                }
            }
            let value = values[degree];
            controls.push([value.x, value.y, value.z, value.w]);
        }
        spans.push(ProofBezierSpan {
            start,
            end,
            controls,
        });
    }
    // A certificate must never omit an active knot interval, however short.
    if spans
        .first()
        .is_none_or(|span| span.start.total_cmp(&knots[degree]).is_ne())
        || spans
            .last()
            .is_none_or(|span| span.end.total_cmp(&knots[count]).is_ne())
        || spans
            .windows(2)
            .any(|pair| pair[0].end.total_cmp(&pair[1].start).is_ne())
    {
        return None;
    }
    Some(spans)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::{InnerSpace, KnotVector, NurbsCurve, ParametricCurve};

    fn p(x: f64, y: f64, z: f64, w: f64) -> Vector4 { Vector4::new(x, y, z, w) }

    #[test]
    fn single_bezier_preserves_homogeneous_controls() {
        let source = vec![
            p(0.0, 0.0, 0.0, 1.0),
            p(2.0, 1.0, 0.0, 2.0),
            p(3.0, 0.0, 0.0, 1.0),
        ];
        let curve = BsplineCurve::new_unchecked(KnotVector::bezier_knot(2), source.clone());
        let spans = proof_bezier_spans(&curve).expect("valid span");
        assert_eq!(spans.len(), 1);
        for (actual, expected) in spans[0].controls.iter().zip(source.iter()) {
            assert_eq!(*actual, [expected.x, expected.y, expected.z, expected.w]);
        }
    }

    #[test]
    fn near_distinct_knots_are_not_grouped() {
        let narrow = 1.0e-9;
        let curve = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5 + narrow, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0, 1.0),
                p(0.25, 0.0, 0.0, 1.0),
                p(0.75, 0.0, 0.0, 1.0),
                p(1.0, 0.0, 0.0, 1.0),
            ],
        );
        let spans = proof_bezier_spans(&curve).expect("short interval retained");
        assert_eq!(spans.len(), 3);
        assert_eq!(spans[1].start.total_cmp(&0.5), std::cmp::Ordering::Equal);
        assert_eq!(
            spans[1].end.total_cmp(&(0.5 + narrow)),
            std::cmp::Ordering::Equal
        );
    }

    #[test]
    fn nonclamped_quadratic_matches_original_at_interior_parameters() {
        let source = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.1, 0.2, 0.5, 0.75, 1.0, 1.1]),
            vec![
                p(0.0, 0.0, 0.0, 1.0),
                p(1.0, 2.0, 0.0, 1.0),
                p(2.0, -1.0, 0.0, 1.0),
                p(3.0, 0.0, 0.0, 1.0),
            ],
        );
        let original = NurbsCurve::new(source.clone());
        let spans = proof_bezier_spans(&source).expect("valid nonclamped");
        assert_eq!(spans.len(), 2);
        for span in spans {
            let bezier = NurbsCurve::new(BsplineCurve::new_unchecked(
                KnotVector::from(
                    vec![span.start; 3]
                        .into_iter()
                        .chain(vec![span.end; 3])
                        .collect::<Vec<_>>(),
                ),
                span.controls
                    .into_iter()
                    .map(|[x, y, z, w]| p(x, y, z, w))
                    .collect(),
            ));
            for index in 1..10 {
                let t = span.start + (span.end - span.start) * f64::from(index) / 10.0;
                assert!((original.evaluate(t) - bezier.evaluate(t)).magnitude() < 1.0e-9);
            }
        }
    }

    // Independent exact-knot de Boor evaluator, for differential tests only.
    // Proof code itself never compares sampled evaluations.
    fn reference_homogeneous(curve: &BsplineCurve<Vector4>, parameter: f64) -> Vector4 {
        let knots = curve.knot_vector();
        let count = curve.control_points().len();
        let degree = knots.len() - count - 1;
        let index = (degree..count)
            .find(|&index| knots[index] <= parameter && parameter < knots[index + 1])
            .expect("test value lies in positive active knot span");
        let mut controls = curve.control_points()[index - degree..=index].to_vec();
        for round in 1..=degree {
            for local in (round..=degree).rev() {
                let left = knots[index - degree + local];
                let right = knots[index + 1 + local - round];
                let alpha = (parameter - left) / (right - left);
                controls[local] = controls[local - 1] * (1.0 - alpha) + controls[local] * alpha;
            }
        }
        controls[degree]
    }

    #[test]
    fn narrow_spans_and_repeated_knots_match_independent_de_boor() {
        let curves = [
            BsplineCurve::new_unchecked(
                KnotVector::from(vec![0.0, 0.0, 0.5, 0.5 + 1.0e-9, 1.0, 1.0]),
                (0..4)
                    .map(|i| p(f64::from(i), 1.0, 0.1, 1.0 + 0.1 * f64::from(i)))
                    .collect(),
            ),
            BsplineCurve::new_unchecked(
                KnotVector::from(vec![
                    0.0,
                    0.0,
                    0.0,
                    0.3,
                    0.3,
                    0.3 + 1.0e-9,
                    0.6,
                    1.0,
                    1.0,
                    1.0,
                ]),
                (0..7)
                    .map(|i| {
                        p(
                            0.2 * f64::from(i),
                            (f64::from(i) - 3.0).powi(2),
                            0.03 * f64::from(i),
                            1.0 + 0.2 * f64::from(i),
                        )
                    })
                    .collect(),
            ),
        ];
        for curve in curves {
            let pieces = proof_bezier_spans(&curve).expect("all positive knot spans retained");
            for span in pieces {
                let degree = span.controls.len() - 1;
                let bezier = BsplineCurve::new_unchecked(
                    KnotVector::from(
                        vec![span.start; degree + 1]
                            .into_iter()
                            .chain(vec![span.end; degree + 1])
                            .collect::<Vec<_>>(),
                    ),
                    span.controls
                        .into_iter()
                        .map(|[x, y, z, w]| p(x, y, z, w))
                        .collect(),
                );
                for fraction in [0.125, 0.375, 0.625, 0.875] {
                    let t = (span.end - span.start).mul_add(fraction, span.start);
                    let reference = reference_homogeneous(&curve, t);
                    let extracted = reference_homogeneous(&bezier, t);
                    let error = (reference - extracted).magnitude();
                    assert!(
                        error <= 5.0e-9 * reference.magnitude().max(1.0),
                        "failed for t={t}, difference={error}"
                    );
                }
            }
        }
    }

    #[test]
    fn invalid_multiplicity_and_weight_refuse() {
        let invalid = BsplineCurve::new_unchecked(
            KnotVector::from(vec![0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0]),
            vec![p(0.0, 0.0, 0.0, 1.0); 5],
        );
        assert!(proof_bezier_spans(&invalid).is_none());
        let bad_weight = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![p(0.0, 0.0, 0.0, 0.0), p(1.0, 0.0, 0.0, 1.0)],
        );
        assert!(proof_bezier_spans(&bad_weight).is_none());
    }
}
