//! Proof-specific rational B-spline knot-span extraction.
//!
//! Algebraic geometry certificates visit EVERY positive input knot interval,
//! even one narrower than the kernel modeling tolerance.
//! Bernstein coefficients are obtained by the symmetric B-spline blossom,
//! without tolerant knot grouping, snapping, or a generic Cut operation.

use crate::prelude::{BsplineCurve, Vector4};

const MAX_PROOF_DEGREE: usize = 64;
const MAX_PROOF_SPANS: usize = 8192;
const MAX_BLOSSOM_STEPS: usize = 32_000_000;

/// Original positive-length knot span and its homogeneous Bernstein controls.
pub struct ProofBezierSpan {
    /// Inclusive start parameter of this positive-length span.
    pub start: f64,
    /// Inclusive end parameter of this positive-length span.
    pub end: f64,
    /// Homogeneous Bernstein control coefficients in parameter order.
    pub controls: Vec<[f64; 4]>,
}

fn finite_homogeneous(point: Vector4) -> bool {
    [point.x, point.y, point.z, point.w]
        .iter()
        .all(|value| value.is_finite())
        && point.w > 0.0
}

// `total_cmp` is useful for sorting IEEE-754 inputs but treats -0.0
// and +0.0 as different. They are the same mathematical knot value, while
// any two distinct finite nonzero floats must remain different spans.
pub fn same_algebraic_knot(left: f64, right: f64) -> bool {
    left.total_cmp(&right).is_eq() || (left.abs().to_bits() == 0 && right.abs().to_bits() == 0)
}

// A single validated layout and blossom recurrence for ordinary and
// outward-rounded coefficient arithmetic. A surface uses the same recurrence
// successively along both axes, including every positive knot interval.
fn span_layout(knots: &[f64], count: usize) -> Option<(usize, Vec<usize>)> {
    if count < 2 || knots.len() <= count + 1 {
        return None;
    }
    let degree = knots.len().checked_sub(count + 1)?;
    if degree == 0
        || degree >= count
        || degree > MAX_PROOF_DEGREE
        || !knots.iter().all(|knot| knot.is_finite())
        || !knots.windows(2).all(|pair| pair[0] <= pair[1])
        || knots[degree] >= knots[count]
    {
        return None;
    }
    let mut multiplicity = 1_usize;
    for pair in knots.windows(2) {
        if same_algebraic_knot(pair[0], pair[1]) {
            multiplicity += 1;
            if multiplicity > degree + 1 {
                return None;
            }
        } else {
            multiplicity = 1;
        }
    }
    let indices: Vec<_> = (degree..count)
        .filter(|&i| knots[i] < knots[i + 1])
        .collect();
    let steps_per_span = degree
        .checked_mul(degree + 1)?
        .checked_mul(degree + 1)?
        .checked_div(2)?;
    if indices.is_empty()
        || indices.len() > MAX_PROOF_SPANS
        || indices.len().checked_mul(steps_per_span)? > MAX_BLOSSOM_STEPS
        || !same_algebraic_knot(knots[*indices.first()?], knots[degree])
        || !same_algebraic_knot(knots[*indices.last()? + 1], knots[count])
        || indices
            .windows(2)
            .any(|pair| !same_algebraic_knot(knots[pair[0] + 1], knots[pair[1]]))
    {
        return None;
    }
    Some((degree, indices))
}

fn blossom<T: Copy>(
    points: &[T],
    knots: &[f64],
    index: usize,
    degree: usize,
    blend: impl Fn(T, T, f64, f64, f64) -> Option<T>,
) -> Option<Vec<T>> {
    let start = knots[index];
    let end = knots[index + 1];
    let source = &points[index - degree..=index];
    // A clamped Bézier chain already has the desired coefficients. This also
    // preserves their exact floating representation in interval proofs.
    if (0..degree).all(|i| {
        same_algebraic_knot(knots[index - i], start)
            && same_algebraic_knot(knots[index + 1 + i], end)
    }) {
        return Some(source.to_vec());
    }
    let mut controls = Vec::with_capacity(degree + 1);
    for right_repeats in 0..=degree {
        let mut values = source.to_vec();
        for round in 1..=degree {
            let parameter = if round <= degree - right_repeats {
                start
            } else {
                end
            };
            for local in (round..=degree).rev() {
                let left = knots[index - degree + local];
                let right = knots[index + 1 + local - round];
                values[local] = blend(values[local - 1], values[local], parameter, left, right)?;
            }
        }
        controls.push(values[degree]);
    }
    Some(controls)
}

/// Extract the full active domain into exactly partitioned Bernstein pieces.
///
/// None indicates invalid knot/control layout, impossible multiplicity,
/// arithmetic overflow/underflow, or exceeded explicit proof-work budget.
/// This numerical step requires a separate arithmetic error bound when used
/// to certify geometric tolerances.
#[must_use = "all active spans must be available before a whole-domain proof"]
pub fn proof_bezier_spans(curve: &BsplineCurve<Vector4>) -> Option<Vec<ProofBezierSpan>> {
    let points = curve.control_points();
    let knots = curve.knot_vector().as_ref();
    if !points.iter().copied().all(finite_homogeneous) {
        return None;
    }
    let (degree, indices) = span_layout(knots, points.len())?;
    indices
        .into_iter()
        .map(|index| {
            let controls = blossom(
                points,
                knots,
                index,
                degree,
                |a, b, parameter, left, right| {
                    let width = right - left;
                    let numerator = parameter - left;
                    if !width.is_finite() || width <= 0.0 || !numerator.is_finite() {
                        return None;
                    }
                    let alpha = numerator / width;
                    if !alpha.is_finite() || !(0.0..=1.0).contains(&alpha) {
                        return None;
                    }
                    let value = a * (1.0 - alpha) + b * alpha;
                    finite_homogeneous(value).then_some(value)
                },
            )?;
            Some(ProofBezierSpan {
                start: knots[index],
                end: knots[index + 1],
                controls: controls.into_iter().map(|v| [v.x, v.y, v.z, v.w]).collect(),
            })
        })
        .collect()
}

use crate::proof_interval::Interval;

/// Four outward-rounded homogeneous Cartesian coefficients (X, Y, Z, W).
pub type IntervalControl = [Interval; 4];

/// Bernstein coefficient enclosures on one positive original knot interval.
#[derive(Clone, Debug)]
pub struct IntervalBezierSpan {
    /// Bernstein controls in parameter order.
    pub values: Vec<IntervalControl>,
    /// Original parameter interval.
    pub parameter: [f64; 2],
}

impl IntervalBezierSpan {
    /// Subdivide at the exact local midpoint by outward-rounded de Casteljau.
    #[must_use = "both subspans together enclose the complete original span"]
    pub fn split(&self) -> Option<[Self; 2]> {
        let mut work = self.values.clone();
        let mut first = vec![*work.first()?];
        let mut last = vec![*work.last()?];
        for count in (1..work.len()).rev() {
            for i in 0..count {
                work[i] = std::array::from_fn(|k| work[i][k].plus(work[i + 1][k]).scale(0.5));
            }
            first.push(work[0]);
            last.push(work[count - 1]);
        }
        last.reverse();
        let middle = self.parameter[0].midpoint(self.parameter[1]);
        Some([
            Self {
                values: first,
                parameter: [self.parameter[0], middle],
            },
            Self {
                values: last,
                parameter: [middle, self.parameter[1]],
            },
        ])
    }
}

pub fn interval_bezier_spans(
    points: &[IntervalControl],
    knots: &[f64],
) -> Option<Vec<IntervalBezierSpan>> {
    if !points
        .iter()
        .all(|p| p.iter().all(|v| v.finite()) && p[3].lo > 0.0)
    {
        return None;
    }
    let (degree, indices) = span_layout(knots, points.len())?;
    indices
        .into_iter()
        .map(|index| {
            let values = blossom(
                points,
                knots,
                index,
                degree,
                |a, b, parameter, left, right| {
                    let left = Interval::point(left);
                    let alpha = Interval::point(parameter)
                        .minus(left)
                        .checked_div(Interval::point(right).minus(left))?
                        .intersection(Interval::range(0.0, 1.0))?;
                    Some(std::array::from_fn(|i| {
                        let hull = a[i].union(b[i]);
                        a[i].plus(b[i].minus(a[i]).times(alpha))
                            .intersection(hull)
                            .unwrap_or(hull)
                    }))
                },
            )?;
            Some(IntervalBezierSpan {
                values,
                parameter: [knots[index], knots[index + 1]],
            })
        })
        .collect()
}

/// Enclose every active Bézier span of the actual homogeneous curve.
#[must_use = "the interval proof requires all active spans"]
pub fn proof_interval_bezier_spans(
    curve: &BsplineCurve<Vector4>,
) -> Option<Vec<IntervalBezierSpan>> {
    let points: Vec<_> = curve
        .control_points()
        .iter()
        .map(|v| [v.x, v.y, v.z, v.w].map(Interval::point))
        .collect();
    interval_bezier_spans(&points, curve.knot_vector().as_ref())
}

/// Positive rational control-hull bounds shared by curve and surface proofs.
#[derive(Clone, Copy, Debug)]
pub struct RationalControlBounds {
    weight: Interval,
    /// Cartesian coordinate hull of the complete Bernstein domain.
    pub position: [Interval; 3],
}

impl RationalControlBounds {
    /// Bound a positive rational control net without evaluating samples.
    #[must_use = "unresolved rational bounds cannot certify a chart"]
    pub fn from_controls<'a>(
        controls: impl IntoIterator<Item = &'a IntervalControl>,
    ) -> Option<Self> {
        let mut controls = controls.into_iter();
        let first = controls.next()?;
        let mut weight = first[3];
        let mut position = [
            first[0].checked_div(weight)?,
            first[1].checked_div(weight)?,
            first[2].checked_div(weight)?,
        ];
        for p in controls {
            weight = weight.union(p[3]);
            for (i, bound) in position.iter_mut().enumerate() {
                *bound = bound.union(p[i].checked_div(p[3])?);
            }
        }
        (weight.lo > 0.0 && weight.finite()).then_some(Self { weight, position })
    }

    /// Bound a Cartesian derivative from adjacent homogeneous coefficients.
    ///
    /// The derivative uses the local unit Bernstein parameter. A positive
    /// original knot-span width changes the scale, never its orientation.
    #[must_use = "the entire derivative enclosure is needed for orientation proofs"]
    pub fn derivative<'a>(
        self,
        pairs: impl IntoIterator<Item = (&'a IntervalControl, &'a IntervalControl)>,
        degree: usize,
    ) -> Option<[Interval; 3]> {
        let degree = f64::from(u32::try_from(degree).ok()?);
        let mut derivative = pairs
            .into_iter()
            .map(|(a, b)| std::array::from_fn::<_, 4, _>(|i| b[i].minus(a[i]).scale(degree)));
        let mut bounds = derivative.next()?;
        for p in derivative {
            for (bound, value) in bounds.iter_mut().zip(p) {
                *bound = bound.union(value);
            }
        }
        let rate = |i: usize| {
            bounds[i]
                .minus(self.position[i].times(bounds[3]))
                .checked_div(self.weight)
        };
        Some([rate(0)?, rate(1)?, rate(2)?])
    }
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
    fn signed_zero_is_one_algebraic_knot_for_multiplicity() {
        // Degree one permits at most two repetitions of any knot value.
        // -0.0 and +0.0 are algebraically the same parameter, despite
        // total_cmp distinguishing their IEEE encodings for sorting.
        let invalid = BsplineCurve::new_unchecked(
            KnotVector::from(vec![-0.0, -0.0, 0.0, 0.0, 1.0, 1.0]),
            vec![p(0.0, 0.0, 0.0, 1.0); 4],
        );
        assert!(proof_bezier_spans(&invalid).is_none());
    }

    #[test]
    fn algebraically_equal_signed_zero_knots_preserve_both_real_spans() {
        // A valid degree-one curve with exactly two interior knots at zero;
        // the mixed signs are a harmless floating-point representation.
        let curve = BsplineCurve::new_unchecked(
            KnotVector::from(vec![-1.0, -1.0, -0.0, 0.0, 1.0, 1.0]),
            vec![
                p(0.0, 0.0, 0.0, 1.0),
                p(0.25, 0.0, 0.0, 1.0),
                p(0.75, 0.0, 0.0, 1.0),
                p(1.0, 0.0, 0.0, 1.0),
            ],
        );
        let spans = proof_bezier_spans(&curve).expect("valid -0/+0 join must be covered");
        assert_eq!(spans.len(), 2);
        assert_eq!(spans[0].start.total_cmp(&-1.0), std::cmp::Ordering::Equal);
        assert_eq!(spans[1].end.total_cmp(&1.0), std::cmp::Ordering::Equal);
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
