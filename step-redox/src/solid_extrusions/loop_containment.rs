//! Geometric role classification for unlabelled planar wires.
//!
//! A positive-weight rational Bernstein span is contained in its Euclidean
//! control hull. Separate every pair of authored wires by subdividing those
//! enclosures, then classify a boundary witness of one wire against the
//! other. The witness cannot cross that boundary. Replacing a span by its
//! chord inside an enclosure disjoint from the witness preserves winding
//! parity; these are certified homotopies, not sampled curve approximations.
//! This classifies roles of the authored closed wires, not their simplicity
//! or occupied solid. Those checks remain part of the native and final
//! complete-boundary proofs.

use super::{RecoveredProfileCurve, profile_contains_with_rotated_ray};
use crate::math2::distance;
use monstertruck_geometry::{
    prelude::{BsplineCurve, KnotVector, Vector4},
    proof_bezier_spans, subdivide_homogeneous_bezier,
};
use num_traits::ToPrimitive;
use std::f64::consts::{FRAC_PI_2, TAU};

type Controls = Vec<[f64; 4]>;

const MAX_DEPTH: usize = 48;
const MAX_WORK: usize = 1_000_000;

#[derive(Clone, Copy)]
struct Bounds {
    lo: [f64; 2],
    hi: [f64; 2],
}

impl Bounds {
    fn include(self, other: Self) -> Self {
        Self {
            lo: std::array::from_fn(|i| self.lo[i].min(other.lo[i])),
            hi: std::array::from_fn(|i| self.hi[i].max(other.hi[i])),
        }
    }

    fn separate(self, other: Self) -> bool {
        (0..2).any(|i| self.hi[i] < other.lo[i] || other.hi[i] < self.lo[i])
    }

    fn excludes(self, point: [f64; 2]) -> bool {
        (0..2).any(|i| point[i] < self.lo[i] || point[i] > self.hi[i])
    }

    fn diameter(self) -> f64 { distance(self.lo, self.hi) }
}

struct Span {
    controls: Controls,
    bounds: Bounds,
    // Retain the original curve's numerical conditioning allowance through
    // every subdivision, rather than recomputing it from smaller hulls.
    resolution: f64,
    depth: usize,
}

fn project(p: &[f64; 4]) -> Option<[f64; 2]> {
    let xy = [p[0] / p[3], p[1] / p[3]];
    (p[3] > 0.0 && xy.iter().all(|v| v.is_finite())).then_some(xy)
}

fn homogeneous(p: [f64; 2], weight: f64) -> [f64; 4] { [p[0] * weight, p[1] * weight, 0.0, weight] }

impl Span {
    fn new(controls: Controls, resolution: f64, depth: usize, tolerance: f64) -> Option<Self> {
        if controls.len() < 2 || depth > MAX_DEPTH {
            return None;
        }
        let roundoff = resolution * (depth + 1).to_f64()?;
        if !roundoff.is_finite() || roundoff >= tolerance {
            return None;
        }
        let pad = tolerance + roundoff;
        let mut bounds = Bounds {
            lo: [f64::INFINITY; 2],
            hi: [f64::NEG_INFINITY; 2],
        };
        for p in &controls {
            let xy = project(p)?;
            bounds = bounds.include(Bounds {
                lo: xy.map(|x| (x - pad).next_down()),
                hi: xy.map(|x| (x + pad).next_up()),
            });
        }
        bounds
            .lo
            .iter()
            .chain(&bounds.hi)
            .all(|x| x.is_finite())
            .then_some(Self {
                controls,
                bounds,
                resolution,
                depth,
            })
    }

    fn endpoints(&self) -> Option<[[f64; 2]; 2]> {
        Some([
            project(self.controls.first()?)?,
            project(self.controls.last()?)?,
        ])
    }

    fn split(&self, tolerance: f64) -> Option<[Self; 2]> {
        let (left, right) = subdivide_homogeneous_bezier(&self.controls)?;
        Some([
            Self::new(left, self.resolution, self.depth + 1, tolerance)?,
            Self::new(right, self.resolution, self.depth + 1, tolerance)?,
        ])
    }
}

/// Use the same numerical-conditioning convention as the full-domain
/// analytic incidence proofs. A large coordinate/weight ratio can only
/// make the enclosure test inconclusive, never admit a closer boundary.
fn resolution(controls: &[[f64; 4]]) -> Option<f64> {
    let mut min_weight = f64::INFINITY;
    let mut max_weight = 0.0_f64;
    let mut scale = 1.0_f64;
    for p in controls {
        let xy = project(p)?;
        min_weight = min_weight.min(p[3]);
        max_weight = max_weight.max(p[3]);
        scale = scale.max(xy[0].abs()).max(xy[1].abs());
    }
    let result = 256.0
        * f64::EPSILON
        * scale
        * (max_weight / min_weight).powi(2)
        * controls.len().to_f64()?.powi(2);
    result.is_finite().then_some(result)
}

fn polynomial_spans(points: &[[f64; 2]]) -> Option<(Vec<Controls>, f64)> {
    if points.len() < 2 {
        return None;
    }
    let controls = points
        .iter()
        .map(|&p| homogeneous(p, 1.0))
        .collect::<Controls>();
    let margin = resolution(&controls)?;
    Some((vec![controls], margin))
}

fn spline_spans(
    points: &[[f64; 2]],
    degree: usize,
    knots: &[f64],
    weights: Option<&[f64]>,
) -> Option<(Vec<Controls>, f64)> {
    let count = points.len();
    if degree == 0
        || degree >= count
        || knots.len() != count.checked_add(degree)?.checked_add(1)?
        || weights
            .is_some_and(|w| w.len() != count || w.iter().any(|w| !w.is_finite() || *w <= 0.0))
    {
        return None;
    }
    // Authored profile endpoints designate the clamped active domain.
    let same_knot = |a: f64, b: f64| {
        a.total_cmp(&b).is_eq() || (a.abs().to_bits() == 0 && b.abs().to_bits() == 0)
    };
    if !knots[..degree].iter().all(|&k| same_knot(k, knots[degree]))
        || !knots[count..].iter().all(|&k| same_knot(k, knots[count]))
    {
        return None;
    }
    let max_weight = weights.map_or(1.0, |w| w.iter().copied().fold(0.0, f64::max));
    let controls = points
        .iter()
        .enumerate()
        .map(|(i, &p)| homogeneous(p, weights.map_or(1.0, |w| w[i] / max_weight)))
        .collect::<Controls>();
    // Include blossom extraction conditioning as well as subdivision.
    // Repeated knots are exact; only positive active intervals introduce
    // interpolation denominators.
    let domain = knots[count] - knots[degree];
    if !domain.is_finite() || domain <= 0.0 {
        return None;
    }
    let min_span = knots
        .windows(2)
        .filter_map(|w| (w[1] > w[0]).then_some((w[1] - w[0]) / domain))
        .reduce(f64::min)?;
    let margin = resolution(&controls)? / min_span;
    let curve = BsplineCurve::new(
        KnotVector::try_from(knots.to_vec()).ok()?,
        controls
            .iter()
            .map(|&[x, y, z, w]| Vector4::new(x, y, z, w))
            .collect(),
    )
    .ok()?;
    Some((
        proof_bezier_spans(&curve)?
            .into_iter()
            .map(|span| span.controls)
            .collect(),
        margin,
    ))
}

fn arc_spans(center: [f64; 2], radius: f64, start: f64, end: f64) -> Option<(Vec<Controls>, f64)> {
    let sweep = end - start;
    if !radius.is_finite() || radius <= 0.0 || !sweep.is_finite() || sweep.abs() > TAU {
        return None;
    }
    let count = (sweep.abs() / FRAC_PI_2).ceil().to_u32()?;
    if count == 0 {
        return None;
    }
    let point = |t: f64| {
        [
            radius.mul_add(t.cos(), center[0]),
            radius.mul_add(t.sin(), center[1]),
        ]
    };
    let mut result = Vec::new();
    let mut margin = 0.0_f64;
    for i in 0..count {
        let a = sweep.mul_add(f64::from(i) / f64::from(count), start);
        let b = if i + 1 == count {
            end
        } else {
            sweep.mul_add(f64::from(i + 1) / f64::from(count), start)
        };
        let middle = a.midpoint(b);
        let w = ((b - a) / 2.0).cos();
        let controls = vec![
            homogeneous(point(a), 1.0),
            [
                center[0].mul_add(w, radius * middle.cos()),
                center[1].mul_add(w, radius * middle.sin()),
                0.0,
                w,
            ],
            homogeneous(point(b), 1.0),
        ];
        margin = margin.max(resolution(&controls)?);
        result.push(controls);
    }
    Some((result, margin))
}

fn curve_spans(curve: &RecoveredProfileCurve, tolerance: f64) -> Option<Vec<Span>> {
    let (pieces, margin) = match curve {
        RecoveredProfileCurve::Line {
            start_mm, end_mm, ..
        } => polynomial_spans(&[*start_mm, *end_mm])?,
        RecoveredProfileCurve::CircleArc {
            center_mm,
            radius_mm,
            start_angle_rad,
            end_angle_rad,
            ..
        } => arc_spans(*center_mm, *radius_mm, *start_angle_rad, *end_angle_rad)?,
        RecoveredProfileCurve::Bezier {
            control_points_mm, ..
        } => polynomial_spans(control_points_mm)?,
        RecoveredProfileCurve::BSpline {
            degree,
            control_points_mm,
            knots,
            weights,
            ..
        } => spline_spans(control_points_mm, *degree, knots, weights.as_deref())?,
    };
    let spans = pieces
        .into_iter()
        .map(|p| Span::new(p, margin, 0, tolerance))
        .collect::<Option<Vec<_>>>()?;
    if distance(spans.first()?.endpoints()?[0], curve.start_point()?) > tolerance
        || distance(spans.last()?.endpoints()?[1], curve.end_point()?) > tolerance
    {
        return None;
    }
    Some(spans)
}

struct Wire {
    spans: Vec<Span>,
    bounds: Bounds,
    witness: [f64; 2],
}

impl Wire {
    fn new(curves: &[RecoveredProfileCurve], tolerance: f64) -> Option<Self> {
        let spans = curves
            .iter()
            .map(|c| curve_spans(c, tolerance))
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();
        let witness = spans.first()?.endpoints()?[0];
        for (i, span) in spans.iter().enumerate() {
            if distance(
                span.endpoints()?[1],
                spans[(i + 1) % spans.len()].endpoints()?[0],
            ) > tolerance
            {
                return None;
            }
        }
        let bounds = spans
            .iter()
            .map(|span| span.bounds)
            .reduce(Bounds::include)?;
        Some(Self {
            spans,
            bounds,
            witness,
        })
    }
}

/// Every recursive pair is covered by Bernstein enclosures of the original
/// spans; touching, crossing and unresolved pairs do not establish separation.
fn separated(a: &Span, b: &Span, tolerance: f64, remaining: &mut usize) -> Option<()> {
    *remaining = remaining.checked_sub(1)?;
    if a.bounds.separate(b.bounds) {
        return Some(());
    }
    // Choose the larger enclosure, avoiding repeated refinement of the
    // smaller operand when a long straight span meets several short arcs.
    if a.bounds.diameter() >= b.bounds.diameter() {
        for child in a.split(tolerance)? {
            separated(&child, b, tolerance, remaining)?;
        }
    } else {
        for child in b.split(tolerance)? {
            separated(a, &child, tolerance, remaining)?;
        }
    }
    Some(())
}

fn append_chords(
    span: &Span,
    witness: [f64; 2],
    tolerance: f64,
    remaining: &mut usize,
    chords: &mut Vec<RecoveredProfileCurve>,
) -> Option<()> {
    *remaining = remaining.checked_sub(1)?;
    if span.bounds.excludes(witness) {
        let [start_mm, end_mm] = span.endpoints()?;
        if distance(start_mm, end_mm) > tolerance {
            chords.push(RecoveredProfileCurve::Line {
                source_edge_ids: Vec::new(),
                start_mm,
                end_mm,
            });
        }
    } else {
        for child in span.split(tolerance)? {
            append_chords(&child, witness, tolerance, remaining, chords)?;
        }
    }
    Some(())
}

fn contains(wire: &Wire, witness: [f64; 2], tolerance: f64, remaining: &mut usize) -> Option<bool> {
    if wire.bounds.excludes(witness) {
        return Some(false);
    }
    let mut chords = Vec::new();
    for span in &wire.spans {
        append_chords(span, witness, tolerance, remaining, &mut chords)?;
    }
    profile_contains_with_rotated_ray(&chords, witness, tolerance)
}

/// Prove one outer boundary enclosing all other (mutually nonnesting) wires.
/// Authored bound order, winding and surface sense play no role.
pub fn prove_loop_roles(loops: &[Vec<RecoveredProfileCurve>], tolerance: f64) -> Option<Vec<bool>> {
    if loops.is_empty() || !tolerance.is_finite() || tolerance <= 0.0 {
        return None;
    }
    let wires = loops
        .iter()
        .map(|ring| Wire::new(ring, tolerance))
        .collect::<Option<Vec<_>>>()?;
    let mut enclosed = vec![0; wires.len()];
    let mut remaining = MAX_WORK;
    for (i, a) in wires.iter().enumerate() {
        for (j, b) in wires.iter().enumerate().skip(i + 1) {
            if a.bounds.separate(b.bounds) {
                continue;
            }
            for first in &a.spans {
                for second in &b.spans {
                    separated(first, second, tolerance, &mut remaining)?;
                }
            }
            let ab = contains(a, b.witness, tolerance, &mut remaining)?;
            let ba = contains(b, a.witness, tolerance, &mut remaining)?;
            if ab && ba {
                return None;
            }
            enclosed[i] += usize::from(ab);
            enclosed[j] += usize::from(ba);
        }
    }
    let roles = enclosed
        .iter()
        .map(|&n| n == wires.len() - 1)
        .collect::<Vec<_>>();
    (roles.iter().filter(|&&v| v).count() == 1
        && enclosed
            .iter()
            .zip(&roles)
            .all(|(&n, &outer)| outer || n == 0))
    .then_some(roles)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn polygon(points: &[[f64; 2]]) -> Vec<RecoveredProfileCurve> {
        points
            .iter()
            .zip(points.iter().cycle().skip(1))
            .map(|(&start_mm, &end_mm)| RecoveredProfileCurve::Line {
                source_edge_ids: vec![],
                start_mm,
                end_mm,
            })
            .collect()
    }
    fn spline_hole(x: f64, y: f64) -> Vec<RecoveredProfileCurve> {
        vec![
            RecoveredProfileCurve::BSpline {
                source_edge_ids: vec![],
                degree: 2,
                control_points_mm: vec![[x, y], [x + 0.25, y - 0.5], [x + 0.5, y]],
                knots: vec![0., 0., 0., 1., 1., 1.],
                weights: Some(vec![1., 0.8, 1.]),
            },
            RecoveredProfileCurve::Line {
                source_edge_ids: vec![],
                start_mm: [x + 0.5, y],
                end_mm: [x, y],
            },
        ]
    }

    #[test]
    fn nonconvex_outer_with_unlabelled_spline_holes_ignores_winding() {
        let outer = polygon(&[
            [0., 0.],
            [4., 0.],
            [4., 4.],
            [2., 4.],
            [2., 2.],
            [1., 2.],
            [1., 4.],
            [0., 4.],
        ]);
        let first = spline_hole(0.1, 1.);
        let second = spline_hole(2.5, 3.);
        assert_eq!(
            prove_loop_roles(&[second.clone(), outer.clone(), first.clone()], 1.0e-8),
            Some(vec![false, true, false])
        );
        let reversed = outer
            .iter()
            .rev()
            .map(RecoveredProfileCurve::reversed)
            .collect::<Vec<_>>();
        assert_eq!(
            prove_loop_roles(&[first, reversed, second], 1.0e-8),
            Some(vec![false, true, false])
        );
    }

    #[test]
    fn separated_curved_outer_and_holes_refuse_touching_crossing_and_nesting() {
        let circle = |radius_mm| {
            vec![RecoveredProfileCurve::CircleArc {
                source_edge_ids: vec![],
                center_mm: [0., 0.],
                radius_mm,
                start_angle_rad: 0.,
                end_angle_rad: TAU,
            }]
        };
        assert_eq!(
            prove_loop_roles(&[circle(1.), circle(3.)], 1.0e-8),
            Some(vec![false, true])
        );
        assert_eq!(
            prove_loop_roles(&[circle(1.), circle(2.), circle(3.)], 1.0e-8),
            None
        );
        assert_eq!(prove_loop_roles(&[circle(2.), circle(2.)], 1.0e-8), None);
        let crossing = polygon(&[[-3., -0.5], [3., -0.5], [3., 0.5], [-3., 0.5]]);
        assert_eq!(prove_loop_roles(&[circle(2.), crossing], 1.0e-8), None);
    }
}
