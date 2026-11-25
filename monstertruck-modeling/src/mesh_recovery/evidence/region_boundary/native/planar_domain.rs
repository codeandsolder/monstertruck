//! Detect certified transverse crossings in planar native trim domains.
//!
//! Carrier incidence and a topological wire do not imply an embedded
//! boundary. In particular, a sparsely witnessed conic section can bow
//! through another edge of the same source patch. Work in a local plane
//! chart using positive-weight rational Bezier hulls. Two curve pieces
//! must cross if each traverses both sides of the other's entire strip:
//! inside the strips' common parallelogram they connect opposite pairs
//! of sides. No tessellation or guessed curve intersection is evidence.
//!
//! This is a positive rejection witness, not a complete simplicity proof.
//! Tangencies, overlapping curves and unresolved near-incidence are not
//! certified by this test.

use crate::{InnerSpace, Point3, Vector2, Vector3, Wire};
use monstertruck_geometry::{
    prelude::TryIntoHomogeneousBsplineCurve, proof_bezier_spans, subdivide_homogeneous_bezier,
};

const MAX_SUBDIVISION: usize = 12;
const MAX_PAIR_VISITS: usize = 131_072;

struct Span {
    controls: Vec<[f64; 4]>,
    bounds: [Vector2; 2],
    ends: [Vector2; 2],
    margin: f64,
}

#[derive(Clone, Copy)]
struct Strip {
    origin: Vector2,
    normal: Vector2,
    range: [f64; 2],
}

fn point(control: [f64; 4]) -> Vector2 {
    Vector2::new(control[0] / control[3], control[1] / control[3])
}

impl Span {
    fn new(controls: Vec<[f64; 4]>, margin: f64) -> Option<Self> {
        let mut lower = Vector2::new(f64::INFINITY, f64::INFINITY);
        let mut upper = Vector2::new(f64::NEG_INFINITY, f64::NEG_INFINITY);
        for &control in &controls {
            if control[3] <= 0.0 || !control.iter().all(|v| v.is_finite()) {
                return None;
            }
            let p = point(control);
            if !p.x.is_finite() || !p.y.is_finite() {
                return None;
            }
            lower.x = lower.x.min(p.x);
            lower.y = lower.y.min(p.y);
            upper.x = upper.x.max(p.x);
            upper.y = upper.y.max(p.y);
        }
        let ends = [point(*controls.first()?), point(*controls.last()?)];
        Some(Self {
            controls,
            bounds: [lower, upper],
            ends,
            margin,
        })
    }

    fn strip(&self) -> Option<Strip> {
        let delta = self.ends[1] - self.ends[0];
        let length = delta.magnitude();
        if !length.is_finite() || length <= self.margin {
            return None;
        }
        let normal = Vector2::new(-delta.y, delta.x) / length;
        let origin = self.ends[0];
        let range = self
            .controls
            .iter()
            .map(|&p| (point(p) - origin).dot(normal))
            .fold([0.0_f64, 0.0_f64], |[lo, hi], value| {
                [lo.min(value), hi.max(value)]
            });
        Some(Strip {
            origin,
            normal,
            range,
        })
    }

    fn overlaps(&self, other: &Self) -> bool {
        let margin = self.margin + other.margin;
        (0..2).all(|axis| {
            self.bounds[0][axis] <= other.bounds[1][axis] + margin
                && other.bounds[0][axis] <= self.bounds[1][axis] + margin
        })
    }

    fn split(&self) -> Option<[Self; 2]> {
        let (left, right) = subdivide_homogeneous_bezier(&self.controls)?;
        Some([
            Self::new(left, self.margin)?,
            Self::new(right, self.margin)?,
        ])
    }

    fn size(&self) -> f64 { (self.bounds[1] - self.bounds[0]).magnitude() }

    // Ordered rational control coordinates imply monotonicity for
    // positive weights, irrespective of the curve's parameter speed.
    fn monotone(&self) -> bool {
        let direction = self.ends[1] - self.ends[0];
        if direction.magnitude() <= self.margin {
            return false;
        }
        self.controls
            .windows(2)
            .all(|pair| (point(pair[1]) - point(pair[0])).dot(direction) >= 0.0)
    }
}

impl Strip {
    fn traversed_by(self, span: &Span, margin: f64) -> bool {
        let [a, b] = span.ends.map(|p| (p - self.origin).dot(self.normal));
        a.min(b) < self.range[0] - margin && a.max(b) > self.range[1] + margin
    }

    fn separates(self, span: &Span, margin: f64) -> bool {
        let [lo, hi] = span
            .controls
            .iter()
            .map(|&p| (point(p) - self.origin).dot(self.normal))
            .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], value| {
                [lo.min(value), hi.max(value)]
            });
        lo > self.range[1] + margin || hi < self.range[0] - margin
    }

    fn width(self) -> f64 { self.range[1] - self.range[0] }
}

fn crossing(a: &Span, b: &Span, depth: usize, remaining: &mut usize) -> bool {
    if *remaining == 0 || !a.overlaps(b) {
        return false;
    }
    *remaining -= 1;
    let margin = a.margin + b.margin;
    let sa = a.strip();
    let sb = b.strip();
    if sa.is_some_and(|s| s.separates(b, margin)) || sb.is_some_and(|s| s.separates(a, margin)) {
        return false;
    }
    if let (Some(sa), Some(sb)) = (sa, sb)
        && sa.traversed_by(b, margin)
        && sb.traversed_by(a, margin)
    {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let width = |span: &Span, strip: Option<Strip>| {
        if span.controls.len() <= 2 {
            0.0
        } else {
            strip.map_or_else(|| span.size(), Strip::width)
        }
    };
    let wa = width(a, sa);
    let wb = width(b, sb);
    if wa.max(wb) <= margin {
        return false;
    }
    let (split, other) = if wa >= wb { (a, b) } else { (b, a) };
    split.split().is_some_and(|pieces| {
        pieces
            .iter()
            .any(|piece| crossing(piece, other, depth - 1, remaining))
    })
}

fn self_crossing(span: &Span, depth: usize, remaining: &mut usize) -> bool {
    if depth == 0
        || *remaining == 0
        || span.monotone()
        || span.controls.len() <= 2
        || span.size() <= span.margin
    {
        return false;
    }
    *remaining -= 1;
    span.split().is_some_and(|[a, b]| {
        crossing(&a, &b, depth - 1, remaining)
            || self_crossing(&a, depth - 1, remaining)
            || self_crossing(&b, depth - 1, remaining)
    })
}

fn projected_spans(
    wires: &[Wire],
    origin: Point3,
    basis: [Vector3; 2],
    tolerance: f64,
) -> Option<Vec<Span>> {
    let mut spans = Vec::new();
    for edge in wires.iter().flatten() {
        let curve = edge.oriented_curve().try_into_homogeneous_bspline_curve()?;
        for span in proof_bezier_spans(&curve)? {
            let maximum = span.controls.iter().map(|p| p[3]).fold(0.0_f64, f64::max);
            let minimum = span
                .controls
                .iter()
                .map(|p| p[3])
                .fold(f64::INFINITY, f64::min);
            let world_scale = span
                .controls
                .iter()
                .map(|p| {
                    let q = Point3::new(p[0] / p[3], p[1] / p[3], p[2] / p[3]);
                    q.x.abs().max(q.y.abs()).max(q.z.abs())
                })
                .fold(
                    origin.x.abs().max(origin.y.abs()).max(origin.z.abs()),
                    f64::max,
                );
            let degree = crate::usize_to_f64(span.controls.len());
            let arithmetic =
                128.0 * f64::EPSILON * world_scale.max(1.0) * (maximum / minimum) * degree * degree;
            if !arithmetic.is_finite() {
                return None;
            }
            let controls = span
                .controls
                .into_iter()
                .map(|p| {
                    let w = p[3] / maximum;
                    let q = Point3::new(p[0] / p[3], p[1] / p[3], p[2] / p[3]) - origin;
                    [q.dot(basis[0]) * w, q.dot(basis[1]) * w, 0.0, w]
                })
                .collect();
            spans.push(Span::new(controls, tolerance.max(arithmetic))?);
        }
    }
    Some(spans)
}

/// Return true only with a transverse crossing witness strictly beyond
/// the given incidence and arithmetic uncertainty.
pub(super) fn has_crossing(
    wires: &[Wire],
    origin: Point3,
    normal: Vector3,
    tolerance: f64,
) -> bool {
    let basis = super::normal_basis(normal).into();
    let Some(mut spans) = projected_spans(wires, origin, basis, tolerance) else {
        return false;
    };
    let mut remaining = MAX_PAIR_VISITS;
    if spans
        .iter()
        .any(|span| self_crossing(span, MAX_SUBDIVISION, &mut remaining))
    {
        return true;
    }
    // Sweep the conservative projected hull bounds. Do not spend
    // quadratic work on pairs with disjoint horizontal intervals.
    spans.sort_by(|a, b| (a.bounds[0].x - a.margin).total_cmp(&(b.bounds[0].x - b.margin)));
    for (i, a) in spans.iter().enumerate() {
        for b in &spans[i + 1..] {
            if remaining == 0 {
                return false;
            }
            if b.bounds[0].x - b.margin > a.bounds[1].x + a.margin {
                break;
            }
            if crossing(a, b, MAX_SUBDIVISION, &mut remaining) {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BsplineCurve, Curve, Edge, EuclideanSpace, KnotVector, NurbsCurve, builder};

    #[test]
    fn curved_planar_trim_crossing_is_not_a_valid_native_domain()
    -> Result<(), Box<dyn std::error::Error>> {
        let [a, b, c] = [
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
        ]
        .map(builder::vertex);
        for (mid, crosses) in [
            (Point3::new(-0.5, 3.0, 0.0), true),
            (Point3::new(-0.5, -0.3, 0.0), false),
        ] {
            let curve = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new(
                KnotVector::bezier_knot(2),
                vec![
                    a.point().to_vec().extend(1.0),
                    mid.to_vec().extend(1.0),
                    b.point().to_vec().extend(1.0),
                ],
            )?));
            let wire = Wire::from(vec![
                Edge::new(&a, &b, curve)?,
                builder::line(&b, &c),
                builder::line(&c, &a),
            ]);
            assert_eq!(
                has_crossing(&[wire], a.point(), Vector3::unit_z(), 1.0e-8),
                crosses
            );
        }
        Ok(())
    }
}
