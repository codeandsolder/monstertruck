//! Positive simplicity certificate for projected rational trim domains.
//!
//! Unlike a crossing witness, a failed separation is inconclusive. The
//! entire curve-to-chord homotopy must preserve disjointness, orientation
//! and nesting. Only endpoint alignment within an explicitly charged seam
//! bound is permitted; no native geometry is mutated.

use super::{Span, normal_basis, point, projected_spans};
use crate::{Edge, InnerSpace, Point3, Vector2, Vector3, Wire};
use robust::{Coord, orient2d};

const MAX_DEPTH: usize = 6;
const MAX_PIECES: usize = 8192;
const MAX_PAIR_TESTS: usize = 262_144;
const MAX_AXIS_TESTS: usize = 1_048_576;

/// Positive simplicity and nesting evidence in an oriented planar chart.
///
/// The proof concerns endpoint-aligned projected rational curves; their
/// discrepancy from the original projected curves is bounded separately.
#[derive(Debug, Clone, Copy)]
pub struct PlanarDomainCertificate {
    /// The unique positively oriented outer loop.
    pub outer_wire: usize,
    /// Number of monotone rational Bézier pieces in the certificate.
    pub pieces: usize,
    /// Maximum endpoint alignment and projection arithmetic uncertainty.
    pub seam_displacement_upper_bound: f64,
}

struct Piece {
    controls: Vec<Vector2>,
    bounds: [Vector2; 2],
    margin: f64,
    wire: usize,
    index: usize,
    count: usize,
}

const fn coord(p: Vector2) -> Coord<f64> { Coord { x: p.x, y: p.y } }

fn turn(a: Vector2, b: Vector2, p: Vector2) -> f64 { orient2d(coord(a), coord(b), coord(p)) }

fn same_point(a: Vector2, b: Vector2) -> bool {
    (a.x - b.x).abs() == 0.0 && (a.y - b.y).abs() == 0.0
}

impl Piece {
    fn new(
        span: &Span,
        ends: [Vector2; 2],
        wire: usize,
        index: usize,
        count: usize,
        depth: usize,
    ) -> Option<Self> {
        let mut controls = span.controls.iter().copied().map(point).collect::<Vec<_>>();
        *controls.first_mut()? = ends[0];
        *controls.last_mut()? = ends[1];
        let bounds = controls.iter().fold([ends[0], ends[0]], |[lo, hi], p| {
            [
                Vector2::new(lo.x.min(p.x), lo.y.min(p.y)),
                Vector2::new(hi.x.max(p.x), hi.y.max(p.y)),
            ]
        });
        let margin = span.margin * crate::usize_to_f64(depth + 1);
        let result = Self {
            controls,
            bounds,
            margin,
            wire,
            index,
            count,
        };
        result.monotone().then_some(result)
    }

    fn monotone(&self) -> bool {
        let Some(&first) = self.controls.first() else {
            return false;
        };
        let Some(&last) = self.controls.last() else {
            return false;
        };
        let direction = last - first;
        if direction.magnitude() <= 2.0 * self.margin {
            return false;
        }
        // Positive rational weights preserve strict scalar control order.
        // Require a margin for each difference, not just the endpoints:
        // a rounded zero alone cannot certify monotonicity.
        (0..2).any(|axis| {
            let sign = direction[axis].signum();
            direction[axis].abs() > 2.0 * self.margin
                && self
                    .controls
                    .windows(2)
                    .all(|pair| sign * (pair[1][axis] - pair[0][axis]) > 2.0 * self.margin)
        }) || self.controls.windows(2).all(|pair| {
            (pair[1] - pair[0]).dot(direction) > 2.0 * self.margin * direction.magnitude()
        })
    }

    fn end(&self) -> Vector2 { self.controls[self.controls.len() - 1] }

    const fn next_is(&self, other: &Self) -> bool {
        self.wire == other.wire && (self.index + 1) % self.count == other.index
    }

    fn overlaps(&self, other: &Self) -> bool {
        let margin = self.margin + other.margin;
        (0..2).all(|axis| {
            self.bounds[0][axis] <= other.bounds[1][axis] + margin
                && other.bounds[0][axis] <= self.bounds[1][axis] + margin
        })
    }

    fn interval(&self, origin: Vector2, axis: Vector2) -> [f64; 2] {
        self.controls
            .iter()
            .map(|&p| (p - origin).dot(axis))
            .fold([f64::INFINITY, f64::NEG_INFINITY], |[lo, hi], v| {
                [lo.min(v), hi.max(v)]
            })
    }

    fn separated(&self, other: &Self, remaining_axes: &mut usize) -> bool {
        if !self.overlaps(other) {
            return true;
        }
        let origin = self.controls[0];
        let margin = self.margin + other.margin;
        let mut separates = |axis: Vector2| {
            let Some(remaining) = remaining_axes.checked_sub(1) else {
                return false;
            };
            *remaining_axes = remaining;
            let len = axis.magnitude();
            if len <= margin || !len.is_finite() {
                return false;
            }
            let axis = axis / len;
            let [lo, hi] = self.interval(origin, axis);
            let [other_lo, other_hi] = other.interval(origin, axis);
            hi + margin < other_lo || other_hi + margin < lo
        };
        // A separating axis of two disjoint convex hulls is normal to
        // a hull edge (or along a degenerate point/segment axis). Every
        // hull edge is among pairs of original control points.
        [self, other].into_iter().any(|piece| {
            let chord = piece.end() - piece.controls[0];
            separates(chord)
                || separates(Vector2::new(-chord.y, chord.x))
                || piece.controls.iter().enumerate().any(|(i, &a)| {
                    piece.controls[i + 1..].iter().any(|&b| {
                        let d = b - a;
                        separates(Vector2::new(-d.y, d.x))
                    })
                })
        })
    }

    fn separated_at_join(&self, next: &Self) -> bool {
        let join = self.end();
        if !same_point(join, next.controls[0]) {
            return false;
        }
        // All non-junction controls must occupy opposite open
        // halfplanes through the same junction. Reflect the preceding
        // rays and find an open hemisphere containing both ray sets.
        let margin = self.margin + next.margin;
        let rays = self.controls[..self.controls.len() - 1]
            .iter()
            .map(|&p| join - p)
            .chain(next.controls[1..].iter().map(|&p| p - join))
            .collect::<Vec<_>>();
        let Some(&first) = rays.first() else {
            return false;
        };
        if first.magnitude() <= margin {
            return false;
        }
        let axis = first.normalize();
        let cross_axis = Vector2::new(-axis.y, axis.x);
        let [lo, hi] = rays
            .iter()
            .map(|&v| v.dot(cross_axis).atan2(v.dot(axis)))
            .fold([0.0_f64, 0.0_f64], |[lo, hi], a| [lo.min(a), hi.max(a)]);
        if hi - lo >= std::f64::consts::PI {
            return false;
        }
        let (sin, cos) = lo.midpoint(hi).sin_cos();
        let separator = axis * cos + cross_axis * sin;
        rays.iter().all(|v| v.dot(separator) > margin)
    }
}

// Define each endpoint-aligned edge in its absolute orientation BEFORE
// subdivision. Every incident face then uses exactly the same reference
// curve (or its inverse), irrespective of the proof subdivision depth.
#[derive(Clone)]
struct AlignedSpan {
    span: Span,
    ends: [Vector2; 2],
}

impl AlignedSpan {
    fn inverse(mut self) -> Self {
        self.span.controls.reverse();
        self.span.ends.reverse();
        self.ends.reverse();
        self
    }

    fn refine(self, depth: usize, result: &mut Vec<Self>) -> Option<()> {
        if result.len() >= MAX_PIECES {
            return None;
        }
        if depth == 0 || self.span.controls.len() <= 2 {
            result.push(self);
        } else {
            let [left, right] = self.span.split()?;
            let middle = left.ends[1];
            Self {
                span: left,
                ends: [self.ends[0], middle],
            }
            .refine(depth - 1, result)?;
            Self {
                span: right,
                ends: [middle, self.ends[1]],
            }
            .refine(depth - 1, result)?;
        }
        Some(())
    }

    fn points(&self) -> Vec<Vector2> {
        let mut points = self
            .span
            .controls
            .iter()
            .copied()
            .map(point)
            .collect::<Vec<_>>();
        // Construction requires at least two controls.
        points[0] = self.ends[0];
        let last = points.len() - 1;
        points[last] = self.ends[1];
        points
    }

    // For P the relative homogeneous numerator, the angular derivative
    // has the sign of det(P,P'). The denominator derivative cancels.
    // Positive bounds on all Bernstein control/derivative pairs suffice.
    fn angular_positive(&self, apex: Vector2, half_plane: Vector2, depth: usize) -> bool {
        let margin = self.span.margin * crate::usize_to_f64(depth + 1);
        let points = self.points();
        if !points
            .iter()
            .all(|p| (*p - apex).dot(half_plane) > 2.0 * margin)
        {
            return false;
        }
        let max_weight = self
            .span
            .controls
            .iter()
            .map(|c| c[3])
            .fold(0.0_f64, f64::max);
        let controls = points
            .iter()
            .zip(&self.span.controls)
            .map(|(&p, c)| (p - apex) * (c[3] / max_weight))
            .collect::<Vec<_>>();
        controls.windows(2).all(|pair| {
            let derivative = pair[1] - pair[0];
            let length = derivative.magnitude();
            controls.iter().all(|&p| {
                let norm = p.magnitude();
                let uncertainty = 2.0_f64.mul_add(margin, 2.0_f64.mul_add(norm, length));
                let error = (64.0 * f64::EPSILON * norm).mul_add(length, margin * uncertainty);
                det(p, derivative) > error
            })
        })
    }
}

fn det(a: Vector2, b: Vector2) -> f64 { a.x.mul_add(b.y, -a.y * b.x) }

fn project(p: Point3, origin: Point3, basis: [Vector3; 2]) -> Vector2 {
    let p = p - origin;
    Vector2::new(p.dot(basis[0]), p.dot(basis[1]))
}

struct ProjectedEdge {
    spans: Vec<AlignedSpan>,
    seam: f64,
    margin: f64,
    reverse: bool,
}

impl ProjectedEdge {
    fn new(edge: &Edge, origin: Point3, basis: [Vector3; 2]) -> Option<Self> {
        let reverse = !edge.orientation();
        let absolute = if reverse {
            edge.inverse()
        } else {
            edge.clone()
        };
        let endpoints =
            [absolute.front().point(), absolute.back().point()].map(|p| project(p, origin, basis));
        let raw = projected_spans(std::iter::once((0, &absolute)), origin, basis, None, 0.0)?
            .into_iter()
            .map(|(_, span)| span)
            .collect::<Vec<_>>();
        let mut seam = 0.0_f64;
        let mut margin = 0.0_f64;
        let spans = raw
            .iter()
            .enumerate()
            .map(|(index, span)| {
                let ends = [
                    if index == 0 {
                        endpoints[0]
                    } else {
                        raw[index - 1].ends[1]
                    },
                    if index + 1 == raw.len() {
                        endpoints[1]
                    } else {
                        span.ends[1]
                    },
                ];
                let gap = (ends[0] - span.ends[0])
                    .magnitude()
                    .max((ends[1] - span.ends[1]).magnitude());
                let bound = 2.0_f64.mul_add(span.margin, gap);
                if !bound.is_finite() {
                    return None;
                }
                seam = seam.max(bound);
                margin = margin.max(span.margin);
                let mut controls = span.controls.clone();
                let align = |control: &mut [f64; 4], endpoint: Vector2| {
                    control[0] = endpoint.x * control[3];
                    control[1] = endpoint.y * control[3];
                };
                align(controls.first_mut()?, ends[0]);
                align(controls.last_mut()?, ends[1]);
                Some(AlignedSpan {
                    span: Span::new(controls, span.margin)?,
                    ends,
                })
            })
            .collect::<Option<Vec<_>>>()?;
        (!spans.is_empty()).then_some(Self {
            spans,
            seam,
            margin,
            reverse,
        })
    }

    fn refined(&self, depth: usize, limit: f64) -> Option<(Vec<AlignedSpan>, f64)> {
        let seam = self
            .margin
            .mul_add(crate::usize_to_f64(depth + 1), self.seam);
        if !seam.is_finite() || seam > limit {
            return None;
        }
        let mut result = Vec::new();
        for span in &self.spans {
            span.clone().refine(depth, &mut result)?;
        }
        if self.reverse {
            result = result.into_iter().rev().map(AlignedSpan::inverse).collect();
        }
        Some((result, seam))
    }
}

struct ProjectedLoop {
    edges: Vec<ProjectedEdge>,
}

impl ProjectedLoop {
    fn new(wire: &Wire, origin: Point3, basis: [Vector3; 2]) -> Option<Self> {
        if wire.is_empty() || !wire.is_closed() || !wire.is_simple() {
            return None;
        }
        let mut edges = Vec::new();
        let mut remaining = MAX_PIECES;
        for edge in wire {
            let projected = ProjectedEdge::new(edge, origin, basis)?;
            remaining = remaining.checked_sub(projected.spans.len())?;
            edges.push(projected);
        }
        Some(Self { edges })
    }

    fn pieces(&self, wire: usize, depth: usize, limit: f64) -> Option<(Vec<Piece>, f64)> {
        let mut refined = Vec::new();
        let mut seam = 0.0_f64;
        for edge in &self.edges {
            let (pieces, gap) = edge.refined(depth, limit)?;
            if refined.len().checked_add(pieces.len())? > MAX_PIECES {
                return None;
            }
            refined.extend(pieces);
            seam = seam.max(gap);
        }
        if refined.len() < 3 {
            return None;
        }
        let count = refined.len();
        let pieces = refined
            .iter()
            .enumerate()
            .map(|(index, span)| Piece::new(&span.span, span.ends, wire, index, count, depth))
            .collect::<Option<Vec<_>>>()?;
        Some((pieces, seam))
    }
}

/// Certify that the endpoint-aligned projection of a rational boundary
/// and two straight radial sides is a positively oriented star fan.
/// The native boundary projection differs by at most the returned seam
/// bound. This uses the same absolute edge alignment as the domain
/// certificate, so a shared projected edge cancels exactly in a union.
pub(in crate::mesh_recovery) fn certify_projected_radial_fan(
    edge: &Edge,
    apex: Point3,
    origin: Point3,
    normal: Vector3,
    seam_limit: f64,
) -> Option<f64> {
    let basis = normal_basis(normal).into();
    let projected = ProjectedEdge::new(edge, origin, basis)?;
    let apex = project(apex, origin, basis);
    let [a, b] =
        [edge.front().point(), edge.back().point()].map(|p| project(p, origin, basis) - apex);
    let margin = projected.margin;
    let error = (64.0 * f64::EPSILON * a.magnitude()).mul_add(
        b.magnitude(),
        margin * (a.magnitude() + b.magnitude() + margin),
    );
    if det(a, b) <= error {
        return None;
    }
    let direction = a.normalize() + b.normalize();
    if direction.magnitude() <= 64.0 * f64::EPSILON {
        return None;
    }
    let half_plane = direction.normalize();
    for depth in 0..=MAX_DEPTH {
        let Some((spans, seam)) = projected.refined(depth, seam_limit) else {
            continue;
        };
        if spans
            .iter()
            .all(|span| span.angular_positive(apex, half_plane, depth))
        {
            return Some(seam);
        }
    }
    None
}

fn disjoint(pieces: &[Piece]) -> bool {
    let mut order = (0..pieces.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| {
        (pieces[a].bounds[0].x - pieces[a].margin)
            .total_cmp(&(pieces[b].bounds[0].x - pieces[b].margin))
    });
    let mut remaining = MAX_PAIR_TESTS;
    let mut remaining_axes = MAX_AXIS_TESTS;
    for (i, &first) in order.iter().enumerate() {
        let a = &pieces[first];
        for &second in &order[i + 1..] {
            let b = &pieces[second];
            if b.bounds[0].x - b.margin > a.bounds[1].x + a.margin {
                break;
            }
            let Some(next) = remaining.checked_sub(1) else {
                return false;
            };
            remaining = next;
            let separated = if a.next_is(b) {
                a.separated_at_join(b)
            } else if b.next_is(a) {
                b.separated_at_join(a)
            } else {
                a.separated(b, &mut remaining_axes)
            };
            if !separated {
                return false;
            }
        }
    }
    true
}

// For a simple polygon a lexicographically extreme vertex is convex.
// Its robust local turn gives the loop orientation without a global
// shoelace sum or cancellation on thin, translated domains.
fn positive_orientation(pieces: &[Piece]) -> Option<bool> {
    let count = pieces.len();
    let (index, corner) = pieces.iter().enumerate().min_by(|(_, a), (_, b)| {
        let a = a.controls[0];
        let b = b.controls[0];
        a.x.partial_cmp(&b.x)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
    })?;
    let sign = turn(
        pieces[(index + count - 1) % count].controls[0],
        corner.controls[0],
        corner.end(),
    );
    (sign.is_finite() && sign.abs() > 0.0).then_some(sign > 0.0)
}

fn contains(loop_pieces: &[Piece], p: Vector2) -> bool {
    // Pairwise hull separation ensures p is not on another loop's chord.
    // Robust orientation resolves the directed ray crossings.
    loop_pieces.iter().fold(0_i32, |winding, piece| {
        let a = piece.controls[0];
        let b = piece.end();
        if a.y <= p.y && b.y > p.y && turn(a, b, p) > 0.0 {
            winding + 1
        } else if a.y > p.y && b.y <= p.y && turn(a, b, p) < 0.0 {
            winding - 1
        } else {
            winding
        }
    }) != 0
}

fn nesting(loops: &[Vec<Piece>]) -> Option<usize> {
    let orientations = loops
        .iter()
        .map(|pieces| positive_orientation(pieces))
        .collect::<Option<Vec<_>>>()?;
    let positive = orientations
        .iter()
        .enumerate()
        .filter_map(|(i, &positive)| positive.then_some(i))
        .collect::<Vec<_>>();
    let [outer] = positive.as_slice() else {
        return None;
    };
    for (i, pieces) in loops.iter().enumerate() {
        if i == *outer {
            continue;
        }
        let witness = pieces.first()?.controls[0];
        if !contains(&loops[*outer], witness)
            || loops
                .iter()
                .enumerate()
                .any(|(j, other)| j != i && j != *outer && contains(other, witness))
        {
            return None;
        }
    }
    Some(*outer)
}

/// Certify a simple, positively oriented projected planar trim domain,
/// with a single outer loop and disjoint negatively oriented holes.
///
/// Positive-weight rational Bézier hulls bound the full homotopy to
/// ordered chords. Distinct pieces have disjoint hulls, except adjacent
/// pieces sharing a junction where strict halfplane separation applies.
/// The chord loops' winding and nesting therefore apply to the complete
/// endpoint-aligned rational loops. Original projected curves differ by
/// at most the reported seam bound; no assertion of exact unsnapped
/// incidence or 3-D plane containment is made. Exhausted proof work,
/// ambiguous geometry, or excessive seam uncertainty returns None.
#[must_use]
pub fn certify_projected_planar_domain(
    wires: &[Wire],
    origin: Point3,
    normal: Vector3,
    seam_limit: f64,
) -> Option<PlanarDomainCertificate> {
    if !seam_limit.is_finite()
        || seam_limit <= 0.0
        || !normal.magnitude2().is_finite()
        || (normal.magnitude2() - 1.0).abs() > 1e-10
        || wires.is_empty()
        || ![origin.x, origin.y, origin.z]
            .into_iter()
            .all(f64::is_finite)
    {
        return None;
    }
    let basis = normal_basis(normal).into();
    let mut projected = Vec::new();
    let mut remaining = MAX_PIECES;
    for wire in wires {
        let loop_projection = ProjectedLoop::new(wire, origin, basis)?;
        remaining = remaining.checked_sub(
            loop_projection
                .edges
                .iter()
                .map(|edge| edge.spans.len())
                .sum(),
        )?;
        projected.push(loop_projection);
    }
    for depth in 0..=MAX_DEPTH {
        let mut loops = Vec::new();
        let mut seam = 0.0_f64;
        let mut count = 0_usize;
        for (i, wire) in projected.iter().enumerate() {
            let Some((pieces, gap)) = wire.pieces(i, depth, seam_limit) else {
                break;
            };
            count = count.checked_add(pieces.len())?;
            if count > MAX_PIECES {
                return None;
            }
            seam = seam.max(gap);
            loops.push(pieces);
        }
        if loops.len() != projected.len() {
            continue;
        }
        let Some(outer_wire) = nesting(&loops) else {
            continue;
        };
        let pieces = loops.into_iter().flatten().collect::<Vec<_>>();
        if disjoint(&pieces) {
            return Some(PlanarDomainCertificate {
                outer_wire,
                pieces: count,
                seam_displacement_upper_bound: seam,
            });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BsplineCurve, Curve, EuclideanSpace, KnotVector, NurbsCurve, builder};

    fn polygon(points: &[(f64, f64)]) -> Wire {
        let vertices = points
            .iter()
            .map(|&(x, y)| builder::vertex(Point3::new(x, y, 0.0)))
            .collect::<Vec<_>>();
        vertices
            .iter()
            .zip(vertices.iter().cycle().skip(1))
            .map(|(a, b)| builder::line(a, b))
            .collect()
    }

    fn check(wires: &[Wire]) -> Option<PlanarDomainCertificate> {
        certify_projected_planar_domain(wires, Point3::new(0.0, 0.0, 0.0), Vector3::unit_z(), 1e-6)
    }

    #[test]
    fn concave_outer_with_hole_and_reordered_wires() {
        let outer = polygon(&[
            (0., 0.),
            (5., 0.),
            (5., 5.),
            (3., 5.),
            (3., 3.),
            (2., 3.),
            (2., 5.),
            (0., 5.),
        ]);
        let hole = polygon(&[(0.5, 0.5), (0.5, 2.), (1.5, 2.), (1.5, 0.5)]);
        let result = check(&[hole, outer]).expect("simple oriented domain");
        assert_eq!(result.outer_wire, 1);
        assert_eq!(result.pieces, 12);
        assert!(result.seam_displacement_upper_bound < 1e-9);
    }

    #[test]
    fn rejects_crossings_touching_reversed_and_nested_holes() {
        let outer = polygon(&[(0., 0.), (5., 0.), (5., 5.), (0., 5.)]);
        let crossing = polygon(&[(0., 0.), (4., 4.), (0., 5.), (5., 0.)]);
        let touching = polygon(&[(0., 1.), (0.5, 2.), (1., 1.)]);
        let outside = polygon(&[(6., 1.), (6., 2.), (7., 2.), (7., 1.)]);
        let hole = polygon(&[(1., 1.), (1., 4.), (4., 4.), (4., 1.)]);
        let nested = polygon(&[(2., 2.), (2., 3.), (3., 3.), (3., 2.)]);
        let reversed = polygon(&[(1., 1.), (4., 1.), (4., 4.), (1., 4.)]);
        assert!(check(&[crossing]).is_none());
        for invalid in [touching, outside, reversed] {
            assert!(check(&[outer.clone(), invalid]).is_none());
        }
        assert!(check(&[outer, hole, nested]).is_none());
    }

    fn curved_boundary(control: Point3) -> crate::Result<Wire> {
        let points = [
            Point3::new(0., 0., 0.),
            Point3::new(4., 0., 0.),
            Point3::new(4., 3., 0.),
            Point3::new(0., 3., 0.),
        ];
        let vertices = points.map(builder::vertex);
        let curve = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(2),
            vec![
                points[0].to_vec().extend(1.),
                control.to_vec().extend(1.),
                points[1].to_vec().extend(1.),
            ],
        )));
        Ok([
            Edge::new(&vertices[0], &vertices[1], curve)?,
            builder::line(&vertices[1], &vertices[2]),
            builder::line(&vertices[2], &vertices[3]),
            builder::line(&vertices[3], &vertices[0]),
        ]
        .into_iter()
        .collect())
    }

    #[test]
    fn rational_hull_subdivision_distinguishes_valid_bow_from_crossing() -> crate::Result<()> {
        // The first control hull intersects the opposite edge, but the
        // actual quadratic stays below it and subdivision proves separation.
        let valid = curved_boundary(Point3::new(2., 4., 0.))?;
        let certificate = check(&[valid]).expect("disjoint curved domain");
        assert!(certificate.pieces > 4);
        let crossing = curved_boundary(Point3::new(2., 7., 0.))?;
        assert!(check(&[crossing]).is_none());
        Ok(())
    }

    #[test]
    fn rational_circle_hole_with_positive_weights() -> crate::Result<()> {
        let outer = polygon(&[(-3., -3.), (3., -3.), (3., 3.), (-3., 3.)]);
        let points = [
            Point3::new(1., 0., 0.),
            Point3::new(0., -1., 0.),
            Point3::new(-1., 0., 0.),
            Point3::new(0., 1., 0.),
        ];
        let vertices = points.map(builder::vertex);
        let edges = (0..4)
            .map(|i| {
                let j = (i + 1) % 4;
                let midpoint =
                    Point3::new(points[i].x + points[j].x, points[i].y + points[j].y, 0.);
                let curve = Curve::NurbsCurve(NurbsCurve::new(BsplineCurve::new_unchecked(
                    KnotVector::bezier_knot(2),
                    vec![
                        points[i].to_vec().extend(1.),
                        midpoint.to_vec().extend(1.) * std::f64::consts::FRAC_1_SQRT_2,
                        points[j].to_vec().extend(1.),
                    ],
                )));
                Edge::new(&vertices[i], &vertices[j], curve)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let hole = edges.into_iter().collect::<Wire>();
        let certificate = check(&[outer, hole]).expect("simple rational circle hole");
        assert_eq!(certificate.outer_wire, 0);
        Ok(())
    }
}
