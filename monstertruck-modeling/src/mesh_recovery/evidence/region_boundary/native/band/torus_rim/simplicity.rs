//! Continuous simplicity of winding-one, piecewise tensor-circle chart paths.
//! Each lifted chord is monotone in both individual angular coordinates,
//! but the complete rim need not admit a common globally monotone shear.
use super::{PathKind, RimPath};
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy)]
struct Span {
    path: RimPath,
    origin: f64,
    t: [f64; 2],
}

// The angular coordinate of a symmetric quadratic rational circle arc.
// Its homogeneous middle control is the unit vector at half the sweep.
fn angle(delta: f64, t: f64) -> f64 {
    let (sin_half, cos_half) = (delta * 0.5).sin_cos();
    let (sin_full, cos_full) = delta.sin_cos();
    let complement = 1.0 - t;
    let middle = 2.0 * t * complement;
    (middle.mul_add(sin_half, t * t * sin_full))
        .atan2(complement.mul_add(complement, middle.mul_add(cos_half, t * t * cos_full)))
}

// Exact angular rate for the same positive-weight quadratic.
fn rate(delta: f64, t: f64) -> f64 {
    let w = (2.0 * t * (1.0 - t)).mul_add((delta * 0.5).cos() - 1.0, 1.0);
    2.0 * (delta * 0.5).sin() / w
}

impl Span {
    fn delta(self) -> [f64; 2] { [self.path.advance, self.path.height[1] - self.path.height[0]] }
    fn point(self, t: f64) -> [f64; 2] {
        let d = self.delta();
        let offset = match self.path.kind {
            PathKind::Circle => d.map(|x| x * t),
            PathKind::Lifted => d.map(|x| angle(x, t)),
        };
        [self.origin + offset[0], self.path.height[0] + offset[1]]
    }
    fn bounds(self) -> [[f64; 2]; 2] {
        let a = self.point(self.t[0]);
        let b = self.point(self.t[1]);
        std::array::from_fn(|i| [a[i].min(b[i]), a[i].max(b[i])])
    }
    fn rates(self) -> [[f64; 2]; 2] {
        self.delta().map(|d| {
            if matches!(self.path.kind, PathKind::Circle) {
                return [d, d];
            }
            let a = rate(d, self.t[0]);
            let b = rate(d, self.t[1]);
            let m = rate(d, 0.5_f64.clamp(self.t[0], self.t[1]));
            [a.min(b).min(m), a.max(b).max(m)]
        })
    }
    const fn split(self) -> [Self; 2] {
        let mid = self.t[0].midpoint(self.t[1]);
        [
            Self {
                t: [self.t[0], mid],
                ..self
            },
            Self {
                t: [mid, self.t[1]],
                ..self
            },
        ]
    }
}

#[derive(Clone, Copy)]
enum Junction {
    EndStart,
    StartEnd,
}
fn separated(a: [[f64; 2]; 2], b: [[f64; 2]; 2], margin: f64) -> bool {
    (0..2).any(|i| a[i][1] + margin < b[i][0] || b[i][1] + margin < a[i][0])
}
fn extent(bounds: [[f64; 2]; 2]) -> f64 {
    (bounds[0][1] - bounds[0][0]).max(bounds[1][1] - bounds[1][0])
}
fn dot_min(rates: [[f64; 2]; 2], direction: [f64; 2]) -> f64 {
    (0..2)
        .map(|i| (rates[i][0] * direction[i]).min(rates[i][1] * direction[i]))
        .sum()
}

// A line through the shared endpoint separates the two local interiors
// when both oriented tangents progress strictly through it. Derivative
// bounds prove this over complete subspans, not just at the endpoint.
fn separated_junction(a: Span, b: Span, junction: Junction, margin: f64) -> bool {
    let (a, b) = match junction {
        Junction::EndStart => (a, b),
        Junction::StartEnd => (b, a),
    };
    let ar = a.rates();
    let br = b.rates();
    let normalize = |v: [f64; 2]| -> Option<[f64; 2]> {
        let len = v[0].hypot(v[1]);
        (len > margin).then(|| v.map(|x| x / len))
    };
    let (Some(u), Some(v)) = (
        normalize(ar.map(|r| r[0].midpoint(r[1]))),
        normalize(br.map(|r| r[0].midpoint(r[1]))),
    ) else {
        return false;
    };
    let direction = [u[0] + v[0], u[1] + v[1]];
    // The local lift accumulates rounding error. Verify that the shared
    // coordinate itself agrees to the arithmetic margin, before relying
    // on its known common topological vertex.
    let end = a.point(a.t[1]);
    let start = b.point(b.t[0]);
    end.iter().zip(start).all(|(&x, y)| (x - y).abs() <= margin)
        && dot_min(ar, direction) > margin
        && dot_min(br, direction) > margin
}

fn disjoint(
    a: Span,
    b: Span,
    junction: Option<Junction>,
    margin: f64,
    depth: usize,
    remaining: &mut usize,
) -> bool {
    if *remaining == 0 {
        return false;
    }
    *remaining -= 1;
    let ab = a.bounds();
    let bb = b.bounds();
    if separated(ab, bb, margin) || junction.is_some_and(|j| separated_junction(a, b, j, margin)) {
        return true;
    }
    if depth == 0 {
        return false;
    }
    if extent(ab) >= extent(bb) {
        let [first, last] = a.split();
        disjoint(
            first,
            b,
            junction.filter(|j| matches!(j, Junction::StartEnd)),
            margin,
            depth - 1,
            remaining,
        ) && disjoint(
            last,
            b,
            junction.filter(|j| matches!(j, Junction::EndStart)),
            margin,
            depth - 1,
            remaining,
        )
    } else {
        let [first, last] = b.split();
        disjoint(
            a,
            first,
            junction.filter(|j| matches!(j, Junction::EndStart)),
            margin,
            depth - 1,
            remaining,
        ) && disjoint(
            a,
            last,
            junction.filter(|j| matches!(j, Junction::StartEnd)),
            margin,
            depth - 1,
            remaining,
        )
    }
}

fn pair_disjoint(a: Span, b: Span, i: usize, j: usize, count: usize, margin: f64) -> bool {
    let ab = a.bounds();
    let bb = b.bounds();
    let shift = ((ab[0][0].midpoint(ab[0][1]) - bb[0][0].midpoint(bb[0][1])) / TAU).round();
    [shift - 1.0, shift, shift + 1.0].into_iter().all(|k| {
        if i == j && k.abs() < 0.5 {
            return true;
        }
        let other = Span {
            origin: k.mul_add(TAU, b.origin),
            ..b
        };
        let junction = if j == i + 1 && k.abs() < 0.5 {
            Some(Junction::EndStart)
        } else if i == 0 && j + 1 == count && (k + 1.0).abs() < 0.5 {
            Some(Junction::StartEnd)
        } else {
            None
        };
        disjoint(a, other, junction, margin, 24, &mut 4096)
    })
}

// A winding-one simple cycle on a cylinder is a boundary of an annulus
// with the disjoint reference ring. Check every potentially overlapping
// period image, including a path against its own nonzero period images.
pub(super) fn simple_periodic_cycle(paths: &[RimPath], resolution: f64) -> bool {
    let mut origin = 0.0;
    let spans = paths
        .iter()
        .map(|&path| {
            let span = Span {
                path,
                origin,
                t: [0.0, 1.0],
            };
            origin += path.advance;
            span
        })
        .collect::<Vec<_>>();
    if spans.len() < 2
        || spans.iter().any(|s| {
            let [a, b] = s.delta();
            !a.is_finite()
                || !b.is_finite()
                || a.abs().max(b.abs()) <= resolution
                || a.abs() > PI + resolution
                || (matches!(s.path.kind, PathKind::Lifted)
                    && a.abs().max(b.abs()) >= PI - resolution)
        })
    {
        return false;
    }
    if (origin - TAU).abs() > resolution
        || paths
            .iter()
            .zip(paths.iter().cycle().skip(1))
            .any(|(a, b)| (a.height[1] - b.height[0]).abs() > resolution)
    {
        return false;
    }
    let margin = resolution.max(256.0 * f64::EPSILON * crate::usize_to_f64(spans.len()) * TAU);
    spans.iter().enumerate().all(|(i, &a)| {
        spans
            .iter()
            .enumerate()
            .skip(i)
            .all(|(j, &b)| pair_disjoint(a, b, i, j, spans.len(), margin))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn paths(points: &[[f64; 2]]) -> Vec<RimPath> {
        points
            .windows(2)
            .map(|p| RimPath {
                advance: p[1][0] - p[0][0],
                height: [p[0][1], p[1][1]],
                kind: PathKind::Lifted,
            })
            .collect()
    }

    #[test]
    fn reentrant_cycle_without_a_common_shear_is_simple() {
        let rim = paths(&[
            [0.0, 1.0],
            [1.0, 1.0],
            [0.8, 1.4],
            [2.0, 1.4],
            [1.8, 1.0],
            [3.0, 1.0],
            [4.0, 1.0],
            [5.0, 1.0],
            [TAU, 1.0],
        ]);
        assert!(!super::super::shear_feasible(&rim, 1.0, 1e-8));
        assert!(!super::super::shear_feasible(&rim, -1.0, 1e-8));
        assert!(simple_periodic_cycle(&rim, 1e-8));
        assert!(super::super::simple_cycle(&rim, 1e-8));
        let crossing = paths(&[
            [0.0, 1.0],
            [1.0, 1.0],
            [0.8, 1.4],
            [2.0, 1.4],
            [0.4, 0.8],
            [3.0, 1.0],
            [4.0, 1.0],
            [5.0, 1.0],
            [TAU, 1.0],
        ]);
        assert!(!simple_periodic_cycle(&crossing, 1e-8));
    }
}
