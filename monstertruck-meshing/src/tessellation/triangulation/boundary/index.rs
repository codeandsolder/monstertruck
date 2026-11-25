//! Immutable UV segment hierarchy shared by winding and sampling-clearance
//! queries. Broad-phase rejection never changes the narrow-phase predicate.

use super::{InnerSpace, Itertools, Point2, SurfacePoint};

#[derive(Clone, Copy, Debug)]
struct Bounds {
    min: Point2,
    max: Point2,
}

impl Bounds {
    fn segment(front: Point2, back: Point2) -> Self {
        Self {
            min: Point2::new(front.x.min(back.x), front.y.min(back.y)),
            max: Point2::new(front.x.max(back.x), front.y.max(back.y)),
        }
    }

    fn union(self, other: Self) -> Self {
        Self {
            min: Point2::new(self.min.x.min(other.min.x), self.min.y.min(other.min.y)),
            max: Point2::new(self.max.x.max(other.max.x), self.max.y.max(other.max.y)),
        }
    }

    fn center(self, axis: usize) -> f64 { self.min[axis].midpoint(self.max[axis]) }

    fn crossing(self, point: Point2) -> bool {
        // A horizontal ray to +U cannot cross a segment strictly left of
        // the query. Keep equality for exact on-edge rejection.
        self.min.y <= point.y && point.y <= self.max.y && point.x <= self.max.x
    }

    fn near(self, point: Point2, slack: Point2) -> bool {
        let intervals: [(f64, f64); 2] = std::array::from_fn(|axis| {
            (
                (self.min[axis] - point[axis]) / slack[axis],
                (self.max[axis] - point[axis]) / slack[axis],
            )
        });
        // If the narrow-phase squared norm or dot product could overflow,
        // it may deliberately reject a nonfinite computation. Do not
        // discard such segments using a finite broad-phase comparison.
        let safe_magnitude = f64::MAX.sqrt() / 4.0;
        if intervals.iter().any(|&(lower, upper)| {
            !lower.is_finite()
                || !upper.is_finite()
                || lower.abs().max(upper.abs()) > safe_magnitude
        }) {
            return true;
        }
        intervals.into_iter().all(|(lower, upper)| {
            // The narrow phase interpolates two normalized endpoints
            // with a fraction clamped to [0,1]. Account conservatively
            // for rounding of subtraction, division, segment difference,
            // multiplication, and addition, even for a long segment.
            let coordinate_scale = lower.abs().max(upper.abs()).max(1.0);
            let roundoff = 8.0 * f64::EPSILON * coordinate_scale;
            lower <= 1.0 + roundoff && upper >= -1.0 - roundoff
        })
    }
}

#[derive(Clone, Copy, Debug)]
struct Segment {
    front: Point2,
    back: Point2,
    bounds: Bounds,
}

impl Segment {
    fn new(front: Point2, back: Point2) -> Self {
        Self {
            front,
            back,
            bounds: Bounds::segment(front, back),
        }
    }

    /// The exact same half-open crossing convention as the linear scan.
    /// None means the query lies on the boundary.
    fn winding(self, point: Point2) -> Option<i64> {
        let coord = |p: Point2| robust::Coord { x: p.x, y: p.y };
        let determinant = robust::orient2d(coord(self.front), coord(self.back), coord(point));
        if determinant == 0.0 && point.x >= self.bounds.min.x && point.x <= self.bounds.max.x {
            return None;
        }
        Some(
            if self.front.y <= point.y && point.y < self.back.y && determinant > 0.0 {
                1
            } else if self.back.y <= point.y && point.y < self.front.y && determinant < 0.0 {
                -1
            } else {
                0
            },
        )
    }

    fn clear(self, point: Point2, slack: Point2) -> bool {
        let scaled = |uv: Point2| {
            Point2::new((uv.x - point.x) / slack.x, (uv.y - point.y) / slack.y)
                - Point2::new(0.0, 0.0)
        };
        let start = scaled(self.front);
        let end = scaled(self.back);
        let segment = end - start;
        let squared = segment.magnitude2();
        let fraction = if squared > 0.0 {
            (-start.dot(segment) / squared).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (start + segment * fraction).magnitude2() > 1.0
    }
}

#[derive(Clone, Debug)]
enum Children {
    Leaf(std::ops::Range<usize>),
    Branch([usize; 2]),
}

#[derive(Clone, Debug)]
struct Node {
    bounds: Bounds,
    children: Children,
}

#[derive(Clone, Debug)]
pub(super) struct BoundaryIndex {
    segments: Vec<Segment>,
    nodes: Vec<Node>,
}

// Small leaves amortize hierarchy traversal while retaining contiguous
// segment storage. This is a performance parameter, not a geometric tolerance.
const LEAF_SEGMENTS: usize = 8;

fn build(segments: &mut [Segment], offset: usize, nodes: &mut Vec<Node>) -> usize {
    let bounds = segments[1..]
        .iter()
        .fold(segments[0].bounds, |bounds, segment| {
            bounds.union(segment.bounds)
        });
    let children = if segments.len() <= LEAF_SEGMENTS {
        Children::Leaf(offset..offset + segments.len())
    } else {
        let axis = usize::from(bounds.max.y - bounds.min.y > bounds.max.x - bounds.min.x);
        let middle = segments.len() / 2;
        segments.select_nth_unstable_by(middle, |a, b| {
            a.bounds
                .center(axis)
                .total_cmp(&b.bounds.center(axis))
                .then_with(|| {
                    a.bounds
                        .center(1 - axis)
                        .total_cmp(&b.bounds.center(1 - axis))
                })
        });
        let (left, right) = segments.split_at_mut(middle);
        Children::Branch([
            build(left, offset, nodes),
            build(right, offset + middle, nodes),
        ])
    };
    let index = nodes.len();
    nodes.push(Node { bounds, children });
    index
}

impl BoundaryIndex {
    pub(super) fn new(loops: &[Vec<SurfacePoint>]) -> Option<Self> {
        if loops
            .iter()
            .flatten()
            .any(|point| !point.x.is_finite() || !point.y.is_finite())
        {
            return None;
        }
        let mut segments = loops
            .iter()
            .flat_map(|wire| {
                wire.iter()
                    .circular_tuple_windows()
                    .map(|(front, back)| Segment::new(front.uv, back.uv))
            })
            .collect::<Vec<_>>();
        let mut nodes = Vec::new();
        if !segments.is_empty() {
            build(&mut segments, 0, &mut nodes);
        }
        Some(Self { segments, nodes })
    }

    fn visit(
        &self,
        index: usize,
        overlaps: &impl Fn(Bounds) -> bool,
        predicate: &mut impl FnMut(Segment) -> bool,
    ) -> bool {
        let node = &self.nodes[index];
        if !overlaps(node.bounds) {
            return true;
        }
        match &node.children {
            Children::Leaf(range) => self.segments[range.clone()]
                .iter()
                .copied()
                .filter(|segment| overlaps(segment.bounds))
                .all(predicate),
            Children::Branch(children) => children
                .iter()
                .all(|&child| self.visit(child, overlaps, predicate)),
        }
    }

    fn all_matching(
        &self,
        overlaps: impl Fn(Bounds) -> bool,
        mut predicate: impl FnMut(Segment) -> bool,
    ) -> bool {
        self.nodes
            .len()
            .checked_sub(1)
            .is_none_or(|root| self.visit(root, &overlaps, &mut predicate))
    }

    pub(super) fn contains(&self, point: Point2) -> bool {
        let mut winding = 0_i64;
        let interior = self.all_matching(
            |bounds| bounds.crossing(point),
            |segment| {
                if let Some(crossings) = segment.winding(point) {
                    winding += crossings;
                    true
                } else {
                    false
                }
            },
        );
        interior && winding > 0
    }

    pub(super) fn clear(&self, point: Point2, slack: Point2) -> bool {
        self.all_matching(
            |bounds| bounds.near(point, slack),
            |segment| segment.clear(point, slack),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::super::Point3;
    use super::*;
    use std::f64::consts::TAU;

    fn wire(points: impl IntoIterator<Item = Point2>) -> Vec<SurfacePoint> {
        points
            .into_iter()
            .map(|uv| SurfacePoint::from((uv, Point3::new(uv.x, uv.y, 0.0))))
            .collect()
    }

    fn ring(count: i32, center: Point2, radii: Point2, reverse: bool) -> Vec<SurfacePoint> {
        wire((0..count).map(|i| {
            let angle = TAU * f64::from(i) / f64::from(count) * if reverse { -1.0 } else { 1.0 };
            let (sin, cos) = angle.sin_cos();
            Point2::new(
                radii.x.mul_add(cos, center.x),
                radii.y.mul_add(sin, center.y),
            )
        }))
    }

    fn linear_contains(index: &BoundaryIndex, point: Point2) -> bool {
        let mut winding = 0_i64;
        for &segment in &index.segments {
            if point.y < segment.bounds.min.y || point.y > segment.bounds.max.y {
                continue;
            }
            let Some(crossings) = segment.winding(point) else {
                return false;
            };
            winding += crossings;
        }
        winding > 0
    }

    fn linear_clear(index: &BoundaryIndex, point: Point2, slack: Point2) -> bool {
        index
            .segments
            .iter()
            .all(|&segment| segment.clear(point, slack))
    }

    #[test]
    fn hierarchy_preserves_winding_and_anisotropic_clearance() {
        for (center, radii) in [
            (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0)),
            (Point2::new(1e8, -1e-7), Point2::new(3.0, 1e-9)),
            (Point2::new(-2.0, 3.0), Point2::new(1e-9, 1e5)),
        ] {
            let loops = vec![
                ring(153, center, radii, false),
                ring(67, center, Point2::new(radii.x * 0.4, radii.y * 0.4), true),
            ];
            let index = BoundaryIndex::new(&loops).expect("finite loops");
            let slack = Point2::new(radii.x * 1e-6, radii.y * 1e-6);
            for row in -35..=35 {
                for col in -35..=35 {
                    let query = Point2::new(
                        radii.x.mul_add(f64::from(col) / 30.0, center.x),
                        radii.y.mul_add(f64::from(row) / 30.0, center.y),
                    );
                    assert_eq!(
                        index.contains(query),
                        linear_contains(&index, query),
                        "{query:?}"
                    );
                    assert_eq!(
                        index.clear(query, slack),
                        linear_clear(&index, query, slack),
                        "{query:?}"
                    );
                }
            }
            for point in loops.iter().flatten() {
                assert!(!index.contains(point.uv));
                assert!(!index.clear(point.uv, slack));
                for factor in [-1.000_001, -1.0, -0.999_999, 0.999_999, 1.0, 1.000_001] {
                    let query = Point2::new(slack.x.mul_add(factor, point.x), point.y);
                    assert_eq!(
                        index.clear(query, slack),
                        linear_clear(&index, query, slack),
                        "{query:?}"
                    );
                }
            }
        }
    }

    #[test]
    #[ignore = "release-only winding and clearance timing"]
    fn benchmark_repeated_uv_queries() {
        use std::{hint::black_box, time::Instant};
        let center = Point2::new(0.0, 0.0);
        let loops = vec![
            ring(8192, center, Point2::new(1.0, 1.0), false),
            ring(4096, center, Point2::new(0.4, 0.4), true),
        ];
        let start = Instant::now();
        let index = BoundaryIndex::new(&loops).expect("finite loops");
        let build = start.elapsed();
        let queries = (-70..=70)
            .flat_map(|i| {
                (-70..=70).map(move |j| Point2::new(f64::from(i) / 60.0, f64::from(j) / 60.0))
            })
            .collect::<Vec<_>>();
        let slack = Point2::new(1e-6, 1e-6);
        let start = Instant::now();
        let reference = queries
            .iter()
            .map(|&p| {
                let p = black_box(p);
                linear_contains(&index, p) && linear_clear(&index, p, slack)
            })
            .collect::<Vec<_>>();
        let linear = start.elapsed();
        let start = Instant::now();
        let actual = queries
            .iter()
            .map(|&p| {
                let p = black_box(p);
                index.contains(p) && index.clear(p, slack)
            })
            .collect::<Vec<_>>();
        let indexed = start.elapsed();
        assert_eq!(actual, reference);
        eprintln!(
            "segments={} queries={} linear={linear:?} index_build={build:?} indexed={indexed:?} speedup={:.2}x",
            index.segments.len(),
            queries.len(),
            linear.as_secs_f64() / (build + indexed).as_secs_f64()
        );
    }
}
