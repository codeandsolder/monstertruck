//! Immutable UV segment hierarchy shared by winding and sampling-clearance
//! queries. Broad-phase rejection never changes the narrow-phase predicate.

use super::{InnerProduct, Itertools, Point2, SurfacePoint};
use crate::spatial_hierarchy::{Bounds, SpatialHierarchy};

type UvBounds = Bounds<2>;

impl UvBounds {
    fn crossing(self, point: Point2) -> bool {
        // A horizontal ray to +U cannot cross a segment strictly left of
        // the query. Keep equality for exact on-edge rejection.
        self.low[1] <= point.y && point.y <= self.high[1] && point.x <= self.high[0]
    }

    fn near(self, point: Point2, slack: Point2) -> bool {
        let intervals: [(f64, f64); 2] = std::array::from_fn(|axis| {
            (
                (self.low[axis] - point[axis]) / slack[axis],
                (self.high[axis] - point[axis]) / slack[axis],
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
    bounds: UvBounds,
}

impl Segment {
    fn new(front: Point2, back: Point2) -> Self {
        Self {
            front,
            back,
            bounds: UvBounds::from_points([[front.x, front.y], [back.x, back.y]]),
        }
    }

    /// The exact same half-open crossing convention as the linear scan.
    /// None means the query lies on the boundary.
    fn winding(self, point: Point2) -> Option<i64> {
        let coord = |p: Point2| robust::Coord { x: p.x, y: p.y };
        let determinant = robust::orient2d(coord(self.front), coord(self.back), coord(point));
        if determinant == 0.0 && point.x >= self.bounds.low[0] && point.x <= self.bounds.high[0] {
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
        let squared = segment.norm_squared();
        let fraction = if squared > 0.0 {
            (-start.dot(&segment) / squared).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (start + segment * fraction).norm_squared() > 1.0
    }
}

#[derive(Clone, Debug)]
pub(super) struct BoundaryIndex {
    segments: Vec<Segment>,
    hierarchy: SpatialHierarchy<2>,
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
        let hierarchy = SpatialHierarchy::new(&mut segments, |segment| segment.bounds);
        Some(Self {
            segments,
            hierarchy,
        })
    }

    fn all_matching(
        &self,
        overlaps: impl Fn(UvBounds) -> bool,
        mut predicate: impl FnMut(Segment) -> bool,
    ) -> bool {
        self.hierarchy.all_matching(&overlaps, |range| {
            self.segments[range]
                .iter()
                .copied()
                .filter(|segment| overlaps(segment.bounds))
                .all(&mut predicate)
        })
    }

    pub(super) fn contains(&self, point: Point2) -> bool {
        let mut winding = 0_i64;
        let interior = self.all_matching(
            |bounds| bounds.crossing(point),
            |segment| {
                segment.winding(point).is_some_and(|crossings| {
                    winding += crossings;
                    true
                })
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
            if point.y < segment.bounds.low[1] || point.y > segment.bounds.high[1] {
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
