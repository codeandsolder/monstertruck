use super::{
    BoundingBox, MetricSpace, Point3, PolygonMesh, distance2_point_triangle, nonpositive_tolerance,
};
use array_macro::array;

#[cfg(test)]
use crate::Vector3;

#[derive(Clone, Debug)]
pub(super) struct HashedPointCloud {
    space: Vec<Vec<Point3>>,
    size: [usize; 3],
    range: [[f64; 2]; 3],
}

impl HashedPointCloud {
    fn new(size: [usize; 3], range: [[f64; 2]; 3]) -> Option<Self> {
        if size.contains(&0) {
            return None;
        }
        let cell_count = size.into_iter().try_fold(1usize, usize::checked_mul)?;
        Some(Self {
            space: vec![Vec::new(); cell_count],
            size,
            range,
        })
    }

    pub(super) fn from_points<'a, I>(points: I, min_cell_size: f64) -> Option<Self>
    where
        I: IntoIterator + Clone,
        I::IntoIter: Iterator<Item = &'a Point3>, {
        if !min_cell_size.is_finite() || min_cell_size <= 0.0 {
            return None;
        }
        let mut bounds = BoundingBox::<Point3>::new();
        let mut len = 0usize;
        for point in points.clone() {
            bounds.push(*point);
            len = len.saturating_add(1);
        }
        if len == 0 {
            return None;
        }

        // Grid density is only a search heuristic. Capping the point-count input
        // at u32 keeps all subsequent integer/float conversions exact while still
        // allowing more cells than can realistically be allocated here.
        let capped_len = u32::try_from(len).unwrap_or(u32::MAX);
        let max_axis_bins = integer_cube_root(capped_len).saturating_add(1);
        let diagonal = bounds.diagonal();
        let size = [
            axis_bin_count(diagonal[0], min_cell_size, max_axis_bins)?,
            axis_bin_count(diagonal[1], min_cell_size, max_axis_bins)?,
            axis_bin_count(diagonal[2], min_cell_size, max_axis_bins)?,
        ];
        let range = array![i => [bounds.min()[i], bounds.max()[i]]; 3];
        let mut result = Self::new(size, range)?;
        for point in points {
            result.push(*point);
        }
        Some(result)
    }

    fn push(&mut self, point: Point3) {
        let index = point.hash(self);
        self[index].push(point);
    }

    // Hash every cell that can contain a point within radius of the query
    // AABB, including adjacent bins at floating-point partition boundaries.
    fn cells_near_bounds(
        &self,
        bounds: &BoundingBox<Point3>,
        radius: f64,
    ) -> Option<([usize; 3], [usize; 3])> {
        if !radius.is_finite()
            || radius < 0.0
            || (0..3).any(|axis| !bounds.min()[axis].is_finite() || !bounds.max()[axis].is_finite())
        {
            return None;
        }
        let lower = array![axis => {
            let limit = bounds.min()[axis] - radius;
            if limit.is_finite() {
                hash_axis(limit.next_down(), self.range[axis], self.size[axis])
            } else {
                0
            }
        }; 3];
        let upper = array![axis => {
            let limit = bounds.max()[axis] + radius;
            if limit.is_finite() {
                hash_axis(limit.next_up(), self.range[axis], self.size[axis])
            } else {
                self.size[axis] - 1
            }
        }; 3];
        Some((lower, upper))
    }

    fn any_within_distance2(
        &self,
        bounds: BoundingBox<Point3>,
        radius: f64,
        distance2: impl Fn(Point3) -> f64,
    ) -> bool {
        let radius2 = radius * radius;
        let Some((lower, upper)) = self.cells_near_bounds(&bounds, radius) else {
            return self
                .space
                .iter()
                .flatten()
                .any(|point| distance2(*point) < radius2);
        };
        for x in lower[0]..=upper[0] {
            for y in lower[1]..=upper[1] {
                for z in lower[2]..=upper[2] {
                    if self[[x, y, z]]
                        .iter()
                        .any(|point| distance2(*point) < radius2)
                    {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Find the actual nearest sample to a bounded query geometry.
    ///
    /// The initial search uses the grid, but arbitrary nearby buckets cannot
    /// establish a global minimum. A finite candidate gives an upper bound:
    /// any better point must lie inside the query AABB expanded by its
    /// distance. We hash those limits and visit every potentially better cell.
    fn minimum_distance2(
        &self,
        bounds: BoundingBox<Point3>,
        distance2: impl Fn(Point3) -> f64,
    ) -> f64 {
        let exhaustive = || {
            self.space
                .iter()
                .flatten()
                .fold(f64::INFINITY, |best, point| best.min(distance2(*point)))
        };
        if (0..3).any(|axis| !bounds.min()[axis].is_finite() || !bounds.max()[axis].is_finite()) {
            return exhaustive();
        }

        let low = bounds.min().hash(self).map(|index| index.saturating_sub(1));
        let high = bounds.max().hash(self);
        let high = array![axis => high[axis].saturating_add(1).min(self.size[axis] - 1); 3];
        let mut best = f64::INFINITY;
        for x in low[0]..=high[0] {
            for y in low[1]..=high[1] {
                for z in low[2]..=high[2] {
                    for point in &self[[x, y, z]] {
                        best = best.min(distance2(*point));
                    }
                }
            }
        }
        if best == 0.0 {
            return 0.0;
        }
        if !best.is_finite() {
            // A sparse grid can leave the initial neighborhood empty.
            // One complete scan is cheaper than seeding an arbitrary distant
            // witness and scanning its expanded bounding region again.
            return exhaustive();
        }

        // Expand slightly for floating-point rounding at bin boundaries.
        // The hash is monotone, so all points that can beat this upper bound
        // reside between the hashed limits on each axis.
        let radius = best.sqrt().next_up();
        let Some((lower, upper)) = self.cells_near_bounds(&bounds, radius) else {
            return exhaustive();
        };
        for x in lower[0]..=upper[0] {
            for y in lower[1]..=upper[1] {
                for z in lower[2]..=upper[2] {
                    if (low[0]..=high[0]).contains(&x)
                        && (low[1]..=high[1]).contains(&y)
                        && (low[2]..=high[2]).contains(&z)
                    {
                        continue;
                    }
                    for point in &self[[x, y, z]] {
                        best = best.min(distance2(*point));
                    }
                }
            }
        }
        best
    }

    pub(super) fn distance2<T: DistanceWithPointCloud + ?Sized>(&self, target: &T) -> f64 {
        target.distance2(self)
    }

    pub(super) fn is_colliding<T: DistanceWithPointCloud + ?Sized>(
        &self,
        target: &T,
        tolerance: f64,
    ) -> bool {
        target.is_colliding(self, tolerance)
    }
}

fn integer_cube_root(value: u32) -> u32 {
    let value = u64::from(value);
    let mut root = 0u32;
    loop {
        let Some(next) = root.checked_add(1) else {
            return root;
        };
        if u64::from(next).pow(3) > value {
            return root;
        }
        root = next;
    }
}

fn axis_bin_count(span: f64, min_cell_size: f64, max_bins: u32) -> Option<usize> {
    if !span.is_finite() || span <= 0.0 {
        return Some(1);
    }
    let bins = if span / f64::from(max_bins) < min_cell_size {
        // Preserve the old 1 + floor(span / min_cell_size) behavior without a
        // float-to-integer cast. max_bins is at most ~1626 for a u32 point count.
        let mut bins = 1u32;
        while bins < max_bins && f64::from(bins) * min_cell_size <= span {
            bins += 1;
        }
        bins
    } else {
        max_bins
    };
    usize::try_from(bins).ok()
}

fn hash_axis(value: f64, range: [f64; 2], bins: usize) -> usize {
    if bins <= 1
        || !value.is_finite()
        || !range[0].is_finite()
        || !range[1].is_finite()
        || range[1] <= range[0]
        || value <= range[0]
    {
        return 0;
    }
    if value >= range[1] {
        return bins - 1;
    }

    // Binary-search the same floor(normalized * bins) partitioning used by the
    // former float-to-usize cast. Grid sizes originate from u32 above, so these
    // conversions remain exact in f64.
    let bins_f64 = f64::from(u32::try_from(bins).unwrap_or(u32::MAX));
    let normalized = (value - range[0]) / (range[1] - range[0]);
    let mut low = 0usize;
    let mut high = bins;
    while low + 1 < high {
        let mid = low + (high - low) / 2;
        let mid_f64 = f64::from(u32::try_from(mid).unwrap_or(u32::MAX));
        if normalized < mid_f64 / bins_f64 {
            high = mid;
        } else {
            low = mid;
        }
    }
    low
}

impl std::ops::Index<[usize; 3]> for HashedPointCloud {
    type Output = Vec<Point3>;

    fn index(&self, indices: [usize; 3]) -> &Self::Output {
        &self.space[(indices[0] * self.size[1] + indices[1]) * self.size[2] + indices[2]]
    }
}

impl std::ops::IndexMut<[usize; 3]> for HashedPointCloud {
    fn index_mut(&mut self, indices: [usize; 3]) -> &mut Self::Output {
        &mut self.space[(indices[0] * self.size[1] + indices[1]) * self.size[2] + indices[2]]
    }
}

trait SpaceHash {
    fn hash(self, space: &HashedPointCloud) -> [usize; 3];
}

impl SpaceHash for Point3 {
    fn hash(self, space: &HashedPointCloud) -> [usize; 3] {
        [
            hash_axis(self[0], space.range[0], space.size[0]),
            hash_axis(self[1], space.range[1], space.size[1]),
            hash_axis(self[2], space.range[2], space.size[2]),
        ]
    }
}

impl SpaceHash for usize {
    fn hash(self, space: &HashedPointCloud) -> [usize; 3] {
        [
            self / (space.size[1] * space.size[2]),
            self % (space.size[1] * space.size[2]) / space.size[2],
            self % space.size[2],
        ]
    }
}

pub(super) trait DistanceWithPointCloud {
    fn distance2(&self, space: &HashedPointCloud) -> f64;

    fn is_colliding(&self, space: &HashedPointCloud, tolerance: f64) -> bool;
}

impl DistanceWithPointCloud for Point3 {
    fn distance2(&self, space: &HashedPointCloud) -> f64 {
        let mut bounds = BoundingBox::<Self>::new();
        bounds.push(*self);
        space.minimum_distance2(bounds, |point| MetricSpace::distance2(*self, point))
    }

    fn is_colliding(&self, space: &HashedPointCloud, tolerance: f64) -> bool {
        nonpositive_tolerance!(tolerance, 0.0);
        let mut bounds = BoundingBox::<Self>::new();
        bounds.push(*self);
        space.any_within_distance2(bounds, tolerance, |point| {
            MetricSpace::distance2(*self, point)
        })
    }
}

impl DistanceWithPointCloud for [Point3; 3] {
    fn distance2(&self, space: &HashedPointCloud) -> f64 {
        let bounds: BoundingBox<Point3> = self.iter().collect();
        space.minimum_distance2(bounds, |point| distance2_point_triangle(point, *self))
    }

    fn is_colliding(&self, space: &HashedPointCloud, tolerance: f64) -> bool {
        nonpositive_tolerance!(tolerance, 0.0);
        let bounds: BoundingBox<Point3> = self.iter().collect();
        space.any_within_distance2(bounds, tolerance, |point| {
            distance2_point_triangle(point, *self)
        })
    }
}

impl DistanceWithPointCloud for &PolygonMesh {
    fn distance2(&self, space: &HashedPointCloud) -> f64 {
        let dist2 = self.faces().triangle_iter().fold(-1.0, |dist2, tri| {
            let tri = array![i => self.positions()[tri[i].pos]; 3];
            f64::max(dist2, tri.distance2(space))
        });
        if dist2 < 0.0 { f64::INFINITY } else { dist2 }
    }
    fn is_colliding(&self, space: &HashedPointCloud, tol: f64) -> bool {
        nonpositive_tolerance!(tol, 0.0);
        self.faces().triangle_iter().any(|tri| {
            let tri = array![i => self.positions()[tri[i].pos]; 3];
            tri.is_colliding(space, tol)
        })
    }
}

#[test]
fn hashed_point_cloud_rejects_empty_or_nonpositive_cell_size() {
    let points = [Point3::new(0.0, 0.0, 0.0)];
    assert!(HashedPointCloud::from_points(&points[..0], 1.0).is_none());
    assert!(HashedPointCloud::from_points(&points, 0.0).is_none());
    assert!(HashedPointCloud::from_points(&points, -1.0).is_none());
    assert!(HashedPointCloud::from_points(&points, f64::NAN).is_none());
    assert!(HashedPointCloud::from_points(&points, f64::INFINITY).is_none());
}

#[test]
fn hashed_point_cloud_handles_zero_span_axes_without_nan_indexing() {
    let points = [
        Point3::new(1.0, 2.0, 3.0),
        Point3::new(1.0, 4.0, 3.0),
        Point3::new(1.0, 6.0, 3.0),
    ];
    let Some(cloud) = HashedPointCloud::from_points(&points, 0.5) else {
        return;
    };
    assert_eq!(points[0].hash(&cloud)[0], 0);
    assert_eq!(points[0].hash(&cloud)[2], 0);
    assert!(cloud.distance2(&points[1]).is_finite());
}

#[test]
fn nearest_point_considers_diagonal_and_nonlocal_grid_cells() {
    // A diagonal neighbor is not one of the six face-adjacent buckets.
    let diagonal = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(51.0, 51.0, 51.0),
        Point3::new(100.0, 100.0, 100.0),
    ];
    let Some(cloud) = HashedPointCloud::from_points(&diagonal, 1.0) else {
        panic!("valid point cloud");
    };
    let query = Point3::new(49.0, 49.0, 49.0);
    assert_eq!(cloud.distance2(&query), 12.0);

    // With five X bins, the middle three can be empty. A nearest-point query
    // there must produce a finite answer even if its immediate neighbors do not.
    let points = (0..64)
        .map(|i| Point3::new(if i % 2 == 0 { 0.0 } else { 100.0 }, 0.0, 0.0))
        .collect::<Vec<_>>();
    let Some(cloud) = HashedPointCloud::from_points(&points, 1.0) else {
        panic!("valid point cloud");
    };
    let query = Point3::new(50.0, 0.0, 0.0);
    assert_eq!(cloud.distance2(&query), 2500.0);

    let triangle = [
        Point3::new(50.0, 0.0, 0.0),
        Point3::new(50.0, 1.0, 0.0),
        Point3::new(50.0, 0.0, 1.0),
    ];
    assert_eq!(cloud.distance2(&triangle), 2500.0);
}

#[test]
fn bounded_nearest_distance_matches_exhaustive_geometry() {
    let points = (0..173)
        .map(|i| {
            Point3::new(
                f64::from((i * 17 % 113) - 56),
                f64::from((i * 29 % 127) - 63),
                f64::from((i * 31 % 97) - 48),
            )
        })
        .collect::<Vec<_>>();
    let Some(cloud) = HashedPointCloud::from_points(&points, 4.0) else {
        panic!("valid point cloud");
    };
    for n in 0..90 {
        let query = Point3::new(
            f64::from((n * 37 % 183) - 91),
            f64::from((n * 41 % 139) - 69),
            f64::from((n * 43 % 149) - 74),
        );
        let expected = points
            .iter()
            .map(|point| MetricSpace::distance2(*point, query))
            .fold(f64::INFINITY, f64::min);
        assert_eq!(cloud.distance2(&query), expected, "point query {n}");

        let triangle = [
            query,
            query + Vector3::new(11.0, 0.0, 0.0),
            query + Vector3::new(0.0, 7.0, 3.0),
        ];
        let expected = points
            .iter()
            .map(|point| distance2_point_triangle(*point, triangle))
            .fold(f64::INFINITY, f64::min);
        assert_eq!(cloud.distance2(&triangle), expected, "triangle query {n}");
        for tolerance in [0.0, 0.1, 2.0, 11.0, 80.0] {
            let expected_point = points
                .iter()
                .any(|point| MetricSpace::distance2(*point, query) < tolerance * tolerance);
            let expected_triangle = points
                .iter()
                .any(|point| distance2_point_triangle(*point, triangle) < tolerance * tolerance);
            assert_eq!(cloud.is_colliding(&query, tolerance), expected_point);
            assert_eq!(cloud.is_colliding(&triangle, tolerance), expected_triangle);
        }
    }
}

#[test]
fn distance_point_triangle_test() {
    let triangle = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    assert_eq!(
        12.45 * 12.45,
        distance2_point_triangle(Point3::new(0.25, 0.25, 12.45), triangle)
    );
    assert_eq!(
        12.45 * 12.45,
        distance2_point_triangle(Point3::new(0.25, 0.25, -12.45), triangle)
    );
    assert_eq!(
        2.0,
        distance2_point_triangle(Point3::new(-1.0, 0.5, -1.0), triangle)
    );
    assert_eq!(
        2.0,
        distance2_point_triangle(Point3::new(0.5, -1.0, 1.0), triangle)
    );
    assert_eq!(
        1.5,
        distance2_point_triangle(Point3::new(1.0, 1.0, 1.0), triangle)
    );
    assert_eq!(
        3.0,
        distance2_point_triangle(Point3::new(-1.0, -1.0, -1.0), triangle)
    );
}

#[cfg(test)]
fn exec_space_division_distance() {
    const SPACE_SIZE: f64 = 100.0;
    const NUM_POINTS: usize = 1000;
    const NUM_TRIANGLES: usize = 10;
    const TRIANGLE_DISPLACEMENT: f64 = 1.0;

    let points = (0..NUM_POINTS)
        .map(|_| {
            Point3::new(
                SPACE_SIZE * rand::random::<f64>() - SPACE_SIZE / 2.0,
                SPACE_SIZE * rand::random::<f64>() - SPACE_SIZE / 2.0,
                SPACE_SIZE * rand::random::<f64>() - SPACE_SIZE / 2.0,
            )
        })
        .collect::<Vec<_>>();
    let triangles = (0..NUM_TRIANGLES)
        .map(|_| {
            let pt = Point3::new(
                SPACE_SIZE * rand::random::<f64>() - SPACE_SIZE / 2.0,
                SPACE_SIZE * rand::random::<f64>() - SPACE_SIZE / 2.0,
                SPACE_SIZE * rand::random::<f64>() - SPACE_SIZE / 2.0,
            );
            [
                pt,
                pt + Vector3::new(
                    TRIANGLE_DISPLACEMENT * 2.0 * rand::random::<f64>() - TRIANGLE_DISPLACEMENT,
                    TRIANGLE_DISPLACEMENT * 2.0 * rand::random::<f64>() - TRIANGLE_DISPLACEMENT,
                    TRIANGLE_DISPLACEMENT * 2.0 * rand::random::<f64>() - TRIANGLE_DISPLACEMENT,
                ),
                pt + Vector3::new(
                    TRIANGLE_DISPLACEMENT * 2.0 * rand::random::<f64>() - TRIANGLE_DISPLACEMENT,
                    TRIANGLE_DISPLACEMENT * 2.0 * rand::random::<f64>() - TRIANGLE_DISPLACEMENT,
                    TRIANGLE_DISPLACEMENT * 2.0 * rand::random::<f64>() - TRIANGLE_DISPLACEMENT,
                ),
            ]
        })
        .collect::<Vec<_>>();
    let Some(hashed) = HashedPointCloud::from_points(&points, 1.0) else {
        return;
    };
    let dist_0 = triangles.iter().fold(f64::INFINITY, |dist, triangle| {
        f64::min(dist, hashed.distance2(triangle))
    });
    let dist_1 = points.iter().fold(f64::INFINITY, |dist2, pt| {
        triangles.iter().fold(dist2, |dist2, triangle| {
            f64::min(dist2, distance2_point_triangle(*pt, *triangle))
        })
    });
    assert_eq!(dist_0, dist_1);
}

#[test]
fn space_division_distance() { (0..10).for_each(|_| exec_space_division_distance()) }
