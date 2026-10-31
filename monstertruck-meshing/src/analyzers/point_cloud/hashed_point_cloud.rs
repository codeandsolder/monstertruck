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

    fn is_colliding(&self, space: &HashedPointCloud, tolerance: f64) -> bool {
        nonpositive_tolerance!(tolerance, 0.0);
        self.distance2(space) < tolerance * tolerance
    }
}

impl DistanceWithPointCloud for Point3 {
    fn distance2(&self, space: &HashedPointCloud) -> f64 {
        let idcs = self.hash(space);
        let closure = |dist2: f64, pt: &Self| f64::min(dist2, MetricSpace::distance2(*self, *pt));
        let mut dist2 = space[idcs].iter().fold(f64::INFINITY, closure);
        if idcs[0] > 0 {
            dist2 = space[[idcs[0] - 1, idcs[1], idcs[2]]]
                .iter()
                .fold(dist2, closure);
        }
        if idcs[0] + 1 < space.size[0] {
            dist2 = space[[idcs[0] + 1, idcs[1], idcs[2]]]
                .iter()
                .fold(dist2, closure);
        }
        if idcs[1] > 0 {
            dist2 = space[[idcs[0], idcs[1] - 1, idcs[2]]]
                .iter()
                .fold(dist2, closure);
        }
        if idcs[1] + 1 < space.size[1] {
            dist2 = space[[idcs[0], idcs[1] + 1, idcs[2]]]
                .iter()
                .fold(dist2, closure);
        }
        if idcs[2] > 0 {
            dist2 = space[[idcs[0], idcs[1], idcs[2] - 1]]
                .iter()
                .fold(dist2, closure);
        }
        if idcs[2] + 1 < space.size[2] {
            dist2 = space[[idcs[0], idcs[1], idcs[2] + 1]]
                .iter()
                .fold(dist2, closure);
        }
        dist2
    }
}

impl DistanceWithPointCloud for [Point3; 3] {
    fn distance2(&self, space: &HashedPointCloud) -> f64 {
        let bdd: BoundingBox<Point3> = self.iter().collect();
        let mut range: [[usize; 3]; 2] = [bdd.min().hash(space), bdd.max().hash(space)];
        range[0][0] = usize::max(range[0][0], 1) - 1;
        range[0][1] = usize::max(range[0][1], 1) - 1;
        range[0][2] = usize::max(range[0][2], 1) - 1;
        range[1][0] = usize::min(range[1][0] + 1, space.size[0] - 1);
        range[1][1] = usize::min(range[1][1] + 1, space.size[1] - 1);
        range[1][2] = usize::min(range[1][2] + 1, space.size[2] - 1);
        let mut dist2 = f64::INFINITY;
        (range[0][0]..=range[1][0]).for_each(|ix| {
            (range[0][1]..=range[1][1]).for_each(|iy| {
                (range[0][2]..=range[1][2]).for_each(|iz| {
                    dist2 = space[[ix, iy, iz]].iter().fold(dist2, |dist2, pt| {
                        f64::min(dist2, distance2_point_triangle(*pt, *self))
                    });
                });
            });
        });
        dist2
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
            tri.distance2(space) < tol * tol
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
