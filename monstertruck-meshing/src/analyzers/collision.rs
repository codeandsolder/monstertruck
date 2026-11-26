#[cfg(test)]
use super::EuclideanSpace;
use super::active_sweep::ActiveSweep;
use super::{
    BoundingBox, InnerSpace, KERNEL_EPSILON_SQUARED, MetricSpace, Origin, Point3, PolygonMesh,
    StandardVertex, Vector3, hash, sweep_coordinate_key, triangle_projection_interval,
};
use array_macro::array;

/// Find collisions between two polygon meshes and extract interference lines.
///
/// # Details
///
/// The algorithm for mesh extraction is as follows.
/// All the triangles are projected toward the axis to form intervals.
/// By looking at the overlap of these intervals, we can narrow down the interfering triangles.
/// The overlap of the intervals is obtained by sorting the endpoints.
///
/// # Remarks
///
/// We see that the surfaces that are in contact are not interfering with each other.
/// Therefore, polygon meshes included in one plane are not interfering.
pub trait Collision {
    /// If `self` and `other` collide, then returns only one interference line.
    /// Otherwise, returns `None`.
    fn collide_with(&self, other: &PolygonMesh) -> Option<(Point3, Point3)>;
    /// Extract all interference lines between `self` and `other`.
    /// # Remarks
    /// The results is not arranged so that included lines make continuous maximal polyline curve.
    fn extract_interference(&self, other: &PolygonMesh) -> Vec<(Point3, Point3)>;
}

impl Collision for PolygonMesh {
    #[inline]
    fn collide_with(&self, other: &PolygonMesh) -> Option<(Point3, Point3)> {
        interference_segments(self, other).next()
    }
    #[inline]
    fn extract_interference(&self, other: &PolygonMesh) -> Vec<(Point3, Point3)> {
        interference_segments(self, other).collect()
    }
}

/// Event rank encodes the original tie order: front before back,
/// then operand number, then original triangle index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct EndPoint {
    coordinate: u64,
    rank: usize,
}

impl EndPoint {
    const fn from_seg(seg: (f64, f64), index: usize, total: usize) -> [Self; 2] {
        [
            Self {
                coordinate: sweep_coordinate_key(seg.0),
                rank: index,
            },
            Self {
                coordinate: sweep_coordinate_key(seg.1),
                rank: total + index,
            },
        ]
    }
}

#[derive(Clone, Copy)]
struct TriangleData {
    points: [Point3; 3],
    bounds: BoundingBox<Point3>,
}

impl TriangleData {
    fn new(points: [Point3; 3]) -> Self {
        Self {
            bounds: points.iter().collect(),
            points,
        }
    }
}

fn sorted_endpoints(triangles: &[Vec<TriangleData>; 2], unit: Vector3) -> Vec<EndPoint> {
    let total = triangles.iter().map(Vec::len).sum();
    let mut events = triangles
        .iter()
        .flatten()
        .enumerate()
        .filter(|(_, triangle)| {
            let [a, b, c] = triangle.points;
            !(b - a).cross(c - a).so_small()
        })
        .flat_map(|(index, triangle)| {
            EndPoint::from_seg(
                triangle_projection_interval(triangle.points, unit, 0.0),
                index,
                total,
            )
        })
        .collect::<Vec<_>>();
    events.sort_unstable();
    events
}

fn colliding_segment_pairs(
    sort_endpoints: Vec<EndPoint>,
    counts: [usize; 2],
) -> impl Iterator<Item = (usize, usize)> {
    let mut current = counts.map(ActiveSweep::new);
    let total = counts.iter().sum::<usize>();
    let mut events = sort_endpoints.into_iter();
    // Drain each front event against the unchanged opposite active set
    // before reading the next event. No per-event allocation or copying.
    let mut pending = None::<(usize, usize, usize)>;
    std::iter::from_fn(move || {
        loop {
            if let Some((operand, index, offset)) = pending.as_mut() {
                if let Some(&other) = current[1 - *operand].indices().get(*offset) {
                    *offset += 1;
                    return Some(if *operand == 0 {
                        (*index, other)
                    } else {
                        (other, *index)
                    });
                }
                pending = None;
            }
            let event = events.next()?;
            let back = event.rank >= total;
            let global = if back { event.rank - total } else { event.rank };
            let operand = usize::from(global >= counts[0]);
            let index = if operand == 0 {
                global
            } else {
                global - counts[0]
            };
            if back {
                current[operand].remove(index);
            } else if current[operand].insert(index) {
                pending = Some((operand, index, 0));
            }
        }
    })
}

fn disjoint_bdbs(first: BoundingBox<Point3>, second: BoundingBox<Point3>) -> bool {
    (0..3).any(|axis| {
        first.max()[axis] < second.min()[axis] || second.max()[axis] < first.min()[axis]
    })
}

fn collide_seg_triangle(seg: [Point3; 2], tri: [Point3; 3]) -> Option<Point3> {
    let ab = tri[1] - tri[0];
    let bc = tri[2] - tri[1];
    let ca = tri[0] - tri[2];
    let nor = ab.cross(ca);
    if nor.so_small() {
        return None;
    }
    let ap = seg[0] - tri[0];
    let aq = seg[1] - tri[0];
    let front_plane_distance = ap.dot(nor);
    let back_plane_distance = aq.dot(nor);
    if front_plane_distance * back_plane_distance > 0.0 {
        return None;
    }
    let h = seg[0]
        + front_plane_distance / (front_plane_distance - back_plane_distance) * (seg[1] - seg[0]);
    if f64::signum(ab.cross(nor).dot(h - tri[0]) + KERNEL_EPSILON_SQUARED)
        + f64::signum(bc.cross(nor).dot(h - tri[1]) + KERNEL_EPSILON_SQUARED)
        + f64::signum(ca.cross(nor).dot(h - tri[2]) + KERNEL_EPSILON_SQUARED)
        >= 2.0
    {
        Some(h)
    } else {
        None
    }
}

fn collide_triangles(
    first_triangle: [Point3; 3],
    second_triangle: [Point3; 3],
) -> Option<(Point3, Point3)> {
    let mut tuple = (None, None);
    for point in [
        collide_seg_triangle([first_triangle[0], first_triangle[1]], second_triangle),
        collide_seg_triangle([first_triangle[1], first_triangle[2]], second_triangle),
        collide_seg_triangle([first_triangle[2], first_triangle[0]], second_triangle),
        collide_seg_triangle([second_triangle[0], second_triangle[1]], first_triangle),
        collide_seg_triangle([second_triangle[1], second_triangle[2]], first_triangle),
        collide_seg_triangle([second_triangle[2], second_triangle[0]], first_triangle),
    ] {
        match tuple {
            (None, _) => tuple.0 = point,
            (Some(_), None) => tuple.1 = point,
            (Some(ref mut first), Some(ref mut second)) => {
                if let Some(point) = point {
                    let first_distance = point.distance2(*first);
                    let second_distance = point.distance2(*second);
                    let endpoint_distance = first.distance2(*second);
                    if endpoint_distance < first_distance {
                        *second = point;
                    } else if endpoint_distance < second_distance {
                        *first = point;
                    }
                }
            }
        }
    }
    match tuple {
        (Some(first), Some(second)) => Some((first, second)),
        _ => None,
    }
}

const fn make_pos_tri(poly: &PolygonMesh, face: [StandardVertex; 3]) -> [Point3; 3] {
    array![i => poly.positions()[face[i].pos]; 3]
}

fn interference_segments(
    poly0: &PolygonMesh,
    poly1: &PolygonMesh,
) -> impl Iterator<Item = (Point3, Point3)> {
    let triangles = [poly0, poly1].map(|mesh| {
        if poly0.positions().is_empty() || poly1.positions().is_empty() {
            Vec::new()
        } else {
            mesh.faces()
                .triangle_iter()
                .map(|face| TriangleData::new(make_pos_tri(mesh, face)))
                .collect::<Vec<_>>()
        }
    });
    let counts = triangles.each_ref().map(Vec::len);
    let events = poly0.positions().first().map_or_else(Vec::new, |&first| {
        sorted_endpoints(&triangles, hash::take_one_unit(first))
    });
    colliding_segment_pairs(events, counts).filter_map(move |(left, right)| {
        let first = &triangles[0][left];
        let second = &triangles[1][right];
        if disjoint_bdbs(first.bounds, second.bounds) {
            None
        } else {
            collide_triangles(first.points, second.points)
        }
    })
}

#[cfg(test)]
fn reference_pairs(mut original: Vec<(f64, bool, usize, usize)>) -> Vec<(usize, usize)> {
    original.sort_by(|a, b| {
        super::compare_sweep_coordinate(a.0, b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| a.3.cmp(&b.3))
    });
    let mut active = [Vec::new(), Vec::new()];
    let mut reference = Vec::new();
    for (_, back, operand, index) in original {
        if back {
            let slot = active[operand]
                .iter()
                .position(|&value| value == index)
                .expect("active test interval");
            active[operand].swap_remove(slot);
        } else {
            active[operand].push(index);
            reference.extend(active[1 - operand].iter().map(|&other| {
                if operand == 0 {
                    (index, other)
                } else {
                    (other, index)
                }
            }));
        }
    }
    reference
}

#[test]
fn streaming_pairs_preserve_exact_event_and_active_order() {
    let counts = [73, 91];
    let total = counts.iter().sum();
    let mut events = Vec::new();
    let mut original = Vec::new();
    for (operand, count) in counts.into_iter().enumerate() {
        for index in 0..count {
            // Repeated, nested, zero-width and coincident intervals.
            let low = f64::from(
                u32::try_from((index * 17 + operand * 5) % 19).expect("small test index"),
            ) - 9.0;
            let high = low + f64::from(u32::try_from((index * 11) % 13).expect("small test index"));
            let global = if operand == 0 {
                index
            } else {
                counts[0] + index
            };
            events.extend(EndPoint::from_seg((low, high), global, total));
            original.extend([(low, false, operand, index), (high, true, operand, index)]);
        }
    }
    let reference = reference_pairs(original);
    events.sort_unstable();
    assert_eq!(
        colliding_segment_pairs(events, counts).collect::<Vec<_>>(),
        reference
    );
}

#[test]
fn cached_collision_matches_original_order_and_narrow_phase() -> anyhow::Result<()> {
    use super::{Faces, StandardAttributes};
    let mesh = |vertical: bool, shift: f64| -> anyhow::Result<PolygonMesh> {
        let mut points = Vec::new();
        let mut faces = Vec::new();
        for row in 0..12 {
            for col in 0..9 {
                let x = f64::from(col) * 0.2 - 0.6;
                let y = f64::from(row) * 0.2 - 0.9;
                let z = shift + f64::from((row + col) % 3) * 0.07;
                let triangle = if vertical {
                    [
                        Point3::new(x, z, y),
                        Point3::new(x + 0.29, z, y),
                        Point3::new(x, z, y + 0.29),
                    ]
                } else {
                    [
                        Point3::new(x, y, z),
                        Point3::new(x + 0.29, y, z),
                        Point3::new(x, y + 0.29, z),
                    ]
                };
                let start = points.len();
                points.extend(triangle);
                faces.push([start, start + 1, start + 2]);
                if (row + col) % 7 == 0 {
                    faces.push([start, start + 1, start + 1]);
                }
            }
        }
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions: points,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?)
    };
    for shift in [-0.7, 0.0, 0.2, 4.0] {
        let horizontal = mesh(false, 0.0)?;
        let vertical = mesh(true, shift)?;
        for (left, right) in [(&horizontal, &vertical), (&vertical, &horizontal)] {
            let triangles = [left, right].map(|mesh| {
                mesh.faces()
                    .triangle_iter()
                    .map(|face| make_pos_tri(mesh, face))
                    .collect::<Vec<_>>()
            });
            let unit = hash::take_one_unit(left.positions()[0]);
            let mut events = Vec::new();
            for (operand, list) in triangles.iter().enumerate() {
                for (index, &triangle) in list.iter().enumerate() {
                    if !(triangle[1] - triangle[0])
                        .cross(triangle[2] - triangle[0])
                        .so_small()
                    {
                        let (low, high) = triangle_projection_interval(triangle, unit, 0.0);
                        events.extend([(low, false, operand, index), (high, true, operand, index)]);
                    }
                }
            }
            let expected = reference_pairs(events)
                .into_iter()
                .filter_map(|(i, j)| {
                    let first = triangles[0][i];
                    let second = triangles[1][j];
                    (!disjoint_bdbs(first.iter().collect(), second.iter().collect()))
                        .then(|| collide_triangles(first, second))
                        .flatten()
                })
                .collect::<Vec<_>>();
            assert_eq!(left.extract_interference(right), expected, "shift={shift}");
            assert_eq!(left.collide_with(right), expected.first().copied());
            if shift == 0.0 {
                assert!(!expected.is_empty());
            }
        }
    }
    Ok(())
}

#[test]
fn collide_triangles_test() {
    let first_triangle = [
        Point3::origin(),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let second_triangle = [
        Point3::new(0.0, 0.0, -1.0),
        Point3::new(-1.0, -1.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
    ];
    assert!(collide_triangles(first_triangle, second_triangle).is_some());
    let first_triangle = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let second_triangle = [
        Point3::new(0.0, 0.0, 0.5),
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
    ];
    assert!(collide_triangles(first_triangle, second_triangle).is_none());
}
