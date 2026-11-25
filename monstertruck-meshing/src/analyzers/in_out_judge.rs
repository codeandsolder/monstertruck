use super::{InnerSpace, Matrix3, Origin, Point3, PolygonMesh, SquareMatrix, Vector3, hash};
use crate::spatial_hierarchy::{Bounds, SpatialHierarchy};
use derive_more::{Deref, DerefMut};
use std::borrow::Borrow;

#[cfg(test)]
use crate::{Faces, StandardAttributes};

#[derive(Clone, Copy, Debug)]
struct Ray {
    origin: Point3,
    direction: Vector3,
}

#[derive(Clone, Copy, Debug, Deref, DerefMut)]
struct Triangle([Point3; 3]);

impl Triangle {
    #[inline]
    fn normal(self) -> Vector3 { (self[1] - self[0]).cross(self[2] - self[0]).normalize() }

    fn is_crossing(self, ray: Ray) -> bool {
        let a = self[0] - self[1];
        let b = self[0] - self[2];
        let mat = Matrix3::from_cols(a, b, ray.direction);
        if mat.determinant().so_small() {
            false
        } else {
            let Some(inv) = mat.invert() else {
                return false;
            };
            let uvt = inv * (self[0] - ray.origin);
            uvt[0] > 0.0 && uvt[1] > 0.0 && uvt[0] + uvt[1] < 1.0 && uvt[2] > 0.0
        }
    }
    fn signed_crossing(self, ray: Ray) -> isize {
        if !self.is_crossing(ray) {
            return 0;
        }
        let direction = self.normal().dot(ray.direction);
        if direction > 0.0 {
            1
        } else if direction < 0.0 {
            -1
        } else {
            0
        }
    }
}

// Preserve the polygon fan of the original crossing implementation, rather
// than substituting a different triangulation for nontriangular faces.
fn triangles(mesh: &PolygonMesh) -> impl Iterator<Item = Triangle> + '_ {
    mesh.face_iter().flat_map(move |face| {
        (2..face.len()).map(move |index| {
            Triangle([
                mesh.positions()[face[0].pos],
                mesh.positions()[face[index - 1].pos],
                mesh.positions()[face[index].pos],
            ])
        })
    })
}

type RayBounds = Bounds<3>;

impl RayBounds {
    // Conservatively round each subtraction and division outward. Do not
    // add a geometric tolerance or exclude tangent/parallel boundary rays:
    // the unchanged triangle intersection test makes the final decision.
    fn intersects(self, ray: Ray) -> bool {
        let mut enter = 0.0_f64;
        let mut exit = f64::INFINITY;
        for axis in 0..3 {
            let origin = ray.origin[axis];
            let direction = ray.direction[axis];
            if direction == 0.0 {
                if origin < self.low[axis] || origin > self.high[axis] {
                    return false;
                }
                continue;
            }
            let low = (self.low[axis] - origin).next_down();
            let high = (self.high[axis] - origin).next_up();
            let (near, far) = if direction > 0.0 {
                (low / direction, high / direction)
            } else {
                (high / direction, low / direction)
            };
            enter = enter.max(near.next_down());
            exit = exit.min(far.next_up());
            if enter > exit {
                return false;
            }
        }
        true
    }
}

#[derive(Clone, Copy, Debug)]
struct OrientedTriangle {
    triangle: Triangle,
    sign: isize,
    bounds: RayBounds,
}

/// Immutable spatial index for repeated signed ray queries.
///
/// Each input mesh contributes its original fan triangles and topological
/// orientation. The index changes only candidate discovery; crossing and
/// winding use the same narrow-phase test as `PolygonMesh`.
#[derive(Clone, Debug)]
pub struct SignedRayIndex {
    triangles: Vec<OrientedTriangle>,
    hierarchy: Option<SpatialHierarchy<3>>,
}

impl SignedRayIndex {
    /// Build an owned index from oriented meshes. A false orientation
    /// subtracts that mesh's original winding contribution.
    #[must_use = "retain the spatial ray index before querying signed mesh crossings"]
    pub fn new<M: Borrow<PolygonMesh>>(meshes: impl IntoIterator<Item = (M, bool)>) -> Self {
        let mut entries = Vec::new();
        let mut finite = true;
        for (mesh, forward) in meshes {
            let sign = if forward { 1 } else { -1 };
            for triangle in triangles(mesh.borrow()) {
                finite &= triangle
                    .iter()
                    .all(|p| (0..3).all(|axis| p[axis].is_finite()));
                entries.push(OrientedTriangle {
                    triangle,
                    sign,
                    bounds: RayBounds::from_points(triangle.0.map(|p| [p.x, p.y, p.z])),
                });
            }
        }
        let hierarchy = finite.then(|| SpatialHierarchy::new(&mut entries, |entry| entry.bounds));
        Self {
            triangles: entries,
            hierarchy,
        }
    }

    fn sum_triangles(entries: &[OrientedTriangle], ray: Ray) -> isize {
        entries
            .iter()
            .map(|entry| entry.sign * entry.triangle.signed_crossing(ray))
            .sum()
    }

    /// Count signed crossings of a positive ray. Nonfinite input uses the
    /// exhaustive narrow phase instead of unsafe spatial rejection.
    #[must_use = "use the signed crossing count when classifying a point against the mesh"]
    pub fn signed_crossing_faces(&self, point: Point3, ray_direction: Vector3) -> isize {
        let ray = Ray {
            origin: point,
            direction: ray_direction,
        };
        let Some(hierarchy) = self.hierarchy.as_ref().filter(|_| {
            (0..3).all(|axis| point[axis].is_finite() && ray_direction[axis].is_finite())
        }) else {
            return Self::sum_triangles(&self.triangles, ray);
        };
        let mut total = 0;
        hierarchy.all_matching(
            |bounds| bounds.intersects(ray),
            |range| {
                total += Self::sum_triangles(&self.triangles[range], ray);
                true
            },
        );
        total
    }
}

/// whether a point is in a domain rounded by a closed polygon.
pub trait IncludingPointInDomain {
    /// Count signed number of faces crossing ray with origin `point` and direction `ray_direction`.
    /// Counter increase if the dot product of the ray and the normal of a face is positive,
    /// and decrease if it is negative.
    ///
    /// # Examples
    /// ```
    /// use monstertruck_meshing::prelude::*;
    /// let simplex = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.0, 0.0),
    ///             Point3::new(1.0, 0.0, 0.0),
    ///             Point3::new(0.0, 1.0, 0.0),
    ///             Point3::new(0.0, 0.0, 1.0),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(vec![
    ///         [0, 2, 1],
    ///         [0, 1, 3],
    ///         [0, 3, 2],
    ///         [1, 2, 3],
    ///     ]),
    /// ).expect("fixture mesh is valid");
    ///
    /// assert_eq!(simplex.signed_crossing_faces(Point3::new(0.1, 0.1, 0.1), Vector3::unit_x()), 1);
    /// assert_eq!(simplex.signed_crossing_faces(Point3::new(-0.1, 0.1, 0.1), Vector3::unit_x()), 0);
    /// ```
    fn signed_crossing_faces(&self, point: Point3, ray_direction: Vector3) -> isize;
    /// whether `point` is in a domain rounded by a closed polygon.
    ///
    /// # Examples
    /// ```
    /// use monstertruck_meshing::prelude::*;
    /// let simplex = PolygonMesh::new(
    ///     StandardAttributes {
    ///         positions: vec![
    ///             Point3::new(0.0, 0.0, 0.0),
    ///             Point3::new(1.0, 0.0, 0.0),
    ///             Point3::new(0.0, 1.0, 0.0),
    ///             Point3::new(0.0, 0.0, 1.0),
    ///         ],
    ///         ..Default::default()
    ///     },
    ///     Faces::from_iter(vec![
    ///         [0, 2, 1],
    ///         [0, 1, 3],
    ///         [0, 3, 2],
    ///         [1, 2, 3],
    ///     ]),
    /// ).expect("fixture mesh is valid");
    ///
    /// assert!(simplex.inside(Point3::new(0.1, 0.1, 0.1)));
    /// assert!(!simplex.inside(Point3::new(-0.1, 0.1, 0.1)));
    /// ```
    fn inside(&self, point: Point3) -> bool;
}

impl IncludingPointInDomain for PolygonMesh {
    fn signed_crossing_faces(&self, point: Point3, ray_direction: Vector3) -> isize {
        let ray = Ray {
            origin: point,
            direction: ray_direction,
        };
        triangles(self)
            .map(|triangle| triangle.signed_crossing(ray))
            .sum()
    }
    fn inside(&self, point: Point3) -> bool {
        let dir = hash::take_one_unit(point);
        self.signed_crossing_faces(point, dir) >= 1
    }
}

#[test]
fn inside100() {
    let positions = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ];
    let faces: Faces = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
        .into_iter()
        .collect();
    let simplex = PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        faces,
    )
    .expect("fixture mesh is valid");

    for _ in 0..100 {
        assert!(simplex.inside(Point3::new(0.1, 0.1, 0.1)));
        assert!(!simplex.inside(Point3::new(-0.1, 0.1, 0.1)));
    }
}

#[test]
fn indexed_winding_matches_exhaustive_polygon_fans() {
    let cube = |origin: Point3, scale: f64| {
        let positions = [
            [0., 0., 0.],
            [1., 0., 0.],
            [1., 1., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 0., 1.],
            [1., 1., 1.],
            [0., 1., 1.],
        ]
        .map(|[x, y, z]| origin + Vector3::new(x, y, z) * scale)
        .to_vec();
        let faces = Faces::from_iter([
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ]);
        PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            faces,
        )
        .expect("quad faces")
    };
    let meshes = (0..12)
        .map(|i| {
            let x = f64::from(i);
            (cube(Point3::new(2. * x, x.sin(), 0.), 1.0), i % 3 != 0)
        })
        .collect::<Vec<_>>();
    let index = SignedRayIndex::new(meshes.iter().map(|(mesh, forward)| (mesh, *forward)));
    assert!(index.hierarchy.is_some());
    let directions = [
        Vector3::unit_x(),
        -Vector3::unit_x(),
        Vector3::unit_y(),
        -Vector3::unit_z(),
        Vector3::new(1., 0.17, -0.29),
        Vector3::new(-0.75, 0.5, 0.21),
        Vector3::new(1e-300, 0., 1.),
        Vector3::new(0., 0., 0.),
    ];
    let origins = (-12..=72)
        .flat_map(|i| {
            let x = f64::from(i) * 0.4;
            [
                Point3::new(x, 0.23, 0.45),
                Point3::new(x, -0.3, 0.),
                Point3::new(x, 0.81, 1.0),
            ]
        })
        .chain([
            Point3::new(1e150, 1e150, 1e150),
            Point3::new(f64::NAN, 0., 0.),
        ]);
    for point in origins {
        for direction in directions {
            let expected = meshes
                .iter()
                .map(|(mesh, forward)| {
                    mesh.signed_crossing_faces(point, direction) * if *forward { 1 } else { -1 }
                })
                .sum::<isize>();
            assert_eq!(
                index.signed_crossing_faces(point, direction),
                expected,
                "{point:?} {direction:?}"
            );
        }
    }
    let empty = SignedRayIndex::new(std::iter::empty::<(&PolygonMesh, bool)>());
    assert_eq!(
        empty.signed_crossing_faces(Point3::new(0., 0., 0.), Vector3::unit_x()),
        0
    );
}

/// Benchmark the complete public index build, including gathering its native
/// triangle records. The same fixture can be run against previous revisions.
#[test]
#[ignore = "release-only signed ray index timing"]
fn benchmark_signed_ray_hierarchy_build() -> anyhow::Result<()> {
    use std::{hint::black_box, time::Instant};
    let mut positions = Vec::new();
    let mut faces = Vec::new();
    for layer in 0..64 {
        for row in 0..64 {
            for col in 0..64 {
                let x = f64::from(col) + f64::from(row).sin() * 0.25;
                let y = f64::from(row) + f64::from(layer).cos() * 0.25;
                let z = f64::from(layer) + f64::from(col).sin() * 0.25;
                let start = positions.len();
                positions.extend([
                    Point3::new(x, y, z),
                    Point3::new(x + 0.4, y, z + 0.2),
                    Point3::new(x, y + 0.4, z + 0.1),
                ]);
                faces.push([start, start + 1, start + 2]);
            }
        }
    }
    // Deterministic permutation makes build measurements independent of
    // trivially sorted input; the actual narrow phase keeps original fans.
    let count = faces.len();
    for i in 0..count {
        faces.swap(
            i,
            (i.wrapping_mul(2654435761_usize).wrapping_add(17)) % count,
        );
    }
    let mesh = PolygonMesh::new(
        StandardAttributes {
            positions,
            ..Default::default()
        },
        Faces::from_iter(faces),
    )?;
    for round in 0..4 {
        let start = Instant::now();
        let index = SignedRayIndex::new([(black_box(&mesh), true)]);
        let build = start.elapsed();
        let start = Instant::now();
        let checksum = (0..128)
            .map(|i| {
                let p = Point3::new(f64::from(i % 64) + 0.1, f64::from(i / 2) + 0.12, -1.0);
                index.signed_crossing_faces(black_box(p), Vector3::new(0.013, 0.029, 1.0))
            })
            .sum::<isize>();
        eprintln!(
            "triangles={count} round={round} build={build:?} queries={:?} checksum={checksum}",
            start.elapsed()
        );
        black_box(index);
    }
    Ok(())
}
