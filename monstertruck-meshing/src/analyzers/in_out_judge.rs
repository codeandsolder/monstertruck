use super::{InnerSpace, Matrix3, Origin, Point3, PolygonMesh, SquareMatrix, Vector3, hash};
use derive_more::{Deref, DerefMut};
use std::borrow::Borrow;
use std::ops::Range;

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

#[derive(Clone, Copy, Debug)]
struct RayBounds {
    low: [f64; 3],
    high: [f64; 3],
}

impl RayBounds {
    fn from_triangle(triangle: Triangle) -> Self {
        Self {
            low: std::array::from_fn(|axis| {
                triangle
                    .iter()
                    .map(|p| p[axis])
                    .fold(f64::INFINITY, f64::min)
            }),
            high: std::array::from_fn(|axis| {
                triangle
                    .iter()
                    .map(|p| p[axis])
                    .fold(f64::NEG_INFINITY, f64::max)
            }),
        }
    }

    fn union(self, other: Self) -> Self {
        Self {
            low: std::array::from_fn(|axis| self.low[axis].min(other.low[axis])),
            high: std::array::from_fn(|axis| self.high[axis].max(other.high[axis])),
        }
    }

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

#[derive(Clone, Debug)]
enum NodeKind {
    Leaf(Range<usize>),
    Branch([usize; 2]),
}

#[derive(Clone, Debug)]
struct Node {
    bounds: RayBounds,
    kind: NodeKind,
}

/// Immutable spatial index for repeated signed ray queries.
///
/// Each input mesh contributes its original fan triangles and topological
/// orientation. The index changes only candidate discovery; crossing and
/// winding use the same narrow-phase test as `PolygonMesh`.
#[derive(Clone, Debug)]
pub struct SignedRayIndex {
    triangles: Vec<OrientedTriangle>,
    nodes: Vec<Node>,
    root: Option<usize>,
}

impl SignedRayIndex {
    /// Build an owned index from oriented meshes. A false orientation
    /// subtracts that mesh's original winding contribution.
    #[must_use]
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
                    bounds: RayBounds::from_triangle(triangle),
                });
            }
        }
        let mut nodes = Vec::new();
        let root = if finite && !entries.is_empty() {
            Self::build(&mut entries, 0, &mut nodes)
        } else {
            None
        };
        Self {
            triangles: entries,
            nodes,
            root,
        }
    }

    fn build(
        entries: &mut [OrientedTriangle],
        offset: usize,
        nodes: &mut Vec<Node>,
    ) -> Option<usize> {
        const LEAF_TRIANGLES: usize = 8;
        let bounds = entries
            .iter()
            .map(|entry| entry.bounds)
            .reduce(RayBounds::union)?;
        let kind = if entries.len() <= LEAF_TRIANGLES {
            NodeKind::Leaf(offset..offset + entries.len())
        } else {
            let axis = (0..3).max_by(|&a, &b| {
                (bounds.high[a] - bounds.low[a]).total_cmp(&(bounds.high[b] - bounds.low[b]))
            })?;
            let middle = entries.len() / 2;
            let center =
                |entry: &OrientedTriangle| entry.bounds.low[axis].midpoint(entry.bounds.high[axis]);
            entries.select_nth_unstable_by(middle, |a, b| center(a).total_cmp(&center(b)));
            let (left, right) = entries.split_at_mut(middle);
            NodeKind::Branch([
                Self::build(left, offset, nodes)?,
                Self::build(right, offset + middle, nodes)?,
            ])
        };
        let id = nodes.len();
        nodes.push(Node { bounds, kind });
        Some(id)
    }

    fn sum_triangles(entries: &[OrientedTriangle], ray: Ray) -> isize {
        entries
            .iter()
            .map(|entry| entry.sign * entry.triangle.signed_crossing(ray))
            .sum()
    }

    fn query(&self, node: usize, ray: Ray) -> isize {
        let node = &self.nodes[node];
        if !node.bounds.intersects(ray) {
            return 0;
        }
        match &node.kind {
            NodeKind::Leaf(range) => Self::sum_triangles(&self.triangles[range.clone()], ray),
            NodeKind::Branch([left, right]) => self.query(*left, ray) + self.query(*right, ray),
        }
    }

    /// Count signed crossings of a positive ray. Nonfinite input uses the
    /// exhaustive narrow phase instead of unsafe spatial rejection.
    #[must_use]
    pub fn signed_crossing_faces(&self, point: Point3, ray_direction: Vector3) -> isize {
        let ray = Ray {
            origin: point,
            direction: ray_direction,
        };
        self.root
            .filter(|_| {
                (0..3).all(|axis| point[axis].is_finite() && ray_direction[axis].is_finite())
            })
            .map_or_else(
                || Self::sum_triangles(&self.triangles, ray),
                |root| self.query(root, ray),
            )
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
    assert!(index.nodes.len() > 1);
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
