//! Reuse the same projected triangle intervals for many point clouds.
//! Every query retains the original sweep's inclusive interval and strict
//! point-to-triangle distance predicates.

use super::super::triangle_projection_interval;
use super::{
    EuclideanSpace, InnerSpace, Origin, Point3, PolygonMesh, Vector3, distance2_point_triangle,
    hash, nonpositive_tolerance,
};

#[derive(Clone, Copy, Debug)]
struct Interval {
    low: f64,
    high: f64,
    subtree_high: f64,
    vertices: [usize; 3],
}

/// Immutable projection index for repeated point-cloud proximity queries
/// against one polygon at a fixed tolerance. It uses the same projected
/// intervals and distance test as the one-shot sweep.
pub struct ProjectionNeighborhood<'a> {
    polygon: &'a PolygonMesh,
    tolerance: f64,
    projection: Option<Vector3>,
    intervals: Vec<Interval>,
}

impl<'a> ProjectionNeighborhood<'a> {
    /// Build an index for the supplied polygon and distance tolerance.
    ///
    /// # Panics
    /// Panics for a negative tolerance, as for the one-shot neighborhood query.
    #[must_use = "retain the projection index for repeated neighborhood queries"]
    pub fn new(polygon: &'a PolygonMesh, tolerance: f64) -> Self {
        nonpositive_tolerance!(tolerance, 0.0);
        let projection = polygon
            .positions()
            .first()
            .copied()
            .map(hash::take_one_unit);
        let mut intervals = Vec::new();
        let projection = projection
            .filter(|unit| tolerance.is_finite() && (0..3).all(|axis| unit[axis].is_finite()));
        let projection = projection.filter(|&unit| {
            for face in polygon.faces().triangle_iter() {
                let vertices = face.map(|vertex| vertex.pos);
                let triangle = vertices.map(|index| polygon.positions()[index]);
                if !triangle
                    .iter()
                    .all(|p| (0..3).all(|axis| p[axis].is_finite()))
                {
                    return false;
                }
                if (triangle[1] - triangle[0])
                    .cross(triangle[2] - triangle[0])
                    .so_small()
                {
                    continue;
                }
                let (low, high) = triangle_projection_interval(triangle, unit, tolerance);
                if !low.is_finite() || !high.is_finite() || low > high {
                    return false;
                }
                intervals.push(Interval {
                    low,
                    high,
                    subtree_high: high,
                    vertices,
                });
            }
            true
        });
        if projection.is_some() {
            intervals.sort_unstable_by(|a, b| a.low.total_cmp(&b.low));
            Self::prepare(&mut intervals);
        } else {
            intervals.clear();
        }
        Self {
            polygon,
            tolerance,
            projection,
            intervals,
        }
    }

    // A median of the low-sorted slice is the implicit root. Its subtree
    // maximum high endpoint proves whether any interval can contain a query.
    fn prepare(intervals: &mut [Interval]) -> f64 {
        let middle = intervals.len() / 2;
        let (left, tail) = intervals.split_at_mut(middle);
        let Some((root, right)) = tail.split_first_mut() else {
            return f64::NEG_INFINITY;
        };
        root.subtree_high = root.high.max(Self::prepare(left)).max(Self::prepare(right));
        root.subtree_high
    }

    fn contains(
        intervals: &[Interval],
        positions: &[Point3],
        coordinate: f64,
        point: Point3,
        tolerance2: f64,
    ) -> bool {
        let middle = intervals.len() / 2;
        let (left, tail) = intervals.split_at(middle);
        let Some((root, right)) = tail.split_first() else {
            return false;
        };
        if coordinate > root.subtree_high || coordinate < intervals[0].low {
            return false;
        }
        if root.low <= coordinate {
            if coordinate <= root.high
                && distance2_point_triangle(point, root.vertices.map(|index| positions[index]))
                    < tolerance2
            {
                return true;
            }
            if Self::contains(right, positions, coordinate, point, tolerance2) {
                return true;
            }
        }
        Self::contains(left, positions, coordinate, point, tolerance2)
    }

    /// Whether every point is strictly within the requested distance of
    /// some nondegenerate polygon triangle.
    #[must_use = "the neighborhood classification is returned by this query"]
    pub fn includes(&self, points: &[Point3]) -> bool {
        let Some(unit) = self.projection else {
            return super::sort_end_points::pointcloud_in_polygon_neighborhood(
                self.polygon,
                points,
                self.tolerance,
            );
        };
        // Preserve the original total-order handling for unusual floating
        // inputs instead of introducing a different NaN/overflow policy.
        if points
            .iter()
            .any(|p| !(0..3).all(|axis| p[axis].is_finite()) || !p.to_vec().dot(unit).is_finite())
        {
            return super::sort_end_points::pointcloud_in_polygon_neighborhood(
                self.polygon,
                points,
                self.tolerance,
            );
        }
        let tolerance2 = self.tolerance * self.tolerance;
        points.iter().all(|&point| {
            Self::contains(
                &self.intervals,
                self.polygon.positions(),
                point.to_vec().dot(unit),
                point,
                tolerance2,
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Faces, StandardAttributes};

    #[test]
    fn projection_index_matches_one_shot_sweep() -> anyhow::Result<()> {
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        for x in 0..12 {
            for y in 0..9 {
                let offset = positions.len();
                let x = f64::from(x) * 0.5;
                let y = f64::from(y) * 0.5;
                positions.extend([
                    Point3::new(x, y, 0.0),
                    Point3::new(x + 0.4, y, 0.0),
                    Point3::new(x, y + 0.4, 0.0),
                ]);
                faces.push([offset, offset + 1, offset + 2]);
            }
        }
        // Degenerate triangle must remain excluded, as in the sweep.
        faces.push([0, 0, 0]);
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        for tolerance in [0.0, 1e-6, 0.05, 0.25, 1.0] {
            let index = ProjectionNeighborhood::new(&mesh, tolerance);
            for x in -3..27 {
                for y in -2..21 {
                    let p = Point3::new(f64::from(x) * 0.25, f64::from(y) * 0.25, tolerance * 0.5);
                    let q = Point3::new(p.x, p.y, tolerance);
                    for points in [&[p][..], &[p, q][..], &[][..]] {
                        assert_eq!(
                            index.includes(points),
                            super::super::sort_end_points::pointcloud_in_polygon_neighborhood(
                                &mesh, points, tolerance
                            ),
                            "{points:?}, {tolerance}"
                        );
                    }
                }
            }
        }
        Ok(())
    }
}
