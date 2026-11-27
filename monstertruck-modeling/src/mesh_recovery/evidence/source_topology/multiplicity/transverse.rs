//! A proper crossing between an actual source boundary segment and the
//! interior of a different exact planar source polygon is incompatible
//! with embedding both fragments in the same shell. No tolerance, projected
//! intersection coordinate or triangulation-dependent diagonal is used.
use super::super::super::super::InputPolygon;
use super::{SourceTransverseCrossing, equivalence::simple_planar_basis};
use crate::Point3;
use robust::{Coord, Coord3D, orient2d, orient3d};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
};

const MAX_SWEEP_COMPARISONS: usize = 1 << 21;

const fn coord(p: Point3) -> Coord3D<f64> {
    Coord3D {
        x: p.x,
        y: p.y,
        z: p.z,
    }
}

fn sign(value: f64) -> Option<Ordering> {
    value.is_finite().then(|| value.partial_cmp(&0.0)).flatten()
}

const fn straddles(a: Ordering, b: Ordering) -> bool {
    !matches!(
        (a, b),
        (Ordering::Less, Ordering::Less) | (Ordering::Greater, Ordering::Greater)
    )
}

const fn opposite(a: Ordering, b: Ordering) -> bool {
    matches!(
        (a, b),
        (Ordering::Less, Ordering::Greater) | (Ordering::Greater, Ordering::Less)
    )
}

// Quotient three-space by the line through the segment, choosing the
// coordinate with the largest direction component. The two remaining
// coordinates have cyclic order. Exact signs are enough for a ray parity
// test in this projection: no ill-conditioned intersection point is made.
struct LineProjection {
    ends: [Point3; 2],
    axis: usize,
    positive: bool,
}

impl LineProjection {
    fn new(ends: [Point3; 2]) -> Option<Self> {
        let axis = (0..3).max_by(|&a, &b| {
            (ends[1][a] - ends[0][a])
                .abs()
                .total_cmp(&(ends[1][b] - ends[0][b]).abs())
        })?;
        let direction = ends[1][axis] - ends[0][axis];
        (direction.is_finite() && direction != 0.0).then_some(Self {
            ends,
            axis,
            positive: direction > 0.0,
        })
    }

    fn side(&self, point: Point3, component: usize) -> Option<Ordering> {
        let p = |p: Point3| Coord {
            x: p[self.axis],
            y: p[(self.axis + component) % 3],
        };
        let [a, b] = self.ends.map(p);
        let order = sign(orient2d(a, b, p(point)))?;
        Some(if self.positive {
            order
        } else {
            order.reverse()
        })
    }

    fn turn(&self, a: Point3, b: Point3) -> Option<Ordering> {
        let [p, q] = self.ends.map(coord);
        let order = sign(orient3d(p, q, coord(a), coord(b)))?;
        // orient3d(p,q,a,b) has the opposite sign to the cyclic projected
        // 2D determinant when the chosen axial direction is positive.
        Some(if self.positive {
            order.reverse()
        } else {
            order
        })
    }

    fn contains_origin(&self, points: impl Iterator<Item = [Point3; 2]>) -> Option<bool> {
        let mut inside = false;
        for [a, b] in points {
            let y = [self.side(a, 2)?, self.side(b, 2)?];
            if !straddles(y[0], y[1]) {
                continue;
            }
            let turn = self.turn(a, b)?;
            if turn == Ordering::Equal {
                // The projected origin lies on an edge (or vertex) only
                // if both coordinate intervals contain it. Exclude all
                // boundary contacts, not merely shared source IDs.
                if straddles(self.side(a, 1)?, self.side(b, 1)?) {
                    return Some(false);
                }
                continue;
            }
            let up = y[1] == Ordering::Greater;
            if (y[0] == Ordering::Greater) != up && (turn == Ordering::Greater) == up {
                inside = !inside;
            }
        }
        Some(inside)
    }
}

fn edges(vertices: &[usize]) -> impl Iterator<Item = [usize; 2]> + '_ {
    vertices
        .iter()
        .copied()
        .zip(vertices.iter().copied().cycle().skip(1))
        .map(Into::into)
}

fn crossing_edge(
    positions: &[Point3],
    source: &InputPolygon,
    target: &InputPolygon,
    basis: [Point3; 3],
) -> Option<[usize; 2]> {
    let [plane_a, plane_b, plane_c] = basis.map(coord);
    for edge in edges(&source.vertices) {
        let ends = edge.map(|id| positions[id]);
        let [Some(start), Some(end)] =
            ends.map(|p| sign(orient3d(plane_a, plane_b, plane_c, coord(p))))
        else {
            continue;
        };
        if !opposite(start, end) {
            continue;
        }
        // Strictly crossing the target's plane makes this parallel
        // projection injective on its planar domain. Origin interior in
        // that projection is equivalent to segment interior hitting face
        // interior, even for concave simple n-gons.
        let Some(projection) = LineProjection::new(ends) else {
            continue;
        };
        if projection
            .contains_origin(edges(&target.vertices).map(|ids| ids.map(|id| positions[id])))
            == Some(true)
        {
            return Some(edge);
        }
    }
    None
}

struct FaceEntry {
    polygon: usize,
    class: usize,
    basis: [Point3; 3],
    low: [f64; 3],
    high: [f64; 3],
}

impl FaceEntry {
    fn new(
        positions: &[Point3],
        polygons: &[InputPolygon],
        polygon: usize,
        class: usize,
    ) -> Option<Self> {
        let source = &polygons[polygon];
        let basis = simple_planar_basis(positions, source)?;
        let mut low = [f64::INFINITY; 3];
        let mut high = [f64::NEG_INFINITY; 3];
        for &id in &source.vertices {
            let p = positions[id];
            for i in 0..3 {
                low[i] = low[i].min(p[i]);
                high[i] = high[i].max(p[i]);
            }
        }
        Some(Self {
            polygon,
            class,
            basis,
            low,
            high,
        })
    }

    fn overlaps(&self, other: &Self) -> bool {
        (0..3).all(|i| self.low[i] <= other.high[i] && other.low[i] <= self.high[i])
    }

    fn crossing(
        &self,
        other: &Self,
        positions: &[Point3],
        polygons: &[InputPolygon],
    ) -> Option<SourceTransverseCrossing> {
        for (a, b) in [(self, other), (other, self)] {
            if let Some(edge) = crossing_edge(
                positions,
                &polygons[a.polygon],
                &polygons[b.polygon],
                b.basis,
            ) {
                return Some(SourceTransverseCrossing {
                    polygons: [a.polygon, b.polygon],
                    edge,
                });
            }
        }
        None
    }
}

fn sweep_axis(entries: &[FaceEntry]) -> usize {
    let span = |i: usize| {
        entries
            .iter()
            .map(|e| e.high[i])
            .fold(f64::NEG_INFINITY, f64::max)
            - entries
                .iter()
                .map(|e| e.low[i])
                .fold(f64::INFINITY, f64::min)
    };
    (0..3)
        .max_by(|&a, &b| span(a).total_cmp(&span(b)))
        .unwrap_or_default()
}

pub(super) fn crossings(
    positions: &[Point3],
    polygons: &[InputPolygon],
    fragments: &[&[usize]],
    classes: &[usize],
    known: &[[usize; 2]],
) -> BTreeMap<[usize; 2], SourceTransverseCrossing> {
    // Equivalent fragments occupy the same oriented domain, even with
    // different triangulations. One representative per class suffices.
    let mut represented = BTreeSet::new();
    let mut entries = fragments
        .iter()
        .zip(classes)
        .filter(|(_, class)| represented.insert(**class))
        .flat_map(|(fragment, &class)| {
            fragment
                .iter()
                .filter_map(move |&id| FaceEntry::new(positions, polygons, id, class))
        })
        .collect::<Vec<_>>();
    let axis = sweep_axis(&entries);
    entries.sort_unstable_by(|a, b| {
        a.low[axis]
            .total_cmp(&b.low[axis])
            .then(a.polygon.cmp(&b.polygon))
    });
    let mut excluded = known.iter().copied().collect::<BTreeSet<_>>();
    let mut witnesses = BTreeMap::new();
    let mut comparisons = 0;
    'sweep: for (i, a) in entries.iter().enumerate() {
        for b in &entries[i + 1..] {
            if b.low[axis] > a.high[axis] {
                break;
            }
            comparisons += 1;
            // Partial discovery remains sound: every collected exclusion
            // is individually certified; uniqueness still requires the
            // independent complete partition solver.
            if comparisons > MAX_SWEEP_COMPARISONS {
                break 'sweep;
            }
            if a.class == b.class || !a.overlaps(b) {
                continue;
            }
            let pair = [a.class.min(b.class), a.class.max(b.class)];
            if excluded.contains(&pair) {
                continue;
            }
            if let Some(witness) = a.crossing(b, positions, polygons) {
                excluded.insert(pair);
                witnesses.insert(pair, witness);
            }
        }
    }
    witnesses
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::input_polygon_from_ids;

    #[test]
    fn projected_parity_handles_concavity_boundaries_and_every_axial_direction() {
        let shape = [
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 1.0],
            [1.0, 1.0],
            [1.0, 3.0],
            [0.0, 3.0],
        ];
        let point = |axis: usize, xy: [f64; 2], height: f64| {
            let mut xyz = [0.0; 3];
            xyz[axis] = height;
            xyz[(axis + 1) % 3] = xy[0];
            xyz[(axis + 2) % 3] = xy[1];
            Point3::new(xyz[0], xyz[1], xyz[2])
        };
        for axis in 0..3 {
            for (xy, expected) in [
                ([0.5, 2.0], true),
                ([2.0, 0.5], true),
                ([2.0, 2.0], false),
                ([1.0, 2.0], false),
                ([1.0, 1.0], false),
                ([3.5, 0.5], false),
                ([-1.0, 0.0], false),
            ] {
                for reverse in [false, true] {
                    let ends = if reverse {
                        [point(axis, xy, 1.0), point(axis, xy, -1.0)]
                    } else {
                        [point(axis, xy, -1.0), point(axis, xy, 1.0)]
                    };
                    let projection = LineProjection::new(ends).expect("nonzero segment");
                    assert_eq!(
                        projection.contains_origin(
                            edges(&[0, 1, 2, 3, 4, 5])
                                .map(|ids| ids.map(|id| point(axis, shape[id], 0.0)))
                        ),
                        Some(expected)
                    );
                }
            }
        }
    }

    #[test]
    fn oblique_lines_use_projection_not_world_axis_coordinates() {
        let shape = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 3.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
        ];
        for (x, y, expected) in [
            (0.5, 2.0, true),
            (2.0, 2.0, false),
            (1.0, 2.0, false),
            (1.0, 1.0, false),
        ] {
            for [dx, dy, dz] in [
                [7.0, -2.0, 1.0],
                [-2.0, 7.0, -1.0],
                [2.0, -1.0, 7.0],
                [-7.0, 2.0, -1.0],
            ] {
                let projection = LineProjection::new([
                    Point3::new(x - dx, y - dy, -dz),
                    Point3::new(x + dx, y + dy, dz),
                ])
                .expect("nonzero segment");
                assert_eq!(
                    projection.contains_origin(
                        edges(&[0, 1, 2, 3, 4, 5]).map(|ids| ids.map(|id| shape[id]))
                    ),
                    Some(expected)
                );
            }
        }
    }

    #[test]
    fn proper_crossings_exclude_only_distinct_fragments_with_certified_faces()
    -> Result<(), Box<dyn std::error::Error>> {
        let positions = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(3.0, 3.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
            Point3::new(1.0, 1.0, -1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(2.0, 1.0, 1.0),
            Point3::new(3.0, 1.0, -1.0),
            Point3::new(3.0, 1.0, 1.0),
            Point3::new(4.0, 1.0, 1.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
        ];
        let polygons = [
            vec![0, 1, 2, 3],
            vec![4, 5, 6],
            vec![7, 8, 9],
            vec![10, 11, 12],
        ]
        .into_iter()
        .enumerate()
        .map(|(id, face)| input_polygon_from_ids(&positions, face, id, crate::KERNEL_EPSILON))
        .collect::<Result<Vec<_>, _>>()?;
        let basis = simple_planar_basis(&positions, &polygons[0]).expect("exact plane");
        assert_eq!(
            crossing_edge(&positions, &polygons[1], &polygons[0], basis),
            Some([4, 5])
        );
        assert_eq!(
            crossing_edge(&positions, &polygons[2], &polygons[0], basis),
            None
        );
        assert_eq!(
            crossing_edge(&positions, &polygons[3], &polygons[0], basis),
            None
        );
        let result = crossings(&positions, &polygons, &[&[0], &[1], &[2]], &[4, 5, 6], &[]);
        assert_eq!(result.len(), 1);
        assert_eq!(result[&[4, 5]].polygons, [1, 0]);
        assert_eq!(result[&[4, 5]].edge, [4, 5]);
        assert!(crossings(&positions, &polygons, &[&[0], &[1]], &[4, 5], &[[4, 5]]).is_empty());
        Ok(())
    }
}
