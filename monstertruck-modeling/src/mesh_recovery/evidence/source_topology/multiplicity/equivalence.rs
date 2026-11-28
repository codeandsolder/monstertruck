//! Exact oriented geometric equivalence of regular source fragments.
//!
//! Two finite oriented planar chains with identical boundary chains differ
//! by a compactly supported closed two-chain in a plane, which is zero.
//! Each face must be simple, exact-coplanar and consistently oriented. Thus
//! the signed coverage, including multiplicity, agrees independent of its
//! triangulation. The original source polygons are never altered.
use super::super::super::super::{InputPolygon, seams};
use crate::{Point3, mesh_recovery::planar_predicates};
use robust::{Coord, Coord3D, orient2d, orient3d};
use std::collections::{BTreeMap, BTreeSet};

pub(super) type Profile = Vec<(usize, bool)>;

pub(super) fn polygon_key(polygon: &InputPolygon) -> Vec<usize> {
    let vertices = &polygon.vertices;
    // Canonical oriented cycle, without allowing orientation reversal.
    let Some(minimum) = vertices.iter().min() else {
        return Vec::new();
    };
    vertices
        .iter()
        .enumerate()
        .filter(|(_, id)| *id == minimum)
        .map(|(i, _)| {
            vertices[i..]
                .iter()
                .chain(&vertices[..i])
                .copied()
                .collect::<Vec<_>>()
        })
        .min()
        .unwrap_or_default()
}

const fn coord3(p: Point3) -> Coord3D<f64> {
    Coord3D {
        x: p.x,
        y: p.y,
        z: p.z,
    }
}
const fn project(p: Point3, axis: usize) -> Coord<f64> {
    // Point3 fields are accessed via a cyclic projection without arithmetic.
    match axis {
        0 => Coord { x: p.y, y: p.z },
        1 => Coord { x: p.z, y: p.x },
        _ => Coord { x: p.x, y: p.y },
    }
}

struct Plane {
    basis: [Point3; 3],
    axis: usize,
}

impl Plane {
    fn from_profile(positions: &[Point3], edges: &[[usize; 2]], profile: &Profile) -> Option<Self> {
        let mut boundary = BTreeSet::new();
        for &(edge, forward) in profile {
            if profile.binary_search(&(edge, !forward)).is_err() {
                boundary.extend(edges[edge]);
            }
        }
        Self::from_points(boundary.into_iter().map(|id| positions[id]))
    }

    fn from_points(mut points: impl Iterator<Item = Point3>) -> Option<Self> {
        let first = points.next()?;
        let second = points.next()?;
        let third = points.find(|&p| !seams::on_line(first, second, p))?;
        let axis = (0..3).max_by(|&a, &b| {
            orient2d(project(first, a), project(second, a), project(third, a))
                .abs()
                .total_cmp(
                    &orient2d(project(first, b), project(second, b), project(third, b)).abs(),
                )
        })?;
        let determinant = orient2d(
            project(first, axis),
            project(second, axis),
            project(third, axis),
        );
        (determinant.is_finite() && determinant != 0.0).then_some(Self {
            basis: [first, second, third],
            axis,
        })
    }

    fn polygon_sign(&self, positions: &[Point3], polygon: &InputPolygon) -> Option<bool> {
        let [a, b, c] = self.basis.map(coord3);
        if polygon
            .vertices
            .iter()
            .any(|&id| orient3d(a, b, c, coord3(positions[id])) != 0.0)
        {
            return None;
        }
        let points = polygon
            .vertices
            .iter()
            .map(|&id| project(positions[id], self.axis))
            .collect::<Vec<_>>();
        planar_predicates::polygon_sign(&points, std::convert::identity)
    }

    fn contains_fragment(
        &self,
        positions: &[Point3],
        polygons: &[InputPolygon],
        fragment: &[usize],
    ) -> bool {
        let mut signs = fragment
            .iter()
            .map(|&id| self.polygon_sign(positions, &polygons[id]));
        let Some(Some(first)) = signs.next() else {
            return false;
        };
        signs.all(|sign| sign == Some(first))
    }
}

// An exactly reversed vertex cycle denotes an opposed copy only when the
// source polygon is a simple, nonzero-area, exactly planar region.
pub(super) fn simple_planar_polygon(positions: &[Point3], polygon: &InputPolygon) -> bool {
    simple_planar_basis(positions, polygon).is_some()
}

pub(super) fn simple_planar_basis(
    positions: &[Point3],
    polygon: &InputPolygon,
) -> Option<[Point3; 3]> {
    let plane = Plane::from_points(polygon.vertices.iter().map(|&id| positions[id]))?;
    plane.polygon_sign(positions, polygon)?;
    Some(plane.basis)
}

const fn take_class(next: &mut usize) -> usize {
    let class = *next;
    *next += 1;
    class
}

// Equivalent classes must have identical oriented incidence profiles,
// including both directions of an edge internal to a fragment. This makes
// quotienting partial search states a sound symmetry reduction.
pub(super) fn classes(
    positions: &[Point3],
    polygons: &[InputPolygon],
    fragments: &[&[usize]],
    edges: &[[usize; 2]],
    profiles: &[Profile],
) -> Vec<usize> {
    let mut buckets = BTreeMap::<&Profile, Vec<usize>>::new();
    for (id, profile) in profiles.iter().enumerate() {
        buckets.entry(profile).or_default().push(id);
    }
    let mut classes = vec![0; fragments.len()];
    let mut next = 0;
    for (profile, ids) in buckets {
        if let [id] = ids.as_slice() {
            classes[*id] = next;
            next += 1;
            continue;
        }
        let mut exact = BTreeMap::<Vec<Vec<usize>>, Vec<usize>>::new();
        for id in ids {
            let mut keys = fragments[id]
                .iter()
                .map(|&p| polygon_key(&polygons[p]))
                .collect::<Vec<_>>();
            keys.sort_unstable();
            exact.entry(keys).or_default().push(id);
        }
        // Identical polygon cycles already prove exact equivalence, even
        // for nonplanar fragments. Only distinct tessellations need the
        // additional planar-chain proof, once per exact geometry class.
        let plane = (exact.len() > 1)
            .then(|| Plane::from_profile(positions, edges, profile))
            .flatten();
        let mut planar_class = None;
        for ids in exact.into_values() {
            let planar = plane
                .as_ref()
                .is_some_and(|p| p.contains_fragment(positions, polygons, fragments[ids[0]]));
            let class = if planar {
                *planar_class.get_or_insert_with(|| take_class(&mut next))
            } else {
                take_class(&mut next)
            };
            for id in ids {
                classes[id] = class;
            }
        }
    }
    classes
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::input_polygon_from_ids;

    #[test]
    fn exact_planar_chain_equivalence_is_independent_of_interior_triangulation()
    -> Result<(), Box<dyn std::error::Error>> {
        let positions = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.5, 0.5, 0.0),
            Point3::new(0.5, 0.5, 1e-9),
        ];
        let faces = [
            vec![0, 1, 2, 3],
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
            vec![0, 1, 5],
            vec![1, 2, 5],
            vec![2, 3, 5],
            vec![3, 0, 5],
        ];
        let polygons = faces
            .into_iter()
            .enumerate()
            .map(|(id, face)| input_polygon_from_ids(&positions, face, id, crate::KERNEL_EPSILON))
            .collect::<Result<Vec<_>, _>>()?;
        let edges = [[0, 1], [1, 2], [2, 3], [0, 3]];
        let profile = vec![(0, true), (1, true), (2, true), (3, false)];
        let result = classes(
            &positions,
            &polygons,
            &[&[0], &[1, 2, 3, 4], &[5, 6, 7, 8]],
            &edges,
            &[profile.clone(), profile.clone(), profile],
        );
        assert_eq!(result[0], result[1]);
        assert_ne!(result[0], result[2]);
        Ok(())
    }
}
