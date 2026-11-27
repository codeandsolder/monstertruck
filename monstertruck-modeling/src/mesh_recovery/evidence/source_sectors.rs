//! Exact local radial sector pairing at a multiply used geometric edge.
//! For a positively oriented body the material sector, seen along the
//! canonical edge direction, begins at a reverse use and ends at the next
//! forward use. The opposite convention applies to negative-volume shells.
//! Pairings are only hypotheses until whole-component closure and signed
//! volume are independently certified.

use super::super::{EdgeKey, EdgeUse, InputPolygon, seams};
use crate::{InnerSpace, Point3, Vector3};
use robust::{Coord3D, orient3d};
use std::cmp::Ordering;

/// An accepted oriented angular-sector pairing at an ambiguous source edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SourceSectorPair {
    /// Canonical geometric edge before topological separation.
    pub edge: [usize; 2],
    /// Original input polygon IDs paired across this edge.
    pub polygons: [usize; 2],
    /// The resolved signed volume of the corresponding component is positive.
    pub positive_volume: bool,
}

const fn coord(point: Point3) -> Coord3D<f64> {
    Coord3D {
        x: point.x,
        y: point.y,
        z: point.z,
    }
}

#[derive(Clone, Copy)]
struct Ray {
    polygon: usize,
    forward: bool,
    witness: Point3,
    flipped: bool,
    // Unit E x ray; its dot product distinguishes coincident and opposite
    // rays when the exact four-point orientation is zero.
    transverse_normal: Vector3,
}

const fn scale(vector: Vector3) -> f64 { vector.x.abs().max(vector.y.abs()).max(vector.z.abs()) }

impl Ray {
    fn from_use(
        positions: &[Point3],
        polygons: &[InputPolygon],
        [start, end]: [Point3; 2],
        (polygon, forward): EdgeUse,
    ) -> Option<Self> {
        let source = polygons.get(polygon)?;
        let (witness, cross) = source
            .vertices
            .iter()
            .map(|&id| {
                let point = positions[id];
                (point, seams::oriented_triangle_area(start, end, point))
            })
            .max_by(|a, b| scale(a.1).total_cmp(&scale(b.1)))?;
        let magnitude = scale(cross);
        if magnitude <= 0.0 || !magnitude.is_finite() {
            return None;
        }
        // A warped source n-gon does not define one exact radial ray. Do
        // not arbitrarily order its vertices by an approximate fit normal.
        if source.vertices.iter().any(|&id| {
            orient3d(
                coord(start),
                coord(end),
                coord(witness),
                coord(positions[id]),
            ) != 0.0
        }) {
            return None;
        }
        let normal = (cross / magnitude).normalize();
        let alignment = normal.dot(source.normal);
        if !alignment.is_finite() || alignment == 0.0 {
            return None;
        }
        // Concave polygons may have a distant sample on the exterior
        // half-plane. The oriented polygon normal fixes the local side
        // of this directed boundary edge, not the arbitrary witness.
        let flipped = (alignment > 0.0) != forward;
        Some(Self {
            polygon,
            forward,
            witness,
            flipped,
            transverse_normal: if flipped { -normal } else { normal },
        })
    }
}

struct Fan {
    edge: [Point3; 2],
}

impl Fan {
    fn turn(&self, left: Ray, right: Ray) -> Option<Ordering> {
        let determinant = orient3d(
            coord(self.edge[0]),
            coord(self.edge[1]),
            coord(left.witness),
            coord(right.witness),
        );
        // orient3d(a,b,c,d) has the opposite sign to E dot ((c-a)x(d-a)).
        let signed = if left.flipped == right.flipped {
            -determinant
        } else {
            determinant
        };
        signed.partial_cmp(&0.0)
    }

    fn same_ray(&self, left: Ray, right: Ray) -> bool {
        self.turn(left, right) == Some(Ordering::Equal)
            && left.transverse_normal.dot(right.transverse_normal) > 0.0
    }

    fn half(&self, reference: Ray, ray: Ray) -> Option<bool> {
        Some(match self.turn(reference, ray)? {
            Ordering::Greater => false,
            Ordering::Less => true,
            Ordering::Equal => reference.transverse_normal.dot(ray.transverse_normal) < 0.0,
        })
    }

    fn ordered_groups(&self, rays: &[Ray]) -> Option<Vec<Vec<Ray>>> {
        let reference = *rays.first()?;
        let mut ordered = rays
            .iter()
            .map(|&ray| Some((self.half(reference, ray)?, ray)))
            .collect::<Option<Vec<_>>>()?;
        let mut valid = true;
        ordered.sort_unstable_by(|a, b| {
            a.0.cmp(&b.0)
                .then_with(|| {
                    self.turn(a.1, b.1).map_or_else(
                        || {
                            valid = false;
                            Ordering::Equal
                        },
                        Ordering::reverse,
                    )
                })
                .then(a.1.polygon.cmp(&b.1.polygon))
        });
        if !valid {
            return None;
        }
        let mut groups: Vec<Vec<Ray>> = Vec::new();
        for (_, ray) in ordered {
            if let Some(group) = groups.last_mut()
                && self.same_ray(group[0], ray)
            {
                // Multiple indistinguishable same-oriented uses do not
                // determine unique body topology. No arbitrary index tie.
                if group.iter().any(|other| other.forward == ray.forward) {
                    return None;
                }
                group.push(ray);
            } else {
                groups.push(vec![ray]);
            }
        }
        Some(groups)
    }
}

fn pair_groups(groups: &[Vec<Ray>], positive: bool) -> Option<Vec<[usize; 2]>> {
    let mut ordered = Vec::new();
    // At coincident opposite rays, the exiting material sector precedes
    // the entering one. Pair only sectors of nonzero angular width.
    let entering = !positive;
    for group in groups {
        ordered.extend(group.iter().copied().filter(|r| r.forward != entering));
        ordered.extend(group.iter().copied().filter(|r| r.forward == entering));
    }
    let mut pairs = Vec::new();
    for (left, right) in ordered.iter().zip(ordered.iter().cycle().skip(1)) {
        if left.forward == right.forward {
            return None;
        }
        if left.forward == entering {
            let mut pair = [left.polygon, right.polygon];
            pair.sort_unstable();
            pairs.push(pair);
        }
    }
    Some(pairs)
}

/// Candidate local pairings for positive and negative volume orientations.
pub(super) fn paired_sectors(
    positions: &[Point3],
    polygons: &[InputPolygon],
    edge: EdgeKey,
    uses: &[EdgeUse],
) -> Option<[Option<Vec<[usize; 2]>>; 2]> {
    if uses.len() <= 2 || uses.iter().filter(|u| u.1).count() * 2 != uses.len() {
        return None;
    }
    let edge = edge.map(|id| positions[id]);
    let rays = uses
        .iter()
        .map(|&usage| Ray::from_use(positions, polygons, edge, usage))
        .collect::<Option<Vec<_>>>()?;
    let fan = Fan { edge };
    let groups = fan.ordered_groups(&rays)?;
    Some([pair_groups(&groups, true), pair_groups(&groups, false)])
}

#[cfg(test)]
mod tests {
    use crate::Point3;
    use crate::mesh_recovery::{analyze_mesh_evidence, recover_analytic_mesh};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_topology::shell::ShellCondition;

    #[test]
    fn face_touching_boxes_preserve_both_oriented_closures()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut positions = Vec::new();
        for z in [0.0, 1.0, 2.0] {
            for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
                positions.push(Point3::new(x, y, z));
            }
        }
        let quads = [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ];
        for positive in [true, false] {
            let faces = [0, 4].into_iter().flat_map(|offset| {
                quads.into_iter().flat_map(move |quad| {
                    let q = quad.map(|v| v + offset);
                    let mut triangles = [[q[0], q[1], q[2]], [q[0], q[2], q[3]]];
                    if !positive {
                        triangles.iter_mut().for_each(|t| t.reverse());
                    }
                    triangles
                })
            });
            let mesh = PolygonMesh::new(
                StandardAttributes {
                    positions: positions.clone(),
                    ..Default::default()
                },
                Faces::from_iter(faces),
            )?;
            let recovered = recover_analytic_mesh(&mesh, Default::default())?;
            let evidence = &recovered.evidence;
            assert_eq!(evidence.non_manifold_edge_count, 0);
            assert_eq!(evidence.inconsistent_winding_edge_count, 0);
            assert_eq!(evidence.provenance.separated_vertices.len(), 4);
            // Four perimeter edges and one coincident triangulation diagonal,
            // each with two accepted material-sector pairings.
            assert_eq!(evidence.provenance.sector_pairs.len(), 10);
            assert!(
                evidence
                    .provenance
                    .sector_pairs
                    .iter()
                    .all(|p| p.positive_volume == positive)
            );
            assert!(evidence.provenance.canceled_sheets.is_empty());
            assert_eq!(evidence.provenance.original_polygon_ids.len(), 24);
            let expected_volume = if positive { 2.0 } else { -2.0 };
            assert!(
                (evidence.closed_source_volume().ok_or("closed source")? - expected_volume).abs()
                    < 1.0e-12
            );
            let shell = recovered.complete_shell().ok_or("two closed box shells")?;
            assert_eq!(shell.shell_condition(), ShellCondition::Closed);
            assert_eq!(shell.connected_components().len(), 2);
            assert_eq!(shell.len(), 12);
        }
        Ok(())
    }

    #[test]
    fn duplicate_same_oriented_rays_are_not_arbitrarily_paired()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]].repeat(2)),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, Default::default())?;
        assert_eq!(evidence.non_manifold_edge_count, 6);
        assert!(evidence.provenance.sector_pairs.is_empty());
        assert!(evidence.provenance.separated_vertices.is_empty());
        Ok(())
    }
}
