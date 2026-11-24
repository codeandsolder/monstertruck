//! Source-backed planar repair proposals and shared point reconciliation.
//! Original evidence stays immutable. An approximation is applied only if
//! every affected source vertex can meet all of its incident carriers.

use super::curve_evidence::measure_path_curves;
use super::{
    RegionBoundaryError, RegionBoundaryIndex, SharedTrimGeometry, SharedTrimOptions,
    SurfaceCarrier, SurfaceRegionPartition, trim_intersections::supplied_circles,
};
use crate::{InnerSpace, Point3, Vector3};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Explicit policy for changing the source geometry. Independent trim-fit
/// tolerances remain separate and are never inflated by this budget.
#[derive(Clone, Copy, Debug)]
pub struct ReconciliationOptions {
    /// Maximum Euclidean displacement of each original source vertex.
    pub max_displacement: f64,
    /// Angular resolution for treating independently fitted planes as parallel.
    pub parallel_angle_tolerance: f64,
    /// Maximum deviation of original oriented planar support normals.
    pub max_normal_angle: f64,
    /// Numerical intersection residual allowed after reconciliation.
    pub intersection_tolerance: f64,
}

/// One accepted independently anchored plane and its original ownership.
#[derive(Debug)]
pub struct ReconciledPlane {
    /// Active region ID after consolidation.
    pub region: usize,
    /// Original support IDs consolidated into this region.
    pub supports: Vec<usize>,
    /// Final offset of the plane along its independently observed normal.
    pub normal_shift: f64,
    /// Maximum original point-to-plane residual over its source vertices.
    pub max_plane_residual: f64,
}

/// A reconciled partition with a separate vertex table; source IDs unchanged.
#[derive(Debug)]
pub struct ReconciledPartition {
    /// Selected carriers and regenerated exact source contact graph.
    pub partition: SurfaceRegionPartition,
    /// Reconciled source vertices indexed by original canonical ID.
    pub positions: Vec<Point3>,
    /// Accepted planar consolidations.
    pub planes: Vec<ReconciledPlane>,
    /// Maximum source-to-reconciled point distance.
    pub max_displacement: f64,
    /// Number of rejected anchored plane candidates.
    pub rejected: usize,
}

#[derive(Clone, Copy)]
struct PlaneSeed {
    origin: Point3,
    normal: Vector3,
    region: usize,
}

impl RegionBoundaryIndex<'_> {
    /// Reconcile planar fallback regions using circular constraints from
    /// independently recovered neighboring carriers. All source ownership
    /// and winding are preserved. The result has its own position table.
    ///
    /// # Errors
    /// Invalid budgets or inconsistent partition data return an error.
    pub fn reconcile_planes(
        &self,
        partition: SurfaceRegionPartition,
        trim_options: SharedTrimOptions,
        options: ReconciliationOptions,
    ) -> Result<ReconciledPartition, RegionBoundaryError> {
        if !options.max_displacement.is_finite()
            || options.max_displacement <= 0.0
            || !options.parallel_angle_tolerance.is_finite()
            || options.parallel_angle_tolerance <= 0.0
            || options.parallel_angle_tolerance >= options.max_normal_angle
            || !options.max_normal_angle.is_finite()
            || !(0.0..std::f64::consts::FRAC_PI_2).contains(&options.max_normal_angle)
            || !options.intersection_tolerance.is_finite()
            || options.intersection_tolerance <= 0.0
            || partition.region_by_support.len() != self.evidence.planar_supports.len()
            || partition.contacts != self.partition_contacts(&partition.region_by_polygon)?
            || self
                .evidence
                .source_planar_support
                .iter()
                .zip(&partition.region_by_polygon)
                .any(|(&support, &region)| {
                    partition.region_by_support.get(support) != Some(&region)
                })
        {
            return Err(RegionBoundaryError::InvalidPartition);
        }
        let seeds = self.plane_seeds(&partition, trim_options, options);
        let mut result = ReconciledPartition {
            partition,
            positions: self.evidence.positions.clone(),
            planes: Vec::new(),
            max_displacement: 0.0,
            rejected: 0,
        };
        let incidence = self.vertex_supports()?;
        for seed in seeds {
            if !matches!(
                result.partition.carriers.get(seed.region),
                Some(Some(SurfaceCarrier::Plane { .. }))
            ) {
                continue;
            }
            let Some(proposal) = self.grow_plane(&result.partition, seed, options) else {
                continue;
            };
            if !self.apply_plane(
                &mut result,
                proposal,
                seed,
                &incidence,
                trim_options,
                options,
            )? {
                result.rejected += 1;
            }
        }
        Ok(result)
    }

    fn plane_seeds(
        &self,
        partition: &SurfaceRegionPartition,
        trim: SharedTrimOptions,
        options: ReconciliationOptions,
    ) -> Vec<PlaneSeed> {
        let mut seeds = Vec::new();
        for path in &partition.contacts.paths {
            let Some(right) = path.regions.1 else {
                continue;
            };
            let Some(a) = partition.carriers.get(path.regions.0).copied().flatten() else {
                continue;
            };
            let Some(b) = partition.carriers.get(right).copied().flatten() else {
                continue;
            };
            let (region, normal, curved, plane) = match (a, b) {
                (SurfaceCarrier::Plane { normal, .. }, other)
                    if !matches!(other, SurfaceCarrier::Plane { .. }) =>
                {
                    (path.regions.0, normal, other, a)
                }
                (other, SurfaceCarrier::Plane { normal, .. })
                    if !matches!(other, SurfaceCarrier::Plane { .. }) =>
                {
                    (right, normal, other, b)
                }
                _ => continue,
            };
            let Some(ids) = path.ordered_vertices() else {
                continue;
            };
            // A source-supported circle, rather than an arbitrary plane
            // passing through the endpoints, anchors this candidate.
            let Some(mut fit) = measure_path_curves(&self.evidence.positions, &ids) else {
                continue;
            };
            for circle in supplied_circles(
                &self.evidence.positions,
                &ids,
                &[curved, plane],
                trim.absolute_tolerance,
            ) {
                fit.circle = Some(circle);
                if plane
                    .measure(&fit)
                    .circle
                    .is_some_and(|r| r.upper_bound <= trim.absolute_tolerance)
                {
                    continue;
                }
                if !curved
                    .measure(&fit)
                    .circle
                    .is_some_and(|r| r.upper_bound <= trim.absolute_tolerance)
                    || circle.max_step_radians
                        > std::f64::consts::FRAC_PI_2 + trim.relative_tolerance
                    || circle.normal.cross(normal).magnitude() > options.max_normal_angle.sin()
                {
                    continue;
                }
                let normal = if circle.normal.dot(normal) < 0.0 {
                    -circle.normal
                } else {
                    circle.normal
                };
                seeds.push(PlaneSeed {
                    origin: circle.center,
                    normal,
                    region,
                });
            }
        }
        seeds
    }

    fn vertex_supports(&self) -> Result<Vec<Vec<usize>>, RegionBoundaryError> {
        let mut result = vec![Vec::new(); self.evidence.positions.len()];
        for (id, support) in self.evidence.planar_supports.iter().enumerate() {
            for &polygon in &support.source_face_indices {
                for &vertex in self
                    .evidence
                    .source_polygon_vertices(polygon)
                    .ok_or(RegionBoundaryError::InvalidEvidence)?
                {
                    result
                        .get_mut(vertex)
                        .ok_or(RegionBoundaryError::InvalidEvidence)?
                        .push(id);
                }
            }
        }
        for row in &mut result {
            row.sort_unstable();
            row.dedup();
        }
        Ok(result)
    }

    fn grow_plane(
        &self,
        partition: &SurfaceRegionPartition,
        seed: PlaneSeed,
        options: ReconciliationOptions,
    ) -> Option<ReconciledPlane> {
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for (support, &region) in partition.region_by_support.iter().enumerate() {
            groups.entry(region).or_default().push(support);
        }
        let mut regions = BTreeSet::new();
        let mut visited = BTreeSet::from([seed.region]);
        let mut queue = VecDeque::from([seed.region]);
        let mut max_plane_residual = 0.0_f64;
        while let Some(region) = queue.pop_front() {
            let Some(SurfaceCarrier::Plane { origin, normal }) = partition.carriers[region] else {
                continue;
            };
            let angle = normal
                .cross(seed.normal)
                .magnitude()
                .atan2(normal.dot(seed.normal));
            if angle > options.max_normal_angle {
                continue;
            }
            // Do not absorb a separately resolved parallel plane, even
            // when its separation lies inside the displacement budget.
            if angle <= options.parallel_angle_tolerance
                && (origin - seed.origin).dot(seed.normal).abs() > options.intersection_tolerance
            {
                continue;
            }
            let supports = &groups[&region];
            let mut residual = 0.0_f64;
            for &id in supports {
                let support = &self.evidence.planar_supports[id];
                if support
                    .normal
                    .cross(seed.normal)
                    .magnitude()
                    .atan2(support.normal.dot(seed.normal))
                    + support.max_normal_angle_radians
                    > options.max_normal_angle
                {
                    residual = f64::INFINITY;
                    break;
                }
                for &polygon in &support.source_face_indices {
                    for &v in self.evidence.source_polygon_vertices(polygon)? {
                        residual = residual.max(
                            (self.evidence.positions[v] - seed.origin)
                                .dot(seed.normal)
                                .abs(),
                        );
                    }
                }
            }
            if residual > options.max_displacement {
                continue;
            }
            max_plane_residual = max_plane_residual.max(residual);
            regions.insert(region);
            for &id in supports {
                for &interface in &self.evidence.support_interface_indices[id] {
                    let [a, b] = self.evidence.support_interfaces[interface].supports;
                    let neighbor = partition.region_by_support[if a == id { b } else { a }];
                    if visited.insert(neighbor) {
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        if regions.len() < 2 || max_plane_residual <= options.intersection_tolerance {
            return None;
        }
        let mut supports = regions
            .iter()
            .flat_map(|id| groups[id].iter().copied())
            .collect::<Vec<_>>();
        supports.sort_unstable();
        let boundary = self.boundary(&supports).ok()?;
        if !boundary.pinched_vertices.is_empty() {
            return None;
        }
        if !self.plane_has_anchors(partition, &regions, seed, options.intersection_tolerance) {
            return None;
        }
        Some(ReconciledPlane {
            region: supports[0],
            supports,
            max_plane_residual,
            normal_shift: 0.0,
        })
    }

    fn plane_has_anchors(
        &self,
        partition: &SurfaceRegionPartition,
        regions: &BTreeSet<usize>,
        seed: PlaneSeed,
        tolerance: f64,
    ) -> bool {
        // Require at least three noncollinear independently anchored
        // samples across the entire connected proposed plane region.
        let mut anchors = BTreeSet::new();
        for path in &partition.contacts.paths {
            let Some(right) = path.regions.1 else {
                continue;
            };
            let other = if regions.contains(&path.regions.0) {
                right
            } else if regions.contains(&right) {
                path.regions.0
            } else {
                continue;
            };
            let Some(carrier) = partition.carriers[other] else {
                continue;
            };
            if matches!(carrier, SurfaceCarrier::Plane { .. }) {
                continue;
            }
            for edge in &path.edges {
                for v in edge.vertices {
                    if (self.evidence.positions[v] - seed.origin)
                        .dot(seed.normal)
                        .abs()
                        <= tolerance
                    {
                        anchors.insert(v);
                    }
                }
            }
        }
        noncollinear(&self.evidence.positions, &anchors, tolerance)
    }

    fn apply_plane(
        &self,
        result: &mut ReconciledPartition,
        mut proposal: ReconciledPlane,
        seed: PlaneSeed,
        incidence: &[Vec<usize>],
        trim: SharedTrimOptions,
        options: ReconciliationOptions,
    ) -> Result<bool, RegionBoundaryError> {
        let region = proposal.region;
        let old_regions = proposal
            .supports
            .iter()
            .map(|&id| result.partition.region_by_support[id])
            .collect::<BTreeSet<_>>();
        let mut vertices = BTreeSet::new();
        for &id in &proposal.supports {
            for &polygon in &self.evidence.planar_supports[id].source_face_indices {
                vertices.extend(
                    self.evidence
                        .source_polygon_vertices(polygon)
                        .ok_or(RegionBoundaryError::InvalidEvidence)?
                        .iter()
                        .copied(),
                );
            }
        }
        let mut constraints = Vec::new();
        for v in vertices {
            let neighbors = incidence[v]
                .iter()
                .map(|&id| result.partition.region_by_support[id])
                .collect::<BTreeSet<_>>();
            let other = neighbors
                .difference(&old_regions)
                .map(|id| {
                    result.partition.carriers[*id].ok_or(RegionBoundaryError::InvalidPartition)
                })
                .collect::<Result<Vec<_>, _>>()?;
            constraints.push(VertexConstraints {
                vertex: v,
                source: self.evidence.positions[v],
                others: other,
            });
        }
        let Some(placement) = place_plane(seed, &constraints, options) else {
            return Ok(false);
        };
        let carrier = placement.carrier;
        let changes = placement.changes;
        proposal.normal_shift = placement.shift;
        proposal.max_plane_residual = constraints
            .iter()
            .map(|c| {
                carrier
                    .measure_segment(c.source, c.source)
                    .map_or(f64::INFINITY, |r| r.witnessed_max)
            })
            .fold(0.0_f64, f64::max);
        let mut positions = result.positions.clone();
        let mut touched = BTreeSet::new();
        let mut affected_supports = BTreeSet::new();
        for (v, point, displacement) in changes {
            positions[v] = point;
            if displacement > options.intersection_tolerance {
                touched.insert(v);
                affected_supports.extend(incidence[v].iter().copied());
            }
        }
        if !preserves_polygon_orientation(self.evidence, &positions, &affected_supports) {
            return Ok(false);
        }
        let mut partition = result.partition.clone();
        for &old in &old_regions {
            partition.carriers[old] = None;
        }
        partition.carriers[region] = Some(carrier);
        for &id in &proposal.supports {
            partition.region_by_support[id] = region;
        }
        for (polygon, &support) in self.evidence.source_planar_support.iter().enumerate() {
            partition.region_by_polygon[polygon] = partition.region_by_support[support];
        }
        partition.contacts = self.partition_contacts(&partition.region_by_polygon)?;
        partition.consolidated_supports =
            partition.carriers.len() - partition.carriers.iter().flatten().count();
        if !preserves_contact_geometry(
            result,
            (&partition, &positions),
            &old_regions,
            &touched,
            region,
            trim,
        ) {
            return Ok(false);
        }
        result.partition = partition;
        result.positions = positions;
        result.max_displacement = result
            .positions
            .iter()
            .zip(&self.evidence.positions)
            .map(|(p, q)| (*p - *q).magnitude())
            .fold(0.0_f64, f64::max);
        result.planes.push(proposal);
        Ok(true)
    }
}

fn polygon_area(points: &[Point3], ids: &[usize]) -> Vector3 {
    let origin = points[ids[0]];
    ids.iter()
        .zip(ids.iter().cycle().skip(1))
        .fold(Vector3::new(0.0, 0.0, 0.0), |area, (&a, &b)| {
            area + (points[a] - origin).cross(points[b] - origin)
        })
}

fn preserves_polygon_orientation(
    evidence: &super::MeshEvidence,
    positions: &[Point3],
    supports: &BTreeSet<usize>,
) -> bool {
    supports.iter().all(|&support| {
        evidence.planar_supports[support]
            .source_face_indices
            .iter()
            .all(|&polygon| {
                let Some(ids) = evidence.source_polygon_vertices(polygon) else {
                    return false;
                };
                let old = polygon_area(&evidence.positions, ids);
                let new = polygon_area(positions, ids);
                new.dot(old) > 0.0
                    && ids
                        .iter()
                        .zip(ids.iter().cycle().skip(1))
                        .all(|(&a, &b)| positions[a] != positions[b])
            })
    })
}

#[derive(Default)]
struct ContactQuality {
    unresolved: BTreeSet<[usize; 2]>,
    conics: BTreeSet<[usize; 2]>,
}

fn contact_quality(
    partition: &SurfaceRegionPartition,
    points: &[Point3],
    regions: &BTreeSet<usize>,
    vertices: &BTreeSet<usize>,
    trim: SharedTrimOptions,
) -> ContactQuality {
    let mut quality = ContactQuality::default();
    for path in &partition.contacts.paths {
        if !regions.contains(&path.regions.0)
            && !path.regions.1.is_some_and(|r| regions.contains(&r))
            && !path
                .edges
                .iter()
                .any(|e| e.vertices.iter().any(|v| vertices.contains(v)))
        {
            continue;
        }
        let carriers = std::iter::once(path.regions.0)
            .chain(path.regions.1)
            .filter_map(|id| partition.carriers[id])
            .collect::<Vec<_>>();
        if let Some(fit) = path.segment_analytic_trims(points, &carriers, trim) {
            for piece in fit.pieces {
                let target = match piece.geometry {
                    SharedTrimGeometry::Unresolved => &mut quality.unresolved,
                    SharedTrimGeometry::Conic(_) => &mut quality.conics,
                    _ => continue,
                };
                for offset in 0..piece.edge_count {
                    let edge = path.edges[(piece.first_edge + offset) % path.edges.len()];
                    target.insert(super::edge_key(edge.vertices[0], edge.vertices[1]).0);
                }
            }
        } else {
            quality.unresolved.extend(
                path.edges
                    .iter()
                    .map(|e| super::edge_key(e.vertices[0], e.vertices[1]).0),
            );
        }
    }
    quality
}

fn preserves_contact_geometry(
    before: &ReconciledPartition,
    after: (&SurfaceRegionPartition, &[Point3]),
    regions: &BTreeSet<usize>,
    vertices: &BTreeSet<usize>,
    new_region: usize,
    trim: SharedTrimOptions,
) -> bool {
    let old = contact_quality(
        &before.partition,
        &before.positions,
        regions,
        vertices,
        trim,
    );
    let new = contact_quality(
        after.0,
        after.1,
        &BTreeSet::from([new_region]),
        vertices,
        trim,
    );
    // Independently circle-anchored, displacement-bounded planar repair
    // may replace an already valid noncircular conic by a simpler circular
    // ring (needed by periodic charts). Never introduce unresolved edges.
    // When unresolved coverage is unchanged, require a strict reduction
    // of conic edges without introducing new conic seams.
    new.unresolved.is_subset(&old.unresolved)
        && (new.unresolved.len() < old.unresolved.len()
            || (new.conics.is_subset(&old.conics) && new.conics.len() < old.conics.len()))
}

struct VertexConstraints {
    vertex: usize,
    source: Point3,
    others: Vec<SurfaceCarrier>,
}

struct PlanePlacement {
    carrier: SurfaceCarrier,
    shift: f64,
    changes: Vec<(usize, Point3, f64)>,
    max_displacement: f64,
    squared_displacement: f64,
}

// Exact interval for a plane to intersect the displacement ball around a
// source vertex while satisfying all *linear* neighboring constraints.
// Nonlinear surfaces are checked by the full iterative solver afterwards.
fn admissible_offsets(
    seed: PlaneSeed,
    constraints: &[VertexConstraints],
    options: ReconciliationOptions,
) -> Option<[f64; 2]> {
    let mut range = [f64::NEG_INFINITY, f64::INFINITY];
    for c in constraints {
        let mut basis = Vec::new();
        let mut offset = Vector3::new(0.0, 0.0, 0.0);
        for &carrier in &c.others {
            let SurfaceCarrier::Plane { origin, normal } = carrier else {
                continue;
            };
            let residual = (c.source + offset - origin).dot(normal);
            let independent = orthogonal_component(normal, &basis);
            let length = independent.magnitude();
            if length > options.parallel_angle_tolerance {
                let direction = independent / length;
                offset -= direction * (residual / length);
                basis.push(direction);
            } else if residual.abs() > options.intersection_tolerance {
                return None;
            }
        }
        let available = options
            .max_displacement
            .mul_add(options.max_displacement, -offset.magnitude2());
        if available < 0.0 {
            return None;
        }
        let center = (c.source + offset - seed.origin).dot(seed.normal);
        let span = available.sqrt() * orthogonal_component(seed.normal, &basis).magnitude();
        // Several fully constrained corner vertices may give a
        // zero-width interval. Preserve the intersection under the
        // declared numerical residual, not bitwise equality of offsets.
        range[0] = range[0].max(center - span - options.intersection_tolerance);
        range[1] = range[1].min(center + span + options.intersection_tolerance);
    }
    (range[0] <= range[1]).then_some(range)
}

fn orthogonal_component(vector: Vector3, basis: &[Vector3]) -> Vector3 {
    basis.iter().fold(vector, |v, &axis| v - axis * v.dot(axis))
}

fn place_plane(
    seed: PlaneSeed,
    constraints: &[VertexConstraints],
    options: ReconciliationOptions,
) -> Option<PlanePlacement> {
    let [lower, upper] = admissible_offsets(seed, constraints, options)?;
    let mut offsets = constraints
        .iter()
        .map(|c| {
            (c.source - seed.origin)
                .dot(seed.normal)
                .clamp(lower, upper)
        })
        .collect::<Vec<_>>();
    offsets.extend([
        lower,
        upper,
        lower.midpoint(upper),
        0.0_f64.clamp(lower, upper),
    ]);
    offsets.sort_unstable_by(f64::total_cmp);
    offsets.dedup_by(|a, b| (*a - *b).abs() <= options.intersection_tolerance * 0.1);
    // Include the minimax offset for the linearized constraint system.
    // A bracket is refined only on feasible evaluations; no untested
    // approximation can be committed.
    let mut best = None;
    for shift in offsets {
        consider_placement(&mut best, seed, constraints, options, shift);
    }
    if let Some(current) = &best {
        let mut a = lower;
        let mut b = upper;
        let mut center = current.shift;
        for _ in 0..12 {
            let x = a.midpoint(center);
            let y = center.midpoint(b);
            consider_placement(&mut best, seed, constraints, options, x);
            consider_placement(&mut best, seed, constraints, options, y);
            let next = best.as_ref()?.shift;
            if next < center {
                b = center;
            } else if next > center {
                a = center;
            } else {
                a = x;
                b = y;
            }
            center = next;
        }
    }
    best
}

fn consider_placement(
    best: &mut Option<PlanePlacement>,
    seed: PlaneSeed,
    constraints: &[VertexConstraints],
    options: ReconciliationOptions,
    shift: f64,
) {
    let carrier = SurfaceCarrier::Plane {
        origin: seed.origin + seed.normal * shift,
        normal: seed.normal,
    };
    let mut changes = Vec::with_capacity(constraints.len());
    let mut maximum = 0.0_f64;
    let mut squared = 0.0_f64;
    for c in constraints {
        let mut carriers = c.others.clone();
        carriers.push(carrier);
        let Some(point) = reconcile_point(c.source, &carriers, options) else {
            return;
        };
        let displacement = (point - c.source).magnitude();
        maximum = maximum.max(displacement);
        squared = displacement.mul_add(displacement, squared);
        if best.as_ref().is_some_and(|b| maximum > b.max_displacement) {
            return;
        }
        changes.push((c.vertex, point, displacement));
    }
    if best.as_ref().is_none_or(|b| {
        maximum
            .total_cmp(&b.max_displacement)
            .then(squared.total_cmp(&b.squared_displacement))
            .is_lt()
    }) {
        *best = Some(PlanePlacement {
            carrier,
            shift,
            changes,
            max_displacement: maximum,
            squared_displacement: squared,
        });
    }
}

fn noncollinear(points: &[Point3], ids: &BTreeSet<usize>, tolerance: f64) -> bool {
    let Some(&first) = ids.first() else {
        return false;
    };
    let origin = points[first];
    let Some(direction) = ids
        .iter()
        .map(|&id| points[id] - origin)
        .max_by(|a, b| a.magnitude2().total_cmp(&b.magnitude2()))
    else {
        return false;
    };
    let span = direction.magnitude();
    span > tolerance
        && ids
            .iter()
            .any(|&id| direction.cross(points[id] - origin).magnitude() > tolerance * span)
}

// Gauss-Newton with a rank-revealing orthogonal constraint basis in R³.
// Near-dependent directions do not add unstable displacement components.
pub(super) fn reconcile_point(
    source: Point3,
    carriers: &[SurfaceCarrier],
    options: ReconciliationOptions,
) -> Option<Point3> {
    let mut point = source;
    for _ in 0..16 {
        let mut basis = Vec::<(Vector3, f64)>::new();
        let mut maximum = 0.0_f64;
        for &carrier in carriers {
            let (normal, error) = carrier.signed_constraint(point)?;
            maximum = maximum.max(error.abs());
            let (orthogonal, rhs) =
                basis
                    .iter()
                    .fold((normal, -error), |(n, r), &(axis, value)| {
                        let projection = n.dot(axis);
                        (n - axis * projection, (-projection).mul_add(value, r))
                    });
            let length = orthogonal.magnitude();
            if length > options.parallel_angle_tolerance {
                basis.push((orthogonal / length, rhs / length));
            }
        }
        if maximum <= options.intersection_tolerance {
            return Some(point);
        }
        let delta = basis
            .iter()
            .fold(Vector3::new(0.0, 0.0, 0.0), |sum, &(axis, value)| {
                sum + axis * value
            });
        point += delta;
        if (point - source).magnitude() > options.max_displacement {
            return None;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> ReconciliationOptions {
        ReconciliationOptions {
            max_displacement: 0.0055,
            max_normal_angle: 0.01,
            parallel_angle_tolerance: 1.0e-8,
            intersection_tolerance: 1.0e-9,
        }
    }

    #[test]
    fn plane_placement_preserves_fully_constrained_bevel_corners()
    -> Result<(), Box<dyn std::error::Error>> {
        let axis = Vector3::unit_x();
        let seed = PlaneSeed {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: axis,
            region: 0,
        };
        let corner = Point3::new(0.005, 3.0, 1.0);
        let constraints = [
            VertexConstraints {
                vertex: 0,
                source: corner,
                others: vec![
                    SurfaceCarrier::Plane {
                        origin: corner,
                        normal: Vector3::unit_y(),
                    },
                    SurfaceCarrier::Plane {
                        origin: corner,
                        normal: Vector3::new(1.0, 0.0, 0.05).normalize(),
                    },
                    SurfaceCarrier::Plane {
                        origin: corner,
                        normal: Vector3::new(0.0, 1.0, 0.05).normalize(),
                    },
                ],
            },
            VertexConstraints {
                vertex: 1,
                source: Point3::new(0.0, 1.0, 0.0),
                others: vec![SurfaceCarrier::LinearRevolution {
                    origin: seed.origin,
                    axis,
                    radius_at_origin: 1.0,
                    slope: 0.0,
                }],
            },
        ];
        let placement = place_plane(seed, &constraints, options()).ok_or("placement")?;
        assert!((placement.shift - 0.005).abs() <= 2.0e-9);
        assert!((placement.changes[0].1 - corner).magnitude() < 1.0e-9);
        assert!((placement.changes[1].1 - Point3::new(0.005, 1.0, 0.0)).magnitude() < 2.0e-9);
        assert!(placement.max_displacement <= options().max_displacement);
        Ok(())
    }

    #[test]
    fn conflicting_fixed_corners_cannot_be_merged() {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let seed = PlaneSeed {
            origin,
            normal: Vector3::unit_x(),
            region: 0,
        };
        let constraints = [origin, Point3::new(0.005, 1.0, 0.0)]
            .into_iter()
            .enumerate()
            .map(|(vertex, source)| VertexConstraints {
                vertex,
                source,
                others: [Vector3::unit_x(), Vector3::unit_y(), Vector3::unit_z()]
                    .into_iter()
                    .map(|normal| SurfaceCarrier::Plane {
                        origin: source,
                        normal,
                    })
                    .collect(),
            })
            .collect::<Vec<_>>();
        assert!(place_plane(seed, &constraints, options()).is_none());
    }

    #[test]
    fn reconcile_point_respects_all_nonlinear_and_linear_constraints()
    -> Result<(), Box<dyn std::error::Error>> {
        let origin = Point3::new(0.0, 0.0, 0.0);
        let point = Point3::new(1.001, 0.0, 0.001);
        let carriers = [
            SurfaceCarrier::CircularRevolution {
                axial_center: origin,
                axis: Vector3::unit_z(),
                radial_center: 0.0,
                profile_radius: 1.0,
            },
            SurfaceCarrier::Plane {
                origin,
                normal: Vector3::unit_z(),
            },
            SurfaceCarrier::Plane {
                origin,
                normal: Vector3::unit_y(),
            },
        ];
        let result = reconcile_point(point, &carriers, options()).ok_or("point")?;
        assert!((result - Point3::new(1.0, 0.0, 0.0)).magnitude() < 1.0e-9);
        Ok(())
    }
}
