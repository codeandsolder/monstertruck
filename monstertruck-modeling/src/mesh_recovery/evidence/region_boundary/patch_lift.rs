//! Explicit region-level fallback against already measured neighboring
//! carriers. This does not mutate the strict analytic recovery result.

use super::{
    RegionBoundaryIndex, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionPartition,
    reconcile::{ReconciliationOptions, reconcile_point},
};
use crate::{
    Curve, Edge, InnerSpace, Vertex, builder,
    mesh_recovery::{
        BoundaryPatchApproximation, BoundaryPatchError, BoundaryPatchOptions, PrescribedSourceEdge,
        approximate_boundary_patch,
    },
};
use std::collections::{BTreeMap, BTreeSet};

/// Independent bounds for a source-backed nonanalytic region replacement.
#[derive(Clone, Copy, Debug)]
pub struct RegionPatchOptions {
    /// Original source-vertex displacement and full-facet deviation budgets.
    pub budget: SurfaceRegionBudget,
    /// Numerical residual allowed at retained analytic intersections.
    pub intersection_tolerance: f64,
    /// Rank resolution for independent retained surface constraints.
    pub angular_resolution: f64,
}

/// A connected approximate replacement with exact retained-carrier edges.
#[derive(Debug)]
pub struct LiftedRegionPatch {
    /// Original active source polygon IDs, in patch provenance order.
    pub source_polygons: Vec<usize>,
    /// Oriented manifold patch with measured deviation.
    pub approximation: BoundaryPatchApproximation,
    /// Shared external edges in original source orientation. These exact
    /// identities must also be used to rebuild the adjacent analytic faces.
    pub boundary_edges: Vec<PrescribedSourceEdge>,
}

/// A requested fallback patch could not be certified.
#[derive(Debug, Clone, thiserror::Error)]
pub enum RegionPatchError {
    /// Invalid or pinched source region union, inconsistent partition or options.
    #[error("invalid source region patch")]
    InvalidRegion,
    /// A source junction cannot satisfy all retained incident surfaces.
    #[error("unresolved patch junction {0}")]
    Junction(usize),
    /// A source edge cannot be lifted on its retained carrier within budget.
    #[error("unresolved patch boundary edge {0:?}")]
    Edge([usize; 2]),
    /// Approximate patch assembly failed.
    #[error(transparent)]
    Assembly(#[from] BoundaryPatchError),
}

impl RegionBoundaryIndex<'_> {
    /// Replace an explicitly selected connected region union by a bounded
    /// source-facet approximation. Retained neighboring surfaces determine
    /// exact rational boundary paths; all incident retained constraints
    /// are reconciled jointly before allocating one vertex per source ID.
    ///
    /// This is an independent approximation transaction. The strict
    /// partition and native graph are unchanged. To attach the patch, the
    /// surrounding analytic faces must reuse its returned boundary edges.
    /// The method never claims to identify the original freeform surface.
    ///
    /// # Errors
    /// Rejects ambiguous incidence, incompatible junctions, unresolved
    /// chart paths, exceeded displacement or chord budgets, and invalid
    /// patch shell topology.
    pub fn approximate_region_patch(
        &self,
        partition: &SurfaceRegionPartition,
        regions: &[usize],
        options: RegionPatchOptions,
    ) -> Result<LiftedRegionPatch, RegionPatchError> {
        let selected = regions.iter().copied().collect::<BTreeSet<_>>();
        self.validate_patch_inputs(partition, &selected, options)?;
        let supports = partition
            .region_by_support
            .iter()
            .enumerate()
            .filter_map(|(id, region)| selected.contains(region).then_some(id))
            .collect::<Vec<_>>();
        let boundary = self
            .boundary(&supports)
            .map_err(|_| RegionPatchError::InvalidRegion)?;
        if !boundary.pinched_vertices.is_empty() {
            return Err(RegionPatchError::InvalidRegion);
        }
        let source_polygons = partition
            .region_by_polygon
            .iter()
            .enumerate()
            .filter_map(|(id, region)| selected.contains(region).then_some(id))
            .collect::<Vec<_>>();
        let polygons = source_polygons
            .iter()
            .map(|&id| self.evidence.source_polygon_vertices(id))
            .collect::<Option<Vec<_>>>()
            .ok_or(RegionPatchError::InvalidRegion)?;
        let boundary_ids = boundary
            .loops
            .iter()
            .flatten()
            .flat_map(|edge| edge.vertices)
            .collect::<BTreeSet<_>>();
        let patch_ids = polygons
            .iter()
            .flat_map(|polygon| polygon.iter().copied())
            .collect::<BTreeSet<_>>();
        let vertices =
            self.patch_junctions(partition, &selected, &patch_ids, &boundary_ids, options)?;
        let mut boundary_edges = Vec::new();
        for edge in boundary.loops.iter().flatten() {
            let [a, b] = edge.vertices;
            let front = &vertices[&a];
            let back = &vertices[&b];
            let retained = edge
                .exterior_polygon
                .map(|id| {
                    partition
                        .region_by_polygon
                        .get(id)
                        .and_then(|&region| partition.carriers.get(region))
                        .copied()
                        .flatten()
                        .ok_or(RegionPatchError::InvalidRegion)
                })
                .transpose()?;
            let native = if let Some(carrier) = retained {
                self.lift_patch_edge(carrier, edge.vertices, front, back, options)?
            } else {
                builder::line(front, back)
            };
            boundary_edges.push(PrescribedSourceEdge {
                vertices: edge.vertices,
                edge: native,
            });
        }
        let owners = source_polygons
            .iter()
            .map(|&id| partition.region_by_polygon[id])
            .collect::<Vec<_>>();
        let approximation = approximate_boundary_patch(
            &self.evidence.positions,
            &polygons,
            &boundary_edges,
            BoundaryPatchOptions {
                owner_by_polygon: Some(&owners),
                ..BoundaryPatchOptions::new(options.budget)
            },
        )?;
        Ok(LiftedRegionPatch {
            source_polygons,
            approximation,
            boundary_edges,
        })
    }

    fn validate_patch_inputs(
        &self,
        partition: &SurfaceRegionPartition,
        selected: &BTreeSet<usize>,
        options: RegionPatchOptions,
    ) -> Result<(), RegionPatchError> {
        if selected.is_empty()
            || selected
                .iter()
                .any(|&id| partition.carriers.get(id).is_none_or(Option::is_none))
            || partition.region_by_polygon.len() != self.evidence.source_planar_support.len()
            || partition.region_by_support.len() != self.evidence.planar_supports.len()
            || !options.budget.vertex.is_finite()
            || options.budget.vertex < 0.0
            || !options.budget.chord.is_finite()
            || options.budget.chord <= 0.0
            || !options.intersection_tolerance.is_finite()
            || options.intersection_tolerance <= 0.0
            || !options.angular_resolution.is_finite()
            || options.angular_resolution <= 0.0
        {
            return Err(RegionPatchError::InvalidRegion);
        }
        if self
            .evidence
            .source_planar_support
            .iter()
            .zip(&partition.region_by_polygon)
            .any(|(&support, &region)| {
                partition.region_by_support.get(support) != Some(&region)
                    || partition.carriers.get(region).is_none_or(Option::is_none)
            })
        {
            return Err(RegionPatchError::InvalidRegion);
        }
        Ok(())
    }

    fn patch_junctions(
        &self,
        partition: &SurfaceRegionPartition,
        selected: &BTreeSet<usize>,
        patch_ids: &BTreeSet<usize>,
        boundary: &BTreeSet<usize>,
        options: RegionPatchOptions,
    ) -> Result<BTreeMap<usize, Vertex>, RegionPatchError> {
        let mut incident = boundary
            .iter()
            .map(|&id| (id, BTreeSet::new()))
            .collect::<BTreeMap<_, _>>();
        // Include point-only incident regions, not just the carrier on the
        // current edge. Reconciliation is global over each source junction.
        for (polygon, &region) in partition.region_by_polygon.iter().enumerate() {
            if selected.contains(&region) {
                continue;
            }
            for id in self
                .evidence
                .source_polygon_vertices(polygon)
                .ok_or(RegionPatchError::InvalidRegion)?
            {
                if let Some(regions) = incident.get_mut(id) {
                    regions.insert(region);
                } else if patch_ids.contains(id) {
                    // A retained region touching the interior of this patch at
                    // one vertex cannot be attached using only boundary edges.
                    return Err(RegionPatchError::InvalidRegion);
                }
            }
        }
        let reconciliation = ReconciliationOptions {
            max_displacement: options.budget.vertex,
            parallel_angle_tolerance: options.angular_resolution,
            intersection_tolerance: options.intersection_tolerance,
            max_normal_angle: options.angular_resolution * 2.0,
        };
        incident
            .into_iter()
            .map(|(id, regions)| {
                let source = *self
                    .evidence
                    .positions
                    .get(id)
                    .ok_or(RegionPatchError::InvalidRegion)?;
                let carriers = regions
                    .into_iter()
                    .map(|region| partition.carriers[region])
                    .collect::<Option<Vec<_>>>()
                    .ok_or(RegionPatchError::InvalidRegion)?;
                let fine = ReconciliationOptions {
                    intersection_tolerance: junction_roundoff(source, &carriers)
                        .min(reconciliation.intersection_tolerance),
                    ..reconciliation
                };
                let point = reconcile_point(source, &carriers, fine)
                    .or_else(|| reconcile_point(source, &carriers, reconciliation))
                    .ok_or(RegionPatchError::Junction(id))?;
                Ok((id, builder::vertex(point)))
            })
            .collect()
    }

    fn lift_patch_edge(
        &self,
        carrier: SurfaceCarrier,
        ids: [usize; 2],
        front: &Vertex,
        back: &Vertex,
        options: RegionPatchOptions,
    ) -> Result<Edge, RegionPatchError> {
        lift_native_source_edge(
            Some(carrier),
            ids.map(|id| self.evidence.positions[id]),
            front,
            back,
            options,
        )
        .map(|(edge, _)| edge)
        .ok_or(RegionPatchError::Edge(ids))
    }
}

// A topologically consistent edge may still have its actual curve
// endpoint slightly outside the exact angular wedge defined by its
// shared native vertex. First seek a solution to arithmetic precision;
// retain the original admissibility tolerance if no such point exists.
pub(super) fn junction_roundoff(source: crate::Point3, carriers: &[SurfaceCarrier]) -> f64 {
    let scale = |p: crate::Point3| p.x.abs().max(p.y.abs()).max(p.z.abs());
    let maximum = carriers
        .iter()
        .fold(scale(source).max(1.0), |bound, carrier| {
            let size = match *carrier {
                SurfaceCarrier::Plane { origin, .. } => scale(origin),
                SurfaceCarrier::LinearRevolution {
                    origin,
                    radius_at_origin,
                    ..
                } => scale(origin).max(radius_at_origin.abs()),
                SurfaceCarrier::CircularRevolution {
                    axial_center,
                    radial_center,
                    profile_radius,
                    ..
                } => scale(axial_center)
                    .max(radial_center.abs())
                    .max(profile_radius.abs()),
            };
            bound.max(size)
        });
    32.0 * f64::EPSILON * maximum
}

/// Construct a single source-segment edge with one shared endpoint
/// identity at either end. The bound is measured against the original
/// source chord, not the displaced native vertex chord.
pub(super) fn lift_native_source_edge(
    carrier: Option<SurfaceCarrier>,
    source: [crate::Point3; 2],
    front: &Vertex,
    back: &Vertex,
    options: RegionPatchOptions,
) -> Option<(Edge, f64)> {
    let movement = (front.point() - source[0])
        .magnitude()
        .max((back.point() - source[1]).magnitude());
    if !movement.is_finite() || movement > options.budget.vertex || front == back {
        return None;
    }
    let (edge, bound) = if let Some(carrier) = carrier {
        let lifted = carrier.lift_source_chord_from_junctions(
            source,
            [front.point(), back.point()],
            options.budget,
            options.intersection_tolerance,
        )?;
        (
            Edge::new(front, back, Curve::NurbsCurve(lifted.curve)).ok()?,
            lifted.chord_deviation_upper_bound,
        )
    } else {
        (builder::line(front, back), movement)
    };
    (bound <= options.budget.chord && edge.is_geometric_consistent()).then_some((edge, bound))
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    use crate::{
        BoundedCurve, ParametricCurve, Point3, Vector3,
        mesh_recovery::{EvidenceOptions, analyze_mesh_evidence},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_topology::shell::ShellCondition;

    type SphericalFixture = (
        crate::mesh_recovery::MeshEvidence,
        SurfaceRegionPartition,
        Vec<usize>,
        usize,
    );

    pub(in crate::mesh_recovery::evidence::region_boundary) fn spherical_fixture()
    -> Result<SphericalFixture, Box<dyn std::error::Error>> {
        let q = std::f64::consts::FRAC_1_SQRT_2;
        let positions = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(q, q, 0.0),
            Point3::new(0.0, q, q),
            Point3::new(q, 0.0, q),
        ];
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter([
                [0, 4, 6],
                [4, 1, 5],
                [6, 5, 2],
                [4, 5, 6],
                [1, 4, 3],
                [4, 0, 3],
                [2, 5, 3],
                [5, 1, 3],
                [0, 6, 3],
                [6, 2, 3],
            ]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let mut supports = evidence.source_planar_support[..4].to_vec();
        supports.sort_unstable();
        supports.dedup();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: 0.0,
            profile_radius: 1.0,
        };
        let proposal = index
            .measure_surface_region(
                &supports,
                carrier,
                SurfaceRegionBudget {
                    vertex: 1.0e-8,
                    chord: 0.25,
                },
            )
            .ok_or("sphere")?;
        let partition = index.select_surface_partition(&[proposal], &evidence, 100)?;
        let sphere_region = partition.region_by_polygon[0];
        assert!(
            partition.region_by_polygon[..4]
                .iter()
                .all(|&id| id == sphere_region)
        );
        let mut fallback = partition.region_by_polygon[4..].to_vec();
        fallback.sort_unstable();
        fallback.dedup();
        assert!(!fallback.contains(&sphere_region));
        Ok((evidence, partition, fallback, sphere_region))
    }

    #[test]
    fn region_patch_lifts_shared_spherical_boundary() -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, partition, fallback, _) = spherical_fixture()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let options = RegionPatchOptions {
            budget: SurfaceRegionBudget {
                vertex: 1.0e-9,
                chord: 0.1,
            },
            intersection_tolerance: 1.0e-10,
            angular_resolution: 1.0e-8,
        };
        let patch = index.approximate_region_patch(&partition, &fallback, options)?;
        assert_eq!(patch.source_polygons, (4..10).collect::<Vec<_>>());
        assert_eq!(patch.boundary_edges.len(), 6);
        assert_eq!(
            patch.approximation.shell.shell_condition(),
            ShellCondition::Oriented
        );
        assert!(patch.approximation.shell.is_geometric_consistent());
        assert!(patch.approximation.deviation_upper_bound <= options.budget.chord);
        for prescribed in &patch.boundary_edges {
            let curve = prescribed.edge.oriented_curve();
            let (start, end) = curve.range_tuple();
            for i in 0..=32 {
                let t = f64::from(i) / 32.0;
                let point = curve.evaluate(t.mul_add(end - start, start));
                assert!(((point - Point3::new(0.0, 0.0, 0.0)).magnitude() - 1.0).abs() < 1.0e-9);
            }
        }
        assert!(
            index
                .approximate_region_patch(
                    &partition,
                    &fallback,
                    RegionPatchOptions {
                        budget: SurfaceRegionBudget {
                            chord: 0.02,
                            ..options.budget
                        },
                        ..options
                    }
                )
                .is_err()
        );
        Ok(())
    }
}
