//! Explore independent curved meridian proposals on the region contact graph.
//! Labels are never used for graph construction, fitting, or scoring.

use super::{PartitionCarriers, partition_contact_metrics};
use monstertruck_modeling::mesh_recovery::{
    CircularMeridianEvidence, MeshEvidence, RegionBoundaryIndex, RegionContactGraph,
    SurfaceCarrier, measure_circular_meridian,
};
use monstertruck_modeling::{InnerProduct, Point3, Vector3};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};

// Exploratory relative accuracy budget, shared across placement and profile
// measurements. Actual residuals are reported, not silently discarded.
const BUDGET: f64 = 1.0e-4;

struct CandidateIndex<'a> {
    evidence: &'a MeshEvidence,
    boundary: &'a RegionBoundaryIndex<'a>,
    carriers: &'a PartitionCarriers,
    supports: BTreeMap<usize, Vec<usize>>,
    neighbors: BTreeMap<usize, BTreeSet<usize>>,
}

struct Candidate {
    regions: Vec<usize>,
    supports: Vec<usize>,
    fit: CircularMeridianEvidence,
    boundary_loops: usize,
    source_polygons: usize,
    pinches: usize,
}

impl Candidate {
    fn accepted(&self) -> bool {
        let tolerance = BUDGET * self.fit.geometric_scale();
        self.fit.resolves_curvature_within(tolerance)
            && self.fit.normal_samples == self.source_polygons
            && self.pinches == 0
    }

    fn summary(&self, evidence: &MeshEvidence) -> Value {
        let polygons = self
            .supports
            .iter()
            .flat_map(|&id| {
                evidence.planar_supports[id]
                    .source_face_indices
                    .iter()
                    .copied()
            })
            .collect::<Vec<_>>();
        json!({
            "regions":self.regions,
            "supports":self.supports.len(),
            "source_polygons":polygons,
            "source_vertices":self.fit.source_vertices,
            "resolved_meridian_samples":self.fit.resolved_meridian_samples,
            "axial_center":[self.fit.axial_center.x,self.fit.axial_center.y,self.fit.axial_center.z],
            "axis":[self.fit.axis.x,self.fit.axis.y,self.fit.axis.z],
            "radial_center":self.fit.radial_center,
            "profile_radius":self.fit.profile_radius,
            "observed_meridian_span":self.fit.observed_meridian_span,
            "linear_profile_residual_lower_bound":self.fit.linear_profile_residual_lower_bound,
            "relative_max_vertex_residual":self.fit.max_vertex_residual/self.fit.profile_radius,
            "relative_rms_vertex_residual":self.fit.rms_vertex_residual/self.fit.profile_radius,
            "max_normal_angle_radians":self.fit.max_normal_angle_radians,
            "mixed_normal_orientation":self.fit.mixed_normal_orientation,
            "normal_samples":self.fit.normal_samples,
            "boundary_loops":self.boundary_loops,
            "pinches":self.pinches,
            "accepted":self.accepted(),
        })
    }
}

impl<'a> CandidateIndex<'a> {
    fn new(
        evidence: &'a MeshEvidence,
        boundary: &'a RegionBoundaryIndex<'a>,
        carriers: &'a PartitionCarriers,
        graph: &RegionContactGraph,
        region_by_support: &[usize],
    ) -> Self {
        let mut supports = BTreeMap::<usize, Vec<usize>>::new();
        for (id, &region) in region_by_support.iter().enumerate() {
            if matches!(
                carriers.carriers[region],
                Some(SurfaceCarrier::LinearRevolution { .. })
            ) {
                supports.entry(region).or_default().push(id);
            }
        }
        let mut neighbors = BTreeMap::<usize, BTreeSet<usize>>::new();
        for path in &graph.paths {
            let (left, Some(right)) = path.regions else {
                continue;
            };
            if !supports.contains_key(&left) || !supports.contains_key(&right) {
                continue;
            }
            let Some(reference) = path
                .edges
                .first()
                .map(|edge| evidence.positions[edge.vertices[0]])
            else {
                continue;
            };
            let Some((origin, axis)) = axis_of(carriers, left) else {
                continue;
            };
            let radius = (reference - origin).cross(&axis).norm();
            if radius <= 0.0 {
                continue;
            }
            let (Some(first), Some(second)) = (carriers.carriers[left], carriers.carriers[right])
            else {
                continue;
            };
            let Some(fit) = first.compare_in_ball(second, reference, 0.0) else {
                continue;
            };
            if fit.direction_angle_radians > BUDGET
                || fit.reference_position_difference > BUDGET * radius
                || fit.reference_profile_difference > BUDGET * radius
            {
                continue;
            }
            neighbors.entry(left).or_default().insert(right);
            neighbors.entry(right).or_default().insert(left);
        }
        Self {
            evidence,
            boundary,
            carriers,
            supports,
            neighbors,
        }
    }

    fn measure(&self, regions: Vec<usize>) -> Option<Candidate> {
        let seed = *regions.first()?;
        let (origin, axis) = axis_of(self.carriers, seed)?;
        let extent = regions
            .iter()
            .filter_map(|&id| self.carriers.extents[id])
            .map(|bounds| bounds.radius_from(origin))
            .fold(0.0_f64, f64::max);
        let mut supports = regions
            .iter()
            .flat_map(|id| self.supports[id].iter().copied())
            .collect::<Vec<_>>();
        supports.sort_unstable();
        // Resolve source meridian stations at the observed angular
        // uncertainty of the input carriers, plus floating-point noise.
        let angular = regions
            .iter()
            .filter_map(|&id| axis_of(self.carriers, id))
            .map(|(_, other)| axis.cross(&other).norm())
            .fold(0.0_f64, f64::max);
        let resolution = extent * (angular + 1.0e-10);
        let fit = measure_circular_meridian(self.evidence, &supports, origin, axis, resolution)?;
        let boundary = self.boundary.boundary(&supports).ok()?;
        Some(Candidate {
            regions,
            supports,
            fit,
            boundary_loops: boundary.loops.len(),
            source_polygons: boundary.source_polygon_count,
            pinches: boundary.pinched_vertices.len(),
        })
    }

    fn seeds(&self) -> BTreeSet<Vec<usize>> {
        let mut triples = BTreeSet::new();
        for (&center, neighbors) in &self.neighbors {
            for &left in neighbors {
                for &right in
                    neighbors.range((std::ops::Bound::Excluded(left), std::ops::Bound::Unbounded))
                {
                    let mut ids = vec![left, center, right];
                    ids.sort_unstable();
                    triples.insert(ids);
                }
            }
        }
        triples
    }

    fn grow(&self, mut candidate: Candidate) -> Candidate {
        loop {
            let frontier = candidate
                .regions
                .iter()
                .filter_map(|id| self.neighbors.get(id))
                .flatten()
                .copied()
                .filter(|id| candidate.regions.binary_search(id).is_err())
                .collect::<BTreeSet<_>>();
            let best = frontier
                .into_iter()
                .filter_map(|id| {
                    let mut regions = candidate.regions.clone();
                    let insertion = regions.partition_point(|other| *other < id);
                    regions.insert(insertion, id);
                    self.measure(regions)
                })
                .filter(Candidate::accepted)
                .max_by(|a, b| {
                    a.source_polygons.cmp(&b.source_polygons).then_with(|| {
                        b.fit
                            .max_vertex_residual
                            .total_cmp(&a.fit.max_vertex_residual)
                    })
                });
            let Some(next) = best else {
                return candidate;
            };
            candidate = next;
        }
    }
}

fn axis_of(carriers: &PartitionCarriers, id: usize) -> Option<(Point3, Vector3)> {
    match carriers.carriers[id]? {
        SurfaceCarrier::LinearRevolution { origin, axis, .. } => Some((origin, axis)),
        SurfaceCarrier::Plane { .. } | SurfaceCarrier::CircularRevolution { .. } => None,
    }
}

/// Compare a deterministic, nonoverlapping subset of the accepted proposals
/// against the input partition. All alternatives remain in the report.
fn applied_partition(
    evidence: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    carriers: &PartitionCarriers,
    graph: &RegionContactGraph,
    region_by_support: &[usize],
    candidates: &BTreeMap<Vec<usize>, Candidate>,
) -> Result<Value, String> {
    let mut candidates = candidates.values().collect::<Vec<_>>();
    candidates.sort_unstable_by(|a, b| {
        b.source_polygons
            .cmp(&a.source_polygons)
            .then_with(|| {
                a.fit
                    .max_vertex_residual
                    .total_cmp(&b.fit.max_vertex_residual)
            })
            .then_with(|| a.regions.cmp(&b.regions))
    });
    let mut selected = BTreeSet::new();
    let mut region_map = region_by_support.to_vec();
    let mut merged_carriers = carriers.carriers.clone();
    let mut merged_regions = 0;
    for candidate in candidates {
        if candidate.regions.iter().any(|id| selected.contains(id)) {
            continue;
        }
        let Some(&region_id) = candidate.regions.first() else {
            continue;
        };
        for &support in &candidate.supports {
            region_map[support] = region_id;
        }
        selected.extend(candidate.regions.iter().copied());
        merged_carriers[region_id] = Some(SurfaceCarrier::from_circular_meridian(&candidate.fit));
        merged_regions += 1;
    }
    let region_by_polygon = evidence
        .source_planar_support
        .iter()
        .map(|&support| region_map[support])
        .collect::<Vec<_>>();
    let mut eliminated_paths = 0;
    let mut eliminated_edges = 0;
    for path in &graph.paths {
        let Some(edge) = path.edges.first() else {
            continue;
        };
        if let Some(other) = edge.exterior_polygon
            && region_by_polygon[edge.interior_polygon] == region_by_polygon[other]
        {
            eliminated_paths += 1;
            eliminated_edges += path.edges.len();
        }
    }
    let merged = boundary
        .partition_contacts(&region_by_polygon)
        .map_err(|error| error.to_string())?;
    let carriers = PartitionCarriers::new(evidence, merged_carriers, &region_by_polygon);
    Ok(json!({
        "selected_input_regions":selected.len(),
        "merged_regions":merged_regions,
        "eliminated_contact_paths":eliminated_paths,
        "eliminated_source_edges":eliminated_edges,
        "contact_graph":partition_contact_metrics(&merged,evidence,Some(&carriers)),
    }))
}

pub fn summarize(
    evidence: &MeshEvidence,
    boundary: &RegionBoundaryIndex<'_>,
    carriers: &PartitionCarriers,
    graph: &RegionContactGraph,
    region_by_support: &[usize],
) -> Result<Value, String> {
    let index = CandidateIndex::new(evidence, boundary, carriers, graph, region_by_support);
    let seeds = index.seeds();
    let mut seed_results = Vec::new();
    let mut grown = BTreeMap::new();
    for seed in &seeds {
        let Some(candidate) = index.measure(seed.clone()) else {
            continue;
        };
        seed_results.push(candidate.summary(evidence));
        if candidate.accepted() {
            let candidate = index.grow(candidate);
            grown.entry(candidate.regions.clone()).or_insert(candidate);
        }
    }
    let applied = applied_partition(
        evidence,
        boundary,
        carriers,
        graph,
        region_by_support,
        &grown,
    )?;
    Ok(json!({
        "relative_budget":BUDGET,
        "connected_coaxial_region_pairs":index.neighbors.values().map(BTreeSet::len).sum::<usize>()/2,
        "connected_three_region_seeds":seeds.len(),
        "measured_seed_count":seed_results.len(),
        "seed_candidates":seed_results,
        "grown_candidates":grown.values().map(|c|c.summary(evidence)).collect::<Vec<_>>(),
        "applied_partition":applied,
    }))
}
