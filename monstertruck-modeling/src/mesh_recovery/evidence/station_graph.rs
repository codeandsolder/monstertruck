//! Axis-conditioned source evidence, independent of the local strip classifier.
//! Repeated observed meridian stations compress redundant azimuth samples;
//! source support connectivity is retained separately from the projection.

use super::super::sample_fit::fit_circle_2d;
use super::meridian::{MeridianCluster, circular_profile_normal, cluster_meridian};
use super::{
    CarrierResidual, CircularMeridianEvidence, MeshEvidence, RegionBoundaryIndex, SurfaceCarrier,
    measure_circular_meridian,
};
use crate::{InnerSpace, Point3, Vector3};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// Projected sample and fitting budgets, in original mesh units.
#[derive(Clone, Copy, Debug)]
pub struct MeridianGraphOptions {
    /// Coordinate resolution for seed-anchored station clustering.
    pub station_resolution: f64,
    /// Absolute source-vertex error allowed when testing a profile.
    pub vertex_tolerance: f64,
    /// Relative source error bound, also limited by observed meridian extent.
    pub relative_tolerance: f64,
    /// Source-normal comparison uncertainty in radians.
    pub normal_tolerance_radians: f64,
}

impl MeridianGraphOptions {
    fn valid(self) -> bool {
        self.station_resolution.is_finite()
            && self.station_resolution > 0.0
            && self.vertex_tolerance.is_finite()
            && self.vertex_tolerance > 0.0
            && self.relative_tolerance.is_finite()
            && self.relative_tolerance > 0.0
            && (0.0..std::f64::consts::FRAC_PI_2).contains(&self.normal_tolerance_radians)
    }
}

/// A repeated axial/radial observation under an independently placed axis.
#[derive(Clone, Debug)]
pub struct MeridianStation {
    /// Mean axial coordinate relative to the index axis origin.
    pub axial: f64,
    /// Mean distance from the index axis.
    pub radius: f64,
    /// Canonical original source-position IDs, never synthetic points.
    pub vertices: Vec<usize>,
    /// Observed nonparallel azimuth displacement, in length units.
    pub transverse_span: f64,
}

/// Original supports independently sampling the same meridian interval.
#[derive(Clone, Debug)]
pub struct MeridianTransition {
    /// Endpoints in canonical station order.
    pub stations: [usize; 2],
    /// Entire original planar supports restricted to these two stations.
    pub supports: Vec<usize>,
    /// Distinct support normals span resolved directions around the axis.
    pub repeated: bool,
}

/// A complete, source-connected circular meridian proposal.
#[derive(Clone, Debug)]
pub struct StationMeridianCandidate {
    /// Sorted original support membership, never projected connectivity alone.
    pub supports: Vec<usize>,
    /// Number of independently repeated stations in its source-connected core.
    pub repeated_stations: usize,
    /// Complete original vertex and oriented normal measurements.
    pub fit: CircularMeridianEvidence,
    /// Continuous approximation bounds on every original source chord.
    pub source_chords: CarrierResidual,
    /// Number of original polygons in the candidate.
    pub source_polygons: usize,
    /// Exact source boundary loop count.
    pub boundary_loops: usize,
}

/// Graph and result counts preserve the distinction between proposals and
/// an eventual globally selected B-rep.
#[derive(Debug)]
pub struct StationMeridianDiscovery {
    /// All connected three-station seed hypotheses.
    pub seed_count: usize,
    /// Distinct connected projected profiles with independent fourth stations.
    pub profile_count: usize,
    /// Complete, source-connected candidates, including competing alternatives.
    pub candidates: Vec<StationMeridianCandidate>,
}

struct SupportSamples {
    vertices: Vec<usize>,
    stations: Vec<usize>,
}

/// Observed meridian station graph at a fixed independently fitted axis.
///
/// A repeated station transition can seed a curved surface without first
/// identifying a conical strip. Projected adjacency is never used as actual
/// three-dimensional source connectivity.
pub struct MeridianStationGraph<'a> {
    evidence: &'a MeshEvidence,
    origin: Point3,
    axis: Vector3,
    options: MeridianGraphOptions,
    samples: Vec<(f64, f64)>,
    stations: Vec<MeridianStation>,
    transitions: Vec<MeridianTransition>,
    supports: Vec<SupportSamples>,
    supports_by_station: Vec<Vec<usize>>,
    adjacency: Vec<Vec<usize>>,
}

impl<'a> MeridianStationGraph<'a> {
    /// Project every welded source vertex once and index whole planar supports.
    /// Returns None for invalid axis, tolerances, or source references.
    #[must_use = "the graph retains original source provenance"]
    pub fn new(
        evidence: &'a MeshEvidence,
        origin: Point3,
        axis: Vector3,
        options: MeridianGraphOptions,
    ) -> Option<Self> {
        if !options.valid()
            || !(0.0..=64.0 * f64::EPSILON).contains(&(axis.magnitude2() - 1.0).abs())
        {
            return None;
        }
        let samples = evidence
            .positions
            .iter()
            .map(|&point| {
                let delta = point - origin;
                (delta.dot(axis), delta.cross(axis).magnitude())
            })
            .collect::<Vec<_>>();
        if samples
            .iter()
            .any(|&(z, r)| !z.is_finite() || !r.is_finite())
        {
            return None;
        }
        let clusters = cluster_meridian(&samples, options.station_resolution);
        let (stations, station_by_vertex) = index_stations(evidence, origin, axis, clusters);
        let (supports, supports_by_station) =
            index_supports(evidence, &station_by_vertex, stations.len())?;
        let transitions = build_transitions(evidence, &stations, &supports, options);
        let mut adjacency = vec![Vec::new(); stations.len()];
        for transition in transitions.iter().filter(|edge| edge.repeated) {
            let [a, b] = transition.stations;
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        Some(Self {
            evidence,
            origin,
            axis,
            options,
            samples,
            stations,
            transitions,
            supports,
            supports_by_station,
            adjacency,
        })
    }

    /// All source stations, including unresolved singleton observations.
    #[must_use = "stations include their original source vertices"]
    pub fn stations(&self) -> &[MeridianStation] { &self.stations }

    /// Observed two-station intervals, including intervals without redundancy.
    #[must_use = "transition evidence must be retained"]
    pub fn transitions(&self) -> &[MeridianTransition] { &self.transitions }

    /// Discover circular profile candidates without any native CAD labels.
    ///
    /// Repeated transitions seed connected station triples; a fourth station
    /// independently tests the profile. Each source region must have a
    /// connected core spanning four resolved repeated stations and is then
    /// extended to any geometrically compatible whole support.
    #[must_use = "candidates must be checked before changing topology"]
    pub fn discover_circular(
        &self,
        boundary: &RegionBoundaryIndex<'_>,
    ) -> StationMeridianDiscovery {
        let seeds = self.seeds();
        let mut profiles = BTreeSet::new();
        let mut candidates = BTreeMap::new();
        for seed in &seeds {
            let Some(profile) = self.fit_station_profile(seed) else {
                continue;
            };
            let stations = self.connected_profile(*seed, profile);
            if stations.len() < 4 || !profiles.insert(stations.clone()) {
                continue;
            }
            let Some(profile) = self.fit_station_profile(&stations) else {
                continue;
            };
            if stations.iter().any(|&id| {
                let station = &self.stations[id];
                profile.error((station.axial, station.radius)) > self.options.vertex_tolerance
            }) {
                continue;
            }
            for candidate in self.source_candidates(&stations, profile, boundary) {
                candidates
                    .entry(candidate.supports.clone())
                    .or_insert(candidate);
            }
        }
        let mut candidates = candidates.into_values().collect::<Vec<_>>();
        candidates.sort_unstable_by(|a, b| {
            b.source_polygons
                .cmp(&a.source_polygons)
                .then_with(|| {
                    a.fit
                        .max_vertex_residual
                        .total_cmp(&b.fit.max_vertex_residual)
                })
                .then_with(|| a.supports.cmp(&b.supports))
        });
        StationMeridianDiscovery {
            seed_count: seeds.len(),
            profile_count: profiles.len(),
            candidates,
        }
    }

    fn seeds(&self) -> BTreeSet<[usize; 3]> {
        let mut seeds = BTreeSet::new();
        for (center, neighbors) in self.adjacency.iter().enumerate() {
            for (i, &left) in neighbors.iter().enumerate() {
                for &right in &neighbors[i + 1..] {
                    let mut ids = [left, center, right];
                    ids.sort_unstable();
                    seeds.insert(ids);
                }
            }
        }
        seeds
    }

    fn fit_station_profile(&self, ids: &[usize]) -> Option<Profile> {
        let first = &self.stations[*ids.first()?];
        let local = ids
            .iter()
            .map(|&id| {
                let station = &self.stations[id];
                (station.axial - first.axial, station.radius - first.radius)
            })
            .collect::<Vec<_>>();
        let ([z, r], radius) = fit_circle_2d(&local)?;
        Some(Profile {
            z: z + first.axial,
            r: r + first.radius,
            radius,
        })
    }

    fn connected_profile(&self, seed: [usize; 3], profile: Profile) -> Vec<usize> {
        let mut visited = BTreeSet::from(seed);
        let mut queue = VecDeque::from(seed);
        while let Some(id) = queue.pop_front() {
            for &neighbor in &self.adjacency[id] {
                let station = &self.stations[neighbor];
                if profile.error((station.axial, station.radius)) <= self.options.vertex_tolerance
                    && visited.insert(neighbor)
                {
                    queue.push_back(neighbor);
                }
            }
        }
        visited.into_iter().collect()
    }

    fn source_candidates(
        &self,
        stations: &[usize],
        profile: Profile,
        boundary: &RegionBoundaryIndex<'_>,
    ) -> Vec<StationMeridianCandidate> {
        let source = stations
            .iter()
            .flat_map(|&id| self.supports_by_station[id].iter().copied())
            .collect::<BTreeSet<_>>();
        let mut eligible = BTreeMap::new();
        for id in source {
            if let Some(orientation) = self.support_orientation(id, profile) {
                eligible.insert(id, orientation);
            }
        }
        // First find source-connected cores, not paths joined only in
        // projected space. A core consists of repeated station intervals.
        let mut cores = eligible
            .iter()
            .filter(|&(&id, orientation)| {
                let ids = &self.supports[id].stations;
                orientation.within_sample_span
                    && ids.len() == 2
                    && ids.iter().all(|id| stations.binary_search(id).is_ok())
                    && self.adjacency[ids[0]].contains(&ids[1])
            })
            .map(|(&id, orientation)| (id, orientation.positive))
            .collect::<BTreeMap<_, _>>();
        let mut candidates = Vec::new();
        let mut completed = BTreeSet::new();
        while let Some((&seed, &sign)) = cores.first_key_value() {
            let core = self.component(seed, sign, &mut cores);
            let repeated = self.local_repeated_station_count(&core);
            if repeated < 4 || completed.contains(&seed) {
                continue;
            }
            let supports = self.grow_component(seed, sign, profile, &mut eligible);
            completed.extend(supports.iter().copied());
            // A broad candidate may acquire a pinched boundary, so keep the
            // independently supported core available as an alternative.
            if let Some(candidate) = self.validate(supports, boundary) {
                candidates.push(candidate);
            } else if let Some(candidate) = self.validate(core, boundary) {
                candidates.push(candidate);
            }
        }
        candidates
    }

    // Repeated transitions elsewhere on the same placed axis are not
    // independent support for this source-connected region.
    fn local_repeated_station_count(&self, core: &[usize]) -> usize {
        let mut pairs = BTreeMap::<[usize; 2], Vec<usize>>::new();
        for &id in core {
            if let [a, b] = self.supports[id].stations.as_slice() {
                pairs.entry([*a, *b]).or_default().push(id);
            }
        }
        let mut adjacency = BTreeMap::<usize, Vec<usize>>::new();
        for ([a, b], supports) in pairs {
            if distinct_normals(
                self.evidence,
                &supports,
                self.options.normal_tolerance_radians,
            ) {
                adjacency.entry(a).or_default().push(b);
                adjacency.entry(b).or_default().push(a);
            }
        }
        let mut seen = BTreeSet::new();
        let mut largest = 0;
        for &seed in adjacency.keys() {
            if !seen.insert(seed) {
                continue;
            }
            let mut queue = VecDeque::from([seed]);
            let mut count = 0;
            while let Some(id) = queue.pop_front() {
                count += 1;
                for &other in &adjacency[&id] {
                    if seen.insert(other) {
                        queue.push_back(other);
                    }
                }
            }
            largest = largest.max(count);
        }
        largest
    }

    // Check whole supports lazily across the actual source frontier, even
    // if their stations were not repeated in the projected seed graph.
    fn grow_component(
        &self,
        seed: usize,
        sign: bool,
        profile: Profile,
        eligible: &mut BTreeMap<usize, OrientationEvidence>,
    ) -> Vec<usize> {
        let mut visited = BTreeSet::from([seed]);
        let mut queue = VecDeque::from([seed]);
        let mut selected = vec![seed];
        while let Some(id) = queue.pop_front() {
            for &interface in &self.evidence.support_interface_indices[id] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let other = if a == id { b } else { a };
                if !visited.insert(other) {
                    continue;
                }
                let orientation = eligible
                    .get(&other)
                    .copied()
                    .or_else(|| self.support_orientation(other, profile));
                if let Some(orientation) = orientation {
                    eligible.insert(other, orientation);
                    if orientation.positive == sign {
                        selected.push(other);
                        queue.push_back(other);
                    }
                }
            }
        }
        selected.sort_unstable();
        selected
    }

    fn component(
        &self,
        seed: usize,
        sign: bool,
        remaining: &mut BTreeMap<usize, bool>,
    ) -> Vec<usize> {
        remaining.remove(&seed);
        let mut queue = VecDeque::from([seed]);
        let mut selected = vec![seed];
        while let Some(id) = queue.pop_front() {
            for &interface in &self.evidence.support_interface_indices[id] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let neighbor = if a == id { b } else { a };
                if remaining.get(&neighbor) == Some(&sign) {
                    remaining.remove(&neighbor);
                    selected.push(neighbor);
                    queue.push_back(neighbor);
                }
            }
        }
        selected.sort_unstable();
        selected
    }

    fn support_orientation(&self, id: usize, profile: Profile) -> Option<OrientationEvidence> {
        let support = &self.evidence.planar_supports[id];
        if self.supports[id].stations.len() < 2
            || self.supports[id]
                .vertices
                .iter()
                .any(|&id| profile.error(self.samples[id]) > self.options.vertex_tolerance)
        {
            return None;
        }
        let center = self.origin + self.axis * profile.z;
        let mut orientation: Option<OrientationEvidence> = None;
        for &polygon in &support.source_face_indices {
            let vertices = self.evidence.source_polygon_vertices(polygon)?;
            let measured = polygon_orientation(
                self.evidence,
                vertices,
                support.normal,
                center,
                self.axis,
                profile.r,
                self.options.normal_tolerance_radians,
            )?;
            if let Some(previous) = &mut orientation {
                if previous.positive != measured.positive {
                    return None;
                }
                previous.within_sample_span &= measured.within_sample_span;
            } else {
                orientation = Some(measured);
            }
        }
        orientation
    }

    fn validate(
        &self,
        supports: Vec<usize>,
        boundary: &RegionBoundaryIndex<'_>,
    ) -> Option<StationMeridianCandidate> {
        let fit = measure_circular_meridian(
            self.evidence,
            &supports,
            self.origin,
            self.axis,
            self.options.station_resolution,
        )?;
        let tolerance = self
            .options
            .vertex_tolerance
            .min(self.options.relative_tolerance * fit.geometric_scale());
        if !fit.resolves_curvature_within(tolerance) {
            return None;
        }
        let profile = Profile {
            z: (fit.axial_center - self.origin).dot(self.axis),
            r: fit.radial_center,
            radius: fit.profile_radius,
        };
        let mut orientation = None;
        let mut cores = BTreeMap::new();
        for &id in &supports {
            let measured = self.support_orientation(id, profile)?;
            if orientation.is_some_and(|other| other != measured.positive) {
                return None;
            }
            orientation = Some(measured.positive);
            let ids = &self.supports[id].stations;
            if measured.within_sample_span
                && ids.len() == 2
                && self.adjacency[ids[0]].contains(&ids[1])
            {
                cores.insert(id, measured.positive);
            }
        }
        let mut repeated_stations = 0;
        while let Some((&seed, &sign)) = cores.first_key_value() {
            let core = self.component(seed, sign, &mut cores);
            repeated_stations = repeated_stations.max(self.local_repeated_station_count(&core));
        }
        if repeated_stations < 4 {
            return None;
        }
        let region = boundary.boundary(&supports).ok()?;
        if !region.pinched_vertices.is_empty() || fit.normal_samples != region.source_polygon_count
        {
            return None;
        }
        let source_chords = SurfaceCarrier::from_circular_meridian(&fit)
            .measure_source_chords(self.evidence, &supports)?;
        Some(StationMeridianCandidate {
            supports,
            source_chords,
            repeated_stations,
            fit,
            source_polygons: region.source_polygon_count,
            boundary_loops: region.loops.len(),
        })
    }
}

#[derive(Clone, Copy)]
struct Profile {
    z: f64,
    r: f64,
    radius: f64,
}

impl Profile {
    fn error(self, (z, r): (f64, f64)) -> f64 {
        ((z - self.z).hypot(r - self.r) - self.radius).abs()
    }
}

fn index_stations(
    evidence: &MeshEvidence,
    origin: Point3,
    axis: Vector3,
    clusters: Vec<MeridianCluster>,
) -> (Vec<MeridianStation>, Vec<usize>) {
    let mut by_vertex = vec![0; evidence.positions.len()];
    let stations = clusters
        .into_iter()
        .enumerate()
        .map(|(id, cluster)| {
            let first = evidence.positions[cluster.members[0]] - origin;
            let transverse = first - axis * first.dot(axis);
            let length = transverse.magnitude();
            let mut span = 0.0_f64;
            for &vertex in &cluster.members {
                by_vertex[vertex] = id;
                let delta = evidence.positions[vertex] - origin;
                if length > 0.0 {
                    span = span.max(
                        (delta - axis * delta.dot(axis))
                            .cross(transverse / length)
                            .magnitude(),
                    );
                }
            }
            MeridianStation {
                axial: cluster.center.0,
                radius: cluster.center.1,
                vertices: cluster.members,
                transverse_span: span,
            }
        })
        .collect();
    (stations, by_vertex)
}

fn index_supports(
    evidence: &MeshEvidence,
    by_vertex: &[usize],
    station_count: usize,
) -> Option<(Vec<SupportSamples>, Vec<Vec<usize>>)> {
    let mut by_station = vec![Vec::new(); station_count];
    let mut supports = Vec::with_capacity(evidence.planar_supports.len());
    for (id, support) in evidence.planar_supports.iter().enumerate() {
        let mut vertices = Vec::new();
        for &polygon in &support.source_face_indices {
            vertices.extend_from_slice(evidence.source_polygon_vertices(polygon)?);
        }
        vertices.sort_unstable();
        vertices.dedup();
        let mut stations = vertices.iter().map(|&v| by_vertex[v]).collect::<Vec<_>>();
        stations.sort_unstable();
        stations.dedup();
        for &station in &stations {
            by_station[station].push(id);
        }
        supports.push(SupportSamples { vertices, stations });
    }
    Some((supports, by_station))
}

fn build_transitions(
    evidence: &MeshEvidence,
    stations: &[MeridianStation],
    supports: &[SupportSamples],
    options: MeridianGraphOptions,
) -> Vec<MeridianTransition> {
    let mut pairs = BTreeMap::<[usize; 2], Vec<usize>>::new();
    for (id, support) in supports.iter().enumerate() {
        if let [a, b] = support.stations.as_slice() {
            pairs.entry([*a, *b]).or_default().push(id);
        }
    }
    pairs
        .into_iter()
        .map(|(ids, supports)| {
            let rings = ids.iter().all(|&id| {
                stations[id].vertices.len() >= 3
                    && stations[id].transverse_span > options.station_resolution
            });
            let repeated =
                rings && distinct_normals(evidence, &supports, options.normal_tolerance_radians);
            MeridianTransition {
                stations: ids,
                supports,
                repeated,
            }
        })
        .collect()
}

fn distinct_normals(evidence: &MeshEvidence, supports: &[usize], angle_tolerance: f64) -> bool {
    let Some(&first) = supports.first() else {
        return false;
    };
    let first = evidence.planar_supports[first].normal;
    supports.iter().any(|&id| {
        first.cross(evidence.planar_supports[id].normal).magnitude() > angle_tolerance.sin()
    })
}

#[derive(Clone, Copy)]
struct OrientationEvidence {
    positive: bool,
    within_sample_span: bool,
}

fn polygon_orientation(
    evidence: &MeshEvidence,
    vertices: &[usize],
    normal: Vector3,
    center: Point3,
    axis: Vector3,
    radial_center: f64,
    angle_tolerance: f64,
) -> Option<OrientationEvidence> {
    let mut centroid = Vector3::new(0.0, 0.0, 0.0);
    for &id in vertices {
        centroid += evidence.positions[id] - center;
    }
    centroid /= crate::usize_to_f64(vertices.len());
    let expected = circular_profile_normal(centroid, axis, radial_center)?;
    let dot = expected.dot(normal);
    if dot.abs() <= angle_tolerance.sin() {
        return None;
    }
    let deviation = expected.cross(normal).magnitude().atan2(dot.abs());
    let mut span = 0.0_f64;
    for &id in vertices {
        if let Some(vertex_normal) =
            circular_profile_normal(evidence.positions[id] - center, axis, radial_center)
        {
            span = span.max(
                expected
                    .cross(vertex_normal)
                    .magnitude()
                    .atan2(expected.dot(vertex_normal)),
            );
        }
    }
    Some(OrientationEvidence {
        positive: dot > 0.0,
        within_sample_span: deviation <= span + angle_tolerance,
    })
}
