//! Axis-conditioned source evidence, independent of the local strip classifier.
//! Repeated observed meridian stations compress redundant azimuth samples;
//! source support connectivity is retained separately from the projection.

use super::super::sample_fit::fit_circle_2d;
use super::meridian::{MeridianClusters, circular_profile_normal, cluster_meridian};
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
    /// Maximum continuous chord residual relative to the local radial
    /// scale of an affine meridian. Sparse faceted polygons may interpolate
    /// the same circle; a local approximation envelope keeps those
    /// alternatives distinguishable from well-sampled surfaces.
    pub max_linear_relative_chord_error: f64,
    /// Maximum continuous chord error relative to the observed circular
    /// meridian scale. A separate explicit sampling envelope is needed
    /// because exact vertex interpolation cannot identify a design surface.
    pub max_circular_relative_chord_error: f64,
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
            && self.max_linear_relative_chord_error.is_finite()
            && self.max_linear_relative_chord_error > 0.0
            && self.max_circular_relative_chord_error.is_finite()
            && self.max_circular_relative_chord_error > 0.0
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
    /// Both endpoints have resolved azimuth sampling; even a single planar
    /// tangent strip can link independently repeated adjacent intervals.
    pub sampled: bool,
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
    /// Absolute continuous chord approximation budget.
    pub chord_tolerance: f64,
    /// Number of original polygons in the candidate.
    pub source_polygons: usize,
    /// Exact source boundary loop count.
    pub boundary_loops: usize,
}

/// Affine meridian candidate at an independently placed axis.
#[derive(Clone, Debug)]
pub struct StationLinearCandidate {
    /// Original whole support membership.
    pub supports: Vec<usize>,
    /// Fitted affine radial profile.
    pub carrier: SurfaceCarrier,
    /// Independently repeated stations inside the actual source core.
    pub repeated_stations: usize,
    /// Original source polygon count.
    pub source_polygons: usize,
    /// Maximum admitted source vertex error.
    pub vertex_tolerance: f64,
    /// Measured maximum source vertex error.
    pub max_vertex_residual: f64,
    /// Radial scale of the observed source domain, for chord approximation.
    pub chord_tolerance: f64,
    /// Continuous source chord measurement.
    pub source_chords: CarrierResidual,
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
    station_by_vertex: Vec<usize>,
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
        Self::build(evidence, origin, axis, options, false)?.1
    }

    // Discovery needs at least two independently sampled stations before
    // either a transition or a mixed-support refinement can exist.
    pub(in crate::mesh_recovery) fn for_discovery(
        evidence: &'a MeshEvidence,
        origin: Point3,
        axis: Vector3,
        options: MeridianGraphOptions,
    ) -> Option<(usize, Option<Self>)> {
        Self::build(evidence, origin, axis, options, true)
    }

    fn build(
        evidence: &'a MeshEvidence,
        origin: Point3,
        axis: Vector3,
        options: MeridianGraphOptions,
        skip_unsampled: bool,
    ) -> Option<(usize, Option<Self>)> {
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
        let count = clusters.clusters.len();
        if skip_unsampled
            && clusters
                .clusters
                .iter()
                .filter(|c| c.members.len() >= 3)
                .take(2)
                .count()
                < 2
        {
            return Some((count, None));
        }
        let (stations, station_by_vertex) = index_stations(evidence, origin, axis, clusters);
        if skip_unsampled
            && stations
                .iter()
                .filter(|station| station_is_sampled(station, options.station_resolution))
                .take(2)
                .count()
                < 2
        {
            return Some((count, None));
        }
        let (supports, supports_by_station) =
            index_supports(evidence, &station_by_vertex, stations.len())?;
        let transitions = build_transitions(evidence, &stations, &supports, options);
        let mut adjacency = vec![Vec::new(); stations.len()];
        for transition in transitions.iter().filter(|edge| edge.sampled) {
            let [a, b] = transition.stations;
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        Some((
            count,
            Some(Self {
                evidence,
                origin,
                axis,
                options,
                samples,
                stations,
                station_by_vertex,
                transitions,
                supports,
                supports_by_station,
                adjacency,
            }),
        ))
    }

    /// All source stations, including unresolved singleton observations.
    #[must_use = "stations include their original source vertices"]
    pub fn stations(&self) -> &[MeridianStation] { &self.stations }

    /// Observed two-station intervals, including intervals without redundancy.
    #[must_use = "transition evidence must be retained"]
    pub fn transitions(&self) -> &[MeridianTransition] { &self.transitions }

    /// Detect original planar groups which hide independently sampled
    /// rotational intervals. This is a refinement request, not a surface
    /// assignment. Only split a mixed group when its original polygons
    /// supply a repeated transition or link independently repeated stations.
    #[must_use = "a refinement preserves source polygons and planar alternatives"]
    pub fn mixed_support_refinements(&self) -> Vec<usize> {
        let mut pairs = BTreeMap::<[usize; 2], Vec<(usize, usize)>>::new();
        let repeated_stations = self
            .transitions
            .iter()
            .filter(|edge| edge.repeated)
            .flat_map(|edge| edge.stations)
            .collect::<BTreeSet<_>>();
        for (id, support) in self.supports.iter().enumerate() {
            if support.stations.len() <= 2
                || self.evidence.planar_supports[id].source_face_indices.len() < 2
            {
                continue;
            }
            for &polygon in &self.evidence.planar_supports[id].source_face_indices {
                let Some(vertices) = self.evidence.source_polygon_vertices(polygon) else {
                    continue;
                };
                let mut stations = vertices
                    .iter()
                    .map(|&vertex| self.station_by_vertex[vertex])
                    .collect::<Vec<_>>();
                stations.sort_unstable();
                stations.dedup();
                if let [a, b] = stations.as_slice()
                    && station_is_sampled(&self.stations[*a], self.options.station_resolution)
                    && station_is_sampled(&self.stations[*b], self.options.station_resolution)
                {
                    pairs.entry([*a, *b]).or_default().push((id, polygon));
                }
            }
        }
        let mut split = BTreeSet::new();
        for (pair, members) in pairs {
            let mut parents = members
                .iter()
                .map(|&(parent, _)| parent)
                .collect::<Vec<_>>();
            parents.sort_unstable();
            parents.dedup();
            if pair.iter().all(|id| repeated_stations.contains(id))
                && self.bridge_has_independent_profile(pair)
            {
                split.extend(parents);
                continue;
            }
            let mut observed = parents.clone();
            if let Ok(id) = self
                .transitions
                .binary_search_by_key(&pair, |edge| edge.stations)
            {
                observed.extend(self.transitions[id].supports.iter().copied());
            }
            if !distinct_normals(
                self.evidence,
                &observed,
                self.options.normal_tolerance_radians,
            ) {
                continue;
            }
            let [a, b] = pair;
            let first = &self.stations[a];
            let second = &self.stations[b];
            let span = second.axial - first.axial;
            if span.abs() <= self.options.station_resolution {
                continue;
            }
            let slope = (second.radius - first.radius) / span;
            let carrier = SurfaceCarrier::LinearRevolution {
                origin: self.origin,
                axis: self.axis,
                radius_at_origin: slope.mul_add(-first.axial, first.radius),
                slope,
            };
            let tolerance = first.radius.min(second.radius).mul_add(
                self.options.max_linear_relative_chord_error,
                self.options.vertex_tolerance,
            );
            for (parent, polygon) in members {
                let Some(vertices) = self.evidence.source_polygon_vertices(polygon) else {
                    continue;
                };
                if vertices
                    .iter()
                    .zip(vertices.iter().cycle().skip(1))
                    .all(|(&a, &b)| {
                        carrier
                            .measure_segment(self.evidence.positions[a], self.evidence.positions[b])
                            .is_some_and(|r| r.upper_bound <= tolerance)
                    })
                {
                    split.insert(parent);
                }
            }
        }
        split.into_iter().collect()
    }

    fn bridge_has_independent_profile(&self, [a, b]: [usize; 2]) -> bool {
        // A missing coplanar interval must connect two independently
        // observed neighboring intervals on one non-linear profile.
        for &left in &self.adjacency[a] {
            if left == b {
                continue;
            }
            let Some(profile) = self.fit_station_profile(&[left, a, b]) else {
                continue;
            };
            if self.adjacency[b].iter().any(|&right| {
                if right == a || right == left {
                    return false;
                }
                let station = &self.stations[right];
                profile.error((station.axial, station.radius)) <= self.options.vertex_tolerance
            }) {
                return true;
            }
        }
        false
    }

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
            let stations = self.connected_profile(seed, profile);
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
        Some(Profile::Circular {
            z: z + first.axial,
            r: r + first.radius,
            radius,
        })
    }

    fn connected_profile(&self, seed: &[usize], profile: Profile) -> Vec<usize> {
        let mut visited = seed.iter().copied().collect::<BTreeSet<_>>();
        let mut queue = seed.iter().copied().collect::<VecDeque<_>>();
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
        self.source_regions(stations, profile, 4)
            .into_iter()
            .filter_map(|(core, supports)| {
                self.validate(supports, boundary)
                    .or_else(|| self.validate(core, boundary))
            })
            .collect()
    }

    fn source_regions(
        &self,
        stations: &[usize],
        profile: Profile,
        min_stations: usize,
    ) -> Vec<(Vec<usize>, Vec<usize>)> {
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
            if repeated < min_stations || completed.contains(&seed) {
                continue;
            }
            let supports = self.grow_component(seed, sign, profile, &mut eligible);
            completed.extend(supports.iter().copied());
            candidates.push((core, supports));
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
        let mut independently_sampled = BTreeSet::new();
        for ([a, b], supports) in pairs {
            if !self.adjacency[a].contains(&b) {
                continue;
            }
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
            if distinct_normals(
                self.evidence,
                &supports,
                self.options.normal_tolerance_radians,
            ) {
                independently_sampled.extend([a, b]);
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
                count += usize::from(independently_sampled.contains(&id));
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
        let mut orientation: Option<OrientationEvidence> = None;
        for &polygon in &support.source_face_indices {
            let vertices = self.evidence.source_polygon_vertices(polygon)?;
            let measured = polygon_orientation(
                self.evidence,
                vertices,
                support.normal,
                self.origin,
                self.axis,
                profile,
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
        let profile = Profile::Circular {
            z: (fit.axial_center - self.origin).dot(self.axis),
            r: fit.radial_center,
            radius: fit.profile_radius,
        };
        let repeated_stations = self.repeated_core(&supports, profile)?;
        if repeated_stations < 4 {
            return None;
        }
        let region = boundary.boundary(&supports).ok()?;
        if !region.pinched_vertices.is_empty() || fit.normal_samples != region.source_polygon_count
        {
            return None;
        }
        let chord_tolerance = fit
            .geometric_scale()
            .mul_add(self.options.max_circular_relative_chord_error, tolerance);
        let source_chords = SurfaceCarrier::from_circular_meridian(&fit).measure_source_chords(
            self.evidence,
            &supports,
            chord_tolerance,
        )?;
        Some(StationMeridianCandidate {
            supports,
            source_chords,
            chord_tolerance,
            repeated_stations,
            fit,
            source_polygons: region.source_polygon_count,
            boundary_loops: region.loops.len(),
        })
    }

    fn repeated_core(&self, supports: &[usize], profile: Profile) -> Option<usize> {
        let mut orientation = None;
        let mut cores = BTreeMap::new();
        for &id in supports {
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
        Some(repeated_stations)
    }

    /// Discover affine meridian profiles using the same independent placed
    /// axis, repeated station transitions and original support connectivity.
    /// Two stations interpolate an affine profile; a third independently
    /// tests it. Both remain explicit alternatives for joint selection.
    #[must_use = "linear candidates are competing source region hypotheses"]
    pub fn discover_linear(
        &self,
        boundary: &RegionBoundaryIndex<'_>,
    ) -> Vec<StationLinearCandidate> {
        let mut profiles = BTreeSet::new();
        let mut candidates = BTreeMap::new();
        for transition in self
            .transitions
            .iter()
            .filter(|transition| transition.repeated)
        {
            let [a, b] = transition.stations;
            let first = &self.stations[a];
            let second = &self.stations[b];
            let span = second.axial - first.axial;
            if span.abs() <= self.options.station_resolution {
                continue;
            }
            let slope = (second.radius - first.radius) / span;
            let profile = Profile::Linear {
                radius: slope.mul_add(-first.axial, first.radius),
                slope,
            };
            let stations = self.connected_profile(&[a, b], profile);
            if !profiles.insert(stations.clone()) {
                continue;
            }
            for (core, supports) in self.source_regions(&stations, profile, 2) {
                if let Some(candidate) = self
                    .validate_linear(supports, profile, boundary)
                    .or_else(|| self.validate_linear(core, profile, boundary))
                {
                    candidates
                        .entry(candidate.supports.clone())
                        .or_insert(candidate);
                }
            }
        }
        candidates.into_values().collect()
    }

    fn validate_linear(
        &self,
        supports: Vec<usize>,
        seed: Profile,
        boundary: &RegionBoundaryIndex<'_>,
    ) -> Option<StationLinearCandidate> {
        // Refit all distinct source meridian observations at the fixed axis,
        // giving each resolved station equal weight.
        let ids = supports
            .iter()
            .flat_map(|&id| self.supports[id].stations.iter().copied())
            .collect::<BTreeSet<_>>();
        let first = &self.stations[*ids.first()?];
        let count = crate::usize_to_f64(ids.len());
        let (sum_z, sum_r) = ids.iter().fold((0.0, 0.0), |(z, r), &id| {
            (
                z + self.stations[id].axial - first.axial,
                r + self.stations[id].radius - first.radius,
            )
        });
        let mean_z = sum_z / count;
        let mean_r = sum_r / count;
        let (zz, zr) = ids.iter().fold((0.0, 0.0), |(zz, zr), &id| {
            let z = self.stations[id].axial - first.axial - mean_z;
            let r = self.stations[id].radius - first.radius - mean_r;
            (z.mul_add(z, zz), z.mul_add(r, zr))
        });
        if zz <= self.options.station_resolution.powi(2) {
            return None;
        }
        let slope = zr / zz;
        let radius = slope.mul_add(-(mean_z + first.axial), mean_r + first.radius);
        let profile = Profile::Linear { radius, slope };
        let mut min_z = f64::INFINITY;
        let mut max_z = f64::NEG_INFINITY;
        let mut min_r = f64::INFINITY;
        let mut max_error = 0.0_f64;
        for &id in &supports {
            for &vertex in &self.supports[id].vertices {
                let (z, r) = self.samples[vertex];
                min_z = min_z.min(z);
                max_z = max_z.max(z);
                min_r = min_r.min(r);
                max_error = max_error.max(profile.error((z, r)));
                if seed.error((z, r)) > self.options.vertex_tolerance {
                    return None;
                }
            }
        }
        let scale = min_r.min(max_z - min_z);
        let tolerance = self
            .options
            .vertex_tolerance
            .min(self.options.relative_tolerance * scale);
        if scale <= 0.0 || max_error > tolerance {
            return None;
        }
        let repeated_stations = self.repeated_core(&supports, profile)?;
        if repeated_stations < 2 {
            return None;
        }
        let carrier = SurfaceCarrier::LinearRevolution {
            origin: self.origin,
            axis: self.axis,
            radius_at_origin: radius,
            slope,
        };
        let region = boundary.boundary(&supports).ok()?;
        if !region.pinched_vertices.is_empty() {
            return None;
        }
        let chord_tolerance =
            min_r.mul_add(self.options.max_linear_relative_chord_error, tolerance);
        let source_chords =
            carrier.measure_source_chords(self.evidence, &supports, chord_tolerance)?;
        Some(StationLinearCandidate {
            supports,
            carrier,
            repeated_stations,
            source_polygons: region.source_polygon_count,
            vertex_tolerance: tolerance,
            max_vertex_residual: max_error,
            chord_tolerance,
            source_chords,
        })
    }
}

#[derive(Clone, Copy)]
enum Profile {
    Circular { z: f64, r: f64, radius: f64 },
    Linear { radius: f64, slope: f64 },
}

impl Profile {
    fn error(self, (z, r): (f64, f64)) -> f64 {
        match self {
            Self::Circular {
                z: center,
                r: radial,
                radius,
            } => ((z - center).hypot(r - radial) - radius).abs(),
            Self::Linear { radius, slope } => (r - slope.mul_add(z, radius)).abs(),
        }
    }

    fn normal(self, delta: Vector3, axis: Vector3) -> Option<Vector3> {
        match self {
            Self::Circular { z, r, .. } => circular_profile_normal(delta - axis * z, axis, r),
            Self::Linear { slope, .. } => {
                let radial = delta - axis * delta.dot(axis);
                (radial.magnitude2() > 0.0).then(|| (radial.normalize() - axis * slope).normalize())
            }
        }
    }
}

fn index_stations(
    evidence: &MeshEvidence,
    origin: Point3,
    axis: Vector3,
    clusters: MeridianClusters,
) -> (Vec<MeridianStation>, Vec<usize>) {
    let mut by_vertex = vec![0; evidence.positions.len()];
    let stations = clusters
        .clusters
        .into_iter()
        .enumerate()
        .map(|(id, cluster)| {
            let members = &clusters.order[cluster.members];
            let first = evidence.positions[members[0]] - origin;
            let transverse = first - axis * first.dot(axis);
            let length = transverse.magnitude();
            let mut span = 0.0_f64;
            for &vertex in members {
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
                vertices: members.to_vec(),
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
            let rings = ids
                .iter()
                .all(|&id| station_is_sampled(&stations[id], options.station_resolution));
            let repeated =
                rings && distinct_normals(evidence, &supports, options.normal_tolerance_radians);
            MeridianTransition {
                stations: ids,
                supports,
                sampled: rings,
                repeated,
            }
        })
        .collect()
}

fn station_is_sampled(station: &MeridianStation, resolution: f64) -> bool {
    station.vertices.len() >= 3 && station.transverse_span > resolution
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
    profile: Profile,
    angle_tolerance: f64,
) -> Option<OrientationEvidence> {
    let mut centroid = Vector3::new(0.0, 0.0, 0.0);
    for &id in vertices {
        centroid += evidence.positions[id] - center;
    }
    centroid /= crate::usize_to_f64(vertices.len());
    let expected = profile.normal(centroid, axis)?;
    let dot = expected.dot(normal);
    if dot.abs() <= angle_tolerance.sin() {
        return None;
    }
    let deviation = expected.cross(normal).magnitude().atan2(dot.abs());
    let mut span = 0.0_f64;
    for &id in vertices {
        if let Some(vertex_normal) = profile.normal(evidence.positions[id] - center, axis) {
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
