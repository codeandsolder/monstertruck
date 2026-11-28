//! Locally grow affine radial profiles from independently observed source
//! circles, including irregular end trims and zero-radius cone tips.

use super::super::super::{carrier_identity::append_distinct_proposals, sample_fit::normal_basis};
use super::super::{SourceCircleObservation, StationLinearCandidate};
use super::{
    MeshEvidence, RegionBoundaryIndex, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionProposal,
};
use crate::{InnerSpace, Point3, Vector3};
use std::collections::{BTreeSet, VecDeque};

#[derive(Clone, Copy)]
struct Profile {
    origin: Point3,
    axis: Vector3,
    radius: f64,
    slope: f64,
}

impl Profile {
    const fn from_carrier(carrier: SurfaceCarrier) -> Option<Self> {
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } = carrier
        else {
            return None;
        };
        Some(Self {
            origin,
            axis,
            radius: radius_at_origin,
            slope,
        })
    }

    fn sample(self, point: Point3) -> (f64, Vector3) {
        let delta = point - self.origin;
        let z = delta.dot(self.axis);
        (z, delta - self.axis * z)
    }

    const fn expected(self, z: f64) -> f64 { self.slope.mul_add(z, self.radius) }

    const fn carrier(self) -> SurfaceCarrier {
        SurfaceCarrier::LinearRevolution {
            origin: self.origin,
            axis: self.axis,
            radius_at_origin: self.radius,
            slope: self.slope,
        }
    }
}

// A linear radial profile admits a bound on the entire filled facet.
// For p = sum(w_i p_i) in its convex hull, r(p) = sum(w_i r_i).
// Triangle inequality bounds rho(p) above. For one unit transverse
// direction u, rho(p) >= u.dot(sum(w_i v_i)). Thus vertex-wise bounds
// against the affine budget e + f*r_i certify every convex combination,
// including apex edges and nonconvex polygons. The unit direction is
// chosen from the observed azimuth span; its inequalities are then tested
// directly, rather than trusting the angle calculation as a certificate.
pub(super) fn facet_bound(
    evidence: &MeshEvidence,
    support: usize,
    carrier: SurfaceCarrier,
    tolerance: f64,
    fraction: f64,
) -> Option<f64> {
    let profile = Profile::from_carrier(carrier)?;
    let (u, v) = normal_basis(profile.axis);
    let mut maximum = 0.0_f64;
    for &polygon in &evidence.planar_supports.get(support)?.source_face_indices {
        let ids = evidence.source_polygon_vertices(polygon)?;
        let mut angles = Vec::with_capacity(ids.len());
        let mut samples = Vec::with_capacity(ids.len());
        for &id in ids {
            let (z, radial) = profile.sample(*evidence.positions.get(id)?);
            let expected = profile.expected(z);
            let rho = radial.magnitude();
            if !expected.is_finite() || expected < -tolerance || (rho - expected).abs() > tolerance
            {
                return None;
            }
            if rho > tolerance {
                angles.push(radial.dot(v).atan2(radial.dot(u)));
            }
            samples.push((radial, expected, rho));
        }
        angles.sort_unstable_by(f64::total_cmp);
        let first = *angles.first()?;
        let last = *angles.last()?;
        let mut gap = first + std::f64::consts::TAU - last;
        let mut start = first;
        for pair in angles.windows(2) {
            if pair[1] - pair[0] > gap {
                gap = pair[1] - pair[0];
                start = pair[1];
            }
        }
        let midpoint = (std::f64::consts::TAU - gap).mul_add(0.5, start);
        let (sin, cos) = midpoint.sin_cos();
        let direction = (u * cos + v * sin).normalize();
        for (radial, expected, rho) in samples {
            let budget = fraction.mul_add(expected, tolerance);
            let bound = (expected - radial.dot(direction))
                .max(rho - expected)
                .max(0.0);
            if !budget.is_finite() || budget <= 0.0 || bound > budget {
                return None;
            }
            maximum = maximum.max(bound / budget);
        }
    }
    Some(maximum)
}

#[derive(Clone, Copy)]
pub(in crate::mesh_recovery) struct RimLinearOptions {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub angular_tolerance: f64,
    pub chord_fraction: f64,
}

struct LocalGrowth<'a> {
    evidence: &'a MeshEvidence,
    support_vertices: Vec<Vec<usize>>,
    // Source edges incident to a vertex, with their original support uses.
    ring_edges: Vec<Vec<(usize, usize)>>,
    options: RimLinearOptions,
}

impl<'a> LocalGrowth<'a> {
    fn new(index: &RegionBoundaryIndex<'a>, options: RimLinearOptions) -> Option<Self> {
        let evidence = index.evidence;
        let mut support_vertices = Vec::with_capacity(evidence.planar_supports.len());
        for support in &evidence.planar_supports {
            let mut vertices = Vec::new();
            for &polygon in &support.source_face_indices {
                vertices.extend_from_slice(evidence.source_polygon_vertices(polygon)?);
            }
            vertices.sort_unstable();
            vertices.dedup();
            support_vertices.push(vertices);
        }
        let mut ring_edges = vec![Vec::new(); evidence.positions.len()];
        for (&[a, b], uses) in &index.incidence {
            if uses.len() > 2 || (uses.len() == 2 && uses[0].1 == uses[1].1) {
                continue;
            }
            for &(polygon, _) in uses {
                let support = evidence.source_planar_support[polygon];
                ring_edges[a].push((b, support));
                ring_edges[b].push((a, support));
            }
        }
        Some(Self {
            evidence,
            support_vertices,
            ring_edges,
            options,
        })
    }

    fn seeds(&self, circle: &SourceCircleObservation) -> BTreeSet<usize> {
        let mut seeds = BTreeSet::new();
        for &vertex in &circle.vertices {
            for &(other, support) in &self.ring_edges[vertex] {
                if circle.vertices.binary_search(&other).is_ok() {
                    seeds.insert(support);
                }
            }
        }
        seeds
    }

    fn orientation(&self, profile: Profile, id: usize, tolerance: f64) -> Option<bool> {
        let support = &self.evidence.planar_supports[id];
        let ids = &self.support_vertices[id];
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for &vertex in ids {
            let (z, radial) = profile.sample(self.evidence.positions[vertex]);
            low = low.min(z);
            high = high.max(z);
            if (radial.magnitude() - profile.expected(z)).abs() > tolerance {
                return None;
            }
        }
        if high - low <= tolerance {
            return None; // a rim disk is not a lateral surface
        }
        let carrier = profile.carrier();
        let mut orientation = None;
        for &polygon in &support.source_face_indices {
            let vertices = self.evidence.source_polygon_vertices(polygon)?;
            let centroid = vertices
                .iter()
                .fold(Vector3::new(0.0, 0.0, 0.0), |sum, &id| {
                    sum + (self.evidence.positions[id] - profile.origin)
                })
                / crate::usize_to_f64(vertices.len());
            let expected = carrier.signed_constraint(profile.origin + centroid)?.0;
            let dot = expected.dot(support.normal);
            let deviation = expected.cross(support.normal).magnitude().atan2(dot.abs());
            let span = vertices
                .iter()
                .filter_map(|&id| {
                    carrier
                        .signed_constraint(self.evidence.positions[id])
                        .map(|(normal, _)| {
                            expected
                                .cross(normal)
                                .magnitude()
                                .atan2(expected.dot(normal))
                        })
                })
                .fold(0.0_f64, f64::max);
            if dot.abs() <= self.options.angular_tolerance.sin()
                || deviation
                    > span + self.options.angular_tolerance + support.max_normal_angle_radians
                || orientation.is_some_and(|sign| sign != (dot > 0.0))
            {
                return None;
            }
            orientation = Some(dot > 0.0);
        }
        facet_bound(
            self.evidence,
            id,
            carrier,
            tolerance,
            self.options.chord_fraction,
        )?;
        orientation
    }

    fn grow(&self, profile: Profile, seed: usize, tolerance: f64) -> Vec<usize> {
        let Some(sign) = self.orientation(profile, seed, tolerance) else {
            return Vec::new();
        };
        let mut tested = BTreeSet::from([seed]);
        let mut pending = VecDeque::from([seed]);
        let mut result = vec![seed];
        while let Some(id) = pending.pop_front() {
            for &interface in &self.evidence.support_interface_indices[id] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let other = if a == id { b } else { a };
                if tested.insert(other) && self.orientation(profile, other, tolerance) == Some(sign)
                {
                    result.push(other);
                    pending.push_back(other);
                }
            }
        }
        result.sort_unstable();
        result
    }

    // A source sliver whose entire filled polygon lies within numerical
    // resolution of a single generator has no resolved lateral normal.
    // The placed profile is fixed independently; its affine generator
    // bounds every convex combination of the original vertices.
    fn unresolved_ruling(&self, profile: Profile, support: usize, tolerance: f64) -> bool {
        let ids = &self.support_vertices[support];
        let tolerance = tolerance.min(crate::KERNEL_EPSILON);
        let Some((_, radial)) = ids
            .iter()
            .map(|&id| profile.sample(self.evidence.positions[id]))
            .max_by(|a, b| a.1.magnitude2().total_cmp(&b.1.magnitude2()))
        else {
            return false;
        };
        if radial.magnitude() <= tolerance {
            return false;
        }
        let direction = radial.normalize();
        ids.iter().all(|&id| {
            let (z, radial) = profile.sample(self.evidence.positions[id]);
            let radius = profile.expected(z);
            radius >= -tolerance && (radial - direction * radius).magnitude() <= tolerance
        })
    }

    // Complete only bounded original fans of unresolved generator slivers.
    // They are not independent orientation witnesses and may not become
    // exterior trims. All original supports and polygon IDs remain owned.
    fn complete_ruling_gaps(
        &self,
        index: &RegionBoundaryIndex<'a>,
        profile: Profile,
        mut supports: Vec<usize>,
        tolerance: f64,
    ) -> Option<Vec<usize>> {
        let mut added = BTreeSet::new();
        loop {
            let boundary = index.boundary(&supports).ok()?;
            if boundary.pinched_vertices.is_empty() {
                return (!added.is_empty()
                    && boundary.loops.iter().flatten().all(|edge| {
                        !added.contains(&self.evidence.source_planar_support[edge.interior_polygon])
                    }))
                .then_some(supports);
            }
            let previous = added.len();
            for gap in index.pinch_gaps(&boundary) {
                if gap
                    .added_support_indices
                    .iter()
                    .all(|&id| self.unresolved_ruling(profile, id, tolerance))
                {
                    added.extend(gap.added_support_indices);
                }
            }
            if added.len() == previous {
                return None;
            }
            supports.extend(added.iter().copied());
            supports.sort_unstable();
            supports.dedup();
        }
    }

    // Four connected original samples test a circle beyond interpolation.
    // Follow only edges actually owned by this proposal. An observed
    // interior section is just as valid as an exterior trim.
    fn has_circle_arc(&self, supports: &[usize], samples: &[usize]) -> bool {
        let mut visited = BTreeSet::new();
        for &seed in samples {
            if !visited.insert(seed) {
                continue;
            }
            let mut count = 0;
            let mut pending = VecDeque::from([seed]);
            while let Some(id) = pending.pop_front() {
                count += 1;
                if count >= 4 {
                    return true;
                }
                for &(other, owner) in &self.ring_edges[id] {
                    if supports.binary_search(&owner).is_ok()
                        && samples.binary_search(&other).is_ok()
                        && visited.insert(other)
                    {
                        pending.push_back(other);
                    }
                }
            }
        }
        false
    }

    // A ruling edge is an observation independent of the circle used
    // to place the axis. It must connect two distinct axial positions
    // at the same azimuth (or the independently fitted cone apex).
    fn has_ruling(&self, profile: Profile, supports: &[usize], tolerance: f64) -> bool {
        supports.iter().any(|&support| {
            self.support_vertices[support].iter().any(|&id| {
                self.ring_edges[id].iter().any(|&(other, owner)| {
                    if owner != support || id >= other {
                        return false;
                    }
                    let (za, va) = profile.sample(self.evidence.positions[id]);
                    let (zb, vb) = profile.sample(self.evidence.positions[other]);
                    if (za - zb).abs() <= tolerance {
                        return false;
                    }
                    let ra = va.magnitude();
                    let rb = vb.magnitude();
                    if ra.min(rb) <= tolerance {
                        return profile.expected(if ra < rb { za } else { zb }).abs() <= tolerance;
                    }
                    // Test azimuth in length units at the larger radius.
                    (va / ra - vb / rb).magnitude() * ra.max(rb) <= tolerance
                })
            })
        })
    }

    fn seed_profile(
        &self,
        circle: &SourceCircleObservation,
        seed: usize,
    ) -> Option<(Profile, usize, f64)> {
        let mut profile = Profile {
            origin: circle.center,
            axis: circle.normal,
            radius: circle.radius,
            slope: 0.0,
        };
        let tolerance = self
            .options
            .absolute_tolerance
            .min(self.options.relative_tolerance * circle.radius);
        let &slope_sample = self.support_vertices[seed].iter().max_by(|&&a, &&b| {
            profile
                .sample(self.evidence.positions[a])
                .0
                .abs()
                .total_cmp(&profile.sample(self.evidence.positions[b]).0.abs())
        })?;
        let (z, radial) = profile.sample(self.evidence.positions[slope_sample]);
        if z.abs() <= tolerance {
            return None;
        }
        profile.slope = (radial.magnitude() - circle.radius) / z;
        Some((profile, slope_sample, tolerance))
    }

    fn proposal(
        &self,
        index: &RegionBoundaryIndex<'a>,
        circle: &SourceCircleObservation,
        seed: usize,
        (profile, slope_sample, tolerance): (Profile, usize, f64),
    ) -> Option<SurfaceRegionProposal<'a>> {
        let supports = self.grow(profile, seed, tolerance);
        if supports.len() < 2 {
            return None;
        }
        let first_normal = self.evidence.planar_supports[supports[0]].normal;
        if !supports.iter().any(|&id| {
            first_normal
                .cross(self.evidence.planar_supports[id].normal)
                .magnitude()
                > self.options.angular_tolerance.sin()
        }) {
            return None;
        }
        let ids = supports
            .iter()
            .flat_map(|&id| self.support_vertices[id].iter().copied())
            .collect::<BTreeSet<_>>();
        // The sampled circle must occur as a connected source arc inside
        // this carrier's ownership, including legitimate interior sections.
        if !self.has_circle_arc(&supports, &circle.vertices) {
            return None;
        }
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        let mut axial = Vec::with_capacity(ids.len());
        for &id in &ids {
            let (z, _) = profile.sample(self.evidence.positions[id]);
            low = low.min(z);
            high = high.max(z);
            axial.push(z);
        }
        let tolerance = tolerance.min(self.options.relative_tolerance * (high - low));
        // The rim fixes the axis and radius; one off-rim point fixes the
        // slope. Do not equate exact interpolation with independent
        // identification. Arbitrary triangulations need not expose a
        // ruling edge: an independently sampled second section or a
        // full-rank set of off-rim observations also suffices.
        let off_rim = ids
            .iter()
            .copied()
            .filter(|id| circle.vertices.binary_search(id).is_err())
            .collect::<Vec<_>>();
        if !self.has_ruling(profile, &supports, tolerance)
            && !has_independent_section(self.evidence, profile, &off_rim, tolerance)
            && !independent_profile_rank(
                self.evidence,
                profile,
                &off_rim
                    .into_iter()
                    .filter(|&id| id != slope_sample)
                    .collect::<Vec<_>>(),
                tolerance,
            )
        {
            return None;
        }
        let proposal = index
            .measure_rim_linear_region(
                &supports,
                profile.carrier(),
                tolerance,
                self.options.chord_fraction,
            )
            .or_else(|| {
                let completed = self.complete_ruling_gaps(index, profile, supports, tolerance)?;
                index.measure_rim_linear_region(
                    &completed,
                    profile.carrier(),
                    tolerance,
                    self.options.chord_fraction,
                )
            })?;
        // Two axial samples interpolate a slope. A third resolved sample
        // tests it independently, even if the off-rim samples are irregular.
        let independently_tested = axial
            .iter()
            .any(|&z| z - low > tolerance && high - z > tolerance);
        Some(if independently_tested {
            proposal.with_independent_profile()
        } else {
            proposal
        })
    }
}

// Two independently resolved cross sections establish the same placed axis
// and affine profile even on a triangulation without explicit ruling edges.
fn has_independent_section(
    evidence: &MeshEvidence,
    profile: Profile,
    ids: &[usize],
    tolerance: f64,
) -> bool {
    let mut samples = ids
        .iter()
        .map(|&id| (profile.sample(evidence.positions[id]).0, id))
        .collect::<Vec<_>>();
    samples.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let mut start = 0;
    while start < samples.len() {
        let z = samples[start].0;
        let end = start + samples[start..].partition_point(|&(value, _)| value - z <= tolerance);
        if z.abs() > tolerance && end - start >= 4 {
            let station = samples[start..end]
                .iter()
                .map(|&(_, id)| id)
                .collect::<Vec<_>>();
            if let Some(fit) =
                super::super::super::sample_fit::fit_circle_carrier(&evidence.positions, &station)
            {
                let delta = fit.center - profile.origin;
                if fit.source_span * fit.plane_span_strength > tolerance
                    && fit.max_plane_residual.max(fit.max_radial_residual) <= tolerance
                    && delta.cross(profile.axis).magnitude() <= tolerance
                    && fit.normal.cross(profile.axis).magnitude() * fit.radius <= tolerance
                    && (fit.radius - profile.expected(delta.dot(profile.axis))).abs() <= tolerance
                {
                    return true;
                }
            }
        }
        start = end;
    }
    false
}

// Jacobian of rho-r-s*z for a general affine radial surface. The six
// local parameters are transverse translations, axis tilts, radius and
// slope. Exclude the seed observations and require an overdetermined,
// full-rank system in its own right; merely counting distinct axial
// coordinates would allow nearly degenerate planar interpolations.
// Normalize length parameters by the sampled extent and use column-pivoted
// modified Gram-Schmidt to measure rank without squaring conditioning.
fn independent_profile_rank(
    evidence: &MeshEvidence,
    profile: Profile,
    ids: &[usize],
    tolerance: f64,
) -> bool {
    const PARAMETERS: usize = 6;
    if ids.len() <= PARAMETERS {
        return false;
    }
    let (u, v) = normal_basis(profile.axis);
    let samples = ids
        .iter()
        .map(|&id| profile.sample(evidence.positions[id]))
        .collect::<Vec<_>>();
    let center = samples.iter().map(|&(z, _)| z).sum::<f64>() / crate::usize_to_f64(samples.len());
    let scale = samples
        .iter()
        .map(|&(z, radial)| (z - center).abs().max(radial.magnitude()))
        .fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale <= tolerance {
        return false;
    }
    let mut columns: [Vec<f64>; PARAMETERS] =
        std::array::from_fn(|_| Vec::with_capacity(samples.len()));
    for (z, radial) in samples {
        let radius = radial.magnitude();
        if radius <= tolerance {
            continue; // the apex has no unique radial derivative
        }
        let direction = radial / radius;
        let axial = (z - center) / scale;
        let tilt = profile.slope.mul_add(radius / scale, axial);
        let a = direction.dot(u);
        let b = direction.dot(v);
        for (column, value) in columns
            .iter_mut()
            .zip([a, b, a * tilt, b * tilt, 1.0, axial])
        {
            column.push(value);
        }
    }
    if columns[0].len() <= PARAMETERS {
        return false;
    }
    let norm = |col: &[f64]| col.iter().map(|v| v * v).sum::<f64>().sqrt();
    let resolution =
        (tolerance / scale).max(64.0 * f64::EPSILON) * crate::usize_to_f64(columns[0].len()).sqrt();
    for i in 0..PARAMETERS {
        let Some(pivot) =
            (i..PARAMETERS).max_by(|&a, &b| norm(&columns[a]).total_cmp(&norm(&columns[b])))
        else {
            return false;
        };
        columns.swap(i, pivot);
        let length = norm(&columns[i]);
        if !length.is_finite() || length <= resolution {
            return false;
        }
        let (head, tail) = columns.split_at_mut(i + 1);
        let basis = &mut head[i];
        for value in basis.iter_mut() {
            *value /= length;
        }
        for column in tail {
            // Reorthogonalize the remaining column against the pivot.
            for _ in 0..2 {
                let projection = column
                    .iter()
                    .zip(basis.iter())
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                for (value, &direction) in column.iter_mut().zip(basis.iter()) {
                    *value = (-projection).mul_add(direction, *value);
                }
            }
        }
    }
    true
}

impl<'a> RegionBoundaryIndex<'a> {
    // A placed station profile has two independently repeated sections.
    // Extend it through irregularly triangulated lateral strips without
    // refitting its axis or radius. As with fixed circular growth, use
    // the declared relative source-fit envelope on the independently
    // observed core scale, not the station-identification resolution.
    // The original continuous chord budget remains a separate bound.
    pub(in crate::mesh_recovery) fn expand_station_linear_regions(
        &self,
        candidates: &[&StationLinearCandidate],
        options: RimLinearOptions,
    ) -> Vec<SurfaceRegionProposal<'a>> {
        if candidates.is_empty() {
            return Vec::new();
        }
        let Some(growth) = LocalGrowth::new(self, options) else {
            return Vec::new();
        };
        candidates
            .iter()
            .filter_map(|candidate| {
                let profile = Profile::from_carrier(candidate.carrier)?;
                let seed = *candidate.supports.first()?;
                let (low, high, radius) = candidate
                    .supports
                    .iter()
                    .flat_map(|&id| &growth.support_vertices[id])
                    .map(|&id| profile.sample(self.evidence.positions[id]))
                    .fold(
                        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY),
                        |(low, high, radius), (z, radial)| {
                            (low.min(z), high.max(z), radius.min(radial.magnitude()))
                        },
                    );
                let tolerance = options.relative_tolerance * radius.min(high - low);
                if !tolerance.is_finite() || tolerance <= 0.0 {
                    return None;
                }
                let supports = growth.grow(profile, seed, tolerance);
                if supports.len() <= candidate.supports.len()
                    || !candidate
                        .supports
                        .iter()
                        .all(|id| supports.binary_search(id).is_ok())
                {
                    return None;
                }
                self.measure_surface_region(
                    &supports,
                    candidate.carrier,
                    SurfaceRegionBudget {
                        vertex: tolerance,
                        chord: candidate.chord_tolerance,
                    },
                )
            })
            .collect()
    }

    pub(in crate::mesh_recovery) fn discover_rim_linear_regions(
        &self,
        circles: &[SourceCircleObservation],
        options: RimLinearOptions,
    ) -> Vec<SurfaceRegionProposal<'a>> {
        let Some(growth) = LocalGrowth::new(self, options) else {
            return Vec::new();
        };
        let Some(&reference) = self.evidence.positions.first() else {
            return Vec::new();
        };
        let extent = self
            .evidence
            .positions
            .iter()
            .map(|&point| (point - reference).magnitude())
            .fold(0.0_f64, f64::max);
        let identity_tolerance = options.absolute_tolerance.min(crate::KERNEL_EPSILON * 0.1);
        let mut candidates = Vec::new();
        for circle in circles {
            let mut completed = Vec::<(f64, Vec<usize>)>::new();
            // The same circle fixes all parameters except slope. Bound
            // profile differences over the entire possible growth domain,
            // not just the seed support, before skipping another seed.
            let axial_extent = (reference - circle.center).dot(circle.normal).abs() + extent;
            for seed in growth.seeds(circle) {
                let Some(parameters) = growth.seed_profile(circle, seed) else {
                    continue;
                };
                let slope = parameters.0.slope;
                if completed.iter().any(|(previous, supports)| {
                    supports.binary_search(&seed).is_ok()
                        && (previous - slope).abs() * axial_extent <= identity_tolerance
                }) {
                    continue;
                }
                if let Some(proposal) = growth.proposal(self, circle, seed, parameters) {
                    completed.push((slope, proposal.supports().to_vec()));
                    candidates.push(proposal);
                }
            }
        }
        let mut proposals = Vec::new();
        append_distinct_proposals(
            self.evidence,
            candidates,
            identity_tolerance,
            &mut proposals,
        );
        proposals
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{
        CircleAxisOptions, EvidenceOptions, analyze_mesh_evidence, discover_circle_axes,
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    fn options() -> RimLinearOptions {
        RimLinearOptions {
            absolute_tolerance: 1.0e-8,
            relative_tolerance: 1.0e-4,
            angular_tolerance: 1.0e-5,
            chord_fraction: 0.25,
        }
    }

    fn discover<'a>(index: &RegionBoundaryIndex<'a>) -> Vec<SurfaceRegionProposal<'a>> {
        let circles = discover_circle_axes(
            index.evidence,
            CircleAxisOptions {
                vertex_tolerance: 1.0e-8,
                relative_tolerance: 1.0e-4,
                axis_position_tolerance: 8.0e-8,
                axis_angle_tolerance_radians: 1.0e-7,
            },
        );
        circles.map_or_else(Vec::new, |c| {
            index.discover_rim_linear_regions(&c.circles, options())
        })
    }

    fn cylinder_mesh(n: usize) -> Result<PolygonMesh, Box<dyn std::error::Error>> {
        let mut positions = Vec::new();
        for row in 0..2 {
            for i in 0..n {
                let a = std::f64::consts::TAU * crate::usize_to_f64(i) / crate::usize_to_f64(n);
                let (sin, cos) = a.sin_cos();
                // No second regular axial ring. The trim plane is oblique.
                let z = if row == 0 {
                    0.0
                } else {
                    2.0 + 0.3 * cos + 0.17 * sin
                };
                positions.push(Point3::new(cos, sin, z));
            }
        }
        let quads = (0..n).flat_map(|i| {
            let next = (i + 1) % n;
            [[i, next, n + next], [i, n + next, n + i]]
        });
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(quads),
        )?)
    }

    #[test]
    fn irregular_trim_cylinder_is_recovered_without_repeated_stations()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = cylinder_mesh(12)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let proposals = discover(&index);
        let p = proposals
            .iter()
            .find(|p| p.source_polygons() == 24)
            .ok_or("missing full cylinder")?;
        assert!(p.independently_tested_profile());
        let SurfaceCarrier::LinearRevolution {
            radius_at_origin,
            axis,
            slope,
            ..
        } = p.carrier()
        else {
            return Err("expected linear carrier".into());
        };
        assert!((radius_at_origin - 1.0).abs() < 1.0e-9);
        assert!(axis.cross(Vector3::unit_z()).magnitude() < 1.0e-9);
        assert!(slope.abs() < 1.0e-9);
        Ok(())
    }

    #[test]
    fn station_linear_core_extends_to_irregular_strip_without_refitting()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = cylinder_mesh(12)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let carrier = SurfaceCarrier::LinearRevolution {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radius_at_origin: 1.0,
            slope: 0.0,
        };
        let mut supports = vec![
            evidence.source_planar_support[0],
            evidence.source_planar_support[2],
        ];
        supports.sort_unstable();
        supports.dedup();
        let chord_tolerance = 0.25;
        let candidate = StationLinearCandidate {
            source_chords: carrier
                .measure_source_chords(&evidence, &supports, chord_tolerance)
                .ok_or("missing chord evidence")?,
            source_polygons: index.boundary(&supports)?.source_polygon_count,
            supports,
            carrier,
            repeated_stations: 2,
            vertex_tolerance: 1.0e-8,
            max_vertex_residual: 0.0,
            chord_tolerance,
        };
        let expanded = index.expand_station_linear_regions(&[&candidate], options());
        assert_eq!(expanded.len(), 1);
        assert_eq!(expanded[0].source_polygons(), 24);
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } = expanded[0].carrier()
        else {
            return Err("expected affine carrier".into());
        };
        assert_eq!(origin, Point3::new(0.0, 0.0, 0.0));
        assert_eq!(axis, Vector3::unit_z());
        assert_eq!(radius_at_origin.to_bits(), 1.0_f64.to_bits());
        assert_eq!(slope.to_bits(), 0.0_f64.to_bits());
        Ok(())
    }

    #[test]
    fn station_growth_distinguishes_source_fit_from_station_resolution()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{MeridianGraphOptions, MeridianStationGraph};
        let n = 12;
        let positions = (0..3)
            .flat_map(|row| {
                (0..n).map(move |i| {
                    let a = std::f64::consts::TAU * crate::usize_to_f64(i) / crate::usize_to_f64(n);
                    let (sin, cos) = a.sin_cos();
                    let (r, z) = if row == 2 {
                        (1.00002, 2.0 + 0.3 * cos + 0.17 * sin)
                    } else {
                        (1.0, crate::usize_to_f64(row))
                    };
                    Point3::new(r * cos, r * sin, z)
                })
            })
            .collect();
        let mut faces = (0..2)
            .flat_map(|row| {
                (0..n).flat_map(move |i| {
                    let a = row * n + i;
                    let b = row * n + (i + 1) % n;
                    [vec![a, b, b + n], vec![a, b + n, a + n]]
                })
            })
            .collect::<Vec<_>>();
        faces.push((0..n).rev().collect());
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let graph = MeridianStationGraph::new(
            &evidence,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            MeridianGraphOptions {
                station_resolution: 1.0e-8,
                vertex_tolerance: 1.0e-8,
                relative_tolerance: 1.0e-4,
                max_linear_relative_chord_error: 0.25,
                max_circular_relative_chord_error: 0.25,
                normal_tolerance_radians: 1.0e-5,
            },
        )
        .ok_or("missing station graph")?;
        let candidates = graph.discover_linear(&index);
        let core = candidates
            .iter()
            .find(|c| c.source_polygons == 2 * n)
            .ok_or("missing independently repeated core")?;
        let expanded = index.expand_station_linear_regions(&[core], options());
        assert_eq!(expanded.len(), 1);
        assert_eq!(expanded[0].source_polygons(), 4 * n);
        // The cap has exact on-cylinder vertices, but no lateral extent.
        assert!(
            !expanded[0]
                .supports()
                .contains(&evidence.source_planar_support[4 * n])
        );
        // Tightening the declared source-fit envelope rejects the same
        // noisy strip without touching the independently measured core.
        assert!(
            index
                .expand_station_linear_regions(
                    &[core],
                    RimLinearOptions {
                        relative_tolerance: 1.0e-6,
                        ..options()
                    }
                )
                .is_empty()
        );
        Ok(())
    }

    // A source triangle on an otherwise common generator can have
    // arbitrary normal after rounding. Its original fan is still known.
    fn cone_with_ruling_sliver(offset: f64) -> Result<PolygonMesh, Box<dyn std::error::Error>> {
        let n = 16;
        let mut positions = vec![Point3::new(0.0, 0.0, 0.0)];
        for i in 0..n {
            let angle = std::f64::consts::TAU * crate::usize_to_f64(i) / crate::usize_to_f64(n);
            let (sin, cos) = angle.sin_cos();
            positions.push(Point3::new(cos, sin, 1.0));
        }
        let mid = positions.len();
        positions.push(Point3::new(0.5 + offset, 0.0, 0.5));
        let mut faces = vec![[0, 1, mid], [0, mid, 2], [mid, 1, 2]];
        faces.extend((1..n).map(|i| [0, i + 1, (i + 1) % n + 1]));
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?)
    }

    #[test]
    fn bounded_ruling_gap_preserves_original_facet_ownership()
    -> Result<(), Box<dyn std::error::Error>> {
        for offset in [1.0e-10, 1.0e-3] {
            let mesh = cone_with_ruling_sliver(offset)?;
            let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
            let index = RegionBoundaryIndex::new(&evidence)?;
            let growth = LocalGrowth::new(&index, options()).ok_or("local growth")?;
            let profile = Profile {
                origin: Point3::new(0.0, 0.0, 0.0),
                axis: Vector3::new(0.0, 0.0, 1.0),
                radius: 0.0,
                slope: 1.0,
            };
            let missing = evidence.source_planar_support[0];
            let supports = (0..evidence.planar_supports.len())
                .filter(|&id| id != missing)
                .collect::<Vec<_>>();
            assert!(!index.boundary(&supports)?.pinched_vertices.is_empty());
            let completed = growth.complete_ruling_gaps(&index, profile, supports, 1.0e-8);
            if offset < 1.0e-8 {
                let completed = completed.ok_or("missing completed fan")?;
                let boundary = index.boundary(&completed)?;
                assert!(boundary.pinched_vertices.is_empty());
                assert_eq!(
                    boundary.source_polygon_count,
                    evidence.source_planar_support.len()
                );
                assert!(
                    index
                        .measure_rim_linear_region(&completed, profile.carrier(), 1.0e-8, 0.25)
                        .is_some()
                );
            } else {
                assert!(completed.is_none());
            }
        }
        Ok(())
    }

    #[test]
    fn interior_circle_can_seed_both_sides_of_a_cone() -> Result<(), Box<dyn std::error::Error>> {
        let n = 16;
        let mut positions = Vec::new();
        for row in 0..3 {
            for j in 0..n {
                let angle = std::f64::consts::TAU * crate::usize_to_f64(j) / crate::usize_to_f64(n);
                let (sin, cos) = angle.sin_cos();
                let z = match row {
                    0 => -0.8 + 0.1 * sin,
                    1 => 0.0,
                    _ => 0.8 + 0.1 * cos,
                };
                let radius = 0.3_f64.mul_add(z, 1.0);
                positions.push(Point3::new(radius * cos, radius * sin, z));
            }
        }
        let triangles = (0..2).flat_map(|row| {
            (0..n).flat_map(move |j| {
                let a = row * n + j;
                let b = row * n + (j + 1) % n;
                [[a, b, b + n], [a, b + n, a + n]]
            })
        });
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(triangles),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let circle = SourceCircleObservation {
            vertices: evidence
                .positions
                .iter()
                .enumerate()
                .filter_map(|(id, p)| (p.z == 0.0).then_some(id))
                .collect(),
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::unit_z(),
            radius: 1.0,
            max_plane_residual: 0.0,
            max_radial_residual: 0.0,
            resolved_curvature_span: 1.0,
        };
        let candidates = index.discover_rim_linear_regions(&[circle], options());
        assert!(
            candidates
                .iter()
                .any(|proposal| proposal.source_polygons() == 4 * n)
        );
        Ok(())
    }

    // Oblique staggered triangles do not contain exact generator edges
    // or repeated off-rim sections. Independently informative scattered
    // samples must still recover a general affine radial profile.
    #[test]
    fn scattered_cone_without_generator_edges_is_recovered()
    -> Result<(), Box<dyn std::error::Error>> {
        let n = 16;
        let mut positions = Vec::new();
        for row in 0..4 {
            let height = crate::usize_to_f64(row);
            for j in 0..n {
                let angle = std::f64::consts::TAU * crate::usize_to_f64(j) / crate::usize_to_f64(n)
                    + height * 0.13;
                let (sin, cos) = angle.sin_cos();
                let z = if row == 0 {
                    0.0
                } else {
                    0.7 * height + 0.05 * sin + 0.03 * cos
                };
                let radius = 1.0 + 0.2 * z;
                positions.push(Point3::new(radius * cos, radius * sin, z));
            }
        }
        let triangles = (0..3).flat_map(|row| {
            (0..n).flat_map(move |j| {
                let a = row * n + j;
                let b = row * n + (j + 1) % n;
                [[a, b, b + n], [a, b + n, a + n]]
            })
        });
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(triangles),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let candidates = discover(&index);
        let proposal = candidates
            .iter()
            .find(|p| p.source_polygons() == 6 * n)
            .ok_or("missing scattered cone")?;
        let SurfaceCarrier::LinearRevolution { axis, slope, .. } = proposal.carrier() else {
            return Err("expected linear carrier".into());
        };
        assert!(axis.cross(Vector3::unit_z()).magnitude() < 1e-9);
        assert!((slope.abs() - 0.2).abs() < 1e-9);
        let profile = Profile::from_carrier(proposal.carrier()).ok_or("profile")?;
        let growth = LocalGrowth::new(&index, options()).ok_or("growth")?;
        assert!(!growth.has_ruling(profile, proposal.supports(), 1e-8));
        Ok(())
    }

    // The seed circle and one more point can always interpolate a slope.
    // A two-point off-rim row supplies neither another resolved circle nor
    // a full-rank independent 3D test of the surface.
    #[test]
    fn interpolation_only_does_not_certify_a_rim_surface() -> Result<(), Box<dyn std::error::Error>>
    {
        let samples = [
            (0.0_f64, 0.0_f64),
            (0.4, 0.0),
            (0.8, 0.0),
            (1.2, 0.0),
            (0.3, 1.0),
            (0.9, 1.0),
        ];
        let positions = samples
            .into_iter()
            .map(|(a, z)| {
                let (sin, cos) = a.sin_cos();
                let r = 0.2_f64.mul_add(z, 1.0);
                Point3::new(r * cos, r * sin, z)
            })
            .collect();
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 4], [1, 5, 4], [1, 2, 5], [2, 3, 5]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let circle = SourceCircleObservation {
            vertices: evidence
                .positions
                .iter()
                .enumerate()
                .filter_map(|(id, p)| (p.z == 0.0).then_some(id))
                .collect(),
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vector3::unit_z(),
            radius: 1.0,
            max_plane_residual: 0.0,
            max_radial_residual: 0.0,
            resolved_curvature_span: 0.1,
        };
        assert!(
            index
                .discover_rim_linear_regions(&[circle], options())
                .is_empty()
        );
        Ok(())
    }

    #[test]
    fn coarse_square_prism_is_not_a_cylinder() -> Result<(), Box<dyn std::error::Error>> {
        let mesh = cylinder_mesh(4)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        assert!(discover(&index).iter().all(|p| p.source_polygons() < 8));
        Ok(())
    }

    #[test]
    fn cone_with_apex_passes_local_radius_approximation() -> Result<(), Box<dyn std::error::Error>>
    {
        let n = 12;
        let mut positions = Vec::new();
        for row in 0..2 {
            let radius = if row == 0 { 1.0 } else { 0.5 };
            for i in 0..n {
                let a = std::f64::consts::TAU * crate::usize_to_f64(i) / crate::usize_to_f64(n);
                let (sin, cos) = a.sin_cos();
                positions.push(Point3::new(radius * cos, radius * sin, 1.0 - radius));
            }
        }
        positions.push(Point3::new(0.0, 0.0, 1.0));
        let triangles = (0..n).flat_map(|i| {
            let next = (i + 1) % n;
            [
                [i, next, n + next],
                [i, n + next, n + i],
                [n + i, n + next, 2 * n],
            ]
        });
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(triangles),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let proposals = discover(&index);
        let p = proposals
            .iter()
            .find(|p| p.source_polygons() == 3 * n)
            .ok_or("missing apex-inclusive cone")?;
        assert!(p.independently_tested_profile());
        assert!(p.relative_chord_bound() < 1.0);
        let recovered = crate::mesh_recovery::recover_analytic_mesh(&mesh, Default::default())?;
        let shell = recovered
            .complete_shell()
            .ok_or("missing complete cone cap")?;
        assert_eq!(shell.len(), 2);
        assert_eq!(recovered.native.regions.len(), 1);
        assert!(
            recovered
                .native
                .regions
                .iter()
                .all(|r| !r.requires_chart_cut)
        );

        Ok(())
    }
}
