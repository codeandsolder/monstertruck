//! Independent *spatial* evidence for rotational source-support tracks.
//!
//! A consistent sequence of normals is not enough to establish a cylinder,
//! cone, or even one continuous CAD face. These fits use original source
//! support-plane offsets and every canonical source polygon vertex, without
//! assuming the identity or primitive family of any original STEP surface.

use super::super::sample_fit::normal_basis;
use super::{CurvatureTrack, MeshEvidence};
use crate::{InnerSpace, Point3, Vector3};
use rustc_hash::FxHashSet as HashSet;

/// Approximate cross-sectional source vertex cluster at one axial station.
///
/// Axial coordinates are relative to the reported fitted axis origin.
/// Grouping is based only on source coordinate and inferred axis uncertainty,
/// not original CAD face boundaries or any assumed cylinder/torus.
#[derive(Clone, Debug)]
pub struct AxialRadialStation {
    /// Mean signed axial coordinate of this sample group.
    pub axial_mean: f64,
    /// Mean source-vertex radius around the candidate axis at this station.
    pub radius_mean: f64,
    /// Largest radius minus smallest radius within this station.
    pub radial_range: f64,
    /// Number of unique original source vertices in the station.
    pub vertex_count: usize,
}

/// Independent source-position measurements for a candidate axisymmetric
/// surface. This is descriptive fit evidence, not a CAD surface certificate.
///
/// Plane offsets are modeled with a fixed transverse axis location and one
/// shared scalar offset. Separately, original source vertices are tested
/// against a best-fit radius that changes linearly with axial coordinate.
/// Cylinders have zero slope; cones have a nonzero slope; the same tests can
/// expose inconsistent planar or toroidal source strips.
#[derive(Clone, Debug)]
pub struct SpatialTrackEvidence {
    /// Fixed axis direction from the independently recovered normal track.
    pub axis: Vector3,
    /// A fitted point on the candidate axis, in original world coordinates.
    pub axis_origin: Point3,
    /// Smaller/larger eigenvalue ratio of the support normal transverse
    /// covariance matrix. Near zero means axis placement is ill-conditioned.
    pub transverse_normal_rank_ratio: f64,
    /// Shared signed offset of source support planes from the fitted axis.
    /// For a faceted cylinder this is the signed polygon apothem, not its
    /// circumradius; for a cone the offset depends on axial reference.
    pub signed_support_plane_offset: f64,
    /// Largest original source support-plane offset error at seed origins.
    pub max_support_plane_offset_residual: f64,
    /// RMS original source support-plane offset error at seed origins.
    pub rms_support_plane_offset_residual: f64,
    /// Minimum normal projection onto the candidate axis over *all* source
    /// planar supports, independent of the fitted source vertex radius.
    pub min_support_normal_axial_projection: f64,
    /// Maximum normal projection onto the same candidate axis.
    pub max_support_normal_axial_projection: f64,
    /// Difference of the extrema of support-normal axial projections. A
    /// perfect cylinder has zero; a true cone also has zero on one common
    /// conical flank, but neighboring toroidal strips may change latitude.
    pub support_normal_axial_projection_spread: f64,
    /// Number of *unique canonical source vertices* sampled across all source
    /// polygons of the track, including non-boundary/interior vertices.
    pub unique_source_vertex_count: usize,
    /// Longitudinal extent of those original source points along the axis.
    pub axial_coordinate_span: f64,
    /// Axial source vertex clusters, each with mean axial position, radius,
    /// within-station radial spread, and provenance-preserving vertex count.
    /// One/two stations cannot establish nonzero meridional curvature.
    pub axial_stations: Vec<AxialRadialStation>,
    /// Relative/axis-uncertainty-derived axial resolution of the clustering.
    pub axial_station_resolution: f64,
    /// Maximum absolute change of mean-radius-versus-axial-position slope
    /// between adjacent station intervals. None with fewer than 3 stations.
    /// This describes meridional curvature, not a torus classification.
    pub max_meridian_slope_change: Option<f64>,
    /// Mean observed radial distance of source vertices from the fitted axis.
    pub mean_sample_radius: f64,
    /// Best-fit radial distance at the reference axis origin.
    pub fitted_sample_radius_at_origin: f64,
    /// Independent linear radius-versus-axial-position slope. None means the
    /// sampled source points do not span enough axial space to determine it.
    pub fitted_sample_radius_slope: Option<f64>,
    /// Maximum absolute vertex radial residual to that independent profile.
    pub max_sample_radial_residual: f64,
    /// RMS vertex radial residual to the independently fitted profile.
    pub rms_sample_radial_residual: f64,
    /// Maximum radial residual divided by mean source-vertex radius.
    pub relative_max_sample_radial_residual: f64,
    /// Axial radius slope suggested by the fitted source support normals;
    /// undefined if their cone crosses its apex within the sampled span.
    pub normal_implied_radius_slope: Option<f64>,
    /// Absolute disagreement between the independent original-vertex
    /// circumradius slope and source-support-plane apothem slope. This is
    /// nonzero even for perfectly tessellated conical walls because the
    /// apothem differs from the circumradius by a facet-angle-dependent term.
    pub radius_slope_disagreement: Option<f64>,
}

struct PlaneFit {
    origin: Point3,
    signed_offset: f64,
    rank_ratio: f64,
    max_residual: f64,
    rms_residual: f64,
}

fn fit_support_plane_placement(
    evidence: &MeshEvidence,
    support_ids: &[usize],
    axis: Vector3,
) -> Option<PlaneFit> {
    let reference = evidence.planar_supports.get(*support_ids.first()?)?.origin;
    let (basis_u, basis_v) = normal_basis(axis);
    let data = support_ids
        .iter()
        .map(|&id| {
            let support = evidence.planar_supports.get(id)?;
            let relative = support.origin - reference;
            Some((
                support.normal.dot(basis_u),
                support.normal.dot(basis_v),
                support.normal.dot(relative),
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    let count = crate::usize_to_f64(data.len());
    if count < 3.0 {
        return None;
    }
    let (mean_x, mean_y, mean_d) = data
        .iter()
        .fold((0.0, 0.0, 0.0), |(sx, sy, sd), &(x, y, d)| {
            (sx + x, sy + y, sd + d)
        });
    let (mean_x, mean_y, mean_d) = (mean_x / count, mean_y / count, mean_d / count);
    let (mut xx, mut xy, mut yy, mut dx, mut dy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y, d) in &data {
        let centered_x = x - mean_x;
        let centered_y = y - mean_y;
        let centered_offset = d - mean_d;
        xx = centered_x.mul_add(centered_x, xx);
        xy = centered_x.mul_add(centered_y, xy);
        yy = centered_y.mul_add(centered_y, yy);
        dx = centered_x.mul_add(centered_offset, dx);
        dy = centered_y.mul_add(centered_offset, dy);
    }
    let det = xx.mul_add(yy, -(xy * xy));
    let trace = xx + yy;
    let largest = trace.midpoint((xx - yy).hypot(2.0 * xy));
    if !det.is_finite()
        || !largest.is_finite()
        || largest <= 0.0
        || det <= f64::EPSILON * largest * largest
    {
        return None;
    }
    let rank_ratio = det / (largest * largest);
    let cx = (dx.mul_add(yy, -(dy * xy))) / det;
    let cy = (dy.mul_add(xx, -(dx * xy))) / det;
    let offset = mean_d - mean_x.mul_add(cx, mean_y * cy);
    let (max_residual, sum_squares) =
        data.iter()
            .fold((0.0_f64, 0.0_f64), |(maximum, sum), &(x, y, d)| {
                let residual = d - x.mul_add(cx, y.mul_add(cy, offset));
                (maximum.max(residual.abs()), residual.mul_add(residual, sum))
            });
    let origin = reference + basis_u * cx + basis_v * cy;
    ([
        origin.x,
        origin.y,
        origin.z,
        offset,
        rank_ratio,
        max_residual,
        sum_squares,
    ]
    .into_iter()
    .all(f64::is_finite))
    .then_some(PlaneFit {
        origin,
        signed_offset: offset,
        rank_ratio,
        max_residual,
        rms_residual: (sum_squares / count).sqrt(),
    })
}

/// Group all original vertex radii by resolved axial coordinate, retaining
/// source sample multiplicity. A rank-two axial profile cannot distinguish a
/// real cone from one secant strip of a torus; do not infer curvature there.
fn group_axial_stations(
    samples: &[(f64, f64)],
    tolerance: f64,
) -> (Vec<AxialRadialStation>, Option<f64>) {
    let mut sorted = samples.to_vec();
    sorted.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let mut stations = Vec::new();
    let mut start = 0;
    while start < sorted.len() {
        let axial_anchor = sorted[start].0;
        let mut end = start;
        let (mut axial_sum, mut radius_sum) = (0.0_f64, 0.0_f64);
        let (mut min_radius, mut max_radius) = (f64::INFINITY, f64::NEG_INFINITY);
        while end < sorted.len() && sorted[end].0 - axial_anchor <= tolerance {
            let (axial, radius) = sorted[end];
            axial_sum += axial;
            radius_sum += radius;
            min_radius = min_radius.min(radius);
            max_radius = max_radius.max(radius);
            end += 1;
        }
        let count = crate::usize_to_f64(end - start);
        stations.push(AxialRadialStation {
            axial_mean: axial_sum / count,
            radius_mean: radius_sum / count,
            radial_range: max_radius - min_radius,
            vertex_count: end - start,
        });
        start = end;
    }
    let slopes = stations
        .windows(2)
        .map(|pair| {
            (pair[1].radius_mean - pair[0].radius_mean) / (pair[1].axial_mean - pair[0].axial_mean)
        })
        .collect::<Vec<_>>();
    let max_meridian_slope_change = (slopes.len() >= 2).then(|| {
        slopes
            .windows(2)
            .map(|pair| (pair[1] - pair[0]).abs())
            .fold(0.0_f64, f64::max)
    });
    (stations, max_meridian_slope_change)
}

struct RadialFit {
    unique_count: usize,
    axial_span: f64,
    axial_stations: Vec<AxialRadialStation>,
    station_resolution: f64,
    max_meridian_slope_change: Option<f64>,
    mean_radius: f64,
    radius_at_origin: f64,
    slope: Option<f64>,
    max_residual: f64,
    rms_residual: f64,
    axial_min: f64,
    axial_max: f64,
    axial_mean: f64,
}

pub(super) fn source_meridian_samples(
    evidence: &MeshEvidence,
    support_ids: &[usize],
    origin: Point3,
    axis: Vector3,
) -> Option<Vec<(f64, f64)>> {
    // Use every source polygon vertex, including those not participating in
    // the support-contact interfaces. Dedup canonical position IDs exactly.
    let mut source_vertices = HashSet::<usize>::default();
    for &support in support_ids {
        for &polygon in &evidence.planar_supports.get(support)?.source_face_indices {
            for &vertex in evidence.source_polygon_vertices(polygon)? {
                source_vertices.insert(vertex);
            }
        }
    }
    if source_vertices.is_empty() {
        return None;
    }
    let mut ids = source_vertices.into_iter().collect::<Vec<_>>();
    ids.sort_unstable();
    let mut data = Vec::with_capacity(ids.len());
    for id in ids {
        let relative = *evidence.positions.get(id)? - origin;
        let axial = relative.dot(axis);
        let radial = (relative - axis * axial).magnitude();
        if !axial.is_finite() || !radial.is_finite() {
            return None;
        }
        data.push((axial, radial));
    }
    Some(data)
}

fn fit_source_vertex_radial_profile(
    evidence: &MeshEvidence,
    support_ids: &[usize],
    origin: Point3,
    axis: Vector3,
    inferred_axis_angular_uncertainty: f64,
) -> Option<RadialFit> {
    let data = source_meridian_samples(evidence, support_ids, origin, axis)?;
    let count = crate::usize_to_f64(data.len());
    let (axial_sum, radius_sum, axial_min, axial_max) = data.iter().fold(
        (0.0, 0.0, f64::INFINITY, f64::NEG_INFINITY),
        |(sum_z, sum_r, min_z, max_z), &(z, r)| (sum_z + z, sum_r + r, min_z.min(z), max_z.max(z)),
    );
    let (mean_z, mean_radius) = (axial_sum / count, radius_sum / count);
    if mean_radius <= 0.0 || !mean_radius.is_finite() {
        return None;
    }
    let (zz, zr) = data.iter().fold((0.0, 0.0), |(zz, zr), &(z, r)| {
        let dz = z - mean_z;
        (dz.mul_add(dz, zz), dz.mul_add(r - mean_radius, zr))
    });
    let axial_span = axial_max - axial_min;
    // A direction uncertainty rotates source transverse extent into the
    // axial coordinate; it is therefore an actual sampling uncertainty,
    // not evidence for many independent meridional cross-sections.
    let station_resolution = (128.0 * f64::EPSILON * mean_radius.max(axial_span))
        .max(inferred_axis_angular_uncertainty.sin().abs() * mean_radius);
    let (axial_stations, max_meridian_slope_change) =
        group_axial_stations(&data, station_resolution);
    let slope = (axial_span > 1.0e-9 * mean_radius
        && zz > f64::EPSILON * count * axial_span * axial_span)
        .then_some(zr / zz);
    let radius_at_origin = slope.unwrap_or(0.0).mul_add(-mean_z, mean_radius);
    let (max_residual, sum_squares) =
        data.iter()
            .fold((0.0_f64, 0.0_f64), |(maximum, sum), &(z, r)| {
                let expected = slope.unwrap_or(0.0).mul_add(z, radius_at_origin);
                let error = r - expected;
                (maximum.max(error.abs()), error.mul_add(error, sum))
            });
    Some(RadialFit {
        unique_count: data.len(),
        axial_span,
        axial_stations,
        station_resolution,
        max_meridian_slope_change,
        mean_radius,
        radius_at_origin,
        slope,
        max_residual,
        rms_residual: (sum_squares / count).sqrt(),
        axial_min,
        axial_max,
        axial_mean: mean_z,
    })
}

/// Measure a rotational surface candidate using source support-plane offsets
/// and *independently* fitted original-vertex radius versus axial coordinate.
///
/// A local normal-axis track is allowed to cross CAD features; this function
/// returns actual source geometry fit residuals, not a yes/no primitive
/// classification. Low fit rank or insufficient axial sample span is
/// explicitly reported, never replaced by an assumed center or cone slope.
/// The support IDs must be strictly increasing (unique), the axis must be
/// a finite unit vector, and angular uncertainty must lie in `[0, pi/2]`.
#[must_use = "the spatial fit contains the measured source residuals"]
pub fn measure_spatial_support_set(
    evidence: &MeshEvidence,
    supports: &[usize],
    axis: Vector3,
    angular_uncertainty: f64,
) -> Option<SpatialTrackEvidence> {
    if supports.windows(2).any(|pair| pair[0] >= pair[1])
        || !(0.0..=64.0 * f64::EPSILON).contains(&(axis.magnitude2() - 1.0).abs())
        || !(0.0..=std::f64::consts::FRAC_PI_2).contains(&angular_uncertainty)
    {
        return None;
    }
    let plane = fit_support_plane_placement(evidence, supports, axis)?;
    let samples = fit_source_vertex_radial_profile(
        evidence,
        supports,
        plane.origin,
        axis,
        angular_uncertainty,
    )?;
    let (projection_sum, projection_min, projection_max) = supports
        .iter()
        .map(|&id| evidence.planar_supports[id].normal.dot(axis))
        .fold(
            (0.0_f64, f64::INFINITY, f64::NEG_INFINITY),
            |(sum, min, max), value| (sum + value, min.min(value), max.max(value)),
        );
    let axial_projection = projection_sum / crate::usize_to_f64(supports.len());
    let transverse = axial_projection
        .mul_add(-axial_projection, 1.0)
        .max(0.0)
        .sqrt();
    // Signed apothem at the mean sample axial coordinate is stable even if the
    // chosen reference frame was close to a cone apex.
    let apothem_at_mean = axial_projection.mul_add(-samples.axial_mean, plane.signed_offset);
    let crosses_apex = axial_projection.mul_add(-samples.axial_min, plane.signed_offset)
        * axial_projection.mul_add(-samples.axial_max, plane.signed_offset)
        < 0.0;
    let normal_slope = (transverse > 1.0e-10
        && apothem_at_mean.abs() > 1.0e-10 * samples.mean_radius
        && !crosses_apex)
        .then_some(-axial_projection / transverse * apothem_at_mean.signum());
    let disagreement = samples
        .slope
        .zip(normal_slope)
        .map(|(observed, implied)| (observed - implied).abs());
    Some(SpatialTrackEvidence {
        axis,
        axis_origin: plane.origin,
        transverse_normal_rank_ratio: plane.rank_ratio,
        signed_support_plane_offset: plane.signed_offset,
        max_support_plane_offset_residual: plane.max_residual,
        rms_support_plane_offset_residual: plane.rms_residual,
        min_support_normal_axial_projection: projection_min,
        max_support_normal_axial_projection: projection_max,
        support_normal_axial_projection_spread: projection_max - projection_min,
        unique_source_vertex_count: samples.unique_count,
        axial_coordinate_span: samples.axial_span,
        axial_stations: samples.axial_stations,
        axial_station_resolution: samples.station_resolution,
        max_meridian_slope_change: samples.max_meridian_slope_change,
        mean_sample_radius: samples.mean_radius,
        fitted_sample_radius_at_origin: samples.radius_at_origin,
        fitted_sample_radius_slope: samples.slope,
        max_sample_radial_residual: samples.max_residual,
        rms_sample_radial_residual: samples.rms_residual,
        relative_max_sample_radial_residual: samples.max_residual / samples.mean_radius,
        normal_implied_radius_slope: normal_slope,
        radius_slope_disagreement: disagreement,
    })
}

/// Independently fit complete source geometry of one connected curvature
/// track; original native STEP face IDs and primitive types are never used.
#[must_use = "the fitted carrier evidence must be checked before use"]
pub fn measure_curvature_spatial(
    evidence: &MeshEvidence,
    track: &CurvatureTrack,
) -> Option<SpatialTrackEvidence> {
    let mut supports = track.support_indices.clone();
    supports.sort_unstable();
    supports.dedup();
    measure_spatial_support_set(
        evidence,
        &supports,
        track.seed_axis,
        track.max_seed_axis_error_radians,
    )
}

#[cfg(test)]
mod tests {
    use super::group_axial_stations;

    #[test]
    fn meridional_profile_needs_three_resolved_axial_stations() {
        let two = [(0.0, 1.0), (2.0, 0.5)];
        let (stations, curvature) = group_axial_stations(&two, 1.0e-12);
        assert_eq!(stations.len(), 2);
        assert!(curvature.is_none());

        let non_linear = [
            (-1.0, 1.0),
            (-1.0 + 1.0e-13, 1.0),
            (0.0, 2.0),
            (0.0 + 1.0e-13, 2.0),
            (1.0, 1.0),
            (1.0 + 1.0e-13, 1.0),
        ];
        let (stations, curvature) = group_axial_stations(&non_linear, 1.0e-11);
        assert_eq!(stations.len(), 3);
        assert_eq!(
            stations
                .iter()
                .map(|station| station.vertex_count)
                .sum::<usize>(),
            6
        );
        assert!((stations[1].radius_mean - 2.0).abs() < 1.0e-12);
        assert!((curvature.unwrap() - 2.0).abs() < 1.0e-11);
    }
}
