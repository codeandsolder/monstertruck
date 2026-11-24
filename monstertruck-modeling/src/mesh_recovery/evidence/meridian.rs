//! Circular profile hypotheses from complete source-connected rotational
//! regions. Three meridian stations interpolate a circle; independently
//! resolved additional stations and all original vertices test that proposal.

use super::super::sample_fit::fit_circle_2d;
use super::{MeshEvidence, spatial::source_meridian_samples};
use crate::{InnerSpace, Point3, Vector3};

/// Independent measured circular meridian of a rotational carrier.
///
/// A positive radial center describes a toroidal profile; zero corresponds
/// to a sphere. The caller chooses budgets and checks source topology.
#[derive(Clone, Copy, Debug)]
pub struct CircularMeridianEvidence {
    /// Candidate point on the axis at the meridian circle center's axial coordinate.
    pub axial_center: Point3,
    /// Fixed unit axis from independently placed neighboring observations.
    pub axis: Vector3,
    /// Radial coordinate of the meridian circle center.
    pub radial_center: f64,
    /// Radius of the fitted meridian circle.
    pub profile_radius: f64,
    /// Diagonal of the observed meridian sample bounding box.
    pub observed_meridian_span: f64,
    /// Lower bound on the maximum residual of any straight-line profile.
    /// Obtained from the minimum width of a triangle of resolved samples.
    pub linear_profile_residual_lower_bound: f64,
    /// Number of distinct resolved meridian clusters used for fitting.
    pub resolved_meridian_samples: usize,
    /// Number of unique original source vertices measured.
    pub source_vertices: usize,
    /// Maximum distance residual of every source vertex to the circle.
    pub max_vertex_residual: f64,
    /// RMS distance residual over all source vertices.
    pub rms_vertex_residual: f64,
    /// Maximum measured angular deviation between observed support normals
    /// and analytical normals at the source polygon centroids, sign independent.
    pub max_normal_angle_radians: f64,
    /// Whether both signs of oriented normal agreement occur among source
    /// polygons with resolved analytic normals.
    pub mixed_normal_orientation: bool,
    /// Number of source polygons with a resolved analytic normal.
    pub normal_samples: usize,
}

/// Fit a circular meridian at a fixed independently recovered axis.
///
/// Source samples are clustered within the supplied length resolution in
/// both axial and radial coordinates. Fit cluster means with equal weight
/// (dense azimuth tessellation must not dominate), then measure every vertex.
/// No STEP labels, surface-family hints, or topology modifications are used.
/// A three-cluster interpolation is reported as a hypothesis, not evidence of
/// independent curvature. Four or more clusters can test the model.
#[must_use = "inspect independent sample count and full-source residuals"]
pub fn measure_circular_meridian(
    evidence: &MeshEvidence,
    support_ids: &[usize],
    origin: Point3,
    axis: Vector3,
    resolution: f64,
) -> Option<CircularMeridianEvidence> {
    if !resolution.is_finite()
        || resolution < 0.0
        || !(0.0..=64.0 * f64::EPSILON).contains(&(axis.magnitude2() - 1.0).abs())
        || support_ids.windows(2).any(|ids| ids[0] >= ids[1])
    {
        return None;
    }
    let mut samples = source_meridian_samples(evidence, support_ids, origin, axis)?;
    let clusters = cluster_meridian(&mut samples, resolution);
    let reference = *clusters.first()?;
    let local = clusters
        .iter()
        .map(|&(z, r)| (z - reference.0, r - reference.1))
        .collect::<Vec<_>>();
    let ([axial, radial], profile_radius) = fit_circle_2d(&local)?;
    let (observed_meridian_span, linear_profile_residual_lower_bound) =
        meridian_extent(&samples, &local);
    let axial_center = axial + reference.0;
    let radial_center = radial + reference.1;
    let (maximum, sum_squares) =
        samples
            .iter()
            .fold((0.0_f64, 0.0_f64), |(maximum, sum), &(z, r)| {
                let error = ((z - axial_center).hypot(r - radial_center) - profile_radius).abs();
                (maximum.max(error), error.mul_add(error, sum))
            });
    let center = origin + axis * axial_center;
    let normal = measure_normals(evidence, support_ids, center, axis, radial_center)?;
    Some(CircularMeridianEvidence {
        axial_center: center,
        axis,
        radial_center,
        profile_radius,
        observed_meridian_span,
        linear_profile_residual_lower_bound,
        resolved_meridian_samples: clusters.len(),
        source_vertices: samples.len(),
        max_vertex_residual: maximum,
        rms_vertex_residual: (sum_squares / crate::usize_to_f64(samples.len())).sqrt(),
        max_normal_angle_radians: normal.max_angle,
        mixed_normal_orientation: normal.positive && normal.negative,
        normal_samples: normal.samples,
    })
}

/// Measure the sampled domain and an independent lower bound for the error
/// of a straight line. Any strip enclosing three points has width at least
/// twice their triangle area divided by its longest side.
fn meridian_extent(samples: &[(f64, f64)], local: &[(f64, f64)]) -> (f64, f64) {
    let (z_min, z_max, r_min, r_max) = samples.iter().fold(
        (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ),
        |(z_min, z_max, r_min, r_max), &(z, r)| {
            (z_min.min(z), z_max.max(z), r_min.min(r), r_max.max(r))
        },
    );
    let span = (z_max - z_min).hypot(r_max - r_min);
    let Some(&(z, r)) = local
        .iter()
        .max_by(|a, b| a.0.hypot(a.1).total_cmp(&b.0.hypot(b.1)))
    else {
        return (span, 0.0);
    };
    let cross = |&(u, v): &(f64, f64)| z.mul_add(v, -r * u).abs();
    let Some(point) = local.iter().max_by(|a, b| cross(a).total_cmp(&cross(b))) else {
        return (span, 0.0);
    };
    let longest = z
        .hypot(r)
        .max(point.0.hypot(point.1))
        .max((point.0 - z).hypot(point.1 - r));
    (
        span,
        if longest > 0.0 {
            cross(point) / (2.0 * longest)
        } else {
            0.0
        },
    )
}

fn cluster_meridian(samples: &mut [(f64, f64)], resolution: f64) -> Vec<(f64, f64)> {
    samples.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let mut result = Vec::new();
    let mut axial_begin = 0;
    while axial_begin < samples.len() {
        let anchor = samples[axial_begin].0;
        let axial_end = axial_begin
            + samples[axial_begin..].partition_point(|sample| sample.0 - anchor <= resolution);
        let group = &mut samples[axial_begin..axial_end];
        group.sort_unstable_by(|a, b| a.1.total_cmp(&b.1));
        let mut radial_begin = 0;
        while radial_begin < group.len() {
            let anchor = group[radial_begin].1;
            let radial_end = radial_begin
                + group[radial_begin..].partition_point(|sample| sample.1 - anchor <= resolution);
            let members = &group[radial_begin..radial_end];
            let (axial, radial) = members
                .iter()
                .fold((0.0, 0.0), |(za, ra), &(z, r)| (za + z, ra + r));
            let count = crate::usize_to_f64(members.len());
            result.push((axial / count, radial / count));
            radial_begin = radial_end;
        }
        axial_begin = axial_end;
    }
    result
}

#[derive(Default)]
struct NormalEvidence {
    max_angle: f64,
    positive: bool,
    negative: bool,
    samples: usize,
}

fn measure_normals(
    evidence: &MeshEvidence,
    supports: &[usize],
    center: Point3,
    axis: Vector3,
    radial_center: f64,
) -> Option<NormalEvidence> {
    let mut result = NormalEvidence::default();
    for &id in supports {
        let support = evidence.planar_supports.get(id)?;
        for &polygon in &support.source_face_indices {
            let vertices = evidence.source_polygon_vertices(polygon)?;
            let mut centroid = Vector3::new(0.0, 0.0, 0.0);
            for &vertex in vertices {
                centroid += *evidence.positions.get(vertex)? - center;
            }
            centroid /= crate::usize_to_f64(vertices.len());
            let axial = centroid.dot(axis);
            let radial = centroid - axis * axial;
            let radius = radial.magnitude();
            if radius <= 0.0 {
                continue;
            }
            let gradient = axis * axial + radial * ((radius - radial_center) / radius);
            if gradient.magnitude2() <= 0.0 {
                continue;
            }
            let normal = gradient.normalize();
            let dot = normal.dot(support.normal);
            result.positive |= dot > 0.0;
            result.negative |= dot < 0.0;
            result.max_angle = result
                .max_angle
                .max(normal.cross(support.normal).magnitude().atan2(dot.abs()));
            result.samples += 1;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{EvidenceOptions, analyze_mesh_evidence};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use std::f64::consts::TAU;

    fn mesh(stations: usize) -> Result<PolygonMesh, Box<dyn std::error::Error>> {
        let ring = 16;
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        for station in 0..stations {
            let angle = 0.2 + 0.3 * crate::usize_to_f64(station);
            let (sine, cosine) = angle.sin_cos();
            let radius = 0.5f64.mul_add(cosine, 2.0);
            let axial = 0.5f64.mul_add(sine, 3.0);
            for sample in 0..ring {
                let (sine, cosine) =
                    (TAU * crate::usize_to_f64(sample) / crate::usize_to_f64(ring)).sin_cos();
                positions.push(Point3::new(radius * cosine, radius * sine, axial));
                if station + 1 < stations {
                    let start = station * ring + sample;
                    let next = station * ring + (sample + 1) % ring;
                    faces.push([start, next, next + ring, start + ring]);
                }
            }
        }
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?)
    }

    #[test]
    fn recover_circular_meridian_from_independent_source_stations()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = mesh(5)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let supports = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        let fit = measure_circular_meridian(
            &evidence,
            &supports,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            1.0e-10,
        )
        .ok_or("missing circular profile")?;
        assert_eq!(fit.resolved_meridian_samples, 5);
        assert_eq!(fit.source_vertices, 80);
        assert!((fit.axial_center.z - 3.0).abs() < 1.0e-10);
        assert!((fit.radial_center - 2.0).abs() < 1.0e-10);
        assert!((fit.profile_radius - 0.5).abs() < 1.0e-10);
        assert!(fit.max_vertex_residual < 1.0e-10);
        assert!(fit.observed_meridian_span > 0.5);
        assert!(fit.linear_profile_residual_lower_bound > 0.02);
        assert!(!fit.mixed_normal_orientation);
        assert_eq!(fit.normal_samples, 64);
        Ok(())
    }

    #[test]
    fn two_meridian_stations_do_not_resolve_curvature() -> Result<(), Box<dyn std::error::Error>> {
        let mesh = mesh(2)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let supports = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        assert!(
            measure_circular_meridian(
                &evidence,
                &supports,
                Point3::new(0.0, 0.0, 0.0),
                Vector3::unit_z(),
                1.0e-10
            )
            .is_none()
        );
        Ok(())
    }
}
