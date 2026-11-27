//! Circular profile hypotheses from complete source-connected rotational
//! regions. Three meridian stations interpolate a circle; independently
//! resolved additional stations and all original vertices test that proposal.

use super::super::sample_fit::fit_circle_2d;
use super::{MeshEvidence, spatial::source_meridian_samples};
use crate::{InnerSpace, Point3, Vector3};
use std::ops::Range;

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

impl CircularMeridianEvidence {
    /// Observed scale limits the accuracy budget even for huge fitted radii.
    #[must_use = "the profile budget is limited by its observed domain"]
    pub const fn geometric_scale(&self) -> f64 {
        self.profile_radius.min(self.observed_meridian_span)
    }

    /// An independent fourth meridian observation and resolved nonlinearity
    /// must both survive the supplied absolute tolerance.
    #[must_use = "unresolved curvature is not a surface certificate"]
    pub fn resolves_curvature_within(&self, tolerance: f64) -> bool {
        tolerance.is_finite()
            && tolerance >= 0.0
            && self.resolved_meridian_samples >= 4
            && self.max_vertex_residual <= tolerance
            && self.linear_profile_residual_lower_bound > tolerance
            && !self.mixed_normal_orientation
    }
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
    let samples = source_meridian_samples(evidence, support_ids, origin, axis)?;
    let projection = cluster_meridian(&samples, resolution);
    let clusters = &projection.clusters;
    let reference = clusters.first()?.center;
    let local = clusters
        .iter()
        .map(|cluster| {
            (
                cluster.center.0 - reference.0,
                cluster.center.1 - reference.1,
            )
        })
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

/// Seed-anchored projected clusters, retaining exact sample provenance.
pub(super) struct MeridianCluster {
    pub center: (f64, f64),
    pub members: Range<usize>,
}

pub(super) struct MeridianClusters {
    pub order: Vec<usize>,
    pub clusters: Vec<MeridianCluster>,
}

pub(super) fn cluster_meridian(samples: &[(f64, f64)], resolution: f64) -> MeridianClusters {
    let mut order = (0..samples.len()).collect::<Vec<_>>();
    // The radial sort below canonicalizes every axial group, including
    // axial ties. Sorting by radius or ID here would repeat that work.
    order.sort_unstable_by(|&a, &b| samples[a].0.total_cmp(&samples[b].0));
    let mut result = Vec::new();
    let mut axial_begin = 0;
    while axial_begin < order.len() {
        let anchor = samples[order[axial_begin]].0;
        let axial_end = axial_begin
            + order[axial_begin..]
                .iter()
                .take_while(|&&id| samples[id].0 - anchor <= resolution)
                .count();
        let group = &mut order[axial_begin..axial_end];
        group.sort_unstable_by(|&a, &b| samples[a].1.total_cmp(&samples[b].1).then(a.cmp(&b)));
        let mut radial_begin = 0;
        while radial_begin < group.len() {
            let anchor = samples[group[radial_begin]].1;
            let radial_end = radial_begin
                + group[radial_begin..]
                    .iter()
                    .take_while(|&&id| samples[id].1 - anchor <= resolution)
                    .count();
            let members = &group[radial_begin..radial_end];
            let reference = samples[members[0]];
            let (axial, radial) = members.iter().fold((0.0, 0.0), |(za, ra), &id| {
                (
                    za + samples[id].0 - reference.0,
                    ra + samples[id].1 - reference.1,
                )
            });
            let count = crate::usize_to_f64(members.len());
            result.push(MeridianCluster {
                center: (reference.0 + axial / count, reference.1 + radial / count),
                members: axial_begin + radial_begin..axial_begin + radial_end,
            });
            radial_begin = radial_end;
        }
        axial_begin = axial_end;
    }
    MeridianClusters {
        order,
        clusters: result,
    }
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
            let Some(normal) = circular_profile_normal(centroid, axis, radial_center) else {
                continue;
            };
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

/// Unit gradient of the signed circular meridian residual, except at a
/// singular axis or profile center.
pub(super) fn circular_profile_normal(
    delta: Vector3,
    axis: Vector3,
    radial_center: f64,
) -> Option<Vector3> {
    let axial = delta.dot(axis);
    let radial = delta - axis * axial;
    let radius = radial.magnitude();
    if radial_center == 0.0 {
        return (delta.magnitude2() > 0.0).then(|| delta.normalize());
    }
    if radius <= 0.0 {
        return None;
    }
    let gradient = axis * axial + radial * ((radius - radial_center) / radius);
    (gradient.magnitude2() > 0.0).then(|| gradient.normalize())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{
        EvidenceOptions, MeridianGraphOptions, MeridianStationGraph, RegionBoundaryIndex,
        analyze_mesh_evidence,
    };
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

    fn graph_options() -> MeridianGraphOptions {
        MeridianGraphOptions {
            station_resolution: 1.0e-9,
            vertex_tolerance: 4.0e-9,
            relative_tolerance: 1.0e-4,
            max_linear_relative_chord_error: 0.25,
            max_circular_relative_chord_error: 1.0 - std::f64::consts::FRAC_1_SQRT_2,
            normal_tolerance_radians: 1.0e-6,
        }
    }

    #[test]
    fn cluster_ranges_preserve_seed_anchored_membership_and_order() {
        let samples = [
            (0.0, 3.0),
            (0.5, 1.0),
            (0.8, 1.2),
            (1.1, 1.0),
            (0.0, 1.0),
            (-0.0, 1.0),
            (2.0, 1.5),
            (2.2, 1.0),
        ];
        let projection = cluster_meridian(&samples, 1.0);
        let members = projection
            .clusters
            .iter()
            .map(|c| &projection.order[c.members.clone()])
            .collect::<Vec<_>>();
        // Neither axial nor radial clustering chains through later samples.
        assert_eq!(members, [&[1, 4, 5, 2][..], &[0], &[3, 6], &[7]]);
        assert!((projection.clusters[0].center.0 - 0.325).abs() < 1.0e-15);
        assert!((projection.clusters[0].center.1 - 1.05).abs() < 1.0e-15);
        assert!(cluster_meridian(&[], 1.0).clusters.is_empty());
    }

    #[test]
    fn unsampled_axis_skips_discovery_but_preserves_diagnostic_station_count()
    -> Result<(), Box<dyn std::error::Error>> {
        let evidence = analyze_mesh_evidence(&mesh(5)?, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let origin = Point3::new(0.123, 0.456, 0.789);
        let axis = Vector3::new(1.0, 2.0, 3.0).normalize();
        let full = MeridianStationGraph::new(&evidence, origin, axis, graph_options())
            .ok_or("full graph")?;
        let (count, fast) =
            MeridianStationGraph::for_discovery(&evidence, origin, axis, graph_options())
                .ok_or("discovery graph")?;
        assert_eq!(count, full.stations().len());
        assert!(fast.is_none());
        assert!(full.mixed_support_refinements().is_empty());
        assert!(full.discover_linear(&boundary).is_empty());
        assert_eq!(full.discover_circular(&boundary).seed_count, 0);
        let (count, fast) = MeridianStationGraph::for_discovery(
            &evidence,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            graph_options(),
        )
        .ok_or("sampled discovery graph")?;
        assert_eq!(count, 5);
        assert!(fast.is_some());
        Ok(())
    }

    #[test]
    fn negative_radial_center_cannot_displace_native_planar_regions()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{SurfaceCarrier, SurfaceRegionBudget};
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        // An exact circular fit on the positive-radial portion of a
        // meridian centered at -2. The kernel's native torus cannot
        // represent this signed branch; the source is still recoverable.
        let ring = 8;
        for station in 0..5_usize {
            let z = crate::usize_to_f64(station) - 2.0;
            let radius = 9.0_f64.mul_add(1.0, -(z * z)).sqrt() - 2.0;
            for sample in 0..ring {
                let (sin, cos) =
                    (TAU * crate::usize_to_f64(sample) / crate::usize_to_f64(ring)).sin_cos();
                positions.push(Point3::new(radius * cos, radius * sin, z));
                if station < 4 {
                    let a = station * ring + sample;
                    let b = station * ring + (sample + 1) % ring;
                    faces.push([a, b, b + ring, a + ring]);
                }
            }
        }
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let supports = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        let carrier = SurfaceCarrier::CircularRevolution {
            axial_center: Point3::new(0.0, 0.0, 0.0),
            axis: Vector3::unit_z(),
            radial_center: -2.0,
            profile_radius: 3.0,
        };
        // All original vertices are on the fitted surface; the generous
        // chord budget cannot be the reason for rejection.
        for &point in &evidence.positions {
            assert!(
                carrier
                    .measure_segment(point, point)
                    .ok_or("residual")?
                    .witnessed_max
                    < 1.0e-10
            );
        }
        assert!(
            boundary
                .measure_surface_region(
                    &supports,
                    carrier,
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        chord: 3.0
                    }
                )
                .is_none()
        );
        Ok(())
    }

    #[test]
    fn station_graph_recovers_complete_curve_without_conical_strips()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = mesh(5)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let graph = MeridianStationGraph::new(
            &evidence,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            graph_options(),
        )
        .ok_or("graph")?;
        assert_eq!(graph.stations().len(), 5);
        assert_eq!(
            graph
                .transitions()
                .iter()
                .filter(|transition| transition.repeated)
                .count(),
            4
        );
        let result = graph.discover_circular(&boundary);
        assert_eq!(result.candidates.len(), 1);
        let candidate = &result.candidates[0];
        assert_eq!(candidate.repeated_stations, 5);
        assert_eq!(candidate.source_polygons, 64);
        assert_eq!(candidate.boundary_loops, 2);
        assert!((candidate.fit.radial_center - 2.0).abs() < 1.0e-10);
        assert!((candidate.fit.profile_radius - 0.5).abs() < 1.0e-10);
        Ok(())
    }

    #[test]
    fn two_station_linear_region_is_recovered_from_independent_circle_axes()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{
            CircleAxisOptions, SurfaceCarrier, SurfaceRegionBudget, discover_circle_axes,
        };
        // Eight azimuth sectors with two different axial radii: a truncated
        // cone, not a cylinder or a planar cap.
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        for station in 0..2 {
            let axial = 2.0 * f64::from(station);
            let radius = 0.25f64.mul_add(axial, 1.0);
            for j in 0..8_usize {
                let (sin, cos) = (TAU * crate::usize_to_f64(j) / 8.0).sin_cos();
                positions.push(Point3::new(radius * cos, radius * sin, axial));
                if station == 0 {
                    let next = (j + 1) % 8;
                    faces.push([j, next, next + 8, j + 8]);
                }
            }
        }
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let found = discover_circle_axes(
            &evidence,
            CircleAxisOptions {
                vertex_tolerance: 1.0e-9,
                relative_tolerance: 1.0e-4,
                axis_position_tolerance: 8.0e-9,
                axis_angle_tolerance_radians: 1.0e-7,
            },
        )
        .ok_or("source circle axes")?;
        let placed = found
            .axes
            .iter()
            .find(|axis| {
                axis.stations.len() >= 2
                    && axis.axis.cross(Vector3::unit_z()).magnitude() < 1.0e-8
                    && axis.origin.x.hypot(axis.origin.y) < 1.0e-8
            })
            .ok_or("axis")?;
        let graph =
            MeridianStationGraph::new(&evidence, placed.origin, placed.axis, graph_options())
                .ok_or("graph")?;
        let candidates = graph.discover_linear(&boundary);
        assert_eq!(candidates.len(), 1);
        let c = &candidates[0];
        assert_eq!(c.source_polygons, 8);
        assert_eq!(c.repeated_stations, 2);
        let SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } = c.carrier
        else {
            return Err("linear carrier".into());
        };
        assert!((slope.abs() - 0.25).abs() < 1.0e-10);
        assert!(
            (radius_at_origin + slope * (Point3::new(0.0, 0.0, 0.0) - origin).dot(axis) - 1.0)
                .abs()
                < 1.0e-10
        );
        let proposal = boundary
            .measure_surface_region(
                &c.supports,
                c.carrier,
                SurfaceRegionBudget {
                    vertex: c.vertex_tolerance,
                    chord: c.chord_tolerance,
                },
            )
            .ok_or("proposal")?;
        assert!(!proposal.independently_tested_profile());
        Ok(())
    }

    #[test]
    fn mixed_coplanar_supports_can_be_refined_without_fragmenting_planar_fallback()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        // Six stations around a circular meridian include two exactly
        // horizontal chord intervals. Their toroidal facets are coplanar
        // with the straight tangent-tube extensions at both ends.
        let ring = 6;
        let sections = 6_usize;
        for section in 0..sections {
            let theta = std::f64::consts::FRAC_PI_2
                * crate::usize_to_f64(section.saturating_sub(1).min(3))
                / 3.0;
            let (st, ct) = theta.sin_cos();
            for sample in 0..ring {
                let (sp, cp) =
                    (TAU * crate::usize_to_f64(sample) / crate::usize_to_f64(ring)).sin_cos();
                let radius = 0.25f64.mul_add(cp, 1.0);
                let mut point = Point3::new(radius * ct, radius * st, 0.25 * sp);
                if section == 0 {
                    point.y -= 2.0;
                }
                if section + 1 == sections {
                    point.x -= 2.0;
                }
                positions.push(point);
                if section + 1 < sections {
                    let next = section * ring + (sample + 1) % ring;
                    let id = section * ring + sample;
                    faces.push([id, id + ring, next + ring, next]);
                }
            }
        }
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let baseline = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let graph = MeridianStationGraph::new(
            &baseline,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            graph_options(),
        )
        .ok_or("graph")?;
        let split = graph.mixed_support_refinements();
        assert!(!split.is_empty());
        let refined = baseline.refine_planar_supports(&split)?;
        assert_eq!(refined.positions, baseline.positions);
        assert_eq!(
            refined.source_polygon_vertex_indices,
            baseline.source_polygon_vertex_indices
        );
        assert!(refined.planar_supports.len() > baseline.planar_supports.len());
        let boundary = RegionBoundaryIndex::new(&refined)?;
        let restored = boundary.select_surface_partition(&[], &baseline, 100)?;
        let original = RegionBoundaryIndex::new(&baseline)?
            .partition_contacts(&baseline.source_planar_support)?;
        assert_eq!(
            restored
                .contacts
                .paths
                .iter()
                .map(|p| p.edges.len())
                .sum::<usize>(),
            original.paths.iter().map(|p| p.edges.len()).sum::<usize>()
        );
        assert_eq!(
            restored.carriers.iter().flatten().count(),
            baseline.planar_supports.len()
        );
        let graph = MeridianStationGraph::new(
            &refined,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            graph_options(),
        )
        .ok_or("refined graph")?;
        let candidates = graph.discover_circular(&boundary);
        let candidate = candidates
            .candidates
            .iter()
            .find(|c| c.source_polygons == 18)
            .ok_or("missing complete elbow")?;
        assert!((candidate.fit.radial_center - 1.0).abs() < 1.0e-10);
        assert!((candidate.fit.profile_radius - 0.25).abs() < 1.0e-10);
        Ok(())
    }

    #[test]
    fn whole_curved_profile_wins_over_individually_exact_linear_strips()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{SurfaceCarrier, SurfaceRegionBudget};
        let mesh = mesh(5)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let graph = MeridianStationGraph::new(
            &evidence,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            graph_options(),
        )
        .ok_or("graph")?;
        let curved = graph.discover_circular(&boundary);
        let c = &curved.candidates[0];
        let scale = c.fit.geometric_scale();
        let mut proposals = vec![
            boundary
                .measure_surface_region(
                    &c.supports,
                    SurfaceCarrier::from_circular_meridian(&c.fit),
                    SurfaceRegionBudget {
                        vertex: 1.0e-8,
                        chord: scale,
                    },
                )
                .ok_or("curved proposal")?,
        ];
        for c in graph.discover_linear(&boundary) {
            if let Some(proposal) = boundary.measure_surface_region(
                &c.supports,
                c.carrier,
                SurfaceRegionBudget {
                    vertex: c.vertex_tolerance,
                    chord: c.chord_tolerance,
                },
            ) {
                proposals.push(proposal);
            }
        }
        assert!(proposals.len() > 1);
        assert!(proposals[0].independently_tested_profile());
        let partition = boundary.select_surface_partition(&proposals, &evidence, 10_000)?;
        assert_eq!(partition.selected_proposals, vec![0]);
        assert_eq!(partition.consolidated_supports, c.supports.len() - 1);
        Ok(())
    }

    #[test]
    fn repeated_source_circles_supply_an_independent_meridian_axis()
    -> Result<(), Box<dyn std::error::Error>> {
        use crate::mesh_recovery::{CircleAxisOptions, discover_circle_axes};
        let mesh = mesh(5)?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let found = discover_circle_axes(
            &evidence,
            CircleAxisOptions {
                vertex_tolerance: 1.0e-9,
                relative_tolerance: 1.0e-4,
                axis_position_tolerance: 8.0e-9,
                axis_angle_tolerance_radians: 1.0e-7,
            },
        )
        .ok_or("circle discovery")?;
        let placed = found
            .axes
            .iter()
            .find(|a| {
                a.stations.len() >= 5
                    && a.axis.cross(Vector3::unit_z()).magnitude() < 1.0e-8
                    && a.origin.x.hypot(a.origin.y) < 1.0e-8
            })
            .ok_or("missing independent axis")?;
        let graph =
            MeridianStationGraph::new(&evidence, placed.origin, placed.axis, graph_options())
                .ok_or("graph")?;
        let result = graph.discover_circular(&boundary);
        assert_eq!(result.candidates.len(), 1);
        assert_eq!(result.candidates[0].source_polygons, 64);
        assert_eq!(result.candidates[0].repeated_stations, 5);
        Ok(())
    }

    /// The projected station graph joins two non-touching azimuth sectors.
    /// Neither physical source component independently resolves four stations.
    #[test]
    fn projected_station_connectivity_does_not_merge_disjoint_source_regions()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut positions = Vec::new();
        let mut faces = Vec::new();
        let ring = 5;
        for sector in 0..2 {
            let base = positions.len();
            for station in 0..3 {
                let angle = 0.2 + 0.3 * crate::usize_to_f64(station + sector * 2);
                let (sine, cosine) = angle.sin_cos();
                let radius = 0.5f64.mul_add(cosine, 2.0);
                let axial = 0.5f64.mul_add(sine, 3.0);
                for sample in 0..ring {
                    let azimuth =
                        2.0 * crate::usize_to_f64(sector) + 0.2 * crate::usize_to_f64(sample);
                    let (sine, cosine) = azimuth.sin_cos();
                    positions.push(Point3::new(radius * cosine, radius * sine, axial));
                    if station < 2 && sample + 1 < ring {
                        let start = base + station * ring + sample;
                        faces.push([start, start + 1, start + 1 + ring, start + ring]);
                    }
                }
            }
        }
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let graph = MeridianStationGraph::new(
            &evidence,
            Point3::new(0.0, 0.0, 0.0),
            Vector3::unit_z(),
            graph_options(),
        )
        .ok_or("graph")?;
        assert_eq!(graph.stations().len(), 5);
        let result = graph.discover_circular(&boundary);
        assert!(result.profile_count > 0);
        assert!(result.candidates.is_empty());
        Ok(())
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
