//! Source-curve evidence for one globally shared contact path.

use super::super::super::{sample_fit::fit_circle_carrier, seams::on_line};
use super::super::MeshEvidence;
use super::RegionContactPath;
use crate::{Point3, Vector3};

/// Straight-line hypothesis through both original contact path endpoints.
#[derive(Clone, Copy, Debug)]
pub struct ContactLineEvidence {
    /// Original source endpoints in directed order.
    pub endpoints: [Point3; 2],
    /// Directed unit vector from the first original vertex to the last.
    pub direction: Vector3,
    /// Original endpoint separation.
    pub chord_length: f64,
    /// Largest perpendicular deviation of any original source vertex.
    pub max_vertex_offset: f64,
    /// Largest backwards projected source edge step.
    pub max_reverse_step: f64,
    /// Adaptive exact collinearity of every original stored source point.
    pub exactly_collinear: bool,
}

/// Circular-arc hypothesis from independently sampled shared trim vertices.
#[derive(Clone, Copy, Debug)]
pub struct ContactCircleEvidence {
    /// Fitted circle center in original coordinates.
    pub center: Point3,
    /// Canonical unit normal of the source-derived candidate plane.
    pub normal: Vector3,
    /// Mean projected vertex radius about the fitted center.
    pub radius: f64,
    /// Normalized area strength of the plane-defining source point triple.
    /// Near zero indicates a poorly resolved source plane, or an
    /// independently supplied analytic circle not fitted to a source triple.
    pub plane_span_strength: f64,
    /// Largest distance of any original vertex from the candidate plane.
    pub max_plane_residual: f64,
    /// Largest radial residual of any projected original vertex.
    pub max_vertex_radial_residual: f64,
    /// Largest radial residual over every projected original chord, including
    /// its interior sagitta, not only the sampled vertices.
    pub max_chord_radial_residual: f64,
    /// Sum of signed principal angular steps in original source winding.
    pub signed_sweep_radians: f64,
    /// Sum of absolute angular steps; larger than absolute signed sweep
    /// when the source path reverses along the candidate circle.
    pub unsigned_sweep_radians: f64,
    /// Largest absolute original segment angular step.
    pub max_step_radians: f64,
    /// Number of source segments with zero projected angular progress.
    pub zero_angular_steps: usize,
}

/// Competing measurable geometric explanations of one shared source path.
#[derive(Clone, Copy, Debug)]
pub struct ContactCurveEvidence {
    /// Number of original vertex samples, excluding the repeated closing ID.
    pub source_vertex_count: usize,
    /// True if the source path ends at its first original canonical vertex.
    pub closed: bool,
    /// Total original polygonal path length.
    pub source_length: f64,
    /// Endpoint line candidate, absent for a closed path.
    pub line: Option<ContactLineEvidence>,
    /// Circular candidate, absent when samples do not resolve a plane.
    pub circle: Option<ContactCircleEvidence>,
}

impl RegionContactPath {
    /// Measure competing line and circle hypotheses once on a shared contact.
    ///
    /// No original CAD labels or primitive selection are used. Neither fit
    /// proves the original CAD curve identity or containment by both surfaces.
    #[must_use = "curve fit residuals and source coverage must be inspected"]
    pub fn measure_curves(&self, evidence: &MeshEvidence) -> Option<ContactCurveEvidence> {
        let ids = self.ordered_vertices()?;
        measure_path_curves(&evidence.positions, &ids)
    }

    pub(super) fn ordered_vertices(&self) -> Option<Vec<usize>> {
        let mut ids = Vec::with_capacity(self.edges.len() + 1);
        let mut previous = self.edges.first()?.vertices[0];
        ids.push(previous);
        for edge in &self.edges {
            if edge.vertices[0] != previous || edge.vertices[0] == edge.vertices[1] {
                return None;
            }
            previous = edge.vertices[1];
            ids.push(previous);
        }
        Some(ids)
    }
}

pub(super) fn measure_path_curves(
    points: &[Point3],
    ids: &[usize],
) -> Option<ContactCurveEvidence> {
    if ids.iter().any(|&id| id >= points.len()) {
        return None;
    }
    let source_length = ids
        .windows(2)
        .map(|pair| (points[pair[1]] - points[pair[0]]).norm())
        .sum::<f64>();
    if !source_length.is_finite() || source_length <= 0.0 {
        return None;
    }
    let closed = ids.first() == ids.last();
    let count = ids.len() - usize::from(closed);
    Some(ContactCurveEvidence {
        source_vertex_count: count,
        closed,
        source_length,
        line: (!closed).then(|| measure_line(points, ids)).flatten(),
        circle: measure_circle(points, ids, count),
    })
}

fn measure_line(points: &[Point3], ids: &[usize]) -> Option<ContactLineEvidence> {
    let start = points[*ids.first()?];
    let end = points[*ids.last()?];
    let delta = end - start;
    let chord_length = delta.norm();
    if !chord_length.is_finite() || chord_length <= 0.0 {
        return None;
    }
    let direction = delta / chord_length;
    let mut max_vertex_offset = 0.0_f64;
    let mut max_reverse_step = 0.0_f64;
    let mut previous = 0.0;
    let mut exactly_collinear = true;
    for &id in ids {
        let point = points[id];
        let relative = point - start;
        let along = relative.dot(&direction);
        max_vertex_offset = max_vertex_offset.max(relative.cross(&direction).norm());
        max_reverse_step = max_reverse_step.max(previous - along);
        previous = along;
        exactly_collinear &= on_line(start, end, point);
    }
    Some(ContactLineEvidence {
        endpoints: [start, end],
        direction,
        chord_length,
        max_vertex_offset,
        max_reverse_step,
        exactly_collinear,
    })
}

fn measure_circle(points: &[Point3], ids: &[usize], count: usize) -> Option<ContactCircleEvidence> {
    if count < 3 {
        return None;
    }
    let samples = ids.get(..count)?;
    let carrier = fit_circle_carrier(points, samples)?;
    measure_placed_circle(
        points,
        ids,
        carrier.center,
        carrier.normal,
        carrier.radius,
        carrier.plane_span_strength,
    )
}

/// Measure the original ordered path against supplied circular geometry.
/// An independent carrier can supply a circle even with only two samples.
pub(super) fn measure_placed_circle(
    points: &[Point3],
    ids: &[usize],
    center: Point3,
    normal: Vector3,
    radius: f64,
    plane_span_strength: f64,
) -> Option<ContactCircleEvidence> {
    if !super::containment::unit(normal) || !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let mut max_plane_residual = 0.0_f64;
    let mut max_radial_residual = 0.0_f64;
    let projected = ids
        .iter()
        .map(|&id| {
            let relative = *points.get(id)? - center;
            let axial = relative.dot(&normal);
            let planar = relative - normal * axial;
            max_plane_residual = max_plane_residual.max(axial.abs());
            max_radial_residual = max_radial_residual.max((planar.norm() - radius).abs());
            Some(planar)
        })
        .collect::<Option<Vec<_>>>()?;
    let mut result = ContactCircleEvidence {
        center,
        normal,
        radius,
        plane_span_strength,
        max_plane_residual,
        max_vertex_radial_residual: max_radial_residual,
        max_chord_radial_residual: max_radial_residual,
        signed_sweep_radians: 0.0,
        unsigned_sweep_radians: 0.0,
        max_step_radians: 0.0,
        zero_angular_steps: 0,
    };
    for pair in projected.windows(2) {
        let a = pair[0];
        let b = pair[1];
        let step = normal.dot(&a.cross(&b)).atan2(a.dot(&b));
        if !step.is_finite() {
            return None;
        }
        let delta = b - a;
        let length_squared = delta.norm_squared();
        if length_squared > 0.0 {
            let closest = a + delta * (-a.dot(&delta) / length_squared).clamp(0.0, 1.0);
            result.max_chord_radial_residual = result
                .max_chord_radial_residual
                .max((closest.norm() - radius).abs());
        }
        result.signed_sweep_radians += step;
        result.unsigned_sweep_radians += step.abs();
        result.max_step_radians = result.max_step_radians.max(step.abs());
        result.zero_angular_steps += usize::from(step == 0.0);
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{EvidenceOptions, RegionBoundaryIndex, analyze_mesh_evidence};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use std::f64::consts::TAU;

    #[test]
    fn shared_circle_measures_chord_error_independently_of_vertex_fit()
    -> Result<(), Box<dyn std::error::Error>> {
        let center = Point3::new(1.0e6, -1.0e6, 5.0);
        let mut positions = vec![center];
        positions.extend((0..16).map(|i| {
            let angle = TAU * f64::from(i) / 16.0;
            center + Vector3::new(angle.cos(), angle.sin(), 0.0) * 2.0
        }));
        let triangles = (1..=16_usize).map(|i| [0, i, i % 16 + 1]);
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(triangles),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let graph = RegionBoundaryIndex::new(&evidence)?.partition_contacts(&[0; 16])?;
        assert_eq!(graph.paths.len(), 1);
        let fit = graph.paths[0]
            .measure_curves(&evidence)
            .ok_or("missing curve evidence")?;
        let circle = fit.circle.ok_or("missing circle")?;
        assert!(fit.closed && fit.line.is_none());
        assert_eq!(fit.source_vertex_count, 16);
        assert!((circle.radius - 2.0).abs() < 1.0e-9);
        assert!((circle.center - center).norm() < 1.0e-9);
        assert!(circle.max_vertex_radial_residual < 1.0e-9);
        assert!(
            (circle.max_chord_radial_residual - 2.0 * (1.0 - (TAU / 32.0).cos())).abs() < 1.0e-9
        );
        assert!((circle.signed_sweep_radians.abs() - TAU).abs() < 1.0e-12);
        Ok(())
    }

    #[test]
    fn shared_collinear_subdivision_does_not_invent_a_circle()
    -> Result<(), Box<dyn std::error::Error>> {
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(2.0, 0.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                    Point3::new(1.0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 3], [1, 2, 3], [1, 0, 4], [2, 1, 4]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let graph = RegionBoundaryIndex::new(&evidence)?.partition_contacts(&[0, 0, 1, 1])?;
        let path = graph
            .paths
            .iter()
            .find(|p| p.regions == (0, Some(1)))
            .ok_or("missing contact")?;
        assert_eq!(path.edges.len(), 2);
        let fit = path
            .measure_curves(&evidence)
            .ok_or("missing curve evidence")?;
        let line = fit.line.ok_or("missing line")?;
        assert!(line.exactly_collinear && fit.circle.is_none());
        assert_eq!(line.chord_length, 2.0);
        assert_eq!(line.max_reverse_step, 0.0);
        Ok(())
    }
}
