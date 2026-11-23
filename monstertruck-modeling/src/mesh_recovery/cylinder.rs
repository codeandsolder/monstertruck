//! Reconstruct likely circular cylinders from the two cap rings of an STL.
//!
//! An STL faceted cylinder is not geometrically identical to its continuous
//! precursor. The recovered analytic candidate carries radial and chord errors.

use super::{
    InputPolygon, MeshRecoveryError, PlanarPatch, PlanarRecoveryOptions, positions_and_polygons,
    recover_planar_mesh,
};
use crate::{Face, InnerSpace, Point3, Shell, Solid, Vector3, builder};
use monstertruck_mesh::PolygonMesh;
use monstertruck_topology::shell::ShellCondition;
use std::collections::HashSet;
use std::f64::consts::{FRAC_PI_2, TAU};

/// Limits for recognizing a closed two-ring sampled cylindrical solid.
#[derive(Clone, Copy, Debug)]
pub struct CylinderRecoveryOptions {
    /// Planar recovery options for locating both cap boundaries.
    pub planar: PlanarRecoveryOptions,
    /// Maximum sample radius deviation from one common cylinder axis.
    pub max_radial_residual: f64,
    /// Upper bound on facet-to-cylinder deviation including chord sagitta.
    pub max_surface_deviation: f64,
    /// Minimum circumference samples per end cap.
    pub min_ring_vertices: usize,
}

impl Default for CylinderRecoveryOptions {
    fn default() -> Self {
        Self {
            planar: PlanarRecoveryOptions::default(),
            max_radial_residual: 1.0e-5,
            max_surface_deviation: 0.05,
            min_ring_vertices: 12,
        }
    }
}

/// Fitted native cylinder and quantitative source evidence.
#[derive(Debug)]
pub struct RecoveredCylinder {
    /// Two planar caps and two analytical revolution faces.
    pub solid: Solid,
    /// Fitted end cap centers in axial order.
    pub cap_centers: [Point3; 2],
    /// Unit axis running from the first cap to the second.
    pub axis: Vector3,
    /// Radius of the fitted continuous cylinder.
    pub radius: f64,
    /// Axial distance between fitted cap centers.
    pub height: f64,
    /// Maximum sample radial fit error.
    pub max_radial_residual: f64,
    /// Bound on deviation from original faceted source surface.
    pub max_surface_deviation: f64,
    /// Original polygon indices for the two end caps.
    pub cap_source_face_indices: [Vec<usize>; 2],
    /// Original polygon indices for the lateral triangulation.
    pub lateral_source_face_indices: Vec<usize>,
}

/// Least-squares circle in a known plane. The relative frame avoids
/// cancellation for large world offsets.
fn fit_circle(points: &[Point3], ids: &[usize], axis: Vector3) -> Option<(Point3, f64)> {
    let start = points[*ids.first()?];
    let delta = points[*ids.get(1)?] - start;
    let projected = delta - axis * delta.dot(axis);
    let norm = projected.magnitude();
    if !norm.is_finite() || norm <= 0.0 {
        return None;
    }
    let xaxis = projected / norm;
    let yaxis = axis.cross(xaxis);
    let samples = ids
        .iter()
        .map(|&id| {
            let delta = points[id] - start;
            (delta.dot(xaxis), delta.dot(yaxis))
        })
        .collect::<Vec<_>>();
    let n = crate::usize_to_f64(samples.len());
    let mx = samples.iter().map(|&(x, _)| x).sum::<f64>() / n;
    let my = samples.iter().map(|&(_, y)| y).sum::<f64>() / n;
    let mqq = samples
        .iter()
        .map(|&(x, y)| x.mul_add(x, y * y))
        .sum::<f64>()
        / n;
    let (mut xx, mut xy, mut yy, mut qx, mut qy) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y) in &samples {
        let dx = x - mx;
        let dy = y - my;
        let dq = x.mul_add(x, y * y) - mqq;
        xx = dx.mul_add(dx, xx);
        xy = dx.mul_add(dy, xy);
        yy = dy.mul_add(dy, yy);
        qx = (dx * 0.5).mul_add(dq, qx);
        qy = (dy * 0.5).mul_add(dq, qy);
    }
    let det = xx.mul_add(yy, -(xy * xy));
    if !det.is_finite() || det <= f64::EPSILON * (xx + yy).powi(2) {
        return None;
    }
    let cx = (qx.mul_add(yy, -(qy * xy))) / det;
    let cy = (qy.mul_add(xx, -(qx * xy))) / det;
    let center = start + xaxis * cx + yaxis * cy;
    let radius = samples
        .iter()
        .map(|&(x, y)| (x - cx).hypot(y - cy))
        .sum::<f64>()
        / n;
    (radius > 0.0 && radius.is_finite()).then_some((center, radius))
}

/// Require a one-turn ordered circumference with no angular reversal.
fn ring_max_step(
    positions: &[Point3],
    ids: &[usize],
    center: Point3,
    axis: Vector3,
) -> Option<f64> {
    let mut total = 0.0;
    let mut max_step: f64 = 0.0;
    let mut sense = false;
    for i in 0..ids.len() {
        let a = positions[ids[i]] - center;
        let b = positions[ids[(i + 1) % ids.len()]] - center;
        let ra = a - axis * a.dot(axis);
        let rb = b - axis * b.dot(axis);
        let step = axis.dot(ra.cross(rb)).atan2(ra.dot(rb));
        if !step.is_finite() || step == 0.0 || step.abs() >= FRAC_PI_2 {
            return None;
        }
        if i == 0 {
            sense = step.is_sign_positive();
        } else if step.is_sign_positive() != sense {
            return None;
        }
        total += step;
        max_step = max_step.max(step.abs());
    }
    ((total.abs() - TAU).abs() <= 1.0e-6).then_some(max_step)
}

/// Construct a true analytical cylinder without retaining polygon seams.
fn build_cylinder_solid(
    cap_center: Point3,
    axis: Vector3,
    radial_start: Vector3,
    radius: f64,
    height: f64,
) -> Result<Solid, MeshRecoveryError> {
    let bottom = builder::vertex(cap_center + radial_start * radius);
    let top = builder::vertex(bottom.point() + axis * height);
    let profile = builder::line(&top, &bottom);
    let mut shell: Shell =
        builder::revolve(&profile, cap_center, axis, builder::SweepAngle::Closed, 2)
            .map_err(|_| MeshRecoveryError::CylinderConstruction)?;
    let boundaries = shell.extract_boundaries();
    if boundaries.len() != 2 {
        return Err(MeshRecoveryError::CylinderConstruction);
    }
    for boundary in boundaries {
        let cap: Face = builder::try_attach_plane([boundary.inverse()])
            .map_err(|_| MeshRecoveryError::CylinderConstruction)?;
        shell.push(cap);
    }
    Solid::new(vec![shell]).map_err(MeshRecoveryError::Topology)
}

/// Shared geometry and acceptance constraints for checking a complete source
/// mesh. Keeping cap and lateral checks in one place prevents them from using
/// inconsistent cylinder axis, radius, or error budgets.
struct CylindricalSource<'a> {
    positions: &'a [Point3],
    caps: [&'a PlanarPatch; 2],
    centers: [Point3; 2],
    axis: Vector3,
    radius: f64,
    height: f64,
    options: CylinderRecoveryOptions,
}

impl CylindricalSource<'_> {
    fn radial_distance(&self, position: Point3) -> f64 {
        let delta = position - self.centers[0];
        (delta - self.axis * delta.dot(self.axis)).magnitude()
    }

    /// Validate source facet geometry and provide original face provenance,
    /// maximum sample error, and conservative faceting error bound.
    fn validate(
        &self,
        polygons: &[InputPolygon],
        mut widest_angle: f64,
    ) -> Option<(Vec<usize>, f64, f64)> {
        let [bottom, upper] = self.caps.map(|cap| &cap.boundary_loops[0]);
        let lower_set: HashSet<usize> = bottom.iter().copied().collect();
        let upper_set: HashSet<usize> = upper.iter().copied().collect();
        if !lower_set.is_disjoint(&upper_set) {
            return None;
        }
        let mut max_residual = bottom
            .iter()
            .chain(upper.iter())
            .map(|&v| (self.radial_distance(self.positions[v]) - self.radius).abs())
            .fold(0.0_f64, f64::max);
        if !max_residual.is_finite() || max_residual > self.options.max_radial_residual {
            return None;
        }
        // Cap samples inside a polygonal disk remain inside the reconstructed
        // circle; explicit axial/radial bounds prevent lost protrusions.
        let mut cap_deviation = 0.0_f64;
        for (cap, center) in self.caps.iter().zip(self.centers) {
            for &face_id in &cap.source_face_indices {
                for &vertex in &polygons[face_id].vertices {
                    let delta = self.positions[vertex] - center;
                    let axial = delta.dot(self.axis).abs();
                    let transverse = (delta - self.axis * delta.dot(self.axis)).magnitude();
                    let overhang = (transverse - self.radius).max(0.0);
                    if axial > self.options.planar.distance_tolerance
                        || overhang > self.options.max_radial_residual
                    {
                        return None;
                    }
                    cap_deviation = cap_deviation.max(axial + overhang);
                }
            }
        }
        let cap_faces: HashSet<usize> = self
            .caps
            .iter()
            .flat_map(|cap| cap.source_face_indices.iter().copied())
            .collect();
        let mut lateral = Vec::new();
        for (face_id, polygon) in polygons.iter().enumerate() {
            if cap_faces.contains(&face_id) {
                continue;
            }
            // Intermediate axial rings are legal, but every lateral vertex
            // must lie between caps and within the same common radial fit.
            for &v in &polygon.vertices {
                let delta = self.positions[v] - self.centers[0];
                let axial = delta.dot(self.axis);
                if axial < -self.options.max_radial_residual
                    || axial > self.height + self.options.max_radial_residual
                {
                    return None;
                }
                max_residual =
                    max_residual.max((self.radial_distance(self.positions[v]) - self.radius).abs());
            }
            if max_residual > self.options.max_radial_residual {
                return None;
            }
            let centroid = polygon
                .vertices
                .iter()
                .map(|&v| self.positions[v] - self.centers[0])
                .fold(Vector3::new(0.0, 0.0, 0.0), |a, p| a + p)
                / crate::usize_to_f64(polygon.vertices.len());
            let outward = centroid - self.axis * centroid.dot(self.axis);
            if polygon.normal.dot(outward) <= 0.0 {
                return None;
            }
            for i in 0..polygon.vertices.len() {
                let a = self.positions[polygon.vertices[i]] - self.centers[0];
                let b = self.positions[polygon.vertices[(i + 1) % polygon.vertices.len()]]
                    - self.centers[0];
                let ra = a - self.axis * a.dot(self.axis);
                let rb = b - self.axis * b.dot(self.axis);
                let angle = self.axis.dot(ra.cross(rb)).atan2(ra.dot(rb)).abs();
                if !angle.is_finite() || angle >= FRAC_PI_2 {
                    return None;
                }
                widest_angle = widest_angle.max(angle);
            }
            lateral.push(face_id);
        }
        if lateral.is_empty() {
            return None;
        }
        // Facet chord sagitta bounds wall and omitted cap ring arcs. Additional
        // cap offset terms cover tilt, off-axis caps and disk overhang.
        let sagitta = self.radius * (1.0 - (widest_angle * 0.5).cos());
        let centers_delta = self.centers[1] - self.centers[0];
        let lateral_deviation = max_residual + sagitta;
        let cap_deviation =
            cap_deviation + sagitta + (centers_delta - self.axis * self.height).magnitude();
        Some((lateral, max_residual, lateral_deviation.max(cap_deviation)))
    }
}

/// Hypothesize a smooth right circular cylinder from an oriented, closed,
/// two-ring mesh, retaining quantitative approximation evidence.
///
/// Both disk caps must be coplanar patches; all side facets must lie between
/// their planes on the recovered radial carrier, including intermediate rings. A sufficiently circular, closed
/// mesh is only *consistent with* a cylinder; no STL can prove its author's
/// original intended surface. The chord error bound makes the inferred
/// smooth-surface replacement an explicit geometric choice.
///
/// Returns `Ok(None)` if the valid mesh does not meet the supplied constraints.
///
/// # Errors
/// Rejects invalid tolerances, malformed input geometry, or native topology
/// construction failure.
pub fn recover_cylindrical_solid(
    mesh: &PolygonMesh,
    options: CylinderRecoveryOptions,
) -> Result<Option<RecoveredCylinder>, MeshRecoveryError> {
    if !options.max_radial_residual.is_finite()
        || options.max_radial_residual < 0.0
        || !options.max_surface_deviation.is_finite()
        || options.max_surface_deviation < 0.0
        || options.min_ring_vertices < 6
    {
        return Err(MeshRecoveryError::InvalidCylinderOptions);
    }
    let planar = recover_planar_mesh(mesh, options.planar)?;
    if planar.shell.shell_condition() != ShellCondition::Closed {
        return Ok(None);
    }
    let caps = planar
        .patches
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            p.boundary_loops.len() == 1 && p.boundary_loops[0].len() >= options.min_ring_vertices
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let [first, second] = caps.as_slice() else {
        return Ok(None);
    };
    let first_cap = &planar.patches[*first];
    let second_cap = &planar.patches[*second];
    let axis = -first_cap.normal;
    if second_cap.normal.dot(axis) < 1.0 - 1.0e-8 {
        return Ok(None);
    }
    let positions = &planar.canonical_positions;
    let bottom = &first_cap.boundary_loops[0];
    let upper = &second_cap.boundary_loops[0];
    let Some((center0, radius0)) = fit_circle(positions, bottom, axis) else {
        return Ok(None);
    };
    let Some((center1, radius1)) = fit_circle(positions, upper, axis) else {
        return Ok(None);
    };
    let span = center1 - center0;
    let height = span.dot(axis);
    if !height.is_finite()
        || height <= 0.0
        || (span - axis * height).magnitude() > options.max_radial_residual
    {
        return Ok(None);
    }

    let radius = f64::midpoint(radius0, radius1);
    let Some(lower_step) = ring_max_step(positions, bottom, center0, axis) else {
        return Ok(None);
    };
    let Some(upper_step) = ring_max_step(positions, upper, center1, axis) else {
        return Ok(None);
    };
    let (original_positions, polygons) =
        positions_and_polygons(mesh, options.planar.distance_tolerance)?;
    debug_assert_eq!(original_positions, *positions);
    let source = CylindricalSource {
        positions,
        caps: [first_cap, second_cap],
        centers: [center0, center1],
        axis,
        radius,
        height,
        options,
    };
    let Some((lateral, max_residual, deviation)) =
        source.validate(&polygons, lower_step.max(upper_step))
    else {
        return Ok(None);
    };
    if deviation > options.max_surface_deviation {
        return Ok(None);
    }
    let radial_start = positions[bottom[0]] - center0;
    let radial_start = radial_start - axis * radial_start.dot(axis);
    let radial_length = radial_start.magnitude();
    if radial_length == 0.0 || !radial_length.is_finite() {
        return Ok(None);
    }
    let solid = build_cylinder_solid(center0, axis, radial_start / radial_length, radius, height)?;
    Ok(Some(RecoveredCylinder {
        solid,
        cap_centers: [center0, center1],
        axis,
        radius,
        height,
        max_radial_residual: max_residual,
        max_surface_deviation: deviation,
        cap_source_face_indices: [
            first_cap.source_face_indices.clone(),
            second_cap.source_face_indices.clone(),
        ],
        lateral_source_face_indices: lateral,
    }))
}

#[cfg(test)]
mod tests;
