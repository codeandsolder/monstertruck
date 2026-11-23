//! Replace a polygonal through-hole with an analytic cylindrical void inside
//! a mixed planar/cylindrical B-rep, preserving the exterior planar patches.
//!
//! The recognition operates on topology and geometry together: two planar
//! annular cap faces, a connected faceted wall between their inner boundaries,
//! and the remaining planar exterior. All face/edge reconstruction uses native
//! Monstertruck topology, with circular trims shared by both carriers.

use super::cylinder::{CylinderRecoveryOptions, build_cylinder_side, fit_circle, ring_max_step};
use super::{
    HashMap, InputPolygon, MeshRecoveryError, PlanarPatch, RecoveredPlanarMesh, edge_incidence,
    edge_key, positions_and_polygons, recover_planar_mesh, seams,
};
use crate::{Face, InnerSpace, KERNEL_EPSILON, Point3, Shell, Solid, Vector3, Wire};
use monstertruck_mesh::PolygonMesh;
use monstertruck_topology::shell::ShellCondition;
use std::collections::{HashSet, VecDeque};

/// Native analytic B-rep recovered from a faceted through-hole solid.
#[derive(Debug)]
pub struct RecoveredCylindricalThroughHole {
    /// Closed native solid retaining the input's exterior planar faces.
    pub solid: Solid,
    /// Circle center on each of the two original planar face carriers.
    pub cap_centers: [Point3; 2],
    /// Unit axis from the first end plane to the second.
    pub axis: Vector3,
    /// Fitted hole radius.
    pub radius: f64,
    /// End-cap separation.
    pub height: f64,
    /// Max radial deviation of all source wall vertices from the fitted carrier.
    pub max_radial_residual: f64,
    /// Bound on faceted wall and trim deviation from the smooth circle.
    pub max_surface_deviation: f64,
    /// Source polygon indices of each planar end-cap annulus.
    pub cap_source_face_indices: [Vec<usize>; 2],
    /// Source polygon indices belonging to the removed faceted hole wall.
    pub cylindrical_source_face_indices: Vec<usize>,
}

/// The wall is a connected portion of the source closed manifold. Validate
/// that every boundary edge of the removed portion is adjacent to exactly
/// one of the two inner cap loops and that both loops are fully covered.
fn wall_has_exact_cap_boundary(
    polygons: &[InputPolygon],
    positions: &[Point3],
    selected: &HashSet<usize>,
    caps: [&PlanarPatch; 2],
) -> Result<bool, MeshRecoveryError> {
    let mut working = polygons.to_vec();
    seams::split_source_t_junctions(positions, &mut working);
    let (incidence, _) = edge_incidence(&working)?;
    let expected = caps
        .iter()
        .flat_map(|cap| {
            let cycle = &cap.boundary_loops[1];
            (0..cycle.len()).map(|i| edge_key(cycle[i], cycle[(i + 1) % cycle.len()]).0)
        })
        .collect::<HashSet<_>>();
    let cap_polygons = caps
        .iter()
        .flat_map(|cap| cap.source_face_indices.iter().copied())
        .collect::<HashSet<_>>();
    let mut found = HashSet::new();
    let mut wall_adjacency = HashMap::<usize, Vec<usize>>::default();
    for (key, uses) in &incidence {
        let adjacent = uses
            .iter()
            .filter_map(|(face, _)| selected.contains(face).then_some(*face))
            .collect::<Vec<_>>();
        match adjacent.as_slice() {
            [a, b] => {
                wall_adjacency.entry(*a).or_default().push(*b);
                wall_adjacency.entry(*b).or_default().push(*a);
            }
            [_wall] => {
                if uses.len() != 2
                    || !uses.iter().any(|(face, _)| cap_polygons.contains(face))
                    || !expected.contains(key)
                {
                    return Ok(false);
                }
                found.insert(*key);
            }
            [] => {}
            _ => return Ok(false),
        }
    }
    if found != expected {
        return Ok(false);
    }
    // A second internally closed faceted component with the same fitted
    // radius must not silently disappear when replacing the source wall.
    let Some(&first) = selected.iter().min() else {
        return Ok(false);
    };
    let mut reached = HashSet::from([first]);
    let mut pending = VecDeque::from([first]);
    while let Some(face) = pending.pop_front() {
        for &neighbor in wall_adjacency.get(&face).into_iter().flatten() {
            if reached.insert(neighbor) {
                pending.push_back(neighbor);
            }
        }
    }
    Ok(reached == *selected)
}

/// Fitted end cap evidence shared by the cylindrical wall and trim builders.
struct HoleCaps {
    patch_indices: [usize; 2],
    centers: [Point3; 2],
    axis: Vector3,
    radius: f64,
    height: f64,
    max_ring_residual: f64,
}

/// Geometric and topological evidence for the faceted wall being replaced.
struct WallEvidence {
    source_faces: HashSet<usize>,
    removed_patches: HashSet<usize>,
    max_radial_residual: f64,
    max_surface_deviation: f64,
}

/// Classify faceted cylindrical wall polygons by their common circle/axis
/// carrier, inward normals, and a complete two-ring manifold boundary.
struct HoleWall<'a> {
    positions: &'a [Point3],
    polygons: &'a [InputPolygon],
    caps: [&'a PlanarPatch; 2],
    centers: [Point3; 2],
    axis: Vector3,
    radius: f64,
    height: f64,
    options: CylinderRecoveryOptions,
}

impl HoleWall<'_> {
    fn candidate_wall(&self, polygon: &InputPolygon, cap_face: bool) -> Option<(f64, f64)> {
        if cap_face {
            return None;
        }
        let mut max_radial = 0.0_f64;
        let mut max_step = 0.0_f64;
        let mut center = Vector3::new(0.0, 0.0, 0.0);
        for (i, &index) in polygon.vertices.iter().enumerate() {
            let delta = self.positions[index] - self.centers[0];
            let axial = delta.dot(self.axis);
            if axial < -self.options.max_radial_residual
                || axial > self.height + self.options.max_radial_residual
            {
                return None;
            }
            let radial = delta - self.axis * axial;
            let radius = radial.magnitude();
            max_radial = max_radial.max((radius - self.radius).abs());
            if max_radial > self.options.max_radial_residual {
                return None;
            }
            let other = self.positions[polygon.vertices[(i + 1) % polygon.vertices.len()]]
                - self.centers[0];
            let next = other - self.axis * other.dot(self.axis);
            let angle = self
                .axis
                .dot(radial.cross(next))
                .atan2(radial.dot(next))
                .abs();
            if !angle.is_finite() || angle >= std::f64::consts::FRAC_PI_2 {
                return None;
            }
            max_step = max_step.max(angle);
            center += delta;
        }
        let centroid = center / crate::usize_to_f64(polygon.vertices.len());
        let radial = centroid - self.axis * centroid.dot(self.axis);
        if polygon.normal.dot(radial) >= 0.0 {
            return None;
        }
        Some((max_radial, max_step))
    }

    /// Every removed polygon must be locally consistent with the hole wall;
    /// global oriented adjacency and cap boundary coverage are checked too.
    fn validate(
        &self,
        patches: &[PlanarPatch],
        max_ring_residual: f64,
    ) -> Result<Option<WallEvidence>, MeshRecoveryError> {
        let cap_source = self
            .caps
            .iter()
            .flat_map(|cap| cap.source_face_indices.iter().copied())
            .collect::<HashSet<_>>();
        let mut selected = HashSet::new();
        let mut max_radial = max_ring_residual;
        let mut max_step = 0.0_f64;
        for (index, polygon) in self.polygons.iter().enumerate() {
            if let Some((radial, step)) = self.candidate_wall(polygon, cap_source.contains(&index))
            {
                selected.insert(index);
                max_radial = max_radial.max(radial);
                max_step = max_step.max(step);
            }
        }
        if selected.is_empty()
            || !wall_has_exact_cap_boundary(self.polygons, self.positions, &selected, self.caps)?
        {
            return Ok(None);
        }
        // Patches are the only reconstructible planar face units. A partial
        // classification would leave a hole in the original shell: refuse it.
        let mut removed = HashSet::new();
        for (index, patch) in patches.iter().enumerate() {
            let part = patch
                .source_face_indices
                .iter()
                .filter(|face| selected.contains(face))
                .count();
            if part > 0 {
                if part != patch.source_face_indices.len() {
                    return Ok(None);
                }
                removed.insert(index);
            }
        }
        let sagitta = self.radius * (1.0 - (max_step * 0.5).cos());
        let radial_error = max_radial + sagitta;
        let cap_offset = (self.centers[1] - self.centers[0] - self.axis * self.height).magnitude();
        let deviation = radial_error + cap_offset;
        if deviation > self.options.max_surface_deviation {
            return Ok(None);
        }
        Ok(Some(WallEvidence {
            source_faces: selected,
            removed_patches: removed,
            max_radial_residual: max_radial,
            max_surface_deviation: deviation,
        }))
    }
}

fn caps_and_axis(mesh: &RecoveredPlanarMesh, options: CylinderRecoveryOptions) -> Option<HoleCaps> {
    let caps = mesh
        .patches
        .iter()
        .enumerate()
        .filter(|(_, patch)| {
            patch.boundary_loops.len() == 2
                && patch.boundary_loops[1].len() >= options.min_ring_vertices
        })
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    let [a, b] = caps.as_slice() else {
        return None;
    };
    let first = &mesh.patches[*a];
    let second = &mesh.patches[*b];
    let axis = -first.normal;
    if second.normal.dot(axis) < 1.0 - options.planar.angular_tolerance_radians {
        return None;
    }
    let positions = &mesh.canonical_positions;
    let loop_a = &first.boundary_loops[1];
    let loop_b = &second.boundary_loops[1];
    let (center_a, radius_a) = fit_circle(positions, loop_a, axis)?;
    let (center_b, radius_b) = fit_circle(positions, loop_b, axis)?;
    let difference = center_b - center_a;
    let height = difference.dot(axis);
    if !height.is_finite()
        || height <= 0.0
        || (difference - axis * height).magnitude() > options.max_radial_residual
    {
        return None;
    }
    let radius = f64::midpoint(radius_a, radius_b);
    let mut max_ring_residual = 0.0_f64;
    for (cycle, center) in [(loop_a, center_a), (loop_b, center_b)] {
        ring_max_step(positions, cycle, center, axis)?;
        let residual = cycle
            .iter()
            .map(|&v| {
                let delta = positions[v] - center;
                let radial = delta - axis * delta.dot(axis);
                (radial.magnitude() - radius).abs()
            })
            .fold(0.0_f64, f64::max);
        if !residual.is_finite() || residual > options.max_radial_residual {
            return None;
        }
        max_ring_residual = max_ring_residual.max(residual);
    }
    Some(HoleCaps {
        patch_indices: [*a, *b],
        centers: [center_a, center_b],
        axis,
        radius,
        height,
        max_ring_residual,
    })
}

fn reconstruct_shell(
    planar: &RecoveredPlanarMesh,
    cap_indices: [usize; 2],
    removed: &HashSet<usize>,
    side: Shell,
    centers: [Point3; 2],
    axis: Vector3,
    height: f64,
) -> Result<Option<Solid>, MeshRecoveryError> {
    let mut circular_wires: [Option<Wire>; 2] = [None, None];
    for wire in side.extract_boundaries() {
        let Some(point) = wire
            .front_vertex()
            .map(monstertruck_topology::Vertex::point)
        else {
            return Ok(None);
        };
        let axial = (point - centers[0]).dot(axis);
        let idx = if axial.abs() <= KERNEL_EPSILON {
            0
        } else if (axial - height).abs() <= KERNEL_EPSILON {
            1
        } else {
            return Ok(None);
        };
        if circular_wires[idx].replace(wire.inverse()).is_some() {
            return Ok(None);
        }
    }
    let [Some(first_wire), Some(second_wire)] = circular_wires else {
        return Ok(None);
    };
    let circular_wires = [first_wire, second_wire];
    let mut result_faces = Vec::new();
    for (index, face) in planar.shell.iter().enumerate() {
        if index == cap_indices[0] || index == cap_indices[1] {
            let cap_number = usize::from(index == cap_indices[1]);
            let outer = face.boundaries()[0].clone();
            let replaced = Face::new(
                vec![outer, circular_wires[cap_number].clone()],
                face.oriented_surface(),
            )?;
            result_faces.push(replaced);
        } else if !removed.contains(&index) {
            result_faces.push(face.clone());
        }
    }
    result_faces.extend(side);
    let solid = Solid::new(vec![result_faces.into_iter().collect()])?;
    if !solid.is_geometric_consistent() {
        return Ok(None);
    }
    Ok(Some(solid))
}

/// Recognize one true cylindrical through-hole in an otherwise planar
/// manifold, reconstructing exact shared circular side/cap trims.
///
/// This routine handles arbitrary outer planar contours, and triangulation
/// of both annular cap surfaces and the curved hole wall. It preserves other
/// planar faces and their shared exterior topology. It makes an explicitly
/// bounded approximation from a faceted source wall to an analytic cylinder.
///
/// The original mesh must form one consistently oriented, closed, connected
/// surface. It must have exactly two end-cap faces with one inner ring each;
/// wall facets with intermediate axial rings are valid. Mixed cylindrical
/// intersections that break this two-cap topology are a separate problem.
///
/// Returns `Ok(None)` for valid source meshes that do not fit.
///
/// # Errors
/// Rejects invalid recovery options, invalid source mesh geometry, and native
/// topology construction failures.
pub fn recover_cylindrical_through_hole(
    mesh: &PolygonMesh,
    options: CylinderRecoveryOptions,
) -> Result<Option<RecoveredCylindricalThroughHole>, MeshRecoveryError> {
    if !options.max_radial_residual.is_finite()
        || options.max_radial_residual < 0.0
        || !options.max_surface_deviation.is_finite()
        || options.max_surface_deviation < 0.0
        || options.min_ring_vertices < 6
    {
        return Err(MeshRecoveryError::InvalidCylinderOptions);
    }
    let planar = recover_planar_mesh(mesh, options.planar)?;
    if planar.shell.shell_condition() != ShellCondition::Closed || !planar.shell.is_connected() {
        return Ok(None);
    }
    let Some(fitted) = caps_and_axis(&planar, options) else {
        return Ok(None);
    };
    let HoleCaps {
        patch_indices: caps,
        centers,
        axis,
        radius,
        height,
        max_ring_residual,
    } = fitted;
    let positions = &planar.canonical_positions;
    let (original, polygons) = positions_and_polygons(mesh, options.planar.distance_tolerance)?;
    debug_assert_eq!(original, *positions);
    let wall = HoleWall {
        positions,
        polygons: &polygons,
        caps: [&planar.patches[caps[0]], &planar.patches[caps[1]]],
        centers,
        axis,
        radius,
        height,
        options,
    };
    let Some(evidence) = wall.validate(&planar.patches, max_ring_residual)? else {
        return Ok(None);
    };
    let radial = positions[planar.patches[caps[0]].boundary_loops[1][0]] - centers[0];
    let radial = radial - axis * radial.dot(axis);
    if radial.magnitude() <= KERNEL_EPSILON {
        return Ok(None);
    }
    let side = build_cylinder_side(centers[0], axis, radial.normalize(), radius, height, false)?;
    let Some(solid) = reconstruct_shell(
        &planar,
        caps,
        &evidence.removed_patches,
        side,
        centers,
        axis,
        height,
    )?
    else {
        return Ok(None);
    };
    let mut source_faces = evidence.source_faces.into_iter().collect::<Vec<_>>();
    source_faces.sort_unstable();
    Ok(Some(RecoveredCylindricalThroughHole {
        solid,
        cap_centers: centers,
        axis,
        radius,
        height,
        max_radial_residual: evidence.max_radial_residual,
        max_surface_deviation: evidence.max_surface_deviation,
        cap_source_face_indices: [
            planar.patches[caps[0]].source_face_indices.clone(),
            planar.patches[caps[1]].source_face_indices.clone(),
        ],
        cylindrical_source_face_indices: source_faces,
    }))
}

#[cfg(test)]
mod tests;
