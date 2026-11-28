//! A curved proposal can interpolate all vertices of a planar chord while
//! claiming a facet that actually belongs to an adjacent observed plane.
//! Correct that ownership only at failed contacts and only when the complete
//! original support lies on the plane at native intersection precision.

use super::{AnalyticRecoveryError, AnalyticRecoveryOptions, NativeAssembly};
use crate::{
    InnerSpace, KERNEL_EPSILON,
    mesh_recovery::{MeshEvidence, NativeTrimOptions, RegionBoundaryIndex, SurfaceCarrier},
};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// One source-backed carrier ownership correction, without moving source points.
#[derive(Clone, Copy, Debug)]
pub struct CoplanarSupportReassignment {
    /// Original planar support identifier.
    pub support: usize,
    /// Fitted region that previously owned the support.
    pub old_region: usize,
    /// Connected planar region in the corrected partition.
    pub new_region: usize,
    /// Maximum vertex-to-plane distance in either source position table.
    pub max_plane_residual: f64,
}

fn plane_residual(
    evidence: &MeshEvidence,
    assembly: &NativeAssembly,
    support: usize,
    target: usize,
    options: NativeTrimOptions,
) -> Option<f64> {
    let partition = &assembly.reconciled.partition;
    if matches!(
        partition.carriers[partition.region_by_support[support]],
        Some(SurfaceCarrier::Plane { .. })
    ) {
        return None;
    }
    let Some(SurfaceCarrier::Plane { origin, normal }) = partition.carriers[target] else {
        return None;
    };
    let plane = &evidence.planar_supports[support];
    let angle = normal
        .cross(plane.normal)
        .magnitude()
        .atan2(normal.dot(plane.normal));
    if angle + plane.max_normal_angle_radians > options.angular_resolution {
        return None;
    }
    let mut residual: f64 = 0.0;
    for &polygon in &plane.source_face_indices {
        for &vertex in evidence.source_polygon_vertices(polygon)? {
            for positions in [&evidence.positions, &assembly.reconciled.positions] {
                let distance = (positions[vertex] - origin).dot(normal).abs();
                if !distance.is_finite() || distance > options.intersection_tolerance {
                    return None;
                }
                residual = residual.max(distance);
            }
        }
    }
    Some(residual)
}

fn proposals(
    evidence: &MeshEvidence,
    assembly: &NativeAssembly,
    options: NativeTrimOptions,
) -> BTreeMap<usize, (usize, f64)> {
    let mut queue = VecDeque::new();
    for (trim, native) in assembly.trims.trims.iter().zip(&assembly.native.trims) {
        if native.is_ok() {
            continue;
        }
        let Some(right) = trim.regions.1 else {
            continue;
        };
        for edge in &trim.source_edges {
            let Some(exterior) = edge.exterior_polygon else {
                continue;
            };
            queue.push_back((evidence.source_planar_support[edge.interior_polygon], right));
            queue.push_back((evidence.source_planar_support[exterior], trim.regions.0));
        }
    }
    let mut visited = BTreeSet::new();
    let mut selected = BTreeMap::new();
    while let Some((support, target)) = queue.pop_front() {
        if !visited.insert((support, target)) {
            continue;
        }
        let Some(residual) = plane_residual(evidence, assembly, support, target, options) else {
            continue;
        };
        if let Some(&(previous, _)) = selected.get(&support) {
            if previous != target {
                // Competing plane ownership is not uniquely witnessed.
                return BTreeMap::new();
            }
            continue;
        }
        selected.insert(support, (target, residual));
        for &interface in &evidence.support_interface_indices[support] {
            let [a, b] = evidence.support_interfaces[interface].supports;
            queue.push_back((if a == support { b } else { a }, target));
        }
    }
    selected
}

pub(super) fn failed_edges(assembly: &NativeAssembly) -> BTreeSet<[usize; 2]> {
    assembly
        .trims
        .trims
        .iter()
        .zip(&assembly.native.trims)
        .filter(|(_, native)| native.is_err())
        .flat_map(|(trim, _)| {
            trim.source_edges.iter().map(|edge| {
                let mut vertices = edge.vertices;
                vertices.sort_unstable();
                vertices
            })
        })
        .collect()
}

pub(super) fn valid_regions(assembly: &NativeAssembly) -> BTreeSet<usize> {
    assembly
        .native
        .regions
        .iter()
        .filter(|region| region.faces.is_ok() && !region.requires_chart_cut)
        .map(|region| region.region)
        .collect()
}

pub(super) fn repair(
    evidence: &MeshEvidence,
    assembly: &mut NativeAssembly,
    options: &AnalyticRecoveryOptions,
) -> Result<Vec<CoplanarSupportReassignment>, AnalyticRecoveryError> {
    let trim = options.trim_options(evidence);
    let native_options = NativeTrimOptions {
        max_vertex_displacement: trim.absolute_tolerance,
        intersection_tolerance: (KERNEL_EPSILON * 0.1).min(trim.absolute_tolerance),
        angular_resolution: options.axis_relative_resolution,
    };
    let selected = proposals(evidence, assembly, native_options);
    if selected.is_empty() {
        return Ok(Vec::new());
    }
    let original = &assembly.reconciled.partition;
    let mut labels = original.region_by_support.clone();
    for (&support, &(target, _)) in &selected {
        labels[support] = target;
    }
    let index =
        RegionBoundaryIndex::new(evidence)?.with_source_chord_budget(options.source_chord_budget);
    let partition = index.repartition_connected(&labels, |label| original.carriers[label])?;
    if options.source_chord_violation(&index, &partition).is_some() {
        return Ok(Vec::new());
    }
    let trims = index.assemble_shared_trims(&partition, &assembly.reconciled.positions, trim)?;
    let native = index.construct_native_trims(
        &partition,
        &trims,
        &assembly.reconciled.positions,
        NativeTrimOptions {
            intersection_tolerance: KERNEL_EPSILON * 0.1,
            ..native_options
        },
    )?;
    // Reconstruct once as a transaction. Every previously successful support
    // must retain valid native ownership, and no new failed edge is allowed.
    let candidate = NativeAssembly {
        reconciled: super::ReconciledPartition {
            partition,
            positions: Vec::new(),
            planes: Vec::new(),
            max_displacement: assembly.reconciled.max_displacement,
            rejected: assembly.reconciled.rejected,
        },
        trims,
        native,
    };
    let before = failed_edges(assembly);
    let after = failed_edges(&candidate);
    let old_valid = valid_regions(assembly);
    let new_valid = valid_regions(&candidate);
    let new_partition = &candidate.reconciled.partition;
    if after.len() >= before.len()
        || !after.is_subset(&before)
        || candidate.trims.ambiguous_edges != assembly.trims.ambiguous_edges
        || original
            .region_by_support
            .iter()
            .zip(&new_partition.region_by_support)
            .any(|(old, new)| old_valid.contains(old) && !new_valid.contains(new))
    {
        return Ok(Vec::new());
    }
    let corrections = selected
        .into_iter()
        .map(
            |(support, (_, max_plane_residual))| CoplanarSupportReassignment {
                support,
                old_region: original.region_by_support[support],
                new_region: new_partition.region_by_support[support],
                max_plane_residual,
            },
        )
        .collect();
    // The source and reconciled coordinates and original reconciliation
    // evidence remain intact. Only the connected ownership changes.
    assembly.reconciled.partition = candidate.reconciled.partition;
    for plane in &mut assembly.reconciled.planes {
        plane.region = assembly.reconciled.partition.region_by_support[plane.supports[0]];
    }
    assembly.trims = candidate.trims;
    assembly.native = candidate.native;
    Ok(corrections)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Point3, Vector3, mesh_recovery::analyze_mesh_evidence};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_topology::shell::ShellCondition;

    #[test]
    fn chord_facet_returns_to_its_observed_plane() -> Result<(), Box<dyn std::error::Error>> {
        let n = 16_usize;
        let mut positions = Vec::new();
        for z in [0.0, 2.0] {
            for i in 0..n {
                let theta = std::f64::consts::TAU * i as f64 / n as f64;
                positions.push(Point3::new(theta.cos(), theta.sin(), z));
            }
        }
        // Faces stores triangles, quadrangles and larger polygons in
        // separate buckets. Keep every side polygon triangular so that
        // the intended source facet IDs are independent of that ordering.
        let mut polygons = (0..n - 1)
            .flat_map(|i| [vec![i, i + 1, i + 1 + n], vec![i, i + 1 + n, i + n]])
            .collect::<Vec<_>>();
        let chord = polygons.len();
        // The closing chord is planar, but the cylinder interpolates all
        // its vertices. Simulate an incorrectly claimed triangle.
        polygons.push(vec![n - 1, 0, n]);
        polygons.push(vec![n - 1, n, 2 * n - 1]);
        polygons.push((0..n).rev().collect());
        polygons.push((n..2 * n).collect());
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(polygons),
        )?;
        let options = AnalyticRecoveryOptions::default();
        let original = analyze_mesh_evidence(&mesh, options.evidence)?;
        let evidence = original.refine_planar_supports(&[original.source_planar_support[chord]])?;
        let wrong = evidence.source_planar_support[chord];
        let plane = evidence.source_planar_support[chord + 1];
        assert_ne!(wrong, plane);
        assert!(
            evidence.planar_supports[wrong]
                .normal
                .dot(evidence.planar_supports[plane].normal)
                > 1.0 - 1.0e-12
        );
        let parent = evidence.source_planar_support[0];
        let mut labels = (0..evidence.planar_supports.len()).collect::<Vec<_>>();
        for &support in &evidence.source_planar_support[..=chord] {
            labels[support] = parent;
        }
        let boundary = RegionBoundaryIndex::new(&evidence)?;
        let partition = boundary.repartition_connected(&labels, |label| {
            if label == parent {
                Some(SurfaceCarrier::LinearRevolution {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    axis: Vector3::unit_z(),
                    radius_at_origin: 1.0,
                    slope: 0.0,
                })
            } else {
                let support = &evidence.planar_supports[label];
                Some(SurfaceCarrier::Plane {
                    origin: support.origin,
                    normal: support.normal,
                })
            }
        })?;
        let mut assembly = super::super::assemble_native(&evidence, partition, &options)?;
        assert!(!failed_edges(&assembly).is_empty());
        let original_positions = assembly.reconciled.positions.clone();
        let corrections = repair(&evidence, &mut assembly, &options)?;
        assert_eq!(corrections.len(), 1);
        assert_eq!(corrections[0].support, wrong);
        assert_eq!(assembly.reconciled.positions, original_positions);
        assert_eq!(
            assembly.reconciled.partition.region_by_support[wrong],
            assembly.reconciled.partition.region_by_support[plane]
        );
        assert!(failed_edges(&assembly).is_empty());
        let mut shell = crate::Shell::new();
        for region in &assembly.native.regions {
            assert!(!region.requires_chart_cut);
            shell.extend(region.faces.as_ref().map_err(|e| *e)?.iter().cloned());
        }
        assert_eq!(shell.shell_condition(), ShellCondition::Closed);
        assert!(shell.is_geometric_consistent());
        Ok(())
    }
}
