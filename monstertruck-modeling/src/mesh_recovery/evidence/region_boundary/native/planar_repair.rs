//! Reconsider ambiguous planar circle fits when native boundaries prove a
//! crossing. A circle through a few coplanar samples is not independent
//! evidence that the source polyline is curved.

use super::super::{OrientedTrim, SharedAnalyticTrim, curve_evidence::measure_path_curves};
use super::{
    NativeTrimError, NativeTrimGraph, NativeTrimOptions, RegionBoundaryIndex, SharedTrimGeometry,
    SharedTrimGraph, SurfaceCarrier, SurfaceRegionPartition,
};
use crate::Point3;
use std::collections::BTreeSet;

fn chord_replacement(
    trim: &SharedAnalyticTrim,
    partition: &SurfaceRegionPartition,
    positions: &[Point3],
    tolerance: f64,
) -> Option<Vec<SharedAnalyticTrim>> {
    if !matches!(trim.geometry, SharedTrimGeometry::Circle(_)) {
        return None;
    }
    let carriers = std::iter::once(trim.regions.0)
        .chain(trim.regions.1)
        .map(|id| partition.carriers.get(id).copied().flatten())
        .collect::<Option<Vec<_>>>()?;
    if !carriers
        .iter()
        .all(|c| matches!(c, SurfaceCarrier::Plane { .. }))
    {
        return None;
    }
    trim.source_edges
        .iter()
        .map(|edge| {
            let fit = measure_path_curves(positions, &edge.vertices)?.line?;
            if !carriers.iter().all(|carrier| {
                carrier
                    .measure_segment(fit.endpoints[0], fit.endpoints[1])
                    .is_some_and(|r| r.upper_bound <= tolerance)
            }) {
                return None;
            }
            Some(SharedAnalyticTrim {
                contact_path: trim.contact_path,
                regions: trim.regions,
                source_edges: vec![*edge],
                geometry: SharedTrimGeometry::Line(fit),
            })
        })
        .collect()
}

// Keep a single authoritative trim table. The two incident faces refer
// to the same replacement segments in opposite traversal directions.
fn replace(graph: &mut SharedTrimGraph, id: usize, replacements: Vec<SharedAnalyticTrim>) {
    let count = replacements.len();
    let first_appended = graph.trims.len();
    let mut iter = replacements.into_iter();
    if let Some(first) = iter.next() {
        graph.trims[id] = first;
        graph.trims.extend(iter);
    }
    let ids = std::iter::once(id)
        .chain(first_appended..first_appended + count.saturating_sub(1))
        .collect::<Vec<_>>();
    for region in &mut graph.regions {
        let Ok(boundary) = &mut region.boundary else {
            continue;
        };
        for cycle in &mut boundary.loops {
            if !cycle.iter().any(|r| r.trim == id) {
                continue;
            }
            let mut next = Vec::with_capacity(cycle.len() + count.saturating_sub(1));
            for reference in cycle.iter() {
                if reference.trim != id {
                    next.push(*reference);
                    continue;
                }
                let refs = ids.iter().map(|&trim| OrientedTrim {
                    trim,
                    reversed: reference.reversed,
                });
                if reference.reversed {
                    next.extend(refs.rev());
                } else {
                    next.extend(refs);
                }
            }
            *cycle = next;
        }
    }
}

fn valid_regions(native: &NativeTrimGraph) -> BTreeSet<usize> {
    native
        .regions
        .iter()
        .filter(|r| r.faces.is_ok() && !r.requires_chart_cut)
        .map(|r| r.region)
        .collect()
}

fn failed_edges(graph: &SharedTrimGraph, native: &NativeTrimGraph) -> BTreeSet<[usize; 2]> {
    graph
        .trims
        .iter()
        .zip(&native.trims)
        .filter(|(_, r)| r.is_err())
        .flat_map(|(trim, _)| {
            trim.source_edges.iter().map(|e| {
                let [a, b] = e.vertices;
                [a.min(b), a.max(b)]
            })
        })
        .collect()
}

impl RegionBoundaryIndex<'_> {
    /// Repair a witnessed invalid planar domain by reconsidering ambiguous
    /// free circular trims. Exact source chords are permitted only on
    /// planar carriers; independently determined curved intersections are
    /// never replaced here.
    pub(in crate::mesh_recovery) fn repair_planar_crossings(
        &self,
        partition: &SurfaceRegionPartition,
        positions: &[Point3],
        graph: &mut SharedTrimGraph,
        native: &mut NativeTrimGraph,
        options: NativeTrimOptions,
    ) -> Result<usize, NativeTrimError> {
        let crossings =
            native.planar_crossings(partition, graph, options.intersection_tolerance)?;
        let ids = crossings
            .iter()
            .flat_map(|w| w.edges.map(|(id, _)| id))
            .collect::<BTreeSet<_>>();
        let mut candidate = graph.clone();
        let mut replaced = 0;
        for id in ids {
            let Some(chords) = chord_replacement(
                &graph.trims[id],
                partition,
                positions,
                options.intersection_tolerance,
            ) else {
                continue;
            };
            replace(&mut candidate, id, chords);
            replaced += 1;
        }
        if replaced == 0 {
            return Ok(0);
        }
        let rebuilt = self.construct_native_trims(partition, &candidate, positions, options)?;
        let before = valid_regions(native);
        let after = valid_regions(&rebuilt);
        if after.len() > before.len()
            && before.is_subset(&after)
            && failed_edges(&candidate, &rebuilt).is_subset(&failed_edges(graph, native))
        {
            *graph = candidate;
            *native = rebuilt;
            return Ok(replaced);
        }
        Ok(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{EvidenceOptions, SharedTrimOptions, analyze_mesh_evidence};
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};

    #[test]
    fn recovers_polygonal_boundary_when_free_arc_crosses_opposite_edge()
    -> Result<(), Box<dyn std::error::Error>> {
        // Four samples on a circle, but the longest chord is much coarser
        // than the other two. The source domain is a narrow polygonal strip.
        let circle_point = |a: f64| Point3::new(4.0 * a.sin(), 4.0 * a.cos(), 0.0);
        let mut positions = [-0.58, -0.475, 0.475, 0.58].map(circle_point).to_vec();
        positions.extend([
            Point3::new(positions[3].x, 3.8, 0.0),
            Point3::new(positions[0].x, 3.8, 0.0),
        ]);
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter([vec![0, 1, 2, 3, 4, 5]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let partition = index.select_surface_partition(&[], &evidence, 100)?;
        let mut graph = index.assemble_shared_trims(
            &partition,
            &evidence.positions,
            SharedTrimOptions {
                absolute_tolerance: 1e-8,
                relative_tolerance: 1e-4,
                max_interval_edges: 128,
            },
        )?;
        let options = NativeTrimOptions {
            max_vertex_displacement: 1e-8,
            intersection_tolerance: 1e-10,
            angular_resolution: 1e-8,
        };
        let mut native =
            index.construct_native_trims(&partition, &graph, &evidence.positions, options)?;
        let witnesses =
            native.planar_crossings(&partition, &graph, options.intersection_tolerance)?;
        assert!(
            !witnesses.is_empty(),
            "expected crossing of unsupported free arc"
        );
        let before = graph
            .trims
            .iter()
            .flat_map(|t| t.source_edges.iter().map(|e| e.vertices))
            .collect::<BTreeSet<_>>();
        assert_eq!(
            index.repair_planar_crossings(
                &partition,
                &evidence.positions,
                &mut graph,
                &mut native,
                options
            )?,
            1
        );
        assert!(
            native
                .regions
                .iter()
                .all(|r| r.faces.is_ok() && !r.requires_chart_cut)
        );
        assert!(native.trims.iter().all(Result::is_ok));
        assert_eq!(
            before,
            graph
                .trims
                .iter()
                .flat_map(|t| t.source_edges.iter().map(|e| e.vertices))
                .collect()
        );
        assert_eq!(
            index.repair_planar_crossings(
                &partition,
                &evidence.positions,
                &mut graph,
                &mut native,
                options
            )?,
            0
        );
        Ok(())
    }
}
