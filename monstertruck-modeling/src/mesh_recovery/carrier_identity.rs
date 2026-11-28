//! Geometry-aware candidate identity independent of source ownership.

use super::{
    MeshEvidence, StationLinearCandidate, StationMeridianCandidate, SurfaceCarrier,
    SurfaceRegionProposal,
};
use crate::Point3;
use std::collections::BTreeMap;

fn source_ball(evidence: &MeshEvidence, supports: &[usize]) -> Option<(Point3, f64)> {
    let reference = evidence.planar_supports.get(*supports.first()?)?.origin;
    let mut extent = 0.0_f64;
    for &support in supports {
        for &polygon in &evidence.planar_supports.get(support)?.source_face_indices {
            for &vertex in evidence.source_polygon_vertices(polygon)? {
                extent = extent.max((*evidence.positions.get(vertex)? - reference).norm());
            }
        }
    }
    Some((reference, extent))
}

fn equivalent(a: SurfaceCarrier, b: SurfaceCarrier, ball: (Point3, f64), tolerance: f64) -> bool {
    a.compare_in_ball(b, ball.0, ball.1)
        .is_some_and(|fit| fit.residual_difference_upper_bound <= tolerance)
}

pub(in crate::mesh_recovery) trait GeometryCandidate {
    fn supports(&self) -> &[usize];
    fn carrier(&self) -> SurfaceCarrier;
    fn error(&self) -> f64;
}

impl GeometryCandidate for StationLinearCandidate {
    fn supports(&self) -> &[usize] { &self.supports }
    fn carrier(&self) -> SurfaceCarrier { self.carrier }
    fn error(&self) -> f64 { self.max_vertex_residual }
}

impl GeometryCandidate for StationMeridianCandidate {
    fn supports(&self) -> &[usize] { &self.supports }
    fn carrier(&self) -> SurfaceCarrier { SurfaceCarrier::from_circular_meridian(&self.fit) }
    fn error(&self) -> f64 { self.fit.max_vertex_residual }
}

struct GeometryClass<T> {
    anchor: SurfaceCarrier,
    candidate: T,
}

pub(in crate::mesh_recovery) struct GeometryBuckets<T> {
    classes: BTreeMap<Vec<usize>, Vec<GeometryClass<T>>>,
}

impl<T> Default for GeometryBuckets<T> {
    fn default() -> Self {
        Self {
            classes: BTreeMap::new(),
        }
    }
}

impl<T: GeometryCandidate> GeometryBuckets<T> {
    pub(in crate::mesh_recovery) fn insert(
        &mut self,
        candidate: T,
        evidence: &MeshEvidence,
        tolerance: f64,
    ) {
        let classes = self
            .classes
            .entry(candidate.supports().to_vec())
            .or_default();
        let carrier = candidate.carrier();
        if !classes.is_empty()
            && let Some(ball) = source_ball(evidence, candidate.supports())
            && let Some(existing) = classes.iter_mut().find(|class| {
                equivalent(class.anchor, carrier, ball, tolerance)
                    && equivalent(class.candidate.carrier(), carrier, ball, tolerance)
            })
        {
            // Fixed anchors prevent transitive near-coincidence drift.
            if candidate.error() < existing.candidate.error() {
                existing.candidate = candidate;
            }
        } else {
            classes.push(GeometryClass {
                anchor: carrier,
                candidate,
            });
        }
    }

    pub(in crate::mesh_recovery) fn into_candidates(self) -> Vec<T> {
        self.classes
            .into_values()
            .flatten()
            .map(|class| class.candidate)
            .collect()
    }
}

pub(in crate::mesh_recovery) fn append_distinct_proposals<'a>(
    evidence: &MeshEvidence,
    candidates: Vec<SurfaceRegionProposal<'a>>,
    tolerance: f64,
    proposals: &mut Vec<SurfaceRegionProposal<'a>>,
) {
    let mut by_support = BTreeMap::<Vec<usize>, Vec<usize>>::new();
    for (id, p) in proposals.iter().enumerate() {
        by_support
            .entry(p.supports().to_vec())
            .or_default()
            .push(id);
    }
    for candidate in candidates {
        let matches = by_support.entry(candidate.supports().to_vec()).or_default();
        if !matches.is_empty()
            && let Some(ball) = source_ball(evidence, candidate.supports())
            && matches.iter().any(|&id| {
                let previous = &proposals[id];
                (!candidate.independently_tested_profile()
                    || previous.independently_tested_profile())
                    && equivalent(previous.carrier(), candidate.carrier(), ball, tolerance)
            })
        {
            continue;
        }
        matches.push(proposals.len());
        proposals.push(candidate);
    }
}
