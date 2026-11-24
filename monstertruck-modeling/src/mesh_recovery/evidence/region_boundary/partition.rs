//! Joint source-support selection. A local surface hypothesis never owns
//! another region's polygons merely because its carrier interpolates them.

use super::super::{MeshEvidence, meridian::cluster_meridian, spatial::source_meridian_samples};
use super::{RegionBoundaryError, RegionBoundaryIndex, RegionContactGraph, SurfaceCarrier};
use std::collections::{BTreeMap, VecDeque};

/// Explicit source approximation budgets, not a claim about the native CAD.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceRegionBudget {
    /// Maximum residual at every original source vertex.
    pub vertex: f64,
    /// Maximum continuous source-chord residual.
    pub chord: f64,
}

/// A measured connected region hypothesis. Construct through the source
/// index so support ownership and admissibility cannot diverge.
#[derive(Clone)]
pub struct SurfaceRegionProposal<'a> {
    evidence: &'a MeshEvidence,
    supports: Vec<usize>,
    carrier: SurfaceCarrier,
    source_polygons: usize,
    relative_vertex_error: f64,
    relative_chord_bound: f64,
    independently_tested_profile: bool,
}

impl SurfaceRegionProposal<'_> {
    /// Sorted, unique original support membership.
    #[must_use = "support membership defines candidate conflicts"]
    pub fn supports(&self) -> &[usize] { &self.supports }

    /// Independently fitted surface carrier.
    #[must_use = "the selected region retains its carrier"]
    pub const fn carrier(&self) -> SurfaceCarrier { self.carrier }

    /// Number of original polygons covered.
    #[must_use = "the count measures source coverage"]
    pub const fn source_polygons(&self) -> usize { self.source_polygons }

    /// True if observed meridian samples independently test the profile
    /// beyond its minimum interpolation data. Not a native face certificate.
    #[must_use = "interpolation alone does not resolve profile ambiguity"]
    pub const fn independently_tested_profile(&self) -> bool { self.independently_tested_profile }

    /// Maximum source-vertex residual divided by its budget.
    #[must_use = "fit error is distinct from chord approximation"]
    pub const fn relative_vertex_error(&self) -> f64 { self.relative_vertex_error }

    /// Maximum continuous chord error bound divided by its budget.
    #[must_use = "the chord error bound is used for admissibility"]
    pub const fn relative_chord_bound(&self) -> f64 { self.relative_chord_bound }
}

/// A complete provisional partition of the original support graph.
/// Singleton planes retain unresolved source geometry exactly; they are
/// not necessarily true original planar CAD faces.
#[derive(Debug)]
pub struct SurfaceRegionPartition {
    /// Input proposal indices selected without any overlap.
    pub selected_proposals: Vec<usize>,
    /// Canonical region ID for each original support.
    pub region_by_support: Vec<usize>,
    /// Canonical region ID for every original source polygon.
    pub region_by_polygon: Vec<usize>,
    /// Carrier at each active region ID; inactive IDs have no carrier.
    pub carriers: Vec<Option<SurfaceCarrier>>,
    /// Exact shared original-source contact paths after selection.
    pub contacts: RegionContactGraph,
    /// Number of original planar-support regions eliminated.
    pub consolidated_supports: usize,
    /// Overlap components for which the configured objective was optimized.
    pub exact_components: usize,
    /// Components where the configured search budget limited optimization.
    pub limited_components: usize,
}

impl<'a> RegionBoundaryIndex<'a> {
    /// Validate one connected complete-support union and its carrier.
    /// No native face labels or family-specific preference enter selection.
    #[must_use = "only admissible measured proposals can enter selection"]
    pub fn measure_surface_region(
        &self,
        supports: &[usize],
        carrier: SurfaceCarrier,
        budget: SurfaceRegionBudget,
    ) -> Option<SurfaceRegionProposal<'a>> {
        if supports.len() < 2
            || supports.windows(2).any(|pair| pair[0] >= pair[1])
            || !budget.vertex.is_finite()
            || budget.vertex <= 0.0
            || !budget.chord.is_finite()
            || budget.chord <= 0.0
            || !self.connected_supports(supports)
        {
            return None;
        }
        let boundary = self.boundary(supports).ok()?;
        if !boundary.pinched_vertices.is_empty() {
            return None;
        }
        let mut vertex_error = 0.0_f64;
        for &support in supports {
            for &polygon in &self.evidence.planar_supports[support].source_face_indices {
                for &vertex in self.evidence.source_polygon_vertices(polygon)? {
                    let point = *self.evidence.positions.get(vertex)?;
                    let residual = carrier.measure_segment(point, point)?.witnessed_max;
                    if residual > budget.vertex {
                        return None;
                    }
                    vertex_error = vertex_error.max(residual);
                }
            }
        }
        let chord = carrier.measure_source_chords(self.evidence, supports, budget.chord)?;
        if chord.upper_bound > budget.chord {
            return None;
        }
        Some(SurfaceRegionProposal {
            evidence: self.evidence,
            supports: supports.to_vec(),
            carrier,
            source_polygons: boundary.source_polygon_count,
            relative_vertex_error: vertex_error / budget.vertex,
            relative_chord_bound: chord.upper_bound / budget.chord,
            independently_tested_profile: profile_is_tested(
                self.evidence,
                supports,
                carrier,
                budget.vertex,
            ),
        })
    }

    fn connected_supports(&self, supports: &[usize]) -> bool {
        if supports
            .iter()
            .any(|&id| id >= self.evidence.planar_supports.len())
        {
            return false;
        }
        let mut seen = vec![false; supports.len()];
        seen[0] = true;
        let mut queue = VecDeque::from([supports[0]]);
        while let Some(id) = queue.pop_front() {
            for &interface in &self.evidence.support_interface_indices[id] {
                let [a, b] = self.evidence.support_interfaces[interface].supports;
                let other = if a == id { b } else { a };
                if let Ok(index) = supports.binary_search(&other)
                    && !seen[index]
                {
                    seen[index] = true;
                    queue.push_back(other);
                }
            }
        }
        seen.into_iter().all(|yes| yes)
    }

    /// Select compatible regions maximizing independently tested profile
    /// support coverage first and consolidation second in each
    /// independent overlap component, subject to a search-node budget.
    /// Unselected supports retain their original planar carrier.
    ///
    /// # Errors
    /// Returns an error for inconsistent source evidence or proposals whose
    /// stored membership is not valid for this index.
    pub fn select_surface_partition(
        &self,
        proposals: &[SurfaceRegionProposal<'_>],
        planar_baseline: &MeshEvidence,
        search_nodes_per_component: usize,
    ) -> Result<SurfaceRegionPartition, RegionBoundaryError> {
        let parent_by_support = planar_parents(self.evidence, planar_baseline)
            .ok_or(RegionBoundaryError::InvalidEvidence)?;
        let count = self.evidence.planar_supports.len();
        let mut by_support = vec![Vec::new(); count];
        for (id, proposal) in proposals.iter().enumerate() {
            if !std::ptr::eq(proposal.evidence, self.evidence) {
                return Err(RegionBoundaryError::InvalidEvidence);
            }
            for &support in &proposal.supports {
                by_support
                    .get_mut(support)
                    .ok_or(RegionBoundaryError::InvalidSupport(support))?
                    .push(id);
            }
        }
        let mut order = (0..proposals.len()).collect::<Vec<_>>();
        order.sort_unstable_by(|&a, &b| {
            let (left, right) = (&proposals[a], &proposals[b]);
            right
                .independently_tested_profile
                .cmp(&left.independently_tested_profile)
                .then_with(|| right.supports.len().cmp(&left.supports.len()))
                .then_with(|| {
                    left.relative_vertex_error
                        .total_cmp(&right.relative_vertex_error)
                })
                .then_with(|| {
                    left.relative_chord_bound
                        .total_cmp(&right.relative_chord_bound)
                })
                .then(a.cmp(&b))
        });
        let mut components = overlap_components(proposals, &by_support, &order);
        let mut selected = Vec::new();
        let mut exact_components = 0;
        let mut limited_components = 0;
        for component in &mut components {
            let (chosen, exact) = solve_component(
                component,
                proposals,
                &by_support,
                search_nodes_per_component,
            );
            selected.extend(chosen);
            if exact {
                exact_components += 1;
            } else {
                limited_components += 1;
            }
        }
        selected.sort_unstable();
        let mut region_by_support = vec![usize::MAX; count];
        let mut carriers = vec![None; count];
        for &id in &selected {
            let proposal = &proposals[id];
            let region = proposal.supports[0];
            for &support in &proposal.supports {
                region_by_support[support] = region;
            }
            carriers[region] = Some(proposal.carrier);
        }
        restore_planar_components(
            self.evidence,
            planar_baseline,
            &parent_by_support,
            &mut region_by_support,
            &mut carriers,
        );
        let consolidated_supports = count - carriers.iter().flatten().count();
        let region_by_polygon = self
            .evidence
            .source_planar_support
            .iter()
            .map(|&support| region_by_support[support])
            .collect::<Vec<_>>();
        let contacts = self.partition_contacts(&region_by_polygon)?;
        Ok(SurfaceRegionPartition {
            selected_proposals: selected,
            region_by_support,
            region_by_polygon,
            carriers,
            contacts,
            consolidated_supports,
            exact_components,
            limited_components,
        })
    }
}

// Baseline plane hypotheses remain valid only on their own source polygons.
// Refinement must never alter source geometry or mix independently observed
// planar parents.
fn planar_parents(evidence: &MeshEvidence, baseline: &MeshEvidence) -> Option<Vec<usize>> {
    if evidence.positions != baseline.positions
        || evidence.source_polygon_vertex_offsets != baseline.source_polygon_vertex_offsets
        || evidence.source_polygon_vertex_indices != baseline.source_polygon_vertex_indices
    {
        return None;
    }
    evidence
        .planar_supports
        .iter()
        .map(|support| {
            let parent = *baseline
                .source_planar_support
                .get(*support.source_face_indices.first()?)?;
            baseline.planar_supports.get(parent)?;
            support
                .source_face_indices
                .iter()
                .all(|&polygon| baseline.source_planar_support.get(polygon) == Some(&parent))
                .then_some(parent)
        })
        .collect()
}

fn restore_planar_components(
    evidence: &MeshEvidence,
    baseline: &MeshEvidence,
    parents: &[usize],
    regions: &mut [usize],
    carriers: &mut [Option<SurfaceCarrier>],
) {
    for seed in 0..regions.len() {
        if regions[seed] != usize::MAX {
            continue;
        }
        let parent = parents[seed];
        regions[seed] = seed;
        carriers[seed] = Some(SurfaceCarrier::from_planar_support(
            &baseline.planar_supports[parent],
        ));
        let mut queue = VecDeque::from([seed]);
        while let Some(id) = queue.pop_front() {
            for &interface in &evidence.support_interface_indices[id] {
                let [a, b] = evidence.support_interfaces[interface].supports;
                let other = if a == id { b } else { a };
                if regions[other] == usize::MAX && parents[other] == parent {
                    regions[other] = seed;
                    queue.push_back(other);
                }
            }
        }
    }
}

fn overlap_components(
    proposals: &[SurfaceRegionProposal<'_>],
    by_support: &[Vec<usize>],
    order: &[usize],
) -> Vec<Vec<usize>> {
    let mut visited = vec![false; proposals.len()];
    let mut visited_supports = vec![false; by_support.len()];
    let mut component_by_id = vec![0; proposals.len()];
    let mut components = Vec::<Vec<usize>>::new();
    for &seed in order {
        if visited[seed] {
            continue;
        }
        let component_id = components.len();
        components.push(Vec::new());
        let mut queue = VecDeque::from([seed]);
        visited[seed] = true;
        while let Some(id) = queue.pop_front() {
            component_by_id[id] = component_id;
            for &support in &proposals[id].supports {
                if visited_supports[support] {
                    continue;
                }
                visited_supports[support] = true;
                for &other in &by_support[support] {
                    if !visited[other] {
                        visited[other] = true;
                        queue.push_back(other);
                    }
                }
            }
        }
    }
    for &id in order {
        components[component_by_id[id]].push(id);
    }
    components
}

fn profile_is_tested(
    evidence: &MeshEvidence,
    supports: &[usize],
    carrier: SurfaceCarrier,
    resolution: f64,
) -> bool {
    let (origin, axis, required) = match carrier {
        SurfaceCarrier::Plane { .. } => return true,
        SurfaceCarrier::LinearRevolution { origin, axis, .. } => (origin, axis, 3),
        SurfaceCarrier::CircularRevolution {
            axial_center, axis, ..
        } => (axial_center, axis, 4),
    };
    let Some(samples) = source_meridian_samples(evidence, supports, origin, axis) else {
        return false;
    };
    // Repeated azimuth samples are one independent meridional observation.
    // Singleton samples from an inaccurate placed axis cannot manufacture
    // an independently tested profile.
    cluster_meridian(&samples, resolution)
        .iter()
        .filter(|station| station.members.len() >= 3)
        .count()
        >= required
}

type Score = [usize; 2];

const fn add_score(a: Score, b: Score) -> Score { [a[0] + b[0], a[1] + b[1]] }

fn solve_component(
    component: &[usize],
    proposals: &[SurfaceRegionProposal<'_>],
    by_support: &[Vec<usize>],
    budget: usize,
) -> (Vec<usize>, bool) {
    if component.len() == 1 {
        return (component.to_vec(), true);
    }
    let mut used = std::collections::BTreeSet::new();
    let greedy = component
        .iter()
        .copied()
        .filter(|&id| {
            let supports = &proposals[id].supports;
            if supports.iter().any(|id| used.contains(id)) {
                return false;
            }
            used.extend(supports.iter().copied());
            true
        })
        .collect::<Vec<_>>();
    let ids = component
        .iter()
        .enumerate()
        .map(|(i, &id)| (id, i))
        .collect::<BTreeMap<_, _>>();
    let weights = component
        .iter()
        .map(|&id| {
            let proposal = &proposals[id];
            let count = proposal.supports.len();
            [
                if proposal.independently_tested_profile {
                    count
                } else {
                    0
                },
                count - 1,
            ]
        })
        .collect::<Vec<_>>();
    let conflicts = component
        .iter()
        .map(|&id| {
            proposals[id]
                .supports
                .iter()
                .flat_map(|&support| &by_support[support])
                .map(|other| ids[other])
                .collect::<std::collections::BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let incumbent = greedy.iter().map(|id| ids[id]).collect::<Vec<_>>();
    let mut search = IndependentSet {
        weights: &weights,
        conflicts: &conflicts,
        best_weight: incumbent
            .iter()
            .fold([0, 0], |score, &id| add_score(score, weights[id])),
        best: incumbent,
        remaining_nodes: budget,
        limited: false,
    };
    search.solve();
    (
        search.best.into_iter().map(|id| component[id]).collect(),
        !search.limited,
    )
}

struct SearchState {
    remaining: Vec<usize>,
    chosen: Vec<usize>,
    score: Score,
}

struct IndependentSet<'a> {
    weights: &'a [Score],
    conflicts: &'a [Vec<usize>],
    best: Vec<usize>,
    best_weight: Score,
    remaining_nodes: usize,
    limited: bool,
}

impl IndependentSet<'_> {
    fn solve(&mut self) {
        let mut stack = vec![SearchState {
            remaining: (0..self.weights.len()).collect(),
            chosen: Vec::new(),
            score: [0, 0],
        }];
        while let Some(state) = stack.pop() {
            if state.score > self.best_weight {
                self.best_weight = state.score;
                self.best.clone_from(&state.chosen);
            }
            if state.remaining.is_empty() {
                continue;
            }
            let bound = state
                .remaining
                .iter()
                .fold(state.score, |score, &id| add_score(score, self.weights[id]));
            if bound <= self.best_weight {
                continue;
            }
            if self.remaining_nodes == 0 {
                self.limited = true;
                break;
            }
            self.remaining_nodes -= 1;
            // Independent choices can be forced into the optimum. For
            // competing choices branch on the strongest source evidence
            // first, then on the densest overlap.
            let Some(&id) = state.remaining.iter().max_by_key(|&&id| {
                (
                    self.weights[id],
                    self.conflicts[id].len(),
                    std::cmp::Reverse(id),
                )
            }) else {
                continue;
            };
            let include = state
                .remaining
                .iter()
                .copied()
                .filter(|other| self.conflicts[id].binary_search(other).is_err())
                .collect::<Vec<_>>();
            if include.len() + 1 != state.remaining.len() {
                stack.push(SearchState {
                    remaining: state
                        .remaining
                        .into_iter()
                        .filter(|&other| other != id)
                        .collect(),
                    chosen: state.chosen.clone(),
                    score: state.score,
                });
            }
            let mut chosen = state.chosen;
            chosen.push(id);
            stack.push(SearchState {
                remaining: include,
                chosen,
                score: add_score(state.score, self.weights[id]),
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn global_selection_beats_largest_first() {
        let weights = [[0, 4], [0, 3], [0, 3]];
        let conflicts = [vec![0, 1, 2], vec![0, 1], vec![0, 2]];
        let mut search = IndependentSet {
            weights: &weights,
            conflicts: &conflicts,
            best: vec![0],
            best_weight: [0, 4],
            remaining_nodes: 100,
            limited: false,
        };
        search.solve();
        search.best.sort_unstable();
        assert_eq!(search.best, vec![1, 2]);
        assert_eq!(search.best_weight, [0, 6]);
        assert!(!search.limited);
    }

    #[test]
    fn independent_profile_evidence_outweighs_interpolation_coverage() {
        let weights = [[3, 2], [0, 50], [0, 50]];
        let conflicts = [vec![0, 1, 2], vec![0, 1], vec![0, 2]];
        let mut search = IndependentSet {
            weights: &weights,
            conflicts: &conflicts,
            best: vec![1, 2],
            best_weight: [0, 100],
            remaining_nodes: 100,
            limited: false,
        };
        search.solve();
        assert_eq!(search.best, vec![0]);
        assert_eq!(search.best_weight, [3, 2]);
        assert!(!search.limited);
    }

    #[test]
    fn search_limit_preserves_feasible_incumbent() {
        let weights = [[0, 4], [0, 3], [0, 3]];
        let conflicts = [vec![0, 1, 2], vec![0, 1], vec![0, 2]];
        let mut search = IndependentSet {
            weights: &weights,
            conflicts: &conflicts,
            best: vec![0],
            best_weight: [0, 4],
            remaining_nodes: 0,
            limited: false,
        };
        search.solve();
        assert_eq!(search.best, vec![0]);
        assert!(search.limited);
    }

    #[test]
    fn large_overlap_component_does_not_have_a_machine_word_limit() {
        let weights = vec![[0, 1]; 80];
        let conflicts = (0..80)
            .map(|_| (0..80).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let mut search = IndependentSet {
            weights: &weights,
            conflicts: &conflicts,
            best: vec![0],
            best_weight: [0, 1],
            remaining_nodes: 1000,
            limited: false,
        };
        search.solve();
        assert_eq!(search.best_weight, [0, 1]);
        assert!(!search.limited);
    }

    #[test]
    fn empty_proposals_preserve_source_partition() -> Result<(), Box<dyn std::error::Error>> {
        use crate::Point3;
        use crate::mesh_recovery::{EvidenceOptions, analyze_mesh_evidence};
        use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
        let mesh = PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(0.0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2], [1, 0, 3]]),
        )?;
        let evidence = analyze_mesh_evidence(&mesh, EvidenceOptions::default())?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let partition = index.select_surface_partition(&[], &evidence, 100)?;
        assert_eq!(partition.region_by_polygon, evidence.source_planar_support);
        assert_eq!(partition.consolidated_supports, 0);
        assert_eq!(
            partition.contacts,
            index.partition_contacts(&evidence.source_planar_support)?
        );
        assert!(partition.selected_proposals.is_empty());
        Ok(())
    }
}
