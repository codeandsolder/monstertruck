//! Explicitly budgeted support-local recovery of failed analytic contacts.
//! A failed shared edge gives two possible source support releases. Search
//! covers of these incidence constraints by increasing mandatory released polygon cost,
//! retaining alternate choices. The native transaction and its continuous
//! source/patch certificates, not the combinatorial cover, decide admission.

use super::{
    HybridRecoveryError, HybridRefinementOptions, HybridRefinementSeed, HybridSource,
    NativeTrimError, NativeTrimGraph, RefinedHybridRecovery, SurfaceRegionPartition,
    refinement::ReleaseContext,
};
use crate::mesh_recovery::{RecoveredAnalyticMesh, SharedTrimGeometry, SharedTrimGraph};
use monstertruck_topology::shell::ShellCondition;
use std::{
    cmp::Reverse,
    collections::{BTreeMap, BTreeSet},
};

mod greedy;

/// Search limits for local releases around failed source contacts.
#[derive(Clone, Copy, Debug)]
pub struct HybridSupportSearchOptions {
    /// Original geometric policies. Here `max_added_polygons` is the total
    /// authorized source polygons, including initial contact support seeds.
    pub refinement: HybridRefinementOptions,
    /// Maximum distinct support selection states, including queued states.
    pub max_states: usize,
    /// Maximum complete native transactions.
    pub max_trials: usize,
}

/// Successful complete shell from a source-local support release search.
#[derive(Debug)]
pub struct HybridSupportSearchRecovery {
    /// Explicit source supports seeded by failed contact incidence.
    pub seed_supports: Vec<usize>,
    /// Total replaced source polygons, including necessary expansion.
    pub released_polygons: usize,
    /// Number of candidate selection states visited.
    pub states: usize,
    /// Number of complete native transactions attempted.
    pub trials: usize,
    /// Whether the starting partition was a competing ownership fit.
    pub ownership_alternative: bool,
    /// Native result with its refined partition and source certificates.
    pub refined: RefinedHybridRecovery,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CoverCandidate {
    cost: usize,
    selected: Vec<usize>,
    seeds: Vec<usize>,
}

type CrossingChoices = Vec<Vec<Vec<usize>>>;

struct CrossingAlternative {
    seeds: Vec<usize>,
    closure: BTreeSet<usize>,
}

struct CrossingSearch {
    choices: Vec<Vec<CrossingAlternative>>,
    frontier: BTreeSet<(usize, Reverse<usize>, CoverCandidate)>,
}

impl CrossingSearch {
    fn missing(&self, selected: &[usize]) -> Option<&[CrossingAlternative]> {
        self.choices
            .iter()
            .find(|alternatives| {
                !alternatives.iter().any(|alternative| {
                    alternative
                        .closure
                        .iter()
                        .all(|id| selected.binary_search(id).is_ok())
                })
            })
            .map(Vec::as_slice)
    }
}

struct CoverSearch {
    choices: Vec<Vec<usize>>,
    candidates: Vec<usize>,
    sizes: Vec<usize>,
    closures: BTreeMap<usize, BTreeSet<usize>>,
    max_polygons: usize,
    max_states: usize,
    frontier: BTreeSet<(usize, Reverse<usize>, CoverCandidate)>,
    crossing_searches: Vec<CrossingSearch>,
    deferred_rejected: Vec<CoverCandidate>,
    visited: BTreeSet<Vec<usize>>,
    states: usize,
    // Reserve one state for a deterministic feasible cover if exact
    // enumeration is too wide to produce a first native transaction.
    greedy_root: Option<CoverCandidate>,
}

impl CoverSearch {
    fn new(
        contacts: BTreeSet<[usize; 2]>,
        closures: BTreeMap<usize, BTreeSet<usize>>,
        sizes: Vec<usize>,
        options: HybridSupportSearchOptions,
    ) -> Result<Self, HybridRecoveryError> {
        if options.max_states == 0
            || options.max_trials == 0
            || closures.iter().any(|(&id, closure)| {
                !closure.contains(&id) || closure.iter().any(|&member| member >= sizes.len())
            })
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        let mut choices = contacts
            .into_iter()
            .map(|pair| {
                pair.into_iter()
                    .filter(|id| closures.contains_key(id))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        if choices.iter().any(Vec::is_empty) {
            return Err(HybridRecoveryError::NoSupportCandidate);
        }
        choices.sort_by(|a, b| a.len().cmp(&b.len()).then(a.cmp(b)));
        choices.dedup();
        let candidates = choices
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let mut search = Self {
            choices,
            candidates,
            sizes,
            closures,
            max_polygons: options.refinement.max_added_polygons,
            max_states: options.max_states,
            frontier: BTreeSet::new(),
            crossing_searches: Vec::new(),
            deferred_rejected: Vec::new(),
            visited: BTreeSet::new(),
            states: 0,
            greedy_root: None,
        };
        // Every single-choice contact is mandatory. Apply these before
        // branching, avoiding exponential state growth from forced pairs.
        let mandatory = search
            .choices
            .iter()
            .filter(|c| c.len() == 1)
            .map(|c| c[0])
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let selected = mandatory
            .iter()
            .flat_map(|id| search.closures[id].iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let cost = selected
            .iter()
            .try_fold(0_usize, |total, &id| total.checked_add(search.sizes[id]))
            .ok_or(HybridRecoveryError::InvalidInput)?;
        if !search.choices.is_empty() {
            let root = CoverCandidate {
                cost,
                selected,
                seeds: mandatory,
            };
            if search.choices.iter().any(|choice| {
                choice
                    .iter()
                    .all(|id| root.selected.binary_search(id).is_err())
            }) {
                search.greedy_root = Some(root.clone());
            }
            search.enqueue(root);
        }
        Ok(search)
    }

    // Allocate a feasible lower-bound charge for each uncovered contact.
    // Any completion contains one of its full mandatory closures. Shared
    // members are required whichever side wins; charge them once. Then
    // charge the smaller remaining capacity of the disjoint alternatives
    // to both sides. Each support's original cost is allocated at most
    // once across all contacts, so the sum cannot exceed completion cost.
    fn lower_bound(&self, cost: usize, selected: &[usize]) -> Option<usize> {
        let mut capacity = self.sizes.clone();
        for &id in selected {
            capacity[id] = 0;
        }
        let mut bound = cost;
        for choice in &self.choices {
            if choice.iter().any(|id| selected.binary_search(id).is_ok()) {
                continue;
            }
            let charge = match choice.as_slice() {
                [id] => self.closures[id].iter().try_fold(0_usize, |sum, &id| {
                    sum.checked_add(std::mem::take(&mut capacity[id]))
                })?,
                [a, b] => Self::charge_pair(&self.closures[a], &self.closures[b], &mut capacity)?,
                _ => return None,
            };
            bound = bound.checked_add(charge)?;
        }
        Some(bound)
    }

    fn charge_pair(
        a: &BTreeSet<usize>,
        b: &BTreeSet<usize>,
        capacity: &mut [usize],
    ) -> Option<usize> {
        let common = a.intersection(b).try_fold(0_usize, |sum, &id| {
            sum.checked_add(std::mem::take(&mut capacity[id]))
        })?;
        let available = |members: &BTreeSet<usize>, other: &BTreeSet<usize>| {
            members
                .difference(other)
                .try_fold(0_usize, |sum, &id| sum.checked_add(capacity[id]))
        };
        let charge = available(a, b)?.min(available(b, a)?);
        for (members, other) in [(a, b), (b, a)] {
            let mut remaining = charge;
            for &id in members.difference(other) {
                let debit = remaining.min(capacity[id]);
                capacity[id] -= debit;
                remaining -= debit;
                if remaining == 0 {
                    break;
                }
            }
        }
        common.checked_add(charge)
    }

    fn state_limit(&self) -> usize { self.max_states - usize::from(self.greedy_root.is_some()) }

    fn enqueue(&mut self, candidate: CoverCandidate) {
        if self.visited.len() >= self.state_limit() || self.visited.contains(&candidate.selected) {
            return;
        }
        if let Some(bound) = self.lower_bound(candidate.cost, &candidate.selected)
            && bound <= self.max_polygons
        {
            self.visited.insert(candidate.selected.clone());
            // Prefer progress toward a complete cover on equal lower bounds.
            self.frontier
                .insert((bound, Reverse(candidate.cost), candidate));
        }
    }

    fn expanded(&self, candidate: &CoverCandidate, supports: &[usize]) -> Option<CoverCandidate> {
        let mut next = candidate.clone();
        for &support in supports {
            if next.selected.binary_search(&support).is_ok() {
                continue;
            }
            for &member in self.closures.get(&support)? {
                if let Err(at) = next.selected.binary_search(&member) {
                    next.cost = next.cost.checked_add(self.sizes[member])?;
                    if next.cost > self.max_polygons {
                        return None;
                    }
                    next.selected.insert(at, member);
                }
            }
            let at = next.seeds.partition_point(|&id| id < support);
            next.seeds.insert(at, support);
        }
        (next.selected != candidate.selected).then_some(next)
    }

    fn extend(&mut self, candidate: &CoverCandidate, support: usize) {
        if self.visited.len() < self.state_limit()
            && let Some(next) = self.expanded(candidate, &[support])
        {
            self.enqueue(next);
        }
    }

    // A lower-bound charge for a disjunction of support sets. Allocate
    // the smallest available capacity of any alternative, then debit at
    // least that amount from each alternative, counting shared debits
    // only once. Thus any cover pays the charge, and a support's capacity
    // cannot be charged again by a different witness.
    fn charge_alternatives(
        alternatives: &[CrossingAlternative],
        capacity: &mut [usize],
    ) -> Option<usize> {
        let charge = alternatives
            .iter()
            .map(|alternative| {
                alternative
                    .closure
                    .iter()
                    .try_fold(0_usize, |sum, &id| sum.checked_add(capacity[id]))
            })
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .min()?;
        let mut debits = BTreeMap::<usize, usize>::new();
        for alternative in alternatives {
            let paid = alternative.closure.iter().try_fold(0_usize, |sum, id| {
                sum.checked_add(debits.get(id).copied().unwrap_or(0))
            })?;
            let mut remaining = charge.saturating_sub(paid);
            for &id in &alternative.closure {
                let debit = remaining.min(capacity[id]);
                capacity[id] -= debit;
                remaining -= debit;
                *debits.entry(id).or_default() += debit;
                if remaining == 0 {
                    break;
                }
            }
        }
        Some(charge)
    }

    fn crossing_bound(&self, state: &CoverCandidate, search: &CrossingSearch) -> Option<usize> {
        let mut capacity = self.sizes.clone();
        for &id in &state.selected {
            capacity[id] = 0;
        }
        search
            .choices
            .iter()
            .try_fold(state.cost, |bound, alternatives| {
                bound.checked_add(Self::charge_alternatives(alternatives, &mut capacity)?)
            })
    }

    fn expand_crossing(&mut self, state: &CoverCandidate, search: &mut CrossingSearch) {
        let Some(alternatives) = search.missing(&state.selected) else {
            return;
        };
        let next = alternatives
            .iter()
            .filter_map(|alternative| self.expanded(state, &alternative.seeds))
            .collect::<Vec<_>>();
        for candidate in next {
            if self.visited.len() >= self.state_limit() {
                break;
            }
            if self.visited.contains(&candidate.selected) {
                continue;
            }
            if let Some(bound) = self.crossing_bound(&candidate, search)
                && bound <= self.max_polygons
            {
                self.visited.insert(candidate.selected.clone());
                search
                    .frontier
                    .insert((bound, Reverse(candidate.cost), candidate));
            }
        }
    }

    // A crossing can be addressed either by replacing its entire failed
    // face or releasing the supports on the opposite side of one of its
    // curved shared trims. Preserve the alternatives across transactions:
    // rejection of the lowest-cost candidate must not discard the frontier.
    fn witnessed_cover(
        &mut self,
        candidate: &CoverCandidate,
        choices: &[Vec<Vec<usize>>],
    ) -> Option<CoverCandidate> {
        let choices = choices
            .iter()
            .map(|alternatives| {
                alternatives
                    .iter()
                    .filter(|seeds| {
                        seeds
                            .iter()
                            .any(|id| candidate.selected.binary_search(id).is_err())
                    })
                    .map(|seeds| {
                        let mut closure = BTreeSet::new();
                        for id in seeds {
                            closure.extend(self.closures.get(id)?.iter().copied());
                        }
                        Some(CrossingAlternative {
                            seeds: seeds.clone(),
                            closure,
                        })
                    })
                    .collect::<Option<Vec<_>>>()
            })
            .collect::<Option<Vec<_>>>()?;
        let mut search = CrossingSearch {
            choices,
            frontier: BTreeSet::new(),
        };
        self.expand_crossing(candidate, &mut search);
        if !search.frontier.is_empty() {
            self.crossing_searches.push(search);
        }
        self.next_crossing_cover()
    }

    fn next_crossing_cover(&mut self) -> Option<CoverCandidate> {
        while let Some(mut search) = self.crossing_searches.pop() {
            while let Some((_, _, candidate)) = search.frontier.pop_first() {
                self.states += 1;
                if search.missing(&candidate.selected).is_none() {
                    if !search.frontier.is_empty() {
                        self.crossing_searches.push(search);
                    }
                    return Some(candidate);
                }
                self.expand_crossing(&candidate, &mut search);
            }
        }
        None
    }

    fn next_cover(&mut self) -> Option<CoverCandidate> {
        if let Some(candidate) = self.next_crossing_cover() {
            return Some(candidate);
        }
        // Broad supersets are a fallback, not evidence for the next local
        // repair. Do not spend the finite state allowance on them while a
        // witnessed repair chain is still producing transactions.
        for candidate in std::mem::take(&mut self.deferred_rejected) {
            self.extend_rejected(&candidate);
        }
        while let Some((_, _, candidate)) = self.frontier.pop_first() {
            self.states += 1;
            let missing = self
                .choices
                .iter()
                .find(|choice| {
                    !choice
                        .iter()
                        .any(|id| candidate.selected.binary_search(id).is_ok())
                })
                .cloned();
            if let Some(missing) = missing {
                for support in missing {
                    self.extend(&candidate, support);
                }
            } else {
                self.greedy_root = None;
                return Some(candidate);
            }
            // A single complete cover needs at most one new seed per
            // constraint. After that much exact search without a cover,
            // interleave one greedy transaction before consuming the
            // remaining frontier and all its repair capacity.
            if self.states >= self.choices.len().saturating_add(1).min(self.state_limit())
                && let Some(candidate) = self.take_greedy_cover()
            {
                return Some(candidate);
            }
        }
        self.take_greedy_cover()
    }

    // A failed geometric transaction can require both sides of a contact
    // to be released. Consider supersets only of observed failed-contact
    // support IDs, never arbitrary unrelated regions.
    fn extend_rejected(&mut self, candidate: &CoverCandidate) {
        for i in 0..self.candidates.len() {
            self.extend(candidate, self.candidates[i]);
            if self.visited.len() >= self.max_states {
                break;
            }
        }
    }
}

fn incomplete_recovery_error(
    seeds: &[usize],
    refined: &RefinedHybridRecovery,
) -> HybridRecoveryError {
    HybridRecoveryError::IncompleteSupportRecovery {
        seed_supports: seeds.to_vec().into_boxed_slice(),
        failed_regions: refined
            .recovery
            .native
            .regions
            .iter()
            .filter_map(|r| r.faces.as_ref().err().map(|&error| (r.region, error)))
            .collect(),
        failed_trims: refined
            .recovery
            .native
            .trims
            .iter()
            .enumerate()
            .filter_map(|(id, t)| t.as_ref().err().map(|&error| (id, error)))
            .collect(),
        chart_cuts: refined
            .recovery
            .native
            .regions
            .iter()
            .filter(|r| r.requires_chart_cut)
            .map(|r| (r.region, refined.parent_by_region[r.region]))
            .collect(),
        chart_obstacles: refined
            .recovery
            .native
            .regions
            .iter()
            .filter(|r| r.requires_chart_cut && !r.chart_obstacles.is_empty())
            .map(|r| (r.region, r.chart_obstacles.clone()))
            .collect(),
    }
}

impl HybridSource<'_> {
    fn failed_support_contacts(&self) -> Result<BTreeSet<[usize; 2]>, HybridRecoveryError> {
        if !self.trims.ambiguous_edges.is_empty()
            || self.trims.trims.len() != self.native.trims.len()
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        let mut contacts = BTreeSet::new();
        for (trim, native) in self.trims.trims.iter().zip(&self.native.trims) {
            if native.is_ok() {
                continue;
            }
            for edge in &trim.source_edges {
                let a = *self
                    .evidence
                    .source_planar_support
                    .get(edge.interior_polygon)
                    .ok_or(HybridRecoveryError::InvalidInput)?;
                let b = *self
                    .evidence
                    .source_planar_support
                    .get(
                        edge.exterior_polygon
                            .ok_or(HybridRecoveryError::InvalidInput)?,
                    )
                    .ok_or(HybridRecoveryError::InvalidInput)?;
                if a == b {
                    return Err(HybridRecoveryError::InvalidInput);
                }
                contacts.insert([a.min(b), a.max(b)]);
            }
        }
        Ok(contacts)
    }

    fn contact_support_closures(
        &self,
        contacts: &BTreeSet<[usize; 2]>,
        permitted: &[usize],
        context: &ReleaseContext<'_>,
    ) -> Result<BTreeMap<usize, BTreeSet<usize>>, HybridRecoveryError> {
        if permitted
            .iter()
            .any(|&id| id >= self.evidence.planar_supports.len())
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        let allowed = permitted.iter().copied().collect::<BTreeSet<_>>();
        let endpoints = contacts.iter().flatten().copied().collect::<BTreeSet<_>>();
        let mut closures = BTreeMap::new();
        let mut rejected = BTreeMap::new();
        for id in endpoints.intersection(&allowed).copied() {
            match context.closure(id) {
                Ok(closure) => {
                    closures.insert(id, closure);
                }
                Err(
                    error @ (HybridRecoveryError::ForbiddenSupport(_)
                    | HybridRecoveryError::SupportLimit { .. }),
                ) => {
                    rejected.insert(id, error);
                }
                Err(error) => return Err(error),
            }
        }
        // A contact with no feasible endpoint cannot be covered even by
        // the union of other seed closures: mandatory closure is transitive.
        for contact in contacts {
            if contact.iter().all(|id| !closures.contains_key(id)) {
                return Err(contact
                    .iter()
                    .find_map(|id| rejected.remove(id))
                    .unwrap_or(HybridRecoveryError::NoSupportCandidate));
            }
        }
        Ok(closures)
    }

    fn crossing_choices(
        &self,
        partition: &SurfaceRegionPartition,
        graph: &SharedTrimGraph,
        native: &NativeTrimGraph,
        tolerance: f64,
    ) -> Result<CrossingChoices, HybridRecoveryError> {
        let crossings = native.planar_crossings(partition, graph, tolerance)?;
        let mut choices = BTreeSet::new();
        for region in native
            .regions
            .iter()
            .filter(|r| matches!(r.faces, Err(NativeTrimError::PlanarTrimCrossing)))
        {
            let whole = partition
                .region_by_support
                .iter()
                .enumerate()
                .filter_map(|(id, &r)| (r == region.region).then_some(id))
                .collect::<Vec<_>>();
            let mut witnessed = false;
            for crossing in crossings.iter().filter(|c| c.region == region.region) {
                witnessed = true;
                let mut alternatives = vec![whole.clone()];
                for (id, _) in crossing.edges {
                    let trim = graph
                        .trims
                        .get(id)
                        .ok_or(HybridRecoveryError::InvalidInput)?;
                    // Exact straight source edges provide no evidence that
                    // their opposite carrier introduced the crossing.
                    if matches!(
                        trim.geometry,
                        SharedTrimGeometry::Line(_) | SharedTrimGeometry::Unresolved
                    ) {
                        continue;
                    }
                    let other = if trim.regions.0 == region.region {
                        trim.regions.1
                    } else if trim.regions.1 == Some(region.region) {
                        Some(trim.regions.0)
                    } else {
                        return Err(HybridRecoveryError::InvalidInput);
                    };
                    let Some(other) = other else {
                        continue;
                    };
                    let mut supports = BTreeSet::new();
                    for edge in &trim.source_edges {
                        for polygon in
                            std::iter::once(edge.interior_polygon).chain(edge.exterior_polygon)
                        {
                            let support = *self
                                .evidence
                                .source_planar_support
                                .get(polygon)
                                .ok_or(HybridRecoveryError::InvalidInput)?;
                            if partition.region_by_support.get(support) == Some(&other) {
                                supports.insert(support);
                            }
                        }
                    }
                    if !supports.is_empty() {
                        alternatives.push(supports.into_iter().collect());
                    }
                }
                alternatives.sort_unstable();
                alternatives.dedup();
                choices.insert(alternatives);
            }
            if !witnessed {
                choices.insert(vec![whole]);
            }
        }
        Ok(choices.into_iter().collect())
    }

    // A localized chart proof failure names the actual shared native
    // intervals that obstruct admission. Try replacing their incident
    // planar supports before considering an entire dependent rim or region.
    // These are repair alternatives, not assertions of transverse crossing.
    fn chart_edge_supports(
        &self,
        partition: &SurfaceRegionPartition,
        graph: &SharedTrimGraph,
        native: &NativeTrimGraph,
        region: usize,
        (id, ordinal): (usize, usize),
    ) -> [Vec<usize>; 2] {
        let Some(trim) = graph.trims.get(id) else {
            return [Vec::new(), Vec::new()];
        };
        // Angular sections and single-carrier lifts use one native edge
        // per original source segment. If a chart already subdivided the
        // trim, fall back to the complete interval rather than guessing.
        let one_to_one = matches!(
            trim.geometry,
            SharedTrimGeometry::ToroidalSection(_)
                | SharedTrimGeometry::LinearSection(_)
                | SharedTrimGeometry::Conic(_)
                | SharedTrimGeometry::Unresolved
        ) && native
            .trims
            .get(id)
            .and_then(|t| t.as_ref().ok())
            .is_some_and(|edges| edges.len() == trim.source_edges.len());
        let source_edges = if one_to_one {
            trim.source_edges
                .get(ordinal)
                .map_or(trim.source_edges.as_slice(), std::slice::from_ref)
        } else {
            trim.source_edges.as_slice()
        };
        let mut result = [BTreeSet::new(), BTreeSet::new()];
        for edge in source_edges {
            for polygon in std::iter::once(edge.interior_polygon).chain(edge.exterior_polygon) {
                let support = self.evidence.source_planar_support[polygon];
                let owner = partition.region_by_support[support];
                result[usize::from(owner != region)].insert(support);
            }
        }
        result.map(|side| side.into_iter().collect())
    }

    fn chart_choices(
        &self,
        partition: &SurfaceRegionPartition,
        graph: &SharedTrimGraph,
        native: &NativeTrimGraph,
    ) -> CrossingChoices {
        let mut choices = BTreeSet::new();
        for region in native.regions.iter().filter(|r| r.requires_chart_cut) {
            let whole = partition
                .region_by_support
                .iter()
                .enumerate()
                .filter_map(|(id, &r)| (r == region.region).then_some(id))
                .collect::<Vec<_>>();
            if whole.is_empty() {
                continue;
            }
            // Every unresolved chart has a whole-source-region fallback.
            // Prefer a small boundary-side change when the failure exposes
            // a localized span or an affected cycle. These are hypotheses,
            // not assertions that a chart proof failed for geometric reasons.
            let mut fallback = vec![whole];
            if let Some(boundary) = graph
                .regions
                .iter()
                .find(|r| r.region == region.region)
                .and_then(|r| r.boundary.as_ref().ok())
            {
                for cycle in &boundary.loops {
                    let mut adjacent = BTreeSet::new();
                    for reference in cycle {
                        let trim = &graph.trims[reference.trim];
                        // Release the other side of a complete variable-rim
                        // path, independent of the analytic carrier family.
                        if matches!(
                            trim.geometry,
                            SharedTrimGeometry::Line(_) | SharedTrimGeometry::Circle(_)
                        ) {
                            continue;
                        }
                        for ordinal in 0..trim.source_edges.len() {
                            let [_, opposite] = self.chart_edge_supports(
                                partition,
                                graph,
                                native,
                                region.region,
                                (reference.trim, ordinal),
                            );
                            adjacent.extend(opposite);
                        }
                    }
                    if !adjacent.is_empty() {
                        fallback.push(adjacent.into_iter().collect());
                    }
                }
                if fallback.len() > 2 {
                    let combined = fallback[1..]
                        .iter()
                        .flatten()
                        .copied()
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect();
                    fallback.push(combined);
                }
            }
            fallback.sort_unstable();
            fallback.dedup();
            let mut localized = false;
            for obstacle in &region.chart_obstacles {
                let Some(edges) = obstacle.edges() else {
                    continue;
                };
                let mut alternatives = fallback.clone();
                for edge in edges {
                    alternatives.extend(
                        self.chart_edge_supports(partition, graph, native, region.region, edge)
                            .into_iter()
                            .filter(|supports| !supports.is_empty()),
                    );
                }
                alternatives.sort_unstable();
                alternatives.dedup();
                localized = true;
                choices.insert(alternatives);
            }
            if !localized {
                choices.insert(fallback);
            }
        }
        choices.into_iter().collect()
    }

    fn domain_choices(
        &self,
        partition: &SurfaceRegionPartition,
        graph: &SharedTrimGraph,
        native: &NativeTrimGraph,
        tolerance: f64,
    ) -> Result<CrossingChoices, HybridRecoveryError> {
        let mut choices = self.crossing_choices(partition, graph, native, tolerance)?;
        choices.extend(self.chart_choices(partition, graph, native));
        Ok(choices)
    }

    fn ensure_crossing_closures(
        search: &mut CoverSearch,
        supports: &[usize],
        context: &ReleaseContext<'_>,
    ) -> bool {
        if supports.iter().map(|&id| search.sizes[id]).sum::<usize>() > search.max_polygons {
            return false;
        }
        for &support in supports {
            if let std::collections::btree_map::Entry::Vacant(entry) =
                search.closures.entry(support)
            {
                let Ok(closure) = context.closure(support) else {
                    return false;
                };
                entry.insert(closure);
            }
        }
        true
    }

    fn crossing_followup(
        search: &mut CoverSearch,
        candidate: &CoverCandidate,
        choices: &mut [Vec<Vec<usize>>],
        context: &ReleaseContext<'_>,
    ) -> Option<CoverCandidate> {
        for alternatives in choices.iter_mut() {
            alternatives.retain(|set| {
                set.iter()
                    .any(|id| candidate.selected.binary_search(id).is_err())
                    && Self::ensure_crossing_closures(search, set, context)
            });
        }
        if choices.iter().any(Vec::is_empty) {
            return None;
        }
        search.witnessed_cover(candidate, choices)
    }

    fn initial_cover_search(
        &self,
        permitted_supports: &[usize],
        context: &ReleaseContext<'_>,
        options: HybridSupportSearchOptions,
    ) -> Result<(CoverSearch, Option<CoverCandidate>), HybridRecoveryError> {
        let contacts = self.failed_support_contacts()?;
        let closures = self.contact_support_closures(&contacts, permitted_supports, context)?;
        let sizes = self
            .evidence
            .planar_supports
            .iter()
            .map(|s| s.source_face_indices.len())
            .collect();
        let mut search = CoverSearch::new(contacts, closures, sizes, options)?;
        let followup = if search.choices.is_empty() {
            Self::crossing_followup(
                &mut search,
                &CoverCandidate {
                    cost: 0,
                    selected: Vec::new(),
                    seeds: Vec::new(),
                },
                &mut self.domain_choices(
                    &self.reconciled.partition,
                    self.trims,
                    self.native,
                    options.refinement.strict.intersection_tolerance,
                )?,
                context,
            )
        } else {
            None
        };
        Ok((search, followup))
    }

    fn search_contacts(
        &self,
        permitted_supports: &[usize],
        options: HybridSupportSearchOptions,
    ) -> Result<HybridSupportSearchRecovery, HybridRecoveryError> {
        let context = ReleaseContext::new(
            self.evidence,
            &self.reconciled.partition,
            permitted_supports,
            options.refinement,
        )?;
        let (mut search, mut followup) =
            self.initial_cover_search(permitted_supports, &context, options)?;
        let mut last = HybridRecoveryError::NoSupportCandidate;
        let mut trials = 0;
        let mut last_seeds = Vec::new();
        while trials < options.max_trials {
            let Some(candidate) = followup.take().or_else(|| search.next_cover()) else {
                break;
            };
            trials += 1;
            last_seeds.clone_from(&candidate.seeds);
            let seed_cost: usize = candidate.seeds.iter().map(|&id| search.sizes[id]).sum();
            let limit = options.refinement.max_added_polygons - seed_cost;
            match self.construct_refined_with_context(
                HybridRefinementSeed::Supports(&candidate.seeds),
                &context,
                limit,
            ) {
                Ok(refined) => {
                    if refined
                        .recovery
                        .complete_shell()
                        .is_some_and(|shell| shell.shell_condition() == ShellCondition::Closed)
                    {
                        return Ok(HybridSupportSearchRecovery {
                            seed_supports: candidate.seeds,
                            released_polygons: seed_cost + refined.added_polygons,
                            states: search.states,
                            trials,
                            ownership_alternative: false,
                            refined,
                        });
                    }
                    if trials < options.max_trials {
                        match self.domain_choices(
                            &refined.partition,
                            &refined.trims,
                            &refined.recovery.native,
                            options.refinement.strict.intersection_tolerance,
                        ) {
                            Ok(mut choices) => {
                                followup = Self::crossing_followup(
                                    &mut search,
                                    &candidate,
                                    &mut choices,
                                    &context,
                                );
                            }
                            Err(error) => {
                                last = error;
                                break;
                            }
                        }
                    }
                    last = incomplete_recovery_error(&candidate.seeds, &refined);
                }
                Err(error) => {
                    if trials < options.max_trials {
                        let mut choices = context.repair_choices(&error, &candidate.selected);
                        if !choices.is_empty() {
                            followup = Self::crossing_followup(
                                &mut search,
                                &candidate,
                                &mut choices,
                                &context,
                            );
                        }
                    }
                    last = error;
                }
            }
            if trials < options.max_trials {
                search.deferred_rejected.push(candidate);
            }
        }
        Err(HybridRecoveryError::SupportSearch {
            states: search.states,
            queued_states: search.visited.len(),
            trials,
            seed_supports: last_seeds.into_boxed_slice(),
            last: Box::new(last),
        })
    }
}

impl RecoveredAnalyticMesh {
    /// Try source support releases at failed shared native trims and
    /// positively witnessed crossing planar domains.
    ///
    /// Every failed contact must be covered by at least one explicitly
    /// permitted incident support. Enumerate smallest mandatory release costs
    /// first, prioritizing newly witnessed crossing follow-ups over other
    /// contact alternatives, then apply the full hybrid certificates. A successful result is a complete closed shell, not
    /// a material or self-intersection certificate. This bounded search
    /// does not assert global minimum total released polygons.
    ///
    /// # Errors
    /// Invalid inputs, open or ambiguous source contact incidence, no
    /// permitted cover, exhausted budgets, or failed geometric transactions.
    pub fn construct_unresolved_contact_hybrid(
        &self,
        permitted_supports: &[usize],
        options: HybridSupportSearchOptions,
    ) -> Result<HybridSupportSearchRecovery, HybridRecoveryError> {
        let primary = HybridSource::primary(self);
        let Some(alternative) = HybridSource::alternative(self) else {
            return primary.search_contacts(permitted_supports, options);
        };
        // Both candidate searches share the caller's total work budget.
        // Reserve half for the strict starting partition, except when
        // only one trial is authorized and the alternate has fewer failures.
        let alternative_trials = options.max_trials.div_ceil(2);
        let alternative_states = options.max_states.div_ceil(2);
        let alternative_options = HybridSupportSearchOptions {
            max_trials: alternative_trials,
            max_states: alternative_states,
            ..options
        };
        match alternative.search_contacts(permitted_supports, alternative_options) {
            Ok(mut result) => {
                result.ownership_alternative = true;
                Ok(result)
            }
            Err(alternative_error) => {
                let (states, queued_states, trials) = match &alternative_error {
                    HybridRecoveryError::SupportSearch {
                        states,
                        queued_states,
                        trials,
                        ..
                    } => (*states, *queued_states, *trials),
                    _ => (0, 0, 0),
                };
                let remaining = HybridSupportSearchOptions {
                    max_trials: options.max_trials - trials,
                    max_states: options.max_states - queued_states,
                    ..options
                };
                if remaining.max_trials == 0 || remaining.max_states == 0 {
                    return Err(alternative_error);
                }
                match primary.search_contacts(permitted_supports, remaining) {
                    Ok(mut result) => {
                        result.trials += trials;
                        result.states += states;
                        Ok(result)
                    }
                    Err(primary) => Err(HybridRecoveryError::OwnershipAlternatives {
                        alternative: Box::new(alternative_error),
                        primary: Box::new(primary),
                    }),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::{
        NativeTrimOptions, RegionPatchOptions, SharedTrimOptions, SurfaceRegionBudget,
    };

    fn options(limit: usize, states: usize) -> HybridSupportSearchOptions {
        HybridSupportSearchOptions {
            refinement: HybridRefinementOptions {
                patch: RegionPatchOptions {
                    budget: SurfaceRegionBudget {
                        vertex: 1e-5,
                        chord: 0.01,
                    },
                    intersection_tolerance: 1e-7,
                    angular_resolution: 1e-8,
                },
                strict: NativeTrimOptions {
                    max_vertex_displacement: 1e-5,
                    intersection_tolerance: 1e-7,
                    angular_resolution: 1e-8,
                },
                trim: SharedTrimOptions {
                    absolute_tolerance: 1e-5,
                    relative_tolerance: 1e-4,
                    max_interval_edges: 128,
                },
                max_added_polygons: limit,
            },
            max_states: states,
            max_trials: 4,
        }
    }

    fn identity_closures(ids: impl IntoIterator<Item = usize>) -> BTreeMap<usize, BTreeSet<usize>> {
        ids.into_iter()
            .map(|id| (id, BTreeSet::from([id])))
            .collect()
    }

    fn next_pair(search: &mut CoverSearch) -> Option<(usize, Vec<usize>)> {
        search.next_cover().map(|c| (c.cost, c.selected))
    }

    #[test]
    fn witnessed_crossing_followup_preserves_budget_and_other_contact_covers()
    -> Result<(), HybridRecoveryError> {
        let mut search = CoverSearch::new(
            BTreeSet::from([[0, 1]]),
            identity_closures(0..4),
            vec![1, 2, 1, 1],
            options(3, 32),
        )?;
        let first = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        let extra = search
            .witnessed_cover(&first, &[vec![vec![2, 3]]])
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(extra.selected, vec![0, 2, 3]);
        assert_eq!(extra.cost, 3);
        assert!(
            search
                .witnessed_cover(&first, &[vec![vec![2, 3]]])
                .is_none()
        );
        assert_eq!(next_pair(&mut search), Some((2, vec![1])));
        let mut crossing_only = CoverSearch::new(
            BTreeSet::new(),
            identity_closures([0, 1]),
            vec![1, 1],
            options(2, 1),
        )?;
        let empty = CoverCandidate {
            cost: 0,
            selected: Vec::new(),
            seeds: Vec::new(),
        };
        assert_eq!(
            crossing_only
                .witnessed_cover(&empty, &[vec![vec![0, 1]]])
                .map(|c| c.selected),
            Some(vec![0, 1])
        );
        assert!(
            crossing_only
                .witnessed_cover(&empty, &[vec![vec![0]]])
                .is_none()
        );
        assert!(crossing_only.next_cover().is_none());
        Ok(())
    }

    #[test]
    fn consecutive_witness_repairs_do_not_queue_unrelated_supersets()
    -> Result<(), HybridRecoveryError> {
        let mut search = CoverSearch::new(
            (1..=5).map(|id| [0, id]).collect(),
            identity_closures(0..8),
            vec![1; 8],
            options(8, 8),
        )?;
        let first = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(first.selected, vec![0]);
        let second = search
            .witnessed_cover(&first, &[vec![vec![6]]])
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        search.deferred_rejected.push(first);
        assert_eq!(second.selected, vec![0, 6]);
        let third = search
            .witnessed_cover(&second, &[vec![vec![7]]])
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        search.deferred_rejected.push(second);
        assert_eq!(third.selected, vec![0, 6, 7]);
        assert!(search.visited.len() < search.max_states);
        // When the witnessed frontier is exhausted, the ordinary rejected
        // supersets become available again, subject to the same budget.
        assert!(search.next_cover().is_some());
        assert!(search.deferred_rejected.is_empty());
        assert!(search.visited.len() <= search.max_states);
        Ok(())
    }

    #[test]
    fn new_obstacle_cannot_be_satisfied_by_a_previous_release() -> Result<(), HybridRecoveryError> {
        let mut search = CoverSearch::new(
            BTreeSet::from([[0, 1]]),
            identity_closures(0..4),
            vec![1, 2, 1, 1],
            options(3, 32),
        )?;
        let first = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(first.selected, vec![0]);
        let repair = search
            .witnessed_cover(&first, &[vec![vec![0], vec![2, 3]]])
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(repair.selected, vec![0, 2, 3]);
        assert!(search.witnessed_cover(&first, &[vec![vec![0]]]).is_none());
        Ok(())
    }

    #[test]
    fn crossing_cover_charges_unions_of_multi_support_alternatives()
    -> Result<(), HybridRecoveryError> {
        let mut search = CoverSearch::new(
            BTreeSet::new(),
            identity_closures(0..6),
            vec![70, 2, 2, 2, 1, 1],
            options(8, 32),
        )?;
        let empty = CoverCandidate {
            cost: 0,
            selected: Vec::new(),
            seeds: Vec::new(),
        };
        // One large failed planar region can instead release three
        // opposite-side supports. Repeated witnesses share their cost.
        let choices = [
            vec![vec![0], vec![1, 2, 3]],
            vec![vec![0], vec![1, 2, 3]],
            vec![vec![4], vec![5]],
        ];
        let result = search
            .witnessed_cover(&empty, &choices)
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(result.cost, 7);
        assert!(result.selected.starts_with(&[1, 2, 3]));
        assert_eq!(result.selected.len(), 4);
        let other = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(other.cost, 7);
        assert_ne!(other.selected, result.selected);
        assert!(search.visited.len() <= search.max_states);
        Ok(())
    }

    #[test]
    fn disjoint_crossing_choices_use_admissible_union_bound() -> Result<(), HybridRecoveryError> {
        let mut search = CoverSearch::new(
            BTreeSet::new(),
            identity_closures(0..96),
            vec![1; 96],
            options(48, 64),
        )?;
        let empty = CoverCandidate {
            cost: 0,
            selected: Vec::new(),
            seeds: Vec::new(),
        };
        let choices = (0..24)
            .map(|i| vec![vec![4 * i, 4 * i + 1], vec![4 * i + 2, 4 * i + 3]])
            .collect::<Vec<_>>();
        let result = search
            .witnessed_cover(&empty, &choices)
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(result.cost, 48);
        assert_eq!(search.states, 24);
        assert!(search.visited.len() <= 64);
        Ok(())
    }

    #[test]
    fn weighted_contact_covers_keep_independent_alternatives() -> Result<(), HybridRecoveryError> {
        let contacts = BTreeSet::from([[0, 1], [1, 2]]);
        let mut search = CoverSearch::new(
            contacts,
            identity_closures([0, 1, 2]),
            vec![2, 3, 2],
            options(7, 32),
        )?;
        let first = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!((first.cost, first.selected.as_slice()), (3, &[1][..]));
        search.extend_rejected(&first);
        assert_eq!(next_pair(&mut search), Some((4, vec![0, 2])));
        assert_eq!(next_pair(&mut search), Some((5, vec![0, 1])));
        assert_eq!(next_pair(&mut search), Some((5, vec![1, 2])));
        Ok(())
    }

    #[test]
    fn greedy_fallback_finds_a_complete_cover_with_one_state() -> Result<(), HybridRecoveryError> {
        let contacts = (0..24).map(|id| [2 * id, 2 * id + 1]).collect();
        let sizes = (0..48).map(|id| 1 + id % 2).collect::<Vec<_>>();
        let mut search =
            CoverSearch::new(contacts, identity_closures(0..48), sizes, options(24, 1))?;
        let candidate = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!(candidate.cost, 24);
        assert_eq!(candidate.selected, (0..48).step_by(2).collect::<Vec<_>>());
        assert_eq!(search.states, 1);
        assert_eq!(search.visited.len(), 1);
        assert!(search.next_cover().is_none());
        Ok(())
    }

    #[test]
    fn greedy_fallback_charges_shared_closure_members_and_respects_budget()
    -> Result<(), HybridRecoveryError> {
        let contacts = BTreeSet::from([[0, 1], [2, 3], [4, 5]]);
        let mut closures = identity_closures(0..6);
        closures.insert(0, BTreeSet::from([0, 2, 4]));
        let sizes = vec![1, 2, 1, 2, 1, 2];
        for (limit, expected) in [(3, Some((3, vec![0, 2, 4]))), (2, None)] {
            let mut search = CoverSearch::new(
                contacts.clone(),
                closures.clone(),
                sizes.clone(),
                options(limit, 1),
            )?;
            assert_eq!(next_pair(&mut search), expected);
            assert!(search.visited.len() <= 1);
        }
        Ok(())
    }

    #[test]
    fn disjoint_contacts_reach_a_cover_with_linear_state_budget() -> Result<(), HybridRecoveryError>
    {
        let contacts = (0..24).map(|i| [2 * i, 2 * i + 1]).collect();
        let mut search = CoverSearch::new(
            contacts,
            identity_closures(0..48),
            vec![1; 48],
            options(24, 64),
        )?;
        let first = next_pair(&mut search);
        assert_eq!(first, Some((24, (0..48).step_by(2).collect())));
        assert_eq!(search.states, 25);
        assert_eq!(search.next_cover().map(|c| c.cost), Some(24));
        Ok(())
    }

    #[test]
    fn weighted_dual_bounds_and_cover_order_match_exhaustive_search()
    -> Result<(), HybridRecoveryError> {
        let contacts = BTreeSet::from([[0, 1], [1, 2], [0, 2], [2, 3], [3, 4]]);
        let sizes = vec![2, 4, 3, 1, 5];
        let mut search = CoverSearch::new(
            contacts.clone(),
            identity_closures(0..5),
            sizes.clone(),
            options(15, 64),
        )?;
        let selections: Vec<_> = (0..32_u32)
            .map(|mask| {
                let selected: Vec<_> = (0..5).filter(|id| mask & (1 << id) != 0).collect();
                let cost = selected.iter().map(|&id| sizes[id]).sum::<usize>();
                (cost, selected)
            })
            .collect();
        let covers: BTreeSet<_> = selections
            .iter()
            .filter(|(_, s)| {
                contacts
                    .iter()
                    .all(|edge| edge.iter().any(|id| s.contains(id)))
            })
            .cloned()
            .collect();
        for (cost, selected) in &selections {
            let minimum = covers
                .iter()
                .find(|(_, s)| selected.iter().all(|id| s.contains(id)));
            assert!(minimum.is_some_and(|(minimum, _)| {
                search
                    .lower_bound(*cost, selected)
                    .is_some_and(|bound| bound <= *minimum)
            }));
        }
        let mut actual = BTreeSet::new();
        let mut previous = 0;
        while let Some(candidate) = search.next_cover() {
            assert!(candidate.cost >= previous);
            previous = candidate.cost;
            search.extend_rejected(&candidate);
            assert!(actual.insert((candidate.cost, candidate.selected)));
        }
        assert_eq!(actual, covers);
        Ok(())
    }

    #[test]
    fn mandatory_closures_share_cost_and_deduplicate_transactions()
    -> Result<(), HybridRecoveryError> {
        // 0 and 1 belong to the same mandatory release component. 2 and 3
        // each require 4; paying for 4 twice would reject the true optimum.
        let closures = BTreeMap::from([
            (0, BTreeSet::from([0, 1])),
            (1, BTreeSet::from([0, 1])),
            (2, BTreeSet::from([2, 4])),
            (3, BTreeSet::from([3, 4])),
        ]);
        let contacts = BTreeSet::from([[0, 2], [1, 3]]);
        let mut search = CoverSearch::new(contacts, closures, vec![3, 3, 1, 1, 2], options(6, 32))?;
        assert_eq!(search.lower_bound(0, &[]), Some(4));
        let first = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!((first.cost, first.selected.as_slice()), (4, &[2, 3, 4][..]));
        search.extend_rejected(&first);
        let second = search
            .next_cover()
            .ok_or(HybridRecoveryError::NoSupportCandidate)?;
        assert_eq!((second.cost, second.selected.as_slice()), (6, &[0, 1][..]));
        search.extend_rejected(&second);
        assert_eq!(search.next_cover(), None);
        Ok(())
    }

    #[test]
    fn cover_search_obeys_permissions_and_state_budget() -> Result<(), HybridRecoveryError> {
        let contacts = BTreeSet::from([[0, 1], [1, 2]]);
        let mut search = CoverSearch::new(
            contacts.clone(),
            identity_closures([0, 2]),
            vec![2, 3, 2],
            options(4, 1),
        )?;
        assert_eq!(next_pair(&mut search), Some((4, vec![0, 2])));
        assert_eq!(search.next_cover(), None);
        assert_eq!(search.states, 1);
        assert!(matches!(
            CoverSearch::new(
                contacts.clone(),
                identity_closures([0]),
                vec![2, 3, 2],
                options(4, 32)
            ),
            Err(HybridRecoveryError::NoSupportCandidate)
        ));
        let mut search = CoverSearch::new(
            contacts,
            identity_closures([0, 2]),
            vec![2, 3, 2],
            options(3, 32),
        )?;
        assert_eq!(search.next_cover(), None);
        Ok(())
    }
}
