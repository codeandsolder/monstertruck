//! Explicitly budgeted support-local recovery of failed analytic contacts.
//! A failed shared edge gives two possible source support releases. Search
//! covers of these incidence constraints by increasing mandatory released polygon cost,
//! retaining alternate choices. The native transaction and its continuous
//! source/patch certificates, not the combinatorial cover, decide admission.

use super::{
    HybridRecoveryError, HybridRefinementOptions, HybridRefinementSeed, RefinedHybridRecovery,
    refinement::support_closure,
};
use crate::mesh_recovery::RecoveredAnalyticMesh;
use monstertruck_topology::shell::ShellCondition;
use std::{
    cmp::Reverse,
    collections::{BTreeMap, BTreeSet},
};

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
    /// Native result with its refined partition and source certificates.
    pub refined: RefinedHybridRecovery,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CoverCandidate {
    cost: usize,
    selected: Vec<usize>,
    seeds: Vec<usize>,
}

struct CoverSearch {
    choices: Vec<Vec<usize>>,
    candidates: Vec<usize>,
    sizes: Vec<usize>,
    closures: BTreeMap<usize, BTreeSet<usize>>,
    max_polygons: usize,
    max_states: usize,
    frontier: BTreeSet<(usize, Reverse<usize>, CoverCandidate)>,
    visited: BTreeSet<Vec<usize>>,
    states: usize,
}

impl CoverSearch {
    fn new(
        contacts: BTreeSet<[usize; 2]>,
        closures: BTreeMap<usize, BTreeSet<usize>>,
        sizes: Vec<usize>,
        options: HybridSupportSearchOptions,
    ) -> Result<Self, HybridRecoveryError> {
        if contacts.is_empty()
            || options.max_states == 0
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
            visited: BTreeSet::new(),
            states: 0,
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
        search.enqueue(CoverCandidate {
            cost,
            selected,
            seeds: mandatory,
        });
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

    fn enqueue(&mut self, candidate: CoverCandidate) {
        if self.visited.len() >= self.max_states || self.visited.contains(&candidate.selected) {
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

    fn extend(&mut self, candidate: &CoverCandidate, support: usize) {
        if candidate.selected.binary_search(&support).is_ok()
            || self.visited.len() >= self.max_states
        {
            return;
        }
        let mut next = candidate.clone();
        for &member in &self.closures[&support] {
            if let Err(at) = next.selected.binary_search(&member) {
                let Some(cost) = next.cost.checked_add(self.sizes[member]) else {
                    return;
                };
                if cost > self.max_polygons {
                    return;
                }
                next.cost = cost;
                next.selected.insert(at, member);
            }
        }
        let at = next.seeds.partition_point(|&id| id < support);
        next.seeds.insert(at, support);
        self.enqueue(next);
    }

    fn next_cover(&mut self) -> Option<CoverCandidate> {
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
                return Some(candidate);
            }
        }
        None
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

impl RecoveredAnalyticMesh {
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
        options: HybridRefinementOptions,
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
            match support_closure(
                &self.evidence,
                &self.reconciled.partition,
                id,
                permitted,
                options,
            ) {
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

    /// Try source support releases at failed shared native trims.
    ///
    /// Every failed contact must be covered by at least one explicitly
    /// permitted incident support. Enumerate smallest mandatory release costs
    /// first, then use the existing witnessed expansion and full hybrid
    /// certificates. A successful result is a complete closed shell, not
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
        let contacts = self.failed_support_contacts()?;
        let closures =
            self.contact_support_closures(&contacts, permitted_supports, options.refinement)?;
        let sizes = self
            .evidence
            .planar_supports
            .iter()
            .map(|s| s.source_face_indices.len())
            .collect();
        let mut search = CoverSearch::new(contacts, closures, sizes, options)?;
        let mut last = HybridRecoveryError::NoSupportCandidate;
        let mut trials = 0;
        while trials < options.max_trials {
            let Some(candidate) = search.next_cover() else {
                break;
            };
            trials += 1;
            let seed_cost: usize = candidate.seeds.iter().map(|&id| search.sizes[id]).sum();
            let refinement = HybridRefinementOptions {
                max_added_polygons: options.refinement.max_added_polygons - seed_cost,
                ..options.refinement
            };
            match self.construct_refined_hybrid(
                HybridRefinementSeed::Supports(&candidate.seeds),
                permitted_supports,
                refinement,
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
                            refined,
                        });
                    }
                    last = HybridRecoveryError::IncompleteSupportRecovery {
                        failed_regions: refined
                            .recovery
                            .native
                            .regions
                            .iter()
                            .filter(|r| r.faces.is_err())
                            .count(),
                        failed_trims: refined
                            .recovery
                            .native
                            .trims
                            .iter()
                            .filter(|t| t.is_err())
                            .count(),
                        chart_cuts: refined
                            .recovery
                            .native
                            .regions
                            .iter()
                            .filter(|r| r.requires_chart_cut)
                            .count(),
                    };
                }
                Err(error) => last = error,
            }
            search.extend_rejected(&candidate);
        }
        Err(HybridRecoveryError::SupportSearch {
            states: search.states,
            trials,
            last: Box::new(last),
        })
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
