//! A feasible first cover for broad, bounded contact searches.
//!
//! This is an incumbent heuristic, not a geometric admission rule or an
//! assertion of minimum release cost. Exact enumeration remains available.

use super::{CoverCandidate, CoverSearch};
use std::{cmp::Ordering, collections::BTreeSet};

struct GreedySeed {
    support: usize,
    contacts: Vec<usize>,
}

#[derive(Clone, Copy)]
struct Gain {
    support: usize,
    contacts: usize,
    cost: usize,
}

impl Gain {
    fn compare(self, other: Self) -> Ordering {
        // Compare coverage per additional original polygon exactly, without
        // division, float rounding or overflow at usize-width operands.
        ((self.contacts as u128) * (other.cost as u128))
            .cmp(&((other.contacts as u128) * (self.cost as u128)))
            .then_with(|| other.cost.cmp(&self.cost))
            .then_with(|| other.support.cmp(&self.support))
    }
}

impl CoverSearch {
    fn greedy_seeds(&self) -> Vec<GreedySeed> {
        let mut incidence = vec![Vec::new(); self.sizes.len()];
        for (contact, choice) in self.choices.iter().enumerate() {
            for &support in choice {
                incidence[support].push(contact);
            }
        }
        self.candidates
            .iter()
            .map(|&support| GreedySeed {
                support,
                contacts: self.closures[&support]
                    .iter()
                    .flat_map(|&member| incidence[member].iter().copied())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
            })
            .collect()
    }

    fn marginal_gain(
        &self,
        seed: &GreedySeed,
        candidate: &CoverCandidate,
        missing: &[bool],
    ) -> Option<Gain> {
        if candidate.selected.binary_search(&seed.support).is_ok() {
            return None;
        }
        let contacts = seed.contacts.iter().filter(|&&id| missing[id]).count();
        if contacts == 0 {
            return None;
        }
        let cost = self.closures[&seed.support]
            .iter()
            .filter(|id| candidate.selected.binary_search(id).is_err())
            .try_fold(0_usize, |sum, &id| sum.checked_add(self.sizes[id]))?;
        (candidate.cost.checked_add(cost)? <= self.max_polygons).then_some(Gain {
            support: seed.support,
            contacts,
            cost,
        })
    }

    fn prune_greedy_seeds(&self, candidate: CoverCandidate) -> Option<CoverCandidate> {
        let mut counts = vec![0_usize; self.sizes.len()];
        for seed in &candidate.seeds {
            for &id in &self.closures[seed] {
                counts[id] += 1;
            }
        }
        let mut seeds = candidate.seeds;
        seeds.retain(|seed| {
            let closure = &self.closures[seed];
            for &id in closure {
                counts[id] -= 1;
            }
            let needed = self
                .choices
                .iter()
                .any(|choice| choice.iter().all(|&id| counts[id] == 0));
            if needed {
                for &id in closure {
                    counts[id] += 1;
                }
            }
            needed
        });
        let selected = counts
            .iter()
            .enumerate()
            .filter_map(|(id, &count)| (count > 0).then_some(id))
            .collect::<Vec<_>>();
        let cost = selected
            .iter()
            .try_fold(0_usize, |sum, &id| sum.checked_add(self.sizes[id]))?;
        Some(CoverCandidate {
            cost,
            selected,
            seeds,
        })
    }

    fn greedy_cover(&self, mut candidate: CoverCandidate) -> Option<CoverCandidate> {
        if candidate.cost > self.max_polygons {
            return None;
        }
        let seeds = self.greedy_seeds();
        let mut missing = self
            .choices
            .iter()
            .map(|choice| {
                choice
                    .iter()
                    .all(|id| candidate.selected.binary_search(id).is_err())
            })
            .collect::<Vec<_>>();
        while missing.iter().any(|&value| value) {
            let best = seeds
                .iter()
                .filter_map(|seed| self.marginal_gain(seed, &candidate, &missing))
                .max_by(|a, b| a.compare(*b))?;
            candidate = self.expanded(&candidate, &[best.support])?;
            for (id, value) in missing.iter_mut().enumerate() {
                if *value
                    && self.choices[id]
                        .iter()
                        .any(|support| candidate.selected.binary_search(support).is_ok())
                {
                    *value = false;
                }
            }
        }
        self.prune_greedy_seeds(candidate)
    }

    pub(super) fn take_greedy_cover(&mut self) -> Option<CoverCandidate> {
        let root = self.greedy_root.take()?;
        let candidate = self.greedy_cover(root)?;
        if self.visited.len() >= self.max_states || !self.visited.insert(candidate.selected.clone())
        {
            return None;
        }
        self.states += 1;
        Some(candidate)
    }
}
