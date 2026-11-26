//! Explicitly budgeted support-local recovery of failed analytic contacts.
//! A failed shared edge gives two possible source support releases. Search
//! covers of these incidence constraints by increasing source polygon cost,
//! retaining alternate choices. The native transaction and its continuous
//! source/patch certificates, not the combinatorial cover, decide admission.

use super::{
    HybridRecoveryError, HybridRefinementOptions, HybridRefinementSeed, RefinedHybridRecovery,
};
use crate::mesh_recovery::RecoveredAnalyticMesh;
use monstertruck_topology::shell::ShellCondition;
use std::collections::BTreeSet;

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

struct CoverSearch {
    choices: Vec<Vec<usize>>,
    candidates: Vec<usize>,
    sizes: Vec<usize>,
    max_polygons: usize,
    max_states: usize,
    frontier: BTreeSet<(usize, Vec<usize>)>,
    visited: BTreeSet<Vec<usize>>,
    states: usize,
}

impl CoverSearch {
    fn new(
        contacts: BTreeSet<[usize; 2]>,
        allowed: &BTreeSet<usize>,
        sizes: Vec<usize>,
        options: HybridSupportSearchOptions,
    ) -> Result<Self, HybridRecoveryError> {
        if contacts.is_empty()
            || options.max_states == 0
            || options.max_trials == 0
            || allowed.iter().any(|&id| id >= sizes.len())
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        let mut choices = contacts
            .into_iter()
            .map(|pair| {
                pair.into_iter()
                    .filter(|id| allowed.contains(id))
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
        let cost = mandatory
            .iter()
            .try_fold(0_usize, |total, &id| total.checked_add(search.sizes[id]))
            .ok_or(HybridRecoveryError::InvalidInput)?;
        search.enqueue(cost, mandatory);
        Ok(search)
    }

    fn enqueue(&mut self, cost: usize, selection: Vec<usize>) {
        if cost <= self.max_polygons
            && self.visited.len() < self.max_states
            && self.visited.insert(selection.clone())
        {
            self.frontier.insert((cost, selection));
        }
    }

    fn extend(&mut self, cost: usize, selected: &[usize], support: usize) {
        let Err(at) = selected.binary_search(&support) else {
            return;
        };
        let Some(cost) = cost.checked_add(self.sizes[support]) else {
            return;
        };
        if cost > self.max_polygons || self.visited.len() >= self.max_states {
            return;
        }
        let mut next = selected.to_vec();
        next.insert(at, support);
        self.enqueue(cost, next);
    }

    fn next_cover(&mut self) -> Option<(usize, Vec<usize>)> {
        while let Some((cost, selected)) = self.frontier.pop_first() {
            self.states += 1;
            let missing = self
                .choices
                .iter()
                .find(|choice| !choice.iter().any(|id| selected.binary_search(id).is_ok()))
                .cloned();
            if let Some(missing) = missing {
                for support in missing {
                    self.extend(cost, &selected, support);
                }
            } else {
                return Some((cost, selected));
            }
        }
        None
    }

    // A failed geometric transaction can require both sides of a contact
    // to be released. Consider supersets only of observed failed-contact
    // support IDs, never arbitrary unrelated regions.
    fn extend_rejected(&mut self, cost: usize, selected: &[usize]) {
        for i in 0..self.candidates.len() {
            self.extend(cost, selected, self.candidates[i]);
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

    /// Try source support releases at failed shared native trims.
    ///
    /// Every failed contact must be covered by at least one explicitly
    /// permitted incident support. Enumerate smallest seed polygon costs
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
        let allowed = permitted_supports.iter().copied().collect::<BTreeSet<_>>();
        let sizes = self
            .evidence
            .planar_supports
            .iter()
            .map(|s| s.source_face_indices.len())
            .collect();
        let mut search = CoverSearch::new(contacts, &allowed, sizes, options)?;
        let mut last = HybridRecoveryError::NoSupportCandidate;
        let mut trials = 0;
        while trials < options.max_trials {
            let Some((cost, seed_supports)) = search.next_cover() else {
                break;
            };
            trials += 1;
            let refinement = HybridRefinementOptions {
                max_added_polygons: options.refinement.max_added_polygons - cost,
                ..options.refinement
            };
            match self.construct_refined_hybrid(
                HybridRefinementSeed::Supports(&seed_supports),
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
                            seed_supports,
                            released_polygons: cost + refined.added_polygons,
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
            search.extend_rejected(cost, &seed_supports);
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

    #[test]
    fn weighted_contact_covers_keep_independent_alternatives() -> Result<(), HybridRecoveryError> {
        let contacts = BTreeSet::from([[0, 1], [1, 2]]);
        let mut search = CoverSearch::new(
            contacts,
            &BTreeSet::from([0, 1, 2]),
            vec![2, 3, 2],
            options(7, 32),
        )?;
        assert_eq!(search.next_cover(), Some((3, vec![1])));
        search.extend_rejected(3, &[1]);
        assert_eq!(search.next_cover(), Some((4, vec![0, 2])));
        assert_eq!(search.next_cover(), Some((5, vec![0, 1])));
        assert_eq!(search.next_cover(), Some((5, vec![1, 2])));
        Ok(())
    }

    #[test]
    fn cover_search_obeys_permissions_and_state_budget() -> Result<(), HybridRecoveryError> {
        let contacts = BTreeSet::from([[0, 1], [1, 2]]);
        let mut search = CoverSearch::new(
            contacts.clone(),
            &BTreeSet::from([0, 2]),
            vec![2, 3, 2],
            options(4, 1),
        )?;
        assert_eq!(search.next_cover(), Some((4, vec![0, 2])));
        assert_eq!(search.next_cover(), None);
        assert_eq!(search.states, 1);
        assert!(matches!(
            CoverSearch::new(
                contacts.clone(),
                &BTreeSet::from([0]),
                vec![2, 3, 2],
                options(4, 32)
            ),
            Err(HybridRecoveryError::NoSupportCandidate)
        ));
        let mut search = CoverSearch::new(
            contacts,
            &BTreeSet::from([0, 2]),
            vec![2, 3, 2],
            options(3, 32),
        )?;
        assert_eq!(search.next_cover(), None);
        Ok(())
    }
}
