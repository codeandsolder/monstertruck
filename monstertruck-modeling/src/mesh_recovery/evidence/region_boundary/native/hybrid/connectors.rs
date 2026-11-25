//! Bounded search over explicitly authorized connector-region releases.
//! A new continuous crossing is a reason to propose another transaction,
//! never a reason to admit a self-crossing retained analytic trim.

use super::{HybridRecovery, HybridRecoveryError, NativeTrimError, NativeTrimOptions};
use crate::mesh_recovery::{RecoveredAnalyticMesh, RegionBoundaryIndex, RegionPatchOptions};
use std::collections::{BTreeMap, BTreeSet};

/// Explicit search and geometric budgets for resolving hybrid crossings.
#[derive(Clone, Copy, Debug)]
pub struct HybridConnectorOptions {
    /// Unchanged source approximation bounds.
    pub patch: RegionPatchOptions,
    /// Original strict native construction policy.
    pub strict: NativeTrimOptions,
    /// Maximum additional source polygons authorized for replacement.
    pub max_added_polygons: usize,
    /// Maximum number of complete native transactions.
    pub max_trials: usize,
}

/// Why an extra original region was proposed for replacement.
#[derive(Clone, Copy, Debug)]
pub struct HybridConnectorRelease {
    /// Original region ID.
    pub region: usize,
    /// Crossing retained trim that prompted the release, or None for
    /// a bounded source-fan connector induced by that release.
    pub triggering_trim: Option<usize>,
}

/// Successful connector search and its complete source selection.
#[derive(Debug)]
pub struct HybridConnectorRecovery {
    /// Sorted original selected region IDs.
    pub selected_regions: Vec<usize>,
    /// Additional connector proposals in trial order.
    pub releases: Vec<HybridConnectorRelease>,
    /// Additional source polygons replaced in the successful transaction.
    pub added_polygons: usize,
    /// Number of attempted native transactions.
    pub trials: usize,
    /// Native result with all ordinary hybrid certificates.
    pub recovery: HybridRecovery,
}

type Selection = (usize, Vec<usize>);
type Releases = Vec<HybridConnectorRelease>;

struct ConnectorSearch<'a> {
    recovery: &'a RecoveredAnalyticMesh,
    index: RegionBoundaryIndex<'a>,
    required: BTreeSet<usize>,
    allowed: Vec<usize>,
    sizes: Vec<usize>,
    options: HybridConnectorOptions,
    frontier: BTreeMap<Selection, Releases>,
    visited: BTreeSet<Vec<usize>>,
}

impl<'a> ConnectorSearch<'a> {
    fn new(
        recovery: &'a RecoveredAnalyticMesh,
        required: &[usize],
        allowed: &[usize],
        options: HybridConnectorOptions,
    ) -> Result<Self, HybridRecoveryError> {
        let partition = &recovery.reconciled.partition;
        if required.is_empty()
            || options.max_trials == 0
            || required
                .iter()
                .chain(allowed)
                .any(|&r| partition.carriers.get(r).is_none_or(Option::is_none))
        {
            return Err(HybridRecoveryError::InvalidInput);
        }
        let mut sizes = vec![0; partition.carriers.len()];
        for &id in &partition.region_by_polygon {
            sizes[id] += 1;
        }
        let required = required.iter().copied().collect::<BTreeSet<_>>();
        let initial = required.iter().copied().collect::<Vec<_>>();
        let allowed = allowed
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self {
            recovery,
            index: RegionBoundaryIndex::new(&recovery.evidence)?,
            required,
            allowed,
            sizes,
            options,
            frontier: BTreeMap::from([((0, initial.clone()), Vec::new())]),
            visited: BTreeSet::from([initial]),
        })
    }

    fn successors(
        &mut self,
        selected: &[usize],
        releases: &[HybridConnectorRelease],
        trim: usize,
    ) -> Result<(), HybridRecoveryError> {
        let contact = self
            .recovery
            .trims
            .trims
            .get(trim)
            .ok_or(HybridRecoveryError::InvalidInput)?;
        let partition = &self.recovery.reconciled.partition;
        for region in std::iter::once(contact.regions.0).chain(contact.regions.1) {
            if self.allowed.binary_search(&region).is_err()
                || selected.binary_search(&region).is_ok()
            {
                continue;
            }
            let mut seeds = selected.to_vec();
            seeds.push(region);
            let next = self
                .index
                .complete_patch_fans(partition, &seeds, &self.allowed)?;
            let cost = next
                .iter()
                .filter(|r| !self.required.contains(r))
                .map(|&r| self.sizes[r])
                .sum::<usize>();
            if cost > self.options.max_added_polygons || !self.visited.insert(next.clone()) {
                continue;
            }
            let mut provenance = releases.to_vec();
            provenance.extend(
                next.iter()
                    .filter(|r| selected.binary_search(r).is_err())
                    .map(|&r| HybridConnectorRelease {
                        region: r,
                        triggering_trim: (r == region).then_some(trim),
                    }),
            );
            self.frontier.insert((cost, next), provenance);
        }
        Ok(())
    }

    fn solve(mut self) -> Result<HybridConnectorRecovery, HybridRecoveryError> {
        let mut last = None;
        let mut trials = 0;
        while trials < self.options.max_trials {
            let Some(((added_polygons, selected_regions), releases)) = self.frontier.pop_first()
            else {
                break;
            };
            trials += 1;
            match self.recovery.construct_hybrid(
                &selected_regions,
                self.options.patch,
                self.options.strict,
            ) {
                Ok(recovery) => {
                    return Ok(HybridConnectorRecovery {
                        selected_regions,
                        releases,
                        added_polygons,
                        trials,
                        recovery,
                    });
                }
                Err(error) => {
                    if let HybridRecoveryError::RetainedTrim {
                        trim,
                        error: NativeTrimError::RevolutionTrimCrossing,
                    } = error
                    {
                        self.successors(&selected_regions, &releases, trim)?;
                    }
                    last = Some(error);
                }
            }
        }
        Err(HybridRecoveryError::ConnectorSearch {
            trials,
            last: Box::new(last.unwrap_or(HybridRecoveryError::InvalidInput)),
        })
    }
}

impl RecoveredAnalyticMesh {
    /// Search explicitly permitted extra regions when a hybrid transaction
    /// creates a witnessed continuous crossing of a retained trim.
    ///
    /// Candidate releases come only from the two sides of the rejected
    /// retained contact and the bounded source-fan closure. Try the
    /// smallest additional polygon ownership first; independent alternatives
    /// are retained within the trial budget. This is a proposal search, not
    /// a proof that releasing a carrier is the only possible repair.
    ///
    /// # Errors
    /// Rejects invalid inputs or exhaustion of the authorized search,
    /// retaining the last geometric rejection. No trial may bypass the
    /// continuous crossing, approximation or retained-region checks.
    pub fn construct_hybrid_with_connectors(
        &self,
        required: &[usize],
        permitted_connectors: &[usize],
        options: HybridConnectorOptions,
    ) -> Result<HybridConnectorRecovery, HybridRecoveryError> {
        ConnectorSearch::new(self, required, permitted_connectors, options)?.solve()
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::super::patch_lift::tests::spherical_fixture;
    use super::*;
    use crate::mesh_recovery::{
        AnalyticRecoveryOptions, ReconciledPartition, SharedTrimOptions, SurfaceRegionBudget,
        recover_analytic_evidence,
    };

    #[test]
    fn connector_search_respects_permission_and_polygon_budget()
    -> Result<(), Box<dyn std::error::Error>> {
        let (evidence, partition, fallback, sphere) = spherical_fixture()?;
        let index = RegionBoundaryIndex::new(&evidence)?;
        let trims = index.assemble_shared_trims(
            &partition,
            &evidence.positions,
            SharedTrimOptions {
                absolute_tolerance: 1.0e-7,
                relative_tolerance: 1.0e-5,
                max_interval_edges: 100,
            },
        )?;
        let strict = NativeTrimOptions {
            max_vertex_displacement: 1.0e-7,
            intersection_tolerance: 1.0e-10,
            angular_resolution: 1.0e-8,
        };
        let native =
            index.construct_native_trims(&partition, &trims, &evidence.positions, strict)?;
        let discovery =
            recover_analytic_evidence(&evidence, AnalyticRecoveryOptions::default())?.discovery;
        let reconciled = ReconciledPartition {
            partition,
            positions: evidence.positions.clone(),
            planes: Vec::new(),
            max_displacement: 0.0,
            rejected: 0,
        };
        let recovery = RecoveredAnalyticMesh {
            evidence,
            discovery,
            reconciled,
            trims,
            native,
        };
        let options = HybridConnectorOptions {
            patch: RegionPatchOptions {
                budget: SurfaceRegionBudget {
                    vertex: 1.0e-9,
                    chord: 0.1,
                },
                intersection_tolerance: 1.0e-10,
                angular_resolution: 1.0e-8,
            },
            strict,
            max_added_polygons: 0,
            max_trials: 1,
        };
        let result = recovery.construct_hybrid_with_connectors(&fallback, &[], options)?;
        assert_eq!(result.trials, 1);
        assert_eq!(result.added_polygons, 0);
        assert!(result.releases.is_empty());
        assert_eq!(result.selected_regions, fallback);
        assert!(result.recovery.complete_shell().is_some());

        let required = &fallback[..1];
        let trim = recovery
            .trims
            .trims
            .iter()
            .position(|contact| {
                contact.regions.1.is_some()
                    && contact.regions.0 != required[0]
                    && contact.regions.1 != Some(required[0])
                    && (contact.regions.0 == sphere || contact.regions.1 == Some(sphere))
            })
            .ok_or("retained contact")?;
        let mut denied = ConnectorSearch::new(&recovery, required, &[], options)?;
        denied.frontier.clear();
        denied.successors(required, &[], trim)?;
        assert!(denied.frontier.is_empty());
        let mut over_budget = ConnectorSearch::new(&recovery, required, &[sphere], options)?;
        over_budget.frontier.clear();
        over_budget.successors(required, &[], trim)?;
        assert!(over_budget.frontier.is_empty());
        let mut permitted = ConnectorSearch::new(
            &recovery,
            required,
            &[sphere],
            HybridConnectorOptions {
                max_added_polygons: 4,
                ..options
            },
        )?;
        permitted.frontier.clear();
        permitted.successors(required, &[], trim)?;
        let ((cost, selection), releases) =
            permitted.frontier.pop_first().ok_or("authorized trial")?;
        assert_eq!(cost, 4);
        assert!(selection.contains(&sphere));
        assert_eq!(releases.len(), 1);
        assert_eq!(releases[0].triggering_trim, Some(trim));
        Ok(())
    }
}
