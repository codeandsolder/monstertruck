//! Closed native mesh assembly from certified independent source chains.
//! Intermediate discovery graphs live only for the duration of a component.

use super::{
    AnalyticRecoveryError, AnalyticRecoveryOptions, HybridPlanarSimplification,
    HybridRecoveryError, HybridRefinementOptions, HybridSupportSearchOptions,
    HybridSupportSearchRecovery, NativeTrimOptions, RecoveredAnalyticMesh, RegionPatchOptions,
    SurfaceRegionBudget, partition_recovery_mesh, recover_analytic_mesh,
};
use crate::{KERNEL_EPSILON, Shell};
use monstertruck_mesh::PolygonMesh;
use monstertruck_topology::shell::ShellCondition;
use std::time::Duration;
use web_time::Instant;

/// Explicit authorization for source-backed approximate patches.
#[derive(Clone, Copy, Debug)]
pub struct NativeHybridOptions {
    /// Physical approximation allowances, in source units.
    pub budget: SurfaceRegionBudget,
    /// Total original polygons allowed to relinquish fitted carrier ownership.
    pub max_released_polygons: usize,
    /// Candidate states allowed for each independent recovery job.
    pub max_states_per_component: usize,
    /// Complete native transactions allowed for each independent job.
    pub max_trials_per_component: usize,
    /// Certify and coalesce coplanar patch faces after a successful search.
    pub simplify_planar_fans: bool,
}

impl NativeHybridOptions {
    fn search_options(
        self,
        recovery: &RecoveredAnalyticMesh,
        analytic: &AnalyticRecoveryOptions,
        remaining: usize,
    ) -> HybridSupportSearchOptions {
        let trim = analytic.trim_options(&recovery.evidence);
        let intersection_tolerance = KERNEL_EPSILON * 0.1;
        let angular_resolution = analytic.axis_relative_resolution;
        HybridSupportSearchOptions {
            refinement: HybridRefinementOptions {
                patch: RegionPatchOptions {
                    budget: self.budget,
                    intersection_tolerance,
                    angular_resolution,
                },
                strict: NativeTrimOptions {
                    max_vertex_displacement: trim.absolute_tolerance,
                    intersection_tolerance,
                    angular_resolution,
                },
                trim,
                max_added_polygons: remaining,
            },
            max_states: self.max_states_per_component,
            max_trials: self.max_trials_per_component,
        }
    }
}

/// Policy for complete closed native reconstruction.
#[derive(Clone, Debug, Default)]
pub struct NativeRecoveryOptions {
    /// Analytic surface, topology and trim discovery policy.
    pub analytic: AnalyticRecoveryOptions,
    /// Explicit hybrid permission. None permits only analytic reconstruction.
    pub hybrid: Option<NativeHybridOptions>,
}

impl NativeRecoveryOptions {
    pub(super) fn analytic_policy(&self) -> Result<AnalyticRecoveryOptions, AnalyticRecoveryError> {
        let mut analytic = self.analytic.clone();
        if let Some(hybrid) = self.hybrid {
            let budget = hybrid.budget;
            if !budget.vertex.is_finite()
                || budget.vertex < 0.0
                || !budget.chord.is_finite()
                || budget.chord <= 0.0
                || hybrid.max_states_per_component == 0
                || hybrid.max_trials_per_component == 0
                || analytic
                    .source_chord_budget
                    .is_some_and(|b| b.to_bits() != budget.chord.to_bits())
            {
                return Err(AnalyticRecoveryError::InvalidOptions);
            }
            analytic.source_chord_budget = Some(budget.chord);
        }
        analytic.validate()?;
        Ok(analytic)
    }
}

/// Borrowed evidence delivered before a component's discovery graph is discarded.
///
/// The graph uses component-local indices; `original_polygons` maps to the
/// original assembly input.
pub struct NativeComponentReport<'a> {
    /// Original source polygons in the component mesh's face order.
    pub original_polygons: &'a [usize],
    /// Source chain closure was certified before discovery.
    pub independently_closed_source: bool,
    /// Analytic evidence, or its explicit construction failure.
    pub analytic: Result<&'a RecoveredAnalyticMesh, &'a AnalyticRecoveryError>,
    /// Hybrid result or error, if a hybrid transaction was attempted.
    pub hybrid: Option<Result<&'a HybridSupportSearchRecovery, &'a HybridRecoveryError>>,
    /// Planar simplification certificate, if successful.
    pub simplification: Option<&'a HybridPlanarSimplification>,
    /// Complete closed and nonsingular shell for this component, if recovered.
    pub shell: Option<&'a Shell>,
    /// Analytic discovery and assembly duration.
    pub analytic_duration: Duration,
    /// Hybrid search, patch assembly and simplification duration.
    pub hybrid_duration: Duration,
}

/// Complete assembly result, without retaining every intermediate graph.
#[derive(Debug)]
pub struct NativeMeshRecovery {
    /// Present only if all jobs produced a complete closed, nonsingular,
    /// geometrically consistent assembly. Not a solid/intersection proof.
    pub shell: Option<Shell>,
    /// Number of independently processed source partitions.
    pub component_count: usize,
    /// Total released source polygons, never above the global allowance.
    pub released_polygons: usize,
}

fn closed_shell(shell: &Shell) -> bool {
    !shell.is_empty()
        && shell.shell_condition() == ShellCondition::Closed
        && shell.singular_vertices().is_empty()
        && shell.is_geometric_consistent()
}

struct Attempt {
    analytic: Result<RecoveredAnalyticMesh, AnalyticRecoveryError>,
    hybrid: Option<Result<HybridSupportSearchRecovery, HybridRecoveryError>>,
    simplification: Option<HybridPlanarSimplification>,
    shell: Option<Shell>,
    analytic_duration: Duration,
    hybrid_duration: Duration,
    released_polygons: usize,
}

impl Attempt {
    fn recover(
        mesh: &PolygonMesh,
        analytic: &AnalyticRecoveryOptions,
        hybrid: Option<NativeHybridOptions>,
        remaining: usize,
    ) -> Self {
        let start = Instant::now();
        let recovery = recover_analytic_mesh(mesh, analytic);
        let analytic_duration = start.elapsed();
        let mut attempt = Self {
            analytic: recovery,
            hybrid: None,
            simplification: None,
            shell: None,
            analytic_duration,
            hybrid_duration: Duration::ZERO,
            released_polygons: 0,
        };
        let Ok(recovery) = &attempt.analytic else {
            return attempt;
        };
        attempt.shell = recovery.complete_shell().filter(closed_shell);
        if attempt.shell.is_some() {
            return attempt;
        }
        let Some(hybrid) = hybrid else {
            return attempt;
        };
        let start = Instant::now();
        let permitted = (0..recovery.evidence.planar_supports.len()).collect::<Vec<_>>();
        let mut result = recovery.construct_unresolved_contact_hybrid(
            &permitted,
            hybrid.search_options(recovery, analytic, remaining),
        );
        if let Ok(search) = &mut result {
            if hybrid.simplify_planar_fans {
                attempt.simplification = search.refined.recovery.simplify_planar_fans();
            }
            attempt.released_polygons = search.released_polygons;
            attempt.shell = search
                .refined
                .recovery
                .complete_shell()
                .filter(closed_shell);
        }
        attempt.hybrid = Some(result);
        attempt.hybrid_duration = start.elapsed();
        attempt
    }
}

/// Recover a closed native assembly from independent source chains.
///
/// The observer runs once per component (including failures) while its
/// source and trim evidence is live. Callers needing only a shell can
/// pass |_| {}. Every component shares the declared release allowance;
/// search work limits are per component. No partial assembly is returned.
///
/// # Errors
/// Invalid recovery policy, invalid source ingestion, or inconsistent
/// release accounting.
pub fn recover_native_mesh(
    mesh: &PolygonMesh,
    options: &NativeRecoveryOptions,
    mut observe: impl FnMut(NativeComponentReport<'_>),
) -> Result<NativeMeshRecovery, AnalyticRecoveryError> {
    let analytic = options.analytic_policy()?;
    let partitions = partition_recovery_mesh(mesh)?;
    let component_count = partitions.len();
    let allowance = options.hybrid.map_or(0, |h| h.max_released_polygons);
    let mut remaining = allowance;
    let mut shell = Shell::new();
    let mut all_closed = true;
    for part in partitions {
        let mut attempt = Attempt::recover(&part.mesh, &analytic, options.hybrid, remaining);
        remaining = remaining
            .checked_sub(attempt.released_polygons)
            .ok_or(AnalyticRecoveryError::InvalidOptions)?;
        observe(NativeComponentReport {
            original_polygons: &part.original_polygons,
            independently_closed_source: part.independently_closed,
            analytic: attempt.analytic.as_ref(),
            hybrid: attempt.hybrid.as_ref().map(Result::as_ref),
            simplification: attempt.simplification.as_ref(),
            shell: attempt.shell.as_ref(),
            analytic_duration: attempt.analytic_duration,
            hybrid_duration: attempt.hybrid_duration,
        });
        if let Some(mut recovered) = attempt.shell.take() {
            shell.append(&mut recovered);
        } else {
            all_closed = false;
        }
    }
    Ok(NativeMeshRecovery {
        shell: (all_closed && closed_shell(&shell)).then_some(shell),
        component_count,
        released_polygons: allowance - remaining,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point3;
    use monstertruck_mesh::{Faces, StandardAttributes};

    fn tetrahedron(offset: f64) -> PolygonMesh {
        PolygonMesh::new_unchecked(
            StandardAttributes {
                positions: vec![
                    Point3::new(offset, 0., 0.),
                    Point3::new(offset + 1., 0., 0.),
                    Point3::new(offset, 1., 0.),
                    Point3::new(offset, 0., 1.),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]),
        )
    }

    #[test]
    fn native_assembly_observes_each_closed_chain() -> Result<(), Box<dyn std::error::Error>> {
        let mut mesh = tetrahedron(0.);
        mesh.merge(tetrahedron(1.));
        let mut observed = Vec::new();
        let result = recover_native_mesh(&mesh, &NativeRecoveryOptions::default(), |part| {
            assert!(part.independently_closed_source);
            assert!(part.analytic.is_ok());
            assert!(part.hybrid.is_none());
            assert!(part.shell.is_some());
            observed.extend_from_slice(part.original_polygons);
        })?;
        assert_eq!(observed, (0..8).collect::<Vec<_>>());
        assert_eq!(result.component_count, 2);
        assert_eq!(result.released_polygons, 0);
        assert_eq!(
            result
                .shell
                .ok_or("closed assembly")?
                .connected_components()
                .len(),
            2
        );
        Ok(())
    }

    #[test]
    fn native_assembly_refuses_partial_export() -> Result<(), Box<dyn std::error::Error>> {
        let mut mesh = tetrahedron(0.);
        mesh.merge(PolygonMesh::new_unchecked(
            StandardAttributes {
                positions: vec![
                    Point3::new(4., 0., 0.),
                    Point3::new(5., 0., 0.),
                    Point3::new(4., 1., 0.),
                ],
                ..Default::default()
            },
            Faces::from_iter([[0, 1, 2]]),
        ));
        let mut closed = Vec::new();
        let result = recover_native_mesh(&mesh, &NativeRecoveryOptions::default(), |part| {
            closed.push(part.shell.is_some());
        })?;
        assert_eq!(closed, [true, false]);
        assert!(result.shell.is_none());
        Ok(())
    }

    #[test]
    fn hybrid_policy_cannot_disagree_with_analytic_source_budget() {
        let options = NativeRecoveryOptions {
            analytic: AnalyticRecoveryOptions {
                source_chord_budget: Some(0.01),
                ..Default::default()
            },
            hybrid: Some(NativeHybridOptions {
                budget: SurfaceRegionBudget {
                    vertex: 0.02,
                    chord: 0.02,
                },
                max_released_polygons: 100,
                max_states_per_component: 10,
                max_trials_per_component: 2,
                simplify_planar_fans: true,
            }),
        };
        assert!(matches!(
            options.analytic_policy(),
            Err(AnalyticRecoveryError::InvalidOptions)
        ));
    }
}
