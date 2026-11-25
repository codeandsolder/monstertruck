//! Label-blind feasibility census of bounded hybrid region replacements.
//! Strict recovery remains unchanged; the global trial can export a complete hybrid shell.
//! Each attempted patch retains all other regions as analytic constraints.

use anyhow::{Context, Result};
use monstertruck_mesh::obj;
use monstertruck_modeling::mesh_recovery::{
    AnalyticRecoveryOptions, HybridConnectorOptions, HybridRecoveryError, HybridRefinementOptions,
    NativeTrimOptions, RecoveredAnalyticMesh, RegionBoundaryIndex, RegionPatchError,
    RegionPatchOptions, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionPartition,
    recover_analytic_mesh,
};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

fn connected_groups(
    partition: &SurfaceRegionPartition,
    candidates: BTreeSet<usize>,
) -> Result<Vec<Vec<usize>>> {
    let mut touches_curved = BTreeSet::new();
    for path in &partition.contacts.paths {
        let (a, Some(b)) = path.regions else {
            continue;
        };
        for (inside, outside) in [(a, b), (b, a)] {
            if candidates.contains(&inside)
                && !candidates.contains(&outside)
                && !matches!(
                    partition.carriers[outside],
                    Some(SurfaceCarrier::Plane { .. })
                )
            {
                touches_curved.insert(inside);
            }
        }
    }
    Ok(partition
        .connected_region_groups(&candidates.into_iter().collect::<Vec<_>>())?
        .into_iter()
        .filter(|group| group.iter().any(|id| touches_curved.contains(id)))
        .collect())
}

struct CandidateSelection {
    groups: Vec<Vec<usize>>,
    initial_groups: Vec<Vec<usize>>,
    required: BTreeSet<usize>,
    polygons_by_region: BTreeMap<usize, usize>,
}

fn candidate_groups(
    recovery: &RecoveredAnalyticMesh,
    index: &RegionBoundaryIndex<'_>,
) -> Result<CandidateSelection> {
    let partition = &recovery.reconciled.partition;
    let failed = recovery
        .native
        .regions
        .iter()
        .filter(|r| {
            r.faces.is_err()
                && matches!(
                    partition.carriers[r.region],
                    Some(SurfaceCarrier::Plane { .. })
                )
        })
        .map(|r| r.region)
        .collect::<BTreeSet<_>>();
    let initial = connected_groups(partition, failed)?;
    let required = initial.iter().flatten().copied().collect::<BTreeSet<_>>();
    let mut polygons_by_region = BTreeMap::<usize, usize>::new();
    for &region in &partition.region_by_polygon {
        *polygons_by_region.entry(region).or_default() += 1;
    }
    // One source facet is the minimal interpolation data for a plane.
    // Do not automatically replace any region backed by multiple facets.
    let connectors = polygons_by_region
        .iter()
        .filter_map(|(&region, &count)| {
            (count == 1
                && matches!(
                    partition.carriers[region],
                    Some(SurfaceCarrier::Plane { .. })
                ))
            .then_some(region)
        })
        .collect::<Vec<_>>();
    let selected = index.complete_patch_fans(
        partition,
        &required.iter().copied().collect::<Vec<_>>(),
        &connectors,
    )?;
    let groups = connected_groups(partition, selected.into_iter().collect())?;
    Ok(CandidateSelection {
        groups,
        initial_groups: initial,
        required,
        polygons_by_region,
    })
}

const fn failure_class(error: RegionPatchError) -> &'static str {
    match error {
        RegionPatchError::InvalidRegion => "region",
        RegionPatchError::Junction(_) => "junction",
        RegionPatchError::Edge(_) => "edge",
        RegionPatchError::Assembly(_) => "assembly",
    }
}

#[derive(Default)]
struct Census {
    groups: Vec<Value>,
    failures: BTreeMap<&'static str, usize>,
    selected: BTreeSet<usize>,
    attempted_polygons: usize,
    connector_polygons: usize,
    accepted_polygons: usize,
    accepted_required_polygons: usize,
    accepted_groups: usize,
}

struct PatchAttempt<'a> {
    index: &'a RegionBoundaryIndex<'a>,
    partition: &'a SurfaceRegionPartition,
    selection: &'a CandidateSelection,
    options: RegionPatchOptions,
}
impl PatchAttempt<'_> {
    fn attempt(&self, group: &[usize], census: &mut Census) -> bool {
        let polygons = group
            .iter()
            .map(|id| self.selection.polygons_by_region[id])
            .sum::<usize>();
        let required = group
            .iter()
            .filter(|id| self.selection.required.contains(id))
            .map(|id| self.selection.polygons_by_region[id])
            .sum::<usize>();
        census.attempted_polygons += polygons;
        census.connector_polygons += polygons - required;
        let result = match self
            .index
            .approximate_region_patch(self.partition, group, self.options)
        {
            Ok(patch) => {
                census.accepted_groups += 1;
                census.accepted_polygons += patch.source_polygons.len();
                census.accepted_required_polygons += required;
                census.selected.extend(group.iter().copied());
                json!({"regions":group,"source_polygons":polygons,"required_polygons":required,
                    "native_faces":patch.approximation.shell.len(),"boundary_edges":patch.boundary_edges.len(),
                    "vertex_displacement":patch.approximation.vertex_displacement,
                    "deviation_bound":patch.approximation.deviation_upper_bound,
                    "seam_deviation_bound":patch.approximation.seam_deviation_upper_bound})
            }
            Err(error) => {
                *census.failures.entry(failure_class(error)).or_default() += 1;
                json!({"regions":group,"source_polygons":polygons,"required_polygons":required,"error":error.to_string()})
            }
        };
        let accepted = result.get("error").is_none();
        census.groups.push(result);
        accepted
    }
}
impl Census {
    fn report(&self) -> Value {
        json!({"attempted_groups":self.groups.len(),"attempted_polygons":self.attempted_polygons,
            "connector_polygons":self.connector_polygons,"accepted_groups":self.accepted_groups,
            "accepted_polygons":self.accepted_polygons,"accepted_required_polygons":self.accepted_required_polygons,
            "failures":self.failures,"groups":self.groups})
    }
}

fn retained_trim_diagnostic(
    recovery: &RecoveredAnalyticMesh,
    error: &HybridRecoveryError,
) -> Value {
    if let HybridRecoveryError::ConnectorSearch { last, .. } = error {
        return retained_trim_diagnostic(recovery, last);
    }
    let HybridRecoveryError::RetainedTrim { trim: id, .. } = error else {
        return Value::Null;
    };
    let Some(trim) = recovery.trims.trims.get(*id) else {
        return Value::Null;
    };
    let endpoints = trim
        .source_edges
        .iter()
        .flat_map(|edge| edge.vertices)
        .collect::<BTreeSet<_>>();
    let partition = &recovery.reconciled.partition;
    json!({"id":id,"regions":[Some(trim.regions.0),trim.regions.1],
      "geometry":format!("{:?}",trim.geometry),
      "carriers":std::iter::once(trim.regions.0).chain(trim.regions.1).map(|r|format!("{:?}",partition.carriers[r])).collect::<Vec<_>>(),
      "edges":trim.source_edges.iter().map(|e|json!({"vertices":e.vertices,"polygons":[Some(e.interior_polygon),e.exterior_polygon]})).collect::<Vec<_>>(),
      "vertices":endpoints.into_iter().map(|id|json!({"id":id,"source":recovery.evidence.positions[id],"reconciled":recovery.reconciled.positions[id],"native":recovery.native.vertices[id].as_ref().map(monstertruck_modeling::Vertex::point)})).collect::<Vec<_>>()})
}

fn global_trial(
    recovery: &RecoveredAnalyticMesh,
    options: RegionPatchOptions,
    selected: &BTreeSet<usize>,
    shell_output: Option<&Path>,
    crossing_limit: Option<usize>,
) -> Result<Value> {
    if selected.is_empty() {
        return Ok(json!({"status":"no_candidates"}));
    }
    let start = Instant::now();
    let analytic = AnalyticRecoveryOptions::default();
    let strict = NativeTrimOptions {
        max_vertex_displacement: analytic.trim_options(&recovery.evidence).absolute_tolerance,
        intersection_tolerance: monstertruck_modeling::KERNEL_EPSILON * 0.1,
        angular_resolution: analytic.axis_relative_resolution,
    };
    let required = selected.iter().copied().collect::<Vec<_>>();
    let attempt = if let Some(limit) = crossing_limit {
        let partition = &recovery.reconciled.partition;
        let mut counts = BTreeMap::<usize, usize>::new();
        for &region in &partition.region_by_polygon {
            *counts.entry(region).or_default() += 1;
        }
        let allowed = counts
            .iter()
            .filter_map(|(&region, &count)| {
                (count == 1
                    && matches!(
                        partition.carriers[region],
                        Some(SurfaceCarrier::Plane { .. })
                    ))
                .then_some(region)
            })
            .collect::<Vec<_>>();
        recovery
            .construct_hybrid_with_connectors(
                &required,
                &allowed,
                HybridConnectorOptions {
                    patch: options,
                    strict,
                    max_added_polygons: limit,
                    max_trials: limit.saturating_mul(4).max(1),
                },
            )
            .map(|result| {
                (
                    result.recovery,
                    result.selected_regions.len(),
                    result.added_polygons,
                    result.trials,
                    result
                        .releases
                        .iter()
                        .map(|r| json!({"region":r.region,"triggering_trim":r.triggering_trim}))
                        .collect::<Vec<_>>(),
                )
            })
    } else {
        recovery
            .construct_hybrid(&required, options, strict)
            .map(|result| (result, required.len(), 0, 1, Vec::new()))
    };
    Ok(match attempt {
        Ok((result, selected_count, added_polygons, trials, releases)) => {
            let before = recovery
                .native
                .regions
                .iter()
                .filter(|r| r.faces.is_err())
                .count();
            let after = result
                .native
                .regions
                .iter()
                .filter(|r| r.faces.is_err())
                .count();
            let old_trims = recovery.native.trims.iter().filter(|r| r.is_err()).count();
            let new_trims = result.native.trims.iter().filter(|r| r.is_err()).count();
            let shell = result.complete_shell();
            let shell_exported = if let (Some(path), Some(shell)) = (shell_output, &shell) {
                let mut writer = BufWriter::new(
                    File::create_new(path)
                        .with_context(|| format!("creating {}", path.display()))?,
                );
                serde_json::to_writer(&mut writer, &shell.compress())?;
                writer.flush()?;
                true
            } else {
                false
            };
            json!({"status":"ok","selected_regions":selected_count,"connector_added_polygons":added_polygons,"connector_trials":trials,"connector_releases":releases,
                "patches":result.patches.len(),"lifted_trims":result.lifted_trims.len(),
                "before_failed_regions":before,"failed_regions":after,
                "before_failed_trims":old_trims,"failed_trims":new_trims,
                "chart_cuts":result.native.regions.iter().filter(|r|r.requires_chart_cut).count(),
                "approximate_polygons":result.patches.iter().map(|p|p.polygons.len()).sum::<usize>(),
                "max_deviation":result.patches.iter().map(|p|p.deviation_upper_bound).fold(0.0_f64,f64::max),
                "max_seam_deviation":result.patches.iter().map(|p|p.seam_deviation_upper_bound).fold(0.0_f64,f64::max),
                "shell_condition":shell.as_ref().map(|s|format!("{:?}",s.shell_condition())),
                "native_faces":shell.as_ref().map(|s|s.len()),"shell_exported":shell_exported,
                "seconds":start.elapsed().as_secs_f64()})
        }
        Err(error) => json!({"status":"rejected","selected_regions":selected.len(),
            "error":error.to_string(),"retained_trim":retained_trim_diagnostic(recovery,&error),"seconds":start.elapsed().as_secs_f64()}),
    })
}

fn census(
    recovery: &RecoveredAnalyticMesh,
    budget: SurfaceRegionBudget,
    shell_output: Option<&Path>,
    crossing_limit: Option<usize>,
) -> Result<Value> {
    let index = RegionBoundaryIndex::new(&recovery.evidence)?;
    let partition = &recovery.reconciled.partition;
    let options = RegionPatchOptions {
        budget,
        intersection_tolerance: monstertruck_modeling::KERNEL_EPSILON * 0.1,
        angular_resolution: 1.0e-8,
    };
    let start = Instant::now();
    let selection = candidate_groups(recovery, &index)?;
    let attempt = PatchAttempt {
        index: &index,
        partition,
        selection: &selection,
        options,
    };
    let mut expanded = Census::default();
    let mut fallback = Census::default();
    let original_by_region = selection
        .initial_groups
        .iter()
        .enumerate()
        .flat_map(|(id, group)| group.iter().map(move |&r| (r, id)))
        .collect::<BTreeMap<_, _>>();
    for group in &selection.groups {
        if attempt.attempt(group, &mut expanded) {
            continue;
        }
        let original = group
            .iter()
            .filter_map(|r| original_by_region.get(r))
            .copied()
            .collect::<BTreeSet<_>>();
        for id in original {
            let subset = &selection.initial_groups[id];
            if subset != group {
                attempt.attempt(subset, &mut fallback);
            }
        }
    }
    let seconds = start.elapsed().as_secs_f64();
    let selected = expanded
        .selected
        .union(&fallback.selected)
        .copied()
        .collect::<BTreeSet<_>>();
    let mut report = expanded.report();
    let fields = report.as_object_mut().context("census object")?;
    fields.extend([
        ("vertex_budget".to_string(), json!(budget.vertex)),
        ("chord_budget".to_string(), json!(budget.chord)),
        ("seconds".to_string(), json!(seconds)),
        ("fallback".to_string(), fallback.report()),
        (
            "combined_accepted_polygons".to_string(),
            json!(expanded.accepted_polygons + fallback.accepted_polygons),
        ),
        (
            "combined_accepted_required_polygons".to_string(),
            json!(expanded.accepted_required_polygons + fallback.accepted_required_polygons),
        ),
        (
            "global_hybrid".to_string(),
            global_trial(recovery, options, &selected, shell_output, crossing_limit)?,
        ),
    ]);
    Ok(report)
}

fn refined_trial(
    recovery: &RecoveredAnalyticMesh,
    budget: SurfaceRegionBudget,
    limit: usize,
    shell_output: Option<&Path>,
) -> Result<Value> {
    let index = RegionBoundaryIndex::new(&recovery.evidence)?;
    let selection = candidate_groups(recovery, &index)?;
    let required = selection
        .groups
        .iter()
        .flatten()
        .copied()
        .collect::<BTreeSet<_>>();
    let partition = &recovery.reconciled.partition;
    // This opt-in diagnostic permits only support-local releases from
    // originally fitted curved regions, not removal of established planes.
    let permitted = partition
        .region_by_support
        .iter()
        .enumerate()
        .filter_map(|(support, &region)| {
            (!matches!(
                partition.carriers[region],
                Some(SurfaceCarrier::Plane { .. })
            ))
            .then_some(support)
        })
        .collect::<Vec<_>>();
    let analytic = AnalyticRecoveryOptions::default();
    let trim = analytic.trim_options(&recovery.evidence);
    let options = HybridRefinementOptions {
        patch: RegionPatchOptions {
            budget,
            intersection_tolerance: monstertruck_modeling::KERNEL_EPSILON * 0.1,
            angular_resolution: analytic.axis_relative_resolution,
        },
        strict: NativeTrimOptions {
            max_vertex_displacement: trim.absolute_tolerance,
            intersection_tolerance: monstertruck_modeling::KERNEL_EPSILON * 0.1,
            angular_resolution: analytic.axis_relative_resolution,
        },
        trim,
        max_added_polygons: limit,
    };
    let start = Instant::now();
    let result = match recovery.construct_refined_hybrid(
        &required.iter().copied().collect::<Vec<_>>(),
        &permitted,
        options,
    ) {
        Ok(refined) => {
            let result = refined.recovery;
            let shell = result.complete_shell();
            let exported = if let (Some(path), Some(shell)) = (shell_output, &shell) {
                let mut writer = BufWriter::new(File::create_new(path)?);
                serde_json::to_writer(&mut writer, &shell.compress())?;
                writer.flush()?;
                true
            } else {
                false
            };
            json!({"status":"ok","required_regions":required.len(),"selected_regions":refined.selected_regions,
                "added_polygons":refined.added_polygons,"conflicts":refined.conflicts.iter().map(|c|json!({
                    "support":c.support,"parent_region":c.parent_region,"edge":c.edge,"distance_lower_witness":c.distance_lower_witness
                })).collect::<Vec<_>>(),
                "patches":result.patches.len(),"lifted_trims":result.lifted_trims.len(),
                "failed_regions":result.native.regions.iter().filter(|r|r.faces.is_err()).map(|r|json!({"region":r.region,"error":r.faces.as_ref().err().map(ToString::to_string)})).collect::<Vec<_>>(),
                "failed_trims":result.native.trims.iter().filter(|r|r.is_err()).count(),
                "chart_cuts":result.native.regions.iter().filter(|r|r.requires_chart_cut).count(),
                "approximate_polygons":result.patches.iter().map(|p|p.polygons.len()).sum::<usize>(),
                "max_deviation":result.patches.iter().map(|p|p.deviation_upper_bound).fold(0.0_f64,f64::max),
                "max_seam_deviation":result.patches.iter().map(|p|p.seam_deviation_upper_bound).fold(0.0_f64,f64::max),
                "shell_condition":shell.as_ref().map(|s|format!("{:?}",s.shell_condition())),
                "native_faces":shell.as_ref().map(|s|s.len()),"shell_exported":exported,
                "seconds":start.elapsed().as_secs_f64()})
        }
        Err(error) => {
            json!({"status":"rejected","error":error.to_string(),"seconds":start.elapsed().as_secs_f64()})
        }
    };
    Ok(result)
}

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let path = PathBuf::from(args.next().context("expected OBJ path")?);
    let read_budget = |arg: Option<std::ffi::OsString>| -> Result<f64> {
        arg.context("expected explicit approximation budget")?
            .to_str()
            .context("non-UTF8 budget")?
            .parse()
            .context("numeric budget")
    };
    let vertex = read_budget(args.next())?;
    let chord = read_budget(args.next())?;
    let mut args = args.peekable();
    let shell_output = args
        .peek()
        .is_some_and(|arg| arg != "--release-conflicts" && arg != "--resolve-crossings")
        .then(|| args.next().map(PathBuf::from))
        .flatten();
    let mode = args
        .next()
        .map(|flag| -> Result<_> {
            anyhow::ensure!(
                flag == "--release-conflicts" || flag == "--resolve-crossings",
                "unknown trial mode"
            );
            let limit = args
                .next()
                .context("expected additional polygon limit")?
                .to_str()
                .context("non-UTF8 limit")?
                .parse::<usize>()
                .context("numeric polygon limit")?;
            Ok((flag, limit))
        })
        .transpose()?;
    anyhow::ensure!(
        args.next().is_none()
            && vertex.is_finite()
            && vertex >= 0.0
            && chord.is_finite()
            && chord > 0.0,
        "usage: mesh_recovery_patch_diagnostics OBJ VERTEX_BUDGET CHORD_BUDGET [SHELL_JSON] [--release-conflicts|--resolve-crossings MAX_ADDED_POLYGONS]"
    );
    let start = Instant::now();
    let mesh = obj::read(BufReader::new(File::open(&path)?))?;
    let recovery = recover_analytic_mesh(&mesh, AnalyticRecoveryOptions::default())?;
    let recovery_seconds = start.elapsed().as_secs_f64();
    let budget = SurfaceRegionBudget { vertex, chord };
    let report = match mode {
        Some((flag, limit)) if flag == "--release-conflicts" => {
            refined_trial(&recovery, budget, limit, shell_output.as_deref())?
        }
        Some((_, limit)) => census(&recovery, budget, shell_output.as_deref(), Some(limit))?,
        None => census(&recovery, budget, shell_output.as_deref(), None)?,
    };
    println!(
        "{}",
        json!({"source":path,"recovery_seconds":recovery_seconds,"hybrid":report})
    );
    Ok(())
}
