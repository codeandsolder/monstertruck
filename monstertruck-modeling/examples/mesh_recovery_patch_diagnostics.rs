//! Label-blind feasibility census of bounded hybrid region replacements.
//! This does not modify the strict recovery or claim complete shell repair.
//! Each attempted patch retains all other regions as analytic constraints.

use anyhow::{Context, Result};
use monstertruck_mesh::obj;
use monstertruck_modeling::mesh_recovery::{
    AnalyticRecoveryOptions, RecoveredAnalyticMesh, RegionBoundaryIndex, RegionPatchError,
    RegionPatchOptions, SurfaceCarrier, SurfaceRegionBudget, SurfaceRegionPartition,
    recover_analytic_mesh,
};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::BufReader,
    path::PathBuf,
    time::Instant,
};

fn connected_groups(
    partition: &SurfaceRegionPartition,
    mut candidates: BTreeSet<usize>,
) -> Vec<Vec<usize>> {
    let mut adjacency = BTreeMap::<usize, BTreeSet<usize>>::new();
    let mut touches_curved = BTreeSet::new();
    for path in &partition.contacts.paths {
        let (a, Some(b)) = path.regions else {
            continue;
        };
        match (candidates.contains(&a), candidates.contains(&b)) {
            (true, true) => {
                adjacency.entry(a).or_default().insert(b);
                adjacency.entry(b).or_default().insert(a);
            }
            (true, false)
                if !matches!(partition.carriers[b], Some(SurfaceCarrier::Plane { .. })) =>
            {
                touches_curved.insert(a);
            }
            (false, true)
                if !matches!(partition.carriers[a], Some(SurfaceCarrier::Plane { .. })) =>
            {
                touches_curved.insert(b);
            }
            _ => {}
        }
    }
    let mut groups = Vec::new();
    while let Some(start) = candidates.pop_first() {
        let mut group = Vec::new();
        let mut queue = vec![start];
        while let Some(region) = queue.pop() {
            group.push(region);
            if let Some(neighbors) = adjacency.get(&region) {
                queue.extend(neighbors.iter().copied().filter(|id| candidates.remove(id)));
            }
        }
        if group.iter().any(|id| touches_curved.contains(id)) {
            group.sort_unstable();
            groups.push(group);
        }
    }
    groups
}

struct CandidateSelection {
    groups: Vec<Vec<usize>>,
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
    let initial = connected_groups(partition, failed);
    let required = initial.into_iter().flatten().collect::<BTreeSet<_>>();
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
    let groups = connected_groups(partition, selected.into_iter().collect());
    Ok(CandidateSelection {
        groups,
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

fn census(recovery: &RecoveredAnalyticMesh, budget: SurfaceRegionBudget) -> Result<Value> {
    let index = RegionBoundaryIndex::new(&recovery.evidence)?;
    let partition = &recovery.reconciled.partition;
    let options = RegionPatchOptions {
        budget,
        intersection_tolerance: monstertruck_modeling::KERNEL_EPSILON * 0.1,
        angular_resolution: 1.0e-8,
    };
    let start = Instant::now();
    let CandidateSelection {
        groups,
        required,
        polygons_by_region,
    } = candidate_groups(recovery, &index)?;
    let mut failures = BTreeMap::new();
    let mut results = Vec::new();
    let mut accepted_polygons = 0;
    let mut attempted_polygons = 0;
    let mut accepted = 0;
    let mut accepted_required_polygons = 0;
    let mut connector_polygons = 0;
    for group in &groups {
        let polygons = group.iter().map(|id| polygons_by_region[id]).sum::<usize>();
        attempted_polygons += polygons;
        let required_polygons = group
            .iter()
            .filter(|id| required.contains(id))
            .map(|id| polygons_by_region[id])
            .sum::<usize>();
        connector_polygons += polygons - required_polygons;
        let result = match index.approximate_region_patch(partition, group, options) {
            Ok(patch) => {
                accepted += 1;
                accepted_polygons += patch.source_polygons.len();
                accepted_required_polygons += required_polygons;
                json!({
                    "regions":group,"source_polygons":polygons,"required_polygons":required_polygons,
                    "native_faces":patch.approximation.shell.len(),
                    "boundary_edges":patch.boundary_edges.len(),
                    "vertex_displacement":patch.approximation.vertex_displacement,
                    "deviation_bound":patch.approximation.deviation_upper_bound,
                })
            }
            Err(error) => {
                *failures.entry(failure_class(error)).or_insert(0_usize) += 1;
                json!({"regions":group,"source_polygons":polygons,"required_polygons":required_polygons,"error":error.to_string()})
            }
        };
        results.push(result);
    }
    Ok(json!({
        "vertex_budget":budget.vertex,"chord_budget":budget.chord,
        "attempted_groups":groups.len(),"attempted_polygons":attempted_polygons,
        "accepted_groups":accepted,"accepted_polygons":accepted_polygons,
        "accepted_required_polygons":accepted_required_polygons,"connector_polygons":connector_polygons,
        "failures":failures,"seconds":start.elapsed().as_secs_f64(),"groups":results
    }))
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
    anyhow::ensure!(
        args.next().is_none()
            && vertex.is_finite()
            && vertex >= 0.0
            && chord.is_finite()
            && chord > 0.0,
        "usage: mesh_recovery_patch_diagnostics OBJ VERTEX_BUDGET CHORD_BUDGET"
    );
    let start = Instant::now();
    let mesh = obj::read(BufReader::new(File::open(&path)?))?;
    let recovery = recover_analytic_mesh(&mesh, AnalyticRecoveryOptions::default())?;
    let recovery_seconds = start.elapsed().as_secs_f64();
    let report = census(&recovery, SurfaceRegionBudget { vertex, chord })?;
    println!(
        "{}",
        json!({"source":path,"recovery_seconds":recovery_seconds,"hybrid":report})
    );
    Ok(())
}
