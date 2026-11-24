//! Minimal, label-blind batch benchmark for the public recovery pipeline.
//! Independent oracle evaluation belongs in the corpus harness, not this binary.
//! Each invocation consumes one OBJ and emits one machine-readable record.
use anyhow::{Context, Result};
use monstertruck_mesh::obj;
use monstertruck_modeling::mesh_recovery::{
    AnalyticRecoveryOptions, RecoveredAnalyticMesh, SharedTrimGeometry, SurfaceCarrier,
    analyze_mesh_evidence, recover_analytic_evidence,
};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
    time::Instant,
};

fn carrier_summary(carrier: SurfaceCarrier) -> Value {
    match carrier {
        SurfaceCarrier::Plane { origin, normal } => json!({
            "kind":"plane","origin":[origin.x,origin.y,origin.z],
            "normal":[normal.x,normal.y,normal.z],
        }),
        SurfaceCarrier::LinearRevolution {
            origin,
            axis,
            radius_at_origin,
            slope,
        } => json!({
            "kind":"linear","origin":[origin.x,origin.y,origin.z],
            "axis":[axis.x,axis.y,axis.z],"radius":radius_at_origin,"slope":slope,
        }),
        SurfaceCarrier::CircularRevolution {
            axial_center,
            axis,
            radial_center,
            profile_radius,
        } => json!({
            "kind":"circular","origin":[axial_center.x,axial_center.y,axial_center.z],
            "axis":[axis.x,axis.y,axis.z],"radial_center":radial_center,"profile_radius":profile_radius,
        }),
    }
}

fn trim_summary(
    trims: &[monstertruck_modeling::mesh_recovery::SharedAnalyticTrim],
    carriers: &[Option<SurfaceCarrier>],
) -> ([usize; 4], BTreeMap<String, usize>) {
    let mut trim_edges = [0_usize; 4];
    let mut unresolved_contact_pairs = BTreeMap::new();
    for trim in trims {
        let kind = match trim.geometry {
            SharedTrimGeometry::Line(_) => 0,
            SharedTrimGeometry::Circle(_) => 1,
            SharedTrimGeometry::Conic(_) => 3,
            SharedTrimGeometry::Unresolved => 2,
        };
        trim_edges[kind] += trim.source_edges.len();
        if kind == 2 {
            let region_kind = |id| match carriers.get(id).copied().flatten() {
                Some(SurfaceCarrier::Plane { .. }) => "plane",
                Some(SurfaceCarrier::LinearRevolution { .. }) => "linear",
                Some(SurfaceCarrier::CircularRevolution { radial_center, .. })
                    if radial_center.abs() <= monstertruck_modeling::KERNEL_EPSILON =>
                {
                    "sphere"
                }
                Some(SurfaceCarrier::CircularRevolution { .. }) => "circular",
                None => "missing",
            };
            let mut pair = [
                region_kind(trim.regions.0),
                trim.regions.1.map_or("open", region_kind),
            ];
            pair.sort_unstable();
            *unresolved_contact_pairs
                .entry(pair.join("+"))
                .or_insert(0_usize) += trim.source_edges.len();
        }
    }
    (trim_edges, unresolved_contact_pairs)
}

// Preserve the exact source samples and provenance of every unresolved trim
// so independent diagnostics need not repeat surface recognition.
fn unresolved_trims(recovery: &RecoveredAnalyticMesh) -> Vec<Value> {
    let provenance = &recovery.evidence.provenance.original_polygon_ids;
    recovery.trims.trims.iter().filter(|trim| matches!(trim.geometry, SharedTrimGeometry::Unresolved)).map(|trim| {
        let edges = trim.source_edges.iter().map(|edge| {
            let vertices = edge.vertices.map(|id| {
                let p = recovery.reconciled.positions[id];
                [p.x, p.y, p.z]
            });
            json!({"vertices":edge.vertices,"points":vertices,"polygons":[Some(provenance[edge.interior_polygon]),edge.exterior_polygon.map(|id|provenance[id])]})
        }).collect::<Vec<_>>();
        json!({"contact_path":trim.contact_path,"regions":trim.regions,"edges":edges})
    }).collect()
}

fn benchmark(path: &Path, evidence_only: bool, normalize: bool) -> Result<Value> {
    let start = Instant::now();
    let mesh = obj::read(BufReader::new(File::open(path)?)).context("read OBJ")?;
    let read_seconds = start.elapsed().as_secs_f64();
    let options = AnalyticRecoveryOptions {
        normalize_coincident_sheets: normalize,
        ..Default::default()
    };
    let evidence = analyze_mesh_evidence(&mesh, options.evidence)?;
    let evidence_seconds = start.elapsed().as_secs_f64() - read_seconds;
    if evidence_only {
        return Ok(json!({
            "status":"ok","source":path,"source_polygons":mesh.face_iter().count(),
            "read_seconds":read_seconds,"evidence_seconds":evidence_seconds,
            "source_supports":evidence.planar_supports.len(),
            "collapsed_source_facets":evidence.provenance.collapsed_facets.iter().map(|f|json!({"polygon":f.polygon,"vertices":f.vertices,"ordered_vertices":f.ordered_vertices})).collect::<Vec<_>>(),
            "original_polygon_ids":evidence.provenance.original_polygon_ids,
            "separated_source_vertices":evidence.provenance.separated_vertices.len(),
            "source_sector_pairs":evidence.provenance.sector_pairs.len(),
            "source_nonmanifold_edges":evidence.non_manifold_edge_count,
            "source_winding_conflicts":evidence.inconsistent_winding_edge_count,
            "source_signed_volume":evidence.closed_source_volume(),
        }));
    }
    let recovery = recover_analytic_evidence(&evidence, options)?;
    let partition = &recovery.reconciled.partition;
    let mut errors = BTreeMap::new();
    let mut failed_regions = 0;
    let mut native_faces = 0;
    let mut chart_cuts = 0;
    let mut region_info = Vec::new();
    for region in &recovery.native.regions {
        let carrier = partition
            .carriers
            .get(region.region)
            .copied()
            .flatten()
            .context("missing carrier for recovered region")?;
        let (faces, error) = match &region.faces {
            Ok(faces) => {
                native_faces += faces.len();
                (faces.len(), None)
            }
            Err(error) => {
                failed_regions += 1;
                *errors.entry(error.to_string()).or_insert(0_usize) += 1;
                (0, Some(error.to_string()))
            }
        };
        chart_cuts += usize::from(region.requires_chart_cut);
        region_info.push(
            json!({"id":region.region,"carrier":carrier_summary(carrier),
            "native_faces":faces,"pending_chart_cut":region.requires_chart_cut,"error":error}),
        );
    }
    let canceled_sheets = recovery
        .evidence
        .provenance
        .canceled_sheets
        .iter()
        .map(|pair| {
            json!({"sides":pair.sides,"resolved_edges":pair.resolved_edges,
            "max_plane_residual":pair.max_plane_residual})
        })
        .collect::<Vec<_>>();
    let (trim_edges, unresolved_contact_pairs) =
        trim_summary(&recovery.trims.trims, &partition.carriers);
    let complete = recovery.complete_shell();
    Ok(json!({
        "status":"ok","source":path,"source_polygons":mesh.face_iter().count(),
        "read_seconds":read_seconds,"evidence_seconds":evidence_seconds,
        "total_seconds":start.elapsed().as_secs_f64(),
        "regions":region_info,"region_by_polygon":partition.region_by_polygon,
        "original_polygon_ids":recovery.evidence.provenance.original_polygon_ids,
        "collapsed_source_facets":recovery.evidence.provenance.collapsed_facets.iter().map(|f|json!({"polygon":f.polygon,"vertices":f.vertices,"ordered_vertices":f.ordered_vertices})).collect::<Vec<_>>(),
        "separated_source_vertices":recovery.evidence.provenance.separated_vertices.iter().map(|v|json!([v.source_vertex,v.separated_vertex])).collect::<Vec<_>>(),
        "source_sector_pairs":recovery.evidence.provenance.sector_pairs.iter().map(|p|json!({"edge":p.edge,"polygons":p.polygons,"positive_volume":p.positive_volume})).collect::<Vec<_>>(),
        "canceled_sheets":canceled_sheets,
        "normalize_coincident_sheets":normalize,
        "failed_regions":failed_regions,"native_faces":native_faces,"pending_chart_cuts":chart_cuts,
        "native_region_errors":errors,
        "native_trim_errors":recovery.native.trims.iter().filter(|trim|trim.is_err()).count(),
        "ambiguous_source_edges":recovery.trims.ambiguous_edges.len(),
        "source_nonmanifold_edges":recovery.evidence.non_manifold_edge_count,
        "source_winding_conflicts":recovery.evidence.inconsistent_winding_edge_count,
        "source_signed_volume":recovery.evidence.closed_source_volume(),
        "line_source_edges":trim_edges[0],"circle_source_edges":trim_edges[1],"conic_source_edges":trim_edges[3],
        "unresolved_source_edges":trim_edges[2],
        "unresolved_contact_pairs":unresolved_contact_pairs,
        "unresolved_trims":unresolved_trims(&recovery),
        "trim_absolute_tolerance":options.trim_options(&recovery.evidence).absolute_tolerance,
        "trim_relative_tolerance":options.trim_options(&recovery.evidence).relative_tolerance,
        "limited_partition_components":partition.limited_components,
        "rim_linear_proposals":recovery.discovery.rim_linear_proposals,
        "compared_native_neighborhoods":recovery.discovery.compared_native_neighborhoods,
        "expanded_native_neighborhoods":recovery.discovery.expanded_native_neighborhoods,
        "shell_condition":complete.as_ref().map(|shell|format!("{:?}",shell.shell_condition())),
        "shell_components":complete.as_ref().map(|shell|shell.connected_components().len()),
    }))
}

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let first = args.next().context("expected one OBJ path")?;
    let evidence_only = first == "--evidence-only";
    let normalize = first == "--normalize-coincident-sheets";
    let path = PathBuf::from(if evidence_only || normalize {
        args.next().context("expected OBJ path after option")?
    } else {
        first
    });
    anyhow::ensure!(
        args.next().is_none(),
        "expected [--evidence-only | --normalize-coincident-sheets] OBJ"
    );
    let result = benchmark(&path, evidence_only, normalize).unwrap_or_else(|error| {
        json!({
            "status":"recovery_error","source":path,"error":format!("{error:#}"),
        })
    });
    println!("{result}");
    Ok(())
}
