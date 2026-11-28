//! Minimal, label-blind batch benchmark for the public recovery pipeline.
//! Independent oracle evaluation belongs in the corpus harness, not this binary.
//! Each invocation consumes one OBJ and emits one machine-readable record.
use anyhow::{Context, Result};
use monstertruck_mesh::obj;
use monstertruck_modeling::mesh_recovery::{
    AnalyticRecoveryOptions, FeatureEdgeKind, MeshEvidence, RecoveredAnalyticMesh,
    SharedAnalyticTrim, SharedTrimGeometry, SurfaceCarrier, analyze_mesh_evidence,
    partition_recovery_mesh, recover_analytic_evidence,
};
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
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
) -> ([usize; 6], BTreeMap<String, usize>) {
    let mut trim_edges = [0_usize; 6];
    let mut unresolved_contact_pairs = BTreeMap::new();
    for trim in trims {
        let kind = match trim.geometry {
            SharedTrimGeometry::Line(_) => 0,
            SharedTrimGeometry::Circle(_) => 1,
            SharedTrimGeometry::Conic(_) => 3,
            SharedTrimGeometry::ToroidalSection(_) => 4,
            SharedTrimGeometry::LinearSection(_) => 5,
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

// Preserve the original, reconciled and native endpoint coordinates and
// original polygon provenance, without rerunning region discovery.
fn trim_edges(recovery: &RecoveredAnalyticMesh, trim: &SharedAnalyticTrim) -> Vec<Value> {
    let provenance = &recovery.evidence.provenance.original_polygon_ids;
    trim.source_edges.iter().map(|edge| {
        let source = edge.vertices.map(|id| {
            let p = recovery.evidence.positions[id];
            [p.x,p.y,p.z]
        });
        let points = edge.vertices.map(|id| {
            let p = recovery.reconciled.positions[id];
            [p.x,p.y,p.z]
        });
        let native = edge.vertices.map(|id| {
            recovery.native.vertices[id].as_ref().map(|vertex| {
                let p = vertex.point();
                [p.x,p.y,p.z]
            })
        });
        json!({"vertices":edge.vertices,"source":source,"points":points,"native":native,
            "polygons":[Some(provenance[edge.interior_polygon]),edge.exterior_polygon.map(|id|provenance[id])]})
    }).collect()
}

// Separate unresolved analytical hypotheses from accepted trim geometries
// that could not be constructed by the native kernel.
fn trim_diagnostics(recovery: &RecoveredAnalyticMesh) -> (Vec<Value>, Vec<Value>) {
    let mut unresolved = Vec::new();
    let mut failed = Vec::new();
    for (id, (trim, native)) in recovery
        .trims
        .trims
        .iter()
        .zip(&recovery.native.trims)
        .enumerate()
    {
        if matches!(trim.geometry, SharedTrimGeometry::Unresolved) {
            unresolved.push(
                json!({"contact_path":trim.contact_path,"regions":trim.regions,
                "edges":trim_edges(recovery,trim)}),
            );
        } else if let Err(error) = native {
            failed.push(
                json!({"trim":id,"contact_path":trim.contact_path,"regions":trim.regions,
                "geometry":format!("{:?}",trim.geometry),"error":error.to_string(),
                "edges":trim_edges(recovery,trim)}),
            );
        }
    }
    (unresolved, failed)
}

// Preserve the complete post-welding incidence of unresolved source
// edges; these defects precede analytical carrier and trim selection.
fn ambiguous_source_incidence(evidence: &MeshEvidence) -> Vec<Value> {
    if evidence.non_manifold_edge_count == 0 && evidence.inconsistent_winding_edge_count == 0 {
        return Vec::new();
    }
    let mut components = vec![usize::MAX; evidence.planar_supports.len()];
    for seed in 0..components.len() {
        if components[seed] != usize::MAX {
            continue;
        }
        components[seed] = seed;
        let mut stack = vec![seed];
        while let Some(support) = stack.pop() {
            for &interface in &evidence.support_interface_indices[support] {
                for &other in &evidence.support_interfaces[interface].supports {
                    if components[other] == usize::MAX {
                        components[other] = seed;
                        stack.push(other);
                    }
                }
            }
        }
    }
    let point = |id: usize| {
        let p = evidence.positions[id];
        [p.x, p.y, p.z]
    };
    evidence
        .edges
        .iter()
        .filter(|edge| {
            matches!(
                edge.kind,
                FeatureEdgeKind::NonManifoldJunction | FeatureEdgeKind::WindingConflict
            )
        })
        .map(|edge| {
            let uses = edge
                .source_face_indices
                .iter()
                .zip(&edge.source_forward_winding)
                .map(|(&polygon, &forward)| {
                    let vertices = evidence
                        .source_polygon_vertices(polygon)
                        .unwrap_or_default();
                    json!({"polygon":evidence.provenance.original_polygon_ids[polygon],
                "forward":forward,"regular_component":components[evidence.source_planar_support[polygon]],"vertices":vertices,
                "points":vertices.iter().copied().map(point).collect::<Vec<_>>()})
                })
                .collect::<Vec<_>>();
            json!({"vertices":edge.vertices,"points":edge.vertices.map(point),"uses":uses})
        })
        .collect()
}

struct RegionSummary {
    errors: BTreeMap<String, usize>,
    failed: usize,
    faces: usize,
    chart_cuts: usize,
    info: Vec<Value>,
}

fn planar_crossing_diagnostics(recovery: &RecoveredAnalyticMesh) -> Result<Vec<Value>> {
    let partition = &recovery.reconciled.partition;
    Ok(recovery.native.planar_crossings(partition,&recovery.trims,monstertruck_modeling::KERNEL_EPSILON*0.1)?
        .into_iter().map(|crossing| {
            let edges = crossing.edges.map(|(trim,ordinal)| {
                let t = &recovery.trims.trims[trim];
                let native = recovery.native.trims[trim].as_ref().ok().and_then(|edges|edges.get(ordinal));
                json!({"trim":trim,"ordinal":ordinal,"regions":t.regions,"geometry":format!("{:?}",t.geometry),
                    "native":native.map(|e|json!({"front":e.front().point(),"back":e.back().point(),"curve":format!("{:?}",e.curve())})),
                    "source_edges":t.source_edges.iter().map(|e| {
                        let supports = [Some(e.interior_polygon),e.exterior_polygon].map(|p|p.map(|id|recovery.evidence.source_planar_support[id]));
                        json!({"vertices":e.vertices,"points":e.vertices.map(|id|recovery.reconciled.positions[id]),
                            "supports":supports,"sizes":supports.map(|s|s.map(|id|recovery.evidence.planar_supports[id].source_face_indices.len()))})
                    }).collect::<Vec<_>>()})
            });
            json!({"region":crossing.region,"edges":edges})
        }).collect())
}

fn region_summary(recovery: &RecoveredAnalyticMesh) -> Result<RegionSummary> {
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
            "native_faces":faces,"pending_chart_cut":region.requires_chart_cut,"chart_obstacles":region.chart_obstacles.iter().map(|obstacle|format!("{obstacle:?}")).collect::<Vec<_>>(),"error":error}),
        );
    }
    Ok(RegionSummary {
        errors,
        failed: failed_regions,
        faces: native_faces,
        chart_cuts,
        info: region_info,
    })
}

fn evidence_summary(
    evidence: &MeshEvidence,
    path: &Path,
    source_polygons: usize,
    read_seconds: f64,
    evidence_seconds: f64,
) -> Value {
    json!({
            "status":"ok","source":path,"source_polygons":source_polygons,
            "read_seconds":read_seconds,"evidence_seconds":evidence_seconds,
            "source_supports":evidence.planar_supports.len(),
            "collapsed_source_facets":evidence.provenance.collapsed_facets.iter().map(|f|json!({"polygon":f.polygon,"vertices":f.vertices,"ordered_vertices":f.ordered_vertices})).collect::<Vec<_>>(),
            "original_polygon_ids":evidence.provenance.original_polygon_ids,
            "separated_source_vertices":evidence.provenance.separated_vertices.len(),
            "source_sector_pairs":evidence.provenance.sector_pairs.len(),
            "source_forced_pairs":evidence.provenance.forced_pairs.len(),
            "source_multiplicity_resolutions":evidence.provenance.multiplicity_resolutions.len(),
            "source_nonmanifold_edges":evidence.non_manifold_edge_count,
            "ambiguous_source_incidence":ambiguous_source_incidence(evidence),
            "source_winding_conflicts":evidence.inconsistent_winding_edge_count,
            "source_signed_volume":evidence.closed_source_volume(),
    })
}

fn write_shell_output(
    path: &Path,
    complete: Option<&monstertruck_modeling::Shell>,
    result: &mut Value,
) -> Result<()> {
    let export_start = Instant::now();
    let exported = if let Some(shell) = complete {
        let mut writer = BufWriter::new(File::create_new(path)?);
        serde_json::to_writer(&mut writer, &shell.compress())?;
        writer.flush()?;
        true
    } else {
        false
    };
    result["shell_exported"] = json!(exported);
    result["shell_export_seconds"] = json!(export_start.elapsed().as_secs_f64());
    Ok(())
}

fn benchmark(
    path: &Path,
    evidence_only: bool,
    normalize: bool,
    joint_fit: bool,
    shell_output: Option<&Path>,
) -> Result<Value> {
    let start = Instant::now();
    let mesh = obj::read(BufReader::new(File::open(path)?)).context("read OBJ")?;
    let read_seconds = start.elapsed().as_secs_f64();
    let mut options = AnalyticRecoveryOptions {
        normalize_coincident_sheets: normalize,
        ..Default::default()
    };
    if !joint_fit {
        options.max_joint_fit_trials = 0;
    }
    let evidence = analyze_mesh_evidence(&mesh, options.evidence)?;
    let evidence_seconds = start.elapsed().as_secs_f64() - read_seconds;
    if evidence_only {
        return Ok(evidence_summary(
            &evidence,
            path,
            mesh.face_iter().count(),
            read_seconds,
            evidence_seconds,
        ));
    }
    let recovery = recover_analytic_evidence(&evidence, &options)?;
    let partition = &recovery.reconciled.partition;
    let regions = region_summary(&recovery)?;
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
    let (unresolved_trims, failed_native_trims) = trim_diagnostics(&recovery);
    let mut result = json!({
        "status":"ok","source":path,"source_polygons":mesh.face_iter().count(),
        "read_seconds":read_seconds,"evidence_seconds":evidence_seconds,
        "total_seconds":start.elapsed().as_secs_f64(),
        "regions":regions.info,"region_by_polygon":partition.region_by_polygon,
        "coplanar_reassignments":recovery.coplanar_reassignments.iter().map(|r|json!({"support":r.support,"old_region":r.old_region,"new_region":r.new_region,"max_plane_residual":r.max_plane_residual})).collect::<Vec<_>>(),
        "joint_fit_rejections":recovery.joint_fit_rejections.iter().map(|r|json!({"regions":r.regions,"method":format!("{:?}",r.method),"rejection":format!("{:?}",r.rejection)})).collect::<Vec<_>>(),
        "joint_fits":recovery.joint_fits.iter().map(|r|json!({"method":format!("{:?}",r.method),"regions":r.regions,"refitted_regions":r.refitted_regions,"repaired_edges":r.repaired_edges,"ownership_additions":r.ownership_additions,"max_displacement":r.max_displacement,"max_allowance":r.max_allowance,"max_additional_carrier_residual":r.max_additional_carrier_residual})).collect::<Vec<_>>(),
        "original_polygon_ids":recovery.evidence.provenance.original_polygon_ids,
        "collapsed_source_facets":recovery.evidence.provenance.collapsed_facets.iter().map(|f|json!({"polygon":f.polygon,"vertices":f.vertices,"ordered_vertices":f.ordered_vertices})).collect::<Vec<_>>(),
        "separated_source_vertices":recovery.evidence.provenance.separated_vertices.iter().map(|v|json!([v.source_vertex,v.separated_vertex])).collect::<Vec<_>>(),
        "source_forced_pairs":recovery.evidence.provenance.forced_pairs.iter().map(|p|json!({"edge":p.edge,"polygons":p.polygons})).collect::<Vec<_>>(),
        "source_multiplicity_resolutions":recovery.evidence.provenance.multiplicity_resolutions.iter().map(|r|json!({"method":format!("{:?}",r.method),"states":r.states,"equivalent_fragments":r.equivalent_fragments,"opposed_faces":r.opposed_faces,"transverse_crossings":r.transverse_crossings.iter().map(|w|json!({"polygons":w.polygons,"edge":w.edge})).collect::<Vec<_>>(),"pairs":r.pairs.iter().map(|p|json!({"edge":p.edge,"polygons":p.polygons})).collect::<Vec<_>>()})).collect::<Vec<_>>(),
        "source_sector_pairs":recovery.evidence.provenance.sector_pairs.iter().map(|p|json!({"edge":p.edge,"polygons":p.polygons,"positive_volume":p.positive_volume})).collect::<Vec<_>>(),
        "canceled_sheets":canceled_sheets,
        "normalize_coincident_sheets":normalize,
        "failed_regions":regions.failed,"native_faces":regions.faces,"pending_chart_cuts":regions.chart_cuts,
        "native_region_errors":regions.errors,
        "native_trim_errors":recovery.native.trims.iter().filter(|trim|trim.is_err()).count(),
        "ambiguous_source_edges":recovery.trims.ambiguous_edges.len(),
        "source_nonmanifold_edges":recovery.evidence.non_manifold_edge_count,
        "ambiguous_source_incidence":ambiguous_source_incidence(&recovery.evidence),
        "source_winding_conflicts":recovery.evidence.inconsistent_winding_edge_count,
        "source_signed_volume":recovery.evidence.closed_source_volume(),
    });
    let trim_summary = json!({
        "line_source_edges":trim_edges[0],"circle_source_edges":trim_edges[1],"conic_source_edges":trim_edges[3],"toroidal_section_source_edges":trim_edges[4],"linear_section_source_edges":trim_edges[5],
        "unresolved_source_edges":trim_edges[2],
        "unresolved_contact_pairs":unresolved_contact_pairs,
        "unresolved_trims":unresolved_trims,
        "failed_native_trims":failed_native_trims,
        "planar_crossings":planar_crossing_diagnostics(&recovery)?,
        "trim_absolute_tolerance":options.trim_options(&recovery.evidence).absolute_tolerance,
        "trim_relative_tolerance":options.trim_options(&recovery.evidence).relative_tolerance,
        "limited_partition_components":partition.limited_components,
        "rim_linear_proposals":recovery.discovery.rim_linear_proposals,
        "expanded_circular_proposals":recovery.discovery.expanded_circular_proposals,"expanded_linear_proposals":recovery.discovery.expanded_linear_proposals,
        "compared_native_neighborhoods":recovery.discovery.compared_native_neighborhoods,
        "expanded_native_neighborhoods":recovery.discovery.expanded_native_neighborhoods,
        "carrier_alternative_trials":recovery.discovery.carrier_alternative_trials,
        "carrier_alternative_acceptances":recovery.discovery.carrier_alternative_acceptances,
        "shell_condition":complete.as_ref().map(|shell|format!("{:?}",shell.shell_condition())),
        "shell_components":complete.as_ref().map(|shell|shell.connected_components().len()),
    });
    result
        .as_object_mut()
        .context("benchmark result is not an object")?
        .extend(
            trim_summary
                .as_object()
                .context("trim summary is not an object")?
                .iter()
                .map(|(key, value)| (key.clone(), value.clone())),
        );
    if let Some(path) = shell_output {
        write_shell_output(path, complete.as_ref(), &mut result)?;
    }
    Ok(result)
}

fn partition_summary(path: &Path) -> Result<Value> {
    let mesh = obj::read(BufReader::new(File::open(path)?)).context("read OBJ")?;
    let start = Instant::now();
    let parts = partition_recovery_mesh(&mesh)?;
    Ok(json!({
        "status": "ok",
        "source": path,
        "partition_seconds": start.elapsed().as_secs_f64(),
        "source_polygons": mesh.face_iter().count(),
        "components": parts.iter().map(|part| json!({
            "source_polygons": part.original_polygons.len(),
            "independently_closed": part.independently_closed,
            "original_polygon_ids": part.original_polygons,
        })).collect::<Vec<_>>(),
    }))
}

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let first = args.next().context("expected one OBJ path")?;
    let partition_only = first == "--partition-only";
    let evidence_only = first == "--evidence-only";
    let normalize = first == "--normalize-coincident-sheets";
    let joint_fit = first != "--no-joint-fit";
    let path = PathBuf::from(
        if partition_only || evidence_only || normalize || !joint_fit {
            args.next().context("expected OBJ path after option")?
        } else {
            first
        },
    );
    let shell_output = args
        .next()
        .map(|flag| -> Result<PathBuf> {
            anyhow::ensure!(
                flag == "--shell-output" && !evidence_only && !partition_only,
                "expected --shell-output after OBJ"
            );
            Ok(PathBuf::from(
                args.next().context("expected shell output path")?,
            ))
        })
        .transpose()?;
    anyhow::ensure!(
        args.next().is_none(),
        "expected [--partition-only | --evidence-only | --normalize-coincident-sheets | --no-joint-fit] OBJ [--shell-output JSON]"
    );
    if partition_only {
        println!("{}", partition_summary(&path)?);
        return Ok(());
    }
    let result = benchmark(
        &path,
        evidence_only,
        normalize,
        joint_fit,
        shell_output.as_deref(),
    )
    .unwrap_or_else(|error| {
        json!({
            "status":"recovery_error","source":path,"error":format!("{error:#}"),
        })
    });
    println!("{result}");
    Ok(())
}
