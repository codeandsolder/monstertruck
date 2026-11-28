//! End-to-end native mesh recovery verification, independent of CAD labels.
//! Volume measurements are numerical diagnostics, not an intersection proof.
//! STEP output is held in memory unless explicitly requested.

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use monstertruck_io::step::{
    load::Table,
    save::{CompleteStepDisplay, StepHeaderDescriptor, StepModel, StepModels},
};
use monstertruck_mesh::{PolygonMesh, obj};
use monstertruck_meshing::{
    analyzers::{CalcVolume, Topology},
    tessellation::{
        MeshableShape, face_drop_count, reset_face_drop_count, stitch_compressed_shell,
    },
};
use monstertruck_modeling::{
    Curve, EuclideanSpace, KERNEL_EPSILON, Point3, Shell, Solid, Surface,
    ToCompressedTrimmedParameterCurves,
    mesh_recovery::{AnalyticRecoveryOptions, recover_analytic_mesh},
};
use monstertruck_topology::{
    compress::CompressedShell, shell::ShellCondition, trimmed::TrimmedShell,
};
use serde_json::{Value, json};
use std::{
    fs::File,
    io::{BufReader, Write},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum InputKind {
    #[default]
    Obj,
    NativeShell,
    Step,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum VerificationMode {
    #[default]
    Full,
    Topology,
    Export,
}

/// Recover native B-rep regions, verify strict tessellation and STEP reload.
#[derive(Parser, Debug)]
struct Args {
    /// Maximum explicitly allowed planar repair displacement in source units.
    #[arg(long, default_value_t = 0.0)]
    max_planar_displacement: f64,
    /// Maximum additional source adjustment in local joint reconstruction.
    #[arg(long, default_value_t = 0.01)]
    max_reconstruction_adjustment: f64,
    /// Known positional uncertainty of source samples, in source units.
    #[arg(long, default_value_t = 0.0)]
    source_position_uncertainty: f64,
    /// Known source tessellation deflection in model units (not kernel tolerance).
    #[arg(long)]
    source_chord_budget: Option<f64>,
    /// Explicitly permit cancellation of proven reverse coincident sheets.
    #[arg(long)]
    normalize_coincident_sheets: bool,
    /// First tessellation tolerance in source units; second is four times finer.
    #[arg(long, default_value_t = 0.01)]
    tolerance: f64,
    /// Input representation: source OBJ, compressed native shell JSON or STEP.
    #[arg(long, value_enum, default_value_t)]
    input: InputKind,
    /// Write the STEP exchange to this path (one input only).
    #[arg(long)]
    step_output: Option<PathBuf>,
    /// Full numerical verification, topology-only reload, or export without reload.
    #[arg(long, value_enum, default_value_t)]
    verification: VerificationMode,
    /// OBJ, compressed native shell JSON or STEP inputs.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,
}

fn numerical_volume(shell: &Shell, tolerance: f64) -> Value {
    let _ = reset_face_drop_count();
    let mesh = shell.compress().triangulation(tolerance);
    match stitch_compressed_shell(&mesh, KERNEL_EPSILON) {
        Ok(mesh) if face_drop_count() == 0 => measure_polygon_volume(mesh, tolerance),
        Ok(_) => {
            json!({"tolerance":tolerance,"error":"empty tessellation face","dropped_faces":face_drop_count()})
        }
        Err(error) => {
            json!({"tolerance":tolerance,"error":error.to_string(),"dropped_faces":face_drop_count()})
        }
    }
}

fn measure_polygon_volume(mut mesh: PolygonMesh, tolerance: f64) -> Value {
    if let Some(&origin) = mesh.positions().first() {
        for point in mesh.positions_mut() {
            *point -= origin.to_vec();
        }
    }
    let condition = mesh.shell_condition();
    let signed_integral = mesh.volume();
    json!({"tolerance":tolerance,"dropped_faces":0,
        "polygons":mesh.face_iter().count(),"mesh_condition":format!("{condition:?}"),
        "signed_integral":signed_integral,
        "closed_signed_volume":(condition==ShellCondition::Closed).then_some(signed_integral)})
}

fn step_roundtrip(shell: &Shell, solid_eligible: bool, args: &Args) -> Result<Value> {
    let step = if solid_eligible {
        let solid = Solid::new(vec![shell.clone()])?;
        let compressed = solid.compress_with_parameter_curves(KERNEL_EPSILON * 100.0);
        CompleteStepDisplay::new(
            StepModel::from(&compressed),
            StepHeaderDescriptor::default(),
        )
        .to_string()
    } else {
        // A disconnected collection is not one connected STEP shell.
        let compressed = shell
            .connected_components()
            .iter()
            .map(|component| component.compress_with_parameter_curves(KERNEL_EPSILON * 100.0))
            .collect::<Vec<_>>();
        let models = compressed.iter().collect::<StepModels<_, _, _>>();
        CompleteStepDisplay::new(models, StepHeaderDescriptor::default()).to_string()
    };
    let table = (args.verification != VerificationMode::Export)
        .then(|| Table::from_step(&step))
        .transpose()?;
    if let Some(path) = args.step_output.as_deref() {
        let mut file = File::create_new(path).context("create STEP output")?;
        file.write_all(step.as_bytes())?;
        file.flush()?;
    }
    let Some(table) = table else {
        return Ok(json!({"bytes":step.len(),"validation":"not_requested","solid_brep":false}));
    };
    inspect_step(
        &table,
        step.len(),
        Some(shell.len()),
        args.tolerance * 0.25,
        args.verification == VerificationMode::Full
            && shell.shell_condition() == ShellCondition::Closed,
    )
}

fn inspect_step(
    table: &Table,
    bytes: usize,
    expected_faces: Option<usize>,
    tolerance: f64,
    tessellate: bool,
) -> Result<Value> {
    let mut restored_faces = 0;
    let mut restored_edges = 0;
    let mut restored_conditions = Vec::new();
    let mut restored_endpoints_consistent = true;
    let mut restored_parameter_curves = 0;
    let mut restored_volumes = Vec::new();
    for native in table.shell.values() {
        let restored = table.to_compressed_trimmed_shell(native)?;
        restored_faces += restored.faces.len();
        restored_edges += restored.edges.len();
        // Tessellate with the reloaded exact face-local trims, not a second
        // parameter projection from the original native boundary.
        if tessellate {
            let _ = reset_face_drop_count();
            let meshed = restored.triangulation(tolerance);
            let dropped = face_drop_count();
            if dropped != 0
                || meshed.faces.iter().any(|face| {
                    face.surface
                        .as_ref()
                        .is_none_or(|mesh| mesh.face_iter().next().is_none())
                })
            {
                restored_volumes.push(json!({"tolerance":tolerance,"dropped_faces":dropped,"error":"incomplete restored tessellation"}));
            } else {
                restored_volumes.push(match stitch_compressed_shell(&meshed, KERNEL_EPSILON) {
                    Ok(mesh) => measure_polygon_volume(mesh, tolerance),
                    Err(error) => json!({"tolerance":tolerance,"error":error.to_string(),"dropped_faces":dropped}),
                });
            }
        }
        let restored = TrimmedShell::try_from(restored)?;
        restored_parameter_curves += restored
            .faces()
            .iter()
            .flat_map(|face| face.trims().iter().flatten())
            .filter(|trim| trim.is_some())
            .count();
        let plain = restored.erase_trims();
        restored_conditions.push(format!("{:?}", plain.shell_condition()));
        restored_endpoints_consistent &=
            plain.edge_iter().all(|edge| edge.is_geometric_consistent());
    }
    Ok(
        json!({"bytes":bytes,"faces":restored_faces,"edges":restored_edges,
        "face_count_preserved":expected_faces.map(|expected| restored_faces == expected),"restored_tessellation_requested":tessellate,"solid_breps":table.manifold_solid_brep.len(),
        "restored_shell_conditions":restored_conditions,"restored_endpoints_consistent":restored_endpoints_consistent,
        "restored_parameter_curves":restored_parameter_curves,"restored_volumes":restored_volumes,
        "planes":table.plane.len(),"cylinders":table.cylindrical_surface.len(),"cones":table.conical_surface.len(),
        "tori":table.toroidal_surface.len(),"spheres":table.spherical_surface.len(),
        "b_spline_surfaces":table.b_spline_surface_with_knots.len()+table.rational_b_spline_surface.len()}),
    )
}

fn verify(path: &Path, args: &Args) -> Result<Value> {
    let mesh = obj::read(BufReader::new(File::open(path)?)).context("read OBJ")?;
    let result = recover_analytic_mesh(
        &mesh,
        AnalyticRecoveryOptions {
            max_planar_displacement: args.max_planar_displacement,
            max_reconstruction_adjustment: args.max_reconstruction_adjustment,
            source_position_uncertainty: args.source_position_uncertainty,
            source_chord_budget: args.source_chord_budget,
            normalize_coincident_sheets: args.normalize_coincident_sheets,
            ..Default::default()
        },
    )?;
    let mut summary = json!({"source":path,"source_polygons":mesh.face_iter().count(),
        "active_polygons":result.evidence.provenance.original_polygon_ids.len(),"canceled_sheets":result.evidence.provenance.canceled_sheets.len(),
        "source_signed_volume":result.evidence.closed_source_volume(),
        "regions":result.native.regions.len(),
        "failed_regions":result.native.regions.iter().filter(|r|r.faces.is_err()).count(),
        "pending_chart_cuts":result.native.regions.iter().filter(|r|r.requires_chart_cut).count(),
        "ambiguous_source_edges":result.trims.ambiguous_edges.len()});
    if let Some(shell) = result.complete_shell() {
        let condition = shell.shell_condition();
        summary["native_faces"] = json!(shell.len());
        summary["shell_condition"] = json!(format!("{condition:?}"));
        summary["shell_components"] = json!(shell.connected_components().len());
        if condition == ShellCondition::Closed && args.verification == VerificationMode::Full {
            let coarse = numerical_volume(&shell, args.tolerance);
            let fine = numerical_volume(&shell, args.tolerance * 0.25);
            let converged = fine["closed_signed_volume"]
                .as_f64()
                .zip(coarse["closed_signed_volume"].as_f64())
                .and_then(|(a, b)| (a.abs() > 0.0).then(|| (a - b).abs() / a.abs()));
            let relative_source_difference = fine["closed_signed_volume"]
                .as_f64()
                .zip(result.evidence.closed_source_volume())
                .and_then(|(a, b)| (b.abs() > 0.0).then(|| (a - b).abs() / b.abs()));
            summary["volume"] = json!({"coarse":coarse,"fine":fine,
                "relative_resolution_difference":converged,"relative_source_difference":relative_source_difference});
        }
        // Do not call a disconnected, open or unmeasurable shell a solid.
        // The numerical checks are still not a global self-intersection proof.
        let solid_eligible = args.verification == VerificationMode::Full
            && condition == ShellCondition::Closed
            && shell.is_connected()
            && result
                .evidence
                .closed_source_volume()
                .is_some_and(|v| v > 0.0)
            && summary["volume"]["fine"]["closed_signed_volume"]
                .as_f64()
                .is_some_and(|v| v > 0.0);
        summary["step"] = match step_roundtrip(&shell, solid_eligible, args) {
            Ok(info) => info,
            Err(error) => json!({"error":error.to_string()}),
        };
    }
    Ok(summary)
}

fn verify_native(path: &Path, args: &Args) -> Result<Value> {
    let compressed: CompressedShell<Point3, Curve, Surface> =
        serde_json::from_reader(BufReader::new(File::open(path)?))
            .context("read compressed native shell")?;
    let shell = Shell::extract(compressed)?;
    let condition = shell.shell_condition();
    let components = shell.connected_components().len();
    let mut summary = json!({"source":path,"native_faces":shell.len(),"shell_condition":format!("{condition:?}"),
        "shell_components":components});
    if args.verification == VerificationMode::Full {
        let coarse = numerical_volume(&shell, args.tolerance);
        let fine = numerical_volume(&shell, args.tolerance * 0.25);
        summary["volume"] = json!({"coarse":coarse,"fine":fine});
    }
    // The native JSON does not certify material classification.
    summary["step"] = step_roundtrip(&shell, false, args)
        .unwrap_or_else(|error| json!({"error":format!("{error:#}")}));
    Ok(summary)
}

fn verify_step_input(path: &Path, args: &Args) -> Result<Value> {
    let step = std::fs::read_to_string(path)?;
    let table = Table::from_step(&step)?;
    let mut summary = inspect_step(
        &table,
        step.len(),
        None,
        args.tolerance,
        args.verification == VerificationMode::Full,
    )?;
    summary["source"] = json!(path);
    Ok(summary)
}

fn main() -> Result<()> {
    let args = Args::parse();
    anyhow::ensure!(
        args.tolerance.is_finite() && args.tolerance >= 4.0 * KERNEL_EPSILON,
        "invalid tessellation tolerance"
    );
    anyhow::ensure!(
        args.step_output.is_none() || args.inputs.len() == 1,
        "STEP output requires one input"
    );
    anyhow::ensure!(
        args.input != InputKind::Step
            || (args.step_output.is_none() && args.verification != VerificationMode::Export),
        "STEP input supports full or topology verification only"
    );
    anyhow::ensure!(
        args.verification != VerificationMode::Export || args.step_output.is_some(),
        "export mode requires STEP output"
    );
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    for path in &args.inputs {
        let report = (match args.input {
            InputKind::Step => verify_step_input(path, &args),
            InputKind::NativeShell => verify_native(path, &args),
            InputKind::Obj => verify(path, &args),
        })
        .unwrap_or_else(|error| json!({"source":path,"error":format!("{error:#}")}));
        writeln!(out, "{report}")?;
        out.flush()?;
    }
    Ok(())
}
