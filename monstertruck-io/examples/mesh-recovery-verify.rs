//! End-to-end native mesh recovery verification, independent of CAD labels.
//! Volume measurements are numerical diagnostics, not an intersection proof.
//! All STEP output is held in memory unless another caller explicitly saves it.

use anyhow::{Context, Result};
use clap::Parser;
use monstertruck_io::step::{
    load::Table,
    save::{CompleteStepDisplay, StepModel},
};
use monstertruck_mesh::{PolygonMesh, obj};
use monstertruck_meshing::{
    analyzers::{CalcVolume, Topology},
    filters::OptimizingFilter,
    tessellation::{
        MeshableShape, MeshedShape, face_drop_count, reset_face_drop_count, shell_to_polygon_strict,
    },
};
use monstertruck_modeling::{
    EuclideanSpace, KERNEL_EPSILON, Shell, Solid, ToCompressedTrimmedParameterCurves,
    mesh_recovery::{AnalyticRecoveryOptions, recover_analytic_mesh},
};
use monstertruck_topology::{shell::ShellCondition, trimmed::TrimmedShell};
use serde_json::{Value, json};
use std::{
    fs::File,
    io::{BufReader, Write},
    path::{Path, PathBuf},
};

/// Recover native B-rep regions, verify strict tessellation and STEP reload.
#[derive(Parser, Debug)]
struct Args {
    /// Maximum explicitly allowed planar repair displacement in source units.
    #[arg(long, default_value_t = 0.0)]
    max_planar_displacement: f64,
    /// Explicitly permit cancellation of proven reverse coincident sheets.
    #[arg(long)]
    normalize_coincident_sheets: bool,
    /// First tessellation tolerance in source units; second is four times finer.
    #[arg(long, default_value_t = 0.01)]
    tolerance: f64,
    /// OBJ inputs.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,
}

fn numerical_volume(shell: &Shell, tolerance: f64) -> Value {
    let _ = reset_face_drop_count();
    match shell_to_polygon_strict(shell, tolerance) {
        Err(error) => {
            json!({"tolerance":tolerance,"error":error.to_string(),"dropped_faces":face_drop_count()})
        }
        Ok(mesh) => {
            let dropped = face_drop_count();
            if dropped != 0 {
                return json!({"tolerance":tolerance,"dropped_faces":dropped,"error":"empty tessellation face"});
            }
            measure_polygon_volume(mesh, tolerance)
        }
    }
}

fn measure_polygon_volume(mut mesh: PolygonMesh, tolerance: f64) -> Value {
    if let Some(&origin) = mesh.positions().first() {
        for point in mesh.positions_mut() {
            *point -= origin.to_vec();
        }
    }
    // Only for the tessellated diagnostic mesh, not the native B-rep.
    // This joins independently emitted face vertices at numerical resolution.
    mesh.put_together_same_attrs(KERNEL_EPSILON);
    let condition = mesh.shell_condition();
    let signed_integral = mesh.volume();
    json!({"tolerance":tolerance,"dropped_faces":0,
        "polygons":mesh.face_iter().count(),"mesh_condition":format!("{condition:?}"),
        "signed_integral":signed_integral,
        "closed_signed_volume":(condition==ShellCondition::Closed).then_some(signed_integral)})
}

fn step_roundtrip(shell: &Shell, solid_eligible: bool, tolerance: f64) -> Result<Value> {
    let step = if solid_eligible {
        let solid = Solid::new(vec![shell.clone()])?;
        let compressed = solid.compress_with_parameter_curves(KERNEL_EPSILON * 100.0);
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string()
    } else {
        let compressed = shell.compress_with_parameter_curves(KERNEL_EPSILON * 100.0);
        CompleteStepDisplay::new(StepModel::from(&compressed), Default::default()).to_string()
    };
    let table = Table::from_step(&step)?;
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
        if shell.shell_condition() == ShellCondition::Closed {
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
                restored_volumes.push(measure_polygon_volume(meshed.to_polygon(), tolerance));
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
        json!({"bytes":step.len(),"faces":restored_faces,"edges":restored_edges,
        "face_count_preserved":restored_faces==shell.len(),"solid_breps":table.manifold_solid_brep.len(),
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
            normalize_coincident_sheets: args.normalize_coincident_sheets,
            ..Default::default()
        },
    )?;
    let mut summary = json!({"source":path,"source_polygons":mesh.face_iter().count(),
        "active_polygons":result.original_polygon_ids.len(),"canceled_sheets":result.canceled_sheets.len(),
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
        if condition == ShellCondition::Closed {
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
        let solid_eligible = condition == ShellCondition::Closed
            && shell.is_connected()
            && result
                .evidence
                .closed_source_volume()
                .is_some_and(|v| v > 0.0)
            && summary["volume"]["fine"]["closed_signed_volume"]
                .as_f64()
                .is_some_and(|v| v > 0.0);
        summary["step"] = match step_roundtrip(&shell, solid_eligible, args.tolerance * 0.25) {
            Ok(info) => info,
            Err(error) => json!({"error":error.to_string()}),
        };
    }
    Ok(summary)
}

fn main() -> Result<()> {
    let args = Args::parse();
    anyhow::ensure!(
        args.tolerance.is_finite() && args.tolerance >= 4.0 * KERNEL_EPSILON,
        "invalid tessellation tolerance"
    );
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    for path in &args.inputs {
        let report = verify(path, &args)
            .unwrap_or_else(|error| json!({"source":path,"error":format!("{error:#}")}));
        writeln!(out, "{report}")?;
        out.flush()?;
    }
    Ok(())
}
