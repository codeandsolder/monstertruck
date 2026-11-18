//! Read-only STEP 3-D edge locus separation, using an incident source plane.
#[cfg(feature = "cad-kernel-monstertruck")]
use anyhow::Context;
use anyhow::{Result, bail};
#[cfg(feature = "cad-kernel-monstertruck")]
use clap::Parser;
#[cfg(feature = "cad-kernel-monstertruck")]
use std::path::PathBuf;

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Parser)]
#[command(
    name = "step-edge-locus-compare",
    about = "Bound a source-to-target full B-spline edge separation via an incident source plane"
)]
struct Cli {
    source: PathBuf,
    source_edge_id: u64,
    source_face_id: u64,
    target: PathBuf,
    target_edge_id: u64,
    #[arg(long, default_value_t = 1.0e-7)]
    endpoint_tolerance_mm: f64,
    #[arg(long, default_value_t = 1.0e-10)]
    target_gap_mm: f64,
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn read_step(path: &PathBuf) -> Result<step_p21::ast::Exchange> {
    let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let source = std::str::from_utf8(&bytes).context("STEP file is not UTF-8")?;
    step_p21::parser::parse(source).context("parse STEP file")
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> Result<()> {
    use step_redox::source_edge_locus::{
        StepEdge, bound_complete_step_edges_separation_by_plane, incident_edge_plane,
    };

    let args = Cli::parse();
    let source = read_step(&args.source)?;
    let target = read_step(&args.target)?;
    let [source_section] = source.data.as_slice() else {
        bail!("source STEP requires exactly one DATA section");
    };
    let [target_section] = target.data.as_slice() else {
        bail!("target STEP requires exactly one DATA section");
    };

    let plane = incident_edge_plane(
        args.source_face_id,
        args.source_edge_id,
        &source_section.entities,
    )?;
    let result = bound_complete_step_edges_separation_by_plane(
        StepEdge {
            id: args.source_edge_id,
            entities: &source_section.entities,
        },
        StepEdge {
            id: args.target_edge_id,
            entities: &target_section.entities,
        },
        plane,
        args.endpoint_tolerance_mm,
        args.target_gap_mm,
    )?;

    serde_json::to_writer_pretty(
        std::io::stdout().lock(),
        &serde_json::json!({
            "source": {"path": args.source, "edge_id": args.source_edge_id,
                "plane_face_id": args.source_face_id},
            "target": {"path": args.target, "edge_id": args.target_edge_id},
            "reference_plane": {"origin_mm": plane.origin_mm, "normal": plane.normal},
            "source_plane_max_deviation_mm": {
                "lower": result.source.lower, "upper": result.source.upper,
                "arithmetic_margin": result.source.arithmetic_margin,
                "subdivisions": result.source.subdivisions
            },
            "target_plane_max_deviation_mm": {
                "lower": result.target.lower, "upper": result.target.upper,
                "arithmetic_margin": result.target.arithmetic_margin,
                "subdivisions": result.target.subdivisions
            },
            "minimum_source_to_target_directed_distance_mm":
                result.minimum_source_to_target_mm
        }),
    )?;
    println!();
    Ok(())
}

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> Result<()> {
    bail!("step-edge-locus-compare requires --features cad-kernel-monstertruck")
}
