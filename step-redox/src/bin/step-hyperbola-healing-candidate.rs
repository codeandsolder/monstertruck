//! Report a separately provenance-tracked analytic replacement for a STEP edge.
//!
//! Read-only: emits JSON evidence and a STEP geometry *fragment*, never
//! silently mutates the original source or claims B-rep equivalence.

#[cfg(feature = "cad-kernel-monstertruck")]
use anyhow::Context;
use anyhow::{Result, bail};
#[cfg(feature = "cad-kernel-monstertruck")]
use clap::Parser;
#[cfg(feature = "cad-kernel-monstertruck")]
use serde::Serialize;
#[cfg(feature = "cad-kernel-monstertruck")]
use std::path::PathBuf;

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Parser)]
#[command(
    name = "step-hyperbola-healing-candidate",
    about = "Read-only proof and native analytic hyperbola proposal for a source STEP edge"
)]
struct Cli {
    input: PathBuf,
    solid_id: u64,
    edge_id: u64,
    plane_face_id: u64,
    cone_face_id: u64,
    /// Millimetres per source STEP coordinate, resolved from the input context.
    #[arg(long)]
    coordinate_scale_mm: f64,
    /// Maximum physical mismatch of each original vertex and healed hyperbola.
    #[arg(long)]
    endpoint_tolerance_mm: f64,
    /// Desired global Bernstein bracket refinement gap.
    #[arg(long)]
    target_gap_mm: f64,
}

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Serialize)]
struct CandidateReport<'a> {
    evidence: &'a step_redox::source_edge_healing::HyperbolaHealingEvidence,
    oriented_hyperbola_parameters: [f64; 2],
    /// A syntactically valid STEP entity fragment, NOT a complete STEP file.
    native_hyperbola_step_fragment: String,
    /// Analytic intersection with corrected plane and original-cone support;
    /// this is a STEP geometry fragment, not a healed topological solid.
    native_intersection_step_fragment: String,
    /// New open STEP edge with original directed vertices and two source coedges.
    /// This is not a complete healed face, shell or solid.
    healed_open_edge: step_redox::healed_edge_topology::HealedEdgeTopologyFragment,
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> Result<()> {
    use monstertruck_io::step::{load::ParsedStepSource, save::StepDisplay};
    use step_redox::{
        native_hyperbola_intersection::native_healed_intersection_curve,
        source_edge_healing::{HyperbolaHealingRequest, propose_plane_cone_hyperbola_healing},
        source_edge_locus::StepEdge,
    };

    let cli = Cli::parse();
    let bytes =
        std::fs::read(&cli.input).with_context(|| format!("read {}", cli.input.display()))?;
    let parsed = ParsedStepSource::parse(&bytes)?;
    let [section] = parsed.exchange.data.as_slice() else {
        bail!("the input STEP must contain exactly one DATA section");
    };
    let candidate = propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
        solid_id: cli.solid_id,
        source: StepEdge {
            id: cli.edge_id,
            entities: &section.entities,
        },
        source_plane_face_id: cli.plane_face_id,
        source_cone_face_id: cli.cone_face_id,
        source_coordinate_scale_mm: cli.coordinate_scale_mm,
        endpoint_tolerance_mm: cli.endpoint_tolerance_mm,
        target_gap_mm: cli.target_gap_mm,
    })?;
    let healed_open_edge = step_redox::healed_edge_topology::healed_edge_with_source_vertices(
        &candidate,
        &section.entities,
        cli.endpoint_tolerance_mm,
    )?;
    let native = candidate.modeling_curve();
    let report = CandidateReport {
        evidence: &candidate.evidence,
        oriented_hyperbola_parameters: [
            candidate.trim.start_parameter,
            candidate.trim.end_parameter,
        ],
        native_hyperbola_step_fragment: StepDisplay::new(&native, 1).to_string(),
        native_intersection_step_fragment: StepDisplay::new(
            &native_healed_intersection_curve(&candidate),
            1,
        )
        .to_string(),
        healed_open_edge,
    };
    serde_json::to_writer_pretty(std::io::stdout().lock(), &report)?;
    println!();
    Ok(())
}

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> Result<()> {
    bail!("step-hyperbola-healing-candidate requires --features cad-kernel-monstertruck")
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;

    #[test]
    fn the_source_scale_and_both_error_budgets_are_explicit() -> Result<()> {
        assert!(
            Cli::try_parse_from([
                "step-hyperbola-healing-candidate",
                "part.step",
                "1",
                "2",
                "3",
                "4",
            ])
            .is_err()
        );
        let cli = Cli::try_parse_from([
            "step-hyperbola-healing-candidate",
            "part.step",
            "1",
            "2",
            "3",
            "4",
            "--coordinate-scale-mm",
            "1",
            "--endpoint-tolerance-mm",
            "0.00000001",
            "--target-gap-mm",
            "0.000000001",
        ])?;
        assert_eq!(cli.solid_id, 1);
        assert_eq!(cli.edge_id, 2);
        assert_eq!(cli.plane_face_id, 3);
        assert_eq!(cli.cone_face_id, 4);
        Ok(())
    }
}
