//! Check source STEP B-spline loci against their own incident planar carriers.
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
    name = "step-plane-edge-audit",
    about = "Audit full source B-spline loci against their own literal STEP PLANE carriers"
)]
struct Cli {
    input: PathBuf,
    solid_id: u64,
    /// Absolute edge-to-face incidence tolerance, in STEP millimeters.
    #[arg(long)]
    tolerance_mm: f64,
    /// Requested geometric gap between lower and upper deviation bounds.
    #[arg(long, default_value_t = 1.0e-10)]
    target_gap_mm: f64,
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> Result<()> {
    let cli = Cli::parse();
    let bytes =
        std::fs::read(&cli.input).with_context(|| format!("read {}", cli.input.display()))?;
    let source = monstertruck_io::step::load::ParsedStepSource::parse(&bytes)?.exchange;
    let [section] = source.data.as_slice() else {
        bail!("source STEP requires exactly one DATA section");
    };
    let coordinate_scale_mm = step_redox::source_edge_locus::source_solid_length_scale_mm(
        cli.solid_id,
        &section.entities,
    )?;
    let audit = step_redox::source_edge_locus::audit_source_solid_planar_bspline_incidences(
        cli.solid_id,
        &section.entities,
        cli.tolerance_mm,
        cli.target_gap_mm,
        coordinate_scale_mm,
    )?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), &audit)?;
    println!();
    Ok(())
}

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> Result<()> {
    bail!("step-plane-edge-audit requires --features cad-kernel-monstertruck")
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;

    #[test]
    fn caller_must_choose_an_explicit_source_tolerance() -> Result<()> {
        assert!(
            Cli::try_parse_from(["step-plane-edge-audit", "fixture.step", "15"]).is_err(),
            "undocumented default units would silently misclassify source geometry"
        );
        let cli = Cli::try_parse_from([
            "step-plane-edge-audit",
            "fixture.step",
            "15",
            "--tolerance-mm",
            "0.00001",
        ])?;
        assert_eq!(cli.solid_id, 15);
        assert_eq!(cli.tolerance_mm, 1.0e-5);
        Ok(())
    }
}
