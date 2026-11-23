//! Inspect source-authored spline caps without changing their topology or geometry.
use anyhow::{Context, Result, ensure};
use clap::Parser;
use std::path::PathBuf;
use step_redox::boundary_proof::audit_source_spline_cap_precision;

#[derive(Parser)]
#[command(
    name = "step-audit-spline-caps",
    about = "Report complete source precision bounds on both caps of each STEP spline strip"
)]
struct Cli {
    input: PathBuf,
    solid_id: u64,
    #[arg(long, default_value_t = 0)]
    data_section: usize,
    /// Fail if at least one source spline-strip face lacks both certified caps.
    #[arg(long)]
    require_complete: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let bytes = std::fs::read(&cli.input)
        .with_context(|| format!("read source STEP {}", cli.input.display()))?;
    let audit = audit_source_spline_cap_precision(&bytes, cli.data_section, cli.solid_id)?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), &audit)?;
    println!();
    if cli.require_complete {
        ensure!(
            audit.strip_faces == audit.fully_certified_strip_faces,
            "{} / {} source spline-strip faces have two precision-certified caps",
            audit.fully_certified_strip_faces,
            audit.strip_faces
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_solid_and_source_path_are_required() -> Result<()> {
        assert!(Cli::try_parse_from(["step-audit-spline-caps", "input.step"]).is_err());
        let parsed = Cli::try_parse_from([
            "step-audit-spline-caps",
            "input.step",
            "100",
            "--data-section",
            "1",
            "--require-complete",
        ])?;
        assert_eq!(parsed.solid_id, 100);
        assert_eq!(parsed.data_section, 1);
        assert!(parsed.require_complete);
        Ok(())
    }
}
