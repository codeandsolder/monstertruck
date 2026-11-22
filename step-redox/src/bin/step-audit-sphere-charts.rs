//! Write a distinct, strictly boundary-verified source for an import audit.

use anyhow::{Context, Result};
use clap::Parser;
use std::{io::Write, path::PathBuf};
use step_redox::boundary_proof::audit_sphere_charts;

#[derive(Parser)]
#[command(about = "Normalize eligible spherical charts without changing source geometry")]
struct Cli {
    input: PathBuf,
    solid_id: u64,
    /// A new output file; never overwrites an existing source or result.
    output: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let input =
        std::fs::read(&cli.input).with_context(|| format!("read {}", cli.input.display()))?;
    let audit = audit_sphere_charts(&input, cli.solid_id)?;
    let mut output = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&cli.output)
        .with_context(|| format!("create {}", cli.output.display()))?;
    output.write_all(audit.step.as_bytes())?;
    output.sync_all()?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), &audit.manifest)?;
    println!();
    Ok(())
}
