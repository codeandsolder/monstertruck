//! Write a detached, provenance-bearing STEP proof view with fully bounded
//! spline-cap curve proposals, never overwriting the authoring source.

#[cfg(feature = "cad-kernel-monstertruck")]
mod enabled {
    use anyhow::{Context, Result, ensure};
    use clap::Parser;
    use std::{fs::OpenOptions, io::Write, path::PathBuf};

    #[derive(Parser)]
    #[command(
        name = "step-source-cap-rechart",
        about = "Rechart source spline caps in a detached STEP proof view; preserve all original geometry and provenance"
    )]
    struct Cli {
        source_step: PathBuf,
        source_solid_id: u64,
        proof_step: PathBuf,
        manifest_json: PathBuf,
        #[arg(long, default_value_t = 0)]
        data_section: usize,
        /// Reject any result lacking a complete, independently certified closed boundary.
        #[arg(long)]
        require_complete_proof: bool,
    }

    pub fn run() -> Result<()> {
        let cli = Cli::parse();
        let original = std::fs::read(&cli.source_step)
            .with_context(|| format!("read original STEP {}", cli.source_step.display()))?;
        let proposal = step_redox::propose_source_cap_rechart(
            &original,
            cli.data_section,
            cli.source_solid_id,
        )?;
        ensure!(
            !cli.require_complete_proof || proposal.manifest.proof_view_boundary_supported,
            "detached source boundary proof is incomplete: {:?}",
            proposal.manifest.proof_view_boundary_reason
        );
        // create_new=true makes source overwrite impossible; no environmental
        // flags can replace the original or overwrite an existing audit.
        let mut proof_file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&cli.proof_step)
            .with_context(|| {
                format!(
                    "create private STEP proof view {}",
                    cli.proof_step.display()
                )
            })?;
        let mut manifest_file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&cli.manifest_json)
            .with_context(|| {
                format!(
                    "create private patch manifest {}",
                    cli.manifest_json.display()
                )
            })?;
        let step = step_redox::write_exchange(&proposal.exchange)?;
        proof_file.write_all(step.as_bytes())?;
        serde_json::to_writer_pretty(&mut manifest_file, &proposal.manifest)?;
        manifest_file.write_all(b"\n")?;
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn requires_both_detached_outputs_and_optional_full_boundary_proof() -> Result<()> {
            assert!(
                Cli::try_parse_from([
                    "step-source-cap-rechart",
                    "original.step",
                    "3877",
                    "proof.step"
                ])
                .is_err()
            );
            let parsed = Cli::try_parse_from([
                "step-source-cap-rechart",
                "original.step",
                "3877",
                "proof.step",
                "manifest.json",
                "--data-section",
                "1",
                "--require-complete-proof",
            ])?;
            assert_eq!(parsed.source_solid_id, 3877);
            assert_eq!(parsed.data_section, 1);
            assert!(parsed.require_complete_proof);
            Ok(())
        }
    }
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> anyhow::Result<()> { enabled::run() }

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> anyhow::Result<()> {
    anyhow::bail!("step-source-cap-rechart requires --features cad-kernel-monstertruck")
}
