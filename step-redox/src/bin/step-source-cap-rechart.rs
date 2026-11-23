//! Write a detached, provenance-bearing STEP proof view with fully bounded
//! spline-cap curve proposals, never overwriting the authoring source.

#[cfg(feature = "cad-kernel-monstertruck")]
mod enabled {
    use anyhow::{Context, Result, ensure};
    use clap::Parser;
    use std::{
        fs::{self, OpenOptions},
        io::Write,
        path::{Path, PathBuf},
    };

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
        let step = step_redox::write_exchange(&proposal.exchange)?;
        let mut manifest = serde_json::to_vec_pretty(&proposal.manifest)?;
        manifest.push(b'\n');
        write_pair_new(
            &cli.proof_step,
            &cli.manifest_json,
            step.as_bytes(),
            &manifest,
        )?;
        Ok(())
    }

    /// Create both outputs without replacing any existing file. If either
    /// creation or write fails, remove only files created in this invocation
    /// rather than leaving an orphaned proof without its provenance manifest.
    fn write_pair_new(proof: &Path, manifest: &Path, step: &[u8], metadata: &[u8]) -> Result<()> {
        let mut proof_file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(proof)
            .with_context(|| format!("create detached STEP {}", proof.display()))?;
        let mut manifest_file = match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(manifest)
        {
            Ok(file) => file,
            Err(error) => {
                drop(proof_file);
                let _ = fs::remove_file(proof);
                return Err(error)
                    .with_context(|| format!("create provenance manifest {}", manifest.display()));
            }
        };
        let result = (|| -> Result<()> {
            proof_file.write_all(step)?;
            manifest_file.write_all(metadata)?;
            proof_file.flush()?;
            manifest_file.flush()?;
            Ok(())
        })();
        if result.is_err() {
            drop(proof_file);
            drop(manifest_file);
            let _ = fs::remove_file(proof);
            let _ = fs::remove_file(manifest);
        }
        result
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn output_pair_is_create_only_and_does_not_leave_orphans() -> Result<()> {
            let temporary = tempfile::tempdir()?;
            let proof = temporary.path().join("view.step");
            let manifest = temporary.path().join("manifest.json");
            fs::write(&manifest, b"existing manifest")?;
            assert!(write_pair_new(&proof, &manifest, b"step", b"manifest").is_err());
            assert!(!proof.exists());
            assert_eq!(fs::read(&manifest)?, b"existing manifest");
            fs::remove_file(&manifest)?;
            write_pair_new(&proof, &manifest, b"step", b"manifest")?;
            assert_eq!(fs::read(&proof)?, b"step");
            assert_eq!(fs::read(&manifest)?, b"manifest");
            assert!(write_pair_new(&proof, &manifest, b"changed", b"changed").is_err());
            assert_eq!(fs::read(&proof)?, b"step");
            Ok(())
        }

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
