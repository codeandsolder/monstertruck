#[cfg(feature = "cad-kernel-monstertruck")]
mod enabled {
    use anyhow::{Context, Result, bail};
    use clap::Parser;
    use std::{fs::OpenOptions, io::Write, path::PathBuf};
    use step_redox::{
        cad_kernel::{CadKernel, monstertruck::MonstertruckKernel},
        cad_model::NodeId,
    };

    #[derive(Parser)]
    #[command(
        name = "step-cad-rebuild",
        about = "Reconstruct a source solid or inspect selected intermediate constructive nodes"
    )]
    struct Cli {
        input: PathBuf,
        output: PathBuf,
        /// Emit these candidate DAG node roots as separate STEP solids.
        /// When specified, OUTPUT must be an existing audit directory.
        #[arg(long, value_delimiter = ',', num_args = 1..)]
        nodes: Vec<usize>,
    }

    pub fn run() -> Result<()> {
        let cli = Cli::parse();
        let bytes =
            std::fs::read(&cli.input).with_context(|| format!("read {}", cli.input.display()))?;
        let kernel = MonstertruckKernel::project_standard();
        let report =
            step_redox::complete_ir::recover_complete_ir_with_kernel_bytes(&bytes, &kernel)?;
        if report.body_count != 1 || report.constructively_recovered_solid_count != 1 {
            bail!(
                "expected exactly one constructively recovered body, got bodies={} constructive={}",
                report.body_count,
                report.constructively_recovered_solid_count
            );
        }
        let fragments = report.selected_fragments().collect::<Vec<_>>();
        let [selected] = fragments.as_slice() else {
            bail!(
                "expected exactly one selected constructive fragment, got {}",
                fragments.len()
            );
        };
        let fragment = &selected.fragment;
        let audit_nodes = !cli.nodes.is_empty();
        if audit_nodes && !cli.output.is_dir() {
            bail!("--nodes requires OUTPUT to be an existing directory");
        }
        let nodes = if audit_nodes {
            cli.nodes
        } else {
            vec![fragment.root.0]
        };
        for id in nodes {
            let root = NodeId(id);
            fragment
                .model
                .node(root)
                .with_context(|| format!("invalid CAD node #{id}"))?;
            let evaluated = kernel.evaluate(&fragment.model, root)?;
            let summary = kernel.summarize(&evaluated);
            if !summary.geometrically_consistent || summary.shells != 1 {
                bail!(
                    "Monstertruck did not rebuild one geometrically consistent solid at node #{id}: {summary:?}"
                );
            }
            let output = if audit_nodes {
                cli.output.join(format!("node-{id}.step"))
            } else {
                cli.output.clone()
            };
            // Existing audit snapshots and the author's STEP are never overwritten.
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&output)
                .with_context(|| format!("create rebuilt STEP {}", output.display()))?;
            file.write_all(kernel.to_step(&evaluated)?.as_bytes())?;
            eprintln!(
                "reconstructed node #{id} => {}: {summary:?}",
                output.display()
            );
        }
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        #[test]
        fn parses_single_output_and_selected_nodes() -> Result<()> {
            let single = Cli::try_parse_from(["step-cad-rebuild", "input.step", "rebuilt.step"])?;
            assert!(single.nodes.is_empty());
            let batch = Cli::try_parse_from([
                "step-cad-rebuild",
                "input.step",
                "audit-dir",
                "--nodes",
                "84,104,109,116",
            ])?;
            assert_eq!(batch.nodes, [84, 104, 109, 116]);
            Ok(())
        }
    }
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> anyhow::Result<()> { enabled::run() }

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> anyhow::Result<()> {
    anyhow::bail!("step-cad-rebuild requires --features cad-kernel-monstertruck")
}
