//! Opt-in source STEP healer. Always writes a *new* destination, and prints
//! a machine-readable provenance manifest to stdout.

#[cfg(feature = "cad-kernel-monstertruck")]
use anyhow::{Context, Result, anyhow, bail};
#[cfg(feature = "cad-kernel-monstertruck")]
use clap::Parser;
#[cfg(feature = "cad-kernel-monstertruck")]
use std::{io::Write, path::PathBuf};

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Parser)]
#[command(
    name = "step-heal-source-hyperbolas",
    about = "Create a separate fully proven healed STEP document, keeping original source intact"
)]
struct Cli {
    input: PathBuf,
    /// Must not exist. This tool never overwrites an existing source or output.
    output: PathBuf,
    /// Repeat `DATA_SECTION:SOLID_ID:EDGE_ID:PLANE_FACE_ID:CONE_FACE_ID`
    #[arg(long, required = true)]
    repair: Vec<String>,
    /// Maximum original STEP vertex to new analytic edge error in mm.
    #[arg(long)]
    endpoint_tolerance_mm: f64,
    /// Desired global original spline-carrier proof gap in mm.
    #[arg(long)]
    target_gap_mm: f64,
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn parse_repair(
    specification: &str,
    endpoint_tolerance_mm: f64,
    target_gap_mm: f64,
) -> Result<step_redox::full_source_hyperbola_healing::FullSourceHyperbolaRepair> {
    use step_redox::full_source_hyperbola_healing::FullSourceHyperbolaRepair;
    let parts = specification.split(':').collect::<Vec<_>>();
    let [section, solid, edge, plane, cone] = parts.as_slice() else {
        bail!("repair must use DATA_SECTION:SOLID:EDGE:PLANE_FACE:CONE_FACE");
    };
    let parse = |s: &str| -> Result<u64> {
        s.parse::<u64>()
            .with_context(|| format!("invalid positive integer source entity ID {s}"))
            .and_then(|id| {
                if id == 0 {
                    bail!("source STEP entity ID cannot be zero");
                }
                Ok(id)
            })
    };
    let data_section = section
        .parse::<usize>()
        .map_err(|_| anyhow!("invalid DATA section ordinal {section}"))?;
    Ok(FullSourceHyperbolaRepair {
        data_section,
        solid_id: parse(solid)?,
        edge_id: parse(edge)?,
        plane_face_id: parse(plane)?,
        cone_face_id: parse(cone)?,
        endpoint_tolerance_mm,
        target_gap_mm,
    })
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> Result<()> {
    use step_redox::full_source_hyperbola_healing::heal_step_document;
    let cli = Cli::parse();
    let mut specs = Vec::with_capacity(cli.repair.len());
    for repair in &cli.repair {
        specs.push(parse_repair(
            repair,
            cli.endpoint_tolerance_mm,
            cli.target_gap_mm,
        )?);
    }
    let original = std::fs::read(&cli.input)
        .with_context(|| format!("read original STEP file {}", cli.input.display()))?;
    let revised = heal_step_document(&original, &specs)?;
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&cli.output)
        .with_context(|| {
            format!(
                "create new, never-overwritten STEP output {}",
                cli.output.display()
            )
        })?;
    file.write_all(revised.step.as_bytes())
        .with_context(|| format!("write new STEP output {}", cli.output.display()))?;
    file.sync_all()
        .with_context(|| format!("sync new STEP output {}", cli.output.display()))?;
    // Data is printed only after the file is durably written.
    serde_json::to_writer_pretty(std::io::stdout().lock(), &revised.repairs)?;
    println!();
    Ok(())
}

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> anyhow::Result<()> {
    anyhow::bail!("step-heal-source-hyperbolas requires --features cad-kernel-monstertruck")
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;

    #[test]
    fn explicit_source_section_and_all_original_entities_are_required() -> Result<()> {
        assert!(parse_repair("0:15385:93481:125603", 1e-8, 1e-9).is_err());
        assert!(parse_repair("0:0:93481:125603:58899", 1e-8, 1e-9).is_err());
        assert!(parse_repair("wrong:15385:93481:125603:58899", 1e-8, 1e-9).is_err());
        let r = parse_repair("0:15385:93481:125603:58899", 1e-8, 1e-9)?;
        assert_eq!(r.data_section, 0);
        assert_eq!(r.solid_id, 15385);
        assert_eq!(r.edge_id, 93481);
        assert_eq!(r.plane_face_id, 125603);
        assert_eq!(r.cone_face_id, 58899);
        Ok(())
    }

    #[test]
    fn cli_requires_new_destination_budgets_and_repair_specification() -> Result<()> {
        assert!(
            Cli::try_parse_from(["step-heal-source-hyperbolas", "old.step", "new.step",]).is_err()
        );
        let cli = Cli::try_parse_from([
            "step-heal-source-hyperbolas",
            "old.step",
            "new.step",
            "--repair",
            "0:15385:93481:125603:58899",
            "--endpoint-tolerance-mm",
            "0.00000001",
            "--target-gap-mm",
            "0.000000001",
        ])?;
        assert_eq!(cli.repair.len(), 1);
        Ok(())
    }
}
