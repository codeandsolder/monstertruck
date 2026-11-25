//! Opt-in source STEP healer. Always writes a new destination and emits
//! machine-readable, source-distinct provenance.

#[cfg(feature = "cad-kernel-monstertruck")]
use anyhow::{Context, Result, bail};
#[cfg(feature = "cad-kernel-monstertruck")]
use clap::{Parser, Subcommand};
#[cfg(feature = "cad-kernel-monstertruck")]
use std::{
    io::Write,
    path::{Path, PathBuf},
};

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Parser)]
#[command(
    name = "step-heal-source-edges",
    about = "Construct a separate provenance-tracked healed STEP document"
)]
struct Cli {
    input: PathBuf,
    /// Must not exist. The original source is never overwritten.
    output: PathBuf,
    /// Maximum source-vertex to analytical edge error in mm.
    #[arg(long, global = true)]
    endpoint_tolerance_mm: f64,
    /// Desired global source spline-carrier proof gap in mm.
    #[arg(long, global = true)]
    target_gap_mm: f64,
    #[command(subcommand)]
    healing: Healing,
}

#[cfg(feature = "cad-kernel-monstertruck")]
#[derive(Subcommand)]
enum Healing {
    /// Exact hyperbola from an axis-parallel plane and a cone.
    Hyperbolas {
        /// Repeat DATA_SECTION:SOLID_ID:EDGE_ID:PLANE_FACE_ID:CONE_FACE_ID.
        #[arg(long, required = true)]
        repair: Vec<String>,
    },
    /// Regular analytical branch of two nonparallel cylinders.
    Cylinders {
        /// Repeat DATA_SECTION:SOLID_ID:EDGE_ID.
        #[arg(long, required = true)]
        repair: Vec<String>,
        /// Maximum whole-locus error of the serialized STEP leader in mm.
        #[arg(long)]
        leader_tolerance_mm: f64,
        /// Maximum number of cubic segments per replacement edge.
        #[arg(long)]
        max_leader_segments: u32,
    },
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn parse_repair<const N: usize>(specification: &str) -> Result<(usize, [u64; N])> {
    let (section, ids) = specification
        .split_once(':')
        .context("repair must begin with DATA_SECTION followed by positive STEP IDs")?;
    let data_section = section
        .parse::<usize>()
        .context("invalid DATA section ordinal")?;
    let ids = ids
        .split(':')
        .map(|s| {
            let id = s
                .parse::<u64>()
                .with_context(|| format!("invalid positive source entity ID {s}"))?;
            if id == 0 {
                bail!("source STEP entity ID cannot be zero");
            }
            Ok(id)
        })
        .collect::<Result<Vec<_>>>()?;
    let ids = ids
        .try_into()
        .map_err(|_: Vec<u64>| anyhow::anyhow!("repair requires {N} source entity IDs"))?;
    Ok((data_section, ids))
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn write_new<M: serde::Serialize>(output: &Path, step: &str, manifest: &M) -> Result<()> {
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)
        .with_context(|| {
            format!(
                "create new, never-overwritten STEP output {}",
                output.display()
            )
        })?;
    file.write_all(step.as_bytes())?;
    file.sync_all()?;
    serde_json::to_writer_pretty(std::io::stdout().lock(), manifest)?;
    println!();
    Ok(())
}

#[cfg(feature = "cad-kernel-monstertruck")]
fn main() -> Result<()> {
    use step_redox::{
        full_source_cylinder_healing::{FullSourceCylinderRepair, heal_cylinder_step_document},
        full_source_hyperbola_healing::{FullSourceHyperbolaRepair, heal_step_document},
    };
    let cli = Cli::parse();
    let original = std::fs::read(&cli.input)
        .with_context(|| format!("read original STEP {}", cli.input.display()))?;
    let endpoint_tolerance_mm = cli.endpoint_tolerance_mm;
    let target_gap_mm = cli.target_gap_mm;
    match cli.healing {
        Healing::Hyperbolas { repair } => {
            let specs = repair
                .iter()
                .map(|s| {
                    let (data_section, [solid_id, edge_id, plane_face_id, cone_face_id]) =
                        parse_repair(s)?;
                    Ok(FullSourceHyperbolaRepair {
                        data_section,
                        solid_id,
                        edge_id,
                        plane_face_id,
                        cone_face_id,
                        endpoint_tolerance_mm,
                        target_gap_mm,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            let revised = heal_step_document(&original, &specs)?;
            write_new(&cli.output, &revised.step, &revised.repairs)
        }
        Healing::Cylinders {
            repair,
            leader_tolerance_mm,
            max_leader_segments,
        } => {
            let specs = repair
                .iter()
                .map(|s| {
                    let (data_section, [solid_id, edge_id]) = parse_repair(s)?;
                    Ok(FullSourceCylinderRepair {
                        data_section,
                        solid_id,
                        edge_id,
                        endpoint_tolerance_mm,
                        target_gap_mm,
                        leader_tolerance_mm,
                        max_leader_segments,
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            let revised = heal_cylinder_step_document(&original, &specs)?;
            write_new(&cli.output, &revised.step, &revised.repairs)
        }
    }
}

#[cfg(not(feature = "cad-kernel-monstertruck"))]
fn main() -> anyhow::Result<()> {
    anyhow::bail!("step-heal-source-edges requires --features cad-kernel-monstertruck")
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;

    #[test]
    fn explicit_source_section_and_all_original_entities_are_required() -> Result<()> {
        assert!(parse_repair::<4>("0:15385:93481:125603").is_err());
        assert!(parse_repair::<4>("0:0:93481:125603:58899").is_err());
        assert!(parse_repair::<4>("wrong:15385:93481:125603:58899").is_err());
        assert_eq!(
            parse_repair::<4>("0:15385:93481:125603:58899")?,
            (0, [15385, 93481, 125603, 58899])
        );
        assert_eq!(parse_repair::<2>("0:171:1049")?, (0, [171, 1049]));
        Ok(())
    }

    #[test]
    fn cli_requires_new_destination_budgets_and_repair_specification() -> Result<()> {
        assert!(Cli::try_parse_from(["step-heal-source-edges", "old.step", "new.step"]).is_err());
        let cli = Cli::try_parse_from([
            "step-heal-source-edges",
            "old.step",
            "new.step",
            "--endpoint-tolerance-mm",
            "0.00000001",
            "--target-gap-mm",
            "0.000000001",
            "cylinders",
            "--repair",
            "0:171:1049",
            "--leader-tolerance-mm",
            "0.0000000002",
            "--max-leader-segments",
            "8192",
        ])?;
        assert!(matches!(cli.healing,Healing::Cylinders{repair,..} if repair.len()==1));
        Ok(())
    }
}
