//! One typed, read-only STEP analysis CLI. Each command emits a single JSON document.
use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use serde::Serialize;
use std::io::Write;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    name = "step-inspect",
    about = "Analyze STEP solids and recovery evidence as JSON"
)]
struct Cli {
    #[command(subcommand)]
    command: Analysis,
}

#[derive(Debug, Subcommand)]
enum Analysis {
    Support {
        input: PathBuf,
    },
    Revolution {
        input: PathBuf,
    },
    RadialSlot {
        input: PathBuf,
    },
    Extrusion {
        input: PathBuf,
    },
    ExtrusionCore {
        input: PathBuf,
    },
    DominantExtrusion {
        input: PathBuf,
    },
    EdgeTreatment {
        input: PathBuf,
    },
    BoundaryFeature {
        input: PathBuf,
    },
    Sweep {
        input: PathBuf,
    },
    HealedExtrusionCore {
        input: PathBuf,
        solid_id: u64,
        cap_a: u64,
        cap_b: u64,
    },
}

impl Analysis {
    fn input(&self) -> &PathBuf {
        match self {
            Self::Support { input }
            | Self::Revolution { input }
            | Self::RadialSlot { input }
            | Self::Extrusion { input }
            | Self::ExtrusionCore { input }
            | Self::DominantExtrusion { input }
            | Self::EdgeTreatment { input }
            | Self::BoundaryFeature { input }
            | Self::Sweep { input }
            | Self::HealedExtrusionCore { input, .. } => input,
        }
    }
}

#[derive(Serialize)]
struct SweepAnalysis {
    round_sweeps: Vec<step_redox::solid_sweeps::RecoveredRoundSweep>,
    revolved_round_tail: Vec<step_redox::solid_sweeps::RecoveredRevolvedRoundTail>,
    open_rectangular: Vec<step_redox::solid_sweeps::RecoveredOpenRectangularSweep>,
}

fn json<T: Serialize>(value: &T, output: &mut impl Write) -> Result<()> {
    serde_json::to_writer(output, value)?;
    Ok(())
}

fn inspect(command: &Analysis, data: &[u8], output: &mut impl Write) -> Result<()> {
    match command {
        Analysis::Support { .. } => json(
            &step_redox::detect_solid_surface_signatures_bytes(data)?,
            output,
        ),
        Analysis::Revolution { .. } => {
            json(&step_redox::detect_solid_revolutions_bytes(data)?, output)
        }
        Analysis::RadialSlot { .. } => json(
            &step_redox::detect_radial_slot_revolutions_bytes(data)?,
            output,
        ),
        Analysis::Extrusion { .. } => {
            json(&step_redox::detect_solid_extrusions_bytes(data)?, output)
        }
        Analysis::ExtrusionCore { .. } => json(
            &step_redox::detect_extrusion_core_candidates_bytes(data)?,
            output,
        ),
        Analysis::DominantExtrusion { .. } => json(
            &step_redox::detect_dominant_extrusion_candidates_bytes(data)?,
            output,
        ),
        Analysis::EdgeTreatment { .. } => json(
            &step_redox::detect_edge_treatment_evidence_bytes(data)?,
            output,
        ),
        Analysis::BoundaryFeature { .. } => {
            json(&step_redox::analyze_boundary_features_bytes(data)?, output)
        }
        Analysis::Sweep { .. } => json(
            &SweepAnalysis {
                round_sweeps: step_redox::detect_round_sweeps_bytes(data)?,
                revolved_round_tail: step_redox::detect_revolved_round_tails_bytes(data)?,
                open_rectangular: step_redox::detect_open_rectangular_sweeps_bytes(data)?,
            },
            output,
        ),
        Analysis::HealedExtrusionCore {
            solid_id,
            cap_a,
            cap_b,
            ..
        } => json(
            &step_redox::recover_treatment_healed_extrusion_core_from_cap_faces_bytes(
                data,
                *solid_id,
                [*cap_a, *cap_b],
            )?,
            output,
        ),
    }
}

fn main() -> Result<()> {
    let command = Cli::parse().command;
    let path = command.input();
    let data = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    inspect(&command, &data, &mut std::io::stdout().lock())
}

#[cfg(test)]
mod tests {
    use super::*;

    const TORUS: &[u8] = include_bytes!("../../validation/fixtures/native_torus.step");

    #[test]
    fn scanner_subcommands_select_the_expected_analysis() -> Result<()> {
        let input = PathBuf::from("part.step");
        let cases = [
            "support",
            "revolution",
            "radial-slot",
            "extrusion",
            "extrusion-core",
            "dominant-extrusion",
            "edge-treatment",
            "boundary-feature",
            "sweep",
        ];
        for name in cases {
            let args = Cli::try_parse_from(["step-inspect", name, "part.step"])?;
            assert_eq!(args.command.input(), &input);
        }
        let healed = Cli::try_parse_from([
            "step-inspect",
            "healed-extrusion-core",
            "part.step",
            "42",
            "11",
            "12",
        ])?;
        assert!(matches!(
            healed.command,
            Analysis::HealedExtrusionCore {
                solid_id: 42,
                cap_a: 11,
                cap_b: 12,
                ..
            }
        ));
        Ok(())
    }

    #[test]
    fn scanner_json_preserves_library_result_schema() -> Result<()> {
        let command = Analysis::Support {
            input: PathBuf::new(),
        };
        let mut output = Vec::new();
        inspect(&command, TORUS, &mut output)?;
        assert_eq!(
            output,
            serde_json::to_vec(&step_redox::detect_solid_surface_signatures_bytes(TORUS)?)?
        );

        let command = Analysis::Revolution {
            input: PathBuf::new(),
        };
        output.clear();
        inspect(&command, TORUS, &mut output)?;
        assert_eq!(
            output,
            serde_json::to_vec(&step_redox::detect_solid_revolutions_bytes(TORUS)?)?
        );
        Ok(())
    }
}
