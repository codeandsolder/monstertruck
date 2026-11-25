use anyhow::Result;
use clap::{Parser, Subcommand};
use std::io::Write as _;
use std::{path::PathBuf, process::ExitCode};

#[derive(Parser)]
#[command(about = "Scan and extract deduplicatable MANIFOLD_SOLID_BREP bodies")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Build, validate and inspect a deduplicated body corpus.
    Corpus {
        #[command(subcommand)]
        action: step_redox::body_corpus::batch::CorpusCommand,
    },
    /// Emit one JSON record per manifold solid; unsupported canonical identities fall back file-locally.
    Scan { input: PathBuf },
    /// Extract selected bodies as standalone STEP files.
    Extract {
        input: PathBuf,
        out_dir: PathBuf,
        /// Selection formatted as `DATA_SECTION:SOLID_ID:FINGERPRINT`.
        #[arg(long = "solid")]
        solids: Vec<String>,
    },
}

fn write_json(value: &impl serde::Serialize) -> Result<()> {
    let mut output = std::io::stdout().lock();
    serde_json::to_writer(&mut output, value)?;
    output.write_all(b"\n")?;
    Ok(())
}

fn main() -> Result<ExitCode> {
    let cli = Cli::parse();
    match cli.command {
        Command::Corpus { action } => return step_redox::body_corpus::batch::run(action),
        Command::Scan { input } => {
            let report = step_redox::body_corpus::scan(&input)?;
            write_json(&report)?;
        }
        Command::Extract {
            input,
            out_dir,
            solids,
        } => {
            let objects = step_redox::body_corpus::extract_selected(&input, &out_dir, &solids)?;
            write_json(&objects)?;
        }
    }
    Ok(ExitCode::SUCCESS)
}
