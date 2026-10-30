//! convert from truck shape json to step file.
//!
//! ### usage
//!
//! ```bash
//! shape-to-step <input shape file> [output shape file]
//! ```

use anyhow::{Context, Result};
use clap::Parser;
use monstertruck_io::step::save;
use monstertruck_modeling::{Curve, Point3, Surface};
use monstertruck_topology::compress::CompressedSolid;
use std::io::Write;

/// convert from truck shape json to step file.
#[derive(Parser, Debug)]
struct Args {
    /// file name of truck shape json
    input_shape_file: String,
    /// step file name
    #[arg(default_value = "output.stp")]
    output_step_file: String,
}

fn main() -> Result<()> {
    let Args {
        input_shape_file,
        output_step_file,
    } = Args::parse();

    let shape_file = std::fs::read(&input_shape_file)
        .with_context(|| format!("failed to read {input_shape_file}"))?;
    let compressed: CompressedSolid<Point3, Curve, Surface> =
        serde_json::from_slice(&shape_file)
            .with_context(|| format!("failed to parse {input_shape_file}"))?;
    let step_string = save::CompleteStepDisplay::new(
        save::StepModel::from(&compressed),
        save::StepHeaderDescriptor {
            organization_system: "shape-to-step".to_owned(),
            ..save::StepHeaderDescriptor::default()
        },
    )
    .to_string();
    let mut step_file = std::fs::File::create(&output_step_file)
        .with_context(|| format!("failed to create {output_step_file}"))?;
    step_file
        .write_all(step_string.as_bytes())
        .with_context(|| format!("failed to write {output_step_file}"))?;
    step_p21::parser::parse(&step_string).context("generated STEP failed self-parse")?;
    Ok(())
}
