mod geometry;
mod paired_volume;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};
use num_traits::ToPrimitive;
use occt_wasm::{OcctKernel, StepImport, wasmtime::Engine};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::{Path, PathBuf},
    process::ExitCode,
};

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Report independent OCCT geometry metrics in the native oracle schema.
    Analyze { input: PathBuf },
    /// Analyze multiple STEP files with one compiled WASM module and isolated instances.
    AnalyzeBatch {
        #[arg(required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,
    },
    /// Compare independent geometry, renderings, and optional strict mesh hashes.
    Compare {
        reference: PathBuf,
        candidate: PathBuf,
        #[arg(long)]
        thresholds: PathBuf,
        #[arg(long, value_enum, default_value = "identity")]
        alignment: geometry::compare::Alignment,
        #[arg(long, value_delimiter = ',', value_enum)]
        views: Vec<geometry::render::View>,
        #[arg(long)]
        out: PathBuf,
    },
    /// Compare saved native or WASM reports, without initializing OCCT.
    CompareMetrics {
        reference: PathBuf,
        candidate: PathBuf,
        #[arg(long)]
        thresholds: PathBuf,
        #[arg(long, value_enum, default_value = "identity")]
        alignment: geometry::compare::Alignment,
    },
    /// Hash the strict world-coordinate tessellation at 1e-9 mm precision.
    MeshHash { input: PathBuf },
    /// Compare signed face integrals at two quadrature tolerances and four references.
    PairedVolume {
        reference: PathBuf,
        candidate: PathBuf,
        #[arg(long)]
        absolute_tolerance: f64,
    },
}

fn load(kernel: &mut OcctKernel, path: &Path) -> Result<(StepImport, usize)> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let import = kernel
        .import_step_bytes(&bytes)
        .with_context(|| format!("import {}", path.display()))?;
    Ok((import, bytes.len()))
}

fn strict_mesh_hash(coordinates: &[f64], quantum: f64) -> Result<String> {
    let scale = quantum.recip();
    ensure!(
        quantum.is_finite() && quantum > 0.0 && scale.is_finite(),
        "invalid mesh hash quantum"
    );
    let triangles = coordinates.as_chunks::<9>();
    ensure!(triangles.1.is_empty(), "incomplete triangle data");
    let mut tris = Vec::with_capacity(triangles.0.len());
    for triangle in triangles.0 {
        let mut points = [[0_i64; 3]; 3];
        for (point, coords) in points.iter_mut().zip(triangle.as_chunks::<3>().0) {
            for (value, &coordinate) in point.iter_mut().zip(coords) {
                // Reject instead of saturating an unrepresentable hash coordinate.
                // Match the reference oracle: multiply by the reciprocal
                // before rounding rather than dividing each coordinate.
                *value = (coordinate * scale)
                    .round_ties_even()
                    .to_i64()
                    .context("mesh coordinate exceeds hash range")?;
            }
        }
        points.sort_unstable();
        tris.push(points);
    }
    tris.sort_unstable();
    let mut hasher = Sha256::new();
    for tri in tris {
        for point in tri {
            for value in point {
                hasher.update(value.to_le_bytes());
            }
        }
    }
    Ok(hex::encode(hasher.finalize()))
}

fn analyze_path(kernel: &mut OcctKernel, input: &Path) -> Result<geometry::GeometryReport> {
    let (import, bytes) = load(kernel, input)?;
    let report = geometry::analyze(kernel, input, bytes, import)?;
    kernel.release(import.shape)?;
    Ok(report)
}

#[derive(Serialize)]
struct MeshComparison {
    reference: String,
    candidate: String,
    pass: bool,
}

#[derive(Serialize)]
struct FullComparison {
    reference: geometry::GeometryReport,
    candidate: geometry::GeometryReport,
    geometry: geometry::compare::Comparison,
    renders: Vec<geometry::render::RenderComparison>,
    #[serde(skip_serializing_if = "Option::is_none")]
    strict_mesh_hash: Option<MeshComparison>,
    pass: bool,
}

fn compare_steps(
    reference: &Path,
    candidate: &Path,
    thresholds: &geometry::compare::Thresholds,
    alignment: geometry::compare::Alignment,
    views: &[geometry::render::View],
    out: &Path,
) -> Result<FullComparison> {
    let mut kernel = OcctKernel::new()?;
    let (a, a_bytes) = load(&mut kernel, reference)?;
    let (b, b_bytes) = load(&mut kernel, candidate)?;
    let ref_report = geometry::analyze(&mut kernel, reference, a_bytes, a)?;
    let cand_report = geometry::analyze(&mut kernel, candidate, b_bytes, b)?;
    let geometry = geometry::compare::compare(&ref_report, &cand_report, alignment, thresholds)?;
    let ref_tris = kernel.tessellate_triangles64(a.shape, 0.03, 0.3)?;
    let cand_tris = kernel.tessellate_triangles64(b.shape, 0.03, 0.3)?;
    let renders = geometry::render::compare(
        &ref_tris,
        &cand_tris,
        geometry.alignment_translation_mm,
        thresholds,
        views,
        &out.join("renders"),
    )?;
    let strict_mesh_hash = if thresholds.strict_mesh_hash {
        let reference =
            strict_mesh_hash(&kernel.tessellate_triangles64(a.shape, 0.002, 0.2)?, 1e-9)?;
        let candidate =
            strict_mesh_hash(&kernel.tessellate_triangles64(b.shape, 0.002, 0.2)?, 1e-9)?;
        Some(MeshComparison {
            pass: reference == candidate,
            reference,
            candidate,
        })
    } else {
        None
    };
    let pass = geometry.metrics_pass
        && renders.iter().all(|r| r.pass)
        && strict_mesh_hash.as_ref().is_none_or(|r| r.pass);
    Ok(FullComparison {
        reference: ref_report,
        candidate: cand_report,
        geometry,
        renders,
        strict_mesh_hash,
        pass,
    })
}

fn main() -> Result<ExitCode> {
    let cli = Cli::parse();
    match cli.command {
        Command::CompareMetrics {
            reference,
            candidate,
            thresholds,
            alignment,
        } => {
            let reference = serde_json::from_reader(fs::File::open(reference)?)?;
            let candidate = serde_json::from_reader(fs::File::open(candidate)?)?;
            let thresholds = serde_json::from_reader(fs::File::open(thresholds)?)?;
            let report =
                geometry::compare::compare(&reference, &candidate, alignment, &thresholds)?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            Ok(if report.metrics_pass {
                ExitCode::SUCCESS
            } else {
                ExitCode::FAILURE
            })
        }
        Command::Analyze { input } => {
            let mut kernel = OcctKernel::new()?;
            let report = analyze_path(&mut kernel, &input)?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            Ok(ExitCode::SUCCESS)
        }
        Command::AnalyzeBatch { inputs } => {
            let engine = Engine::new(&OcctKernel::engine_config())?;
            let module = OcctKernel::compile(&engine)?;
            let reports = inputs
                .iter()
                .map(|path| {
                    let mut kernel = OcctKernel::from_module(&module)?;
                    analyze_path(&mut kernel, path)
                })
                .collect::<Result<Vec<_>>>()?;
            println!("{}", serde_json::to_string_pretty(&reports)?);
            Ok(ExitCode::SUCCESS)
        }
        Command::Compare {
            reference,
            candidate,
            thresholds,
            alignment,
            views,
            out,
        } => {
            let thresholds = serde_json::from_reader(fs::File::open(thresholds)?)?;
            let report =
                compare_steps(&reference, &candidate, &thresholds, alignment, &views, &out)?;
            let json = serde_json::to_string_pretty(&report)?;
            fs::write(out.join("report.json"), &json)?;
            println!("{json}");
            Ok(if report.pass {
                ExitCode::SUCCESS
            } else {
                ExitCode::FAILURE
            })
        }
        Command::MeshHash { input } => {
            let mut kernel = OcctKernel::new()?;
            let (import, _) = load(&mut kernel, &input)?;
            let triangles = kernel.tessellate_triangles64(import.shape, 0.002, 0.2)?;
            println!("{}", strict_mesh_hash(&triangles, 1e-9)?);
            Ok(ExitCode::SUCCESS)
        }
        Command::PairedVolume {
            reference,
            candidate,
            absolute_tolerance,
        } => {
            let mut kernel = OcctKernel::new()?;
            let (a, _) = load(&mut kernel, &reference)?;
            let (b, _) = load(&mut kernel, &candidate)?;
            let report =
                paired_volume::paired_volume(&mut kernel, a.shape, b.shape, absolute_tolerance)?;
            println!("{}", serde_json::to_string_pretty(&report)?);
            Ok(if report.within_absolute_tolerance {
                ExitCode::SUCCESS
            } else {
                ExitCode::FAILURE
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hash_is_order_independent_but_retains_picometer_coordinates() -> Result<()> {
        let triangles = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let reversed = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert_eq!(
            strict_mesh_hash(&triangles, 1e-9)?,
            strict_mesh_hash(&reversed, 1e-9)?
        );
        let mut changed = triangles;
        changed[3] += 2e-9;
        assert_ne!(
            strict_mesh_hash(&triangles, 1e-9)?,
            strict_mesh_hash(&changed, 1e-9)?
        );
        assert!(strict_mesh_hash(&[f64::INFINITY; 9], 1e-9).is_err());
        assert!(strict_mesh_hash(&[0.0; 8], 1e-9).is_err());
        assert!(strict_mesh_hash(&triangles, 0.0).is_err());
        Ok(())
    }
}
