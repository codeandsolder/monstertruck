mod paired_volume;

use anyhow::{Context, Result, ensure};
use clap::{Parser, Subcommand};
use num_traits::ToPrimitive;
use occt_wasm::{OcctKernel, ShapeHandle};
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
    /// Report independent OCCT geometry metrics and exact triangle hash.
    Analyze { input: PathBuf },
    /// Compare signed face integrals at two quadrature tolerances and four references.
    PairedVolume {
        reference: PathBuf,
        candidate: PathBuf,
        #[arg(long)]
        absolute_tolerance: f64,
    },
}

fn load(kernel: &mut OcctKernel, path: &Path) -> Result<ShapeHandle> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    kernel
        .import_step_bytes(&bytes)
        .with_context(|| format!("import {}", path.display()))
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
        for (point, coords) in points.iter_mut().zip(triangle.chunks_exact(3)) {
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
    Ok(format!("{:x}", hasher.finalize()))
}

fn analyze(kernel: &mut OcctKernel, path: &Path) -> Result<()> {
    let shape = load(kernel, path)?;
    let bbox = kernel.get_bounding_box_loose(shape, true)?;
    let volume = kernel.get_volume(shape)?;
    let area = kernel.get_surface_area(shape)?;
    let faces = kernel.sub_shape_occurrence_count(shape, "face")?;
    let edges = kernel.sub_shape_occurrence_count(shape, "edge")?;
    let vertices = kernel.sub_shape_occurrence_count(shape, "vertex")?;
    let triangles = kernel.tessellate_triangles64(shape, 0.002, 0.2)?;
    let mesh_hash = strict_mesh_hash(&triangles, 1e-9)?;
    println!("path={}", path.display());
    println!(
        "bbox={:.12},{:.12},{:.12},{:.12},{:.12},{:.12}",
        bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z
    );
    println!("volume={volume:.15}\narea={area:.15}");
    println!("faces={faces} edges={edges} vertices={vertices}");
    println!("triangles={}\nmesh_hash={mesh_hash}", triangles.len() / 9);
    Ok(())
}

fn main() -> Result<ExitCode> {
    let cli = Cli::parse();
    let mut kernel = OcctKernel::new()?;
    match cli.command {
        Command::Analyze { input } => {
            analyze(&mut kernel, &input)?;
            Ok(ExitCode::SUCCESS)
        }
        Command::PairedVolume {
            reference,
            candidate,
            absolute_tolerance,
        } => {
            let a = load(&mut kernel, &reference)?;
            let b = load(&mut kernel, &candidate)?;
            let report = paired_volume::paired_volume(&mut kernel, a, b, absolute_tolerance)?;
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
