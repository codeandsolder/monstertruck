//! Per-solid Boolean residual diagnostics. Pairing follows the release
//! inventory matcher, rather than silently truncating a sorted zip.
use super::{Solid, compare::match_solids};
use anyhow::{Context, Result, ensure};
use clap::Args;
use occt_wasm::{OcctKernel, ShapeHandle};
use serde::Serialize;
use std::{path::PathBuf, process::ExitCode};

#[derive(Debug, Args)]
pub struct BooleanArgs {
    reference: PathBuf,
    candidate: PathBuf,
    /// Absolute tolerance in cubic millimetres; no relative tolerance fallback.
    #[arg(long)]
    absolute_tolerance: f64,
}

struct Inventory {
    handles: Vec<ShapeHandle>,
    solids: Vec<Solid>,
}

impl Inventory {
    fn measure(kernel: &mut OcctKernel, input: ShapeHandle) -> Result<Self> {
        let handles = kernel.get_sub_shape_occurrences(input, "solid")?;
        let solids = handles
            .iter()
            .map(|&solid| {
                ensure!(kernel.is_valid(solid)?, "invalid BRep solid");
                Solid::measure(kernel, solid)
            })
            .collect::<Result<_>>()?;
        Ok(Self { handles, solids })
    }
}

#[derive(Debug, Serialize)]
struct Residual {
    reference_index: usize,
    candidate_index: usize,
    reference_volume: f64,
    candidate_volume: f64,
    common_volume: f64,
    reference_minus_candidate: f64,
    candidate_minus_reference: f64,
    common_missing_volume: f64,
    pass: bool,
}

#[derive(Debug, Serialize)]
struct Report {
    absolute_tolerance: f64,
    reference_solids: usize,
    candidate_solids: usize,
    /// Matching is deterministic metric matching, not optimal material matching.
    pairs: Vec<Residual>,
    #[serde(skip_serializing_if = "Option::is_none")]
    unmatched_reason: Option<&'static str>,
    pass: bool,
}

pub(crate) fn result_volume(kernel: &mut OcctKernel, shape: ShapeHandle) -> Result<f64> {
    let volume = kernel.get_b_rep_properties(shape)?.volume;
    kernel.release(shape)?;
    ensure!(volume.is_finite(), "invalid Boolean result volume");
    Ok(volume)
}

fn compare(
    kernel: &mut OcctKernel,
    a: ShapeHandle,
    b: ShapeHandle,
    tolerance: f64,
) -> Result<Report> {
    ensure!(
        tolerance.is_finite() && tolerance > 0.0,
        "volume tolerance must be finite and positive"
    );
    let reference = Inventory::measure(kernel, a)?;
    let candidate = Inventory::measure(kernel, b)?;
    let mut report = Report {
        absolute_tolerance: tolerance,
        reference_solids: reference.solids.len(),
        candidate_solids: candidate.solids.len(),
        pairs: Vec::new(),
        unmatched_reason: None,
        pass: false,
    };
    if report.reference_solids != report.candidate_solids || reference.solids.is_empty() {
        report.unmatched_reason =
            Some("nonempty equal-size solid inventories are required for per-solid diagnostics");
        return Ok(report);
    }
    for (i, m) in match_solids(&reference.solids, &candidate.solids, [0.0; 3])? {
        let j = m.index;
        let a = reference.handles[i];
        let b = candidate.handles[j];
        let va = reference.solids[i].volume;
        let vb = candidate.solids[j].volume;
        let common = kernel
            .common(a, b)
            .with_context(|| format!("common of solid pair {i}/{j}"))?;
        let common_volume = result_volume(kernel, common)?;
        let ab = kernel
            .cut(a, b)
            .with_context(|| format!("reference minus candidate solid pair {i}/{j}"))?;
        let a_minus_b = result_volume(kernel, ab)?;
        let ba = kernel
            .cut(b, a)
            .with_context(|| format!("candidate minus reference solid pair {i}/{j}"))?;
        let b_minus_a = result_volume(kernel, ba)?;
        let missing = (va - common_volume).abs().max((vb - common_volume).abs());
        ensure!(missing.is_finite(), "nonfinite common volume residual");
        report.pairs.push(Residual {
            reference_index: i,
            candidate_index: j,
            reference_volume: va,
            candidate_volume: vb,
            common_volume,
            reference_minus_candidate: a_minus_b,
            candidate_minus_reference: b_minus_a,
            common_missing_volume: missing,
            pass: [missing, a_minus_b.abs(), b_minus_a.abs()]
                .into_iter()
                .all(|v| v <= tolerance),
        });
    }
    report.pairs.sort_by_key(|pair| pair.reference_index);
    report.pass = report.pairs.iter().all(|pair| pair.pass);
    Ok(report)
}

pub fn run(args: &BooleanArgs) -> Result<ExitCode> {
    let mut kernel = OcctKernel::new()?;
    let (a, _) = crate::load(&mut kernel, &args.reference)?;
    let (b, _) = crate::load(&mut kernel, &args.candidate)?;
    let report = compare(&mut kernel, a.shape, b.shape, args.absolute_tolerance)?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(if report.pass {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn directed_cuts_common_and_reordered_inventory() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let small = kernel.make_box(1.0, 2.0, 3.0)?;
        let large = kernel.make_box_from_corners(5.0, 0.0, 0.0, 7.0, 2.0, 3.0)?;
        let a = kernel.make_compound(&[small, large])?;
        let b = kernel.make_compound(&[large, small])?;
        let report = compare(&mut kernel, a, b, 1e-9)?;
        assert!(report.pass);
        assert_eq!(report.pairs.len(), 2);
        assert_eq!(report.pairs[0].candidate_index, 1);
        assert_eq!(report.pairs[1].candidate_index, 0);
        let moved = kernel.translate(small, 0.5, 0.0, 0.0)?;
        let report = compare(&mut kernel, small, moved, 1e-9)?;
        assert!(!report.pass);
        let pair = report.pairs.first().context("no paired solid")?;
        for value in [
            pair.common_volume,
            pair.reference_minus_candidate,
            pair.candidate_minus_reference,
            pair.common_missing_volume,
        ] {
            assert!((value - 3.0).abs() < 1e-9);
        }
        let different_count = compare(&mut kernel, small, a, 1e-9)?;
        assert!(!different_count.pass);
        assert!(different_count.unmatched_reason.is_some());
        assert!(compare(&mut kernel, small, small, f64::NAN).is_err());
        Ok(())
    }
}
