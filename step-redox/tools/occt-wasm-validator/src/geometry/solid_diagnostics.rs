//! Independent single-solid diagnostics. Import entity IDs from OCCT's parsed
//! STEP model, never textual order or an approximate shape correspondence.
use crate::{
    geometry::{self, boolean::result_volume},
    paired_volume,
};
use anyhow::{Context, Result, ensure};
use clap::Args;
use occt_wasm::{OcctKernel, ShapeHandle};
use serde::Serialize;
use std::{
    fs,
    path::{Path, PathBuf},
    process::ExitCode,
};

const ABS_GEOM_TOL: f64 = 1e-9;
const ABS_VOLUME_TOL: f64 = 1e-10;
const ABS_AREA_TOL: f64 = 1e-9;
const REL_TOL: f64 = 1e-10;
const ADAPTIVE_EPSILON: f64 = 1e-8;

#[derive(Debug, Args)]
pub struct SolidArgs {
    reference: PathBuf,
    candidate: PathBuf,
    /// Exact textual STEP #identifier, not model rank, for the reference shape.
    #[arg(long)]
    reference_entity_id: Option<u32>,
    /// Exact textual STEP #identifier, not model rank, for the candidate shape.
    #[arg(long)]
    candidate_entity_id: Option<u32>,
}

struct Selected {
    shape: ShapeHandle,
    entity_id: Option<u32>,
    rank: Option<u32>,
}

fn load_selected(kernel: &mut OcctKernel, path: &Path, id: Option<u32>) -> Result<Selected> {
    let bytes = fs::read(path).with_context(|| format!("read {}", path.display()))?;
    let (shape, rank) = match id {
        Some(entity_id) => {
            let imported = kernel
                .import_step_entity_bytes(&bytes, entity_id)
                .with_context(|| {
                    format!("transfer STEP entity #{entity_id} from {}", path.display())
                })?;
            (imported.shape, Some(imported.rank))
        }
        None => (
            kernel
                .import_step_bytes(&bytes)
                .with_context(|| format!("import {}", path.display()))?
                .shape,
            None,
        ),
    };
    ensure!(
        kernel.sub_shape_occurrence_count(shape, "solid")? == 1,
        "{}: selected shape must contain exactly one solid",
        path.display()
    );
    Ok(Selected {
        shape,
        entity_id: id,
        rank,
    })
}

#[derive(Debug, Serialize)]
struct Measured {
    #[serde(flatten)]
    metrics: geometry::Solid,
    center: [f64; 3],
    bbox: [f64; 6],
    brep_valid: bool,
}

impl Measured {
    fn new(kernel: &mut OcctKernel, shape: ShapeHandle) -> Result<Self> {
        let metrics = geometry::Solid::measure(kernel, shape)?;
        let props = kernel.get_b_rep_properties(shape)?;
        let bounds = geometry::Bounds::new(&kernel.get_bounding_box_loose(shape, true)?)?;
        ensure!(
            props.center_of_mass.iter().all(|x| x.is_finite()),
            "nonfinite center of mass"
        );
        Ok(Self {
            metrics,
            center: props.center_of_mass,
            bbox: bounds.bbox,
            brep_valid: kernel.is_valid(shape)?,
        })
    }
}

#[derive(Debug, Serialize)]
struct Integral {
    volume: f64,
    /// OCCT returns a relative estimate; it is not a certified numerical bound.
    estimated_error: f64,
}

#[derive(Debug, Serialize)]
struct AdaptiveAgreement {
    a: Integral,
    b: Integral,
    delta: f64,
    consistency_tolerance: f64,
    consistent: bool,
}

fn volume_tolerance(a: f64, b: f64) -> f64 { ABS_VOLUME_TOL.max(REL_TOL * a.abs().max(b.abs())) }

fn adaptive_agreement(
    kernel: &mut OcctKernel,
    a: ShapeHandle,
    b: ShapeHandle,
) -> Result<AdaptiveAgreement> {
    let source = kernel.get_adaptive_volume(a, ADAPTIVE_EPSILON, true, true)?;
    let candidate = kernel.get_adaptive_volume(b, ADAPTIVE_EPSILON, true, true)?;
    let delta = (source.volume - candidate.volume).abs();
    // Preserve the existing native-oracle allowance for both the returned
    // relative estimate and a potentially absolute interpretation.
    let allowance = source
        .relative_error
        .max(source.relative_error * source.volume.abs())
        + candidate
            .relative_error
            .max(candidate.relative_error * candidate.volume.abs());
    let tolerance = volume_tolerance(source.volume, candidate.volume) + allowance;
    ensure!(
        [delta, allowance, tolerance]
            .into_iter()
            .all(|x| x.is_finite()),
        "nonfinite adaptive integration comparison"
    );
    Ok(AdaptiveAgreement {
        a: Integral {
            volume: source.volume,
            estimated_error: source.relative_error,
        },
        b: Integral {
            volume: candidate.volume,
            estimated_error: candidate.relative_error,
        },
        delta,
        consistency_tolerance: tolerance,
        consistent: delta <= tolerance,
    })
}

#[derive(Debug, Serialize)]
struct Errors {
    volume: f64,
    area: f64,
    center: f64,
    bbox_report_only: f64,
}

fn max_abs_difference<const N: usize>(a: [f64; N], b: [f64; N]) -> Result<f64> {
    let maximum = a
        .into_iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    ensure!(maximum.is_finite(), "nonfinite geometric difference");
    Ok(maximum)
}

#[derive(Debug, Serialize)]
struct BooleanResidual {
    a_minus_b: Option<f64>,
    b_minus_a: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    a_minus_b_error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    b_minus_a_error: Option<String>,
}

fn directed_cut(kernel: &mut OcctKernel, a: ShapeHandle, b: ShapeHandle) -> Result<f64> {
    let shape = kernel.cut(a, b)?;
    ensure!(
        kernel.is_valid(shape)?,
        "Boolean cut produced an invalid BRep"
    );
    Ok(result_volume(kernel, shape)?.abs())
}

fn boolean_residual(kernel: &mut OcctKernel, a: ShapeHandle, b: ShapeHandle) -> BooleanResidual {
    let mut run = |lhs, rhs| match directed_cut(kernel, lhs, rhs) {
        Ok(volume) => (Some(volume), None),
        Err(error) => (None, Some(format!("{error:#}"))),
    };
    let (a_minus_b, a_minus_b_error) = run(a, b);
    let (b_minus_a, b_minus_a_error) = run(b, a);
    BooleanResidual {
        a_minus_b,
        b_minus_a,
        a_minus_b_error,
        b_minus_a_error,
    }
}

#[derive(Debug, Serialize)]
#[serde(untagged)]
enum PairedDiagnostic {
    Measured(paired_volume::PairedVolume),
    Failed {
        within_absolute_tolerance: bool,
        error: String,
    },
}

#[derive(Debug, Serialize)]
struct Tolerances {
    geometry_abs: f64,
    volume_abs: f64,
    boolean_volume_effective: f64,
    area_abs: f64,
    relative: f64,
}

#[derive(Debug, Serialize)]
struct Report {
    reference_entity_id: Option<u32>,
    reference_entity_rank: Option<u32>,
    candidate_entity_id: Option<u32>,
    candidate_entity_rank: Option<u32>,
    a: Measured,
    b: Measured,
    geometry_equal: bool,
    adaptive_volume: AdaptiveAgreement,
    paired_absolute_volume: PairedDiagnostic,
    errors: Errors,
    boolean_residual_volume: BooleanResidual,
    tolerances: Tolerances,
}

fn compare(kernel: &mut OcctKernel, a: Selected, b: Selected) -> Result<Report> {
    let measured_a = Measured::new(kernel, a.shape)?;
    let measured_b = Measured::new(kernel, b.shape)?;
    let volume_error = (measured_a.metrics.volume - measured_b.metrics.volume).abs();
    let area_error = (measured_a.metrics.area - measured_b.metrics.area).abs();
    ensure!(
        volume_error.is_finite() && area_error.is_finite(),
        "nonfinite solid properties"
    );
    let errors = Errors {
        volume: volume_error,
        area: area_error,
        center: max_abs_difference(measured_a.center, measured_b.center)?,
        bbox_report_only: max_abs_difference(measured_a.bbox, measured_b.bbox)?,
    };
    // This deliberately does not require matching face counts: valid
    // subdivided surfaces can represent exactly the same material.
    let adaptive = adaptive_agreement(kernel, a.shape, b.shape)?;
    let residual = boolean_residual(kernel, a.shape, b.shape);
    let boolean_tolerance = volume_tolerance(measured_a.metrics.volume, measured_b.metrics.volume);
    let geometry_equal = adaptive.consistent
        && measured_a.brep_valid
        && measured_b.brep_valid
        && residual
            .a_minus_b
            .is_some_and(|value| value <= boolean_tolerance)
        && residual
            .b_minus_a
            .is_some_and(|value| value <= boolean_tolerance);
    // Paired absolute integration is a stricter independent diagnostic,
    // never a replacement for the established adaptive/Boolean decision.
    let paired = match paired_volume::paired_volume(kernel, a.shape, b.shape, ABS_VOLUME_TOL) {
        Ok(result) => PairedDiagnostic::Measured(result),
        Err(error) => PairedDiagnostic::Failed {
            within_absolute_tolerance: false,
            error: format!("{error:#}"),
        },
    };
    Ok(Report {
        reference_entity_id: a.entity_id,
        reference_entity_rank: a.rank,
        candidate_entity_id: b.entity_id,
        candidate_entity_rank: b.rank,
        a: measured_a,
        b: measured_b,
        geometry_equal,
        adaptive_volume: adaptive,
        paired_absolute_volume: paired,
        errors,
        boolean_residual_volume: residual,
        tolerances: Tolerances {
            geometry_abs: ABS_GEOM_TOL,
            volume_abs: ABS_VOLUME_TOL,
            boolean_volume_effective: boolean_tolerance,
            area_abs: ABS_AREA_TOL,
            relative: REL_TOL,
        },
    })
}

pub fn run(args: &SolidArgs) -> Result<ExitCode> {
    let mut kernel = OcctKernel::new()?;
    let a = load_selected(&mut kernel, &args.reference, args.reference_entity_id)?;
    let b = load_selected(&mut kernel, &args.candidate, args.candidate_entity_id)?;
    let report = compare(&mut kernel, a, b)?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(if report.geometry_equal {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn selected(shape: ShapeHandle) -> Selected {
        Selected {
            shape,
            entity_id: None,
            rank: None,
        }
    }

    #[test]
    fn independent_geometry_gates_and_report_only_paired_volume() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let whole = kernel.make_box(2., 3., 4.)?;
        let equal = compare(&mut kernel, selected(whole), selected(whole))?;
        assert!(equal.geometry_equal);
        assert!(equal.adaptive_volume.consistent);
        assert_eq!(equal.boolean_residual_volume.a_minus_b, Some(0.));
        assert!(matches!(
            equal.paired_absolute_volume,
            PairedDiagnostic::Measured(_)
        ));
        let left = kernel.make_box(1., 3., 4.)?;
        let right = kernel.translate(left, 1., 0., 0.)?;
        let subdivided = kernel.fuse(left, right)?;
        let split = compare(&mut kernel, selected(whole), selected(subdivided))?;
        assert!(split.geometry_equal);
        assert!(split.a.metrics.topology.faces != split.b.metrics.topology.faces);
        let shifted = kernel.translate(whole, 0.5, 0., 0.)?;
        let negative = compare(&mut kernel, selected(whole), selected(shifted))?;
        assert!(!negative.geometry_equal);
        assert!(negative.adaptive_volume.consistent); // same volume but wrong material placement
        assert!(
            (negative
                .boolean_residual_volume
                .a_minus_b
                .context("missing cut")?
                - 6.)
                .abs()
                < 1e-8
        );
        let changed = kernel.make_box(2., 3., 5.)?;
        let invalid = compare(&mut kernel, selected(whole), selected(changed))?;
        assert!(!invalid.adaptive_volume.consistent);
        assert!(!invalid.geometry_equal);
        assert!(invalid.errors.volume > 0.);
        Ok(())
    }

    #[test]
    fn source_entity_identity_and_single_solid_precondition() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let a = kernel.make_box(2., 3., 4.)?;
        let b = kernel.make_box(1., 2., 3.)?;
        let compound = kernel.make_compound(&[a, b])?;
        let dir = tempfile::tempdir()?;
        let source = dir.path().join("compound.step");
        let candidate = dir.path().join("box.step");
        let text = kernel.export_step(compound)?;
        fs::write(&source, &text)?;
        fs::write(&candidate, kernel.export_step(a)?)?;
        assert!(load_selected(&mut kernel, &source, None).is_err());
        let mut chosen = None;
        for record in text.split(';') {
            if let Some((label, body)) = record.split_once('=')
                && body.trim_start().starts_with("MANIFOLD_SOLID_BREP(")
                && let Some(id) = label.trim().strip_prefix('#')
            {
                let id: u32 = id.parse()?;
                let imported = load_selected(&mut kernel, &source, Some(id))?;
                if (kernel.get_b_rep_properties(imported.shape)?.volume - 24.).abs() < 1e-9 {
                    chosen = Some(id);
                }
            }
        }
        let source_id = chosen.context("missing requested solid entity")?;
        let selected_source = load_selected(&mut kernel, &source, Some(source_id))?;
        let selected_target = load_selected(&mut kernel, &candidate, None)?;
        assert!(selected_source.rank.is_some());
        let comparison = compare(&mut kernel, selected_source, selected_target)?;
        assert!(comparison.geometry_equal);
        assert_eq!(comparison.reference_entity_id, Some(source_id));
        assert!(load_selected(&mut kernel, &source, Some(u32::MAX)).is_err());
        let _ = PathBuf::new();
        Ok(())
    }
}
