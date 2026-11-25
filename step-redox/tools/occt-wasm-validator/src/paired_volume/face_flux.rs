//! Report the actual oriented face occurrences that contribute to the
//! independent signed volume integral, without assuming equal face counts.

use super::{EPSILONS, accurate_sum, finite_difference};
use anyhow::{Context, Result, bail, ensure};
use occt_wasm::{OcctKernel, ShapeHandle};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

mod geodesic;
mod rational_bezier;
mod substitution;

#[derive(Debug, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
enum AnalyticSphere {
    Available { flux: geodesic::GeodesicFlux },
    Unsupported { reason: String },
}

#[derive(Debug, Serialize)]
struct FaceFlux {
    shell_occurrence: usize,
    face_occurrence: usize,
    orientation: String,
    surface_kind: String,
    area: f64,
    bounds: [f64; 6],
    uv_bounds: [f64; 4],
    signed_volume: [f64; 4],
    error_estimate: [f64; 4],
    coarse_signed_volume: [f64; 4],
    coarse_error_estimate: [f64; 4],
    #[serde(skip_serializing_if = "Option::is_none")]
    analytic_sphere: Option<AnalyticSphere>,
}

#[derive(Debug, Serialize)]
struct FamilyFlux {
    surface_kind: String,
    a_faces: usize,
    b_faces: usize,
    a: [f64; 4],
    b: [f64; 4],
    delta: [f64; 4],
    coarse_delta: [f64; 4],
}

#[derive(Debug, Serialize)]
pub struct FaceFluxReport {
    references: [[f64; 3]; 4],
    quadrature_eps: [f64; 2],
    a: Vec<FaceFlux>,
    b: Vec<FaceFlux>,
    surface_families: Vec<FamilyFlux>,
    #[serde(skip_serializing_if = "Option::is_none")]
    analytic_substitution: Option<substitution::Substitution>,
}

fn array<const N: usize>(values: Vec<f64>, label: &str) -> Result<[f64; N]> {
    let values: [f64; N] = values
        .try_into()
        .map_err(|_| anyhow::anyhow!("invalid {label} length"))?;
    ensure!(values.iter().all(|x| x.is_finite()), "nonfinite {label}");
    Ok(values)
}

fn face_metadata(
    kernel: &mut OcctKernel,
    face: ShapeHandle,
    shell_occurrence: usize,
    face_occurrence: usize,
    orientation: String,
) -> Result<FaceFlux> {
    let bounds = kernel.get_optimal_bounds(face, false, false)?;
    let bounds = [
        bounds.min.x,
        bounds.min.y,
        bounds.min.z,
        bounds.max.x,
        bounds.max.y,
        bounds.max.z,
    ];
    ensure!(
        bounds.iter().all(|x| x.is_finite()) && (0..3).all(|i| bounds[i] <= bounds[i + 3]),
        "invalid face bounds"
    );
    let area = kernel.get_surface_area(face)?;
    ensure!(area.is_finite() && area >= 0.0, "invalid face area");
    let surface_kind = kernel.surface_type(face)?;
    let analytic_sphere =
        (surface_kind == "sphere").then(|| match geodesic::evaluate(kernel, face, &orientation) {
            Ok(flux) => AnalyticSphere::Available { flux },
            Err(error) => AnalyticSphere::Unsupported {
                reason: format!("{error:#}"),
            },
        });
    Ok(FaceFlux {
        shell_occurrence,
        face_occurrence,
        orientation,
        surface_kind,
        area,
        bounds,
        uv_bounds: array(kernel.uv_bounds(face)?, "face UV bounds")?,
        signed_volume: [0.0; 4],
        error_estimate: [0.0; 4],
        coarse_signed_volume: [0.0; 4],
        coarse_error_estimate: [0.0; 4],
        analytic_sphere,
    })
}

/// Match the solid/shell/face occurrence traversal in the kernel's
/// signed-face integration routine, including its orientation filter.
/// A count mismatch is an error rather than a guessed face correspondence.
fn face_occurrences(kernel: &mut OcctKernel, shape: ShapeHandle) -> Result<Vec<FaceFlux>> {
    let solids = kernel.get_sub_shape_occurrences(shape, "solid")?;
    let [solid] = solids.as_slice() else {
        bail!("face flux requires exactly one solid occurrence");
    };
    let mut result = Vec::new();
    for (shell_index, shell) in kernel
        .get_sub_shape_occurrences(*solid, "shell")?
        .into_iter()
        .enumerate()
    {
        for (face_index, face) in kernel
            .get_sub_shape_occurrences(shell, "face")?
            .into_iter()
            .enumerate()
        {
            let orientation = kernel.shape_orientation(face)?;
            match orientation.as_str() {
                "forward" | "reversed" => result.push(face_metadata(
                    kernel,
                    face,
                    shell_index,
                    face_index,
                    orientation,
                )?),
                "internal" | "external" => (),
                _ => bail!("unexpected face orientation {orientation}"),
            }
        }
    }
    ensure!(!result.is_empty(), "no oriented face occurrences");
    Ok(result)
}

fn record_terms(faces: &mut [FaceFlux], index: usize, coarse: &[f64], fine: &[f64]) -> Result<()> {
    let coarse = coarse.as_chunks::<2>().0;
    let fine = fine.as_chunks::<2>().0;
    ensure!(
        coarse.len() == faces.len() && fine.len() == faces.len(),
        "face integral counts ({}, {}) do not match oriented occurrences {}",
        coarse.len(),
        fine.len(),
        faces.len()
    );
    for ((face, &[coarse_mass, coarse_error]), &[mass, error]) in
        faces.iter_mut().zip(coarse).zip(fine)
    {
        face.signed_volume[index] = mass;
        face.error_estimate[index] = error;
        face.coarse_signed_volume[index] = coarse_mass;
        face.coarse_error_estimate[index] = coarse_error;
    }
    Ok(())
}

fn group(faces: &[FaceFlux]) -> BTreeMap<&str, Vec<&FaceFlux>> {
    let mut groups = BTreeMap::<&str, Vec<&FaceFlux>>::new();
    for face in faces {
        groups.entry(&face.surface_kind).or_default().push(face);
    }
    groups
}

fn totals(faces: &[&FaceFlux], select: impl Fn(&FaceFlux) -> &[f64; 4]) -> Result<[f64; 4]> {
    let values = (0..4)
        .map(|i| accurate_sum(faces.iter().map(|face| select(face)[i])))
        .collect::<Result<Vec<_>>>()?;
    array(values, "surface-family signed volumes")
}

fn families(a: &[FaceFlux], b: &[FaceFlux]) -> Result<Vec<FamilyFlux>> {
    let a = group(a);
    let b = group(b);
    let kinds = a.keys().chain(b.keys()).copied().collect::<BTreeSet<_>>();
    kinds
        .into_iter()
        .map(|kind| {
            let left = a.get(kind).map_or(&[][..], Vec::as_slice);
            let right = b.get(kind).map_or(&[][..], Vec::as_slice);
            let a = totals(left, |face| &face.signed_volume)?;
            let b = totals(right, |face| &face.signed_volume)?;
            let delta = difference(a, b)?;
            let coarse_delta = difference(
                totals(left, |face| &face.coarse_signed_volume)?,
                totals(right, |face| &face.coarse_signed_volume)?,
            )?;
            Ok(FamilyFlux {
                surface_kind: kind.to_owned(),
                a_faces: left.len(),
                b_faces: right.len(),
                a,
                b,
                delta,
                coarse_delta,
            })
        })
        .collect()
}

fn difference(a: [f64; 4], b: [f64; 4]) -> Result<[f64; 4]> {
    let values = a
        .into_iter()
        .zip(b)
        .map(|(a, b)| finite_difference(a, b))
        .collect::<Result<Vec<_>>>()?;
    array(values, "surface-family deltas")
}

/// Record the exact terms used by paired volume, not a second integration
/// with a potentially different shape state or duplicate expensive calls.
pub(super) struct FaceFluxAccumulator {
    a: Vec<FaceFlux>,
    b: Vec<FaceFlux>,
}

impl FaceFluxAccumulator {
    pub(super) fn new(kernel: &mut OcctKernel, a: ShapeHandle, b: ShapeHandle) -> Result<Self> {
        Ok(Self {
            a: face_occurrences(kernel, a).context("reference face flux")?,
            b: face_occurrences(kernel, b).context("candidate face flux")?,
        })
    }

    pub(super) fn record(
        &mut self,
        index: usize,
        coarse: [&[f64]; 2],
        fine: [&[f64]; 2],
    ) -> Result<()> {
        record_terms(&mut self.a, index, coarse[0], fine[0])?;
        record_terms(&mut self.b, index, coarse[1], fine[1])
    }

    pub(super) fn finish(mut self, references: [[f64; 3]; 4]) -> Result<FaceFluxReport> {
        for face in self.a.iter_mut().chain(&mut self.b) {
            if let Some(AnalyticSphere::Available { flux }) = &mut face.analytic_sphere {
                flux.set_references(references);
            }
        }
        let surface_families = families(&self.a, &self.b)?;
        let analytic_substitution = substitution::summarize(&self.a, &self.b)?;
        Ok(FaceFluxReport {
            references,
            quadrature_eps: EPSILONS,
            a: self.a,
            b: self.b,
            surface_families,
            analytic_substitution,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An analytic box face's signed volume is (n·(x-reference)) A/3.
    /// Checking its actual sign and magnitude, rather than only a total or
    /// axis dependence, detects swapped occurrences of parallel faces.
    fn check_box_face(
        face: &FaceFlux,
        references: [[f64; 3]; 4],
        low: [f64; 3],
        high: [f64; 3],
        sense: f64,
    ) -> Result<()> {
        let axis = (0..3)
            .find(|&i| (face.bounds[i] - face.bounds[i + 3]).abs() <= 1e-10)
            .context("missing constant face coordinate")?;
        let position = face.bounds[axis];
        let outward = if (position - low[axis]).abs() <= 1e-10 {
            -sense
        } else {
            ensure!(
                (position - high[axis]).abs() <= 1e-10,
                "unexpected box face location"
            );
            sense
        };
        for (index, reference) in references.into_iter().enumerate() {
            let expected = outward * face.area * (position - reference[axis]) / 3.0;
            assert!((face.signed_volume[index] - expected).abs() <= 1e-10);
            assert!((face.coarse_signed_volume[index] - expected).abs() <= 1e-10);
        }
        Ok(())
    }

    #[test]
    fn face_terms_follow_oriented_box_occurrences() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let whole = kernel.make_box(2.0, 3.0, 4.0)?;
        let left = kernel.make_box(1.0, 3.0, 4.0)?;
        let right = kernel.translate(left, 1.0, 0.0, 0.0)?;
        let split = kernel.fuse(left, right)?;
        let paired = super::super::paired_volume(&mut kernel, whole, split, 1e-10, true)?;
        assert!(paired.within_absolute_tolerance);
        let report = paired.face_flux.context("missing requested face terms")?;
        assert_eq!(report.a.len(), 6);
        assert!(report.b.len() > report.a.len());
        let [plane] = report.surface_families.as_slice() else {
            bail!("box should have only planar faces");
        };
        assert_eq!(plane.a_faces, 6);
        assert_eq!(plane.b_faces, report.b.len());
        for (a, b) in plane.a.iter().zip(plane.b) {
            assert!((*a - 24.0).abs() <= 1e-10);
            assert!((b - 24.0).abs() <= 1e-10);
        }
        for face in &report.a {
            check_box_face(face, report.references, [0.0; 3], [2.0, 3.0, 4.0], 1.0)?;
        }

        // A closed solid with a cavity has multiple shells and negative
        // inner-face fluxes. The cumulative orientation, occurrence order and
        // total must remain consistent, not be deduplicated by shape identity.
        let outer = kernel.make_box(3.0, 3.0, 3.0)?;
        let inner = kernel.make_box_from_corners(1.0, 1.0, 1.0, 2.0, 2.0, 2.0)?;
        let hollow = kernel.cut(outer, inner)?;
        let paired = super::super::paired_volume(&mut kernel, hollow, hollow, 1e-10, true)?;
        let report = paired.face_flux.context("missing cavity face terms")?;
        assert_eq!(report.a.len(), 12);
        assert_eq!(
            report
                .a
                .iter()
                .map(|face| face.shell_occurrence)
                .collect::<BTreeSet<_>>()
                .len(),
            2
        );
        for face in &report.a {
            let cavity = face
                .bounds
                .iter()
                .all(|value| *value >= 1.0 - 1e-10 && *value <= 2.0 + 1e-10);
            let (low, high, sense) = if cavity {
                ([1.0; 3], [2.0; 3], -1.0)
            } else {
                ([0.0; 3], [3.0; 3], 1.0)
            };
            check_box_face(face, report.references, low, high, sense)?;
        }
        assert!(
            paired
                .samples
                .iter()
                .all(|sample| (sample.a.volume - 26.0).abs() <= 1e-10)
        );
        Ok(())
    }
}
