//! Report an alternative mixed analytic/quadrature integral without changing
//! the strict raw quadrature gate. Only independently recognized spherical
//! faces are substituted; every unsupported face keeps its measured terms.

use super::{AnalyticSphere, FaceFlux, SpatialSphere, accurate_sum, finite_difference};
use anyhow::Result;
use serde::Serialize;

#[derive(Debug, Serialize)]
struct Coverage {
    substituted_faces: usize,
    unsupported_spherical_faces: usize,
}

#[derive(Clone, Copy)]
pub(super) enum Method {
    Geodesic,
    Spatial,
}

impl Method {
    const fn label(self) -> &'static str {
        match self {
            Self::Geodesic => {
                "diagnostic analytic geodesic spherical flux plus original quadrature"
            }
            Self::Spatial => "diagnostic spatial spherical boundary flux plus original quadrature",
        }
    }

    fn available(self, face: &FaceFlux) -> bool {
        match self {
            Self::Geodesic => {
                matches!(face.analytic_sphere, Some(AnalyticSphere::Available { .. }))
            }
            Self::Spatial => matches!(face.spatial_sphere, Some(SpatialSphere::Available { .. })),
        }
    }
}

impl Coverage {
    fn from_faces(faces: &[FaceFlux], method: Method) -> Self {
        let mut result = Self {
            substituted_faces: 0,
            unsupported_spherical_faces: 0,
        };
        for face in faces {
            if method.available(face) {
                result.substituted_faces += 1;
            } else if face.surface_kind == "sphere" {
                result.unsupported_spherical_faces += 1;
            }
        }
        result
    }
}

#[derive(Debug, Serialize)]
struct Sample {
    a: f64,
    b: f64,
    delta: f64,
    coarse_delta: f64,
    a_quadrature_drift: f64,
    b_quadrature_drift: f64,
}

#[derive(Debug, Serialize)]
pub(super) struct Substitution {
    method: &'static str,
    a_coverage: Coverage,
    b_coverage: Coverage,
    samples: [Sample; 4],
    max_absolute_delta: f64,
    max_quadrature_drift: f64,
    a_reference_spread: f64,
    b_reference_spread: f64,
}

fn terms(face: &FaceFlux, index: usize, method: Method) -> [f64; 2] {
    match method {
        Method::Geodesic => {
            if let Some(AnalyticSphere::Available { flux }) = &face.analytic_sphere {
                return [flux.volumes()[index]; 2];
            }
        }
        Method::Spatial => {
            if let Some(SpatialSphere::Available { flux }) = &face.spatial_sphere {
                return flux.volumes(index);
            }
        }
    }
    [face.coarse_signed_volume[index], face.signed_volume[index]]
}

fn sample(a: &[FaceFlux], b: &[FaceFlux], index: usize, method: Method) -> Result<Sample> {
    let selected_a = a
        .iter()
        .map(|face| terms(face, index, method))
        .collect::<Vec<_>>();
    let selected_b = b
        .iter()
        .map(|face| terms(face, index, method))
        .collect::<Vec<_>>();
    let paired = |rank| {
        // Sum the signed terms of both shapes directly. Subtracting two
        // rounded whole-solid volumes first would discard low-order bits.
        accurate_sum(
            selected_a
                .iter()
                .map(|v: &[f64; 2]| v[rank])
                .chain(selected_b.iter().map(|v| -v[rank])),
        )
    };
    let drift =
        |terms: &[[f64; 2]]| accurate_sum(terms.iter().flat_map(|&[coarse, fine]| [fine, -coarse]));
    Ok(Sample {
        a: accurate_sum(selected_a.iter().map(|v| v[1]))?,
        b: accurate_sum(selected_b.iter().map(|v| v[1]))?,
        delta: paired(1)?,
        coarse_delta: paired(0)?,
        a_quadrature_drift: drift(&selected_a)?.abs(),
        b_quadrature_drift: drift(&selected_b)?.abs(),
    })
}

fn spread(samples: &[Sample; 4], select: impl Fn(&Sample) -> f64) -> Result<f64> {
    let low = samples.iter().map(&select).fold(f64::INFINITY, f64::min);
    let high = samples.iter().map(select).fold(f64::NEG_INFINITY, f64::max);
    finite_difference(high, low)
}

pub(super) fn summarize(
    a: &[FaceFlux],
    b: &[FaceFlux],
    method: Method,
) -> Result<Option<Substitution>> {
    let a_coverage = Coverage::from_faces(a, method);
    let b_coverage = Coverage::from_faces(b, method);
    if a_coverage.substituted_faces + b_coverage.substituted_faces == 0 {
        return Ok(None);
    }
    let samples = (0..4)
        .map(|index| sample(a, b, index, method))
        .collect::<Result<Vec<_>>>()?;
    let samples: [Sample; 4] = samples
        .try_into()
        .map_err(|_| anyhow::anyhow!("invalid substitution sample count"))?;
    Ok(Some(Substitution {
        method: method.label(),
        a_coverage,
        b_coverage,
        max_absolute_delta: samples
            .iter()
            .map(|s| s.delta.abs().max(s.coarse_delta.abs()))
            .fold(0.0, f64::max),
        max_quadrature_drift: samples
            .iter()
            .map(|s| s.a_quadrature_drift.max(s.b_quadrature_drift))
            .fold(0.0, f64::max),
        a_reference_spread: spread(&samples, |s| s.a)?,
        b_reference_spread: spread(&samples, |s| s.b)?,
        samples,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn measured(value: f64) -> FaceFlux {
        FaceFlux {
            shell_occurrence: 0,
            face_occurrence: 0,
            orientation: "forward".to_owned(),
            surface_kind: "plane".to_owned(),
            area: 0.0,
            bounds: [0.0; 6],
            uv_bounds: [0.0; 4],
            signed_volume: [value; 4],
            error_estimate: [0.0; 4],
            coarse_signed_volume: [value; 4],
            coarse_error_estimate: [0.0; 4],
            analytic_sphere: None,
            spatial_sphere: None,
        }
    }

    #[test]
    fn direct_signed_difference_preserves_small_volume_changes() -> Result<()> {
        let a = [measured(1e16), measured(1.0)];
        let b = [measured(1e16)];
        assert!(summarize(&a, &b, Method::Geodesic)?.is_none());
        let sample = sample(&a, &b, 0, Method::Geodesic)?;
        assert_eq!(sample.a, sample.b);
        assert_eq!(sample.delta, 1.0);
        assert_eq!(sample.coarse_delta, 1.0);
        Ok(())
    }
}
