use anyhow::{Context, Result, ensure};
use occt_wasm::{BoundingBox, OcctKernel, ShapeHandle};
use serde::Serialize;

mod face_flux;
use face_flux::{FaceFluxAccumulator, FaceFluxReport};

const EPSILONS: [f64; 2] = [1e-12, 1e-13];

#[derive(Debug, Serialize)]
pub struct Integral {
    pub volume: f64,
    pub summed_face_error_estimate: f64,
}

#[derive(Debug, Serialize)]
pub struct ReferenceSample {
    pub reference: [f64; 3],
    pub a: Integral,
    pub b: Integral,
    pub delta: f64,
    pub coarse_delta: f64,
    pub a_quadrature_drift: f64,
    pub b_quadrature_drift: f64,
}

#[derive(Debug, Serialize)]
struct ReferenceSpread {
    a: f64,
    b: f64,
}

#[derive(Debug, Serialize)]
pub struct PairedVolume {
    method: &'static str,
    quadrature_eps: [f64; 2],
    absolute_tolerance: f64,
    max_absolute_delta: f64,
    max_quadrature_drift: f64,
    pub within_absolute_tolerance: bool,
    reference_spread: ReferenceSpread,
    pub samples: Vec<ReferenceSample>,
    #[serde(skip_serializing_if = "Option::is_none")]
    face_flux: Option<FaceFluxReport>,
}

// Distill a non-overlapping floating-point expansion using error-free
// FastTwoSum, with a final round-to-nearest-even correction. Unlike an
// ordinary fold, this preserves small signed face terms under cancellation
// and makes summation insensitive to traversal order.
fn accurate_sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let mut partials = Vec::<f64>::new();
    for mut x in values {
        ensure!(x.is_finite(), "nonfinite integration term");
        let mut retained = 0;
        for index in 0..partials.len() {
            let mut y = partials[index];
            if x.abs() < y.abs() {
                std::mem::swap(&mut x, &mut y);
            }
            let hi = x + y;
            ensure!(hi.is_finite(), "integration sum overflow");
            let lo = y - (hi - x);
            if lo != 0.0 {
                partials[retained] = lo;
                retained += 1;
            }
            x = hi;
        }
        partials.truncate(retained);
        if x != 0.0 {
            partials.push(x);
        }
    }
    let mut hi = partials.pop().unwrap_or(0.0);
    while let Some(y) = partials.pop() {
        let x = hi;
        hi = x + y;
        ensure!(hi.is_finite(), "integration sum overflow");
        let lo = y - (hi - x);
        if lo != 0.0 {
            if partials
                .last()
                .is_some_and(|&next| (lo > 0.0 && next > 0.0) || (lo < 0.0 && next < 0.0))
            {
                let doubled = lo * 2.0;
                let corrected = hi + doubled;
                if (corrected - hi).total_cmp(&doubled).is_eq() {
                    hi = corrected;
                }
            }
            break;
        }
    }
    ensure!(hi.is_finite(), "integration sum overflow");
    Ok(hi)
}

fn signed_terms(
    kernel: &mut OcctKernel,
    shape: ShapeHandle,
    point: [f64; 3],
    epsilon: f64,
) -> Result<Vec<f64>> {
    let [x, y, z] = point;
    let terms = kernel.get_signed_face_volume_integrals(shape, x, y, z, epsilon)?;
    let pairs = terms.as_chunks::<2>();
    ensure!(
        !pairs.0.is_empty() && pairs.1.is_empty(),
        "invalid face integration buffer"
    );
    ensure!(
        pairs
            .0
            .iter()
            .all(|[mass, error]| mass.is_finite() && error.is_finite() && *error >= 0.0),
        "invalid signed face integration"
    );
    Ok(terms)
}

fn integrate(terms: &[f64]) -> Result<Integral> {
    let pairs = terms.as_chunks::<2>().0;
    Ok(Integral {
        volume: accurate_sum(pairs.iter().map(|pair| pair[0]))?,
        summed_face_error_estimate: accurate_sum(pairs.iter().map(|pair| pair[1]))?,
    })
}

fn references(a: BoundingBox, b: BoundingBox) -> Result<[[f64; 3]; 4]> {
    let lo_a = [a.min.x, a.min.y, a.min.z];
    let hi_a = [a.max.x, a.max.y, a.max.z];
    let lo_b = [b.min.x, b.min.y, b.min.z];
    let hi_b = [b.max.x, b.max.y, b.max.z];
    for (lo, hi) in [(lo_a, hi_a), (lo_b, hi_b)] {
        ensure!(
            lo.iter()
                .zip(hi)
                .all(|(&x, y)| x.is_finite() && y.is_finite() && x <= y),
            "invalid finite bounds"
        );
    }
    let lo: [f64; 3] = std::array::from_fn(|i| lo_a[i].min(lo_b[i]));
    let hi: [f64; 3] = std::array::from_fn(|i| hi_a[i].max(hi_b[i]));
    let [cx, cy, cz] = std::array::from_fn(|i| f64::midpoint(lo[i], hi[i]));
    Ok([
        [cx, cy, cz],
        [hi[0], cy, cz],
        [cx, hi[1], cz],
        [cx, cy, hi[2]],
    ])
}

fn finite_difference(a: f64, b: f64) -> Result<f64> {
    let delta = a - b;
    ensure!(delta.is_finite(), "nonfinite volume difference");
    Ok(delta)
}

fn spread(samples: &[ReferenceSample], select: impl Fn(&ReferenceSample) -> f64) -> Result<f64> {
    let (first, rest) = samples.split_first().context("no volume samples")?;
    let first = select(first);
    let (low, high) = rest
        .iter()
        .map(select)
        .fold((first, first), |(lo, hi), value| {
            (lo.min(value), hi.max(value))
        });
    finite_difference(high, low)
}

pub fn paired_volume(
    kernel: &mut OcctKernel,
    a: ShapeHandle,
    b: ShapeHandle,
    tolerance: f64,
    include_face_flux: bool,
) -> Result<PairedVolume> {
    ensure!(
        tolerance.is_finite() && tolerance > 0.0,
        "volume tolerance must be finite and positive"
    );
    let bbox_a = kernel.get_optimal_bounds(a, false, false)?;
    let bbox_b = kernel.get_optimal_bounds(b, false, false)?;
    let references = references(bbox_a, bbox_b)?;
    let mut flux = include_face_flux
        .then(|| FaceFluxAccumulator::new(kernel, a, b))
        .transpose()?;
    let mut samples = Vec::with_capacity(references.len());
    let mut max_absolute_delta = 0.0_f64;
    let mut max_quadrature_drift = 0.0_f64;
    for (index, reference) in references.into_iter().enumerate() {
        let coarse_a = signed_terms(kernel, a, reference, EPSILONS[0])?;
        let coarse_b = signed_terms(kernel, b, reference, EPSILONS[0])?;
        let fine_a = signed_terms(kernel, a, reference, EPSILONS[1])?;
        let fine_b = signed_terms(kernel, b, reference, EPSILONS[1])?;
        if let Some(flux) = &mut flux {
            flux.record(index, [&coarse_a, &coarse_b], [&fine_a, &fine_b])?;
        }
        let coarse_a = integrate(&coarse_a)?;
        let coarse_b = integrate(&coarse_b)?;
        let fine_a = integrate(&fine_a)?;
        let fine_b = integrate(&fine_b)?;
        let delta = finite_difference(fine_a.volume, fine_b.volume)?;
        let coarse_delta = finite_difference(coarse_a.volume, coarse_b.volume)?;
        let a_quadrature_drift = finite_difference(fine_a.volume, coarse_a.volume)?.abs();
        let b_quadrature_drift = finite_difference(fine_b.volume, coarse_b.volume)?.abs();
        max_absolute_delta = max_absolute_delta.max(delta.abs()).max(coarse_delta.abs());
        max_quadrature_drift = max_quadrature_drift
            .max(a_quadrature_drift)
            .max(b_quadrature_drift);
        samples.push(ReferenceSample {
            reference,
            a: fine_a,
            b: fine_b,
            delta,
            coarse_delta,
            a_quadrature_drift,
            b_quadrature_drift,
        });
    }
    Ok(PairedVolume {
        method: "shared-reference signed Gauss-Kronrod with spline spans",
        quadrature_eps: EPSILONS,
        absolute_tolerance: tolerance,
        max_absolute_delta,
        max_quadrature_drift,
        within_absolute_tolerance: max_absolute_delta <= tolerance
            && max_quadrature_drift <= tolerance,
        reference_spread: ReferenceSpread {
            a: spread(&samples, |sample| sample.a.volume)?,
            b: spread(&samples, |sample| sample.b.volume)?,
        },
        samples,
        face_flux: flux.map(|flux| flux.finish(references)).transpose()?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sum_retains_small_terms_and_rounds_ties() -> Result<()> {
        for terms in [[1e16, 1.0, -1e16], [-1e16, 1e16, 1.0], [1.0, -1e16, 1e16]] {
            assert!(accurate_sum(terms)?.total_cmp(&1.0).is_eq());
        }
        let half_ulp = f64::EPSILON / 2.0;
        assert!(accurate_sum([1.0, half_ulp])?.total_cmp(&1.0).is_eq());
        assert!(
            accurate_sum([1.0, half_ulp, half_ulp / 2.0])?
                .total_cmp(&(1.0 + f64::EPSILON))
                .is_eq()
        );
        assert!(accurate_sum([f64::NAN]).is_err());
        Ok(())
    }

    #[test]
    fn shared_reference_gate_preserves_signed_volume() -> Result<()> {
        let mut kernel = crate::test_kernel()?;
        let whole = kernel.make_box(2.0, 3.0, 4.0)?;
        let left = kernel.make_box(1.0, 3.0, 4.0)?;
        let right = kernel.translate(left, 1.0, 0.0, 0.0)?;
        let split = kernel.fuse(left, right)?;
        assert!(
            kernel.sub_shape_occurrence_count(split, "face")?
                > kernel.sub_shape_occurrence_count(whole, "face")?
        );
        let tolerance = 1e-10;
        let forward = paired_volume(&mut kernel, whole, split, tolerance, false)?;
        let reverse = paired_volume(&mut kernel, split, whole, tolerance, false)?;
        assert!(forward.within_absolute_tolerance && reverse.within_absolute_tolerance);
        for (a, b) in forward.samples.iter().zip(&reverse.samples) {
            assert!(
                a.reference
                    .iter()
                    .zip(b.reference)
                    .all(|(&x, y)| x.total_cmp(&y).is_eq())
            );
            assert!((a.a.volume - 24.0).abs() <= tolerance);
            assert!((a.b.volume - 24.0).abs() <= tolerance);
            assert!((a.delta + b.delta).abs() <= tolerance);
        }
        for change in [1e-6, 1e-10] {
            let changed = kernel.make_box(2.0, 3.0, 4.0 + change)?;
            let report = paired_volume(&mut kernel, whole, changed, tolerance, false)?;
            assert!(!report.within_absolute_tolerance);
            for sample in report.samples {
                assert!((sample.delta + 6.0 * change).abs() <= tolerance);
            }
        }
        let outer = kernel.make_box(3.0, 3.0, 3.0)?;
        let cavity = kernel.make_box_from_corners(1.0, 1.0, 1.0, 2.0, 2.0, 2.0)?;
        let hollow = kernel.cut(outer, cavity)?;
        let report = paired_volume(&mut kernel, hollow, hollow, tolerance, false)?;
        assert!(report.within_absolute_tolerance);
        assert!(
            report
                .samples
                .iter()
                .all(|sample| (sample.a.volume - 26.0).abs() <= tolerance)
        );
        let faces = kernel.get_sub_shape_occurrences(whole, "face")?;
        let face = *faces.first().context("box has no faces")?;
        assert!(paired_volume(&mut kernel, whole, face, tolerance, false).is_err());
        assert!(paired_volume(&mut kernel, whole, split, 0.0, false).is_err());
        Ok(())
    }
}
