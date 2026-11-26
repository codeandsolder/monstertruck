//! Embedded Gauss--Kronrod quadrature for independent spatial boundary
//! integrals. Error estimates are empirical, not interval bounds.

use anyhow::{Result, ensure};

const NODES: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_5,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_5,
    0.0,
];
const KRONROD: [f64; 8] = [
    0.022_935_322_010_529_22,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_2,
    0.140_653_259_715_525_9,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_4,
    0.204_432_940_075_298_9,
    0.209_482_141_084_727_8,
];
const GAUSS: [f64; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_7,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

#[derive(Debug)]
pub(super) struct Integral<const N: usize> {
    pub(super) value: [f64; N],
    pub(super) estimated_error: f64,
    pub(super) component_errors: [f64; N],
}

fn rule<const N: usize>(
    integrand: &impl Fn(f64) -> Result<[f64; N]>,
    start: f64,
    end: f64,
) -> Result<Integral<N>> {
    let center = f64::midpoint(start, end);
    let half_width = (end - start) * 0.5;
    let central = integrand(center)?;
    let mut fine = central.map(|v| v * KRONROD[7]);
    let mut coarse = central.map(|v| v * GAUSS[3]);
    for i in 0..7 {
        let offset = NODES[i] * half_width;
        let left = integrand(center - offset)?;
        let right = integrand(center + offset)?;
        for j in 0..N {
            let pair = left[j] + right[j];
            fine[j] = KRONROD[i].mul_add(pair, fine[j]);
            if i % 2 == 1 {
                coarse[j] = GAUSS[i / 2].mul_add(pair, coarse[j]);
            }
        }
    }
    let value = fine.map(|v| v * half_width);
    let component_errors = std::array::from_fn(|i| ((fine[i] - coarse[i]) * half_width).abs());
    let estimated_error = component_errors.into_iter().fold(0.0, f64::max);
    ensure!(
        value.iter().all(|v| v.is_finite()) && estimated_error.is_finite(),
        "nonfinite spatial boundary integral"
    );
    Ok(Integral {
        value,
        estimated_error,
        component_errors,
    })
}

fn adaptive<const N: usize>(
    integrand: &impl Fn(f64) -> Result<[f64; N]>,
    start: f64,
    end: f64,
    tolerance: f64,
    depth: u8,
) -> Result<Integral<N>> {
    let integral = rule(integrand, start, end)?;
    // Do not recursively split on floating-point roundoff alone.
    let scale = integral.value.iter().map(|v| v.abs()).fold(0.0, f64::max);
    if integral.estimated_error <= tolerance.max(32.0 * f64::EPSILON * scale) {
        return Ok(integral);
    }
    ensure!(depth < 20, "spatial boundary quadrature did not converge");
    let center = f64::midpoint(start, end);
    let left = adaptive(integrand, start, center, tolerance * 0.5, depth + 1)?;
    let right = adaptive(integrand, center, end, tolerance * 0.5, depth + 1)?;
    Ok(Integral {
        value: std::array::from_fn(|i| left.value[i] + right.value[i]),
        estimated_error: left.estimated_error + right.estimated_error,
        component_errors: std::array::from_fn(|i| {
            left.component_errors[i] + right.component_errors[i]
        }),
    })
}

pub(super) fn integrate<const N: usize>(
    integrand: &impl Fn(f64) -> Result<[f64; N]>,
    tolerance: f64,
) -> Result<Integral<N>> {
    adaptive(integrand, 0.0, 1.0, tolerance, 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_nonpolynomial_integrand() -> Result<()> {
        let integral = integrate(&|x| Ok([(1.0 + x * x).recip(), x.exp()]), 1e-14)?;
        assert!((integral.value[0] - std::f64::consts::FRAC_PI_4).abs() < 1e-14);
        assert!((integral.value[1] - 1.0_f64.exp_m1()).abs() < 1e-14);
        Ok(())
    }
}
