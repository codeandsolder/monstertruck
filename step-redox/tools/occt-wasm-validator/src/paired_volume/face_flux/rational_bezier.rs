//! Local homogeneous knot refinement and Bernstein product coefficients.
//! Independent of the reconstruction kernel.

use anyhow::{Result, ensure};
use occt_wasm::NurbsCurveData;

pub(super) type Homogeneous = [f64; 4];

struct Curve {
    degree: usize,
    knots: Vec<f64>,
    poles: Vec<Homogeneous>,
}

impl Curve {
    fn insert(&mut self, u: f64) -> Result<()> {
        let p = self.degree;
        let k = self.knots.partition_point(|&knot| knot <= u) - 1;
        let s = self.knots.iter().filter(|&&knot| knot == u).count();
        ensure!(
            k >= p && k < self.poles.len() && s < p,
            "invalid knot refinement"
        );
        let mut poles = vec![[0.0; 4]; self.poles.len() + 1];
        poles[..=k - p].copy_from_slice(&self.poles[..=k - p]);
        poles[k - s + 1..].copy_from_slice(&self.poles[k - s..]);
        for (i, pole) in poles.iter_mut().enumerate().take(k - s + 1).skip(k - p + 1) {
            let alpha = (u - self.knots[i]) / (self.knots[i + p] - self.knots[i]);
            ensure!((0.0..=1.0).contains(&alpha), "invalid refinement blend");
            *pole = std::array::from_fn(|j| {
                self.poles[i][j].mul_add(alpha, self.poles[i - 1][j] * (1.0 - alpha))
            });
        }
        self.poles = poles;
        self.knots.insert(k + 1, u);
        Ok(())
    }

    fn refine(&mut self, u: f64) -> Result<()> {
        while self.knots.iter().filter(|&&knot| knot == u).count() < self.degree {
            self.insert(u)?;
        }
        Ok(())
    }
}

/// Extract all Bezier spans of the actual finite trim (not its whole carrier).
/// Clamping and positive weights make Bernstein convex-hull bounds applicable.
pub(super) fn spans(
    data: &NurbsCurveData,
    center: [f64; 3],
    radius: f64,
    mut first: f64,
    mut last: f64,
) -> Result<Vec<Vec<Homogeneous>>> {
    let p = usize::try_from(data.degree)?;
    ensure!(p > 0 && !data.periodic, "periodic or constant spline");
    ensure!(
        data.knots.len() >= 2 && data.knots.len() == data.multiplicities.len(),
        "invalid knot data"
    );
    ensure!(
        data.knots.iter().all(|x| x.is_finite()) && data.knots.windows(2).all(|x| x[0] < x[1]),
        "invalid knot values"
    );
    let multiplicities = data
        .multiplicities
        .iter()
        .map(|&m| usize::try_from(m))
        .collect::<Result<Vec<_>, _>>()?;
    ensure!(
        multiplicities[0] == p + 1 && multiplicities[multiplicities.len() - 1] == p + 1,
        "unclamped spline"
    );
    ensure!(
        multiplicities[1..multiplicities.len() - 1]
            .iter()
            .all(|&m| m > 0 && m <= p),
        "discontinuous spline"
    );
    let [lo, hi] = [data.knots[0], data.knots[data.knots.len() - 1]];
    let param_tol = 64.0 * f64::EPSILON * (hi.abs().max(lo.abs()).max(hi - lo));
    if (first - lo).abs() <= param_tol {
        first = lo;
    }
    if (last - hi).abs() <= param_tol {
        last = hi;
    }
    ensure!(
        first >= lo && last <= hi && first < last,
        "trim outside spline domain"
    );
    let (points, remainder) = data.poles.as_chunks::<3>();
    ensure!(
        remainder.is_empty() && points.len() > p,
        "invalid pole data"
    );
    ensure!(
        !data.rational || data.weights.len() == points.len(),
        "invalid weights"
    );
    let max_weight = if data.rational {
        ensure!(
            data.weights.iter().all(|w| w.is_finite() && *w > 0.0),
            "nonpositive weights"
        );
        data.weights.iter().copied().fold(0.0, f64::max)
    } else {
        1.0
    };
    let poles = points
        .iter()
        .enumerate()
        .map(|(i, point)| {
            let w = if data.rational {
                data.weights[i] / max_weight
            } else {
                1.0
            };
            [
                (point[0] - center[0]) / radius * w,
                (point[1] - center[1]) / radius * w,
                (point[2] - center[2]) / radius * w,
                w,
            ]
        })
        .collect::<Vec<_>>();
    ensure!(
        poles
            .iter()
            .all(|pole| pole.iter().all(|x| x.is_finite()) && pole[3] > 0.0),
        "invalid homogeneous poles"
    );
    let knots = data
        .knots
        .iter()
        .zip(multiplicities)
        .flat_map(|(&k, m)| std::iter::repeat_n(k, m))
        .collect::<Vec<_>>();
    ensure!(
        knots.len() == poles.len() + p + 1,
        "invalid B-spline dimensions"
    );
    let mut curve = Curve {
        degree: p,
        knots,
        poles,
    };
    for &u in &data.knots[1..data.knots.len() - 1] {
        curve.refine(u)?;
    }
    curve.refine(first)?;
    curve.refine(last)?;
    let spans = (p..curve.poles.len())
        .filter(|&k| {
            curve.knots[k] < curve.knots[k + 1]
                && curve.knots[k] >= first
                && curve.knots[k + 1] <= last
        })
        .map(|k| curve.poles[k - p..=k].to_vec())
        .collect::<Vec<_>>();
    ensure!(!spans.is_empty(), "empty refined trim");
    Ok(spans)
}

fn binomial(n: usize, k: usize) -> f64 {
    let k = k.min(n - k);
    (0..k).fold(1.0, |value, i| value * (n - i) as f64 / (i + 1) as f64)
}

/// Bernstein coefficients of a bilinear product, allowing vector/scalar
/// products without changing the polynomial representation.
pub(super) fn product(
    a: &[Homogeneous],
    b: &[Homogeneous],
    bilinear: impl Fn(Homogeneous, Homogeneous) -> f64,
) -> Result<Vec<f64>> {
    ensure!(!a.is_empty() && !b.is_empty(), "empty Bernstein polynomial");
    let p = a.len() - 1;
    let q = b.len() - 1;
    let mut result = vec![0.0; p + q + 1];
    for (i, &left) in a.iter().enumerate() {
        for (j, &right) in b.iter().enumerate() {
            let weight = binomial(p, i) * binomial(q, j) / binomial(p + q, i + j);
            result[i + j] = weight.mul_add(bilinear(left, right), result[i + j]);
        }
    }
    ensure!(
        result.iter().all(|x| x.is_finite()),
        "nonfinite Bernstein product"
    );
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refine_nonuniform_quadratic_preserves_spans_and_trim() -> Result<()> {
        let data = NurbsCurveData {
            degree: 2,
            rational: false,
            periodic: false,
            knots: vec![0.0, 0.25, 1.0],
            multiplicities: vec![3, 1, 3],
            poles: vec![
                0.0, 0.0, 0.0, 0.125, 0.0, 0.0, 0.625, 0.0, 0.0, 1.0, 0.0, 0.0,
            ],
            weights: vec![],
        };
        let result = spans(&data, [0.0; 3], 1.0, 0.125, 0.75)?;
        assert_eq!(result.len(), 2);
        assert_eq!(result[0][0][0], 0.125);
        assert_eq!(result[0][2][0], 0.25);
        assert_eq!(result[1][0][0], 0.25);
        assert_eq!(result[1][2][0], 0.75);
        Ok(())
    }
}
