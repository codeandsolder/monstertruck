//! Refine already discovered curve parameters to numerical convergence.
//!
//! A modeling/welding tolerance is an acceptance bound, not a sensible
//! stopping criterion for a Newton solve. Retain the best iterate and do
//! not change either input curve to make an approximate vertex fit.
use monstertruck_geometry::prelude::*;

const MAX_REFINEMENTS: usize = 32;
const MAX_BACKTRACKS: usize = 16;

fn improve<const N: usize>(
    mut parameters: [f64; N],
    ranges: [(f64, f64); N],
    residual: impl Fn([f64; N]) -> Vector3,
    direction: impl Fn([f64; N], Vector3) -> Option<[f64; N]>,
) -> Option<[f64; N]> {
    if parameters.iter().any(|value| !value.is_finite())
        || ranges
            .iter()
            .any(|&(lo, hi)| !lo.is_finite() || !hi.is_finite() || lo >= hi)
    {
        return None;
    }
    for (parameter, (lo, hi)) in parameters.iter_mut().zip(ranges) {
        *parameter = parameter.clamp(lo, hi);
    }
    let mut gap = residual(parameters);
    let mut error = gap.magnitude2();
    if !error.is_finite() {
        return None;
    }
    for _ in 0..MAX_REFINEMENTS {
        if error == 0. {
            break;
        }
        let Some(step) = direction(parameters, gap) else {
            break;
        };
        if step.iter().any(|value| !value.is_finite()) {
            break;
        }
        let mut scale = 1.;
        let mut improved = false;
        for _ in 0..MAX_BACKTRACKS {
            let trial = std::array::from_fn(|i| {
                step[i]
                    .mul_add(scale, parameters[i])
                    .clamp(ranges[i].0, ranges[i].1)
            });
            let trial_gap = residual(trial);
            let trial_error = trial_gap.magnitude2();
            if trial_error < error {
                parameters = trial;
                gap = trial_gap;
                error = trial_error;
                improved = true;
                break;
            }
            scale *= 0.5;
        }
        if !improved {
            break;
        }
    }
    Some(parameters)
}

/// Improve the seed of a point-to-curve inverse. Gauss–Newton minimizes
/// the full spatial residual, with backtracking for curved or clamped steps.
pub(super) fn point_parameter<C: BoundedCurve<Point = Point3, Vector = Vector3>>(
    curve: &C,
    point: Point3,
    seed: f64,
) -> Option<f64> {
    let [parameter] = improve(
        [seed],
        [curve.range_tuple()],
        |[t]| curve.evaluate(t) - point,
        |[t], residual| {
            let tangent = curve.derivative(t);
            let norm = tangent.magnitude2();
            (norm > 0.).then(|| [-tangent.dot(residual) / norm])
        },
    )?;
    Some(parameter)
}

/// Refine a crossing using a two-parameter least-squares Newton step. Nearly
/// parallel tangents may make the solve singular: keep the best iterate
/// and let the caller decide whether its actual residual is acceptable.
pub(super) fn crossing_parameters<C: BoundedCurve<Point = Point3, Vector = Vector3>>(
    lhs: &C,
    rhs: &C,
    seeds: [f64; 2],
) -> Option<[f64; 2]> {
    improve(
        seeds,
        [lhs.range_tuple(), rhs.range_tuple()],
        |[u, v]| lhs.evaluate(u) - rhs.evaluate(v),
        |[u, v], gap| {
            let a = lhs.derivative(u);
            let b = rhs.derivative(v);
            let aa = a.magnitude2();
            let bb = b.magnitude2();
            let ab = -a.dot(b);
            let determinant = aa.mul_add(bb, -(ab * ab));
            if !determinant.is_finite() || determinant <= f64::EPSILON * aa * bb {
                return None;
            }
            let ga = -a.dot(gap);
            let gb = b.dot(gap);
            Some([
                ga.mul_add(bb, -(gb * ab)) / determinant,
                aa.mul_add(gb, -(ab * ga)) / determinant,
            ])
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossing_refines_an_already_weldable_cubic_seed() {
        let cubic = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(3),
            vec![
                Point3::new(0., 0., 0.),
                Point3::new(1., 0., 0.),
                Point3::new(2., 1., 0.),
                Point3::new(3., 2., 0.),
            ],
        );
        let line = BsplineCurve::new_unchecked(
            KnotVector::bezier_knot(1),
            vec![Point3::new(1., -1., 0.), Point3::new(1., 2., 0.)],
        );
        let seed = [1. / 3. + 1.0e-7, 1. / 3. + 8. / 81.];
        let initial = cubic.evaluate(seed[0]) - line.evaluate(seed[1]);
        assert!(initial.magnitude() < 1.0e-6);
        assert!(initial.magnitude() > 1.0e-8);
        let [u, v] = crossing_parameters(&cubic, &line, seed).unwrap();
        assert!((u - 1. / 3.).abs() < 1.0e-14);
        assert!(cubic.evaluate(u).distance(line.evaluate(v)) < 1.0e-14);
        let t = point_parameter(&cubic, cubic.evaluate(u), seed[0]).unwrap();
        assert!((t - u).abs() < 1.0e-14);
    }

    #[test]
    fn clamped_and_singular_steps_keep_best_finite_iterate() {
        let a = Line(Point3::new(0., 0., 0.), Point3::new(1., 0., 0.));
        let b = Line(Point3::new(0., 1., 0.), Point3::new(1., 1., 0.));
        let [u, v] = crossing_parameters(&a, &b, [0.5, 0.5]).unwrap();
        assert!((a.evaluate(u).distance(b.evaluate(v)) - 1.).abs() < f64::EPSILON);
        let end = point_parameter(&a, Point3::new(2., 0., 0.), 0.5).unwrap();
        assert!(a.evaluate(end).distance(Point3::new(1., 0., 0.)) < f64::EPSILON);
    }
}
