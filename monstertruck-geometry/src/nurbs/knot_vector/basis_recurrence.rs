//! Common recurrence for one B-spline basis window or a complete derivative jet.

use super::{BasisWindow, KnotVector, inv_or_zero_strict};
use smallvec::SmallVec;

type Values = SmallVec<[f64; 32]>;

pub(super) struct BasisRecurrence<'a> {
    knots: &'a KnotVector,
    degree: usize,
    parameter: f64,
    last: usize,
    span: usize,
    start: usize,
}

impl<'a> BasisRecurrence<'a> {
    pub(super) fn new(knots: &'a KnotVector, degree: usize, parameter: f64) -> Self {
        let last = knots.len() - 1;
        // Extend the first/last nonempty span outside the knot range, exactly
        // as in the single-window evaluator.
        let index = knots
            .floor(parameter)
            .unwrap_or_else(|| knots.strict_multiplicity(0).saturating_sub(1));
        let span = if index == last {
            last.saturating_sub(knots.strict_multiplicity(last))
        } else {
            index
        };
        Self {
            knots,
            degree,
            parameter,
            last,
            span,
            start: span.saturating_sub(degree),
        }
    }

    fn initial(&self) -> Values {
        // The trailing sentinel is read by the last recurrence iteration.
        let mut values = smallvec::smallvec![0.0; self.degree + 2];
        values[self.span - self.start] = 1.0;
        values
    }

    #[inline]
    fn polynomial_step(&self, rank: usize, values: &mut [f64]) {
        let base = self.span.saturating_sub(rank);
        let delta = self.knots[base + rank] - self.knots[base];
        let mut left = inv_or_zero_strict(delta) * (self.parameter - self.knots[base]);
        for index in base..=usize::min(self.span, self.last - rank - 1) {
            let delta = self.knots[index + rank + 1] - self.knots[index + 1];
            let right = inv_or_zero_strict(delta) * (self.knots[index + rank + 1] - self.parameter);
            let offset = index - self.start;
            values[offset] = right.mul_add(values[offset + 1], left * values[offset]);
            left = 1.0 - right;
        }
    }

    #[inline]
    fn derivative_step(&self, rank: usize, values: &mut [f64]) {
        let base = self.span.saturating_sub(rank);
        let delta = self.knots[base + rank] - self.knots[base];
        let mut left = inv_or_zero_strict(delta);
        for index in base..=usize::min(self.span, self.last - rank - 1) {
            let delta = self.knots[index + rank + 1] - self.knots[index + 1];
            let right = inv_or_zero_strict(delta);
            let offset = index - self.start;
            values[offset] = right.mul_add(-values[offset + 1], left * values[offset])
                * crate::usize_to_f64(rank);
            left = right;
        }
    }

    fn finish(&self, mut values: Values) -> BasisWindow {
        let count = self.last - self.degree;
        values.truncate(usize::min(self.degree + 1, count - self.start));
        BasisWindow::new(self.start, values, count)
    }

    pub(super) fn evaluate(&self, derivative: usize) -> BasisWindow {
        let mut values = self.initial();
        let polynomial_rank = self.degree - derivative;
        for rank in 1..=polynomial_rank {
            self.polynomial_step(rank, &mut values);
        }
        for rank in polynomial_rank + 1..=self.degree {
            self.derivative_step(rank, &mut values);
        }
        self.finish(values)
    }

    pub(super) fn jet<const N: usize>(&self) -> [BasisWindow; N] {
        // Derivative rank r shares the first degree-r polynomial steps of
        // every lower-order derivative. Visit ranks in descending order and
        // extend this prefix in place. The remaining derivative recurrence
        // uses a private copy, so no arithmetic or summation order changes.
        let mut windows = std::array::from_fn(|_| BasisWindow::empty(self.last - self.degree));
        let mut prefix = self.initial();
        let mut completed = 0;
        for (derivative, window) in windows.iter_mut().enumerate().take(self.degree + 1).rev() {
            let polynomial_rank = self.degree - derivative;
            for rank in completed + 1..=polynomial_rank {
                self.polynomial_step(rank, &mut prefix);
            }
            completed = polynomial_rank;
            let mut values = if derivative == 0 {
                std::mem::take(&mut prefix)
            } else {
                prefix.clone()
            };
            for rank in polynomial_rank + 1..=self.degree {
                self.derivative_step(rank, &mut values);
            }
            *window = self.finish(values);
        }
        windows
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_recurrence_matches_individual_windows_bit_for_bit() {
        for (degree, knots) in [
            (0, KnotVector::uniform_knot(0, 4)),
            (1, KnotVector::uniform_knot(1, 4)),
            (3, KnotVector::uniform_knot(3, 4)),
            (
                3,
                KnotVector::from(vec![0., 0., 0., 0., 0.25, 0.25, 0.75, 1., 1., 1., 1.]),
            ),
            (2, KnotVector::from(vec![0., 0.25, 0.5, 0.75, 1.])),
            (40, KnotVector::bezier_knot(40)),
        ] {
            for parameter in [-0.1, 0., 0.125, 0.25, 0.5, 0.75, 0.999, 1., 1.1] {
                let jet = knots.bspline_basis_jet::<5>(degree, parameter);
                for (rank, window) in jet.iter().enumerate() {
                    let single = knots.bspline_basis_functions_unchecked(degree, rank, parameter);
                    assert_eq!(window.start_index(), single.start_index());
                    assert_eq!(window.total_length(), single.total_length());
                    assert!(
                        window
                            .values()
                            .iter()
                            .map(|v| v.to_bits())
                            .eq(single.values().iter().map(|v| v.to_bits()))
                    );
                }
            }
        }
    }
}
