//! Common recurrence for one B-spline basis window or a complete derivative jet.

use super::{BasisWindow, KnotVector, inv_or_zero_strict};
use crate::nurbs::BasisSlice;
use smallvec::SmallVec;

type Values = SmallVec<[f64; 32]>;

/// All nonzero derivative ranks in one coefficient buffer. Each row retains
/// its recurrence sentinel, while active slices expose only valid basis values.
pub(in crate::nurbs) struct BasisJet<const N: usize> {
    coefficients: Values,
    stride: usize,
    active_ranks: usize,
    start: usize,
    len: usize,
}

impl<const N: usize> BasisJet<N> {
    pub(in crate::nurbs) fn windows(&self) -> [BasisSlice<'_>; N] {
        std::array::from_fn(|rank| {
            if rank >= self.active_ranks {
                BasisSlice {
                    start: 0,
                    values: &[],
                }
            } else {
                let offset = rank * self.stride;
                BasisSlice {
                    start: self.start,
                    values: &self.coefficients[offset..offset + self.len],
                }
            }
        })
    }
}

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

    pub(super) fn jet<const N: usize>(&self) -> BasisJet<N> {
        // Derivative rank r shares the first degree-r polynomial steps of
        // every lower-order derivative. Each row becomes the next prefix
        // before its own derivative steps overwrite it. This allocates one
        // coefficient block even when a high-degree jet spills to the heap.
        let stride = self.degree + 2;
        let active_ranks = N.min(self.degree + 1);
        let mut coefficients: Values = smallvec::smallvec![0.0; stride * active_ranks];
        let mut completed = 0;
        if active_ranks != 0 {
            coefficients[(active_ranks - 1) * stride + self.span - self.start] = 1.0;
        }
        for derivative in (0..active_ranks).rev() {
            let offset = derivative * stride;
            let end = offset + stride;
            let polynomial_rank = self.degree - derivative;
            for rank in completed + 1..=polynomial_rank {
                self.polynomial_step(rank, &mut coefficients[offset..end]);
            }
            completed = polynomial_rank;
            if derivative != 0 {
                coefficients.copy_within(offset..end, offset - stride);
            }
            for rank in polynomial_rank + 1..=self.degree {
                self.derivative_step(rank, &mut coefficients[offset..end]);
            }
        }
        BasisJet {
            coefficients,
            stride,
            active_ranks,
            start: self.start,
            len: (self.degree + 1).min(self.last - self.degree - self.start),
        }
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
                for (rank, window) in jet.windows().iter().enumerate() {
                    let single = knots.bspline_basis_functions_unchecked(degree, rank, parameter);
                    assert_eq!(window.start, single.start_index());
                    assert!(
                        window
                            .values
                            .iter()
                            .map(|v| v.to_bits())
                            .eq(single.values().iter().map(|v| v.to_bits()))
                    );
                }
            }
        }
    }
}
