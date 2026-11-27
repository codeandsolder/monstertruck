//! Embedded Gauss--Kronrod quadrature for independent spatial boundary
//! integrals. Error estimates are empirical, not interval bounds.

use anyhow::{Context, Result, ensure};
use std::{cmp::Ordering, collections::BinaryHeap};

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
    absolute_integrals: [f64; N],
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
    let mut absolute = central.map(|v| v.abs() * KRONROD[7]);
    for i in 0..7 {
        let offset = NODES[i] * half_width;
        let left = integrand(center - offset)?;
        let right = integrand(center + offset)?;
        for j in 0..N {
            let pair = left[j] + right[j];
            fine[j] = KRONROD[i].mul_add(pair, fine[j]);
            absolute[j] = KRONROD[i].mul_add(left[j].abs() + right[j].abs(), absolute[j]);
            if i % 2 == 1 {
                coarse[j] = GAUSS[i / 2].mul_add(pair, coarse[j]);
            }
        }
    }
    let value = fine.map(|v| v * half_width);
    let component_errors = std::array::from_fn(|i| ((fine[i] - coarse[i]) * half_width).abs());
    let estimated_error = component_errors.into_iter().fold(0.0, f64::max);
    let integral = Integral {
        value,
        estimated_error,
        component_errors,
        absolute_integrals: absolute.map(|v| v * half_width),
    };
    integral.check_finite()?;
    Ok(integral)
}

const MAX_DEPTH: u8 = 20;

impl<const N: usize> Integral<N> {
    fn check_finite(&self) -> Result<()> {
        ensure!(
            self.estimated_error.is_finite()
                && self
                    .value
                    .iter()
                    .chain(&self.component_errors)
                    .chain(&self.absolute_integrals)
                    .all(|v| v.is_finite()),
            "nonfinite spatial boundary integral"
        );
        Ok(())
    }

    fn combine(left: &Self, right: &Self) -> Self {
        Self {
            value: std::array::from_fn(|i| left.value[i] + right.value[i]),
            estimated_error: left.estimated_error + right.estimated_error,
            component_errors: std::array::from_fn(|i| {
                left.component_errors[i] + right.component_errors[i]
            }),
            absolute_integrals: std::array::from_fn(|i| {
                left.absolute_integrals[i] + right.absolute_integrals[i]
            }),
        }
    }

    fn unconverged_component(&self, tolerance: f64) -> Option<usize> {
        (0..N)
            .filter_map(|i| {
                let target = tolerance.max(32.0 * f64::EPSILON * self.absolute_integrals[i]);
                let error = self.component_errors[i];
                (error > target).then_some((i, error / target))
            })
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(i, _)| i)
    }
}

struct Interval<const N: usize> {
    start: f64,
    end: f64,
    depth: u8,
    parent: Option<usize>,
    children: Option<[usize; 2]>,
    integral: Integral<N>,
}

#[derive(Clone, Copy)]
struct Priority {
    error: f64,
    index: usize,
}

impl PartialEq for Priority {
    fn eq(&self, other: &Self) -> bool { self.cmp(other) == Ordering::Equal }
}
impl Eq for Priority {}
impl PartialOrd for Priority {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}
impl Ord for Priority {
    fn cmp(&self, other: &Self) -> Ordering {
        self.error
            .total_cmp(&other.error)
            .then_with(|| self.index.cmp(&other.index))
    }
}

struct Partition<const N: usize> {
    nodes: Vec<Interval<N>>,
    queues: [BinaryHeap<Priority>; N],
}

impl<const N: usize> Partition<N> {
    fn insert(&mut self, interval: Interval<N>) -> usize {
        let index = self.nodes.len();
        if interval.depth < MAX_DEPTH {
            for (queue, &error) in self
                .queues
                .iter_mut()
                .zip(&interval.integral.component_errors)
            {
                queue.push(Priority { error, index });
            }
        }
        self.nodes.push(interval);
        index
    }

    fn refine(
        &mut self,
        component: usize,
        integrand: &impl Fn(f64) -> Result<[f64; N]>,
    ) -> Result<bool> {
        let index = loop {
            let Some(priority) = self.queues[component].pop() else {
                return Ok(false);
            };
            if self.nodes[priority.index].children.is_none() {
                if priority.error == 0.0 {
                    return Ok(false);
                }
                break priority.index;
            }
        };
        let node = &self.nodes[index];
        let center = f64::midpoint(node.start, node.end);
        let depth = node.depth + 1;
        let start = node.start;
        let end = node.end;
        let left = self.insert(Interval {
            start,
            end: center,
            depth,
            parent: Some(index),
            children: None,
            integral: rule(integrand, start, center)?,
        });
        let right = self.insert(Interval {
            start: center,
            end,
            depth,
            parent: Some(index),
            children: None,
            integral: rule(integrand, center, end)?,
        });
        self.nodes[index].children = Some([left, right]);
        // Recompute only the ancestors. Pairwise sums avoid accumulated
        // subtract/add drift as intervals are replaced by their children.
        let mut current = Some(index);
        while let Some(i) = current {
            let [a, b] = self.nodes[i]
                .children
                .context("missing quadrature children")?;
            self.nodes[i].integral =
                Integral::combine(&self.nodes[a].integral, &self.nodes[b].integral);
            self.nodes[i].integral.check_finite()?;
            current = self.nodes[i].parent;
        }
        Ok(true)
    }
}

pub(super) fn integrate<const N: usize>(
    integrand: &impl Fn(f64) -> Result<[f64; N]>,
    tolerance: f64,
) -> Result<Integral<N>> {
    ensure!(
        tolerance.is_finite() && tolerance > 0.0,
        "invalid quadrature tolerance"
    );
    let root = rule(integrand, 0.0, 1.0)?;
    if root.unconverged_component(tolerance).is_none() {
        return Ok(root);
    }
    let mut partition = Partition {
        nodes: Vec::new(),
        queues: std::array::from_fn(|_| BinaryHeap::new()),
    };
    partition.insert(Interval {
        start: 0.0,
        end: 1.0,
        depth: 0,
        parent: None,
        children: None,
        integral: root,
    });
    // The requested error applies to the complete integral, not to an
    // equal share of it on every dyadic leaf. Refine the largest remaining
    // contributor to the worst unmet component. A depth-limited leaf keeps
    // contributing its measured error and never vetoes the whole integral
    // when the sum of all leaf estimates already meets the global budget.
    while let Some(component) = partition.nodes[0].integral.unconverged_component(tolerance) {
        ensure!(
            partition.refine(component, integrand)?,
            "spatial boundary quadrature did not converge: component {component}, total error {:e}, target {tolerance:e}",
            partition.nodes[0].integral.component_errors[component]
        );
    }
    Ok(partition.nodes.swap_remove(0).integral)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn global_budget_handles_localized_error_without_tightening_each_leaf() -> Result<()> {
        let split = 0.123_456_7;
        let tolerance = 1e-7;
        let integral = integrate(
            &|x| Ok([if x < split { 1.0 } else { 0.0 }, x * x]),
            tolerance,
        )?;
        assert!((integral.value[0] - split).abs() <= tolerance);
        assert!((integral.value[1] - 1.0 / 3.0).abs() <= tolerance);
        assert!(
            integral
                .component_errors
                .iter()
                .all(|&error| error <= tolerance)
        );
        Ok(())
    }

    #[test]
    fn resolves_nonpolynomial_integrand() -> Result<()> {
        let integral = integrate(&|x| Ok([(1.0 + x * x).recip(), x.exp()]), 1e-14)?;
        assert!((integral.value[0] - std::f64::consts::FRAC_PI_4).abs() < 1e-14);
        assert!((integral.value[1] - 1.0_f64.exp_m1()).abs() < 1e-14);
        Ok(())
    }
}
