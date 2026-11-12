//! Deterministic nearest-seed index for a cached rational surface sample grid.
//!
//! The index only accelerates the discrete seed search. Newton refinement
//! continues to use the exact rational surface and the original solver.

use super::{Bounded, MetricSpace};

#[derive(Clone, Debug)]
struct Node {
    sample: usize,
    axis: usize,
    left: Option<usize>,
    right: Option<usize>,
    min: [f64; 4],
    max: [f64; 4],
}

/// Immutable FP64 spatial index. The cached sample-order index is retained,
/// so exact distance ties choose the original first sample.
#[derive(Clone, Debug)]
pub(super) struct SeedSpatialIndex {
    nodes: Vec<Node>,
    root: Option<usize>,
    dimensions: usize,
}

impl SeedSpatialIndex {
    /// Construct a tree for all finite, dimensional sample positions.
    pub(super) fn build<P>(samples: &[(f64, f64, P)]) -> Option<Self>
    where P: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + Copy {
        let dimensions = samples.first()?.2.dimension();
        if dimensions == 0
            || dimensions > 4
            || samples
                .iter()
                .any(|(_, _, p)| (0..dimensions).any(|axis| !p[axis].is_finite()))
        {
            return None;
        }
        let mut indices = (0..samples.len()).collect::<Vec<_>>();
        let mut tree = Self {
            nodes: Vec::with_capacity(samples.len()),
            root: None,
            dimensions,
        };
        tree.root = tree.partition(samples, &mut indices);
        Some(tree)
    }

    fn partition<P>(&mut self, samples: &[(f64, f64, P)], indices: &mut [usize]) -> Option<usize>
    where P: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + Copy {
        if indices.is_empty() {
            return None;
        }
        let mut min = [f64::INFINITY; 4];
        let mut max = [f64::NEG_INFINITY; 4];
        for &index in indices.iter() {
            for axis in 0..self.dimensions {
                let position = samples[index].2[axis];
                min[axis] = min[axis].min(position);
                max[axis] = max[axis].max(position);
            }
        }
        let split_axis = (0..self.dimensions)
            .max_by(|&a, &b| (max[a] - min[a]).total_cmp(&(max[b] - min[b])))?;
        let mid = indices.len() / 2;
        indices.select_nth_unstable_by(mid, |&a, &b| {
            samples[a].2[split_axis]
                .total_cmp(&samples[b].2[split_axis])
                .then_with(|| a.cmp(&b))
        });
        let sample = indices[mid];
        let node_id = self.nodes.len();
        self.nodes.push(Node {
            sample,
            axis: split_axis,
            left: None,
            right: None,
            min,
            max,
        });
        let (left, right_with_mid) = indices.split_at_mut(mid);
        let (_, right) = right_with_mid.split_first_mut()?;
        let left_child = self.partition(samples, left);
        let right_child = self.partition(samples, right);
        self.nodes[node_id].left = left_child;
        self.nodes[node_id].right = right_child;
        Some(node_id)
    }

    /// The original strict-minimum index. `None` means the caller must use
    /// the linear scan (nonfinite coordinate or every distance nonfinite).
    pub(super) fn nearest<P>(&self, samples: &[(f64, f64, P)], target: P) -> Option<usize>
    where P: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + Copy {
        if (0..self.dimensions).any(|axis| !target[axis].is_finite()) {
            return None;
        }
        let mut state = (f64::INFINITY, None);
        if let Some(root) = self.root {
            self.visit(root, samples, target, &mut state);
        }
        state.1
    }

    fn visit<P>(
        &self,
        index: usize,
        samples: &[(f64, f64, P)],
        target: P,
        best: &mut (f64, Option<usize>),
    ) where
        P: Bounded<Scalar = f64> + MetricSpace<Metric = f64> + Copy,
    {
        let node = &self.nodes[index];
        let distance = samples[node.sample].2.distance2(target);
        if distance < best.0 || (distance == best.0 && best.1.is_some_and(|i| node.sample < i)) {
            *best = (distance, Some(node.sample));
        }
        let go_left_first = target[node.axis] < samples[node.sample].2[node.axis];
        let near = if go_left_first { node.left } else { node.right };
        let far = if go_left_first { node.right } else { node.left };
        if let Some(child) = near {
            self.visit(child, samples, target, best);
        }
        if let Some(child) = far {
            // An AABB provides a geometrical lower bound. Pad the bound by
            // floating-point coordinate rounding; extra traversal is safe.
            // If any bound overflows or is not finite, never prune it.
            let bounds = &self.nodes[child];
            let mut lower = 0.0;
            let mut scale = 1.0_f64;
            for axis in 0..self.dimensions {
                let q = target[axis];
                let gap = if q < bounds.min[axis] {
                    bounds.min[axis] - q
                } else if q > bounds.max[axis] {
                    q - bounds.max[axis]
                } else {
                    0.0
                };
                lower += gap * gap;
                scale = scale
                    .max(q.abs())
                    .max(bounds.min[axis].abs())
                    .max(bounds.max[axis].abs());
            }
            let roundoff = 128.0 * f64::EPSILON * scale * scale;
            if !(lower.is_finite() && roundoff.is_finite()) || lower <= best.0 + roundoff {
                self.visit(child, samples, target, best);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::{Point2, Point3};

    fn exhaustive<P: MetricSpace<Metric = f64> + Copy>(
        samples: &[(f64, f64, P)],
        point: P,
    ) -> Option<usize> {
        samples
            .iter()
            .enumerate()
            .fold((f64::INFINITY, None), |best, (i, sample)| {
                let dist = sample.2.distance2(point);
                if dist < best.0 { (dist, Some(i)) } else { best }
            })
            .1
    }

    #[test]
    fn exact_first_minimum_for_curved_grids_and_near_ties_in_two_and_three_dimensions() {
        let sample3 = (0..=50)
            .flat_map(|u| {
                (0..=50).map(move |v| {
                    let angle = 0.12 + 1.5 * f64::from(v) / 50.0;
                    let p = Point3::new(
                        67.1195 + 1.27 * f64::from(u) / 50.0,
                        5.36457 + 0.6731 * angle.cos(),
                        1.67943 + 0.6731 * angle.sin(),
                    );
                    (0.0, 0.0, p)
                })
            })
            .collect::<Vec<_>>();
        let tree = SeedSpatialIndex::build(&sample3).expect("finite three-dimensional grid");
        for index in (0..sample3.len()).step_by(7) {
            let a = sample3[index].2;
            let b = sample3[(index + 51).min(sample3.len() - 1)].2;
            for delta in [-1.0e-10, 0.0, 1.0e-10] {
                let point = Point3::new((a.x + b.x) / 2.0 + delta, a.y, a.z);
                assert_eq!(tree.nearest(&sample3, point), exhaustive(&sample3, point));
            }
        }
        let sample2 = vec![
            (0.0, 0.0, Point2::new(1.0, 0.0)),
            (1.0, 0.0, Point2::new(1.0, 0.0)),
            (0.0, 1.0, Point2::new(2.0, 1.0)),
        ];
        let tree2 = SeedSpatialIndex::build(&sample2).expect("finite two-dimensional grid");
        assert_eq!(tree2.nearest(&sample2, Point2::new(1.0, 0.0)), Some(0));
        assert_eq!(tree2.nearest(&sample2, Point2::new(2.0, 1.0)), Some(2));
        assert_eq!(tree2.nearest(&sample2, Point2::new(f64::NAN, 0.0)), None);
    }
}
