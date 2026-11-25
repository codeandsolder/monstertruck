//! Immutable, dimension-independent bounding hierarchy. Build by partitioning
//! contiguous compact bounds/centroid records, then gather the geometry once
//! into contiguous leaf order.
//! Callers own all narrow-phase semantics and nonfinite-input handling.

use std::ops::Range;

#[derive(Clone, Copy, Debug)]
pub struct Bounds<const N: usize> {
    pub(crate) low: [f64; N],
    pub(crate) high: [f64; N],
}

impl<const N: usize> Bounds<N> {
    pub(crate) fn from_points<const M: usize>(points: [[f64; N]; M]) -> Self {
        let mut low = [f64::INFINITY; N];
        let mut high = [f64::NEG_INFINITY; N];
        for point in points {
            for axis in 0..N {
                low[axis] = low[axis].min(point[axis]);
                high[axis] = high[axis].max(point[axis]);
            }
        }
        Self { low, high }
    }

    fn union(self, other: Self) -> Self {
        Self {
            low: std::array::from_fn(|axis| self.low[axis].min(other.low[axis])),
            high: std::array::from_fn(|axis| self.high[axis].max(other.high[axis])),
        }
    }

    fn center(self) -> [f64; N] {
        std::array::from_fn(|axis| self.low[axis].midpoint(self.high[axis]))
    }

    fn longest_axis(self) -> usize {
        let span = |axis: usize| self.high[axis] - self.low[axis];
        (1..N).fold(0, |best, axis| {
            if span(axis).total_cmp(&span(best)).is_gt() {
                axis
            } else {
                best
            }
        })
    }
}

#[derive(Clone, Debug)]
enum Children {
    Leaf(Range<usize>),
    Branch([usize; 2]),
}

#[derive(Clone, Debug)]
struct Node<const N: usize> {
    bounds: Bounds<N>,
    children: Children,
}

#[derive(Clone, Debug)]
pub struct SpatialHierarchy<const N: usize> {
    nodes: Vec<Node<N>>,
}

const LEAF_SIZE: usize = 8;

/// Keep bounds, sort keys and source index adjacent during partitioning.
/// Indirectly gathering bounds and centroids through a permuted index list
/// causes random memory accesses at every recursive level.
#[derive(Clone, Copy)]
struct BuildRecord<const N: usize> {
    bounds: Bounds<N>,
    center: [f64; N],
    index: usize,
}

fn build<const N: usize>(
    records: &mut [BuildRecord<N>],
    offset: usize,
    nodes: &mut Vec<Node<N>>,
) -> usize {
    let bounds = records[1..]
        .iter()
        .fold(records[0].bounds, |bounds, record| {
            bounds.union(record.bounds)
        });
    let children = if records.len() <= LEAF_SIZE {
        Children::Leaf(offset..offset + records.len())
    } else {
        let axis = bounds.longest_axis();
        let middle = records.len() / 2;
        records.select_nth_unstable_by(middle, |a, b| {
            a.center[axis]
                .total_cmp(&b.center[axis])
                .then_with(|| a.index.cmp(&b.index))
        });
        let (left, right) = records.split_at_mut(middle);
        Children::Branch([
            build(left, offset, nodes),
            build(right, offset + middle, nodes),
        ])
    };
    let id = nodes.len();
    nodes.push(Node { bounds, children });
    id
}

impl<const N: usize> SpatialHierarchy<N> {
    /// Rearrange entries into contiguous leaves. No geometry is moved during
    /// recursive median selection; compact records retain bounds and
    /// precomputed centroid keys in contiguous storage.
    pub(crate) fn new<T: Copy>(entries: &mut Vec<T>, bounds: impl Fn(T) -> Bounds<N>) -> Self {
        if entries.is_empty() {
            return Self { nodes: Vec::new() };
        }
        if entries.len() <= LEAF_SIZE {
            let node = Node {
                bounds: entries[1..]
                    .iter()
                    .copied()
                    .fold(bounds(entries[0]), |acc, entry| acc.union(bounds(entry))),
                children: Children::Leaf(0..entries.len()),
            };
            return Self { nodes: vec![node] };
        }
        let mut records = entries
            .iter()
            .copied()
            .enumerate()
            .map(|(index, entry)| {
                let bounds = bounds(entry);
                BuildRecord {
                    bounds,
                    center: bounds.center(),
                    index,
                }
            })
            .collect::<Vec<_>>();
        let mut nodes = Vec::new();
        build(&mut records, 0, &mut nodes);
        *entries = records
            .into_iter()
            .map(|record| entries[record.index])
            .collect();
        Self { nodes }
    }

    fn visit(
        &self,
        index: usize,
        overlaps: &impl Fn(Bounds<N>) -> bool,
        predicate: &mut impl FnMut(Range<usize>) -> bool,
    ) -> bool {
        let node = &self.nodes[index];
        if !overlaps(node.bounds) {
            return true;
        }
        match &node.children {
            Children::Leaf(range) => predicate(range.clone()),
            Children::Branch(children) => children
                .iter()
                .all(|&child| self.visit(child, overlaps, predicate)),
        }
    }

    /// Visit all potentially overlapping contiguous leaf ranges; false
    /// from the caller terminates the query immediately.
    pub(crate) fn all_matching(
        &self,
        overlaps: impl Fn(Bounds<N>) -> bool,
        mut predicate: impl FnMut(Range<usize>) -> bool,
    ) -> bool {
        self.nodes
            .len()
            .checked_sub(1)
            .is_none_or(|root| self.visit(root, &overlaps, &mut predicate))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_entry_has_exactly_one_leaf_and_pruning_preserves_overlaps() {
        let mut entries = (0_u32..257)
            .map(|i| {
                let x = f64::from(i % 17);
                let y = f64::from(i / 17);
                (i, Bounds::from_points([[x, y], [x + 0.3, y + 0.4]]))
            })
            .collect::<Vec<_>>();
        let tree = SpatialHierarchy::new(&mut entries, |entry| entry.1);
        let mut visited = vec![0; entries.len()];
        assert!(tree.all_matching(
            |_| true,
            |range| {
                for i in range {
                    visited[entries[i].0 as usize] += 1;
                }
                true
            }
        ));
        assert!(visited.iter().all(|&count| count == 1));
        for x in 0..18 {
            for y in 0..17 {
                let point = [f64::from(x) + 0.1, f64::from(y) + 0.2];
                let overlap = |bounds: Bounds<2>| {
                    (0..2).all(|axis| {
                        bounds.low[axis] <= point[axis] && point[axis] <= bounds.high[axis]
                    })
                };
                let expected = entries.iter().filter(|entry| overlap(entry.1)).count();
                let mut actual = 0;
                assert!(tree.all_matching(overlap, |range| {
                    actual += entries[range]
                        .iter()
                        .filter(|entry| overlap(entry.1))
                        .count();
                    true
                }));
                assert_eq!(actual, expected);
            }
        }
        assert!(!tree.all_matching(|_| true, |_| false));
    }
}
