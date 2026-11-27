//! Incidence quotient over interchangeable source fragments. A class
//! represents several distinct copies of the same oriented geometry and
//! has the identical set of oriented uses at every ambiguous edge.
use super::search::Constraint;
use std::collections::{BTreeMap, BTreeSet};

pub(super) struct ClassGraph {
    pub(super) labels: Vec<usize>,
    pub(super) class_of: Vec<usize>,
    pub(super) instances: Vec<Vec<usize>>,
    pub(super) edges: Vec<[Vec<usize>; 2]>,
    pub(super) profiles: Vec<Vec<(usize, usize)>>,
    pub(super) excluded: Vec<BTreeSet<usize>>,
}

impl ClassGraph {
    pub(super) fn new(
        constraints: &[Constraint],
        classes: &[usize],
        exclusions: &[[usize; 2]],
    ) -> Option<Self> {
        let mut grouped = BTreeMap::<usize, Vec<usize>>::new();
        for (id, &class) in classes.iter().enumerate() {
            grouped.entry(class).or_default().push(id);
        }
        let (labels, instances): (Vec<_>, Vec<_>) = grouped.into_iter().unzip();
        let mut class_of = vec![0; classes.len()];
        for (class, ids) in instances.iter().enumerate() {
            for &id in ids {
                class_of[id] = class;
            }
        }
        let mut profiles = vec![Vec::new(); instances.len()];
        let mut edges = Vec::with_capacity(constraints.len());
        for (edge, constraint) in constraints.iter().enumerate() {
            if constraint.forward.len() != constraint.reverse.len() || constraint.forward.len() < 2
            {
                return None;
            }
            let mut sides = [Vec::new(), Vec::new()];
            for (direction, side) in [&constraint.forward, &constraint.reverse]
                .into_iter()
                .enumerate()
            {
                let mut seen = BTreeSet::new();
                let mut counts = BTreeMap::<usize, usize>::new();
                for &fragment in side {
                    if !seen.insert(fragment) {
                        return None;
                    }
                    *counts.entry(*class_of.get(fragment)?).or_default() += 1;
                }
                for (class, count) in counts {
                    // No class may be only partly present on a side.
                    if count != instances[class].len() {
                        return None;
                    }
                    profiles[class].push((edge, direction));
                    sides[direction].push(class);
                }
            }
            edges.push(sides);
        }
        let mut excluded = vec![BTreeSet::new(); instances.len()];
        for &[a, b] in exclusions {
            let a = labels.binary_search(&a).ok()?;
            let b = labels.binary_search(&b).ok()?;
            if a == b {
                return None;
            }
            excluded[a].insert(b);
            excluded[b].insert(a);
        }
        Some(Self {
            labels,
            class_of,
            instances,
            edges,
            profiles,
            excluded,
        })
    }

    pub(super) fn anchors(&self) -> impl Iterator<Item = usize> + '_ {
        self.instances
            .iter()
            .enumerate()
            .filter_map(|(id, instances)| (instances.len() == 1).then_some(id))
    }
}
