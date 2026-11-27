//! Bounded incidence partition search modulo exact fragment equivalence.
use super::class_graph::ClassGraph;
use std::collections::{BTreeMap, BTreeSet};

// The bitset represents a local ambiguous cluster, never the whole mesh.
// Exhaustion is inconclusive; no partial or first-found partition is accepted.
pub(super) const MAX_FRAGMENTS: usize = u128::BITS as usize;
const MAX_STATES: usize = 8192;
const MAX_WORK: usize = 1 << 20;

pub(super) struct Constraint {
    pub(super) forward: Vec<usize>,
    pub(super) reverse: Vec<usize>,
}

#[derive(Clone)]
struct State {
    roots: Vec<usize>,
    members: Vec<u128>,
    conflicts: Vec<u128>,
}

impl State {
    fn new(count: usize, constraints: &[Constraint]) -> Option<Self> {
        if count == 0 || count > MAX_FRAGMENTS {
            return None;
        }
        let members = (0..count).map(|id| 1u128 << id).collect::<Vec<_>>();
        let mut conflicts = vec![0; count];
        for constraint in constraints {
            if constraint.forward.len() != constraint.reverse.len() || constraint.forward.len() < 2
            {
                return None;
            }
            for side in [&constraint.forward, &constraint.reverse] {
                let mut mask = 0;
                for &id in side {
                    let bit = *members.get(id)?;
                    if mask & bit != 0 {
                        return None;
                    }
                    mask |= bit;
                }
                for &id in side {
                    conflicts[id] |= mask & !members[id];
                }
            }
        }
        Some(Self {
            roots: (0..count).collect(),
            members,
            conflicts,
        })
    }

    fn compatible(&self, a: usize, b: usize) -> bool { self.conflicts[a] & self.members[b] == 0 }

    fn join(&self, [a, b]: [usize; 2]) -> Self {
        let mut state = self.clone();
        let (keep, remove) = (a.min(b), a.max(b));
        let mut members = state.members[remove];
        while members != 0 {
            let id = members.trailing_zeros() as usize;
            state.roots[id] = keep;
            members &= members - 1;
        }
        state.members[keep] |= state.members[remove];
        state.conflicts[keep] |= state.conflicts[remove];
        state.members[remove] = 0;
        state.conflicts[remove] = 0;
        state
    }

    fn signature(&self, classes: &[usize]) -> Vec<Vec<usize>> {
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for (&root, &class) in self.roots.iter().zip(classes) {
            groups.entry(root).or_default().push(class);
        }
        let mut signature = groups
            .into_values()
            .map(|mut ids| {
                ids.sort_unstable();
                ids
            })
            .collect::<Vec<_>>();
        signature.sort_unstable();
        signature
    }
}

enum Branch {
    Complete,
    Dead,
    Joins(Vec<[usize; 2]>),
}

struct Search<'a> {
    constraints: &'a [Constraint],
    classes: &'a [usize],
    unique_classes: BTreeSet<usize>,
    visited: BTreeSet<Vec<Vec<usize>>>,
    solution: Option<(Vec<usize>, Vec<Vec<usize>>)>,
    remaining: usize,
}

impl Search<'_> {
    fn branch(&mut self, state: &State) -> Option<Branch> {
        let mut best = None::<Vec<[usize; 2]>>;
        for constraint in self.constraints {
            for (side, other) in [
                (&constraint.forward, &constraint.reverse),
                (&constraint.reverse, &constraint.forward),
            ] {
                for &id in side {
                    let root = state.roots[id];
                    self.remaining = self.remaining.checked_sub(other.len())?;
                    if other.iter().any(|&opposite| state.roots[opposite] == root) {
                        continue;
                    }
                    let joins = other
                        .iter()
                        .filter_map(|&opposite| {
                            let target = state.roots[opposite];
                            state.compatible(root, target).then_some([root, target])
                        })
                        .collect::<Vec<_>>();
                    if joins.is_empty() {
                        return Some(Branch::Dead);
                    }
                    if best.as_ref().is_none_or(|old| joins.len() < old.len()) {
                        best = Some(joins);
                    }
                }
            }
        }
        Some(best.map_or(Branch::Complete, Branch::Joins))
    }

    fn visit(&mut self, state: &State) -> Option<()> {
        self.remaining = self.remaining.checked_sub(self.classes.len())?;
        let signature = state.signature(self.classes);
        if self.visited.contains(&signature) {
            return Some(());
        }
        if self.visited.len() >= MAX_STATES {
            return None;
        }
        self.visited.insert(signature.clone());
        match self.branch(state)? {
            Branch::Dead => Some(()),
            Branch::Joins(joins) => {
                for pair in joins {
                    self.visit(&state.join(pair))?;
                }
                Some(())
            }
            Branch::Complete => {
                // A completely duplicated shell has no independent source
                // witness for its body multiplicity. Distinct noninterchangeable
                // fragments must anchor every resulting shell.
                if signature
                    .iter()
                    .any(|group| !group.iter().any(|id| self.unique_classes.contains(id)))
                {
                    return None;
                }
                if let Some((_, old)) = &self.solution {
                    if *old != signature {
                        return None;
                    }
                } else {
                    self.solution = Some((state.roots.clone(), signature));
                }
                Some(())
            }
        }
    }
}

pub(super) struct Solution {
    pub(super) roots: Vec<usize>,
    pub(super) states: usize,
}

pub(super) fn solve(constraints: &[Constraint], classes: &[usize]) -> Option<Solution> {
    let state = State::new(classes.len(), constraints)?;
    let graph = ClassGraph::new(constraints, classes)?;
    let unique_classes = graph.anchors().map(|id| graph.labels[id]).collect();
    let mut search = Search {
        constraints,
        classes,
        unique_classes,
        visited: BTreeSet::new(),
        solution: None,
        remaining: MAX_WORK,
    };
    search.visit(&state)?;
    Some(Solution {
        roots: search.solution?.0,
        states: search.visited.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn copies_are_interchangeable_but_distinct_caps_are_not() {
        let edge = Constraint {
            forward: vec![0, 1],
            reverse: vec![2, 3],
        };
        let result = solve(&[edge], &[0, 1, 2, 2]).expect("one geometric partition");
        assert_ne!(result.roots[0], result.roots[1]);
        assert_ne!(result.roots[2], result.roots[3]);
        assert_eq!(result.states, 3);
        let edge = Constraint {
            forward: vec![0, 1, 2],
            reverse: vec![3, 4, 5],
        };
        assert!(solve(&[edge], &[0, 1, 2, 3, 4, 4]).is_none());
        let edge = Constraint {
            forward: vec![0, 1],
            reverse: vec![2, 3],
        };
        assert!(solve(&[edge], &[0, 0, 1, 1]).is_none());
    }

    #[test]
    fn incompatible_global_pairing_rejects_even_with_balanced_local_edges() {
        let edges = [
            Constraint {
                forward: vec![0, 1],
                reverse: vec![2, 3],
            },
            Constraint {
                forward: vec![0, 2],
                reverse: vec![1, 3],
            },
            Constraint {
                forward: vec![0, 3],
                reverse: vec![1, 2],
            },
            Constraint {
                forward: vec![0, 4],
                reverse: vec![5, 6],
            },
        ];
        // The first four fragments cannot join any other of those four
        // without taking two same-side uses at some edge. The duplicate
        // pair 5/6 has an identical profile but cannot resolve that conflict.
        assert!(solve(&edges, &[0, 1, 2, 3, 4, 5, 5]).is_none());
    }

    #[test]
    fn global_constraints_can_determine_unique_partition_without_copies() {
        // Both local matchings of every edge are initially feasible.
        // Nevertheless only one complete partition satisfies all four
        // edges. Local necessary-join propagation alone cannot start.
        let edges = [
            Constraint {
                forward: vec![2, 4],
                reverse: vec![3, 6],
            },
            Constraint {
                forward: vec![0, 2],
                reverse: vec![1, 6],
            },
            Constraint {
                forward: vec![1, 6],
                reverse: vec![2, 5],
            },
            Constraint {
                forward: vec![1, 3],
                reverse: vec![4, 5],
            },
        ];
        let result = solve(&edges, &[0, 1, 2, 3, 4, 5, 6]).expect("unique global partition");
        assert_eq!(
            oracle_signature(&result.roots, &[0, 1, 2, 3, 4, 5, 6]),
            vec![vec![0, 1, 5], vec![2, 6], vec![3, 4]]
        );
    }

    fn partitions(count: usize) -> Vec<Vec<usize>> {
        let mut result = vec![vec![0]];
        for _ in 1..count {
            result = result
                .into_iter()
                .flat_map(|group| {
                    let next = group.iter().copied().max().unwrap_or_default() + 1;
                    (0..=next).map(move |label| {
                        let mut extended = group.clone();
                        extended.push(label);
                        extended
                    })
                })
                .collect();
        }
        result
    }

    // Independent exhaustive oracle: partitions rather than sequences
    // of joins. Reject disconnected groups even when their edge counts
    // would otherwise happen to balance.
    fn oracle_partition(partition: &[usize], constraints: &[Constraint]) -> bool {
        let mut adjacency = vec![Vec::new(); partition.len()];
        for constraint in constraints {
            let mut sides = BTreeMap::<usize, [Vec<usize>; 2]>::new();
            for (direction, side) in [&constraint.forward, &constraint.reverse]
                .into_iter()
                .enumerate()
            {
                for &id in side {
                    sides.entry(partition[id]).or_default()[direction].push(id);
                }
            }
            for [forward, reverse] in sides.into_values() {
                let ([a], [b]) = (forward.as_slice(), reverse.as_slice()) else {
                    return false;
                };
                adjacency[*a].push(*b);
                adjacency[*b].push(*a);
            }
        }
        let mut seen = BTreeSet::new();
        let mut groups = BTreeSet::new();
        for (seed, &part) in partition.iter().enumerate() {
            if seen.contains(&seed) {
                continue;
            }
            if !groups.insert(part) {
                return false;
            }
            let mut stack = vec![seed];
            while let Some(id) = stack.pop() {
                if seen.insert(id) {
                    stack.extend(adjacency[id].iter().copied());
                }
            }
        }
        true
    }

    fn oracle_signature(partition: &[usize], classes: &[usize]) -> Vec<Vec<usize>> {
        let mut groups = BTreeMap::<usize, Vec<usize>>::new();
        for (&root, &class) in partition.iter().zip(classes) {
            groups.entry(root).or_default().push(class);
        }
        let mut result = groups
            .into_values()
            .map(|mut group| {
                group.sort_unstable();
                group
            })
            .collect::<Vec<_>>();
        result.sort_unstable();
        result
    }

    #[test]
    fn quotient_search_agrees_with_exhaustive_small_partitions() {
        let count = 5;
        let candidates = (0..1usize << count)
            .filter(|mask| mask.count_ones() == 2)
            .collect::<Vec<_>>();
        let edges = candidates
            .iter()
            .flat_map(|&forward| {
                candidates
                    .iter()
                    .filter(move |&&reverse| forward & reverse == 0 && forward < reverse)
                    .map(move |&reverse| Constraint {
                        forward: (0..count).filter(|&id| forward & (1 << id) != 0).collect(),
                        reverse: (0..count).filter(|&id| reverse & (1 << id) != 0).collect(),
                    })
            })
            .collect::<Vec<_>>();
        let partitions = partitions(count);
        let mut compared = 0;
        for i in 0..edges.len() {
            for j in i + 1..edges.len() {
                let constraints = [
                    Constraint {
                        forward: edges[i].forward.clone(),
                        reverse: edges[i].reverse.clone(),
                    },
                    Constraint {
                        forward: edges[j].forward.clone(),
                        reverse: edges[j].reverse.clone(),
                    },
                ];
                let profiles = (0..count)
                    .map(|id| {
                        constraints
                            .iter()
                            .map(|e| (e.forward.contains(&id), e.reverse.contains(&id)))
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
                if profiles.iter().any(|p| p.iter().all(|&(f, r)| !f && !r)) {
                    continue;
                }
                // Every repeated class is a true incidence symmetry.
                let mut class_ids = BTreeMap::new();
                let classes = profiles
                    .iter()
                    .map(|p| {
                        let next = class_ids.len();
                        *class_ids.entry(p).or_insert(next)
                    })
                    .collect::<Vec<_>>();
                let unique = classes
                    .iter()
                    .copied()
                    .filter(|id| classes.iter().filter(|&other| other == id).count() == 1)
                    .collect::<BTreeSet<_>>();
                let solutions = partitions
                    .iter()
                    .filter(|p| oracle_partition(p, &constraints))
                    .map(|p| oracle_signature(p, &classes))
                    .collect::<BTreeSet<_>>();
                let expected = if solutions.len() == 1 {
                    solutions
                        .first()
                        .filter(|s| s.iter().all(|g| g.iter().any(|id| unique.contains(id))))
                } else {
                    None
                };
                let actual =
                    solve(&constraints, &classes).map(|s| oracle_signature(&s.roots, &classes));
                assert_eq!(actual.as_ref(), expected);
                compared += 1;
            }
        }
        assert!(compared > 10);
    }
}
