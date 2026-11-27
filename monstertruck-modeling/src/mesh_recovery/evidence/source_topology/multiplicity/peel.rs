//! Prove anchored closed shells by necessary joins on exact geometry
//! classes. Unlike instance-level matching, one forced opposite class can
//! have many interchangeable copies. Each closed anchored chain is
//! mandatory in every feasible partition. If the chains consume exactly
//! the supplied multiplicities, their geometric partition is unique.
use super::{class_graph::ClassGraph, search::Solution};
use std::collections::BTreeSet;

// Bound the total candidate/profile checks even for highly ambiguous input.
const MAX_WORK: usize = 1 << 22;

struct Closure<'a> {
    graph: &'a ClassGraph,
    chosen: Vec<usize>,
    used: Vec<[Option<usize>; 2]>,
    pending: BTreeSet<(usize, usize)>,
    forbidden: Vec<bool>,
}

enum ClosureResult {
    Closed(Vec<usize>),
    Unresolved,
}

enum Choice {
    Forced(usize),
    Ambiguous,
    Impossible,
}

impl<'a> Closure<'a> {
    fn new(graph: &'a ClassGraph) -> Self {
        Self {
            graph,
            chosen: Vec::new(),
            used: vec![[None, None]; graph.edges.len()],
            pending: BTreeSet::new(),
            forbidden: vec![false; graph.instances.len()],
        }
    }

    fn add(&mut self, class: usize, remaining: &mut usize) -> Option<bool> {
        let profile = &self.graph.profiles[class];
        *remaining = remaining.checked_sub(
            profile
                .len()
                .checked_add(self.graph.excluded[class].len())?,
        )?;
        if self.forbidden[class]
            || profile
                .iter()
                .any(|&(edge, side)| self.used[edge][side].is_some())
        {
            return Some(false);
        }
        self.chosen.push(class);
        for &other in &self.graph.excluded[class] {
            self.forbidden[other] = true;
        }
        for &(edge, side) in profile {
            self.used[edge][side] = Some(class);
            self.pending.remove(&(edge, side));
            if self.used[edge][1 - side].is_none() {
                self.pending.insert((edge, 1 - side));
            }
        }
        Some(true)
    }

    fn choice(&self, edge: usize, side: usize, remaining: &mut usize) -> Option<Choice> {
        let mut compatible = None;
        for &class in &self.graph.edges[edge][side] {
            let profile = &self.graph.profiles[class];
            *remaining = remaining.checked_sub(profile.len() + 1)?;
            if self.forbidden[class] || profile.iter().any(|&(e, s)| self.used[e][s].is_some()) {
                continue;
            }
            if compatible.replace(class).is_some() {
                return Some(Choice::Ambiguous);
            }
        }
        Some(compatible.map_or(Choice::Impossible, Choice::Forced))
    }

    fn close(
        mut self,
        anchor: usize,
        remaining: &mut usize,
        states: &mut usize,
    ) -> Option<ClosureResult> {
        if !self.add(anchor, remaining)? {
            return Some(ClosureResult::Unresolved);
        }
        *states += 1;
        while !self.pending.is_empty() {
            let mut forced = None;
            for &(edge, side) in &self.pending {
                match self.choice(edge, side, remaining)? {
                    Choice::Impossible => return Some(ClosureResult::Unresolved),
                    Choice::Ambiguous => {}
                    Choice::Forced(class) => {
                        forced = Some(class);
                        break;
                    }
                }
            }
            let Some(class) = forced else {
                return Some(ClosureResult::Unresolved);
            };
            if !self.add(class, remaining)? {
                return Some(ClosureResult::Unresolved);
            }
            *states += 1;
        }
        self.chosen.sort_unstable();
        Some(ClosureResult::Closed(self.chosen))
    }
}

fn allocate_instances(graph: &ClassGraph, shells: &[Vec<usize>]) -> Option<Vec<usize>> {
    let total = graph.instances.iter().map(Vec::len).sum();
    let mut roots = vec![usize::MAX; total];
    let mut assigned = vec![0; graph.instances.len()];
    for shell in shells {
        let mut ids = Vec::with_capacity(shell.len());
        for &class in shell {
            ids.push(*graph.instances[class].get(assigned[class])?);
            assigned[class] += 1;
        }
        let root = *ids.iter().min()?;
        for id in ids {
            roots[id] = root;
        }
    }
    roots.iter().all(|&id| id != usize::MAX).then_some(roots)
}

pub(super) fn solve(graph: &ClassGraph) -> Option<Solution> {
    let mut remaining = MAX_WORK;
    let mut states = 0;
    let mut consumed = vec![0; graph.instances.len()];
    let mut shells = Vec::new();
    for anchor in graph.anchors() {
        if consumed[anchor] != 0 {
            continue;
        }
        let ClosureResult::Closed(shell) =
            Closure::new(graph).close(anchor, &mut remaining, &mut states)?
        else {
            continue;
        };
        for &class in &shell {
            consumed[class] += 1;
            if consumed[class] > graph.instances[class].len() {
                return None;
            }
        }
        shells.push(shell);
    }
    if consumed
        .iter()
        .zip(&graph.instances)
        .any(|(&n, ids)| n != ids.len())
    {
        return None;
    }
    Some(Solution {
        roots: allocate_instances(graph, &shells)?,
        states,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh_recovery::evidence::source_topology::multiplicity::search::Constraint;

    #[test]
    fn forced_class_closure_scales_past_instance_search_limit() {
        let count = 80;
        let constraint = Constraint {
            forward: (0..count).collect(),
            reverse: (count..count * 2).collect(),
        };
        let classes = (0..count)
            .chain(std::iter::repeat_n(count, count))
            .collect::<Vec<_>>();
        let graph = ClassGraph::new(&[constraint], &classes, &[]).expect("valid class graph");
        let result = solve(&graph).expect("all anchored shells forced");
        assert_eq!(result.states, count * 2);
        for id in 0..count {
            assert_eq!(result.roots[id], id);
            assert_eq!(result.roots[id + count], id);
        }
    }

    #[test]
    fn opposed_face_witnesses_make_ambiguous_closure_necessary() {
        let constraints = [Constraint {
            forward: vec![0, 1],
            reverse: vec![2, 3],
        }];
        let classes = [10, 20, 30, 40];
        let initial = ClassGraph::new(&constraints, &classes, &[]).expect("valid class graph");
        assert!(solve(&initial).is_none());
        let graph = ClassGraph::new(&constraints, &classes, &[[10, 30]]).expect("valid exclusions");
        let solution = solve(&graph).expect("necessary compatible partition");
        assert_eq!(solution.roots, [0, 1, 1, 0]);
        let searched = super::super::search::solve(&constraints, &classes, &graph, &[])
            .expect("same exhaustive partition");
        assert_eq!(searched.roots, solution.roots);
        let impossible = ClassGraph::new(&constraints, &classes, &[[10, 30], [10, 40]])
            .expect("valid exclusions");
        assert!(solve(&impossible).is_none());
        assert!(super::super::search::solve(&constraints, &classes, &impossible, &[]).is_none());
    }

    #[test]
    fn forced_closure_never_invents_anchors_or_copies() {
        let constraint = Constraint {
            forward: vec![0, 1],
            reverse: vec![2, 3],
        };
        for classes in [[0, 0, 1, 1], [0, 1, 2, 3]] {
            let graph = ClassGraph::new(std::slice::from_ref(&constraint), &classes, &[])
                .expect("valid class graph");
            assert!(solve(&graph).is_none());
        }
    }
}
