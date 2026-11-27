//! Necessary connectivity constraints for independently closed oriented
//! components. A pairing is forced only if every local perfect matching
//! consistent with the already known component incidence contains it.
//! This is a necessary constraint under the closed orientable shell
//! hypothesis, not an assertion that a globally feasible partition exists.

use super::super::super::EdgeIncidence;
use rustc_hash::FxHashSet as HashSet;
use std::collections::{BTreeSet, VecDeque};

/// A necessary join certified by component incidence and perfect matching,
/// accepted only as part of a closed, nonzero-volume source component.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SourceForcedPair {
    /// Canonical edge IDs at the time of the proof.
    pub edge: [usize; 2],
    /// Original polygon IDs joined across the edge.
    pub polygons: [usize; 2],
}

/// Additional necessary incompatibilities between exact geometry classes.
#[derive(Clone, Copy)]
pub(super) struct GeometryExclusions<'a> {
    pub(super) classes: &'a [usize],
    pub(super) excluded: &'a [BTreeSet<usize>],
}

struct ComponentGeometry {
    classes: Vec<HashSet<usize>>,
    excluded: Vec<HashSet<usize>>,
}

impl ComponentGeometry {
    fn new(input: GeometryExclusions<'_>) -> Self {
        Self {
            classes: input
                .classes
                .iter()
                .map(|&class| std::iter::once(class).collect())
                .collect(),
            excluded: input
                .classes
                .iter()
                .map(|&class| input.excluded[class].iter().copied().collect())
                .collect(),
        }
    }

    fn compatible(&self, a: usize, b: usize) -> bool {
        self.classes[a].is_disjoint(&self.excluded[b])
    }

    fn join(&mut self, a: usize, b: usize) {
        let classes = std::mem::take(&mut self.classes[b]);
        self.classes[a].extend(classes);
        let excluded = std::mem::take(&mut self.excluded[b]);
        self.excluded[a].extend(excluded);
    }
}

struct Components {
    parent: Vec<usize>,
    size: Vec<usize>,
    uses: Vec<HashSet<(usize, bool)>>,
    valid: Vec<bool>,
    geometry: Option<ComponentGeometry>,
}

impl Components {
    fn new(count: usize, geometry: Option<GeometryExclusions<'_>>) -> Self {
        Self {
            parent: (0..count).collect(),
            size: vec![1; count],
            uses: vec![HashSet::default(); count],
            valid: vec![true; count],
            geometry: geometry.map(ComponentGeometry::new),
        }
    }

    fn root(&mut self, mut id: usize) -> usize {
        while self.parent[id] != id {
            self.parent[id] = self.parent[self.parent[id]];
            id = self.parent[id];
        }
        id
    }

    fn compatible(&self, a: usize, b: usize) -> bool {
        self.valid[a]
            && self.valid[b]
            && (a == b
                || (self.uses[a].is_disjoint(&self.uses[b])
                    && self.geometry.as_ref().is_none_or(|g| g.compatible(a, b))))
    }

    fn join(&mut self, mut a: usize, mut b: usize) -> usize {
        if a == b {
            return a;
        }
        if self.size[a] < self.size[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.valid[a] &= self.valid[b];
        if let Some(geometry) = &mut self.geometry {
            self.valid[a] &= geometry.compatible(a, b);
            geometry.join(a, b);
        }
        self.parent[b] = a;
        self.size[a] += self.size[b];
        let uses = std::mem::take(&mut self.uses[b]);
        self.uses[a].extend(uses);
        a
    }
}

// Find a perfect matching by augmenting paths. No enumeration of matchings
// and no assumption about the number of coincident edge uses.
fn perfect_matching(adjacency: &[Vec<usize>]) -> Option<(Vec<usize>, Vec<usize>)> {
    let n = adjacency.len();
    let mut left = vec![None; n];
    let mut right = vec![None::<usize>; n];
    for seed in 0..n {
        let mut parent = vec![None; n];
        let mut seen = vec![false; n];
        seen[seed] = true;
        let mut queue = VecDeque::from([seed]);
        let mut free = None;
        'search: while let Some(i) = queue.pop_front() {
            for &j in &adjacency[i] {
                if parent[j].is_some() {
                    continue;
                }
                parent[j] = Some(i);
                if let Some(next) = right[j] {
                    if !seen[next] {
                        seen[next] = true;
                        queue.push_back(next);
                    }
                } else {
                    free = Some(j);
                    break 'search;
                }
            }
        }
        let mut j = free?;
        loop {
            let i = parent[j]?;
            right[j] = Some(i);
            let previous = left[i].replace(j);
            if let Some(previous) = previous {
                j = previous;
            } else {
                break;
            }
        }
    }
    Some((
        left.into_iter().collect::<Option<Vec<_>>>()?,
        right.into_iter().collect::<Option<Vec<_>>>()?,
    ))
}

// Relative to any perfect matching, another matching can omit its i-th
// pair iff an alternating cycle contains that pair. Contract the matched
// edges and search for a nontrivial directed cycle through i.
fn forced_partner(adjacency: &[Vec<usize>], inverse: &[usize], i: usize) -> bool {
    let mut seen = vec![false; adjacency.len()];
    let mut queue = VecDeque::new();
    seen[i] = true;
    for &j in &adjacency[i] {
        let next = inverse[j];
        if next != i && !seen[next] {
            seen[next] = true;
            queue.push_back(next);
        }
    }
    while let Some(current) = queue.pop_front() {
        for &j in &adjacency[current] {
            let next = inverse[j];
            if next == i {
                return false;
            }
            if !seen[next] {
                seen[next] = true;
                queue.push_back(next);
            }
        }
    }
    true
}

fn forced_join(
    components: &mut Components,
    uses: &[(usize, bool)],
) -> Option<([usize; 2], [usize; 2])> {
    let (mut forward, mut reverse): (Vec<_>, Vec<_>) =
        uses.iter().copied().partition(|&(_, forward)| forward);
    if forward.len() != reverse.len() || forward.len() < 2 {
        return None;
    }
    forward.sort_unstable();
    reverse.sort_unstable();
    let left = forward
        .iter()
        .map(|&(p, _)| components.root(p))
        .collect::<Vec<_>>();
    let right = reverse
        .iter()
        .map(|&(p, _)| components.root(p))
        .collect::<Vec<_>>();
    let adjacency = left
        .iter()
        .map(|&a| {
            right
                .iter()
                .enumerate()
                .filter_map(|(j, &b)| components.compatible(a, b).then_some(j))
                .collect()
        })
        .collect::<Vec<Vec<usize>>>();
    let (matching, inverse) = perfect_matching(&adjacency)?;
    for (i, &j) in matching.iter().enumerate() {
        if left[i] != right[j] && forced_partner(&adjacency, &inverse, i) {
            return Some(([left[i], right[j]], [forward[i].0, reverse[j].0]));
        }
    }
    None
}

pub(super) fn necessary_joins(
    count: usize,
    incidence: &EdgeIncidence,
    geometry: Option<GeometryExclusions<'_>>,
) -> Vec<SourceForcedPair> {
    let mut ambiguous = incidence
        .iter()
        .filter(|(_, uses)| !matches!(uses.as_slice(), [a, b] if a.1 != b.1))
        .collect::<Vec<_>>();
    if !ambiguous.iter().any(|(_, uses)| uses.len() > 2) {
        return Vec::new();
    }
    ambiguous.sort_unstable_by_key(|&(&edge, _)| edge);
    let mut components = Components::new(count, geometry);
    for uses in incidence.values() {
        if let [a, b] = uses.as_slice()
            && a.1 != b.1
        {
            let root_a = components.root(a.0);
            let root_b = components.root(b.0);
            components.join(root_a, root_b);
        }
    }
    for (id, (_, uses)) in ambiguous.iter().enumerate() {
        for &(polygon, forward) in *uses {
            let root = components.root(polygon);
            if !components.uses[root].insert((id, forward)) {
                components.valid[root] = false;
            }
        }
    }
    let mut queue = (0..ambiguous.len()).collect::<VecDeque<_>>();
    let mut queued = vec![true; ambiguous.len()];
    let mut result = Vec::new();
    while let Some(id) = queue.pop_front() {
        queued[id] = false;
        let (&edge, uses) = ambiguous[id];
        if let Some(([a, b], polygons)) = forced_join(&mut components, uses) {
            let root = components.join(a, b);
            result.push(SourceForcedPair { edge, polygons });
            // Only edge constraints incident to the changed component can
            // have changed. Reevaluate those, including the current edge.
            let mut affected = components.uses[root]
                .iter()
                .map(|&(id, _)| id)
                .collect::<Vec<_>>();
            affected.sort_unstable();
            for id in affected {
                if !queued[id] {
                    queued[id] = true;
                    queue.push_back(id);
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_certifies_only_necessary_pairs() {
        // Two alternating cycles connected through a mandatory middle
        // match. No row or column has degree one.
        let graph = vec![
            vec![0, 1],
            vec![0, 1, 2],
            vec![2, 3],
            vec![3, 4],
            vec![3, 4],
        ];
        let (matching, inverse) = perfect_matching(&graph).expect("perfect matching");
        assert_eq!(matching[2], 2);
        assert_eq!(
            (0..5)
                .map(|i| forced_partner(&graph, &inverse, i))
                .collect::<Vec<_>>(),
            vec![false, false, true, false, false]
        );
        assert!(perfect_matching(&[vec![0], vec![0]]).is_none());
    }

    #[test]
    fn exhaustive_three_by_three_matching_certificates() {
        // Enumerate every graph of order three and compare against all
        // six permutations. This oracle is test-only.
        let permutations = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        for bits in 0..(1 << 9) {
            let graph = (0..3)
                .map(|i| {
                    (0..3)
                        .filter(|&j| bits & (1 << (3 * i + j)) != 0)
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let feasible = permutations
                .iter()
                .filter(|p| (0..3).all(|i| graph[i].contains(&p[i])))
                .collect::<Vec<_>>();
            match perfect_matching(&graph) {
                Some((matching, inverse)) => {
                    assert!(!feasible.is_empty());
                    for (i, &j) in matching.iter().enumerate() {
                        assert_eq!(
                            forced_partner(&graph, &inverse, i),
                            feasible.iter().all(|p| p[i] == j)
                        );
                    }
                }
                None => assert!(feasible.is_empty()),
            }
        }
    }

    #[test]
    fn symmetric_incidence_remains_ambiguous() {
        let mut incidence = EdgeIncidence::default();
        incidence.insert([0, 1], vec![(0, true), (1, true), (2, false), (3, false)]);
        assert!(necessary_joins(4, &incidence, None).is_empty());
    }

    #[test]
    fn global_component_conflicts_force_joins() {
        let mut incidence = EdgeIncidence::default();
        // Regular edge joins polygons 0/1, 2/3. A separate ambiguous
        // edge excludes cross-component pairings of 0 with 3 and 2 with 1.
        incidence.insert([0, 1], vec![(0, true), (1, false)]);
        incidence.insert([2, 3], vec![(2, true), (3, false)]);
        incidence.insert([4, 5], vec![(0, true), (2, false), (4, true), (5, false)]);
        incidence.insert([6, 7], vec![(1, true), (3, true), (4, false), (5, false)]);
        // Components 0/1 and 2/3 share a forward use at [6,7];
        // component 4 cannot join 5 because both use [6,7] in reverse.
        // Thus the groups {0,1,5} and {2,3,4} are forced.
        let result = necessary_joins(6, &incidence, None);
        assert_eq!(result.len(), 2);
        assert!(result.iter().any(|p| p.polygons == [0, 5]));
        assert!(
            result
                .iter()
                .any(|p| p.polygons == [4, 2] || p.polygons == [3, 4])
        );
    }
}
