//! The edge-incidence proof alone does not exclude a pinched vertex.
//! Every selected source component must also have one circular link at
//! each vertex, before it is separated from the unresolved source pool.
use super::InputPolygon;
use std::collections::{BTreeMap, BTreeSet};

fn cycle(edges: &[[usize; 2]]) -> bool {
    let mut neighbors = BTreeMap::<usize, Vec<usize>>::new();
    for &[a, b] in edges {
        if a == b {
            return false;
        }
        neighbors.entry(a).or_default().push(b);
        neighbors.entry(b).or_default().push(a);
    }
    if neighbors.values().any(|n| n.len() != 2) {
        return false;
    }
    let Some((&start, _)) = neighbors.first_key_value() else {
        return false;
    };
    let mut seen = BTreeSet::new();
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        if seen.insert(node) {
            stack.extend(neighbors[&node].iter().copied());
        }
    }
    seen.len() == neighbors.len()
}

pub(super) fn manifold(
    polygons: &[InputPolygon],
    components: &[Option<usize>],
    ids: impl IntoIterator<Item = usize>,
) -> bool {
    let mut links = BTreeMap::<(usize, usize), Vec<[usize; 2]>>::new();
    for id in ids {
        let Some(component) = components[id] else {
            return false;
        };
        let vertices = &polygons[id].vertices;
        let count = vertices.len();
        if count < 3 || vertices.iter().copied().collect::<BTreeSet<_>>().len() != count {
            return false;
        }
        for (i, &vertex) in vertices.iter().enumerate() {
            links
                .entry((component, vertex))
                .or_default()
                .push([vertices[(i + count - 1) % count], vertices[(i + 1) % count]]);
        }
    }
    links.values().all(|edges| cycle(edges))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn link_must_be_one_cycle_not_two_touching_components() {
        assert!(cycle(&[[0, 1], [1, 2], [2, 0]]));
        assert!(!cycle(&[[0, 1], [1, 2], [2, 0], [3, 4], [4, 5], [5, 3]]));
        assert!(!cycle(&[[0, 1], [1, 2]]));
        assert!(cycle(&[[0, 1], [0, 1]]));
    }
}
