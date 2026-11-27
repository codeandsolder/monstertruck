//! Contract globally necessary matches before bounded geometric search.
use super::super::{EdgeIncidence, connected_membership, forced};
use super::{
    class_graph::ClassGraph,
    search::{Constraint, Solution},
};
use std::collections::BTreeSet;

pub(super) fn necessary_joins(
    constraints: &[Constraint],
    edges: &[[usize; 2]],
    graph: &ClassGraph,
) -> Vec<[usize; 2]> {
    let incidence: EdgeIncidence = edges
        .iter()
        .copied()
        .zip(constraints.iter().map(|c| {
            c.forward
                .iter()
                .map(|&id| (id, true))
                .chain(c.reverse.iter().map(|&id| (id, false)))
                .collect()
        }))
        .collect();
    forced::necessary_joins(
        graph.class_of.len(),
        &incidence,
        Some(forced::GeometryExclusions {
            classes: &graph.class_of,
            excluded: &graph.excluded,
        }),
    )
    .into_iter()
    .map(|pair| pair.polygons)
    .collect()
}

pub(super) fn complete(
    constraints: &[Constraint],
    graph: &ClassGraph,
    joins: &[[usize; 2]],
) -> Option<Solution> {
    let (roots, count) = connected_membership(
        graph.class_of.len(),
        &EdgeIncidence::default(),
        joins.iter().copied(),
    );
    // Every required join is already necessary. The partition is unique
    // only if all exposed uses are now paired within its connected groups.
    for constraint in constraints {
        let forward = constraint
            .forward
            .iter()
            .map(|&id| roots[id])
            .collect::<BTreeSet<_>>();
        let reverse = constraint
            .reverse
            .iter()
            .map(|&id| roots[id])
            .collect::<BTreeSet<_>>();
        if forward.len() != constraint.forward.len()
            || reverse.len() != constraint.reverse.len()
            || forward != reverse
        {
            return None;
        }
    }
    let anchored = graph
        .anchors()
        .map(|id| roots[graph.instances[id][0]])
        .collect::<BTreeSet<_>>();
    if anchored.len() != count {
        return None;
    }
    // Validate all geometric exclusions on the resulting partition.
    let mut classes = vec![BTreeSet::new(); count];
    for (&root, &class) in roots.iter().zip(&graph.class_of) {
        if !classes[root].is_disjoint(&graph.excluded[class]) {
            return None;
        }
        classes[root].insert(class);
    }
    Some(Solution {
        roots,
        states: joins.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::super::{peel, search};
    use super::*;

    #[test]
    fn global_matching_closes_when_anchor_local_propagation_stalls() {
        let constraints = [
            Constraint {
                forward: vec![0, 3],
                reverse: vec![2, 4],
            },
            Constraint {
                forward: vec![4, 0],
                reverse: vec![5, 1],
            },
            Constraint {
                forward: vec![4, 0],
                reverse: vec![3, 2],
            },
        ];
        let classes = [0, 1, 2, 3, 4, 5];
        let graph =
            ClassGraph::new(&constraints, &classes, &[[1, 3], [2, 3]]).expect("valid exclusions");
        assert!(peel::solve(&graph).is_none());
        let joins = necessary_joins(&constraints, &[[0, 1], [1, 2], [2, 3]], &graph);
        let result = complete(&constraints, &graph, &joins).expect("global closure");
        assert_eq!(result.roots, [0, 0, 0, 1, 1, 1]);
        let exhaustive =
            search::solve(&constraints, &classes, &graph, &[]).expect("unique partition");
        for (i, &root) in result.roots.iter().enumerate() {
            for (j, &other) in result.roots.iter().enumerate() {
                assert_eq!(root == other, exhaustive.roots[i] == exhaustive.roots[j]);
            }
        }
    }
}
