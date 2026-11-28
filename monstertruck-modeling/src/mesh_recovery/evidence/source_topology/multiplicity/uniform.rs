//! Preserve the integral multiplicity of a closed oriented geometric chain.
//! Each exposed edge has precisely one interchangeable class in each direction.
//! The incidence quotient determines every join, up to permutation of copies.
//! This is evidence for repeated source coverage, not for distinct material
//! bodies or the absence of accidental duplicate input records.

use super::{class_graph::ClassGraph, search::Solution};
use crate::mesh_recovery::evidence::source_topology::{EdgeIncidence, connected_membership};

pub(super) fn solve(graph: &ClassGraph) -> Option<Solution> {
    let joins = graph
        .edges
        .iter()
        .map(|[a, b]| {
            let ([a], [b]) = (a.as_slice(), b.as_slice()) else {
                return None;
            };
            Some([*a, *b])
        })
        .collect::<Option<Vec<_>>>()?;
    let (membership, count) =
        connected_membership(graph.instances.len(), &EdgeIncidence::default(), joins);
    let mut groups = vec![Vec::new(); count];
    for (class, group) in membership.iter().copied().enumerate() {
        groups[group].push(class);
    }
    let mut roots = vec![usize::MAX; graph.class_of.len()];
    for (group_id, group) in groups.into_iter().enumerate() {
        let multiplicity = graph.instances[*group.first()?].len();
        if multiplicity < 2
            || group.iter().any(|&class| {
                graph.instances[class].len() != multiplicity
                    || graph.excluded[class]
                        .iter()
                        .any(|&other| membership[other] == group_id)
            })
        {
            return None;
        }
        // Every instance is interchangeable with every other member of its
        // class. Choose a deterministic bijection across necessary joins;
        // the resulting oriented chain is independent of that permutation.
        for copy in 0..multiplicity {
            let root = group
                .iter()
                .map(|&class| graph.instances[class][copy])
                .min()?;
            for &class in &group {
                roots[graph.instances[class][copy]] = root;
            }
        }
    }
    roots
        .iter()
        .all(|&id| id != usize::MAX)
        .then_some(Solution {
            roots,
            states: graph.edges.len(),
        })
}

#[cfg(test)]
mod tests {
    use super::super::search::Constraint;
    use super::*;

    #[test]
    fn uniform_coverage_has_a_unique_geometric_quotient() -> Result<(), &'static str> {
        let constraints = [
            Constraint {
                forward: vec![0, 1, 2],
                reverse: vec![3, 4, 5],
            },
            Constraint {
                forward: vec![3, 4, 5],
                reverse: vec![6, 7, 8],
            },
            Constraint {
                forward: vec![6, 7, 8],
                reverse: vec![0, 1, 2],
            },
        ];
        let classes = [10, 10, 10, 20, 20, 20, 30, 30, 30];
        let graph = ClassGraph::new(&constraints, &classes, &[]).ok_or("class graph")?;
        assert_eq!(
            solve(&graph).ok_or("uniform quotient")?.roots,
            [0, 1, 2, 0, 1, 2, 0, 1, 2]
        );
        let graph = ClassGraph::new(&constraints, &[10, 10, 11, 20, 20, 20, 30, 30, 30], &[])
            .ok_or("ambiguous graph")?;
        assert!(solve(&graph).is_none());
        Ok(())
    }
}
