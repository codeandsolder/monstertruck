//! Split connected planar supports at non-manifold face-domain vertices.
//! Source incidence stays intact: only ownership of already present
//! polygons changes. A multi-source flood from the disjoint sectors of
//! one pinched vertex cuts the original support into connected children.

use super::super::{EdgeIncidence, InputPolygon, edge_key};
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};

fn opposite(incidence: &EdgeIncidence, polygon: usize, a: usize, b: usize) -> Option<usize> {
    let (key, direction) = edge_key(a, b);
    match incidence.get(&key)?.as_slice() {
        [first, second] if *first == (polygon, direction) && second.1 != direction => {
            Some(second.0)
        }
        [first, second] if *second == (polygon, direction) && first.1 != direction => Some(first.0),
        _ => None,
    }
}

fn pinch(
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    members: &HashSet<usize>,
) -> Option<usize> {
    let mut outgoing = HashSet::new();
    let mut pinches = BTreeSet::new();
    for &id in members {
        let vertices = &polygons[id].vertices;
        for (&a, &b) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
            let (key, direction) = edge_key(a, b);
            // Do not conceal ambiguous source edges by manufacturing
            // planar cuts. Ordinary native diagnostics still report them.
            match incidence.get(&key)?.as_slice() {
                [only] if *only == (id, direction) => (),
                [first, second] if first.1 != second.1 => {
                    if opposite(incidence, id, a, b).is_some_and(|other| members.contains(&other)) {
                        continue;
                    }
                }
                _ => return None,
            }
            if !outgoing.insert(a) {
                pinches.insert(a);
            }
        }
    }
    pinches.iter().nth(pinches.len() / 2).copied()
}

fn sectors(
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    members: &HashSet<usize>,
    vertex: usize,
) -> Vec<Vec<usize>> {
    let mut remaining = members
        .iter()
        .copied()
        .filter(|&id| polygons[id].vertices.contains(&vertex))
        .collect::<BTreeSet<_>>();
    let mut groups = Vec::new();
    while let Some(seed) = remaining.pop_first() {
        let mut group = vec![seed];
        let mut cursor = 0;
        while let Some(&id) = group.get(cursor) {
            let vertices = &polygons[id].vertices;
            for (&a, &b) in vertices.iter().zip(vertices.iter().cycle().skip(1)) {
                if (a == vertex || b == vertex)
                    && let Some(other) = opposite(incidence, id, a, b)
                    && remaining.remove(&other)
                {
                    group.push(other);
                }
            }
            cursor += 1;
        }
        groups.push(group);
    }
    groups
}

fn flood(
    members: &HashSet<usize>,
    neighbors: &[Vec<usize>],
    mut seeds: Vec<Vec<usize>>,
) -> Option<Vec<Vec<usize>>> {
    let mut labels = HashMap::with_capacity(members.len());
    let mut queue = VecDeque::new();
    for (label, group) in seeds.iter().enumerate() {
        for &id in group {
            labels.insert(id, label);
            queue.push_back((id, label));
        }
    }
    while let Some((id, label)) = queue.pop_front() {
        for &other in &neighbors[id] {
            if members.contains(&other)
                && let std::collections::hash_map::Entry::Vacant(slot) = labels.entry(other)
            {
                slot.insert(label);
                seeds[label].push(other);
                queue.push_back((other, label));
            }
        }
    }
    (labels.len() == members.len()).then_some(seeds)
}

pub(super) fn split(
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    neighbors: &[Vec<usize>],
    group: Vec<usize>,
) -> Vec<Vec<usize>> {
    let mut pending = vec![group];
    let mut result = Vec::new();
    while let Some(mut group) = pending.pop() {
        let members = group.iter().copied().collect::<HashSet<_>>();
        let children = pinch(polygons, incidence, &members)
            .map(|vertex| sectors(polygons, incidence, &members, vertex))
            .filter(|seeds| seeds.len() > 1)
            .and_then(|seeds| flood(&members, neighbors, seeds));
        if let Some(children) = children {
            // Every child contains a different nonempty sector, so the
            // recursion strictly decreases original polygon membership.
            pending.extend(children);
        } else {
            group.sort_unstable();
            result.push(group);
        }
    }
    result.sort_unstable_by_key(|group| group[0]);
    result
}
