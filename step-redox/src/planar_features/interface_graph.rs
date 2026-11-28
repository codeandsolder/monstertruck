//! Candidate ownership from complete source interface cycles. Cutting only
//! interface edges, rather than deleting every face which might be a host,
//! keeps potential nested hosts inside their actual feature component.

use super::{HostBoundLookup, ShellTopology};
use std::collections::{BTreeSet, HashMap, HashSet, VecDeque, hash_map::Entry};

type Owners = HashMap<u64, HashSet<u64>>;

fn cut_owners(lookups: &HashMap<u64, HostBoundLookup>) -> Owners {
    let mut owners = Owners::new();
    for (&host, bounds) in lookups {
        for (edges, matches) in bounds {
            if matches.len() == 1 {
                for &edge in edges {
                    owners.entry(edge).or_default().insert(host);
                }
            }
        }
    }
    owners
}

struct Partition {
    components: Vec<HashSet<u64>>,
    membership: HashMap<u64, usize>,
}

impl Partition {
    fn new(topology: &ShellTopology, owners: &Owners) -> Self {
        let mut membership = HashMap::new();
        let mut components = Vec::new();
        for &seed in &topology.faces {
            if membership.contains_key(&seed) {
                continue;
            }
            let index = components.len();
            membership.insert(seed, index);
            let mut queue = VecDeque::from([seed]);
            let mut component = HashSet::from([seed]);
            while let Some(face) = queue.pop_front() {
                for edge in topology.face_edges.get(&face).into_iter().flatten() {
                    if owners.contains_key(edge) {
                        continue;
                    }
                    for &neighbor in topology.edge_faces.get(edge).into_iter().flatten() {
                        if let Entry::Vacant(entry) = membership.entry(neighbor) {
                            entry.insert(index);
                            component.insert(neighbor);
                            queue.push_back(neighbor);
                        }
                    }
                }
            }
            components.push(component);
        }
        Self {
            components,
            membership,
        }
    }

    fn dependencies(
        &self,
        topology: &ShellTopology,
        owners: &Owners,
    ) -> (Vec<BTreeSet<usize>>, BTreeSet<usize>) {
        let mut required = vec![BTreeSet::new(); self.components.len()];
        let mut seeds = BTreeSet::new();
        for (edge, hosts) in owners {
            let Some([a, b]) = topology.edge_faces.get(edge).map(Vec::as_slice) else {
                continue;
            };
            let (Some(&ca), Some(&cb)) = (self.membership.get(a), self.membership.get(b)) else {
                continue;
            };
            if ca == cb {
                continue;
            }
            let (a_owns, b_owns) = (hosts.contains(a), hosts.contains(b));
            // An outside host can seal an interface by deleting its inner
            // bound. If a candidate contains the only owner of an edge, it
            // must also contain the opposite side; otherwise that edge is
            // left with no retained host to heal. Two-sided ownership is not
            // forced in either direction and is checked at candidate proof.
            if a_owns {
                seeds.insert(cb);
                if !b_owns {
                    required[ca].insert(cb);
                }
            }
            if b_owns {
                seeds.insert(ca);
                if !a_owns {
                    required[cb].insert(ca);
                }
            }
        }
        (required, seeds)
    }
}

fn dependency_closure(seed: usize, required: &[BTreeSet<usize>]) -> BTreeSet<usize> {
    let mut closure = BTreeSet::from([seed]);
    let mut queue = VecDeque::from([seed]);
    while let Some(component) = queue.pop_front() {
        for &next in &required[component] {
            if closure.insert(next) {
                queue.push_back(next);
            }
        }
    }
    closure
}

fn retained_hosts(
    topology: &ShellTopology,
    owners: &Owners,
    component: &HashSet<u64>,
) -> Option<HashSet<u64>> {
    let mut hosts = HashSet::new();
    for face in component {
        for edge in topology.face_edges.get(face)? {
            let adjacent = topology.edge_faces.get(edge)?;
            if adjacent.len() != 2 {
                return None;
            }
            for &outside in adjacent.iter().filter(|id| !component.contains(id)) {
                if !owners.get(edge).is_some_and(|ids| ids.contains(&outside)) {
                    return None;
                }
                hosts.insert(outside);
            }
        }
    }
    (!hosts.is_empty()).then_some(hosts)
}

/// A feature candidate is the minimal closure of a witnessed cut-side
/// component under source ownership dependencies. This discovers arbitrary
/// subsets of retained hosts without trying the power set of candidate faces.
pub(super) fn components_for_interfaces(
    topology: &ShellTopology,
    lookups: &HashMap<u64, HostBoundLookup>,
) -> Vec<(HashSet<u64>, HashSet<u64>)> {
    let owners = cut_owners(lookups);
    if owners.is_empty() {
        return Vec::new();
    }
    let partition = Partition::new(topology, &owners);
    let (required, seeds) = partition.dependencies(topology, &owners);
    let mut seen = BTreeSet::new();
    seeds
        .into_iter()
        .filter_map(|seed| {
            let closure = dependency_closure(seed, &required);
            if !seen.insert(closure.clone()) || closure.len() == partition.components.len() {
                return None;
            }
            let component = closure
                .into_iter()
                .flat_map(|i| partition.components[i].iter().copied())
                .collect::<HashSet<_>>();
            let hosts = retained_hosts(topology, &owners, &component)?;
            Some((component, hosts))
        })
        .collect()
}

#[cfg(test)]
mod tests;
