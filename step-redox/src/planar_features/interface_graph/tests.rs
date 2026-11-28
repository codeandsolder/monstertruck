use super::*;
use crate::planar_features::matched_host_bounds;

fn topology(faces: Vec<u64>, edges: &[(u64, [u64; 2])]) -> ShellTopology {
    let edge_faces = edges
        .iter()
        .map(|(edge, faces)| (*edge, faces.to_vec()))
        .collect();
    let mut face_edges = HashMap::<u64, HashSet<u64>>::new();
    for &(edge, attached) in edges {
        for face in attached {
            face_edges.entry(face).or_default().insert(edge);
        }
    }
    ShellTopology {
        faces,
        face_edges,
        edge_faces,
        adjacency: HashMap::new(),
    }
}

fn bound(edge: u64) -> (Vec<u64>, Vec<(u64, u64, String)>) {
    (vec![edge], vec![(edge + 100, edge + 200, "T".into())])
}

#[test]
fn proper_host_subset_preserves_nested_feature_ownership() {
    // The feature joins two retained hosts, contains a third eligible
    // host, and has a nested feature on it. An unrelated feature joins
    // the second retained host.
    let shell = topology(
        vec![1, 2, 3, 4, 5, 6, 7],
        &[
            (10, [1, 2]),
            (11, [3, 4]),
            (12, [5, 6]),
            (20, [1, 3]),
            (21, [2, 4]),
            (22, [4, 5]),
            (23, [2, 7]),
        ],
    );
    let lookups = HashMap::from([
        (1, HashMap::from([bound(20)])),
        (2, HashMap::from([bound(21), bound(23)])),
        (4, HashMap::from([bound(22)])),
    ]);
    let candidates = components_for_interfaces(&shell, &lookups);
    assert_eq!(candidates.len(), 3);
    assert!(candidates.contains(&(HashSet::from([3, 4, 5, 6]), HashSet::from([1, 2]))));
    assert!(candidates.contains(&(HashSet::from([5, 6]), HashSet::from([4]))));
    assert!(candidates.contains(&(HashSet::from([7]), HashSet::from([2]))));
}

#[test]
fn multiple_complete_bounds_on_same_host_are_not_one_motif() {
    let lookup = HashMap::from([
        (vec![1, 2], vec![(101, 201, "T".into())]),
        (vec![3, 4], vec![(102, 202, "F".into())]),
    ]);
    let result = matched_host_bounds(10, &HashSet::from([1, 2, 3, 4]), &lookup)
        .expect("all exposed edges are two complete host bounds");
    assert_eq!(result.len(), 2);
    assert!(matched_host_bounds(10, &HashSet::from([1, 2, 3]), &lookup).is_none());
    assert!(matched_host_bounds(10, &HashSet::from([1, 2, 3, 4, 5]), &lookup).is_none());
}

#[test]
fn seam_uses_do_not_create_false_nonmanifold_boundary() {
    let shell = topology(vec![1, 2], &[(10, [1, 2]), (11, [2, 2])]);
    let lookups = HashMap::from([(1, HashMap::from([bound(10)]))]);
    assert_eq!(
        components_for_interfaces(&shell, &lookups),
        vec![(HashSet::from([2]), HashSet::from([1]))]
    );
}
