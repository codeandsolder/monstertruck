//! Connected boundary evidence between seed-anchored planar facet supports.
//!
//! An interface is a connected component of source mesh edges between two
//! observed local supports, not an asserted analytic CAD trim. A curved chord
//! patch can yield the same topology as a designed planar region.

use super::{
    super::{EdgeIncidence, InputPolygon},
    canonical_unit,
};
use crate::{InnerProduct, Point3, Vector3};
use rustc_hash::FxHashMap as HashMap;
use std::collections::VecDeque;

/// Graph shape of one source-contact component between two local supports.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SupportBoundaryTopology {
    /// Exactly two degree-one endpoints, with all other vertices degree two.
    OpenPath,
    /// All vertices have degree two, with no dangling endpoints.
    ClosedLoop,
    /// Source edges branch or otherwise do not form a single simple path/cycle.
    Branched,
}

/// One connected component of source edges joining two local planar supports.
///
/// A component may bend, form a loop, or branch. Seed-line measurements are
/// descriptive and do not certify that it belongs to a designed straight edge.
#[derive(Clone, Debug)]
pub struct SupportInterface {
    /// Canonically sorted local support IDs.
    pub supports: [usize; 2],
    /// Canonical original source-position index pairs, one per source edge.
    pub source_edges: Vec<[usize; 2]>,
    /// Original source polygon IDs on either side, in support-ID order.
    /// This vector is aligned with the original source edge list.
    pub source_polygon_pairs: Vec<[usize; 2]>,
    /// Exact edge graph classification of this connected component.
    pub topology: SupportBoundaryTopology,
    /// Sum of source polygonal boundary segment lengths.
    pub boundary_length: f64,
    /// Minimum local source facet-normal jump, in radians.
    pub min_turn_radians: f64,
    /// Maximum local source facet-normal jump, in radians.
    pub max_turn_radians: f64,
    /// Source edge-length-weighted average unsigned normal jump.
    pub mean_turn_radians: f64,
    /// Canonical undirected direction of the longest source segment.
    pub seed_direction: Vector3,
    /// Largest endpoint distance from the longest segment's supporting line.
    /// Zero means exact collinearity only up to floating-point arithmetic.
    pub max_seed_line_offset: f64,
    /// Largest undirected segment angle from the longest segment, in radians.
    pub max_seed_direction_spread_radians: f64,
}

#[derive(Clone, Copy)]
struct SourceContact {
    edge: [usize; 2],
    polygons: [usize; 2],
    length: f64,
    turn: f64,
}

fn source_interface_component(
    positions: &[Point3],
    supports: [usize; 2],
    contacts: &[SourceContact],
) -> Option<SupportInterface> {
    let mut source_edges = Vec::with_capacity(contacts.len());
    let mut source_polygon_pairs = Vec::with_capacity(contacts.len());
    let mut boundary_length = 0.0;
    let mut angular_sum = 0.0;
    let mut min_turn_radians = f64::INFINITY;
    let mut max_turn_radians = 0.0_f64;
    let mut vertex_degree = HashMap::<usize, usize>::default();
    let longest = contacts.iter().max_by(|a, b| {
        a.length
            .total_cmp(&b.length)
            .then_with(|| b.edge.cmp(&a.edge))
    })?;
    let seed_origin = positions[longest.edge[0]];
    let seed_direction = canonical_unit(positions[longest.edge[1]] - positions[longest.edge[0]]);
    let mut max_seed_line_offset = 0.0_f64;
    let mut max_seed_direction_spread_radians = 0.0_f64;
    for contact in contacts {
        source_edges.push(contact.edge);
        source_polygon_pairs.push(contact.polygons);
        boundary_length += contact.length;
        angular_sum = contact.turn.mul_add(contact.length, angular_sum);
        min_turn_radians = min_turn_radians.min(contact.turn);
        max_turn_radians = max_turn_radians.max(contact.turn);
        for &vertex in &contact.edge {
            *vertex_degree.entry(vertex).or_default() += 1;
            let offset = positions[vertex] - seed_origin;
            max_seed_line_offset = max_seed_line_offset.max(offset.cross(&seed_direction).norm());
        }
        let direction = (positions[contact.edge[1]] - positions[contact.edge[0]]).normalize();
        let spread = seed_direction
            .cross(&direction)
            .norm()
            .atan2(seed_direction.dot(&direction).abs());
        max_seed_direction_spread_radians = max_seed_direction_spread_radians.max(spread);
    }
    let endpoint_count = vertex_degree
        .values()
        .filter(|&&degree| degree == 1)
        .count();
    let topology = if vertex_degree.values().all(|&degree| degree == 2) {
        SupportBoundaryTopology::ClosedLoop
    } else if endpoint_count == 2 && vertex_degree.values().all(|&degree| degree <= 2) {
        SupportBoundaryTopology::OpenPath
    } else {
        SupportBoundaryTopology::Branched
    };
    Some(SupportInterface {
        supports,
        source_edges,
        source_polygon_pairs,
        topology,
        boundary_length,
        min_turn_radians,
        max_turn_radians,
        mean_turn_radians: angular_sum / boundary_length,
        seed_direction,
        max_seed_line_offset,
        max_seed_direction_spread_radians,
    })
}

/// Build independently connected interfaces per unordered support pair.
///
/// Each source edge is used at most once; malformed or ambiguous source
/// incidence never generates an apparent two-sided interface.
pub(super) fn build_support_interfaces(
    positions: &[Point3],
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    membership: &[usize],
    region_count: usize,
) -> (Vec<SupportInterface>, Vec<Vec<usize>>) {
    let mut pair_contacts = HashMap::<[usize; 2], Vec<SourceContact>>::default();
    for (&edge, uses) in incidence {
        let [(a, a_forward), (b, b_forward)] = uses.as_slice() else {
            continue;
        };
        if a_forward == b_forward {
            continue;
        }
        let first = membership[*a];
        let second = membership[*b];
        if first == second {
            continue;
        }
        let (pair, polygons_ordered) = if first < second {
            ([first, second], [*a, *b])
        } else {
            ([second, first], [*b, *a])
        };
        let direction = positions[edge[1]] - positions[edge[0]];
        let turn = polygons[*a]
            .normal
            .cross(&polygons[*b].normal)
            .norm()
            .atan2(polygons[*a].normal.dot(&polygons[*b].normal));
        pair_contacts.entry(pair).or_default().push(SourceContact {
            edge,
            polygons: polygons_ordered,
            length: direction.norm(),
            turn,
        });
    }
    let mut keys = pair_contacts.keys().copied().collect::<Vec<_>>();
    keys.sort_unstable();
    let mut interfaces = Vec::new();
    let mut region_interfaces = vec![Vec::new(); region_count];
    let mut append_component = |pair, component: &[SourceContact]| {
        if let Some(interface) = source_interface_component(positions, pair, component) {
            let id = interfaces.len();
            interfaces.push(interface);
            region_interfaces[pair[0]].push(id);
            region_interfaces[pair[1]].push(id);
        }
    };
    for pair in keys {
        let Some(mut edges) = pair_contacts.remove(&pair) else {
            continue;
        };
        edges.sort_unstable_by_key(|contact| contact.edge);
        // A single source edge is already a complete connected component.
        // This dominates real CAD tessellations and avoids a per-edge graph.
        if edges.len() == 1 {
            append_component(pair, &edges);
            continue;
        }
        let mut vertex_edges = HashMap::<usize, Vec<usize>>::default();
        for (id, contact) in edges.iter().enumerate() {
            for &vertex in &contact.edge {
                vertex_edges.entry(vertex).or_default().push(id);
            }
        }
        let mut visited = vec![false; edges.len()];
        for first in 0..edges.len() {
            if visited[first] {
                continue;
            }
            visited[first] = true;
            let mut queue = VecDeque::from([first]);
            let mut component = Vec::new();
            while let Some(current) = queue.pop_front() {
                component.push(edges[current]);
                for &vertex in &edges[current].edge {
                    for &next in &vertex_edges[&vertex] {
                        if !visited[next] {
                            visited[next] = true;
                            queue.push_back(next);
                        }
                    }
                }
            }
            component.sort_unstable_by_key(|contact| contact.edge);
            append_component(pair, &component);
        }
    }
    (interfaces, region_interfaces)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contact(points: &[Point3], a: usize, b: usize) -> SourceContact {
        let edge = [a.min(b), a.max(b)];
        SourceContact {
            edge,
            polygons: [0, 1],
            length: (points[a] - points[b]).norm(),
            turn: 0.125,
        }
    }

    #[test]
    fn classifies_connected_closed_polyline_boundary_as_a_loop() {
        let positions = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let contacts = [
            contact(&positions, 0, 1),
            contact(&positions, 1, 2),
            contact(&positions, 2, 3),
            contact(&positions, 3, 0),
        ];
        let result = source_interface_component(&positions, [1, 2], &contacts).unwrap();
        assert_eq!(result.topology, SupportBoundaryTopology::ClosedLoop);
        assert_eq!(result.source_edges.len(), 4);
        assert!((result.boundary_length - 4.0).abs() < 1.0e-12);
        assert!((result.mean_turn_radians - 0.125).abs() < 1.0e-12);
        assert!((result.max_seed_line_offset - 1.0).abs() < 1.0e-12);
        assert!(
            (result.max_seed_direction_spread_radians - std::f64::consts::FRAC_PI_2).abs()
                < 1.0e-12
        );
    }

    #[test]
    fn classifies_nonmanifold_contact_graph_as_branched_not_as_a_trim_wire() {
        let positions = [
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let contacts = [
            contact(&positions, 0, 1),
            contact(&positions, 1, 2),
            contact(&positions, 1, 3),
        ];
        let result = source_interface_component(&positions, [1, 2], &contacts).unwrap();
        assert_eq!(result.topology, SupportBoundaryTopology::Branched);
        assert_eq!(result.source_edges.len(), 3);
        assert_eq!(result.source_polygon_pairs.len(), 3);
    }

    #[test]
    fn empty_contact_collection_is_not_a_topological_boundary() {
        assert!(source_interface_component(&[], [0, 1], &[]).is_none());
    }
}
