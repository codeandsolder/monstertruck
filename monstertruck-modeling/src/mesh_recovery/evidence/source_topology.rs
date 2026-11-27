//! Preserve independently closed oriented source components at geometric contacts.
//! Exact geometric welding is evidence for coincident positions, not for
//! identifying the topological vertices of distinct closed surface chains.

use super::super::{EdgeIncidence, InputPolygon, source_edge_incidence};
use super::{
    SourceProvenance,
    source_sectors::{SourceSectorPair, paired_sectors},
};
use crate::Point3;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::collections::VecDeque;

mod forced;
mod multiplicity;
pub use forced::SourceForcedPair;
pub use multiplicity::{
    SourceMultiplicityMethod, SourceMultiplicityPair, SourceMultiplicityResolution,
};

/// A distinct topological copy of an exactly coincident source position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeparatedSourceVertex {
    /// Original exact-welded vertex ID in the position table.
    pub source_vertex: usize,
    /// New topological vertex ID, with unchanged coordinates.
    pub separated_vertex: usize,
}

// Outward-rounded interval arithmetic for the oriented volume integral.
// A component with zero/unresolved volume (e.g. an attached double sheet)
// must not be isolated as a material body. Inward-oriented cavity shells
// are legitimate, so the sign itself is not constrained.
#[derive(Clone, Copy)]
struct Interval {
    low: f64,
    high: f64,
}

impl Interval {
    const ZERO: Self = Self {
        low: 0.0,
        high: 0.0,
    };
    const UNKNOWN: Self = Self {
        low: f64::NEG_INFINITY,
        high: f64::INFINITY,
    };

    fn difference(a: f64, b: f64) -> Self {
        let d = a - b;
        Self {
            low: d.next_down(),
            high: d.next_up(),
        }
    }

    fn add(self, other: Self) -> Self {
        Self {
            low: (self.low + other.low).next_down(),
            high: (self.high + other.high).next_up(),
        }
    }

    fn subtract(self, other: Self) -> Self {
        Self {
            low: (self.low - other.high).next_down(),
            high: (self.high - other.low).next_up(),
        }
    }

    fn multiply(self, other: Self) -> Self {
        let products = [
            self.low * other.low,
            self.low * other.high,
            self.high * other.low,
            self.high * other.high,
        ];
        if products.iter().any(|x| !x.is_finite()) {
            return Self::UNKNOWN;
        }
        Self {
            low: products
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min)
                .next_down(),
            high: products
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        }
    }

    fn resolves_sign(self, positive: Option<bool>) -> bool {
        self.low.is_finite()
            && self.high.is_finite()
            && ((positive != Some(false) && self.low > 0.0)
                || (positive != Some(true) && self.high < 0.0))
    }
}

#[derive(Clone, Copy)]
struct ComponentVolume {
    origin: Point3,
    six_volume: Interval,
}

impl ComponentVolume {
    fn add_polygon(&mut self, positions: &[Point3], vertices: &[usize]) {
        let local = |id: usize| {
            std::array::from_fn::<_, 3, _>(|i| {
                Interval::difference(positions[id][i], self.origin[i])
            })
        };
        let first = local(vertices[0]);
        for pair in vertices[1..].windows(2) {
            let second = local(pair[0]);
            let third = local(pair[1]);
            let mut determinant = Interval::ZERO;
            for (i, coordinate) in first.into_iter().enumerate() {
                let j = (i + 1) % 3;
                let k = (i + 2) % 3;
                determinant = determinant.add(
                    coordinate.multiply(
                        second[j]
                            .multiply(third[k])
                            .subtract(second[k].multiply(third[j])),
                    ),
                );
            }
            self.six_volume = self.six_volume.add(determinant);
        }
    }
}

// Follow only regular oriented two-use edges to identify independent
// source-connected components. A component is closed when every edge,
// including those shared geometrically with other components, has exactly
// two opposite uses within that component. No pairing is inferred between
// components merely because their geometry coincides.
fn connected_membership(
    polygon_count: usize,
    incidence: &EdgeIncidence,
    pairs: impl IntoIterator<Item = [usize; 2]>,
) -> (Vec<usize>, usize) {
    let mut neighbors = vec![Vec::new(); polygon_count];
    for uses in incidence.values() {
        if let [a, b] = uses.as_slice()
            && a.1 != b.1
        {
            neighbors[a.0].push(b.0);
            neighbors[b.0].push(a.0);
        }
    }
    for [a, b] in pairs {
        neighbors[a].push(b);
        neighbors[b].push(a);
    }
    let mut membership = vec![usize::MAX; polygon_count];
    let mut count = 0;
    for seed in 0..polygon_count {
        if membership[seed] != usize::MAX {
            continue;
        }
        let id = count;
        count += 1;
        membership[seed] = id;
        let mut queue = VecDeque::from([seed]);
        while let Some(current) = queue.pop_front() {
            for &neighbor in &neighbors[current] {
                if membership[neighbor] == usize::MAX {
                    membership[neighbor] = id;
                    queue.push_back(neighbor);
                }
            }
        }
    }
    (membership, count)
}

fn closed_components(
    positions: &[Point3],
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    pairs: impl IntoIterator<Item = [usize; 2]>,
    positive: Option<bool>,
) -> Vec<Option<usize>> {
    let (membership, count) = connected_membership(polygons.len(), incidence, pairs);
    let mut closed = vec![true; count];
    let mut uses_by_component = Vec::new();
    for uses in incidence.values() {
        if matches!(uses.as_slice(), [a, b] if a.1 != b.1) {
            continue;
        }
        uses_by_component.clear();
        uses_by_component.extend(
            uses.iter()
                .map(|&(polygon, forward)| (membership[polygon], forward)),
        );
        uses_by_component.sort_unstable();
        for group in uses_by_component.chunk_by(|a, b| a.0 == b.0) {
            if !matches!(group, [a, b] if a.1 != b.1) {
                closed[group[0].0] = false;
            }
        }
    }
    let mut volumes = vec![None; closed.len()];
    for (polygon, &id) in polygons.iter().zip(&membership) {
        if closed[id] {
            volumes[id]
                .get_or_insert(ComponentVolume {
                    origin: positions[polygon.vertices[0]],
                    six_volume: Interval::ZERO,
                })
                .add_polygon(positions, &polygon.vertices);
        }
    }
    for (closed, volume) in closed.iter_mut().zip(volumes) {
        *closed &= volume.is_some_and(|v| v.six_volume.resolves_sign(positive));
    }
    membership
        .into_iter()
        .map(|id| closed[id].then_some(id))
        .collect()
}

/// Separate exactly coincident vertex IDs belonging to distinct independently
/// closed components. Keep all unresolved/open polygons in one common pool,
/// and leave their incidence untouched. The result flags the closed polygons
/// that must not be subjected to geometric T-junction search against that pool.
pub(super) fn separate_closed_components(
    positions: &mut Vec<Point3>,
    polygons: &mut [InputPolygon],
    provenance: &mut SourceProvenance,
) -> Vec<bool> {
    let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
    let components = closed_components(positions, polygons, &incidence, [], None);
    separate_components(
        positions,
        polygons,
        &components,
        &mut provenance.separated_vertices,
    );
    let mut closed = components.iter().map(Option::is_some).collect::<Vec<_>>();
    if closed.contains(&false) {
        let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
        let selected =
            select_sector_components(positions, polygons, &incidence, &closed, provenance);
        separate_components(
            positions,
            polygons,
            &selected,
            &mut provenance.separated_vertices,
        );
        for (flag, component) in closed.iter_mut().zip(selected) {
            *flag |= component.is_some();
        }
    }
    if closed.contains(&false) {
        let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
        let pairs = forced::necessary_joins(polygons.len(), &incidence, None);
        if !pairs.is_empty() {
            let mut selected = closed_components(
                positions,
                polygons,
                &incidence,
                pairs.iter().map(|pair| pair.polygons),
                None,
            );
            for (component, &already_closed) in selected.iter_mut().zip(&closed) {
                if already_closed {
                    *component = None;
                }
            }
            for mut pair in pairs {
                let [a, b] = pair.polygons;
                if selected[a].is_some() && selected[a] == selected[b] {
                    pair.polygons = pair.polygons.map(|id| provenance.original_polygon_ids[id]);
                    provenance.forced_pairs.push(pair);
                }
            }
            separate_components(
                positions,
                polygons,
                &selected,
                &mut provenance.separated_vertices,
            );
            for (flag, component) in closed.iter_mut().zip(selected) {
                *flag |= component.is_some();
            }
        }
    }
    if closed.contains(&false) {
        let incidence = source_edge_incidence(polygons.iter().map(|p| p.vertices.as_slice()));
        let selected =
            multiplicity::select_components(positions, polygons, &incidence, &closed, provenance);
        separate_components(
            positions,
            polygons,
            &selected,
            &mut provenance.separated_vertices,
        );
        for (flag, component) in closed.iter_mut().zip(selected) {
            *flag |= component.is_some();
        }
    }
    closed
}

fn sector_candidates(
    positions: &[Point3],
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
) -> [Vec<SourceSectorPair>; 2] {
    let mut candidates = [Vec::new(), Vec::new()];
    for (&edge, uses) in incidence {
        if let Some(pairings) = paired_sectors(positions, polygons, edge, uses) {
            for (index, pairs) in pairings.into_iter().enumerate() {
                if let Some(pairs) = pairs {
                    candidates[index].extend(pairs.into_iter().map(|polygons| SourceSectorPair {
                        edge,
                        polygons,
                        positive_volume: index == 0,
                    }));
                }
            }
        }
    }
    for pairs in &mut candidates {
        pairs.sort_unstable_by_key(|p| (p.edge, p.polygons));
    }
    candidates
}

fn select_sector_components(
    positions: &[Point3],
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    closed: &[bool],
    provenance: &mut SourceProvenance,
) -> Vec<Option<usize>> {
    let pairs = sector_candidates(positions, polygons, incidence);
    if pairs.iter().all(Vec::is_empty) {
        return vec![None; polygons.len()];
    }
    let memberships = std::array::from_fn::<_, 2, _>(|index| {
        closed_components(
            positions,
            polygons,
            incidence,
            pairs[index].iter().map(|pair| pair.polygons),
            Some(index == 0),
        )
    });
    // A source polygon cannot belong to two incompatible sign-certified
    // component explanations. Reject every overlapping component as a
    // whole, not individual polygons that would break its closed chain.
    let mut conflicts = [HashSet::default(), HashSet::default()];
    for (id, &already_closed) in closed.iter().enumerate() {
        if !already_closed && let (Some(a), Some(b)) = (memberships[0][id], memberships[1][id]) {
            conflicts[0].insert(a);
            conflicts[1].insert(b);
        }
    }
    let mut group_ids = HashMap::default();
    let mut selected = vec![None; polygons.len()];
    for (id, &already_closed) in closed.iter().enumerate() {
        if already_closed {
            continue;
        }
        for (index, membership) in memberships.iter().enumerate() {
            if let Some(component) = membership[id]
                && !conflicts[index].contains(&component)
            {
                let next = group_ids.len();
                selected[id] = Some(*group_ids.entry((index, component)).or_insert(next));
            }
        }
    }
    for (index, pairs) in pairs.into_iter().enumerate() {
        for mut pair in pairs {
            let [a, b] = pair.polygons;
            if let Some(component) = memberships[index][a]
                && group_ids
                    .get(&(index, component))
                    .is_some_and(|&group| selected[a] == Some(group) && selected[b] == Some(group))
            {
                pair.polygons = pair.polygons.map(|id| provenance.original_polygon_ids[id]);
                provenance.sector_pairs.push(pair);
            }
        }
    }
    selected
}

fn separate_components(
    positions: &mut Vec<Point3>,
    polygons: &mut [InputPolygon],
    components: &[Option<usize>],
    separated_vertices: &mut Vec<SeparatedSourceVertex>,
) {
    let mut owner = vec![None; positions.len()];
    // Reserve original IDs for the unresolved pool. Never change its edges.
    for (polygon, component) in polygons.iter().zip(components) {
        if component.is_none() {
            for &vertex in &polygon.vertices {
                owner[vertex] = Some(usize::MAX);
            }
        }
    }
    let mut copies = HashMap::default();
    for (polygon, component) in polygons.iter_mut().zip(components) {
        let Some(component) = *component else {
            continue;
        };
        for vertex in &mut polygon.vertices {
            let original = *vertex;
            match owner[original] {
                None => owner[original] = Some(component),
                Some(id) if id == component => {}
                Some(_) => {
                    *vertex = *copies.entry((original, component)).or_insert_with(|| {
                        let separated_vertex = positions.len();
                        positions.push(positions[original]);
                        separated_vertices.push(SeparatedSourceVertex {
                            source_vertex: original,
                            separated_vertex,
                        });
                        separated_vertex
                    });
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::Point3;
    use crate::mesh_recovery::{
        AnalyticRecoveryOptions, analyze_mesh_evidence, recover_analytic_mesh,
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_topology::shell::ShellCondition;

    #[test]
    fn independent_closed_tetrahedra_keep_distinct_coincident_edges()
    -> Result<(), Box<dyn std::error::Error>> {
        // The second tetrahedron touches the first along either its full
        // edge or an exactly coincident subsegment with an extra endpoint.
        // Neither is a T-junction connecting their surface chains.
        for partial in [false, true] {
            let start = usize::from(partial) * 4;
            let mesh = PolygonMesh::new(
                StandardAttributes {
                    positions: vec![
                        Point3::new(0.0, 0.0, 0.0),
                        Point3::new(2.0, 0.0, 0.0),
                        Point3::new(0.0, 1.0, 0.0),
                        Point3::new(0.0, 0.0, 1.0),
                        Point3::new(1.0, 0.0, 0.0),
                        Point3::new(1.0, -1.0, 0.0),
                        Point3::new(1.0, 0.0, -1.0),
                    ],
                    ..Default::default()
                },
                Faces::from_iter([
                    [0, 2, 1],
                    [0, 1, 3],
                    [0, 3, 2],
                    [1, 2, 3],
                    [start, 5, 1],
                    [start, 1, 6],
                    [start, 6, 5],
                    [1, 5, 6],
                ]),
            )?;
            let evidence = analyze_mesh_evidence(&mesh, Default::default())?;
            assert_eq!(evidence.non_manifold_edge_count, 0);
            assert_eq!(evidence.inconsistent_winding_edge_count, 0);
            assert_eq!(
                evidence.provenance.separated_vertices.len(),
                if partial { 1 } else { 2 }
            );
            for copy in &evidence.provenance.separated_vertices {
                assert_eq!(
                    evidence.positions[copy.source_vertex],
                    evidence.positions[copy.separated_vertex]
                );
            }
            assert!(
                evidence
                    .source_polygon_vertex_offsets
                    .windows(2)
                    .all(|p| p[1] - p[0] == 3)
            );
            assert!(evidence.closed_source_volume().is_some());
            let recovered = recover_analytic_mesh(&mesh, AnalyticRecoveryOptions::default())?;
            let shell = recovered
                .complete_shell()
                .ok_or("two separate closed shells")?;
            assert_eq!(shell.shell_condition(), ShellCondition::Closed);
            assert_eq!(shell.connected_components().len(), 2);
            assert_eq!(shell.len(), 8);
            assert!(recovered.evidence.provenance.canceled_sheets.is_empty());
        }
        Ok(())
    }
}
