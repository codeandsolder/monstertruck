//! Recover source connectivity when global oriented incidence determines
//! a unique partition, modulo interchangeable exact geometric fragments.
//! Preserve all polygon multiplicity and exact coordinates. Do not select
//! a geometrically inequivalent partition. A uniform oriented quotient
//! preserves repeated coverage without asserting material body identity.
use super::{
    EdgeIncidence, InputPolygon, SourceProvenance, closed_components, connected_membership,
};
use crate::Point3;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

mod class_graph;
mod equivalence;
mod exclusions;
mod peel;
mod propagation;
mod search;
mod transverse;
mod uniform;
mod vertex_links;
use equivalence::Profile;
use search::Constraint;

/// One topological pairing selected modulo exact fragment equivalence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SourceMultiplicityPair {
    /// Canonical edge before source component separation.
    pub edge: [usize; 2],
    /// Original polygon IDs assigned to the same shell across the edge.
    pub polygons: [usize; 2],
}

/// Proof used to establish the unique geometric partition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceMultiplicityMethod {
    /// A unique closed geometric chain repeated with uniform integer coverage.
    /// Does not distinguish coincident bodies from duplicated source records.
    UniformOrientedCoverage,
    /// All anchored shells were necessarily closed on exact class incidence.
    ForcedClassClosure,
    /// Global necessary joins under exact geometric incompatibility constraints.
    GlobalForcedMatching,
    /// Bounded exhaustive search of distinct geometric partitions.
    QuotientSearch,
}

/// A proper crossing of an original polygon edge through another face's interior.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SourceTransverseCrossing {
    /// Original polygon IDs, with the crossing segment belonging to the first.
    pub polygons: [usize; 2],
    /// Exact-welded endpoint IDs before source component separation.
    pub edge: [usize; 2],
}

/// A source-incidence proof that found one geometric partition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SourceMultiplicityResolution {
    /// Geometry/incident partition proof method.
    pub method: SourceMultiplicityMethod,
    /// Number of explored geometric states or forced closure steps.
    pub states: usize,
    /// Interchangeable source fragments, grouped by geometric equivalence.
    /// Each fragment retains every original polygon ID.
    pub equivalent_fragments: Vec<Vec<Vec<usize>>>,
    /// Selected oriented edge uses, preserving their source multiplicity.
    pub pairs: Vec<SourceMultiplicityPair>,
    /// Original polygon pairs witnessing exact opposed-face exclusions.
    pub opposed_faces: Vec<[usize; 2]>,
    /// Proper transverse crossings excluding same-shell assignments.
    pub transverse_crossings: Vec<SourceTransverseCrossing>,
}

#[derive(Clone, Copy)]
struct SourceUse {
    fragment: usize,
    polygon: usize,
    forward: bool,
}

struct AmbiguousEdge {
    edge: [usize; 2],
    uses: Vec<SourceUse>,
}

struct Cluster {
    groups: Vec<usize>,
    edges: Vec<usize>,
}

fn clusters(groups: usize, edges: &[AmbiguousEdge]) -> Vec<Cluster> {
    let mut incident = vec![Vec::new(); groups];
    for (id, edge) in edges.iter().enumerate() {
        for source in &edge.uses {
            incident[source.fragment].push(id);
        }
    }
    let mut seen = vec![false; edges.len()];
    let mut clusters = Vec::new();
    for seed in 0..edges.len() {
        if seen[seed] {
            continue;
        }
        seen[seed] = true;
        let mut queue = VecDeque::from([seed]);
        let mut cluster_edges = Vec::new();
        let mut members = BTreeSet::new();
        while let Some(id) = queue.pop_front() {
            cluster_edges.push(id);
            for source in &edges[id].uses {
                if members.insert(source.fragment) {
                    for &next in &incident[source.fragment] {
                        if !seen[next] {
                            seen[next] = true;
                            queue.push_back(next);
                        }
                    }
                }
            }
        }
        cluster_edges.sort_unstable();
        clusters.push(Cluster {
            groups: members.into_iter().collect(),
            edges: cluster_edges,
        });
    }
    clusters
}

struct ClusterInput<'a> {
    fragments: Vec<&'a [usize]>,
    edges: Vec<[usize; 2]>,
    constraints: Vec<Constraint>,
    profiles: Vec<Profile>,
    uses: Vec<Vec<SourceUse>>,
}

impl<'a> ClusterInput<'a> {
    fn new(
        cluster: &Cluster,
        fragments: &'a [Vec<usize>],
        edges: &[AmbiguousEdge],
    ) -> Option<Self> {
        let fragments = cluster
            .groups
            .iter()
            .map(|&id| fragments[id].as_slice())
            .collect::<Vec<_>>();
        let mut profiles = vec![Vec::new(); fragments.len()];
        let mut constraints = Vec::new();
        let mut edge_keys = Vec::new();
        let mut uses = Vec::new();
        for (local, &id) in cluster.edges.iter().enumerate() {
            let edge = &edges[id];
            let mut constraint = Constraint {
                forward: Vec::new(),
                reverse: Vec::new(),
            };
            let mut local_uses = Vec::new();
            for &source in &edge.uses {
                let component = cluster.groups.binary_search(&source.fragment).ok()?;
                profiles[component].push((local, source.forward));
                local_uses.push(SourceUse {
                    fragment: component,
                    ..source
                });
                if source.forward {
                    constraint.forward.push(component);
                } else {
                    constraint.reverse.push(component);
                }
            }
            if constraint.forward.len() != constraint.reverse.len() || constraint.forward.len() < 2
            {
                return None;
            }
            constraints.push(constraint);
            edge_keys.push(edge.edge);
            uses.push(local_uses);
        }
        for profile in &mut profiles {
            profile.sort_unstable();
            if profile.windows(2).any(|p| p[0] == p[1]) {
                return None;
            }
        }
        Some(Self {
            fragments,
            edges: edge_keys,
            constraints,
            profiles,
            uses,
        })
    }

    fn solve(
        &self,
        classes: &[usize],
        graph: &class_graph::ClassGraph,
    ) -> Option<(SourceMultiplicityMethod, search::Solution)> {
        if let Some(solution) = uniform::solve(graph) {
            Some((SourceMultiplicityMethod::UniformOrientedCoverage, solution))
        } else if let Some(solution) = peel::solve(graph) {
            Some((SourceMultiplicityMethod::ForcedClassClosure, solution))
        } else {
            let joins = propagation::necessary_joins(&self.constraints, &self.edges, graph);
            if let Some(solution) = propagation::complete(&self.constraints, graph, &joins) {
                Some((SourceMultiplicityMethod::GlobalForcedMatching, solution))
            } else {
                Some((
                    SourceMultiplicityMethod::QuotientSearch,
                    search::solve(&self.constraints, classes, graph, &joins)?,
                ))
            }
        }
    }

    fn resolve(
        &self,
        positions: &[Point3],
        polygons: &[InputPolygon],
    ) -> Option<SourceMultiplicityResolution> {
        let classes = equivalence::classes(
            positions,
            polygons,
            &self.fragments,
            &self.edges,
            &self.profiles,
        );
        let graph = class_graph::ClassGraph::new(&self.constraints, &classes, &[])?;
        let initial = self.solve(&classes, &graph);
        let mut opposed_faces = Vec::new();
        let mut transverse_crossings = Vec::new();
        let (method, solution) = if let Some(result) = initial {
            result
        } else {
            let exclusions =
                exclusions::opposed_faces(positions, polygons, &self.fragments, &classes);
            let mut labels = exclusions.iter().map(|e| e.classes).collect::<Vec<_>>();
            opposed_faces.extend(exclusions.into_iter().map(|e| e.polygons));
            let opposed = if labels.is_empty() {
                None
            } else {
                let graph = class_graph::ClassGraph::new(&self.constraints, &classes, &labels)?;
                self.solve(&classes, &graph)
            };
            if let Some(result) = opposed {
                result
            } else {
                let crossings =
                    transverse::crossings(positions, polygons, &self.fragments, &classes, &labels);
                if crossings.is_empty() {
                    return None;
                }
                for (pair, witness) in crossings {
                    labels.push(pair);
                    transverse_crossings.push(witness);
                }
                let graph = class_graph::ClassGraph::new(&self.constraints, &classes, &labels)?;
                self.solve(&classes, &graph)?
            }
        };
        let mut pairs = Vec::new();
        for (&edge, uses) in self.edges.iter().zip(&self.uses) {
            for source in uses {
                if !source.forward {
                    continue;
                }
                let root = solution.roots[source.fragment];
                let other = uses
                    .iter()
                    .find(|other| !other.forward && solution.roots[other.fragment] == root)?;
                pairs.push(SourceMultiplicityPair {
                    edge,
                    polygons: [source.polygon, other.polygon],
                });
            }
        }
        let mut equivalent = BTreeMap::<usize, Vec<Vec<usize>>>::new();
        for (&class, fragment) in classes.iter().zip(&self.fragments) {
            equivalent.entry(class).or_default().push(fragment.to_vec());
        }
        let equivalent_fragments = equivalent
            .into_values()
            .filter(|group| group.len() > 1)
            .collect();
        Some(SourceMultiplicityResolution {
            method,
            states: solution.states,
            equivalent_fragments,
            pairs,
            opposed_faces,
            transverse_crossings,
        })
    }
}

fn original_ids(resolution: &mut SourceMultiplicityResolution, original: &[usize]) {
    for witness in &mut resolution.transverse_crossings {
        witness.polygons = witness.polygons.map(|id| original[id]);
    }
    for witness in &mut resolution.opposed_faces {
        *witness = witness.map(|id| original[id]);
    }
    for pair in &mut resolution.pairs {
        pair.polygons = pair.polygons.map(|id| original[id]);
    }
    for fragment in resolution.equivalent_fragments.iter_mut().flatten() {
        for id in fragment {
            *id = original[*id];
        }
    }
}

pub(super) fn select_components(
    positions: &[Point3],
    polygons: &[InputPolygon],
    incidence: &EdgeIncidence,
    closed: &[bool],
    provenance: &mut SourceProvenance,
) -> Vec<Option<usize>> {
    if !incidence.values().any(|uses| uses.len() > 2) {
        return vec![None; polygons.len()];
    }
    let (membership, count) = connected_membership(polygons.len(), incidence, []);
    let mut fragments = vec![Vec::new(); count];
    for (polygon, &group) in membership.iter().enumerate() {
        fragments[group].push(polygon);
    }
    let mut edges = incidence
        .iter()
        .filter(|(_, uses)| !matches!(uses.as_slice(),[a,b] if a.1!=b.1))
        .map(|(&edge, uses)| AmbiguousEdge {
            edge,
            uses: uses
                .iter()
                .map(|&(polygon, forward)| SourceUse {
                    fragment: membership[polygon],
                    polygon,
                    forward,
                })
                .collect(),
        })
        .collect::<Vec<_>>();
    edges.sort_unstable_by_key(|e| e.edge);
    let mut resolutions = Vec::new();
    for cluster in clusters(count, &edges) {
        if let Some(input) = ClusterInput::new(&cluster, &fragments, &edges)
            && let Some(resolution) = input.resolve(positions, polygons)
        {
            resolutions.push((cluster, resolution));
        }
    }
    if resolutions.is_empty() {
        return vec![None; polygons.len()];
    }
    let candidates = closed_components(
        positions,
        polygons,
        incidence,
        resolutions
            .iter()
            .flat_map(|(_, r)| r.pairs.iter().map(|p| p.polygons)),
        None,
    );
    let mut selected = vec![None; polygons.len()];
    for (cluster, mut resolution) in resolutions {
        let ids = cluster
            .groups
            .iter()
            .flat_map(|&group| fragments[group].iter().copied());
        if ids.clone().any(|id| closed[id] || candidates[id].is_none())
            || !vertex_links::manifold(polygons, &candidates, ids.clone())
        {
            continue;
        }
        for id in ids {
            selected[id] = candidates[id];
        }
        original_ids(&mut resolution, &provenance.original_polygon_ids);
        provenance.multiplicity_resolutions.push(resolution);
    }
    selected
}

#[cfg(test)]
mod tests {
    use crate::{
        Point3,
        mesh_recovery::{AnalyticRecoveryOptions, analyze_mesh_evidence, recover_analytic_mesh},
    };
    use monstertruck_mesh::{Faces, PolygonMesh, StandardAttributes};
    use monstertruck_topology::shell::ShellCondition;

    fn pyramids(
        alternate: bool,
        reversed: bool,
        shuffle: bool,
    ) -> Result<PolygonMesh, Box<dyn std::error::Error>> {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.25, 0.25, 1.0),
            Point3::new(0.75, 0.75, 2.0),
        ];
        let mut faces = vec![
            [0, 2, 1],
            [0, 3, 2],
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ];
        faces.extend(if alternate {
            [[0, 3, 1], [1, 3, 2]]
        } else {
            [[0, 2, 1], [0, 3, 2]]
        });
        faces.extend([[0, 1, 5], [1, 2, 5], [2, 3, 5], [3, 0, 5]]);
        if reversed {
            for face in &mut faces {
                face.swap(1, 2);
            }
        }
        if shuffle {
            faces.rotate_left(5);
            faces.reverse();
        }
        Ok(PolygonMesh::new(
            StandardAttributes {
                positions,
                ..Default::default()
            },
            Faces::from_iter(faces),
        )?)
    }

    #[test]
    fn repeated_complete_oriented_source_chains_preserve_every_polygon()
    -> Result<(), Box<dyn std::error::Error>> {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        let single = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        for reversed in [false, true] {
            for copies in [2, 3, 5] {
                let mut faces = Vec::new();
                for (i, face) in single.iter().enumerate() {
                    for copy in 0..copies {
                        let mut cycle = *face;
                        cycle.rotate_left((i + copy) % 3);
                        if reversed {
                            cycle.reverse();
                        }
                        faces.push(cycle);
                    }
                }
                faces.rotate_left(3);
                let mesh = PolygonMesh::new(
                    StandardAttributes {
                        positions: positions.clone(),
                        ..Default::default()
                    },
                    Faces::from_iter(faces),
                )?;
                let evidence = analyze_mesh_evidence(&mesh, Default::default())?;
                assert_eq!(evidence.non_manifold_edge_count, 0);
                assert_eq!(evidence.provenance.original_polygon_ids.len(), copies * 4);
                let cert = &evidence.provenance.multiplicity_resolutions[0];
                assert_eq!(
                    cert.method,
                    super::SourceMultiplicityMethod::UniformOrientedCoverage
                );
                assert_eq!(cert.pairs.len(), copies * 6);
                assert!(evidence.provenance.canceled_sheets.is_empty());
                let volume = evidence.closed_source_volume().ok_or("closed source")?;
                let expected = f64::from(u32::try_from(copies)?) / 6.0;
                assert!((volume - if reversed { -expected } else { expected }).abs() < 1e-12);
                let recovered = recover_analytic_mesh(&mesh, &AnalyticRecoveryOptions::default())?;
                let shell = recovered
                    .complete_shell()
                    .ok_or("closed repeated coverage")?;
                assert_eq!(shell.connected_components().len(), copies);
                assert_eq!(shell.len(), copies * 4);
            }
        }
        Ok(())
    }

    #[test]
    fn coincident_oriented_caps_keep_distinct_body_cores() -> Result<(), Box<dyn std::error::Error>>
    {
        for alternate in [false, true] {
            for reversed in [false, true] {
                for shuffle in [false, true] {
                    let mesh = pyramids(alternate, reversed, shuffle)?;
                    let evidence = analyze_mesh_evidence(&mesh, Default::default())?;
                    assert_eq!(evidence.non_manifold_edge_count, 0);
                    assert_eq!(evidence.provenance.original_polygon_ids.len(), 12);
                    assert_eq!(evidence.provenance.multiplicity_resolutions.len(), 1);
                    assert_eq!(evidence.provenance.separated_vertices.len(), 4);
                    assert!(evidence.provenance.canceled_sheets.is_empty());
                    let cert = &evidence.provenance.multiplicity_resolutions[0];
                    assert_eq!(
                        cert.equivalent_fragments.len(),
                        if alternate { 1 } else { 2 }
                    );
                    assert_eq!(cert.pairs.len(), if alternate { 8 } else { 10 });
                    let volume = evidence.closed_source_volume().ok_or("closed source")?;
                    assert!((volume - if reversed { -1.0 } else { 1.0 }).abs() < 1e-12);
                    let recovered =
                        recover_analytic_mesh(&mesh, &AnalyticRecoveryOptions::default())?;
                    let shell = recovered
                        .complete_shell()
                        .ok_or("closed recovered pyramids")?;
                    assert_eq!(shell.shell_condition(), ShellCondition::Closed);
                    assert_eq!(shell.connected_components().len(), 2);
                    assert_eq!(shell.len(), 10);
                }
            }
        }
        Ok(())
    }
}
