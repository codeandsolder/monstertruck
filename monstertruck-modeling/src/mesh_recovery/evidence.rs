//! Primitive-independent geometric evidence from polygon mesh incidence.
//!
//! Source incidence and collinear chains are topology/coordinate facts.
//! Dihedral selection and parallel groups are hypotheses with explicit
//! thresholds; they must not alter the source vertices or edge graph.

use super::{InputPolygon, MeshRecoveryError, edge_incidence, positions_and_polygons, seams};
use crate::{InnerSpace, KERNEL_EPSILON, Point3, Vector3};
use monstertruck_mesh::PolygonMesh;
use std::collections::{HashMap, HashSet};

/// Thresholds only for feature discovery and direction clustering.
#[derive(Clone, Copy, Debug)]
pub struct EvidenceOptions {
    /// Source edge is retained when its normal discontinuity exceeds this
    /// angle in radians. Coplanar triangulation edges are suppressed.
    pub min_dihedral_angle_radians: f64,
    /// Undirected straight chains differing by at most this angle in radians
    /// may join a candidate parallel direction group.
    pub max_parallel_angle_radians: f64,
}

impl Default for EvidenceOptions {
    fn default() -> Self {
        Self {
            min_dihedral_angle_radians: 1.0e-6,
            max_parallel_angle_radians: 1.0e-6,
        }
    }
}

/// Reason a source segment is geometrically interesting.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FeatureEdgeKind {
    /// Source surface boundary (one incident polygon).
    Boundary,
    /// Two polygons have measurably different normals. This is not yet
    /// evidence of a designed sharp edge: it can be smooth tessellation.
    Dihedral,
}

/// Source edge with canonical positions and original polygon provenance.
#[derive(Clone, Debug)]
pub struct FeatureEdge {
    /// Undirected canonical source position indices.
    pub vertices: [usize; 2],
    /// Source polygon indices on each side of the edge.
    pub source_face_indices: Vec<usize>,
    /// Boundary or measurable dihedral.
    pub kind: FeatureEdgeKind,
    /// Dihedral angle in radians; none for source boundaries.
    pub dihedral_angle_radians: Option<f64>,
}

/// Maximal exact straight-through feature chain.
#[derive(Clone, Debug)]
pub struct CollinearChain {
    /// Ordered canonical vertices, retaining all original intermediate points.
    pub vertices: Vec<usize>,
    /// Indices in the feature-edge array between successive chain vertices.
    pub edge_indices: Vec<usize>,
    /// Category shared across the chain.
    pub kind: FeatureEdgeKind,
}

/// Relative line placement in one candidate parallel bundle's seed frame.
#[derive(Clone, Debug)]
pub struct ParallelPlacement {
    /// Index in [`MeshEvidence::chains`].
    pub chain_index: usize,
    /// Perpendicular displacement from seed line to this chain's first point.
    /// For approximately parallel directions this is a pointwise offset,
    /// not a globally constant separation.
    pub transverse_offset: Vector3,
    /// Magnitude of `transverse_offset` in mesh coordinate units.
    pub transverse_distance: f64,
    /// Projected start/end extent in the seed line's local frame.
    pub longitudinal_interval: [f64; 2],
    /// Angular difference from the *fixed seed direction*, in radians.
    pub angular_residual_radians: f64,
}

/// Angular hypothesis with spatial evidence for each straight-line member.
///
/// Every candidate is checked against the seed, never against an arbitrary
/// neighboring member; tolerance chains cannot drift through orientation
/// space and incorrectly merge nonparallel carriers.
#[derive(Clone, Debug)]
pub struct ParallelBundle {
    /// Fixed, canonical undirected unit seed axis.
    pub direction: Vector3,
    /// First point on the seed chain, providing a local reference frame.
    pub origin: Point3,
    /// Source chains and their relative placement measurements.
    pub placements: Vec<ParallelPlacement>,
}

impl ParallelBundle {
    /// Transverse separation of two members in the bundle's seed plane.
    /// For nearly parallel candidates, it is only a local measurement.
    #[must_use]
    pub fn transverse_separation(&self, first: usize, second: usize) -> Option<f64> {
        let a = self.placements.get(first)?.transverse_offset;
        let b = self.placements.get(second)?.transverse_offset;
        Some((a - b).magnitude())
    }
}

/// Relationship between two straight chain support lines, independent of a
/// high-level surface or CAD feature hypothesis.
#[derive(Clone, Debug)]
pub struct ChainPairEvidence {
    /// Undirected angle in radians, from parallel (zero) to perpendicular.
    pub angle_radians: f64,
    /// Adaptive exact 3D orientation predicate: all four stored source
    /// endpoints belong to a common supporting plane.
    pub exactly_coplanar: bool,
    /// Both source endpoints of the second chain lie on the first infinite
    /// line, according to adaptive exact 2D projected orientation predicates.
    pub exactly_coincident_carriers: bool,
    /// Closest distance between infinite support lines (not finite segments).
    /// For nearly parallel lines this may occur far outside either segment.
    pub supporting_line_distance: f64,
    /// Separation of the second line's first point from the first line;
    /// preserves a useful local spatial measurement for near-parallel lines.
    pub reference_transverse_distance: f64,
}

/// Reusable primitive-free geometric evidence for a source mesh.
#[derive(Debug)]
pub struct MeshEvidence {
    /// Exactly welded source positions.
    pub positions: Vec<Point3>,
    /// Boundary and dihedral source edges in canonical index order.
    pub edges: Vec<FeatureEdge>,
    /// Maximal exact collinear straight-through chains.
    pub chains: Vec<CollinearChain>,
    /// Approximate direction groups with measured spatial relationships.
    pub parallel_bundles: Vec<ParallelBundle>,
    /// Locally coplanar triangulation edges not included as features.
    pub suppressed_interior_edge_count: usize,
}

impl MeshEvidence {
    /// Compare the infinite support lines of any two recovered feature
    /// chains, computing geometric relations lazily rather than materializing
    /// a quadratic matrix of all chain pairs.
    #[must_use]
    pub fn chain_pair(&self, first: usize, second: usize) -> Option<ChainPairEvidence> {
        let a = self.chains.get(first)?;
        let b = self.chains.get(second)?;
        let [a0, a1] = [
            self.positions[*a.vertices.first()?],
            self.positions[*a.vertices.last()?],
        ];
        let [b0, b1] = [
            self.positions[*b.vertices.first()?],
            self.positions[*b.vertices.last()?],
        ];
        let u = (a1 - a0).normalize();
        let v = (b1 - b0).normalize();
        let cross = u.cross(v);
        let sine = cross.magnitude();
        let angle = sine.atan2(u.dot(v).abs());
        let delta = b0 - a0;
        let reference_distance = delta.cross(u).magnitude();
        let supporting_distance = if sine == 0.0 {
            reference_distance
        } else {
            delta.dot(cross).abs() / sine
        };
        let coord = |p: Point3| robust::Coord3D {
            x: p.x,
            y: p.y,
            z: p.z,
        };
        Some(ChainPairEvidence {
            angle_radians: angle,
            exactly_coplanar: robust::orient3d(coord(a0), coord(a1), coord(b0), coord(b1)) == 0.0,
            exactly_coincident_carriers: seams::on_line(a0, a1, b0) && seams::on_line(a0, a1, b1),
            supporting_line_distance: supporting_distance,
            reference_transverse_distance: reference_distance,
        })
    }
}

fn select_edges(
    polygons: &[InputPolygon],
    options: EvidenceOptions,
) -> Result<(Vec<FeatureEdge>, usize), MeshRecoveryError> {
    let (incidence, _) = edge_incidence(polygons)?;
    let mut keys = incidence.keys().copied().collect::<Vec<_>>();
    keys.sort_unstable();
    let mut edges = Vec::new();
    let mut suppressed = 0;
    for key in keys {
        let uses = &incidence[&key];
        let (kind, angle) = match uses.as_slice() {
            [(_face, _)] => (Some(FeatureEdgeKind::Boundary), None),
            [(left, _), (right, _)] => {
                let a = polygons[*left].normal;
                let b = polygons[*right].normal;
                let angle = a.cross(b).magnitude().atan2(a.dot(b));
                if angle > options.min_dihedral_angle_radians {
                    (Some(FeatureEdgeKind::Dihedral), Some(angle))
                } else {
                    suppressed += 1;
                    (None, None)
                }
            }
            _ => return Err(MeshRecoveryError::NonManifoldEdge(key)),
        };
        if let Some(kind) = kind {
            edges.push(FeatureEdge {
                vertices: key,
                source_face_indices: uses.iter().map(|(face, _)| *face).collect(),
                kind,
                dihedral_angle_radians: angle,
            });
        }
    }
    Ok((edges, suppressed))
}

/// Require a degree-two feature vertex, same edge category, and the vertex
/// strictly between both other endpoints on one exact source line.
fn straight_through(positions: &[Point3], edges: &[FeatureEdge]) -> (Vec<Vec<usize>>, Vec<bool>) {
    let mut incident = vec![Vec::new(); positions.len()];
    for (id, edge) in edges.iter().enumerate() {
        for &vertex in &edge.vertices {
            incident[vertex].push(id);
        }
    }
    let through = incident
        .iter()
        .enumerate()
        .map(|(vertex, neighbors)| {
            let [a, b] = neighbors.as_slice() else {
                return false;
            };
            if edges[*a].kind != edges[*b].kind {
                return false;
            }
            let other = |edge: &FeatureEdge| {
                if edge.vertices[0] == vertex {
                    edge.vertices[1]
                } else {
                    edge.vertices[0]
                }
            };
            seams::on_segment(
                positions[other(&edges[*a])],
                positions[other(&edges[*b])],
                positions[vertex],
            )
            .is_some()
        })
        .collect();
    (incident, through)
}

fn trace_chains(positions: &[Point3], edges: &[FeatureEdge]) -> Vec<CollinearChain> {
    let (incident, through) = straight_through(positions, edges);
    let mut visited = vec![false; edges.len()];
    let mut chains = Vec::new();
    for (first_id, edge) in edges.iter().enumerate() {
        if visited[first_id] {
            continue;
        }
        let [a, b] = edge.vertices;
        let start = if through[a] && !through[b] { b } else { a };
        let mut vertices = vec![start];
        let mut ids = Vec::new();
        let mut current = start;
        let mut next_id = first_id;
        loop {
            if visited[next_id] {
                break;
            }
            visited[next_id] = true;
            ids.push(next_id);
            let segment = edges[next_id].vertices;
            current = if segment[0] == current {
                segment[1]
            } else {
                segment[0]
            };
            vertices.push(current);
            if !through[current] {
                break;
            }
            let Some(&candidate) = incident[current]
                .iter()
                .find(|&&id| id != next_id && !visited[id])
            else {
                break;
            };
            next_id = candidate;
        }
        // Certify the entire chord, not only locally straight junctions.
        let first = positions[vertices[0]];
        let last = positions[vertices[vertices.len() - 1]];
        let exact = first != last
            && vertices[1..vertices.len() - 1]
                .iter()
                .all(|&v| seams::on_segment(first, last, positions[v]).is_some());
        if exact {
            chains.push(CollinearChain {
                vertices,
                edge_indices: ids,
                kind: edge.kind,
            });
        } else {
            for id in ids {
                chains.push(CollinearChain {
                    vertices: edges[id].vertices.to_vec(),
                    edge_indices: vec![id],
                    kind: edges[id].kind,
                });
            }
        }
    }
    chains
}

fn canonical_unit(vector: Vector3) -> Vector3 {
    let axis = (0..3)
        .max_by(|&a, &b| vector[a].abs().total_cmp(&vector[b].abs()))
        .unwrap_or(0);
    let directed = if vector[axis].is_sign_negative() {
        -vector
    } else {
        vector
    };
    directed.normalize()
}

// Angular-neighbor grid bounds the number of bundle seeds inspected per
// chain, even when a complex curved mesh has many unrelated directions.
// Register both +/- unit direction to avoid hemisphere-boundary artifacts.
fn direction_cell(direction: Vector3, width: f64) -> [f64; 3] {
    // Minimum cell width keeps indices below 2^53: integer cell values and
    // their immediate neighbors are exactly representable, without casts.
    [direction.x, direction.y, direction.z].map(|component| ((component + 1.0) / width).floor())
}

fn parallel_bundles(
    positions: &[Point3],
    chains: &[CollinearChain],
    options: EvidenceOptions,
) -> Vec<ParallelBundle> {
    let mut candidates = chains
        .iter()
        .enumerate()
        .map(|(id, chain)| {
            let delta =
                positions[chain.vertices[chain.vertices.len() - 1]] - positions[chain.vertices[0]];
            (id, delta.magnitude2(), canonical_unit(delta))
        })
        .collect::<Vec<_>>();
    candidates.sort_unstable_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    let width = (2.0 * (options.max_parallel_angle_radians * 0.5).sin()).max(f64::EPSILON * 8.0);
    let mut bins = HashMap::<[u64; 3], Vec<usize>>::new();
    let mut bundles: Vec<ParallelBundle> = Vec::new();
    for (chain_index, _, direction) in candidates {
        let start = positions[chains[chain_index].vertices[0]];
        let end = positions[chains[chain_index].vertices[chains[chain_index].vertices.len() - 1]];
        let mut possible = HashSet::new();
        for sign in [direction, -direction] {
            let cell = direction_cell(sign, width);
            for x in [-1.0, 0.0, 1.0] {
                for y in [-1.0, 0.0, 1.0] {
                    for z in [-1.0, 0.0, 1.0] {
                        let neighbor = [
                            (cell[0] + x).to_bits(),
                            (cell[1] + y).to_bits(),
                            (cell[2] + z).to_bits(),
                        ];
                        if let Some(ids) = bins.get(&neighbor) {
                            possible.extend(ids.iter().copied());
                        }
                    }
                }
            }
        }
        let closest = possible
            .into_iter()
            .filter_map(|id| {
                let dot = bundles[id].direction.dot(direction);
                let angle = bundles[id]
                    .direction
                    .cross(direction)
                    .magnitude()
                    .atan2(dot.abs());
                (angle <= options.max_parallel_angle_radians).then_some((id, angle))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        let index = if let Some((id, _)) = closest {
            id
        } else {
            let id = bundles.len();
            bundles.push(ParallelBundle {
                direction,
                origin: start,
                placements: Vec::new(),
            });
            for sign in [direction, -direction] {
                bins.entry(direction_cell(sign, width).map(f64::to_bits))
                    .or_default()
                    .push(id);
            }
            id
        };
        let group = &mut bundles[index];
        let delta = start - group.origin;
        let offset = delta - group.direction * delta.dot(group.direction);
        let begin = delta.dot(group.direction);
        let finish = (end - group.origin).dot(group.direction);
        let angular = group
            .direction
            .cross(direction)
            .magnitude()
            .atan2(group.direction.dot(direction).abs());
        group.placements.push(ParallelPlacement {
            chain_index,
            transverse_offset: offset,
            transverse_distance: offset.magnitude(),
            longitudinal_interval: [begin.min(finish), begin.max(finish)],
            angular_residual_radians: angular,
        });
    }
    bundles
}

/// Extract maximal exact straight feature chains and candidate parallel
/// spatial arrangements without presupposing any analytic CAD primitive.
///
/// Source coordinates are exactly welded; T-junctions are repaired using
/// adaptive exact predicates; coplanar triangulation edges are suppressed.
/// Other source edges retain polygon adjacency and measured dihedrals.
/// Chain collinearity is exact in the stored coordinates, while direction
/// groups are only bounded-angular hypotheses.
///
/// # Errors
/// Returns [`MeshRecoveryError`] for invalid thresholds or invalid source
/// geometry, orientation or manifold edge incidence.
pub fn analyze_mesh_evidence(
    mesh: &PolygonMesh,
    options: EvidenceOptions,
) -> Result<MeshEvidence, MeshRecoveryError> {
    if !options.min_dihedral_angle_radians.is_finite()
        || !(0.0..std::f64::consts::PI).contains(&options.min_dihedral_angle_radians)
        || !options.max_parallel_angle_radians.is_finite()
        || !(0.0..std::f64::consts::FRAC_PI_2).contains(&options.max_parallel_angle_radians)
    {
        return Err(MeshRecoveryError::InvalidEvidenceOptions);
    }
    let (positions, mut polygons) = positions_and_polygons(mesh, KERNEL_EPSILON)?;
    seams::split_source_t_junctions(&positions, &mut polygons);
    let (edges, suppressed_interior_edge_count) = select_edges(&polygons, options)?;
    let chains = trace_chains(&positions, &edges);
    let parallel_bundles = parallel_bundles(&positions, &chains, options);
    Ok(MeshEvidence {
        positions,
        edges,
        chains,
        parallel_bundles,
        suppressed_interior_edge_count,
    })
}

#[cfg(test)]
mod tests;
