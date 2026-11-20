//! Equal-radius three-way rolling-ball junctions between planar carriers.
//! The sharp vertex is the intersection of the three carrier planes, not a
//! projection of an arbitrary terminal point from one cylindrical strip.

use super::{
    CurveSupport, DIR_TOL, EdgeTreatmentEvidence, EdgeTreatmentVirtualEdge, EndpointEdit,
    GEOM_TOL_MM, KernelFaceSplice, KernelNewEdge, KernelPatchCurve, KernelPatchEdgeRef,
    KernelPatchEdgeUse, KernelPatchVertexRef, KernelTopologyPatch, ReducedFace, ReducedSolid,
    SurfaceSupport, TreatmentTopology, edge_face_users, face, normalize_retained_edge_rewrites,
    topology_edges, topology_face,
};
use crate::{
    brep::{PlaneSupport, SphereSupport},
    cad_ir::EdgeSelector,
    edge_treatments::{EdgeTreatmentEvidenceStrength, EdgeTreatmentKind},
    math3::{add, cross, distance, dot, mul, normalize, point_to_unit_line_distance, sub},
};
use anyhow::{Context, Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub struct SphericalJunction {
    pub face_id: u64,
    pub cylinders: [u64; 3],
    pub point: [f64; 3],
}

fn plane(solid: &ReducedSolid, id: u64) -> Result<PlaneSupport> {
    let SurfaceSupport::Plane(plane) = face(solid, id)?.surface else {
        anyhow::bail!("junction carrier #{id} is not planar");
    };
    Ok(plane)
}

fn triple_intersection(planes: [PlaneSupport; 3]) -> Result<[f64; 3]> {
    let [a, b, c] = planes;
    let [an, bn, cn] = [a.normal, b.normal, c.normal].map(|normal| normalize(normal, DIR_TOL));
    let (an, bn, cn) = (
        an.context("degenerate plane")?,
        bn.context("degenerate plane")?,
        cn.context("degenerate plane")?,
    );
    let determinant = dot(an, cross(bn, cn));
    ensure!(
        determinant.abs() > DIR_TOL,
        "junction planes do not have a unique intersection"
    );
    // Work relative to a carrier origin, avoiding cancellation for remote parts.
    let point = add(
        a.origin_mm,
        mul(
            add(
                mul(cross(cn, an), dot(bn, sub(b.origin_mm, a.origin_mm))),
                mul(cross(an, bn), dot(cn, sub(c.origin_mm, a.origin_mm))),
            ),
            determinant.recip(),
        ),
    );
    ensure!(
        planes
            .iter()
            .all(|p| dot(p.normal, sub(point, p.origin_mm)).abs() <= GEOM_TOL_MM),
        "sharp junction does not lie on all carrier planes"
    );
    Ok(point)
}

fn is_round(evidence: &EdgeTreatmentEvidence) -> bool {
    evidence.kind == EdgeTreatmentKind::CylindricalFillet
        && evidence.strength == EdgeTreatmentEvidenceStrength::AnalyticTangency
        && evidence.carrier_pair_unique
        && matches!(
            evidence.virtual_edge,
            EdgeTreatmentVirtualEdge::Segment { .. }
        )
}

type EvidenceIndex<'a> = BTreeMap<u64, &'a EdgeTreatmentEvidence>;

fn certify_sphere(
    solid: &ReducedSolid,
    sphere: &ReducedFace,
    support: SphereSupport,
    candidates: &EvidenceIndex<'_>,
) -> Result<SphericalJunction> {
    let [wire] = sphere.loops.as_slice() else {
        anyhow::bail!("sphere junction has multiple boundaries");
    };
    let mut neighbors = BTreeSet::new();
    for edge in &wire.edges {
        let users = edge_face_users(solid, edge.edge_id);
        ensure!(
            users.len() == 2 && users.contains(&sphere.id),
            "nonmanifold sphere boundary"
        );
        let neighbor = *users
            .iter()
            .find(|&&id| id != sphere.id)
            .context("sphere boundary lacks second face")?;
        ensure!(
            candidates.contains_key(&neighbor),
            "sphere boundary is not an exact round interface"
        );
        let SurfaceSupport::Cylinder(cylinder) = face(solid, neighbor)?.surface else {
            anyhow::bail!("junction neighbor is not cylindrical");
        };
        let axis = normalize(cylinder.axis, DIR_TOL).context("degenerate cylinder axis")?;
        let CurveSupport::Circle(circle) = edge.support else {
            anyhow::bail!("sphere/cylinder seam is not a circle");
        };
        ensure!(
            (cylinder.radius_mm - support.radius_mm).abs() <= GEOM_TOL_MM
                && (circle.radius_mm - support.radius_mm).abs() <= GEOM_TOL_MM
                && distance(circle.center_mm, support.center_mm) <= GEOM_TOL_MM
                && dot(circle.normal, axis).abs() >= 1.0 - DIR_TOL
                && point_to_unit_line_distance(support.center_mm, cylinder.axis_origin_mm, axis)
                    <= GEOM_TOL_MM,
            "sphere/cylinder seam does not certify a common rolling ball"
        );
        neighbors.insert(neighbor);
    }
    let cylinders: [u64; 3] = neighbors
        .into_iter()
        .collect::<Vec<_>>()
        .try_into()
        .map_err(|_| anyhow::anyhow!("sphere does not join exactly three rounds"))?;
    let mut pairs = BTreeSet::new();
    let mut carrier_ids = BTreeSet::new();
    for id in cylinders {
        let mut pair = candidates[&id]
            .carriers
            .each_ref()
            .map(|carrier| carrier.face_id);
        pair.sort_unstable();
        pairs.insert(pair);
        carrier_ids.extend(pair);
    }
    ensure!(
        pairs.len() == 3 && carrier_ids.len() == 3,
        "sphere does not join three distinct pairs of three planes"
    );
    let planes = carrier_ids
        .into_iter()
        .map(|id| plane(solid, id))
        .collect::<Result<Vec<_>>>()?;
    let planes: [PlaneSupport; 3] = planes
        .try_into()
        .map_err(|_| anyhow::anyhow!("missing junction planes"))?;
    ensure!(
        planes.iter().all(
            |p| (dot(p.normal, sub(support.center_mm, p.origin_mm)).abs() - support.radius_mm)
                .abs()
                <= GEOM_TOL_MM
        ),
        "sphere is not tangent to every sharp carrier plane"
    );
    let point = triple_intersection(planes)?;
    Ok(SphericalJunction {
        face_id: sphere.id,
        cylinders,
        point,
    })
}

/// A connected component of exact three-way rolling-ball junctions and
/// their cylindrical seams. A component may have other planar endpoints.
pub struct RoundNetwork<'a> {
    pub rounds: Vec<&'a EdgeTreatmentEvidence>,
    pub junctions: Vec<SphericalJunction>,
}

pub fn round_networks<'a>(
    solid: &ReducedSolid,
    evidence: &'a [EdgeTreatmentEvidence],
) -> Vec<RoundNetwork<'a>> {
    let candidates = evidence
        .iter()
        .filter(|e| is_round(e))
        .map(|e| (e.treatment_face_id, e))
        .collect::<EvidenceIndex<'_>>();
    let junctions = solid
        .faces()
        .filter_map(|face| {
            let SurfaceSupport::Sphere(support) = face.surface else {
                return None;
            };
            certify_sphere(solid, face, support, &candidates).ok()
        })
        .collect::<Vec<_>>();
    let mut by_round = BTreeMap::<u64, Vec<usize>>::new();
    for (i, junction) in junctions.iter().enumerate() {
        for cylinder in junction.cylinders {
            by_round.entry(cylinder).or_default().push(i);
        }
    }
    let mut seen = BTreeSet::new();
    let mut junctions = junctions.into_iter().map(Some).collect::<Vec<_>>();
    let mut networks = Vec::new();
    for start in 0..junctions.len() {
        if !seen.insert(start) {
            continue;
        }
        let mut stack = vec![start];
        let mut members = Vec::new();
        let mut rounds = BTreeSet::new();
        while let Some(index) = stack.pop() {
            let Some(junction) = junctions[index].take() else {
                continue;
            };
            for id in junction.cylinders {
                rounds.insert(id);
                for &next in &by_round[&id] {
                    if seen.insert(next) {
                        stack.push(next);
                    }
                }
            }
            members.push(junction);
        }
        networks.push(RoundNetwork {
            rounds: rounds.into_iter().map(|id| candidates[&id]).collect(),
            junctions: members,
        });
    }
    networks
}

struct Terminal {
    face_id: u64,
    edges: BTreeSet<u64>,
    vertices: BTreeSet<u64>,
    point: [f64; 3],
}

fn terminal_runs(
    solid: &ReducedSolid,
    item: &EdgeTreatmentEvidence,
    junctions: &BTreeMap<u64, [f64; 3]>,
) -> Result<Vec<Terminal>> {
    let EdgeTreatmentVirtualEdge::Segment { start_mm, end_mm } = item.virtual_edge else {
        anyhow::bail!("round lacks a sharp segment");
    };
    let source = face(solid, item.treatment_face_id)?;
    let interfaces = item
        .carriers
        .iter()
        .flat_map(|c| c.shared_edge_ids.iter().copied())
        .collect::<BTreeSet<_>>();
    let mut terminal = BTreeMap::<u64, (BTreeSet<u64>, BTreeSet<u64>)>::new();
    for edge in source
        .loop_edges()
        .filter(|e| !interfaces.contains(&e.edge_id))
    {
        let users = edge_face_users(solid, edge.edge_id);
        ensure!(
            users.len() == 2 && users.contains(&source.id),
            "nonmanifold round terminal"
        );
        let id = *users
            .iter()
            .find(|&&id| id != source.id)
            .context("no terminal neighbor")?;
        let (edges, vertices) = terminal.entry(id).or_default();
        edges.insert(edge.edge_id);
        vertices.extend([edge.start_vertex, edge.end_vertex]);
    }
    ensure!(
        terminal.len() == 2,
        "round must have two distinct terminal runs"
    );
    terminal
        .into_iter()
        .map(|(id, (edges, vertices))| {
            let point = if let Some(point) = junctions.get(&id) {
                *point
            } else {
                let plane = plane(solid, id)?;
                // The cap is a surviving plane, not another inferred feature.
                let ends = [start_mm, end_mm]
                    .into_iter()
                    .filter(|p| dot(plane.normal, sub(*p, plane.origin_mm)).abs() <= GEOM_TOL_MM)
                    .collect::<Vec<_>>();
                let [point] = ends.as_slice() else {
                    anyhow::bail!("round terminal does not have one sharp endpoint on its cap");
                };
                *point
            };
            let axis = normalize(sub(end_mm, start_mm), DIR_TOL).context("zero inverse segment")?;
            ensure!(
                point_to_unit_line_distance(point, start_mm, axis) <= GEOM_TOL_MM,
                "sharp junction not on inverse line"
            );
            Ok(Terminal {
                face_id: id,
                edges,
                vertices,
                point,
            })
        })
        .collect()
}

#[derive(Default)]
struct VertexMap {
    points: Vec<[f64; 3]>,
    old: BTreeMap<u64, usize>,
    sphere: BTreeMap<u64, usize>,
}

impl VertexMap {
    fn terminal(&mut self, terminal: &Terminal, is_sphere: bool) -> Result<usize> {
        let existing = terminal
            .vertices
            .iter()
            .filter_map(|id| self.old.get(id).copied())
            .collect::<BTreeSet<_>>();
        ensure!(
            existing.len() <= 1,
            "round terminal has conflicting sharp vertex owners"
        );
        let index = if let Some(index) = self
            .sphere
            .get(&terminal.face_id)
            .copied()
            .filter(|_| is_sphere)
        {
            ensure!(
                existing.iter().all(|other| *other == index),
                "sphere joins different sharp vertices"
            );
            index
        } else if let Some(index) = existing.into_iter().next() {
            index
        } else {
            let index = self.points.len();
            self.points.push(terminal.point);
            index
        };
        ensure!(
            distance(self.points[index], terminal.point) <= GEOM_TOL_MM,
            "inconsistent shared terminal point"
        );
        if is_sphere {
            self.sphere.insert(terminal.face_id, index);
        }
        for &id in &terminal.vertices {
            self.old.insert(id, index);
        }
        Ok(index)
    }
}

pub struct PlannedRoundNetwork {
    pub patch: KernelTopologyPatch,
    pub selectors: Vec<(u64, EdgeSelector)>,
    pub junction_faces: Vec<u64>,
}

fn carrier_splice(
    solid: &ReducedSolid,
    carrier: &crate::edge_treatments::EdgeTreatmentCarrier,
    map: &VertexMap,
    indices: [usize; 2],
    edge_index: usize,
    section: usize,
    provenance: &impl TreatmentTopology,
) -> Result<KernelFaceSplice> {
    let face = face(solid, carrier.face_id)?;
    let shared = carrier
        .shared_edge_ids
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let (previous, next) = super::carrier_run_neighbors(face, &shared)?;
    let first = map
        .old
        .get(&previous.end_vertex)
        .copied()
        .context("carrier start not on round terminal")?;
    let last = map
        .old
        .get(&next.start_vertex)
        .copied()
        .context("carrier end not on round terminal")?;
    ensure!(
        [first, last] == indices || [last, first] == indices,
        "carrier run does not join both sharp vertices"
    );
    for index in indices {
        let p = plane(solid, face.id)?;
        ensure!(
            dot(p.normal, sub(map.points[index], p.origin_mm)).abs() <= GEOM_TOL_MM,
            "sharp rail outside carrier plane"
        );
    }
    Ok(KernelFaceSplice {
        face: topology_face(provenance, section, face.id)?,
        remove_edges: topology_edges(provenance, section, shared)?,
        replacement: vec![KernelPatchEdgeUse {
            edge: KernelPatchEdgeRef::New(edge_index),
            orientation: first == indices[0],
        }],
    })
}

fn retained_edits(
    solid: &ReducedSolid,
    map: &VertexMap,
    removed_edges: &BTreeSet<u64>,
    removed_faces: &BTreeSet<u64>,
) -> Result<Vec<EndpointEdit>> {
    let mut edits = Vec::new();
    for face in solid.faces().filter(|f| !removed_faces.contains(&f.id)) {
        for edge in face
            .loop_edges()
            .filter(|e| !removed_edges.contains(&e.edge_id))
        {
            for vertex in [edge.start_vertex, edge.end_vertex] {
                let Some(index) = map.old.get(&vertex).copied() else {
                    continue;
                };
                let CurveSupport::Line(line) = edge.support else {
                    anyhow::bail!("junction would move a retained nonlinear edge");
                };
                let point = map.points[index];
                ensure!(
                    point_to_unit_line_distance(
                        point,
                        line.origin_mm,
                        normalize(line.direction, DIR_TOL).context("degenerate retained line")?
                    ) <= GEOM_TOL_MM,
                    "sharp junction is outside surviving line support"
                );
                edits.push(EndpointEdit {
                    carrier_face_id: face.id,
                    edge: edge.clone(),
                    old_vertex: vertex,
                    new_vertex: index,
                    new_point: point,
                    replacement_curve: None,
                });
            }
        }
    }
    Ok(edits)
}

/// Invert a complete connected rolling-ball network. Every sphere is removed
/// with all three adjacent rounds; all newly shortened or extended retained
/// line endpoints retain their original geometric support. No coordinate weld.
pub fn plan_round_network(
    solid: &ReducedSolid,
    network: &RoundNetwork<'_>,
    section: usize,
    provenance: &impl TreatmentTopology,
) -> Result<PlannedRoundNetwork> {
    let junctions = network
        .junctions
        .iter()
        .map(|j| (j.face_id, j.point))
        .collect::<BTreeMap<_, _>>();
    let junction_faces = junctions.keys().copied().collect::<Vec<_>>();
    let removed_faces = network
        .rounds
        .iter()
        .map(|e| e.treatment_face_id)
        .chain(junction_faces.iter().copied())
        .collect::<BTreeSet<_>>();
    let mut removed_edges = BTreeSet::new();
    for &id in &removed_faces {
        removed_edges.extend(super::unique_face_edges(face(solid, id)?)?);
    }
    let mut map = VertexMap::default();
    let mut new_edges = Vec::new();
    let mut splices = Vec::new();
    let mut selectors = Vec::new();
    for item in &network.rounds {
        let terminals = terminal_runs(solid, item, &junctions)?;
        let [first, last] = terminals.as_slice() else {
            anyhow::bail!("round has wrong number of terminal runs");
        };
        let indices = [
            map.terminal(first, junctions.contains_key(&first.face_id))?,
            map.terminal(last, junctions.contains_key(&last.face_id))?,
        ];
        let [start_mm, end_mm] = indices.map(|i| map.points[i]);
        ensure!(
            distance(start_mm, end_mm) > GEOM_TOL_MM,
            "inverse sharp edge degenerates"
        );
        let edge_index = new_edges.len();
        new_edges.push(KernelNewEdge {
            vertices: (
                KernelPatchVertexRef::New(indices[0]),
                KernelPatchVertexRef::New(indices[1]),
            ),
            curve: KernelPatchCurve::Line { start_mm, end_mm },
        });
        selectors.push((
            item.treatment_face_id,
            EdgeSelector::Segment { start_mm, end_mm },
        ));
        // Every original terminal run belongs either to a removed sphere or
        // a surviving planar cap. A cap run collapses to one sharp vertex.
        for terminal in &terminals {
            if !removed_faces.contains(&terminal.face_id) {
                ensure!(
                    super::edges_form_one_cyclic_run(
                        face(solid, terminal.face_id)?,
                        &terminal.edges
                    ),
                    "terminal cap run not contiguous"
                );
                splices.push(KernelFaceSplice {
                    face: topology_face(provenance, section, terminal.face_id)?,
                    remove_edges: topology_edges(
                        provenance,
                        section,
                        terminal.edges.iter().copied(),
                    )?,
                    replacement: Vec::new(),
                });
            }
        }
        for carrier in &item.carriers {
            ensure!(
                !removed_faces.contains(&carrier.face_id),
                "round depends on another removed carrier"
            );
            splices.push(carrier_splice(
                solid, carrier, &map, indices, edge_index, section, provenance,
            )?);
        }
    }
    // Whole network must assign every shared old corner consistently.
    let edits = retained_edits(solid, &map, &removed_edges, &removed_faces)?;
    let patch = KernelTopologyPatch {
        remove_faces: removed_faces
            .iter()
            .map(|&id| topology_face(provenance, section, id))
            .collect::<Result<_>>()?,
        new_vertices: map.points,
        new_edges,
        rewrite_edges: normalize_retained_edge_rewrites(&edits, section, provenance)?,
        face_splices: splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(provenance, section, removed_edges)?,
    };
    Ok(PlannedRoundNetwork {
        patch,
        selectors,
        junction_faces,
    })
}
