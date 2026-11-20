//! Invert a closed perimeter of terminal-connected analytic treatments.
//!
//! The sharp junction is owned by its original shared terminal edge, rather
//! than inferred independently from neighbors that are themselves removed.
//! Each of the two original carrier runs must certify that same junction.

use super::{
    CircleFrame, DIR_TOL, EndpointEdit, GEOM_TOL_MM, carrier_run_neighbors, circle_carrier_run,
    circle_carrier_supports_inverse, circle_line_endpoint_edit, edge_face_users, face,
    normalize_retained_edge_rewrites, open_inverse_circle_arc, same_virtual_circle, topology_edges,
    topology_face, unique_face_edges,
};
use crate::{
    brep::{CurveSupport, OrientedEdgeUse, SurfaceSupport},
    edge_treatments::{
        EdgeTreatmentEvidence, EdgeTreatmentEvidenceStrength, EdgeTreatmentVirtualEdge,
    },
    kernel_source_load::KernelTopologyProvenance,
    kernel_topology::{
        KernelFaceSplice, KernelNewEdge, KernelPatchCurve, KernelPatchEdgeRef, KernelPatchEdgeUse,
        KernelPatchVertexRef, KernelTopologyPatch,
    },
    math3::{distance, dot, norm, normalize, point_to_unit_line_distance, sub},
    reduced_brep::{ReducedFace, ReducedSolid},
};
use anyhow::{Context, Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

type EvidenceGroup<'a> = Vec<&'a EdgeTreatmentEvidence>;

fn terminal_edges(solid: &ReducedSolid, evidence: &EdgeTreatmentEvidence) -> Result<BTreeSet<u64>> {
    let edges = unique_face_edges(face(solid, evidence.treatment_face_id)?)?;
    let interfaces = evidence
        .carriers
        .iter()
        .flat_map(|c| c.shared_edge_ids.iter().copied())
        .collect::<BTreeSet<_>>();
    ensure!(
        interfaces.is_subset(&edges),
        "treatment has foreign carrier edges"
    );
    Ok(edges.difference(&interfaces).copied().collect())
}

/// Identify closed terminal-connected components independently of source ID
/// order. Carrier/treatment adjacency is a different relation: a face used
/// as support cannot be removed in the same inverse layer.
pub fn closed_components<'a>(
    solid: &ReducedSolid,
    evidence: &[&'a EdgeTreatmentEvidence],
) -> (Vec<EvidenceGroup<'a>>, EvidenceGroup<'a>) {
    let terminals = evidence
        .iter()
        .map(|e| terminal_edges(solid, e).ok())
        .collect::<Vec<_>>();
    let mut owners = BTreeMap::<u64, Vec<usize>>::new();
    for (i, edges) in terminals.iter().enumerate() {
        if let Some(edges) = edges {
            for edge in edges {
                owners.entry(*edge).or_default().push(i);
            }
        }
    }
    let mut seen = BTreeSet::new();
    let mut closed = Vec::new();
    let mut other = Vec::new();
    for i in 0..evidence.len() {
        if !seen.insert(i) {
            continue;
        }
        let mut stack = vec![i];
        let mut members = vec![i];
        while let Some(current) = stack.pop() {
            if let Some(edges) = &terminals[current] {
                for neighbor in edges.iter().flat_map(|id| &owners[id]) {
                    if seen.insert(*neighbor) {
                        stack.push(*neighbor);
                        members.push(*neighbor);
                    }
                }
            }
        }
        members.sort_unstable();
        let group = members.iter().map(|&i| evidence[i]).collect::<Vec<_>>();
        let group_ids = group
            .iter()
            .map(|e| e.treatment_face_id)
            .collect::<BTreeSet<_>>();
        let distinct_curves = group
            .iter()
            .skip(1)
            .any(|e| !same_virtual_circle(&group[0].virtual_edge, &e.virtual_edge));
        let eligible = group.len() > 1
            && distinct_curves
            && members.iter().all(|&i| {
                terminals[i].as_ref().is_some_and(|edges| {
                    edges.len() == 2 && edges.iter().all(|id| owners[id].len() == 2)
                })
            })
            && group.iter().all(|e| {
                e.carrier_pair_unique
                    && e.strength == EdgeTreatmentEvidenceStrength::AnalyticTangency
                    && e.carriers.iter().all(|c| !group_ids.contains(&c.face_id))
            });
        if eligible {
            closed.push(group);
        } else {
            other.extend(group);
        }
    }
    (closed, other)
}

#[derive(Clone, Copy)]
struct Junction {
    terminal: u64,
    point: [f64; 3],
}

struct CarrierRun {
    face_id: u64,
    face_forward: bool,
    shared: BTreeSet<u64>,
    previous: OrientedEdgeUse,
    next: OrientedEdgeUse,
    ends: [Junction; 2],
    circular_increasing: Option<bool>,
}

fn terminal_at_vertex(terminals: &[&OrientedEdgeUse], vertex: u64) -> Result<u64> {
    let mut matching = terminals
        .iter()
        .filter(|e| e.start_vertex == vertex || e.end_vertex == vertex);
    let first = matching
        .next()
        .context("carrier endpoint has no original terminal edge")?;
    ensure!(
        matching.next().is_none(),
        "carrier endpoint has ambiguous terminal edges"
    );
    Ok(first.edge_id)
}

fn segment_junction(
    carrier: &ReducedFace,
    start: [f64; 3],
    end: [f64; 3],
    source: [f64; 3],
) -> Result<[f64; 3]> {
    let delta = sub(end, start);
    let length = norm(delta);
    let axis = normalize(delta, GEOM_TOL_MM).context("degenerate sharp segment")?;
    let SurfaceSupport::Plane(plane) = &carrier.surface else {
        anyhow::bail!("linear coupled treatment needs a planar carrier");
    };
    ensure!(
        [start, end]
            .into_iter()
            .all(|p| dot(sub(p, plane.origin_mm), plane.normal).abs() <= GEOM_TOL_MM),
        "sharp segment is not on original carrier #{}",
        carrier.id
    );
    let projected = dot(sub(source, start), axis);
    if projected.abs() <= GEOM_TOL_MM {
        Ok(start)
    } else if (projected - length).abs() <= GEOM_TOL_MM {
        Ok(end)
    } else {
        anyhow::bail!("source terminal does not project to either sharp segment endpoint")
    }
}

fn carrier_run(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    carrier: &crate::edge_treatments::EdgeTreatmentCarrier,
    terminals: &[&OrientedEdgeUse],
) -> Result<CarrierRun> {
    let face = face(solid, carrier.face_id)?;
    let shared = carrier
        .shared_edge_ids
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    ensure!(
        !shared.is_empty() && shared.len() == carrier.shared_edge_ids.len(),
        "coupled carrier has invalid shared edge use"
    );
    let (previous, next) = carrier_run_neighbors(face, &shared)?;
    let sources = [
        (previous.end_vertex, previous.end_mm),
        (next.start_vertex, next.start_mm),
    ];
    let terminal_ids = sources.map(|(v, _)| terminal_at_vertex(terminals, v));
    let [a, b] = terminal_ids;
    let ids = [a?, b?];
    ensure!(ids[0] != ids[1], "carrier has identical terminal ends");
    let (points, circular_increasing) = match evidence.virtual_edge {
        EdgeTreatmentVirtualEdge::Segment { start_mm, end_mm } => {
            let axis =
                normalize(sub(end_mm, start_mm), GEOM_TOL_MM).context("degenerate sharp line")?;
            for edge in face.loop_edges().filter(|e| shared.contains(&e.edge_id)) {
                let CurveSupport::Line(line) = edge.support else {
                    anyhow::bail!("linear treatment has non-linear interface");
                };
                ensure!(
                    normalize(line.direction, DIR_TOL)
                        .is_some_and(|d| dot(d, axis).abs() >= 1.0 - DIR_TOL),
                    "linear treatment interface does not follow sharp line"
                );
            }
            (
                [
                    segment_junction(face, start_mm, end_mm, sources[0].1)?,
                    segment_junction(face, start_mm, end_mm, sources[1].1)?,
                ],
                None,
            )
        }
        EdgeTreatmentVirtualEdge::Circle { .. } => {
            let frame = CircleFrame::from_edge(&evidence.virtual_edge)?;
            ensure!(
                circle_carrier_supports_inverse(face, frame),
                "original circular carrier lacks sharp support"
            );
            let run = circle_carrier_run(face, shared.clone(), frame)?;
            (
                [
                    frame.point(frame.angle(sources[0].1)?),
                    frame.point(frame.angle(sources[1].1)?),
                ],
                Some(run.sweep > 0.0),
            )
        }
    };
    ensure!(
        distance(points[0], points[1]) > GEOM_TOL_MM,
        "collapsed coupled inverse edge"
    );
    Ok(CarrierRun {
        face_id: face.id,
        face_forward: face.surface_forward,
        shared,
        previous,
        next,
        ends: std::array::from_fn(|i| Junction {
            terminal: ids[i],
            point: points[i],
        }),
        circular_increasing,
    })
}

struct DirectedChain {
    points: Vec<[f64; 3]>,
    edges: Vec<KernelNewEdge>,
    terminals: [u64; 2],
}

struct Treatment {
    id: u64,
    runs: [CarrierRun; 2],
    geometry: EdgeTreatmentVirtualEdge,
    all_edges: BTreeSet<u64>,
}

impl Treatment {
    fn new(solid: &ReducedSolid, evidence: &EdgeTreatmentEvidence) -> Result<Self> {
        let source = face(solid, evidence.treatment_face_id)?;
        let all_edges = unique_face_edges(source)?;
        let terminal_ids = terminal_edges(solid, evidence)?;
        let terminals = source
            .loop_edges()
            .filter(|e| terminal_ids.contains(&e.edge_id))
            .collect::<Vec<_>>();
        ensure!(
            terminals.len() == 2,
            "closed coupled treatment requires two terminals"
        );
        let [a, b] = &evidence.carriers;
        let runs = [
            carrier_run(solid, evidence, a, &terminals)?,
            carrier_run(solid, evidence, b, &terminals)?,
        ];
        ensure!(
            runs[0].face_id != runs[1].face_id && runs[0].shared.is_disjoint(&runs[1].shared),
            "coupled treatment does not have distinct carrier interfaces"
        );
        for end in runs[0].ends {
            ensure!(
                runs[1]
                    .ends
                    .iter()
                    .any(|other| other.terminal == end.terminal
                        && distance(other.point, end.point) <= GEOM_TOL_MM),
                "coupled treatment #{} carrier endpoints disagree at terminal #{}",
                evidence.treatment_face_id,
                end.terminal
            );
        }
        Ok(Self {
            id: evidence.treatment_face_id,
            runs,
            geometry: evidence.virtual_edge.clone(),
            all_edges,
        })
    }

    /// Produce one directed chain; boundary wrappers use it forward or
    /// reversed according to original terminal ownership, never proximity.
    fn chain(&self) -> Result<DirectedChain> {
        let [a, b] = self.runs[0].ends;
        match &self.geometry {
            EdgeTreatmentVirtualEdge::Segment { .. } => Ok(DirectedChain {
                points: vec![a.point, b.point],
                edges: vec![KernelNewEdge {
                    vertices: (KernelPatchVertexRef::New(0), KernelPatchVertexRef::New(1)),
                    curve: KernelPatchCurve::Line {
                        start_mm: a.point,
                        end_mm: b.point,
                    },
                }],
                terminals: [a.terminal, b.terminal],
            }),
            EdgeTreatmentVirtualEdge::Circle { .. } => {
                let increasing = self.runs[0]
                    .circular_increasing
                    .context("missing circular direction")?;
                let (points, edges) = open_inverse_circle_arc(
                    CircleFrame::from_edge(&self.geometry)?,
                    a.point,
                    b.point,
                    increasing,
                )?;
                Ok(DirectedChain {
                    points,
                    edges,
                    terminals: if increasing {
                        [a.terminal, b.terminal]
                    } else {
                        [b.terminal, a.terminal]
                    },
                })
            }
        }
    }
}

struct TerminalJoint {
    vertex: usize,
    owners: BTreeSet<u64>,
}

struct Perimeter {
    vertices: Vec<[f64; 3]>,
    joints: BTreeMap<u64, TerminalJoint>,
    new_edges: Vec<KernelNewEdge>,
    splices: Vec<KernelFaceSplice>,
    edits: Vec<EndpointEdit>,
    interfaces: BTreeMap<u64, (u64, u64)>,
}

impl Perimeter {
    fn new(solid: &ReducedSolid, treatments: &[Treatment]) -> Result<Self> {
        let mut result = Self {
            vertices: Vec::new(),
            joints: BTreeMap::new(),
            new_edges: Vec::new(),
            splices: Vec::new(),
            edits: Vec::new(),
            interfaces: BTreeMap::new(),
        };
        for treatment in treatments {
            for run in &treatment.runs {
                for &edge in &run.shared {
                    ensure!(
                        result
                            .interfaces
                            .insert(edge, (treatment.id, run.face_id))
                            .is_none(),
                        "coupled source edge has multiple treatment owners"
                    );
                    ensure!(
                        edge_face_users(solid, edge) == BTreeSet::from([treatment.id, run.face_id]),
                        "coupled carrier interface is not manifold"
                    );
                }
            }
            for end in treatment.runs[0].ends {
                let joint = result.joints.entry(end.terminal).or_insert_with(|| {
                    let vertex = result.vertices.len();
                    result.vertices.push(end.point);
                    TerminalJoint {
                        vertex,
                        owners: BTreeSet::new(),
                    }
                });
                ensure!(
                    distance(result.vertices[joint.vertex], end.point) <= GEOM_TOL_MM,
                    "sharp treatment junction disagrees at source terminal #{}",
                    end.terminal
                );
                ensure!(
                    joint.owners.insert(treatment.id),
                    "duplicate treatment junction owner"
                );
            }
        }
        for (&terminal, joint) in &result.joints {
            ensure!(
                joint.owners.len() == 2 && edge_face_users(solid, terminal) == joint.owners,
                "source terminal #{terminal} is not a closed two-treatment junction"
            );
        }
        Ok(result)
    }

    fn add_chain(&mut self, treatment: &Treatment) -> Result<(usize, usize, [u64; 2])> {
        let DirectedChain {
            points,
            edges,
            terminals,
        } = treatment.chain()?;
        let mut vertices = Vec::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            let vertex = if i == 0 {
                self.joints[&terminals[0]].vertex
            } else if i + 1 == points.len() {
                self.joints[&terminals[1]].vertex
            } else {
                let vertex = self.vertices.len();
                self.vertices.push(*p);
                vertex
            };
            ensure!(
                distance(*p, self.vertices[vertex]) <= GEOM_TOL_MM,
                "generated chain endpoint disagrees with terminal"
            );
            vertices.push(vertex);
        }
        let from = self.new_edges.len();
        let remap = |reference| match reference {
            KernelPatchVertexRef::New(i) => Ok(KernelPatchVertexRef::New(vertices[i])),
            KernelPatchVertexRef::Existing(_) => {
                anyhow::bail!("inverse chain has unexpected existing vertex")
            }
        };
        for mut edge in edges {
            edge.vertices = (remap(edge.vertices.0)?, remap(edge.vertices.1)?);
            self.new_edges.push(edge);
        }
        Ok((from, self.new_edges.len(), terminals))
    }

    fn edit_neighbor(
        &mut self,
        treatment: &Treatment,
        run: &CarrierRun,
        neighbor: &OrientedEdgeUse,
        source_vertex: u64,
        end: Junction,
    ) -> Result<()> {
        let joint = &self.joints[&end.terminal];
        if let Some(&(owner, carrier)) = self.interfaces.get(&neighbor.edge_id) {
            ensure!(
                owner != treatment.id && joint.owners.contains(&owner) && carrier == run.face_id,
                "adjacent removed carrier edge is not owned by the same terminal junction"
            );
            return Ok(());
        }
        let CurveSupport::Line(line) = neighbor.support else {
            anyhow::bail!(
                "coupled retained edge #{} is not straight",
                neighbor.edge_id
            );
        };
        let axis = normalize(line.direction, DIR_TOL).context("invalid retained line direction")?;
        let point = self.vertices[joint.vertex];
        ensure!(
            point_to_unit_line_distance(point, line.origin_mm, axis) <= GEOM_TOL_MM,
            "authored retained line #{} does not reach coupled junction",
            neighbor.edge_id
        );
        self.edits.push(circle_line_endpoint_edit(
            run.face_id,
            neighbor.clone(),
            source_vertex,
            joint.vertex,
            point,
        )?);
        Ok(())
    }

    fn add_treatment(
        &mut self,
        treatment: &Treatment,
        section: usize,
        provenance: &KernelTopologyProvenance,
    ) -> Result<()> {
        let (from, to, terminals) = self.add_chain(treatment)?;
        let mut orientations = [false; 2];
        for (i, run) in treatment.runs.iter().enumerate() {
            let order = run.ends.map(|end| end.terminal);
            let forward = order == terminals;
            ensure!(
                forward || order == [terminals[1], terminals[0]],
                "carrier has unrelated chain endpoints"
            );
            ensure!(
                run.circular_increasing
                    .is_none_or(|increasing| increasing == forward),
                "circular carrier traversal contradicts chain direction"
            );
            orientations[i] = forward == run.face_forward;
            let replacement = (from..to)
                .map(|index| KernelPatchEdgeUse {
                    edge: KernelPatchEdgeRef::New(if forward {
                        index
                    } else {
                        to - 1 - (index - from)
                    }),
                    orientation: forward,
                })
                .collect();
            self.splices.push(KernelFaceSplice {
                face: topology_face(provenance, section, run.face_id)?,
                remove_edges: topology_edges(provenance, section, run.shared.iter().copied())?,
                replacement,
            });
            self.edit_neighbor(
                treatment,
                run,
                &run.previous,
                run.previous.end_vertex,
                run.ends[0],
            )?;
            self.edit_neighbor(
                treatment,
                run,
                &run.next,
                run.next.start_vertex,
                run.ends[1],
            )?;
        }
        ensure!(
            orientations[0] != orientations[1],
            "coupled sharp edge lacks opposing shell uses"
        );
        Ok(())
    }
}

pub fn plan_closed_perimeter(
    solid: &ReducedSolid,
    group: &[&EdgeTreatmentEvidence],
    section: usize,
    provenance: &KernelTopologyProvenance,
) -> Result<KernelTopologyPatch> {
    let treatments = group
        .iter()
        .map(|e| {
            Treatment::new(solid, e)
                .with_context(|| format!("coupled treatment #{}", e.treatment_face_id))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut perimeter = Perimeter::new(solid, &treatments)?;
    for treatment in &treatments {
        perimeter.add_treatment(treatment, section, provenance)?;
    }
    // Internal source terminal edges belong to two removed faces.
    let removed_edges = treatments
        .iter()
        .flat_map(|t| t.all_edges.iter().copied())
        .collect::<BTreeSet<_>>();
    Ok(KernelTopologyPatch {
        remove_faces: treatments
            .iter()
            .map(|t| topology_face(provenance, section, t.id))
            .collect::<Result<_>>()?,
        new_vertices: perimeter.vertices,
        new_edges: perimeter.new_edges,
        rewrite_edges: normalize_retained_edge_rewrites(&perimeter.edits, section, provenance)?,
        face_splices: perimeter.splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(provenance, section, removed_edges)?,
    })
}

#[cfg(all(test, feature = "cad-kernel-monstertruck"))]
mod tests {
    use super::*;
    use crate::cad_kernel::monstertruck::MonstertruckKernel;
    use crate::kernel_source_load::SourceStepLoadInspector;
    use crate::kernel_topology::{SourceTopologyEditor, plan_source_edge_treatment_patches};
    use std::io::{Cursor, Read};

    #[test]
    fn real_te_closed_perimeter_is_one_native_topology_transaction() -> Result<()> {
        // Existing bundled, SHA-pinned manufacturer source. No synthetic
        // source geometry or manually adjusted candidate parameters.
        let archive = include_bytes!("../../validation/fixtures/bundled/te-T2111092004-000.zip");
        let mut zip = zip::ZipArchive::new(Cursor::new(archive))?;
        let mut bytes = Vec::new();
        zip.by_index(0)?.read_to_end(&mut bytes)?;
        use sha2::{Digest, Sha256};
        assert_eq!(
            hex::encode(Sha256::digest(&bytes)),
            "b485ad005657eec408ee345f77d68848d54c42940892a9e40aae30d2d18d95cf"
        );
        let kernel = MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(&bytes)?.exchange,
        )?;
        let inspected = kernel.inspect_source_document(&document)?;
        let loaded = inspected
            .solids
            .iter()
            .find(|s| s.source.id == 9354)
            .context("TE source solid absent")?;
        let provenance = loaded
            .topology_provenance
            .as_ref()
            .context("TE source has no stable topology")?;
        let entities = &document.as_ref().data[loaded.source.data_section].entities;
        let index = crate::step_graph::build_index(entities);
        let all = crate::edge_treatments::detect_edge_treatment_evidence_for_solids_with_index(
            entities,
            &index,
            &std::collections::HashSet::from([9354]),
        );
        let reduced =
            ReducedSolid::from_step(9354, entities, &index).context("missing reduced source")?;
        let refs = all.iter().collect::<Vec<_>>();
        let (groups, other) = closed_components(&reduced, &refs);
        assert_eq!(groups.len(), 3);
        assert!(groups.iter().all(|group| group.len() == 8));
        assert_eq!(other.len(), all.len() - 24);
        assert!(
            groups
                .iter()
                .any(|group| group.iter().any(|e| e.treatment_face_id == 7893))
        );
        let mut selected = groups.into_iter().flatten().cloned().collect::<Vec<_>>();
        // Prove all three perimeters in one transaction, independently of
        // candidate ordering. Larger-radius candidate carriers must survive.
        selected.reverse();
        let patch = plan_source_edge_treatment_patches(
            document.as_ref(),
            loaded.source,
            &selected,
            provenance,
        )?;
        assert_eq!(patch.remove_faces.len(), 24);
        let removed = patch.remove_faces.iter().copied().collect::<BTreeSet<_>>();
        let stable_faces = provenance
            .faces
            .iter()
            .map(|f| f.topology_id)
            .filter(|id| !removed.contains(id))
            .collect::<BTreeSet<_>>();
        let edited = kernel.edit_source_topology(&document, loaded.source, &[patch])?;
        let summary = kernel.inspect_edited_source(&edited)?;
        assert_eq!(summary.boundary_shells, 1);
        assert_eq!(summary.faces, 285);
        assert_eq!(
            summary
                .retained_source_topology
                .faces
                .iter()
                .map(|f| f.topology_id)
                .collect::<BTreeSet<_>>(),
            stable_faces
        );
        Ok(())
    }
}
