//! Invert a straight edge treatment, retaining bounded neighboring edges.
//! At a terminal shoulder a tangent contact is also an existing corner of
//! a third face. Preserve that corner and connect it to the sharp edge in
//! the intersection of the carrier and terminal planes.
use super::{
    CurveSupport, DIR_TOL, EdgeTreatmentEvidence, EndpointEdit, GEOM_TOL_MM, KernelFaceSplice,
    KernelNewEdge, KernelPatchCurve, KernelPatchEdgeRef, KernelPatchEdgeUse, KernelPatchVertexRef,
    KernelTopologyPatch, OrientedEdgeUse, ReducedFace, ReducedSolid, SurfaceSupport,
    TreatmentTopology, carrier_run_neighbors, curved_segment, distance, dot, edge_face_users,
    edges_form_one_cyclic_run, face, normalize, normalize_retained_edge_rewrites,
    point_to_unit_line_distance, sub, terminal_run_edges, topology_edges, topology_face,
    topology_vertex, unique_face_edges,
};
use anyhow::{Context, Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

struct SharpSegment<'a> {
    solid: &'a ReducedSolid,
    points: [[f64; 3]; 2],
    terminals: BTreeMap<u64, BTreeSet<u64>>,
}

struct Contact {
    edit: EndpointEdit,
    terminal_face: Option<u64>,
}

struct CarrierPlan {
    splice: KernelFaceSplice,
    contacts: [Contact; 2],
}

struct Connector {
    terminal_face: u64,
    old_vertex: u64,
    new_vertex: usize,
    edge_index: usize,
}

impl SharpSegment<'_> {
    /// A mismatch with a retained line is admissible only when the actual
    /// treatment terminal meets that line's endpoint, both endpoints of the
    /// new connector lie on both incident planes, and the retained line is
    /// shared with a different face. No extrapolation of that line is allowed.
    fn shoulder_terminal(
        &self,
        carrier: &ReducedFace,
        edge: &OrientedEdgeUse,
        old_vertex: u64,
        old_point: [f64; 3],
        target: [f64; 3],
    ) -> Option<u64> {
        let SurfaceSupport::Plane(carrier_plane) = carrier.surface else {
            return None;
        };
        let carrier_normal = normalize(carrier_plane.normal, DIR_TOL)?;
        let axis = normalize(sub(self.points[1], self.points[0]), DIR_TOL)?;
        let retained_users = edge_face_users(self.solid, edge.edge_id);
        if retained_users.len() != 2 || !retained_users.contains(&carrier.id) {
            return None;
        }
        if distance(old_point, target) <= GEOM_TOL_MM
            || [old_point, target].into_iter().any(|point| {
                dot(carrier_normal, sub(point, carrier_plane.origin_mm)).abs() > GEOM_TOL_MM
            })
        {
            return None;
        }
        let mut candidates = self.terminals.iter().filter_map(|(&id, edges)| {
            if retained_users.contains(&id) {
                return None;
            }
            let terminal = face(self.solid, id).ok()?;
            let SurfaceSupport::Plane(plane) = terminal.surface else {
                return None;
            };
            let normal = normalize(plane.normal, DIR_TOL)?;
            if dot(normal, axis).abs() <= DIR_TOL
                || [old_point, target]
                    .into_iter()
                    .any(|point| dot(normal, sub(point, plane.origin_mm)).abs() > GEOM_TOL_MM)
            {
                return None;
            }
            let (before, after) = carrier_run_neighbors(terminal, edges).ok()?;
            [before.end_vertex, after.start_vertex]
                .contains(&old_vertex)
                .then_some(id)
        });
        let selected = candidates.next()?;
        candidates.next().is_none().then_some(selected)
    }

    fn contact(
        &self,
        carrier: &ReducedFace,
        edge: &OrientedEdgeUse,
        previous: bool,
        target: usize,
    ) -> Option<Contact> {
        let CurveSupport::Line(line) = edge.support else {
            return None;
        };
        let direction = normalize(line.direction, DIR_TOL)?;
        let (old_vertex, old_point) = if previous {
            (edge.end_vertex, edge.end_mm)
        } else {
            (edge.start_vertex, edge.start_mm)
        };
        let new_point = self.points[target];
        let terminal_face =
            if point_to_unit_line_distance(new_point, line.origin_mm, direction) <= GEOM_TOL_MM {
                None
            } else {
                Some(self.shoulder_terminal(carrier, edge, old_vertex, old_point, new_point)?)
            };
        Some(Contact {
            edit: EndpointEdit {
                carrier_face_id: carrier.id,
                old_vertex,
                edge: edge.clone(),
                new_vertex: target,
                new_point,
                replacement_curve: None,
            },
            terminal_face,
        })
    }

    fn carrier(
        &self,
        carrier: &ReducedFace,
        shared: &BTreeSet<u64>,
        section: usize,
        provenance: &impl TreatmentTopology,
    ) -> Result<CarrierPlan> {
        let (previous, next) = carrier_run_neighbors(carrier, shared)?;
        let mut candidates = [[0, 1], [1, 0]].into_iter().filter_map(|[a, b]| {
            Some([
                self.contact(carrier, &previous, true, a)?,
                self.contact(carrier, &next, false, b)?,
            ])
        });
        let contacts = candidates.next().with_context(|| {
            format!(
                "proven inverse segment does not meet retained neighbors on carrier face #{}",
                carrier.id
            )
        })?;
        ensure!(
            candidates.next().is_none(),
            "ambiguous inverse segment direction on carrier face #{}",
            carrier.id
        );
        Ok(CarrierPlan {
            splice: KernelFaceSplice {
                face: topology_face(provenance, section, carrier.id)?,
                remove_edges: topology_edges(provenance, section, shared.iter().copied())?,
                replacement: vec![KernelPatchEdgeUse {
                    edge: KernelPatchEdgeRef::New(0),
                    orientation: contacts[0].edit.new_vertex == 0,
                }],
            },
            contacts,
        })
    }

    fn terminal_splices(
        &self,
        connectors: &[Connector],
        section: usize,
        provenance: &impl TreatmentTopology,
    ) -> Result<Vec<KernelFaceSplice>> {
        self.terminals
            .iter()
            .map(|(&id, edges)| {
                let terminal = face(self.solid, id)?;
                ensure!(
                    edges_form_one_cyclic_run(terminal, edges),
                    "terminal edges {edges:?} are not one contiguous run on face #{id}"
                );
                let local = connectors
                    .iter()
                    .filter(|connector| connector.terminal_face == id)
                    .collect::<Vec<_>>();
                let mut replacement = Vec::new();
                if !local.is_empty() {
                    let (before, after) = carrier_run_neighbors(terminal, edges)?;
                    let mut shared_target = None;
                    let mut consumed = 0;
                    for (vertex, orientation) in
                        [(before.end_vertex, true), (after.start_vertex, false)]
                    {
                        let matches = local
                            .iter()
                            .filter(|connector| connector.old_vertex == vertex)
                            .collect::<Vec<_>>();
                        let Some(connector) = matches.first() else {
                            continue;
                        };
                        ensure!(
                            matches.len() == 1,
                            "multiple connectors at terminal vertex #{vertex}"
                        );
                        ensure!(
                            shared_target.is_none_or(|target| target == connector.new_vertex),
                            "terminal connectors disagree on sharp vertex"
                        );
                        shared_target = Some(connector.new_vertex);
                        consumed += 1;
                        replacement.push(KernelPatchEdgeUse {
                            edge: KernelPatchEdgeRef::New(connector.edge_index),
                            orientation,
                        });
                    }
                    ensure!(
                        consumed == local.len(),
                        "connector does not belong to either end of terminal run on face #{id}"
                    );
                }
                Ok(KernelFaceSplice {
                    face: topology_face(provenance, section, id)?,
                    remove_edges: topology_edges(provenance, section, edges.iter().copied())?,
                    replacement,
                })
            })
            .collect()
    }
}

struct SegmentEdits {
    edges: Vec<KernelNewEdge>,
    rewrites: Vec<EndpointEdit>,
    splices: Vec<KernelFaceSplice>,
    connectors: Vec<Connector>,
}

impl SegmentEdits {
    fn new(points: [[f64; 3]; 2]) -> Self {
        Self {
            edges: vec![KernelNewEdge {
                vertices: (KernelPatchVertexRef::New(0), KernelPatchVertexRef::New(1)),
                curve: KernelPatchCurve::Line {
                    start_mm: points[0],
                    end_mm: points[1],
                },
            }],
            rewrites: Vec::new(),
            splices: Vec::new(),
            connectors: Vec::new(),
        }
    }

    fn add_carrier(
        &mut self,
        plan: CarrierPlan,
        section: usize,
        provenance: &impl TreatmentTopology,
    ) -> Result<()> {
        let mut splice = plan.splice;
        for (i, contact) in plan.contacts.into_iter().enumerate() {
            let edit = contact.edit;
            let Some(terminal_face) = contact.terminal_face else {
                self.rewrites.push(edit);
                continue;
            };
            let old_point = if i == 0 {
                edit.edge.end_mm
            } else {
                edit.edge.start_mm
            };
            let edge_index = self.edges.len();
            self.edges.push(KernelNewEdge {
                vertices: (
                    KernelPatchVertexRef::Existing(topology_vertex(
                        provenance,
                        section,
                        edit.old_vertex,
                    )?),
                    KernelPatchVertexRef::New(edit.new_vertex),
                ),
                curve: KernelPatchCurve::Line {
                    start_mm: old_point,
                    end_mm: edit.new_point,
                },
            });
            self.connectors.push(Connector {
                terminal_face,
                old_vertex: edit.old_vertex,
                new_vertex: edit.new_vertex,
                edge_index,
            });
            let use_ = KernelPatchEdgeUse {
                edge: KernelPatchEdgeRef::New(edge_index),
                orientation: i == 0,
            };
            if i == 0 {
                splice.replacement.insert(0, use_);
            } else {
                splice.replacement.push(use_);
            }
        }
        self.splices.push(splice);
        Ok(())
    }
}

pub(super) fn plan_segment_patch(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    section: usize,
    provenance: &impl TreatmentTopology,
    virtual_start: [f64; 3],
    virtual_end: [f64; 3],
) -> Result<KernelTopologyPatch> {
    ensure!(
        distance(virtual_start, virtual_end) > GEOM_TOL_MM,
        "inverse segment for treatment face #{} is degenerate",
        evidence.treatment_face_id
    );
    if let Some(patch) = curved_segment::plan_curved_segment_patch(
        solid,
        evidence,
        section,
        provenance,
        virtual_start,
        virtual_end,
    )? {
        return Ok(patch);
    }
    let treatment = face(solid, evidence.treatment_face_id)?;
    let treatment_edges = unique_face_edges(treatment)?;
    let mut interfaces = BTreeSet::new();
    for carrier in &evidence.carriers {
        for &id in &carrier.shared_edge_ids {
            ensure!(
                interfaces.insert(id),
                "treatment face #{} assigns source edge #{id} to both carriers",
                evidence.treatment_face_id
            );
        }
    }
    ensure!(
        interfaces.is_subset(&treatment_edges),
        "carrier interface contains an edge absent from treatment face #{}",
        evidence.treatment_face_id
    );
    let terminal_edges = treatment_edges
        .difference(&interfaces)
        .copied()
        .collect::<BTreeSet<_>>();
    ensure!(
        !terminal_edges.is_empty(),
        "segment treatment face #{} has no terminal boundary",
        evidence.treatment_face_id
    );
    let bounds = SharpSegment {
        solid,
        points: [virtual_start, virtual_end],
        terminals: terminal_run_edges(solid, evidence, &terminal_edges)?,
    };
    let mut edits = SegmentEdits::new(bounds.points);
    for carrier in &evidence.carriers {
        let shared = carrier
            .shared_edge_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        ensure!(
            !shared.is_empty() && shared.len() == carrier.shared_edge_ids.len(),
            "carrier face #{} has empty or repeated treatment edges",
            carrier.face_id
        );
        let plan = bounds.carrier(face(solid, carrier.face_id)?, &shared, section, provenance)?;
        edits.add_carrier(plan, section, provenance)?;
    }
    edits
        .splices
        .extend(bounds.terminal_splices(&edits.connectors, section, provenance)?);
    Ok(KernelTopologyPatch {
        remove_faces: vec![topology_face(
            provenance,
            section,
            evidence.treatment_face_id,
        )?],
        new_vertices: bounds.points.to_vec(),
        new_edges: edits.edges,
        rewrite_edges: normalize_retained_edge_rewrites(&edits.rewrites, section, provenance)?,
        face_splices: edits.splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(provenance, section, treatment_edges)?,
    })
}
