//! Separate source-document repairs for regular two-cylinder intersections.
//!
//! The existing face carriers, topology, coedge orientation and vertices
//! are retained. The new finite STEP leader is a bounded approximation of
//! the analytical branch, associated with both original cylinder surfaces.
//! The original spline records remain available for comparison and the
//! returned document is not a certificate of original-source equivalence.

use crate::{
    brep,
    source_edge_healing::{
        CylinderHealingEvidence, CylinderHealingRequest, propose_cylinder_intersection_healing,
    },
    source_edge_locus::StepEdge,
    step_entities::{enumeration_bool, simple_parameters_mut},
    step_graph::{build_index, entity_id, entity_ref_value, simple_record, visit_entity_refs},
    step_io::write_exchange,
};
use anyhow::{Context, Result, bail};
use monstertruck_geometry::prelude::{InnerSpace, Matrix4, Point3, Transformed};
use monstertruck_io::step::{load::ParsedStepSource, save::StepDisplay};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use step_p21::ast::{DataSection, EntityInstance, Name, Parameter};

/// Explicit request to heal a particular source edge in its own solid.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct FullSourceCylinderRepair {
    pub data_section: usize,
    pub solid_id: u64,
    pub edge_id: u64,
    pub endpoint_tolerance_mm: f64,
    pub target_gap_mm: f64,
    pub leader_tolerance_mm: f64,
    pub max_leader_segments: u32,
}

/// Source-distinct geometry and topology provenance for one replacement.
#[derive(Clone, Debug, Serialize)]
pub struct FullSourceCylinderRepairEvidence {
    pub data_section: usize,
    pub original: CylinderHealingEvidence,
    pub original_surface_ids: [u64; 2],
    pub intersection_curve_id: u64,
    pub leader_curve_id: u64,
    pub new_geometry_entity_ids: Vec<u64>,
    pub leader_segments: u32,
    pub leader_degree: usize,
    pub leader_control_count: usize,
    pub leader_error_upper_mm: f64,
    pub leader_arithmetic_margin_mm: f64,
    pub serialized_control_error_upper_mm: f64,
    pub source_to_leader_separation_mm: [f64; 2],
}

/// Separately serialized document with its explicit repair manifest.
#[derive(Clone, Debug)]
pub struct FullSourceCylinderHealing {
    pub step: String,
    pub repairs: Vec<FullSourceCylinderRepairEvidence>,
}

fn verify_scope(
    source: &DataSection,
    index: &HashMap<u64, usize>,
    evidence: &CylinderHealingEvidence,
) -> Result<()> {
    if index.len() != source.entities.len() {
        bail!("original source section contains duplicate entity IDs");
    }
    // Reject a modification of a topological edge shared with an unrelated
    // face or solid. Reusing the *geometry* by other entities is harmless:
    // the old geometry records are never edited or removed.
    let mut parents = HashMap::<u64, Vec<u64>>::new();
    for entity in &source.entities {
        visit_entity_refs(entity, &mut |child| {
            parents.entry(child).or_default().push(entity_id(entity));
        });
    }
    let mut visited = HashSet::new();
    let mut pending = vec![evidence.original_edge_id];
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        if let Some(record) = index
            .get(&id)
            .and_then(|&i| source.entities.get(i))
            .and_then(simple_record)
        {
            let outside = match record.name.as_str() {
                "ADVANCED_FACE" | "FACE_SURFACE" => {
                    id != evidence.primary_face_id && id != evidence.secondary_face_id
                }
                "MANIFOLD_SOLID_BREP" | "BREP_WITH_VOIDS" => id != evidence.original_solid_id,
                _ => false,
            };
            if outside {
                bail!("source edge is also referenced by unrelated topology #{id}");
            }
        }
        if let Some(referencing) = parents.get(&id) {
            pending.extend(referencing);
        }
    }
    Ok(())
}

struct PreparedRepair {
    evidence: CylinderHealingEvidence,
    surfaces: [u64; 2],
    leader: monstertruck_geometry::analytic_cylinder_intersection::IntersectionLeader,
}

fn prepare(source: &DataSection, request: FullSourceCylinderRepair) -> Result<PreparedRepair> {
    let candidate = propose_cylinder_intersection_healing(CylinderHealingRequest {
        solid_id: request.solid_id,
        source: StepEdge {
            id: request.edge_id,
            entities: &source.entities,
        },
        endpoint_tolerance_mm: request.endpoint_tolerance_mm,
        target_gap_mm: request.target_gap_mm,
    })?;
    let index = build_index(&source.entities);
    verify_scope(source, &index, &candidate.evidence)?;
    let surface = |face| {
        brep::face_surface(face, &source.entities, &index)
            .context("original cylinder surface disappeared")
    };
    let surfaces = [
        surface(candidate.evidence.primary_face_id)?,
        surface(candidate.evidence.secondary_face_id)?,
    ];
    let leader = candidate
        .arc
        .bounded_leader(
            request.leader_tolerance_mm * 0.95,
            request.max_leader_segments,
        )
        .context("analytical intersection cannot meet the requested finite STEP leader budget")?;
    if candidate.evidence.locus_separation_mm[0] - leader.maximum_error
        <= candidate.evidence.source_representation_tolerance_mm
    {
        bail!("finite leader would not retain the independently proven source separation");
    }
    Ok(PreparedRepair {
        evidence: candidate.evidence,
        surfaces,
        leader,
    })
}

struct NativeEntities {
    entities: Vec<EntityInstance>,
    intersection: u64,
    leader: u64,
    control_error_mm: f64,
}

fn serialized_control_error(
    prepared: &PreparedRepair,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    leader: u64,
) -> Result<f64> {
    let brep::CurveSupport::BSpline(spline) = brep::curve_support(leader, entities, index) else {
        bail!("serialized leader does not resolve to its polynomial spline");
    };
    let physical = &prepared.leader.curve;
    if spline.weights.is_some()
        || spline.degree != physical.degree()
        || spline.control_points_mm.len() != physical.control_points().len()
        || !spline
            .knots
            .iter()
            .zip(physical.knot_vector().iter())
            .all(|(a, b)| a.total_cmp(b).is_eq())
        || spline.knots.len() != physical.knot_vector().len()
    {
        bail!("serialized leader changed its parameterization");
    }
    let scale = prepared.evidence.source_coordinate_scale_mm;
    let mut error = 0.0_f64;
    for (&point, &expected) in spline
        .control_points_mm
        .iter()
        .zip(physical.control_points())
    {
        let actual = Point3::new(point[0] * scale, point[1] * scale, point[2] * scale);
        let coordinate_scale = [
            actual.x.abs(),
            actual.y.abs(),
            actual.z.abs(),
            expected.x.abs(),
            expected.y.abs(),
            expected.z.abs(),
            1.0_f64,
        ]
        .into_iter()
        .fold(1.0_f64, f64::max);
        let deviation =
            (64.0 * f64::EPSILON).mul_add(coordinate_scale, (actual - expected).magnitude());
        if !deviation.is_finite() {
            bail!("nonfinite serialized leader coordinates");
        }
        error = error.max(deviation);
    }
    Ok(error)
}

fn native_entities(section: &DataSection, prepared: &PreparedRepair) -> Result<NativeEntities> {
    let high = section.entities.iter().map(entity_id).max().unwrap_or(0);
    let intersection = high
        .checked_add(1)
        .context("STEP entity namespace exhausted")?;
    let leader_id = intersection
        .checked_add(1)
        .context("STEP entity namespace exhausted")?;
    let mut leader = prepared.leader.curve.clone();
    let inverse_scale = prepared.evidence.source_coordinate_scale_mm.recip();
    if !inverse_scale.is_finite() || inverse_scale <= 0.0 {
        bail!("invalid coordinate conversion for finite STEP leader");
    }
    leader.transform_by(Matrix4::from_scale(inverse_scale));
    let first = usize::try_from(leader_id)?;
    let [primary, secondary] = prepared.surfaces;
    let fragment = format!(
        "#{intersection}=INTERSECTION_CURVE('',#{leader_id},(#{primary},#{secondary}),.CURVE_3D.);\n{}",
        StepDisplay::new(&leader, first),
    );
    let document = format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('bounded analytical intersection'),'2;1');\n\
        FILE_NAME('geometry','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{fragment}\nENDSEC;\nEND-ISO-10303-21;\n"
    );
    let mut exchange = step_p21::parser::parse(&document)?;
    let local = exchange
        .data
        .pop()
        .context("generated STEP leader has no DATA section")?;
    let index = build_index(&local.entities);
    if index.len() != local.entities.len()
        || !index.contains_key(&intersection)
        || !index.contains_key(&leader_id)
        || index.keys().any(|&id| id <= high)
    {
        bail!("generated intersection has invalid local entity IDs");
    }
    for entity in &local.entities {
        let mut valid = true;
        visit_entity_refs(entity, &mut |id| {
            valid &= index.contains_key(&id) || prepared.surfaces.contains(&id);
        });
        if !valid {
            bail!("generated intersection contains dangling references");
        }
    }
    let control_error_mm = serialized_control_error(prepared, &local.entities, &index, leader_id)?;
    Ok(NativeEntities {
        entities: local.entities,
        intersection,
        leader: leader_id,
        control_error_mm,
    })
}

fn apply(
    section: &mut DataSection,
    request: FullSourceCylinderRepair,
    prepared: PreparedRepair,
) -> Result<FullSourceCylinderRepairEvidence> {
    let NativeEntities {
        entities: extra,
        intersection,
        leader,
        control_error_mm,
    } = native_entities(section, &prepared)?;
    let maximum_error = prepared.leader.maximum_error + control_error_mm;
    if !maximum_error.is_finite()
        || maximum_error > request.leader_tolerance_mm
        || prepared.evidence.locus_separation_mm[0] - maximum_error
            <= prepared.evidence.source_representation_tolerance_mm
    {
        bail!("serialized leader exceeds its physical error or separation budget");
    }
    let index = build_index(&section.entities);
    let original = &prepared.evidence;
    let edge = simple_parameters_mut(section, &index, request.edge_id, "EDGE_CURVE")?;
    if edge.len() != 5
        || entity_ref_value(&edge[1]) != Some(original.source_start_vertex_id)
        || entity_ref_value(&edge[2]) != Some(original.source_end_vertex_id)
        || entity_ref_value(&edge[3]) != Some(original.original_basis_curve_id)
        || enumeration_bool(&edge[4]) != Some(original.source_curve_same_sense)
    {
        bail!("source edge topology changed during repair");
    }
    // The normalized leader parameter follows the original spline basis
    // direction even when the analytical sweep is negative. Consequently
    // the original EDGE_CURVE same_sense and all coedges remain unchanged.
    edge[3] = Parameter::Ref(Name::Entity(intersection));
    let ids = extra.iter().map(entity_id).collect();
    section.entities.extend(extra);
    Ok(FullSourceCylinderRepairEvidence {
        data_section: request.data_section,
        original_surface_ids: prepared.surfaces,
        intersection_curve_id: intersection,
        leader_curve_id: leader,
        new_geometry_entity_ids: ids,
        leader_segments: prepared.leader.segments,
        leader_degree: prepared.leader.curve.degree(),
        leader_control_count: prepared.leader.curve.control_points().len(),
        leader_error_upper_mm: maximum_error,
        leader_arithmetic_margin_mm: prepared.leader.arithmetic_margin,
        serialized_control_error_upper_mm: control_error_mm,
        source_to_leader_separation_mm: [
            (original.locus_separation_mm[0] - maximum_error).max(0.0),
            original.locus_separation_mm[1] + maximum_error,
        ],
        original: prepared.evidence,
    })
}

/// Construct a new source document from explicitly requested edge repairs.
///
/// Every candidate is established from the unmodified source, and all
/// existing entity records except the nominated `EDGE_CURVE` geometry
/// references remain intact. No new body is automatically admitted.
///
/// # Errors
/// Rejects malformed source geometry/topology, duplicate repairs, ambiguous
/// shared topology, insufficient separation, and unachievable STEP error budgets.
pub fn heal_cylinder_step_document(
    bytes: &[u8],
    requests: &[FullSourceCylinderRepair],
) -> Result<FullSourceCylinderHealing> {
    if requests.is_empty() {
        bail!("no requested cylinder intersection repairs");
    }
    let original = ParsedStepSource::parse(bytes)?.exchange;
    crate::step_io::require_supported_sections(&original)?;
    let mut revised = original.clone();
    let mut seen = HashSet::new();
    let mut manifest = Vec::with_capacity(requests.len());
    for &request in requests {
        if !seen.insert((request.data_section, request.edge_id)) {
            bail!("duplicate cylinder intersection repair");
        }
        let source = original
            .data
            .get(request.data_section)
            .context("source DATA section absent")?;
        let prepared = prepare(source, request)?;
        let destination = revised
            .data
            .get_mut(request.data_section)
            .context("destination DATA section absent")?;
        manifest.push(apply(destination, request, prepared)?);
    }
    let step = write_exchange(&revised)?;
    let roundtrip = ParsedStepSource::parse(step.as_bytes())?.exchange;
    if roundtrip.data != revised.data {
        bail!("serialized repaired source changed DATA entity values");
    }
    Ok(FullSourceCylinderHealing {
        step,
        repairs: manifest,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::{BoundedCurve, ParametricCurve};
    use monstertruck_io::step::load::EdgeCurveHolder;
    use step_p21::tables::EntityTable;

    const SOURCE: &[u8] =
        include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");

    fn request(edge_id: u64) -> FullSourceCylinderRepair {
        FullSourceCylinderRepair {
            data_section: 0,
            solid_id: 171,
            edge_id,
            endpoint_tolerance_mm: 1e-8,
            target_gap_mm: 1e-9,
            leader_tolerance_mm: 2e-10,
            max_leader_segments: 8192,
        }
    }

    #[test]
    fn original_cylinder_faces_and_directed_topology_survive_repaired_step() -> Result<()> {
        let original = ParsedStepSource::parse(SOURCE)?.resolve();
        let revised = heal_cylinder_step_document(SOURCE, &[request(1049), request(1063)])?;
        let parsed = ParsedStepSource::parse(revised.step.as_bytes())?.resolve();
        let source = &original.exchange().data[0].entities;
        let target = &parsed.exchange().data[0].entities;
        let index = build_index(target);
        for entity in source {
            let id = entity_id(entity);
            let other = &target[index[&id]];
            if id != 1049 && id != 1063 {
                assert_eq!(entity, other, "source entity #{id} was altered");
            }
        }
        let table = parsed.table(0).context("typed DATA section missing")?;
        assert!(table.entity_report.is_empty());
        for repair in revised.repairs {
            let edge = repair.original.original_edge_id;
            let holder =
                <monstertruck_io::step::load::Table as EntityTable<EdgeCurveHolder>>::get_owned(
                    table, edge,
                )?;
            assert_eq!(holder.same_sense, repair.original.source_curve_same_sense);
            let request = request(edge);
            let candidate = propose_cylinder_intersection_healing(CylinderHealingRequest {
                solid_id: request.solid_id,
                source: StepEdge {
                    id: edge,
                    entities: source,
                },
                endpoint_tolerance_mm: request.endpoint_tolerance_mm,
                target_gap_mm: request.target_gap_mm,
            })?;
            let curve = brep::curve_support(repair.leader_curve_id, target, &index);
            let brep::CurveSupport::BSpline(spline) = curve else {
                bail!("finite leader did not survive serialization");
            };
            let leader = spline.rational_curve().context("invalid imported leader")?;
            let (lo, hi) = leader.range_tuple();
            for i in 0..=512 {
                let t = f64::from(i) / 512.0;
                let reference = candidate
                    .arc
                    .evaluate(candidate.arc.sweep() * t)
                    .context("bad analytical reference")?;
                assert!(
                    (leader.evaluate((hi - lo).mul_add(t, lo)) - reference).magnitude()
                        <= repair.leader_error_upper_mm
                );
            }
            assert_eq!(repair.leader_degree, spline.degree);
            assert_eq!(repair.leader_control_count, spline.control_points_mm.len());
            eprintln!(
                "edge #{edge}: degree {}, {} segments, {} controls, bounded error {} mm",
                repair.leader_degree,
                repair.leader_segments,
                repair.leader_control_count,
                repair.leader_error_upper_mm
            );
        }
        let proof = crate::boundary_proof::compare_analytic_boundaries(171, target, 171, target);
        assert!(
            proof.equivalent,
            "healed self proof: {:?}",
            proof.first.reason
        );
        let original_proof =
            crate::boundary_proof::compare_analytic_boundaries(171, source, 171, target);
        assert!(!original_proof.equivalent);
        assert!(!original_proof.first.supported);
        let initial = original.import_solid(0, 171)?;
        let healed = parsed.import_solid(0, 171)?;
        assert_eq!(initial.provenance, healed.provenance);
        assert_eq!(initial.shell_reports, healed.shell_reports);
        Ok(())
    }
}
