//! Opt-in reconstruction of STEP solids with separately healed exact edges.
//!
//! Original input bytes are never modified. New geometry is inserted under
//! fresh source-section-local entity IDs; source splines, old planes,
//! original face-bound topology, product metadata and other DATA sections
//! remain available for provenance and independent comparison.

use crate::{
    brep,
    native_hyperbola_intersection::native_healed_intersection_curve,
    source_edge_healing::{HyperbolaHealingRequest, propose_plane_cone_hyperbola_healing},
    source_edge_locus::{StepEdge, source_solid_length_scale_mm},
    step_entities::{closure_from, named_entity_parameters},
    step_graph::{
        build_index, entity_id, entity_ref_value, rewrite_entity_refs, visit_entity_refs,
    },
    step_io::write_exchange,
};
use anyhow::{Context, Result, anyhow, bail};
use monstertruck_geometry::prelude::{Matrix4, Transformed};
use monstertruck_io::step::{load::ParsedStepSource, save::StepDisplay};
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use step_p21::ast::{DataSection, EntityInstance, Name, Parameter};

/// A requested independently proven replacement of one source edge.
#[derive(Clone, Copy, Debug, Serialize)]
pub struct FullSourceHyperbolaRepair {
    pub data_section: usize,
    pub solid_id: u64,
    pub edge_id: u64,
    pub plane_face_id: u64,
    pub cone_face_id: u64,
    /// Physical maximum source-vertex to candidate 3D edge mismatch.
    pub endpoint_tolerance_mm: f64,
    /// Physical desired gap of whole-curve analytic departure bounds.
    pub target_gap_mm: f64,
}

/// A provenance record for one repair in one source DATA section.
#[derive(Clone, Debug, Serialize)]
pub struct FullSourceHyperbolaRepairEvidence {
    pub data_section: usize,
    pub original_solid_id: u64,
    pub original_edge_id: u64,
    pub original_basis_curve_id: u64,
    pub original_plane_face_id: u64,
    pub original_cone_face_id: u64,
    pub original_plane_surface_id: u64,
    pub original_cone_surface_id: u64,
    pub new_intersection_curve_id: u64,
    pub new_plane_surface_id: u64,
    pub new_geometry_entity_ids: Vec<u64>,
    pub source_length_scale_mm: f64,
    pub original_curve_same_sense: bool,
    pub healed_curve_same_sense: bool,
    pub source_to_healed_locus_separation_lower_mm: f64,
    pub healed_source_plane_residual_upper_mm: f64,
}

/// Newly written STEP document plus machine-readable, source-distinct proof.
#[derive(Clone, Debug)]
pub struct FullSourceHyperbolaHealing {
    pub step: String,
    pub repairs: Vec<FullSourceHyperbolaRepairEvidence>,
}

fn simple_parameters_mut<'a>(
    section: &'a mut DataSection,
    index: &HashMap<u64, usize>,
    id: u64,
    expected: &str,
) -> Result<&'a mut [Parameter]> {
    let position = *index
        .get(&id)
        .ok_or_else(|| anyhow!("missing source entity #{id}"))?;
    let entity = section
        .entities
        .get_mut(position)
        .ok_or_else(|| anyhow!("invalid source entity index for #{id}"))?;
    let EntityInstance::Simple { id: actual, record } = entity else {
        bail!("source entity #{id} is complex, expected simple {expected}");
    };
    if *actual != id || record.name != expected {
        bail!("source entity #{id} has unexpected type or stale ID");
    }
    let Parameter::List(parameters) = &mut record.parameter else {
        bail!("source entity #{id} lacks a parameter list");
    };
    Ok(parameters)
}

fn reference(parameter: &Parameter) -> Result<u64> {
    entity_ref_value(parameter)
        .ok_or_else(|| anyhow!("expected an actual local STEP entity reference"))
}

fn named_local_id(
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    id: u64,
    expected: &str,
) -> Result<()> {
    named_entity_parameters(id, expected, entities, index)
        .ok_or_else(|| anyhow!("expected native geometry #{id} of type {expected}"))?;
    Ok(())
}

fn analytic_geometry_with_local_refs(
    candidate: &crate::source_edge_healing::HyperbolaHealingCandidate,
    source_scale_mm: f64,
) -> Result<(Vec<EntityInstance>, u64, u64, u64)> {
    if !source_scale_mm.is_finite() || source_scale_mm <= 0.0 {
        bail!("invalid original source coordinate scale");
    }
    // The analytic candidate is evaluated in physical millimetres.
    // STEP geometry must instead use the specific source DATA section unit.
    let inverse_scale = 1.0 / source_scale_mm;
    if !inverse_scale.is_finite() || inverse_scale <= 0.0 {
        bail!("source unit conversion is not representable by a finite transform");
    }
    let mut geometry = native_healed_intersection_curve(candidate);
    geometry.transform_by(Matrix4::from_scale(inverse_scale));
    let mut fragment = String::new();
    std::fmt::write(
        &mut fragment,
        format_args!("{}", StepDisplay::new(&geometry, 1)),
    )
    .map_err(|_| anyhow!("the analytic repaired support cannot be represented in STEP"))?;
    let text = format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('analytic healed geometry'),'2;1');\n\
FILE_NAME('geometry','',(''),(''),'','','');\n\
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{fragment}\
ENDSEC;\nEND-ISO-10303-21;\n"
    );
    let parsed = step_p21::parser::parse(&text)?;
    let [section] = parsed.data.as_slice() else {
        bail!("native analytic geometry must contain exactly one DATA section");
    };
    let index = build_index(&section.entities);
    if index.len() != section.entities.len() {
        bail!("native generated geometry has duplicate entity IDs");
    }
    // The StepFormat writer guarantees the first local entity is the
    // intersection, but all surface associations are discovered from
    // the typed STEP entity graph, not from hardcoded local ID offsets.
    let intersection_ids = section
        .entities
        .iter()
        .filter_map(|entity| {
            let EntityInstance::Simple { id, record } = entity else {
                return None;
            };
            (record.name == "INTERSECTION_CURVE").then_some(*id)
        })
        .collect::<Vec<_>>();
    let [intersection_id] = intersection_ids.as_slice() else {
        bail!("native geometry must contain exactly one intersection curve");
    };
    let parameters = named_entity_parameters(
        *intersection_id,
        "INTERSECTION_CURVE",
        &section.entities,
        &index,
    )
    .ok_or_else(|| anyhow!("native intersection entity missing"))?;
    let Some(Parameter::List(associated)) = parameters.get(2) else {
        bail!("native intersection lacks source surface associations");
    };
    let [plane, cone] = associated.as_slice() else {
        bail!("native analytic intersection must have exactly two supports");
    };
    let plane_id = reference(plane)?;
    let cone_id = reference(cone)?;
    named_local_id(&section.entities, &index, plane_id, "PLANE")?;
    named_local_id(&section.entities, &index, cone_id, "CONICAL_SURFACE")?;
    Ok((
        section.entities.clone(),
        *intersection_id,
        plane_id,
        cone_id,
    ))
}

fn allocate_native_entities(
    section: &DataSection,
    native: Vec<EntityInstance>,
    intersection_id: u64,
    plane_id: u64,
    cone_id: u64,
    original_cone_id: u64,
) -> Result<(Vec<EntityInstance>, u64, u64)> {
    let mut high = section.entities.iter().map(entity_id).max().unwrap_or(0);
    let mut aliases = HashMap::new();
    let index = build_index(&native);
    if index.len() != native.len() {
        bail!("ambiguous native geometry IDs");
    }
    for entity in &native {
        high = high
            .checked_add(1)
            .ok_or_else(|| anyhow!("STEP entity namespace exhausted"))?;
        if aliases.insert(entity_id(entity), high).is_some() {
            bail!("new analytic entity has duplicate source ID");
        }
    }
    let new_intersection = *aliases
        .get(&intersection_id)
        .ok_or_else(|| anyhow!("native intersection ID was not allocated"))?;
    let new_plane = *aliases
        .get(&plane_id)
        .ok_or_else(|| anyhow!("native corrected plane ID was not allocated"))?;
    let new_cone = *aliases
        .get(&cone_id)
        .ok_or_else(|| anyhow!("native cone ID was not allocated"))?;
    let mut remapped = Vec::with_capacity(native.len());
    for mut entity in native {
        let old_id = entity_id(&entity);
        let replacement = aliases[&old_id];
        match &mut entity {
            EntityInstance::Simple { id, .. } | EntityInstance::Complex { id, .. } => {
                *id = replacement;
            }
        }
        rewrite_entity_refs(&mut entity, &aliases);
        remapped.push(entity);
    }
    let mut remapped_section = DataSection {
        meta: vec![],
        entities: remapped,
    };
    let map = build_index(&remapped_section.entities);
    let parameters = simple_parameters_mut(
        &mut remapped_section,
        &map,
        new_intersection,
        "INTERSECTION_CURVE",
    )?;
    let Some(Parameter::List(associated)) = parameters.get_mut(2) else {
        bail!("new intersection has lost its supporting surface list");
    };
    let [plane, cone] = associated.as_mut_slice() else {
        bail!("new intersection no longer has exactly two supports");
    };
    if reference(plane)? != new_plane || reference(cone)? != new_cone {
        bail!("native surface association order no longer matches analytic source");
    }
    *cone = Parameter::Ref(Name::Entity(original_cone_id));

    // The extra generated conical carrier is disconnected after linking the
    // original source cone. Remove *all* now-unreachable generated support
    // records and their dependent frame/point nodes, retaining the rest
    // by an exact referenced-entity closure from the analytic master.
    let local_map = build_index(&remapped_section.entities);
    let reachable = closure_from(new_intersection, &remapped_section.entities, &local_map);
    let retained: Vec<_> = remapped_section
        .entities
        .into_iter()
        .filter(|entity| reachable.contains(&entity_id(entity)))
        .collect();
    if retained.is_empty() || !retained.iter().any(|e| entity_id(e) == new_plane) {
        bail!("closed source geometry association did not retain its corrected plane");
    }
    let retained_ids = retained.iter().map(entity_id).collect::<HashSet<_>>();
    let original_ids = section
        .entities
        .iter()
        .map(entity_id)
        .collect::<HashSet<_>>();
    for entity in &retained {
        let mut bad_ref = None;
        visit_entity_refs(entity, &mut |id| {
            if !retained_ids.contains(&id) && !original_ids.contains(&id) {
                bad_ref = Some(id);
            }
        });
        if let Some(bad_ref) = bad_ref {
            bail!("new analytic source geometry has dangling reference #{bad_ref}");
        }
    }
    Ok((retained, new_intersection, new_plane))
}

/// The two original topology records are the only existing records edited.
#[derive(Clone, Copy)]
struct TopologyRewrite {
    edge_id: u64,
    plane_face_id: u64,
    original_basis_curve_id: u64,
    new_intersection_id: u64,
    original_plane_surface_id: u64,
    new_plane_surface_id: u64,
    original_same_sense: bool,
    healed_same_sense: bool,
}

fn update_existing_source_topology(
    section: &mut DataSection,
    index: &HashMap<u64, usize>,
    rewrite: TopologyRewrite,
) -> Result<()> {
    let edge = simple_parameters_mut(section, index, rewrite.edge_id, "EDGE_CURVE")?;
    if edge.len() != 5
        || reference(&edge[3])? != rewrite.original_basis_curve_id
        || !matches!(&edge[4], Parameter::Enumeration(s) if
            (s == "T") == rewrite.original_same_sense)
    {
        bail!("source edge changed during analytical repair");
    }
    edge[3] = Parameter::Ref(Name::Entity(rewrite.new_intersection_id));
    edge[4] = Parameter::Enumeration(if rewrite.healed_same_sense { "T" } else { "F" }.into());
    let face = simple_parameters_mut(section, index, rewrite.plane_face_id, "ADVANCED_FACE")?;
    if face.len() != 4 || reference(&face[2])? != rewrite.original_plane_surface_id {
        bail!("source plane face changed during analytical repair");
    }
    face[2] = Parameter::Ref(Name::Entity(rewrite.new_plane_surface_id));
    Ok(())
}

/// A proof prepared exclusively from the original input, before any edits.
struct PreparedSourceRepair {
    candidate: crate::source_edge_healing::HyperbolaHealingCandidate,
    source_scale_mm: f64,
    original_plane: u64,
    original_cone: u64,
    native: Vec<EntityInstance>,
    intersection_id: u64,
    plane_id: u64,
    cone_id: u64,
}

fn prepare_source_repair(
    source: &DataSection,
    request: FullSourceHyperbolaRepair,
) -> Result<PreparedSourceRepair> {
    let index = build_index(&source.entities);
    if index.len() != source.entities.len() {
        bail!("original source section has duplicated entity IDs");
    }
    let scale = source_solid_length_scale_mm(request.solid_id, &source.entities)
        .with_context(|| format!("source scale for solid #{}", request.solid_id))?;
    let candidate = propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
        solid_id: request.solid_id,
        source: StepEdge {
            id: request.edge_id,
            entities: &source.entities,
        },
        source_plane_face_id: request.plane_face_id,
        source_cone_face_id: request.cone_face_id,
        source_coordinate_scale_mm: scale,
        endpoint_tolerance_mm: request.endpoint_tolerance_mm,
        target_gap_mm: request.target_gap_mm,
    })?;
    crate::healed_edge_topology::healed_edge_with_source_vertices(
        &candidate,
        &source.entities,
        request.endpoint_tolerance_mm,
    )?;
    let original_plane = brep::face_surface(request.plane_face_id, &source.entities, &index)
        .ok_or_else(|| anyhow!("original source face has no plane"))?;
    let original_cone = brep::face_surface(request.cone_face_id, &source.entities, &index)
        .ok_or_else(|| anyhow!("original source face has no cone"))?;
    named_local_id(&source.entities, &index, original_plane, "PLANE")?;
    named_local_id(&source.entities, &index, original_cone, "CONICAL_SURFACE")?;
    let (native, intersection_id, plane_id, cone_id) =
        analytic_geometry_with_local_refs(&candidate, scale)?;
    Ok(PreparedSourceRepair {
        candidate,
        source_scale_mm: scale,
        original_plane,
        original_cone,
        native,
        intersection_id,
        plane_id,
        cone_id,
    })
}

fn apply_prepared_repair(
    changed: &mut DataSection,
    request: FullSourceHyperbolaRepair,
    prepared: PreparedSourceRepair,
) -> Result<FullSourceHyperbolaRepairEvidence> {
    let PreparedSourceRepair {
        candidate,
        source_scale_mm: scale,
        original_plane,
        original_cone,
        native,
        intersection_id,
        plane_id,
        cone_id,
    } = prepared;
    let (extra, new_intersection, new_plane) = allocate_native_entities(
        changed,
        native,
        intersection_id,
        plane_id,
        cone_id,
        original_cone,
    )?;
    let original_index = build_index(&changed.entities);
    update_existing_source_topology(
        changed,
        &original_index,
        TopologyRewrite {
            edge_id: request.edge_id,
            plane_face_id: request.plane_face_id,
            original_basis_curve_id: candidate.evidence.original_basis_curve_id,
            new_intersection_id: new_intersection,
            original_plane_surface_id: original_plane,
            new_plane_surface_id: new_plane,
            original_same_sense: candidate.evidence.source_curve_same_sense,
            healed_same_sense: candidate.trim.forward(),
        },
    )?;
    let new_geometry_entity_ids = extra.iter().map(entity_id).collect::<Vec<_>>();
    changed.entities.extend(extra);
    let current_index = build_index(&changed.entities);
    if current_index.len() != changed.entities.len() {
        bail!("repaired source geometry has collided entity IDs");
    }
    // Verify the original solid's manifold face set is unchanged.
    let now_faces = brep::solid_face_ids(request.solid_id, &changed.entities, &current_index)
        .ok_or_else(|| anyhow!("healed source solid lost its original face topology"))?;
    if !now_faces.contains(&request.plane_face_id) || !now_faces.contains(&request.cone_face_id) {
        bail!("analytic healed source is no longer incident on requested solid");
    }
    Ok(FullSourceHyperbolaRepairEvidence {
        data_section: request.data_section,
        original_solid_id: request.solid_id,
        original_edge_id: request.edge_id,
        original_basis_curve_id: candidate.evidence.original_basis_curve_id,
        original_plane_face_id: request.plane_face_id,
        original_cone_face_id: request.cone_face_id,
        original_plane_surface_id: original_plane,
        original_cone_surface_id: original_cone,
        new_intersection_curve_id: new_intersection,
        new_plane_surface_id: new_plane,
        new_geometry_entity_ids,
        source_length_scale_mm: scale,
        original_curve_same_sense: candidate.evidence.source_curve_same_sense,
        healed_curve_same_sense: candidate.trim.forward(),
        source_to_healed_locus_separation_lower_mm: candidate
            .evidence
            .minimum_source_to_healed_locus_separation_mm,
        healed_source_plane_residual_upper_mm: candidate
            .evidence
            .healed_source_plane_max_residual_mm,
    })
}

/// Construct a separately healed STEP document without overwriting its source.
///
/// Every requested edge is proved against the unmodified source before any
/// mutation. New analytic curves and corrected planes are inserted into the
/// same DATA section, while original edge IDs, vertices, coedge directions,
/// surrounding shell/solid topology, old splines and metadata are retained.
///
/// # Errors
/// Rejects empty repair lists, duplicate/conflicting requested replacements,
/// invalid or unsupported source solids, source unit/uncertainty mismatch,
/// insufficient globally certified original-spline departure, invalid new
/// analytic carrier geometry and malformed references.
pub fn heal_step_document(
    bytes: &[u8],
    requests: &[FullSourceHyperbolaRepair],
) -> Result<FullSourceHyperbolaHealing> {
    if requests.is_empty() {
        bail!("no requested analytic source repairs");
    }
    let parsed = ParsedStepSource::parse(bytes)?;
    crate::step_io::require_supported_sections(&parsed.exchange)?;
    let mut revised = parsed.exchange.clone();
    let mut seen_edges = HashSet::new();
    let mut seen_faces = HashSet::new();
    let mut manifest = Vec::with_capacity(requests.len());

    for request in requests {
        if !seen_edges.insert((request.data_section, request.edge_id)) {
            bail!("duplicate analytic repair for the same source edge");
        }
        if !seen_faces.insert((request.data_section, request.plane_face_id)) {
            bail!("more than one corrected analytic plane requested on one face");
        }
        // Proof always uses the unmodified source, never earlier repairs.
        let source = parsed
            .exchange
            .data
            .get(request.data_section)
            .ok_or_else(|| anyhow!("requested original STEP DATA section absent"))?;
        let prepared = prepare_source_repair(source, *request)?;
        let changed = revised
            .data
            .get_mut(request.data_section)
            .ok_or_else(|| anyhow!("new document source DATA section absent"))?;
        manifest.push(apply_prepared_repair(changed, *request, prepared)?);
    }
    let text = write_exchange(&revised)?;
    let reparsed = ParsedStepSource::parse(text.as_bytes())?;
    if reparsed.exchange.data.len() != parsed.exchange.data.len() {
        bail!("serialized analytic repair changed DATA section count");
    }
    Ok(FullSourceHyperbolaHealing {
        step: text,
        repairs: manifest,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_entities::vertex_point;
    use monstertruck_geometry::prelude::{InnerSpace, ParametricCurve};
    use monstertruck_io::step::load::{EdgeCurveHolder, Table};
    use step_p21::tables::EntityTable;

    fn fixture() -> &'static [u8] {
        include_bytes!("../validation/fixtures/production_radial_slot.step")
    }

    fn request(edge_id: u64) -> FullSourceHyperbolaRepair {
        let (plane_face_id, cone_face_id) = match edge_id {
            93481 => (125603, 58899),
            120654 => (75450, 36791),
            _ => unreachable!("fixture has only two deliberately inconsistent source cubics"),
        };
        FullSourceHyperbolaRepair {
            data_section: 0,
            solid_id: 15385,
            edge_id,
            plane_face_id,
            cone_face_id,
            endpoint_tolerance_mm: 1.0e-8,
            target_gap_mm: 1.0e-9,
        }
    }

    /// Compare the original and healed versions of one STEP topology entity.
    fn original_and_repaired<'a>(
        id: u64,
        kind: &str,
        original: &'a [EntityInstance],
        original_index: &HashMap<u64, usize>,
        revised: &'a [EntityInstance],
        revised_index: &HashMap<u64, usize>,
    ) -> Result<(&'a [Parameter], &'a [Parameter])> {
        let original = named_entity_parameters(id, kind, original, original_index)
            .ok_or_else(|| anyhow!("source #{id} is not a {kind}"))?;
        let repaired = named_entity_parameters(id, kind, revised, revised_index)
            .ok_or_else(|| anyhow!("healed #{id} is not a {kind}"))?;
        Ok((original, repaired))
    }

    #[test]
    fn both_independently_proved_repairs_keep_real_original_topology_and_vertices() -> Result<()> {
        let original = ParsedStepSource::parse(fixture())?;
        let [original_section] = original.exchange.data.as_slice() else {
            bail!("expected one source DATA section");
        };
        let original_index = build_index(&original_section.entities);
        for repairs in [
            vec![request(93481)],
            vec![request(120654)],
            vec![request(93481), request(120654)],
            vec![request(120654), request(93481)],
        ] {
            let result = heal_step_document(fixture(), &repairs)?;
            assert_eq!(result.repairs.len(), repairs.len());
            let revised = ParsedStepSource::parse(result.step.as_bytes())?;
            let [section] = revised.exchange.data.as_slice() else {
                bail!("output lost original DATA section");
            };
            assert_eq!(revised.exchange.header, original.exchange.header);
            assert_eq!(section.meta, original_section.meta);
            let index = build_index(&section.entities);
            assert_eq!(index.len(), section.entities.len());
            let repaired_ids = repairs.iter().map(|r| r.edge_id).collect::<HashSet<_>>();
            let repaired_faces = repairs
                .iter()
                .map(|r| r.plane_face_id)
                .collect::<HashSet<_>>();
            for entity in &original_section.entities {
                let id = entity_id(entity);
                let current = &section.entities[index[&id]];
                if !repaired_ids.contains(&id) && !repaired_faces.contains(&id) {
                    assert_eq!(current, entity, "source #{id} unexpectedly changed");
                }
            }
            let table = Table::from_step(&result.step)?;
            for detail in &result.repairs {
                let (old_edge, new_edge) = original_and_repaired(
                    detail.original_edge_id,
                    "EDGE_CURVE",
                    &original_section.entities,
                    &original_index,
                    &section.entities,
                    &index,
                )?;
                assert_eq!(old_edge[1], new_edge[1]);
                assert_eq!(old_edge[2], new_edge[2]);
                assert_eq!(reference(&old_edge[3])?, detail.original_basis_curve_id);
                assert_eq!(reference(&new_edge[3])?, detail.new_intersection_curve_id);
                assert_eq!(new_edge[4], Parameter::Enumeration("F".to_owned()));
                assert!(detail.original_curve_same_sense);
                assert!(!detail.healed_curve_same_sense);
                let (old_plane, new_plane) = original_and_repaired(
                    detail.original_plane_face_id,
                    "ADVANCED_FACE",
                    &original_section.entities,
                    &original_index,
                    &section.entities,
                    &index,
                )?;
                assert_eq!(old_plane[1], new_plane[1]);
                assert_eq!(old_plane[3], new_plane[3]);
                assert_eq!(reference(&new_plane[2])?, detail.new_plane_surface_id);
                assert_eq!(reference(&old_plane[2])?, detail.original_plane_surface_id);
                let new_int = named_entity_parameters(
                    detail.new_intersection_curve_id,
                    "INTERSECTION_CURVE",
                    &section.entities,
                    &index,
                )
                .ok_or_else(|| anyhow!("native analytic intersection missing"))?;
                let Some(Parameter::List(associated)) = new_int.get(2) else {
                    bail!("new intersection has no typed associated surfaces");
                };
                assert_eq!(associated.len(), 2);
                assert_eq!(reference(&associated[0])?, detail.new_plane_surface_id);
                assert_eq!(reference(&associated[1])?, detail.original_cone_surface_id);
                assert!(
                    detail
                        .new_geometry_entity_ids
                        .contains(&detail.new_intersection_curve_id)
                );
                assert!(
                    !detail
                        .new_geometry_entity_ids
                        .contains(&detail.original_cone_surface_id)
                );
                let owned = <Table as EntityTable<EdgeCurveHolder>>::get_owned(
                    &table,
                    detail.original_edge_id,
                )?;
                assert!(!owned.same_sense);
                let curve = owned
                    .parse_curve3d()
                    .map_err(|error| anyhow!("real source repaired native edge: {error}"))?;
                let source_start_id = reference(&old_edge[1])?;
                let source_end_id = reference(&old_edge[2])?;
                let start =
                    vertex_point(source_start_id, &original_section.entities, &original_index)
                        .ok_or_else(|| anyhow!("source start vertex missing"))?;
                let end = vertex_point(source_end_id, &original_section.entities, &original_index)
                    .ok_or_else(|| anyhow!("source end vertex missing"))?;
                let (a, b) = curve.parameter_range();
                let (std::ops::Bound::Included(a), std::ops::Bound::Included(b)) = (a, b) else {
                    bail!("source healed native edge has no finite parameter bounds");
                };
                let point = |xyz: [f64; 3]| {
                    monstertruck_geometry::prelude::Point3::new(xyz[0], xyz[1], xyz[2])
                };
                assert!((curve.evaluate(a) - point(start)).magnitude() < 1.0e-8);
                assert!((curve.evaluate(b) - point(end)).magnitude() < 1.0e-8);
                assert!(detail.source_to_healed_locus_separation_lower_mm > 1.0e-5);
            }
        }
        Ok(())
    }

    #[test]
    fn changed_and_conflicting_requests_fail_without_returning_a_partial_document() -> Result<()> {
        assert!(heal_step_document(fixture(), &[]).is_err());
        assert!(heal_step_document(fixture(), &[request(93481), request(93481)]).is_err());
        let mut wrong_face = request(93481);
        wrong_face.plane_face_id = 75450;
        assert!(heal_step_document(fixture(), &[wrong_face]).is_err());
        let mut invalid = request(93481);
        invalid.endpoint_tolerance_mm = f64::NAN;
        assert!(heal_step_document(fixture(), &[invalid]).is_err());
        invalid = request(93481);
        invalid.data_section = 1;
        assert!(heal_step_document(fixture(), &[invalid]).is_err());
        Ok(())
    }

    #[test]
    fn geometry_placement_scales_to_source_coordinates_without_changing_angular_units() -> Result<()>
    {
        let input = ParsedStepSource::parse(fixture())?;
        let source = &input.exchange.data[0].entities;
        let candidate = propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
            solid_id: 15385,
            source: StepEdge {
                id: 93481,
                entities: source,
            },
            source_plane_face_id: 125603,
            source_cone_face_id: 58899,
            source_coordinate_scale_mm: 1.0,
            endpoint_tolerance_mm: 1.0e-8,
            target_gap_mm: 1.0e-9,
        })?;
        assert!(
            analytic_geometry_with_local_refs(&candidate, f64::MIN_POSITIVE / 16.0).is_err(),
            "non-finite inverse STEP length transform must fail without panicking"
        );
        let (mm, root_mm, _, _) = analytic_geometry_with_local_refs(&candidate, 1.0)?;
        let (cm, root_cm, _, _) = analytic_geometry_with_local_refs(&candidate, 10.0)?;
        let mm_table = Table::from_data_section(&DataSection {
            meta: vec![],
            entities: mm,
        });
        let cm_table = Table::from_data_section(&DataSection {
            meta: vec![],
            entities: cm,
        });
        let mm_holder =
            <Table as EntityTable<monstertruck_io::step::load::SurfaceCurveHolder>>::get_owned(
                &mm_table, root_mm,
            )?;
        let cm_holder =
            <Table as EntityTable<monstertruck_io::step::load::SurfaceCurveHolder>>::get_owned(
                &cm_table, root_cm,
            )?;
        let curve_mm = monstertruck_io::step::load::step_geometry::Curve3D::try_from(&mm_holder)
            .map_err(|err| anyhow!("mm original analytic geometry conversion: {err}"))?;
        let curve_cm = monstertruck_io::step::load::step_geometry::Curve3D::try_from(&cm_holder)
            .map_err(|err| anyhow!("cm original analytic geometry conversion: {err}"))?;
        let millimetres = curve_mm.evaluate(0.5);
        let centimetres = curve_cm.evaluate(0.5);
        let leader_point = |curve: &monstertruck_io::step::load::step_geometry::Curve3D|
            -> Result<monstertruck_geometry::prelude::Point3> {
            let monstertruck_io::step::load::step_geometry::Curve3D::IntersectionCurve(geometry) = curve
            else {
                bail!("native source STEP lost its two-surface intersection");
            };
            Ok(geometry.leader().evaluate(0.5))
        };
        let mm_leader = leader_point(&curve_mm)?;
        let cm_leader = leader_point(&curve_cm)?;
        let expected = candidate
            .analytic
            .evaluate(0.5)
            .ok_or_else(|| anyhow!("original hyperbola out of domain"))?;
        let discrepancy = (millimetres - centimetres * 10.0).magnitude();
        let leader_discrepancy = (mm_leader - cm_leader * 10.0).magnitude();
        assert!((mm_leader - expected).magnitude() < 1.0e-9);
        assert!(leader_discrepancy < 1.0e-9);
        assert!(
            discrepancy < 1.0e-9,
            "different source unit placements mm={millimetres:?} cm={centimetres:?} discrepancy={discrepancy}"
        );
        Ok(())
    }
}
