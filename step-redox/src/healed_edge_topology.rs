//! Finite STEP edge geometry and two opposing source coedges.
//!
//! The original manufacturer spline remains untouched. This creates a
//! separate open topology fragment: two vertices, one edge, two coedge uses.
//! No face loops, manifold shell or restored solid are asserted.

use crate::{
    native_hyperbola_intersection::native_healed_intersection_curve,
    source_edge_healing::{HyperbolaHealingCandidate, unique_original_face_edge_use},
    step_entities::{enumeration_bool, named_entity_parameters, vertex_point},
    step_graph::{build_index, entity_id, entity_ref_value},
    step_io::format_real,
};
use anyhow::{Result, anyhow, bail};
use monstertruck_geometry::prelude::{InnerSpace, Point3};
use monstertruck_io::step::save::StepDisplay;
use serde::Serialize;
use step_p21::ast::EntityInstance;

/// Fresh local entity IDs, with original IDs retained in healing evidence.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
pub struct HealedEdgeEntityIds {
    pub intersection_curve: u64,
    pub start_point: u64,
    pub start_vertex: u64,
    pub end_point: u64,
    pub end_vertex: u64,
    pub edge_curve: u64,
    pub plane_oriented_edge: u64,
    pub cone_oriented_edge: u64,
}

/// A new finite, source-vertex-anchored open edge, not a closed solid.
#[derive(Clone, Debug, Serialize)]
pub struct HealedEdgeTopologyFragment {
    pub ids: HealedEdgeEntityIds,
    /// The native hyperbola sense is independent of the source spline sense.
    pub healed_curve_same_sense: bool,
    pub start_residual_mm: f64,
    pub end_residual_mm: f64,
    pub native_step_fragment: String,
}

fn original_vertex_point(
    id: u64,
    source: &[EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
    scale_mm: f64,
) -> Result<Point3> {
    let xyz =
        vertex_point(id, source, index).ok_or_else(|| anyhow!("original vertex #{id} missing"))?;
    let p = Point3::new(xyz[0] * scale_mm, xyz[1] * scale_mm, xyz[2] * scale_mm);
    if [p.x, p.y, p.z].into_iter().any(|x| !x.is_finite()) {
        bail!("source vertex coordinates are not finite");
    }
    Ok(p)
}

fn original_coedge_orientation(
    id: u64,
    edge: u64,
    source: &[EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
) -> Result<bool> {
    let params = named_entity_parameters(id, "ORIENTED_EDGE", source, index)
        .ok_or_else(|| anyhow!("original coedge #{id} missing"))?;
    if params.len() != 5 || params.get(3).and_then(entity_ref_value) != Some(edge) {
        bail!("original coedge does not reference its nominated edge");
    }
    params
        .get(4)
        .and_then(enumeration_bool)
        .ok_or_else(|| anyhow!("original coedge lacks an orientation"))
}

fn next_id(id: u64) -> Result<u64> {
    id.checked_add(1)
        .ok_or_else(|| anyhow!("STEP entity ID overflow"))
}

#[derive(Clone, Copy)]
struct VerifiedSourceVertices {
    start: Point3,
    end: Point3,
    plane_orientation: bool,
    cone_orientation: bool,
    start_residual_mm: f64,
    end_residual_mm: f64,
}

fn verify_source_vertices(
    candidate: &HyperbolaHealingCandidate,
    source: &[EntityInstance],
    endpoint_tolerance_mm: f64,
) -> Result<VerifiedSourceVertices> {
    let evidence = &candidate.evidence;
    let scale = evidence.source_coordinate_scale_mm;
    if !scale.is_finite()
        || scale <= 0.0
        || !endpoint_tolerance_mm.is_finite()
        || endpoint_tolerance_mm <= 0.0
    {
        bail!("source unit scale and endpoint error budget must be positive finite");
    }
    let declared_scale =
        crate::source_edge_locus::source_solid_length_scale_mm(evidence.original_solid_id, source)?;
    let scale_resolution = declared_scale.max(scale) * 128.0 * f64::EPSILON;
    if (declared_scale - scale).abs() > scale_resolution {
        bail!("source STEP length context disagrees with healing candidate coordinates");
    }
    let index = build_index(source);
    if index.len() != source.len() {
        bail!("original source contains ambiguous duplicate STEP entity IDs");
    }
    let original = named_entity_parameters(evidence.original_edge_id, "EDGE_CURVE", source, &index)
        .ok_or_else(|| anyhow!("original edge was not found"))?;
    if original.len() != 5
        || original.get(1).and_then(entity_ref_value) != Some(evidence.source_start_vertex_id)
        || original.get(2).and_then(entity_ref_value) != Some(evidence.source_end_vertex_id)
        || original.get(3).and_then(entity_ref_value) != Some(evidence.original_basis_curve_id)
        || original.get(4).and_then(enumeration_bool) != Some(evidence.source_curve_same_sense)
    {
        bail!("original edge changed since construction of healing candidate");
    }
    let source_plane_coedge = unique_original_face_edge_use(
        evidence.source_plane_face_id,
        evidence.original_edge_id,
        source,
        &index,
    )?;
    let source_cone_coedge = unique_original_face_edge_use(
        evidence.source_cone_face_id,
        evidence.original_edge_id,
        source,
        &index,
    )?;
    if source_plane_coedge != evidence.source_plane_coedge
        || source_cone_coedge != evidence.source_cone_coedge
        || source_plane_coedge.parameter_forward == source_cone_coedge.parameter_forward
    {
        bail!("the original neighboring face-loop incidences changed after the healing proof");
    }
    let plane_orientation = original_coedge_orientation(
        evidence.source_plane_coedge.oriented_edge_id,
        evidence.original_edge_id,
        source,
        &index,
    )?;
    let cone_orientation = original_coedge_orientation(
        evidence.source_cone_coedge.oriented_edge_id,
        evidence.original_edge_id,
        source,
        &index,
    )?;
    let start = original_vertex_point(evidence.source_start_vertex_id, source, &index, scale)?;
    let end = original_vertex_point(evidence.source_end_vertex_id, source, &index, scale)?;
    let new_start = candidate
        .analytic
        .evaluate(candidate.trim.start_parameter)
        .ok_or_else(|| anyhow!("healed curve lacks a valid starting point"))?;
    let new_end = candidate
        .analytic
        .evaluate(candidate.trim.end_parameter)
        .ok_or_else(|| anyhow!("healed curve lacks a valid ending point"))?;
    let start_residual_mm = (new_start - start).magnitude();
    let end_residual_mm = (new_end - end).magnitude();
    if !start_residual_mm.is_finite()
        || !end_residual_mm.is_finite()
        || start_residual_mm > endpoint_tolerance_mm
        || end_residual_mm > endpoint_tolerance_mm
    {
        bail!("original source vertices do not lie within healed endpoint budget");
    }

    Ok(VerifiedSourceVertices {
        start,
        end,
        plane_orientation,
        cone_orientation,
        start_residual_mm,
        end_residual_mm,
    })
}

fn render_open_edge(
    candidate: &HyperbolaHealingCandidate,
    vertices: VerifiedSourceVertices,
) -> Result<HealedEdgeTopologyFragment> {
    let VerifiedSourceVertices {
        start,
        end,
        plane_orientation,
        cone_orientation,
        start_residual_mm,
        end_residual_mm,
    } = vertices;
    let t0 = candidate.trim.start_parameter;
    let t1 = candidate.trim.end_parameter;
    let parameter_resolution = 32.0 * f64::EPSILON * t0.abs().max(t1.abs()).max(1.0);
    if !t0.is_finite() || !t1.is_finite() || (t0 - t1).abs() <= parameter_resolution {
        bail!("healed curve requires distinct finite native trim parameters");
    }
    let healed_curve_same_sense = t0 < t1;
    let geometry = StepDisplay::new(&native_healed_intersection_curve(candidate), 1).to_string();
    let exchange = step_p21::parser::parse(&format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('open edge'),'2;1');\n\
FILE_NAME('geometry','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n\
ENDSEC;\nDATA;\n{geometry}\nENDSEC;\nEND-ISO-10303-21;\n"
    ))?;
    let [section] = exchange.data.as_slice() else {
        bail!("expected exactly one STEP geometry DATA section");
    };
    let max_id = section
        .entities
        .iter()
        .map(entity_id)
        .max()
        .ok_or_else(|| anyhow!("native intersection has no STEP geometry records"))?;
    let start_point = next_id(max_id)?;
    let start_vertex = next_id(start_point)?;
    let end_point = next_id(start_vertex)?;
    let end_vertex = next_id(end_point)?;
    let edge_curve = next_id(end_vertex)?;
    let plane_oriented_edge = next_id(edge_curve)?;
    let cone_oriented_edge = next_id(plane_oriented_edge)?;
    let ids = HealedEdgeEntityIds {
        intersection_curve: 1,
        start_point,
        start_vertex,
        end_point,
        end_vertex,
        edge_curve,
        plane_oriented_edge,
        cone_oriented_edge,
    };
    let flag = |b: bool| if b { ".T." } else { ".F." };
    let xyz = |p: Point3| {
        format!(
            "({},{},{})",
            format_real(p.x),
            format_real(p.y),
            format_real(p.z)
        )
    };
    let native_step_fragment = format!(
        "{geometry}\
#{start_point}=CARTESIAN_POINT('original start',{});\n\
#{start_vertex}=VERTEX_POINT('source start',#{start_point});\n\
#{end_point}=CARTESIAN_POINT('original end',{});\n\
#{end_vertex}=VERTEX_POINT('source end',#{end_point});\n\
#{edge_curve}=EDGE_CURVE('new analytic edge',#{start_vertex},#{end_vertex},#1,{});\n\
#{plane_oriented_edge}=ORIENTED_EDGE('plane incidence',*,*,#{edge_curve},{});\n\
#{cone_oriented_edge}=ORIENTED_EDGE('cone incidence',*,*,#{edge_curve},{});\n",
        xyz(start),
        xyz(end),
        flag(healed_curve_same_sense),
        flag(plane_orientation),
        flag(cone_orientation),
    );
    Ok(HealedEdgeTopologyFragment {
        ids,
        healed_curve_same_sense,
        start_residual_mm,
        end_residual_mm,
        native_step_fragment,
    })
}

/// Emit an open edge with physically original source vertex order and
/// source-authored opposing coedge orientations. Every output coordinate
/// is in millimetres, even when source unit scale is different.
///
/// The new curve sense follows native hyperbola trim direction; for the
/// two manufacturer radial edges it is false despite their original
/// cubic edge sense being true.
///
/// # Errors
/// Refuses invalid source provenance, changed original topology, endpoint
/// departure above the explicit tolerance or invalid native geometry.
pub fn healed_edge_with_source_vertices(
    candidate: &HyperbolaHealingCandidate,
    source: &[EntityInstance],
    endpoint_tolerance_mm: f64,
) -> Result<HealedEdgeTopologyFragment> {
    let verified = verify_source_vertices(candidate, source, endpoint_tolerance_mm)?;
    render_open_edge(candidate, verified)
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Result, bail};
    use monstertruck_geometry::prelude::ParametricCurve;
    use monstertruck_io::step::load::{EdgeCurveHolder, OrientedEdgeHolder, Table};
    use step_p21::{ast::Parameter, tables::EntityTable};

    fn original_manufacturer() -> Result<Vec<(HyperbolaHealingCandidate, Vec<EntityInstance>)>> {
        let text = std::str::from_utf8(include_bytes!(
            "../validation/fixtures/production_radial_slot.step"
        ))?;
        let exchange = step_p21::parser::parse(text)?;
        let [section] = exchange.data.as_slice() else {
            bail!("manufacturer source has no single DATA section");
        };
        let source = &section.entities;
        [(93481, 125603, 58899), (120654, 75450, 36791)]
            .into_iter()
            .map(|(edge_id, plane, cone)| {
                let candidate = crate::source_edge_healing::propose_plane_cone_hyperbola_healing(
                    crate::source_edge_healing::HyperbolaHealingRequest {
                        solid_id: 15385,
                        source: crate::source_edge_locus::StepEdge {
                            id: edge_id,
                            entities: source,
                        },
                        source_plane_face_id: plane,
                        source_cone_face_id: cone,
                        source_coordinate_scale_mm: 1.0,
                        endpoint_tolerance_mm: 1.0e-8,
                        target_gap_mm: 1.0e-9,
                    },
                )?;
                Ok((candidate, source.clone()))
            })
            .collect()
    }

    fn document(fragment: &str) -> String {
        format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('new open edge'),'2;1');\n\
FILE_NAME('edge','2026-11-18',(''),(''),'','','');\n\
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{fragment}\
ENDSEC;\nEND-ISO-10303-21;\n"
        )
    }

    #[test]
    fn real_source_vertices_and_opposing_coedges_survive_bounded_edge_import() -> Result<()> {
        for (candidate, source) in original_manufacturer()? {
            let healed = healed_edge_with_source_vertices(&candidate, &source, 1.0e-8)?;
            assert!(candidate.evidence.source_curve_same_sense);
            assert!(!healed.healed_curve_same_sense);
            assert!(healed.start_residual_mm < 1.0e-11);
            assert!(healed.end_residual_mm < 1.0e-11);
            let full = document(&healed.native_step_fragment);
            let exchange = step_p21::parser::parse(&full)?;
            assert_eq!(exchange.data.len(), 1);
            let table = Table::from_step(&full)?;
            assert!(table.entity_report.is_empty());
            assert_eq!(table.edge_curve.len(), 1);
            assert_eq!(table.vertex_point.len(), 2);
            assert_eq!(table.oriented_edge.len(), 2);
            assert_eq!(table.surface_curve.len(), 1);
            assert_eq!(table.pcurve.len(), 0);
            assert!(table.dummy.is_empty());
            let plane_use = <Table as EntityTable<OrientedEdgeHolder>>::get_owned(
                &table,
                healed.ids.plane_oriented_edge,
            )?;
            let cone_use = <Table as EntityTable<OrientedEdgeHolder>>::get_owned(
                &table,
                healed.ids.cone_oriented_edge,
            )?;
            assert!(!plane_use.orientation);
            assert!(cone_use.orientation);
            assert_eq!(
                plane_use.edge_element.same_sense,
                healed.healed_curve_same_sense
            );
            assert_eq!(
                cone_use.edge_element.same_sense,
                healed.healed_curve_same_sense
            );

            let owned =
                <Table as EntityTable<EdgeCurveHolder>>::get_owned(&table, healed.ids.edge_curve)?;
            assert_eq!(owned.same_sense, healed.healed_curve_same_sense);
            let curve = owned
                .parse_curve3d()
                .map_err(|error| anyhow!("typed bounded edge import: {error}"))?;
            let idx = build_index(&source);
            let source_start = original_vertex_point(
                candidate.evidence.source_start_vertex_id,
                &source,
                &idx,
                candidate.evidence.source_coordinate_scale_mm,
            )?;
            let source_end = original_vertex_point(
                candidate.evidence.source_end_vertex_id,
                &source,
                &idx,
                candidate.evidence.source_coordinate_scale_mm,
            )?;
            let (low, high) = curve.parameter_range();
            let start_t = match low {
                std::ops::Bound::Included(t) => t,
                _ => bail!("native edge lost its finite parameter start"),
            };
            let end_t = match high {
                std::ops::Bound::Included(t) => t,
                _ => bail!("native edge lost its finite parameter end"),
            };
            let first = curve.evaluate(start_t);
            let last = curve.evaluate(end_t);
            assert!((first - source_start).magnitude() < 1.0e-8);
            assert!((last - source_end).magnitude() < 1.0e-8);
        }
        Ok(())
    }

    #[test]
    fn tampering_original_spline_sense_and_endpoint_budget_are_rejected() -> Result<()> {
        let mut candidates = original_manufacturer()?;
        let (candidate, source) = candidates.remove(0);
        assert!(healed_edge_with_source_vertices(&candidate, &source, 0.0).is_err());
        assert!(healed_edge_with_source_vertices(&candidate, &[], 1.0e-8).is_err());
        let mut wrong_units = candidate.clone();
        wrong_units.evidence.source_coordinate_scale_mm = 10.0;
        let unit_error = healed_edge_with_source_vertices(&wrong_units, &source, 1.0e-8)
            .expect_err("original STEP millimetre context must not be overridden");
        assert!(unit_error.to_string().contains("length context"));
        let mut tampered = source;
        let mut altered = false;
        for instance in &mut tampered {
            if let EntityInstance::Simple { id, record } = instance
                && *id == candidate.evidence.original_edge_id
                && let Parameter::List(parameters) = &mut record.parameter
                && let Some(Parameter::Enumeration(sense)) = parameters.get_mut(4)
            {
                *sense = "F".to_owned();
                altered = true;
                break;
            }
        }
        assert!(altered);
        assert!(healed_edge_with_source_vertices(&candidate, &tampered, 1.0e-8).is_err());
        Ok(())
    }

    #[test]
    fn changed_original_face_bound_orientation_is_rejected() -> Result<()> {
        let (candidate, source) = original_manufacturer()?.remove(0);
        let index = build_index(&source);
        let loops =
            crate::brep::face_loops(candidate.evidence.source_plane_face_id, &source, &index)
                .ok_or_else(|| anyhow!("source plane has no face loops"))?;
        let bound_id = loops
            .iter()
            .find(|loop_| {
                loop_
                    .edges
                    .iter()
                    .any(|edge| edge.edge_id == candidate.evidence.original_edge_id)
            })
            .ok_or_else(|| anyhow!("source edge absent from source face loops"))?
            .bound_id;
        let mut changed = source;
        let mut edited = false;
        for entity in &mut changed {
            if let EntityInstance::Simple { id, record } = entity
                && *id == bound_id
                && let Parameter::List(params) = &mut record.parameter
                && let Some(Parameter::Enumeration(orientation)) = params.get_mut(2)
            {
                *orientation = if orientation == "T" { "F" } else { "T" }.to_owned();
                edited = true;
                break;
            }
        }
        assert!(edited);
        assert!(healed_edge_with_source_vertices(&candidate, &changed, 1.0e-8).is_err());
        Ok(())
    }
}
