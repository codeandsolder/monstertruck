//! Detached typed proof-view proposals from author-precision-certified STEP caps.
//!
//! The original source and every authored vertex, face, surface and curve are
//! immutable. New private 3-D curves reuse the original surface poles; only
//! copied `EDGE_CURVE` geometric references may change. The resulting whole
//! boundary diagnostic is evidence, never an automatic recovery admission.

use super::{bspline_curve_with_knots, parse_surface, rational_bspline_curve};
use crate::{
    boundary_proof::{analytic_boundary_diagnostic, audit_source_spline_cap_precision},
    brep::{self, BSplineSupport, CurveSupport},
    math3::distance,
    step_entities::{cartesian_point, named_entity_parameters},
    step_graph::{build_index, entity_id, entity_ref_value, simple_record_mut},
};
use anyhow::{Context, Result, bail, ensure};
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use step_p21::ast::{EntityInstance, Exchange, Parameter};

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SourceCapPrecisionPatch {
    pub solid_id: u64,
    pub face_id: u64,
    pub surface_id: u64,
    pub edge_id: u64,
    pub original_curve_id: u64,
    pub private_curve_id: u64,
    pub surface_row: usize,
    pub source_to_private_locus_upper_mm: f64,
    pub reversed: bool,
}

#[derive(Debug, Serialize)]
pub struct SourceCapRechartManifest {
    pub solid_id: u64,
    pub data_section: usize,
    pub author_uncertainty_mm: f64,
    pub source_boundary_supported: bool,
    pub cap_recharted_boundary_supported: bool,
    pub proof_view_boundary_supported: bool,
    pub proof_view_boundary_reason: Option<String>,
    pub reconciled_face_senses: Vec<crate::boundary_proof::StripFaceSenseRepair>,
    pub patches: Vec<SourceCapPrecisionPatch>,
}

pub struct SourceCapRechartProposal {
    pub exchange: Exchange,
    pub manifest: SourceCapRechartManifest,
}

fn reverse_knot_data(
    multiplicities: Parameter,
    knot_values: Parameter,
) -> Option<(Parameter, Parameter)> {
    let (Parameter::List(mut mult), Parameter::List(knots)) = (multiplicities, knot_values) else {
        return None;
    };
    let lo = crate::step_entities::number(knots.first()?)?;
    let hi = crate::step_entities::number(knots.last()?)?;
    if !lo.is_finite() || !hi.is_finite() || lo >= hi {
        return None;
    }
    let reflected = knots
        .iter()
        .rev()
        .map(|k| {
            let value = (lo + hi) - crate::step_entities::number(k)?;
            value.is_finite().then_some(Parameter::Real(value))
        })
        .collect::<Option<Vec<_>>>()?;
    mult.reverse();
    Some((Parameter::List(mult), Parameter::List(reflected)))
}

fn private_surface_row_curve(
    surface_id: u64,
    row: usize,
    original: &BSplineSupport,
    bound_coordinate: f64,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    id: u64,
) -> Option<(EntityInstance, bool)> {
    let surface = parse_surface(surface_id, entities, index)?;
    let rows = surface.control_points.len();
    let columns = surface.control_points.first()?.len();
    if row >= 2 {
        return None;
    }
    let (swapped, degree, closed, mut multiplicities, mut knots, count) =
        if surface.v_degree == 1 && columns == 2 {
            (
                false,
                surface.u_degree,
                surface.u_closed,
                surface.u_multiplicities,
                surface.u_knots,
                rows,
            )
        } else if surface.u_degree == 1 && rows == 2 {
            (
                true,
                surface.v_degree,
                surface.v_closed,
                surface.v_multiplicities,
                surface.v_knots,
                columns,
            )
        } else {
            return None;
        };
    let mut ids = (0..count)
        .map(|i| {
            if swapped {
                surface.control_points[row][i]
            } else {
                surface.control_points[i][row]
            }
        })
        .collect::<Vec<_>>();
    let mut weights = surface.weights.as_ref().map(|net| {
        (0..count)
            .map(|i| {
                if swapped {
                    net[row][i].clone()
                } else {
                    net[i][row].clone()
                }
            })
            .collect::<Vec<_>>()
    });
    let first = cartesian_point(*ids.first()?, entities, index)?;
    let last = cartesian_point(*ids.last()?, entities, index)?;
    let source_start = *original.control_points_mm.first()?;
    let source_end = *original.control_points_mm.last()?;
    let direct = distance(first, source_start).max(distance(last, source_end));
    let backwards = distance(first, source_end).max(distance(last, source_start));
    // The sense must be unique under the complete original cap locus bound.
    if direct.min(backwards) > bound_coordinate
        || (direct <= bound_coordinate && backwards <= bound_coordinate)
    {
        return None;
    }
    let reversed = backwards < direct;
    if reversed {
        ids.reverse();
        if let Some(weights) = &mut weights {
            weights.reverse();
        }
        (multiplicities, knots) = reverse_knot_data(multiplicities, knots)?;
    }
    let name = surface.name;
    let self_intersect = surface.self_intersect;
    let spec = surface.knot_spec;
    let curve = if let Some(weights) = weights {
        rational_bspline_curve(
            id,
            degree,
            &ids,
            &weights,
            closed,
            self_intersect,
            multiplicities,
            knots,
            spec,
            name,
        )
    } else {
        bspline_curve_with_knots(
            id,
            degree,
            &ids,
            closed,
            self_intersect,
            multiplicities,
            knots,
            spec,
            name,
        )
    };
    Some((curve, reversed))
}

fn identify_rechart_cap_candidates(
    audit: &crate::boundary_proof::SplineCapPrecisionAudit,
    source: &[EntityInstance],
    index: &HashMap<u64, usize>,
    scale: f64,
) -> Result<BTreeMap<u64, (u64, u64, usize, f64)>> {
    // Coalesce complete cap claims by EDGE_CURVE, never by its (possibly
    // shared) geometric curve ID. Different surface-row owners are ambiguous.
    let mut candidates = BTreeMap::<u64, (u64, u64, usize, f64)>::new();
    for cap in &audit.cap_witnesses {
        let params = named_entity_parameters(cap.source_edge_id, "EDGE_CURVE", source, &index)
            .context("source EDGE_CURVE missing")?;
        let curve_id = entity_ref_value(params.get(3).context("missing edge geometry")?)
            .context("edge geometry is not an entity reference")?;
        let original_curve = source
            .get(*index.get(&curve_id).context("missing curve ID")?)
            .context("unresolvable original curve")?;
        ensure!(
            match original_curve {
                EntityInstance::Simple { record, .. } => record.name == "B_SPLINE_CURVE_WITH_KNOTS",
                EntityInstance::Complex { subsuper, .. } => subsuper
                    .0
                    .iter()
                    .any(|r| r.name == "RATIONAL_B_SPLINE_CURVE"),
            },
            "edge #{} has an unsupported shared wrapper or geometry",
            cap.source_edge_id
        );
        let CurveSupport::BSpline(current) = brep::curve_support(curve_id, source, &index) else {
            bail!(
                "edge #{} has an unsupported original B-spline",
                cap.source_edge_id
            );
        };
        let row = super::boundary_strip_profiles(cap.surface_id, source, &index)
            .and_then(|(rows, _)| rows.get(cap.surface_row).cloned())
            .context("source surface no longer has the witnessed cap row")?;
        if current == row {
            continue;
        }
        // A surface and a curve share the *literal* source knot data and
        // control points even when two independent normalized-knot parsers
        // differ in the last floating-point bit. Compare the typed STEP
        // authoring record before treating it as a geometric change.
        if private_surface_row_curve(
            cap.surface_id,
            cap.surface_row,
            &current,
            cap.bound_mm / scale,
            source,
            &index,
            curve_id,
        )
        .is_some_and(|(literal, _)| literal == *original_curve)
        {
            continue;
        }
        if let Some((other_face, other_surface, other_row, _)) = candidates.get(&cap.source_edge_id)
        {
            ensure!(
                *other_surface == cap.surface_id && *other_row == cap.surface_row,
                "ambiguous surface rows for edge #{} (faces #{} and #{})",
                cap.source_edge_id,
                other_face,
                cap.face_id
            );
            continue;
        }
        candidates.insert(
            cap.source_edge_id,
            (cap.face_id, cap.surface_id, cap.surface_row, cap.bound_mm),
        );
    }
    Ok(candidates)
}

fn construct_private_cap_curves(
    source: &[EntityInstance],
    index: &HashMap<u64, usize>,
    candidates: &BTreeMap<u64, (u64, u64, usize, f64)>,
    scale: f64,
    solid_id: u64,
) -> Result<(Vec<SourceCapPrecisionPatch>, Vec<EntityInstance>)> {
    let mut next_id = source
        .iter()
        .map(entity_id)
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .context("STEP entity IDs exhausted")?;
    let mut generated = Vec::new();
    let mut patches = Vec::new();
    for (&edge_id, &(face_id, surface_id, row, bound_mm)) in candidates {
        let params = named_entity_parameters(edge_id, "EDGE_CURVE", source, &index)
            .context("candidate edge absent")?;
        let original_curve_id =
            entity_ref_value(params.get(3).context("candidate geometry absent")?)
                .context("candidate geometry reference invalid")?;
        let CurveSupport::BSpline(current) = brep::curve_support(original_curve_id, source, &index)
        else {
            bail!("candidate edge #{edge_id} B-spline unavailable");
        };
        let (curve, reversed) = private_surface_row_curve(
            surface_id,
            row,
            &current,
            bound_mm / scale,
            source,
            &index,
            next_id,
        )
        .context("cannot author an unambiguous private surface-row curve")?;
        generated.push(curve);
        patches.push(SourceCapPrecisionPatch {
            solid_id,
            face_id,
            surface_id,
            edge_id,
            original_curve_id,
            private_curve_id: next_id,
            surface_row: row,
            source_to_private_locus_upper_mm: bound_mm,
            reversed,
        });
        next_id = next_id.checked_add(1).context("STEP curve ID overflow")?;
    }
    Ok((patches, generated))
}

fn apply_private_cap_patches(
    entities: &mut Vec<EntityInstance>,
    index: &HashMap<u64, usize>,
    patches: &[SourceCapPrecisionPatch],
    generated: Vec<EntityInstance>,
    scale: f64,
) -> Result<()> {
    for patch in patches {
        let record = simple_record_mut(
            entities
                .get_mut(*index.get(&patch.edge_id).context("candidate edge lost")?)
                .context("candidate edge out of bounds")?,
        )
        .context("candidate edge not a simple record")?;
        ensure!(record.name == "EDGE_CURVE", "candidate edge type altered");
        let Parameter::List(params) = &mut record.parameter else {
            bail!("malformed edge parameters")
        };
        ensure!(
            entity_ref_value(params.get(3).context("no original curve")?)
                == Some(patch.original_curve_id),
            "proof view edge geometry no longer equals original"
        );
        params[3] = crate::step_entities::entity_ref(patch.private_curve_id);
    }
    entities.extend(generated);
    let proof_index = build_index(entities);
    ensure!(
        proof_index.len() == entities.len(),
        "proof-view entity IDs collide"
    );
    for patch in patches {
        let CurveSupport::BSpline(private) =
            brep::curve_support(patch.private_curve_id, entities, &proof_index)
        else {
            bail!(
                "generated B-spline #{} is not parseable",
                patch.private_curve_id
            );
        };
        // Validate *literal full-curve authoring* against the retained
        // original surface row, not independently normalized float knot
        // vectors (which can differ by one final ULP despite equal STEP
        // knot parameters and hence exact mathematical NURBS identity).
        let (expected_entity, expected_reverse) = private_surface_row_curve(
            patch.surface_id,
            patch.surface_row,
            &private,
            patch.source_to_private_locus_upper_mm / scale,
            entities,
            &proof_index,
            patch.private_curve_id,
        )
        .context("reconstructed private curve has ambiguous source-row orientation")?;
        ensure!(
            expected_reverse == patch.reversed,
            "proof curve #{} changed direction",
            patch.private_curve_id
        );
        let actual = entities
            .get(
                *proof_index
                    .get(&patch.private_curve_id)
                    .context("new curve reference missing")?,
            )
            .context("new curve not found")?;
        ensure!(
            *actual == expected_entity,
            "edge #{} private curve differs from the literal surface-row STEP data",
            patch.edge_id
        );
    }
    Ok(())
}

/// Create a complete detached source AST with only private cap-curve proposals.
///
/// # Errors
/// Refuses ambiguous shared edges, malformed STEP, missing source precision,
/// incomplete original face trims, or unsupported geometry wrappers. An
/// unproved final whole-solid boundary is returned as a diagnostic, not admitted.
pub fn propose_source_cap_rechart(
    source_step: &[u8],
    data_section: usize,
    solid_id: u64,
) -> Result<SourceCapRechartProposal> {
    let audit = audit_source_spline_cap_precision(source_step, data_section, solid_id)?;
    ensure!(
        audit.fully_certified_strip_faces == audit.strip_faces,
        "{} / {} source strip face trims are precision-proven",
        audit.fully_certified_strip_faces,
        audit.strip_faces
    );
    let mut exchange = monstertruck_io::step::load::ParsedStepSource::parse(source_step)?.exchange;
    crate::step_io::require_supported_sections(&exchange)?;
    let source = &exchange
        .data
        .get(data_section)
        .context("missing DATA section")?
        .entities;
    let index = build_index(source);
    ensure!(
        index.len() == source.len(),
        "source STEP has repeated entity IDs"
    );
    let original_boundary = analytic_boundary_diagnostic(solid_id, source, &index);
    let scale = crate::source_edge_locus::source_solid_length_scale_mm(solid_id, source)?;
    let candidates = identify_rechart_cap_candidates(&audit, source, &index, scale)?;
    let (patches, generated) =
        construct_private_cap_curves(source, &index, &candidates, scale, solid_id)?;
    let entities = &mut exchange.data[data_section].entities;
    apply_private_cap_patches(entities, &index, &patches, generated, scale)?;
    let proof_index = build_index(entities);
    let cap_boundary = analytic_boundary_diagnostic(solid_id, entities, &proof_index);
    let reconciliation = if cap_boundary.supported {
        None
    } else {
        crate::boundary_proof::reconcile_strip_face_senses(solid_id, entities, &proof_index)
    };
    let (reconciled_face_senses, proof_boundary) = if let Some(reconciled) = reconciliation {
        let repairs = reconciled.repairs;
        let diagnostic = reconciled.diagnostic;
        *entities = reconciled.entities;
        (repairs, diagnostic)
    } else {
        (Vec::new(), cap_boundary.clone())
    };
    Ok(SourceCapRechartProposal {
        exchange,
        manifest: SourceCapRechartManifest {
            solid_id,
            data_section,
            author_uncertainty_mm: audit.authored_uncertainty_mm,
            source_boundary_supported: original_boundary.supported,
            cap_recharted_boundary_supported: cap_boundary.supported,
            proof_view_boundary_supported: proof_boundary.supported,
            proof_view_boundary_reason: proof_boundary.reason,
            reconciled_face_senses,
            patches,
        },
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_graph::{build_index, entity_id};
    use monstertruck_io::step::load::ParsedStepSource;

    const ORIGINAL: &[u8] =
        include_bytes!("../../validation/fixtures/native_drafted_round_network_0b5a3544.step");

    #[test]
    fn original_manufacturer_private_proposal_preserves_every_source_entity() -> Result<()> {
        let source = ParsedStepSource::parse(ORIGINAL)?.exchange;
        let original = &source.data[0].entities;
        let proposal = propose_source_cap_rechart(ORIGINAL, 0, 3877)?;
        let updated = &proposal.exchange.data[0].entities;
        let manifest = &proposal.manifest;
        eprintln!(
            "proposed private source cap replacements: {}, original proof={} rechart proof={} reason={:?}",
            manifest.patches.len(),
            manifest.source_boundary_supported,
            manifest.proof_view_boundary_supported,
            manifest.proof_view_boundary_reason
        );
        assert!(manifest.patches.len() >= 10);
        assert!(manifest.proof_view_boundary_supported);
        assert!(!manifest.source_boundary_supported);
        assert!(!manifest.cap_recharted_boundary_supported);
        assert_eq!(manifest.reconciled_face_senses.len(), 10);
        assert_eq!(updated.len(), original.len() + manifest.patches.len());
        let original_index = build_index(original);
        let updated_index = build_index(updated);
        let by_edge = manifest
            .patches
            .iter()
            .map(|p| (p.edge_id, p))
            .collect::<BTreeMap<_, _>>();
        for item in original {
            let id = entity_id(item);
            let changed = &updated[*updated_index.get(&id).context("lost source entity")?];
            if let Some(plan) = by_edge.get(&id) {
                let a = named_entity_parameters(id, "EDGE_CURVE", original, &original_index)
                    .context("original EDGE_CURVE missing")?;
                let b = named_entity_parameters(id, "EDGE_CURVE", updated, &updated_index)
                    .context("rechart EDGE_CURVE missing")?;
                for i in [0, 1, 2, 4] {
                    assert_eq!(a[i], b[i], "source EDGE_CURVE metadata changed");
                }
                assert_eq!(entity_ref_value(&a[3]), Some(plan.original_curve_id));
                assert_eq!(entity_ref_value(&b[3]), Some(plan.private_curve_id));
            } else if let Some(sense) = manifest
                .reconciled_face_senses
                .iter()
                .find(|sense| sense.face_id == id)
            {
                let a = crate::step_graph::simple_record(item).context("original face missing")?;
                let b =
                    crate::step_graph::simple_record(changed).context("recharted face missing")?;
                assert_eq!(a.name, b.name);
                let (Parameter::List(a), Parameter::List(b)) = (&a.parameter, &b.parameter) else {
                    bail!("face parameters malformed");
                };
                assert_eq!(a.len(), 4);
                assert_eq!(b.len(), 4);
                assert_eq!(&a[..3], &b[..3], "non-sense face fields must stay exact");
                assert_eq!(
                    a[3],
                    Parameter::Enumeration(
                        if sense.original_same_sense { "T" } else { "F" }.to_string()
                    )
                );
                assert_eq!(
                    b[3],
                    Parameter::Enumeration(
                        if sense.reconciled_same_sense {
                            "T"
                        } else {
                            "F"
                        }
                        .to_string()
                    )
                );
            } else {
                assert_eq!(item, changed, "original STEP entity #{id} was mutated");
            }
        }
        let step = crate::step_io::write_exchange(&proposal.exchange)?;
        let parsed = ParsedStepSource::parse(step.as_bytes())?.exchange;
        assert_eq!(parsed.data[0].entities.len(), updated.len());
        let again = propose_source_cap_rechart(step.as_bytes(), 0, 3877)?;
        assert!(
            again.manifest.patches.is_empty(),
            "proof-only patch must be idempotent"
        );
        Ok(())
    }

    #[test]
    fn cannot_borrow_missing_author_precision_from_exported_source() -> Result<()> {
        let source = std::str::from_utf8(ORIGINAL)?;
        let missing = source.replacen(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT",
            "NO_UNCERTAINTY_CONTEXT",
            1,
        );
        assert!(propose_source_cap_rechart(missing.as_bytes(), 0, 3877).is_err());
        Ok(())
    }
}
