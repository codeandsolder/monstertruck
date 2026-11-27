//! Explicit source-only chart normalization for independent STEP import audits.
//!
//! A sphere is invariant under a right-handed orthonormal change of its chart.
//! Only directly face-referenced spheres with authored spatial circle trims
//! are eligible. Never rewrite pcurves, other users of a surface, or any
//! occupied geometry to make an external comparator pass.

use super::{DIRECTION_Q, compare_analytic_boundaries, strict_closed_shell_faces};
use crate::brep::{self, CircleSupport, CurveSupport, SphereSupport, SurfaceSupport};
use crate::math3::{add, canonical_direction, cross, dot, mul, norm, normalize, sub};
use crate::step_entities::{entity_ref, named_entity_parameters, push_simple};
use crate::step_graph::{
    ReferenceGraph, build_index, entity_id, entity_ref_value, simple_record_mut,
};
use crate::step_io::{require_supported_sections, write_exchange};
use anyhow::{Context, Result, ensure};
use monstertruck_io::step::load::ParsedStepSource;
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap};
use step_p21::ast::{EntityInstance, Parameter};

/// A sphere chart change, with unchanged original 3D trim and carrier data.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SphereChartChange {
    pub surface_id: u64,
    pub face_ids: Vec<u64>,
    pub original_placement_id: u64,
    pub new_placement_id: u64,
    pub original_axis: [f64; 3],
    pub original_reference: [f64; 3],
    pub new_axis: [f64; 3],
    pub new_reference: [f64; 3],
    pub latitude_edge_ids: Vec<u64>,
    pub meridian_edge_ids: Vec<u64>,
    /// Great-circle-only trimming has no intrinsic latitude/meridian split.
    pub great_circle_edge_ids: Vec<u64>,
}

/// A provenance manifest tied to the exact original and revised bytes.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SphereChartManifest {
    pub solid_id: u64,
    pub original_sha256: String,
    pub recharted_sha256: String,
    pub changes: Vec<SphereChartChange>,
}

/// A separate, fully boundary-verified STEP document for an import audit.
#[derive(Debug, Clone)]
pub struct SphereChartAudit {
    pub step: String,
    pub manifest: SphereChartManifest,
}

struct ChartPlan {
    change: SphereChartChange,
    center_id: u64,
}

fn circles_for_faces(
    faces: &[u64],
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<BTreeMap<u64, CircleSupport>> {
    let mut circles = BTreeMap::new();
    for &face in faces {
        for boundary in brep::face_loops(face, entities, index)? {
            for edge in boundary.edges {
                // A CIRCLE has no surface-chart dependency. In particular,
                // do not unwrap a SURFACE_CURVE or a SEAM_CURVE here.
                named_entity_parameters(edge.curve_id, "CIRCLE", entities, index)?;
                let CurveSupport::Circle(circle) = edge.support else {
                    return None;
                };
                circles.insert(edge.edge_id, circle);
            }
        }
    }
    Some(circles)
}

struct AlignedChart {
    axis: [f64; 3],
    reference: [f64; 3],
    parallels: Vec<u64>,
    meridians: Vec<u64>,
    great_circles: Vec<u64>,
}

fn aligned_chart(
    sphere: SphereSupport,
    circles: &BTreeMap<u64, CircleSupport>,
) -> Option<AlignedChart> {
    let mut parallels = Vec::new();
    let mut meridians = Vec::new();
    for (&edge, &circle) in circles {
        if !super::bounded::circle_on_sphere(circle, sphere.center_mm, sphere.radius_mm) {
            return None;
        }
        let offset = norm(sub(circle.center_mm, sphere.center_mm));
        if offset <= super::ANALYTIC_BOUNDARY_PROOF_TOL_MM {
            meridians.push(edge);
        } else {
            parallels.push(edge);
        }
    }
    // A non-great-circle section anchors the latitude axis. Great circles
    // cannot independently distinguish a latitude from a longitude.
    let axis = canonical_direction(circles[parallels.first()?].normal);
    if norm(cross(axis, sphere.axis)) <= DIRECTION_Q
        || parallels
            .iter()
            .any(|edge| norm(cross(axis, circles[edge].normal)) > DIRECTION_Q)
        || meridians
            .iter()
            .any(|edge| dot(axis, circles[edge].normal).abs() > DIRECTION_Q)
    {
        return None;
    }
    let reference = if let Some(first) = meridians.first() {
        let normal = circles[first].normal;
        if dot(sphere.x_direction, axis).abs() <= DIRECTION_Q
            && dot(sphere.x_direction, normal).abs() <= DIRECTION_Q
        {
            normalize(
                sub(sphere.x_direction, mul(axis, dot(sphere.x_direction, axis))),
                DIRECTION_Q,
            )?
        } else {
            canonical_direction(normalize(cross(normal, axis), DIRECTION_Q)?)
        }
    } else {
        // Latitude-only trimming also has a uniquely anchored polar axis.
        // Its azimuth is arbitrary; reuse a known in-plane circle reference.
        circles[parallels.first()?].x_direction
    };
    Some(AlignedChart {
        axis,
        reference,
        parallels,
        meridians,
        great_circles: Vec::new(),
    })
}

/// A chart is not a physical sphere change. For all-great-circle trims,
/// choose an alternate well-conditioned polar axis whose poles are uniformly
/// away from each authored great-circle plane. Pick the azimuth chart so the
/// entire *minor-arc polygon* lies in an open hemisphere far from its seam.
/// The source and recharted solids still undergo full strict boundary proof.
fn all_great_circle_chart(
    sphere: SphereSupport,
    faces: &[u64],
    circles: &BTreeMap<u64, CircleSupport>,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
) -> Option<AlignedChart> {
    let radius = sphere.radius_mm;
    let tolerance = super::ANALYTIC_BOUNDARY_PROOF_TOL_MM;
    if circles.len() < 3 || !radius.is_finite() || radius <= tolerance {
        return None;
    }
    let mut vertices = Vec::new();
    let normals = circles
        .values()
        .map(|circle| circle.normal)
        .collect::<Vec<_>>();
    for &face in faces {
        let [wire] = brep::face_loops(face, entities, index)?.as_slice() else {
            return None;
        };
        if !wire.outer || wire.edges.len() < 3 {
            return None;
        }
        for edge in &wire.edges {
            let circle = circles.get(&edge.edge_id)?;
            if norm(sub(circle.center_mm, sphere.center_mm)) > tolerance
                || (circle.radius_mm - radius).abs() > tolerance
            {
                return None;
            }
            let y = normalize(cross(circle.normal, circle.x_direction), DIRECTION_Q)?;
            let endpoints = [edge.start_mm, edge.end_mm];
            let mut angles = [0.0; 2];
            for (i, point) in endpoints.into_iter().enumerate() {
                let relative = sub(point, sphere.center_mm);
                if (norm(relative) - radius).abs() > tolerance
                    || dot(relative, circle.normal).abs() > tolerance
                {
                    return None;
                }
                let radial = normalize(relative, DIRECTION_Q)?;
                vertices.push(radial);
                angles[i] = dot(relative, y).atan2(dot(relative, circle.x_direction));
            }
            let sweep = if edge.parameter_forward {
                (angles[1] - angles[0]).rem_euclid(std::f64::consts::TAU)
            } else {
                (angles[0] - angles[1]).rem_euclid(std::f64::consts::TAU)
            };
            // A long arc can cross the longitude seam despite its endpoints
            // being safely positioned. A complete all-great closed sphere is
            // not a small-polygon chart-normalization candidate.
            if !sweep.is_finite()
                || sweep <= DIRECTION_Q
                || sweep >= std::f64::consts::PI - DIRECTION_Q
            {
                return None;
            }
        }
    }
    let centroid = normalize(vertices.iter().copied().fold([0.; 3], add), DIRECTION_Q)?;
    // Find the maximally independent triple of authored great-circle planes.
    // This is deterministic, O(m^3), and independent of face enumeration order.
    let mut best_det = 0.0_f64;
    let mut triple = None;
    for (i, &a) in normals.iter().enumerate() {
        for (j, &b) in normals.iter().enumerate().skip(i + 1) {
            for (k, &c) in normals.iter().enumerate().skip(j + 1) {
                let det = dot(a, cross(b, c)).abs();
                if det > best_det {
                    best_det = det;
                    triple = Some((i, j, k));
                }
            }
        }
    }
    if best_det <= DIRECTION_Q {
        return None;
    }
    let (i, j, k) = triple?;
    let mut best = None::<(f64, [f64; 3], [f64; 3])>;
    for signs in 0_u8..8 {
        let a = if signs & 1 == 0 {
            normals[i]
        } else {
            mul(normals[i], -1.)
        };
        let b = if signs & 2 == 0 {
            normals[j]
        } else {
            mul(normals[j], -1.)
        };
        let c = if signs & 4 == 0 {
            normals[k]
        } else {
            mul(normals[k], -1.)
        };
        let Some(axis) = normalize(add(add(a, b), c), DIRECTION_Q) else {
            continue;
        };
        if norm(cross(axis, sphere.axis)) <= DIRECTION_Q {
            continue;
        }
        let pole_margin = normals
            .iter()
            .map(|n| dot(*n, axis).abs())
            .fold(f64::INFINITY, f64::min);
        if pole_margin <= 0.01 {
            continue;
        }
        let Some(reference) = normalize(sub(centroid, mul(axis, dot(centroid, axis))), DIRECTION_Q)
        else {
            continue;
        };
        let seam_margin = vertices
            .iter()
            .map(|v| dot(*v, reference))
            .fold(f64::INFINITY, f64::min);
        // Both endpoints of every minor great-circle arc lie in this positive
        // hemisphere; by geodesic convexity, the entire arc stays there.
        if seam_margin <= 0.01 {
            continue;
        }
        let score = pole_margin.min(seam_margin);
        if best.as_ref().is_none_or(|(prior, _, _)| score > *prior) {
            best = Some((score, axis, reference));
        }
    }
    let (_, axis, reference) = best?;
    Some(AlignedChart {
        axis,
        reference,
        parallels: Vec::new(),
        meridians: Vec::new(),
        great_circles: circles.keys().copied().collect(),
    })
}

fn plan_chart(
    surface_id: u64,
    faces: Vec<u64>,
    entities: &[EntityInstance],
    index: &HashMap<u64, usize>,
    graph: &ReferenceGraph,
) -> Option<ChartPlan> {
    // The only legal incoming reference here is a direct geometric carrier
    // of one of the selected solid's faces. A pcurve or a derived surface
    // referencing the same sphere makes the chart immutable.
    let users = graph.inbound().get(&surface_id)?;
    if users.len() != faces.len() || users.iter().any(|id| !faces.contains(id)) {
        return None;
    }
    let sphere_params = named_entity_parameters(surface_id, "SPHERICAL_SURFACE", entities, index)?;
    let original_placement_id = entity_ref_value(sphere_params.get(1)?)?;
    let placement =
        named_entity_parameters(original_placement_id, "AXIS2_PLACEMENT_3D", entities, index)?;
    let center_id = entity_ref_value(placement.get(1)?)?;
    let SurfaceSupport::Sphere(sphere) = brep::surface_support(surface_id, entities, index) else {
        return None;
    };
    let circles = circles_for_faces(&faces, entities, index)?;
    let chart = aligned_chart(sphere, &circles)
        .or_else(|| all_great_circle_chart(sphere, &faces, &circles, entities, index))?;
    Some(ChartPlan {
        center_id,
        change: SphereChartChange {
            surface_id,
            face_ids: faces,
            original_placement_id,
            new_placement_id: 0,
            original_axis: sphere.axis,
            original_reference: sphere.x_direction,
            new_axis: chart.axis,
            new_reference: chart.reference,
            latitude_edge_ids: chart.parallels,
            meridian_edge_ids: chart.meridians,
            great_circle_edge_ids: chart.great_circles,
        },
    })
}

fn push_direction(entities: &mut Vec<EntityInstance>, next: &mut u64, direction: [f64; 3]) -> u64 {
    push_simple(
        entities,
        next,
        "DIRECTION",
        vec![
            Parameter::String(String::new()),
            Parameter::List(direction.into_iter().map(Parameter::Real).collect()),
        ],
    )
}

fn apply_chart(
    entities: &mut Vec<EntityInstance>,
    index: &HashMap<u64, usize>,
    next: &mut u64,
    plan: ChartPlan,
) -> Result<SphereChartChange> {
    let ChartPlan {
        mut change,
        center_id,
    } = plan;
    let axis = push_direction(entities, next, change.new_axis);
    let reference = push_direction(entities, next, change.new_reference);
    let placement = push_simple(
        entities,
        next,
        "AXIS2_PLACEMENT_3D",
        vec![
            Parameter::String(String::new()),
            entity_ref(center_id),
            entity_ref(axis),
            entity_ref(reference),
        ],
    );
    let record = simple_record_mut(&mut entities[index[&change.surface_id]])
        .context("sphere record changed")?;
    let Parameter::List(params) = &mut record.parameter else {
        anyhow::bail!("sphere parameters changed");
    };
    ensure!(
        entity_ref_value(&params[1]) == Some(change.original_placement_id),
        "sphere placement changed"
    );
    params[1] = entity_ref(placement);
    change.new_placement_id = placement;
    Ok(change)
}

/// Build a separate chart-normalized source for one closed manifold solid.
///
/// The output retains all original geometry, apart from changing the chart
/// reference of provably eligible spheres. The complete oriented boundary
/// must compare equal after serialization and parsing, using the original
/// strict boundary proof. No external comparator result enters selection.
///
/// # Errors
/// Returns an error for an unsupported STEP exchange, no eligible chart,
/// ID exhaustion, serialization drift, or failed strict boundary proof.
pub fn audit_sphere_charts(input: &[u8], solid_id: u64) -> Result<SphereChartAudit> {
    let original = ParsedStepSource::parse(input)?.exchange;
    require_supported_sections(&original)?;
    ensure!(
        original.data.len() == 1 && original.data[0].meta.is_empty(),
        "sphere chart audit requires one ordinary DATA section"
    );
    let entities = &original.data[0].entities;
    let index = build_index(entities);
    ensure!(index.len() == entities.len(), "duplicate source entity ID");
    let faces =
        strict_closed_shell_faces(solid_id, entities, &index).map_err(anyhow::Error::msg)?;
    let mut by_sphere = BTreeMap::<u64, Vec<u64>>::new();
    for face in faces {
        let surface = brep::face_surface(face, entities, &index).context("missing face carrier")?;
        if named_entity_parameters(surface, "SPHERICAL_SURFACE", entities, &index).is_some() {
            by_sphere.entry(surface).or_default().push(face);
        }
    }
    let graph = ReferenceGraph::new(entities);
    let plans = by_sphere
        .into_iter()
        .filter_map(|(surface, faces)| plan_chart(surface, faces, entities, &index, &graph))
        .collect::<Vec<_>>();
    ensure!(
        !plans.is_empty(),
        "no independently eligible misaligned spherical charts"
    );
    let mut next = entities
        .iter()
        .map(entity_id)
        .max()
        .unwrap_or_default()
        .checked_add(1)
        .context("source entity ID exhaustion")?;
    let added = u64::try_from(plans.len())?
        .checked_mul(3)
        .context("chart entity count overflow")?;
    ensure!(
        next.checked_add(added).is_some(),
        "source entity ID exhaustion"
    );
    let mut revised = original.clone();
    let changes = plans
        .into_iter()
        .map(|plan| apply_chart(&mut revised.data[0].entities, &index, &mut next, plan))
        .collect::<Result<Vec<_>>>()?;
    let step = write_exchange(&revised)?;
    let reparsed = ParsedStepSource::parse(step.as_bytes())?.exchange;
    ensure!(
        reparsed == revised,
        "recharted STEP serialization changed the entity graph"
    );
    let proof =
        compare_analytic_boundaries(solid_id, entities, solid_id, &reparsed.data[0].entities);
    ensure!(
        proof.equivalent,
        "recharted source failed strict boundary equivalence: {:?} / {:?}",
        proof.first.reason,
        proof.second.reason
    );
    let manifest = SphereChartManifest {
        solid_id,
        original_sha256: hex::encode(Sha256::digest(input)),
        recharted_sha256: hex::encode(Sha256::digest(step.as_bytes())),
        changes,
    };
    Ok(SphereChartAudit { step, manifest })
}

#[cfg(test)]
mod tests {
    use super::*;
    use step_p21::ast::Exchange;

    const SOURCE: &[u8] =
        include_bytes!("../../validation/fixtures/native_relief_blind_636f2bd5.step");

    fn sphere_ids(entities: &[EntityInstance]) -> Vec<u64> {
        entities
            .iter()
            .filter_map(|entity| match entity {
                EntityInstance::Simple { id, record } if record.name == "SPHERICAL_SURFACE" => {
                    Some(*id)
                }
                _ => None,
            })
            .collect()
    }

    fn restored_source(
        mut recharted: Exchange,
        original_count: usize,
        changes: &[SphereChartChange],
    ) -> Result<Exchange> {
        let entities = &mut recharted.data[0].entities;
        let index = build_index(entities);
        for change in changes {
            let sphere = simple_record_mut(&mut entities[index[&change.surface_id]])
                .context("recharted sphere")?;
            let Parameter::List(params) = &mut sphere.parameter else {
                anyhow::bail!("sphere parameters");
            };
            assert_eq!(params[1], entity_ref(change.new_placement_id));
            params[1] = entity_ref(change.original_placement_id);
        }
        entities.truncate(original_count);
        Ok(recharted)
    }

    #[test]
    fn source_sphere_chart_audit_preserves_every_original_entity_except_chart_refs() -> Result<()> {
        let original = ParsedStepSource::parse(SOURCE)?.exchange;
        let original_entities = &original.data[0].entities;
        let audit = audit_sphere_charts(SOURCE, 1673)?;
        let revised = ParsedStepSource::parse(audit.step.as_bytes())?.exchange;
        assert_eq!(
            audit.manifest.changes.len(),
            sphere_ids(original_entities).len()
        );
        assert_eq!(audit.manifest.changes.len(), 8);
        assert_eq!(revised.data[0].entities.len(), original_entities.len() + 24);
        assert_eq!(
            audit.manifest.original_sha256,
            hex::encode(Sha256::digest(SOURCE))
        );
        assert_eq!(
            audit.manifest.recharted_sha256,
            hex::encode(Sha256::digest(audit.step.as_bytes()))
        );
        assert_eq!(
            restored_source(revised, original_entities.len(), &audit.manifest.changes)?,
            original
        );
        assert!(audit_sphere_charts(audit.step.as_bytes(), 1673).is_err());
        Ok(())
    }

    #[test]
    fn source_sphere_chart_audit_respects_other_chart_users() -> Result<()> {
        let mut source = ParsedStepSource::parse(SOURCE)?.exchange;
        let entities = &mut source.data[0].entities;
        let mut next = entities
            .iter()
            .map(entity_id)
            .max()
            .context("source entity IDs")?
            + 1;
        let sphere = sphere_ids(entities)[0];
        // An incoming pcurve reference is already sufficient to disallow a
        // chart change, even if the referenced curve representation is absent.
        push_simple(
            entities,
            &mut next,
            "PCURVE",
            vec![
                Parameter::String(String::new()),
                entity_ref(sphere),
                Parameter::NotProvided,
            ],
        );
        let bytes = write_exchange(&source)?;
        let audit = audit_sphere_charts(bytes.as_bytes(), 1673)?;
        assert_eq!(audit.manifest.changes.len(), 7);
        assert!(
            audit
                .manifest
                .changes
                .iter()
                .all(|change| change.surface_id != sphere)
        );
        let revised = ParsedStepSource::parse(audit.step.as_bytes())?.exchange;
        assert_eq!(
            restored_source(
                revised,
                source.data[0].entities.len(),
                &audit.manifest.changes
            )?,
            source
        );
        Ok(())
    }

    #[test]
    fn source_sphere_chart_audit_rejects_invalid_closed_shell() -> Result<()> {
        let mut source = ParsedStepSource::parse(SOURCE)?.exchange;
        let entities = &mut source.data[0].entities;
        let index = build_index(entities);
        let solid = named_entity_parameters(1673, "MANIFOLD_SOLID_BREP", entities, &index)
            .context("source solid")?;
        let shell_id = entity_ref_value(&solid[1]).context("source shell")?;
        let shell = simple_record_mut(&mut entities[index[&shell_id]]).context("shell record")?;
        let Parameter::List(params) = &mut shell.parameter else {
            anyhow::bail!("shell params");
        };
        let Parameter::List(faces) = &mut params[1] else {
            anyhow::bail!("shell faces");
        };
        faces.pop().context("nonempty shell")?;
        assert!(audit_sphere_charts(write_exchange(&source)?.as_bytes(), 1673).is_err());
        Ok(())
    }
}
