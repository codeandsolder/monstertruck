//! Construct an explicitly separate regular intersection of two authored
//! cylindrical carriers incident on the same defective source edge.

use super::{
    SourceEdgeTopology, SourceFaceEdgeUse, point3, source_edge_topology,
    unique_original_face_edge_use,
};
use crate::brep::{self, CurveSupport, SurfaceSupport};
use crate::source_edge_locus::{StepEdge, source_solid_length_scale_mm};
use crate::step_entities::named_entity_parameters;
use crate::step_graph::build_index;
use anyhow::{Context, Result, bail};
use monstertruck_geometry::analytic_cylinder_intersection::{
    BoundedSourceIntersection, CylinderCarrier, CylinderIntersectionArc,
};
use monstertruck_geometry::prelude::{InnerProduct, NurbsCurve, Point3, Vector3, Vector4};
use serde::Serialize;

#[derive(Clone, Copy)]
pub struct CylinderHealingRequest<'a> {
    pub solid_id: u64,
    pub source: StepEdge<'a>,
    pub endpoint_tolerance_mm: f64,
    pub target_gap_mm: f64,
}

/// The original edge remains unchanged. The new curve is an analytical
/// intersection, not an original-source boundary equivalence certificate.
#[derive(Clone, Debug, Serialize)]
pub struct CylinderHealingEvidence {
    pub original_solid_id: u64,
    pub original_edge_id: u64,
    pub original_basis_curve_id: u64,
    pub source_start_vertex_id: u64,
    pub source_end_vertex_id: u64,
    pub source_curve_same_sense: bool,
    pub primary_face_id: u64,
    pub secondary_face_id: u64,
    pub primary_coedge: SourceFaceEdgeUse,
    pub secondary_coedge: SourceFaceEdgeUse,
    pub source_coordinate_scale_mm: f64,
    pub source_representation_tolerance_mm: f64,
    pub primary_deviation_mm: [f64; 2],
    pub secondary_deviation_mm: [f64; 2],
    pub locus_separation_mm: [f64; 2],
    pub original_endpoint_residual_mm: [f64; 2],
    pub healed_endpoint_residual_mm: [f64; 2],
    pub angular_sweep_rad: f64,
    pub minimum_transverse_root_mm: f64,
    pub status: &'static str,
}

#[derive(Clone, Debug)]
pub struct CylinderHealingCandidate {
    pub evidence: CylinderHealingEvidence,
    /// Angular parameter zero to sweep follows the original basis curve,
    /// independently of the `EDGE_CURVE` same-sense flag in the evidence.
    pub arc: CylinderIntersectionArc,
}

#[derive(Clone, Copy)]
struct IncidentCylinder {
    face_id: u64,
    coedge: SourceFaceEdgeUse,
    carrier: CylinderCarrier,
}

fn incident_cylinders(
    request: CylinderHealingRequest<'_>,
    index: &std::collections::HashMap<u64, usize>,
    scale: f64,
) -> Result<[IncidentCylinder; 2]> {
    let source = request.source;
    let faces = brep::solid_face_ids(request.solid_id, source.entities, index)
        .context("source solid does not resolve to a manifold face set")?;
    let mut incident = Vec::new();
    for face_id in faces {
        let loops = brep::face_loops(face_id, source.entities, index)
            .context("cannot establish complete edge incidence in source solid")?;
        if !loops
            .iter()
            .flat_map(|boundary| &boundary.edges)
            .any(|edge| edge.edge_id == source.id)
        {
            continue;
        }
        let coedge = unique_original_face_edge_use(face_id, source.id, source.entities, index)?;
        let surface_id = brep::face_surface(face_id, source.entities, index)
            .context("incident face does not resolve to a surface")?;
        if named_entity_parameters(surface_id, "CYLINDRICAL_SURFACE", source.entities, index)
            .is_none()
        {
            bail!("incident source face is not an authored CYLINDRICAL_SURFACE");
        }
        let SurfaceSupport::Cylinder(cylinder) =
            brep::surface_support(surface_id, source.entities, index)
        else {
            bail!("unsupported authored cylindrical carrier");
        };
        let carrier = CylinderCarrier::new(
            point3(cylinder.axis_origin_mm, scale),
            Vector3::from(cylinder.axis),
            cylinder.radius_mm * scale,
        )
        .context("invalid cylindrical carrier")?;
        incident.push(IncidentCylinder {
            face_id,
            coedge,
            carrier,
        });
    }
    let [first, second]: [IncidentCylinder; 2] = incident.try_into().map_err(|_| {
        anyhow::anyhow!("source edge must have exactly two incident cylindrical faces")
    })?;
    if first.face_id == second.face_id
        || first.coedge.parameter_forward == second.coedge.parameter_forward
    {
        bail!("original neighboring faces do not oppose the same edge traversal");
    }
    Ok([first, second])
}

fn directed_vertex_residuals(
    mut endpoints: [Point3; 2],
    topology: SourceEdgeTopology,
    scale: f64,
    tolerance: f64,
) -> Result<[f64; 2]> {
    if !topology.curve_same_sense {
        endpoints.reverse();
    }
    let residual = [
        (endpoints[0] - point3(topology.start, scale)).norm(),
        (endpoints[1] - point3(topology.end, scale)).norm(),
    ];
    if !residual
        .iter()
        .all(|value| value.is_finite() && *value <= tolerance)
        || topology.start_vertex_id == topology.end_vertex_id
        || (point3(topology.start, scale) - point3(topology.end, scale)).norm() <= tolerance
    {
        bail!("curve does not match distinct authored directed vertices");
    }
    Ok(residual)
}

fn bounded_candidate(
    first: IncidentCylinder,
    second: IncidentCylinder,
    curve: &NurbsCurve<Vector4>,
    tolerance: f64,
    target_gap: f64,
) -> Result<(
    BoundedSourceIntersection,
    IncidentCylinder,
    IncidentCylinder,
)> {
    let mut candidates =
        [(first, second), (second, first)]
            .into_iter()
            .filter_map(|(primary, secondary)| {
                BoundedSourceIntersection::from_spline(
                    primary.carrier,
                    secondary.carrier,
                    curve,
                    tolerance,
                    target_gap,
                )
                .map(|bound| (bound, primary, secondary))
            });
    candidates
        .next()
        .map(|first| {
            candidates.fold(first, |best, candidate| {
                if candidate.0.separation_upper < best.0.separation_upper {
                    candidate
                } else {
                    best
                }
            })
        })
        .context("no supported source-guided regular cylinder intersection")
}

/// Propose a new regular analytical edge only for a proven source defect.
///
/// Resolve source units, uncertainty and both carriers from the original
/// solid. A full-domain source angular certificate chooses the branch; the
/// two authored vertices must independently match the directed endpoints.
/// The returned bounds describe edge loci, not repaired faces or a solid.
///
/// # Errors
/// Returns an error for unsupported topology/geometry, ambiguous angular
/// coverage, insufficient separation evidence or invalid error budgets.
pub fn propose_cylinder_intersection_healing(
    request: CylinderHealingRequest<'_>,
) -> Result<CylinderHealingCandidate> {
    let source = request.source;
    let index = build_index(source.entities);
    let scale = source_solid_length_scale_mm(request.solid_id, source.entities)?;
    let uncertainty = crate::solid_revolutions::declared_source_uncertainty_by_representation_item(
        source.entities,
        &index,
    )
    .get(&request.solid_id)
    .copied()
    .context("source solid has no authored uncertainty")?;
    let topology = source_edge_topology(source, &index)?;
    let [first, second] = incident_cylinders(request, &index, scale)?;
    let CurveSupport::BSpline(mut spline) =
        brep::curve_support(topology.basis_curve_id, source.entities, &index)
    else {
        bail!("source edge has no supported B-spline basis");
    };
    for point in &mut spline.control_points_mm {
        for coordinate in point {
            *coordinate *= scale;
        }
    }
    let curve = spline.rational_curve().context("invalid source spline")?;
    let tolerance = request.endpoint_tolerance_mm;
    let original_endpoint_residual_mm = directed_vertex_residuals(
        [
            point3(
                *spline.control_points_mm.first().context("empty spline")?,
                1.0,
            ),
            point3(
                *spline.control_points_mm.last().context("empty spline")?,
                1.0,
            ),
        ],
        topology,
        scale,
        tolerance,
    )?;

    let (bound, primary, secondary) =
        bounded_candidate(first, second, &curve, tolerance, request.target_gap_mm)?;
    if bound.separation_lower <= uncertainty {
        bail!("source defect does not exceed the authored uncertainty");
    }
    let arc = bound.proposal.arc;
    let residual = directed_vertex_residuals(
        [
            arc.evaluate(0.0).context("invalid analytic start")?,
            arc.evaluate(arc.sweep()).context("invalid analytic end")?,
        ],
        topology,
        scale,
        tolerance,
    )?;
    Ok(CylinderHealingCandidate {
        evidence: CylinderHealingEvidence {
            original_solid_id: request.solid_id,
            original_edge_id: source.id,
            original_basis_curve_id: topology.basis_curve_id,
            source_start_vertex_id: topology.start_vertex_id,
            source_end_vertex_id: topology.end_vertex_id,
            source_curve_same_sense: topology.curve_same_sense,
            primary_face_id: primary.face_id,
            secondary_face_id: secondary.face_id,
            primary_coedge: primary.coedge,
            secondary_coedge: secondary.coedge,
            source_coordinate_scale_mm: scale,
            source_representation_tolerance_mm: uncertainty,
            primary_deviation_mm: [bound.primary_deviation.lower, bound.primary_deviation.upper],
            secondary_deviation_mm: [
                bound.secondary_deviation.lower,
                bound.secondary_deviation.upper,
            ],
            locus_separation_mm: [bound.separation_lower, bound.separation_upper],
            original_endpoint_residual_mm,
            healed_endpoint_residual_mm: residual,
            angular_sweep_rad: arc.sweep(),
            minimum_transverse_root_mm: arc.minimum_root(),
            status: "separate_analytic_edge_requires_face_and_topology_reconstruction",
        },
        arc,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_crossing_cylinder_defects_generate_bounded_separate_edges() -> Result<()> {
        let source = std::str::from_utf8(include_bytes!(
            "../../validation/fixtures/native_mixed_treatment_4a6356.step"
        ))?;
        let exchange = step_p21::parser::parse(source)?;
        let [section] = exchange.data.as_slice() else {
            bail!("fixture must have one DATA section");
        };
        for (edge, other_face) in [(1049, 308), (1063, 898)] {
            let request = CylinderHealingRequest {
                solid_id: 171,
                source: StepEdge {
                    id: edge,
                    entities: &section.entities,
                },
                endpoint_tolerance_mm: 1e-8,
                target_gap_mm: 1e-9,
            };
            let candidate = propose_cylinder_intersection_healing(request)?;
            let evidence = &candidate.evidence;
            eprintln!("{evidence:?}");
            assert_eq!(evidence.primary_face_id, 680);
            assert_eq!(evidence.secondary_face_id, other_face);
            assert!(evidence.locus_separation_mm[0] > 1e-4);
            assert!(evidence.locus_separation_mm[1] < 1e-3);
            assert!(
                evidence
                    .healed_endpoint_residual_mm
                    .iter()
                    .all(|r| *r < 1e-8)
            );
            assert!(
                propose_cylinder_intersection_healing(CylinderHealingRequest {
                    solid_id: 1,
                    ..request
                })
                .is_err()
            );
        }
        Ok(())
    }
}
