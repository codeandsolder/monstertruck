//! Opt-in construction of an analytically healed source STEP edge.
//!
//! The original STEP edge is not modified. A replacement is a separate
//! provenance-bearing candidate and is *never* source-boundary admitted merely
//! because its two analytic incident surfaces intersect exactly.

use crate::brep::{self, SurfaceSupport};
use crate::source_edge_locus::{StepEdge, complete_bspline_edge, incident_edge_plane};
use crate::step_entities::{enumeration_bool, named_entity_parameters, vertex_point};
use crate::step_graph::{build_index, entity_ref_value};
use anyhow::{Context, Result, anyhow, bail};
use monstertruck_geometry::analytic_cone_deviation::{
    ConeDeviationBounds, bound_nurbs_cone_deviation,
};
use monstertruck_geometry::analytic_plane_cone_hyperbola::{
    AxisParallelPlaneConeHyperbola, OrientedHyperbolaTrim,
};
use monstertruck_geometry::analytic_plane_deviation::{
    PlaneDeviationBounds, bound_nurbs_plane_deviation,
};
use monstertruck_geometry::prelude::{
    Matrix4, Point3, Processor, TrimmedCurve, UnitHyperbola, Vector3,
};
use serde::Serialize;

/// Why a replacement is *not* equivalent to an unmodified manufacturer edge.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct HyperbolaHealingEvidence {
    pub original_solid_id: u64,
    pub original_edge_id: u64,
    pub original_basis_curve_id: u64,
    pub source_start_vertex_id: u64,
    pub source_end_vertex_id: u64,
    pub source_curve_same_sense: bool,
    pub source_plane_face_id: u64,
    pub source_cone_face_id: u64,
    pub source_coordinate_scale_mm: f64,
    /// The source STEP representation context is a separate quantity; this
    /// geometric witness does not authorize any increased model tolerance.
    pub source_representation_tolerance_mm: f64,
    pub source_plane_max_deviation_lower_mm: f64,
    pub source_plane_max_deviation_upper_mm: f64,
    pub source_cone_max_deviation_lower_mm: f64,
    pub source_cone_max_deviation_upper_mm: f64,
    /// Strictly geometric lower witness to directed source-to-healed locus
    /// separation: the hyperbola is on its exact analytic cone, while the
    /// source-plane witness subtracts the healed finite trim's global
    /// residual to the original (slightly nonparallel) source plane.
    /// It does not certify an upper Hausdorff distance or matching PCURVEs.
    pub minimum_source_to_healed_locus_separation_mm: f64,
    pub healed_source_start_residual_mm: f64,
    pub healed_source_end_residual_mm: f64,
    pub cone_axis_plane_normal_dot: f64,
    /// Whole-healed-trim upper bound against the original authored plane,
    /// including the accepted source-axis angular conditioning residual.
    pub healed_source_plane_max_residual_mm: f64,
    /// Permanent explicit distinction between an original source edge and a
    /// newly constructed analytic edge requiring topology/PCURVE checks.
    pub status: &'static str,
}

/// Verified source STEP solid and carrier identities with explicit budgets.
///
/// Source uncertainty is resolved from the referenced
/// `GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT`; it cannot be overridden.
#[derive(Clone, Copy)]
pub struct HyperbolaHealingRequest<'a> {
    pub solid_id: u64,
    pub source: StepEdge<'a>,
    pub source_plane_face_id: u64,
    pub source_cone_face_id: u64,
    /// The source file coordinate-to-millimetre scale, resolved by the caller.
    pub source_coordinate_scale_mm: f64,
    pub endpoint_tolerance_mm: f64,
    pub target_gap_mm: f64,
}

/// A constructive geometry proposal with explicit original-edge identity.
///
/// Its source basis curve must be retained by the caller alongside this
/// new edge. This is read-only and makes no changes to input STEP entities.
#[derive(Clone, Debug)]
pub struct HyperbolaHealingCandidate {
    pub evidence: HyperbolaHealingEvidence,
    pub analytic: AxisParallelPlaneConeHyperbola,
    pub trim: OrientedHyperbolaTrim,
}

impl HyperbolaHealingCandidate {
    #[must_use]
    pub fn modeling_curve(&self) -> Processor<TrimmedCurve<UnitHyperbola<Point3>>, Matrix4> {
        self.analytic.modeling_curve(self.trim)
    }
}

fn point3(point: [f64; 3], scale: f64) -> Point3 {
    Point3::new(point[0] * scale, point[1] * scale, point[2] * scale)
}

#[derive(Clone, Copy)]
struct SourceEdgeTopology {
    start_vertex_id: u64,
    end_vertex_id: u64,
    basis_curve_id: u64,
    curve_same_sense: bool,
    start: [f64; 3],
    end: [f64; 3],
}

fn source_edge_topology(
    source: StepEdge<'_>,
    index: &std::collections::HashMap<u64, usize>,
) -> Result<SourceEdgeTopology> {
    let params = named_entity_parameters(source.id, "EDGE_CURVE", source.entities, index)
        .with_context(|| format!("missing source edge #{}", source.id))?;
    let source_start_vertex_id = params
        .get(1)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow!("source edge has no starting vertex"))?;
    let source_end_vertex_id = params
        .get(2)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow!("source edge has no ending vertex"))?;
    let original_basis_curve_id = params
        .get(3)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow!("source edge has no basis curve"))?;
    let source_curve_same_sense = params
        .get(4)
        .and_then(enumeration_bool)
        .ok_or_else(|| anyhow!("source edge has no valid curve sense"))?;
    let start = vertex_point(source_start_vertex_id, source.entities, index)
        .ok_or_else(|| anyhow!("source edge has no starting vertex coordinates"))?;
    let end = vertex_point(source_end_vertex_id, source.entities, index)
        .ok_or_else(|| anyhow!("source edge has no ending vertex coordinates"))?;

    Ok(SourceEdgeTopology {
        start_vertex_id: source_start_vertex_id,
        end_vertex_id: source_end_vertex_id,
        basis_curve_id: original_basis_curve_id,
        curve_same_sense: source_curve_same_sense,
        start,
        end,
    })
}

fn source_analytic_carriers(
    source: StepEdge<'_>,
    plane_face: u64,
    cone_face: u64,
    index: &std::collections::HashMap<u64, usize>,
) -> Result<(crate::source_edge_locus::ReferencePlane, brep::ConeSupport)> {
    let plane = incident_edge_plane(plane_face, source.id, source.entities)?;
    let plane_surface_id = brep::face_surface(plane_face, source.entities, index)
        .ok_or_else(|| anyhow!("source face has no plane surface"))?;
    if named_entity_parameters(plane_surface_id, "PLANE", source.entities, index).is_none() {
        bail!("source face must have a literal authored STEP PLANE");
    }

    let cone_loops = brep::face_loops(cone_face, source.entities, index)
        .ok_or_else(|| anyhow!("source conical face has no usable edge loops"))?;
    if !cone_loops
        .iter()
        .flat_map(|loop_| &loop_.edges)
        .any(|edge| edge.edge_id == source.id)
    {
        bail!("source edge is not incident on the specified conical face");
    }
    let cone_surface_id = brep::face_surface(cone_face, source.entities, index)
        .ok_or_else(|| anyhow!("source conical face has no carrier"))?;
    if named_entity_parameters(cone_surface_id, "CONICAL_SURFACE", source.entities, index).is_none()
    {
        bail!("source conical face must have a literal authored STEP CONICAL_SURFACE");
    }
    let SurfaceSupport::Cone(cone) = brep::surface_support(cone_surface_id, source.entities, index)
    else {
        bail!("source conical carrier is numerically unsupported");
    };
    Ok((plane, cone))
}

fn verified_solid_uncertainty_mm(
    request: HyperbolaHealingRequest<'_>,
    index: &std::collections::HashMap<u64, usize>,
) -> Result<f64> {
    let faces = brep::solid_face_ids(request.solid_id, request.source.entities, index)
        .ok_or_else(|| anyhow!("source solid has no valid closed shell topology"))?;
    if !faces.contains(&request.source_plane_face_id)
        || !faces.contains(&request.source_cone_face_id)
    {
        bail!("both incident surfaces must be faces of the given source solid");
    }
    crate::solid_revolutions::source_tolerance_by_representation_item(
        request.source.entities,
        index,
    )
    .get(&request.solid_id)
    .copied()
    .ok_or_else(|| anyhow!("source solid has no independently declared usable STEP uncertainty"))
}

fn source_carrier_deviations(
    source: StepEdge<'_>,
    plane: crate::source_edge_locus::ReferencePlane,
    cone: brep::ConeSupport,
    scale_mm: f64,
    endpoint_tolerance_mm: f64,
    target_gap_mm: f64,
) -> Result<(PlaneDeviationBounds, ConeDeviationBounds)> {
    let source_curve = complete_bspline_edge(source, endpoint_tolerance_mm / scale_mm)?;
    let raw_target_gap = target_gap_mm / scale_mm;
    let raw_plane_origin = point3(plane.origin_mm, 1.0);
    let raw_normal = Vector3::new(plane.normal[0], plane.normal[1], plane.normal[2]);
    let raw_cone_origin = point3(cone.reference_origin_mm, 1.0);
    let raw_axis = Vector3::new(cone.axis[0], cone.axis[1], cone.axis[2]);
    let plane =
        bound_nurbs_plane_deviation(&source_curve, raw_plane_origin, raw_normal, raw_target_gap)
            .ok_or_else(|| anyhow!("source curve's complete plane bound unsupported"))?;
    let cone = bound_nurbs_cone_deviation(
        &source_curve,
        raw_cone_origin,
        raw_axis,
        cone.reference_radius_mm,
        cone.semi_angle_rad.tan(),
        raw_target_gap,
    )
    .ok_or_else(|| anyhow!("source curve's complete cone bound unsupported"))?;
    Ok((plane, cone))
}

/// Derive the exact hyperbolic carrier and oriented endpoints for a *new* edge.
///
/// The plane and conical face must both be literally authored in STEP and
/// topologically incident on the same original `EDGE_CURVE`. The source cubic
/// must have a global full-domain monotone-trim certificate (not merely
/// matching endpoints). Source-vs-carrier Bernstein bounds are evaluated
/// across every positive spline knot interval. The maximum of the two
/// independent source-to-carrier separation witnesses must exceed the
/// original representation uncertainty to justify a repair candidate.
///
/// The output stays in millimetres for the native kernel. The original STEP
/// document may use other units, hence the caller must provide its previously
/// resolved coordinate scale explicitly. No unsupported unit is assumed to be
/// millimetres. A successful result is not whole-BRep equivalence or admission.
///
/// # Errors
/// Returns an error for invalid or unverified source geometry, nonincident
/// faces, unsupported conic regimes, malformed units/tolerances, or vertices
/// not lying on the constructed hyperbola to the explicit endpoint budget.
pub fn propose_plane_cone_hyperbola_healing(
    request: HyperbolaHealingRequest<'_>,
) -> Result<HyperbolaHealingCandidate> {
    let HyperbolaHealingRequest {
        solid_id,
        source,
        source_plane_face_id,
        source_cone_face_id,
        source_coordinate_scale_mm,
        endpoint_tolerance_mm,
        target_gap_mm,
    } = request;
    if ![
        source_coordinate_scale_mm,
        endpoint_tolerance_mm,
        target_gap_mm,
    ]
    .into_iter()
    .all(|value| value.is_finite() && value > 0.0)
    {
        bail!("source units and tolerance budgets must all be positive and finite");
    }
    if source_plane_face_id == source_cone_face_id {
        bail!("plane and cone must be distinct original STEP faces");
    }
    let index = build_index(source.entities);
    let source_representation_tolerance_mm = verified_solid_uncertainty_mm(request, &index)?;
    let topology = source_edge_topology(source, &index)?;

    let (plane, cone) =
        source_analytic_carriers(source, source_plane_face_id, source_cone_face_id, &index)?;
    let slope = cone.semi_angle_rad.tan();
    let raw_normal = Vector3::new(plane.normal[0], plane.normal[1], plane.normal[2]);
    let raw_axis = Vector3::new(cone.axis[0], cone.axis[1], cone.axis[2]);
    let (plane_bound, cone_bound) = source_carrier_deviations(
        source,
        plane,
        cone,
        source_coordinate_scale_mm,
        endpoint_tolerance_mm,
        target_gap_mm,
    )?;
    let scale = source_coordinate_scale_mm;
    let analytic = AxisParallelPlaneConeHyperbola::from_carriers(
        point3(plane.origin_mm, scale),
        raw_normal,
        point3(cone.reference_origin_mm, scale),
        raw_axis,
        cone.reference_radius_mm * scale,
        slope,
    )
    .ok_or_else(|| anyhow!("source plane/cone intersection is not a supported hyperbola"))?;
    let trim = analytic
        .oriented_trim(
            point3(topology.start, scale),
            point3(topology.end, scale),
            endpoint_tolerance_mm,
        )
        .ok_or_else(|| anyhow!("original directed vertices do not trim the analytic hyperbola"))?;
    let healed_source_plane_max_residual_mm = analytic
        .source_plane_residual_upper(point3(plane.origin_mm, scale), raw_normal, trim)
        .ok_or_else(|| anyhow!("healed edge plane residual cannot be bounded"))?;
    if healed_source_plane_max_residual_mm > endpoint_tolerance_mm {
        bail!("healed edge does not lie on the original authored plane to the required budget");
    }
    // The constructed hyperbola lies exactly on its analytic cone in the
    // mathematical model, but the input plane's tiny angular misalignment
    // gives the finite healed trim a nonzero maximum distance to the
    // *original* plane. Subtract that whole-trim upper bound before applying
    // the reverse triangle inequality to the source-plane lower witness.
    // Otherwise the reported directed source→healed locus distance could be
    // overstated, even when the endpoint discrepancy is extremely small.
    let plane_separation = plane_bound
        .lower
        .mul_add(scale, -healed_source_plane_max_residual_mm)
        .max(0.0);
    let cone_separation = cone_bound.lower * scale;
    let minimum_source_to_healed_locus_separation_mm = plane_separation.max(cone_separation);
    if minimum_source_to_healed_locus_separation_mm <= source_representation_tolerance_mm {
        bail!("source is not proven geometrically distinct from the intended healing");
    }
    Ok(HyperbolaHealingCandidate {
        evidence: HyperbolaHealingEvidence {
            original_solid_id: solid_id,
            original_edge_id: source.id,
            original_basis_curve_id: topology.basis_curve_id,
            source_start_vertex_id: topology.start_vertex_id,
            source_end_vertex_id: topology.end_vertex_id,
            source_curve_same_sense: topology.curve_same_sense,
            source_plane_face_id,
            source_cone_face_id,
            source_coordinate_scale_mm: scale,
            source_representation_tolerance_mm,
            source_plane_max_deviation_lower_mm: plane_bound.lower * scale,
            source_plane_max_deviation_upper_mm: plane_bound.upper * scale,
            source_cone_max_deviation_lower_mm: cone_bound.lower * scale,
            source_cone_max_deviation_upper_mm: cone_bound.upper * scale,
            minimum_source_to_healed_locus_separation_mm,
            healed_source_start_residual_mm: trim.start_residual,
            healed_source_end_residual_mm: trim.end_residual,
            cone_axis_plane_normal_dot: analytic.normal_axis_misalignment,
            healed_source_plane_max_residual_mm,
            status: "new_analytic_edge_requires_explicit_healing_and_brep_proof",
        },
        analytic,
        trim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::{InnerSpace, ParametricCurve};
    use monstertruck_io::step::save::StepDisplay;

    #[test]
    fn original_manufacturer_slot_edges_generate_new_oriented_hyperbolas() -> Result<()> {
        let source = std::str::from_utf8(include_bytes!(
            "../validation/fixtures/production_radial_slot.step"
        ))?;
        let exchange = step_p21::parser::parse(source)?;
        let [section] = exchange.data.as_slice() else {
            bail!("fixture must have exactly one STEP DATA section");
        };
        for (edge_id, plane_face, cone_face, minimum) in [
            (93481, 125603, 58899, 1.57e-4),
            (120654, 75450, 36791, 8.38e-5),
        ] {
            let candidate = propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
                solid_id: 15385,
                source: StepEdge {
                    id: edge_id,
                    entities: &section.entities,
                },
                source_plane_face_id: plane_face,
                source_cone_face_id: cone_face,
                source_coordinate_scale_mm: 1.0,
                endpoint_tolerance_mm: 1.0e-8,
                target_gap_mm: 1.0e-9,
            })?;
            let evidence = &candidate.evidence;
            assert_eq!(evidence.original_solid_id, 15385);
            assert_eq!(evidence.original_edge_id, edge_id);
            assert_eq!(evidence.source_representation_tolerance_mm, 1.0e-5);
            assert!(evidence.minimum_source_to_healed_locus_separation_mm > minimum);
            assert!(!candidate.trim.forward());
            assert!(evidence.healed_source_start_residual_mm < 1.0e-11);
            assert!(evidence.healed_source_end_residual_mm < 1.0e-11);
            assert!(evidence.healed_source_plane_max_residual_mm < 1.0e-8);
            let original_plane_separation = (evidence.source_plane_max_deviation_lower_mm
                - evidence.healed_source_plane_max_residual_mm)
                .max(0.0);
            let certified_separation =
                original_plane_separation.max(evidence.source_cone_max_deviation_lower_mm);
            assert!(
                (evidence.minimum_source_to_healed_locus_separation_mm - certified_separation)
                    .abs()
                    < 1.0e-16,
                "directed separation must account for original-plane residual on the healed trim"
            );
            // A new analytic curve, retaining original topological direction.
            let curve = candidate.modeling_curve();
            let first = curve.evaluate(candidate.trim.start_parameter);
            let last = curve.evaluate(candidate.trim.end_parameter);
            let actual_first = candidate
                .analytic
                .evaluate(candidate.trim.start_parameter)
                .ok_or_else(|| anyhow!("invalid analytic starting vertex"))?;
            let actual_last = candidate
                .analytic
                .evaluate(candidate.trim.end_parameter)
                .ok_or_else(|| anyhow!("invalid analytic ending vertex"))?;
            assert!((first - actual_first).magnitude() < 1.0e-10);
            assert!((last - actual_last).magnitude() < 1.0e-10);
            let serialized = StepDisplay::new(&curve, 900).to_string();
            assert!(serialized.contains("HYPERBOLA("));
            assert!(serialized.contains("AXIS2_PLACEMENT_3D("));
            assert!(!serialized.contains("B_SPLINE_CURVE"));
        }
        Ok(())
    }

    #[test]
    fn unrelated_faces_wrong_solid_and_tampered_uncertainty_are_rejected() -> Result<()> {
        let source = std::str::from_utf8(include_bytes!(
            "../validation/fixtures/production_radial_slot.step"
        ))?;
        let exchange = step_p21::parser::parse(source)?;
        let request = HyperbolaHealingRequest {
            solid_id: 15385,
            source: StepEdge {
                id: 93481,
                entities: &exchange.data[0].entities,
            },
            source_plane_face_id: 125603,
            source_cone_face_id: 58899,
            source_coordinate_scale_mm: 1.0,
            endpoint_tolerance_mm: 1.0e-8,
            target_gap_mm: 1.0e-9,
        };
        assert!(
            propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
                source_plane_face_id: 75450,
                ..request
            })
            .is_err()
        );
        assert!(
            propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
                solid_id: 999_999,
                ..request
            })
            .is_err()
        );
        assert!(
            propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
                source_coordinate_scale_mm: 0.0,
                ..request
            })
            .is_err()
        );
        let tampered = source.replace("LENGTH_MEASURE(0.00001)", "LENGTH_MEASURE(0.001)");
        assert_ne!(source, tampered);
        let untrusted = step_p21::parser::parse(&tampered)?;
        assert!(
            propose_plane_cone_hyperbola_healing(HyperbolaHealingRequest {
                source: StepEdge {
                    id: 93481,
                    entities: &untrusted.data[0].entities
                },
                ..request
            })
            .is_err()
        );
        Ok(())
    }
}
