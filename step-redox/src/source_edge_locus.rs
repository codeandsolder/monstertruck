//! Whole-domain STEP edge-locus comparison via a separating analytic plane.
//!
//! The input `EDGE_CURVE` must cover its entire referenced B-spline domain.
//! STEP permits topological trims of longer basis curves: those are rejected
//! here rather than treating an off-trim excursion as a real edge witness.
//! The geometric comparison is independent of knot/parameter correspondence.

use crate::brep::{self, CurveSupport};
use crate::step_entities::{named_entity_parameters, vertex_point};
use crate::step_graph::{build_index, entity_ref_value};
use anyhow::{Context, Result, anyhow, bail};
use monstertruck_geometry::analytic_plane_deviation::{
    PlaneDirectedSeparationBounds, bound_nurbs_curves_separation_by_plane,
};
use monstertruck_geometry::prelude::{
    BoundedCurve, BsplineCurve, KnotVector, NurbsCurve, ParametricCurve, Point3, Vector3, Vector4,
};
use step_p21::ast::EntityInstance;

/// Raw STEP edge identity in its own DATA section.
#[derive(Clone, Copy)]
pub struct StepEdge<'a> {
    pub id: u64,
    pub entities: &'a [EntityInstance],
}

/// One independently specified Euclidean plane.
#[derive(Clone, Copy)]
pub struct ReferencePlane {
    pub origin_mm: [f64; 3],
    pub normal: [f64; 3],
}

fn point_distance(point: Point3, vertex: [f64; 3]) -> f64 {
    (point.x - vertex[0])
        .hypot(point.y - vertex[1])
        .hypot(point.z - vertex[2])
}

/// Sufficient global proof that the edge cannot revisit an endpoint.
///
/// A nondecreasing control polygon projection gives a nondecreasing rational
/// B-spline projection by the positive-weight normalized B-spline basis'
/// total positivity. Requiring a strict, arithmetic-separated increase for
/// every control and continuity at internal knots also rules out loops and
/// ambiguous repeated visits to either endpoint. Nonmonotone curves may be
/// perfectly valid; this conservative *proof* simply refuses their trims.
fn has_strict_chord_progress(
    control_points: &[[f64; 3]],
    weights: &[f64],
    degree: usize,
    knots: &[f64],
    first: Point3,
    last: Point3,
) -> bool {
    let count = control_points.len();
    let Some(order) = degree.checked_add(1) else {
        return false;
    };
    if count < order || weights.len() != count || knots.len() != count + order {
        return false;
    }
    // Requiring clamped end knots makes both full-domain endpoints identical
    // to the first/last Euclidean controls, rather than merely close to them.
    if knots[..order]
        .iter()
        .any(|knot| knot.total_cmp(&knots[0]).is_ne())
        || knots[knots.len() - order..]
            .iter()
            .any(|knot| knot.total_cmp(&knots[knots.len() - 1]).is_ne())
    {
        return false;
    }
    let chord = [last.x - first.x, last.y - first.y, last.z - first.z];
    let len2 = chord[0].mul_add(chord[0], chord[1].mul_add(chord[1], chord[2] * chord[2]));
    if !len2.is_finite() || len2 <= 0.0 {
        return false;
    }
    // A full internal multiplicity of degree + 1 can conceal a geometric
    // discontinuity. Endpoint matching must not authorize that union of loci.
    let domain_start = knots[degree];
    let domain_end = knots[count];
    if knots.windows(degree + 1).any(|window| {
        window[0].total_cmp(&window[degree]).is_eq()
            && window[0] > domain_start
            && window[0] < domain_end
    }) {
        return false;
    }
    let scale = control_points
        .iter()
        .flatten()
        .chain([&first.x, &first.y, &first.z, &last.x, &last.y, &last.z])
        .fold(1.0_f64, |scale, value| scale.max(value.abs()));
    let weight_min = weights.iter().copied().fold(f64::INFINITY, f64::min);
    let weight_max = weights.iter().copied().fold(0.0_f64, f64::max);
    if !weight_min.is_finite() || weight_min <= 0.0 || !weight_max.is_finite() {
        return false;
    }
    let floor = 256.0 * f64::EPSILON * scale * (weight_max / weight_min);
    if !floor.is_finite() {
        return false;
    }
    let inv_len = len2.sqrt().recip();
    let mut previous = None;
    for control in control_points {
        let relative = [
            control[0] - first.x,
            control[1] - first.y,
            control[2] - first.z,
        ];
        let progress = chord[0].mul_add(
            relative[0],
            chord[1].mul_add(relative[1], chord[2] * relative[2]),
        ) * inv_len;
        if !progress.is_finite() || previous.is_some_and(|prev| progress <= prev + floor) {
            return false;
        }
        previous = Some(progress);
    }
    true
}

/// Validate that the geometric B-spline's complete active parameter range
/// agrees with an `EDGE_CURVE`'s two topological endpoints.
///
/// The STEP B-spline parser validates weights and knots. Endpoint matching
/// in either orientation establishes full-range coverage for an ordinary
/// open edge; unresolved cyclic or self-intersecting topologies need a
/// separate trim/branch proof, never inferred from this metric.
pub(crate) fn complete_bspline_edge(
    step: StepEdge<'_>,
    endpoint_tolerance_mm: f64,
) -> Result<NurbsCurve<Vector4>> {
    let index = build_index(step.entities);
    complete_bspline_edge_with_index(step, &index, endpoint_tolerance_mm)
}

fn complete_bspline_edge_with_index(
    step: StepEdge<'_>,
    index: &std::collections::HashMap<u64, usize>,
    endpoint_tolerance_mm: f64,
) -> Result<NurbsCurve<Vector4>> {
    if !endpoint_tolerance_mm.is_finite() || endpoint_tolerance_mm <= 0.0 {
        bail!("edge endpoint tolerance must be positive and finite");
    }
    let params = named_entity_parameters(step.id, "EDGE_CURVE", step.entities, index)
        .with_context(|| format!("EDGE_CURVE #{} is missing or malformed", step.id))?;
    let start_id = params
        .get(1)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow!("edge #{} has no starting vertex", step.id))?;
    let end_id = params
        .get(2)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow!("edge #{} has no ending vertex", step.id))?;
    let curve_id = params
        .get(3)
        .and_then(entity_ref_value)
        .ok_or_else(|| anyhow!("edge #{} has no geometric curve", step.id))?;
    let start = vertex_point(start_id, step.entities, index)
        .ok_or_else(|| anyhow!("edge #{} has invalid starting vertex", step.id))?;
    let end = vertex_point(end_id, step.entities, index)
        .ok_or_else(|| anyhow!("edge #{} has invalid ending vertex", step.id))?;
    let CurveSupport::BSpline(spline) = brep::curve_support(curve_id, step.entities, index) else {
        bail!("edge #{} has no supported full B-spline geometry", step.id);
    };
    let count = spline.control_points_mm.len();
    let weights = spline.weights.unwrap_or_else(|| vec![1.0; count]);
    let first_control = spline
        .control_points_mm
        .first()
        .copied()
        .ok_or_else(|| anyhow!("edge #{} has no B-spline controls", step.id))?;
    let last_control = spline
        .control_points_mm
        .last()
        .copied()
        .ok_or_else(|| anyhow!("edge #{} has no B-spline controls", step.id))?;
    if !has_strict_chord_progress(
        &spline.control_points_mm,
        &weights,
        spline.degree,
        &spline.knots,
        Point3::new(first_control[0], first_control[1], first_control[2]),
        Point3::new(last_control[0], last_control[1], last_control[2]),
    ) {
        bail!(
            "edge #{} lacks an unambiguous full-domain monotone trim certificate",
            step.id
        );
    }
    let controls = spline
        .control_points_mm
        .into_iter()
        .zip(weights)
        .map(|(point, weight)| {
            Vector4::new(
                point[0] * weight,
                point[1] * weight,
                point[2] * weight,
                weight,
            )
        })
        .collect::<Vec<_>>();
    let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::from(spline.knots),
        controls,
    ));
    let (first, last) = curve.range_tuple();
    if !first.is_finite() || !last.is_finite() || first >= last {
        bail!("edge #{} has invalid active B-spline range", step.id);
    }
    let first_point = curve.evaluate(first);
    let last_point = curve.evaluate(last);
    if point_distance(first_point, [last_point.x, last_point.y, last_point.z])
        <= endpoint_tolerance_mm
    {
        bail!(
            "edge #{} has ambiguous coincident geometric endpoints",
            step.id
        );
    }
    let forward = point_distance(first_point, start) <= endpoint_tolerance_mm
        && point_distance(last_point, end) <= endpoint_tolerance_mm;
    let backward = point_distance(first_point, end) <= endpoint_tolerance_mm
        && point_distance(last_point, start) <= endpoint_tolerance_mm;
    if !forward && !backward {
        bail!(
            "edge #{} trims a subrange of its B-spline geometry",
            step.id
        );
    }
    Ok(curve)
}

/// Resolve a genuine planar carrier incident on a specified source edge.
///
/// A matching plane elsewhere in the same document is insufficient: this
/// validates the oriented topological edge's actual face-loop membership.
///
/// # Errors
/// Returns an error for malformed or nonincident face/edge topology, or a
/// missing or nonplanar source surface.
pub fn incident_edge_plane(
    face_id: u64,
    edge_id: u64,
    entities: &[EntityInstance],
) -> Result<ReferencePlane> {
    let index = build_index(entities);
    let loops = brep::face_loops(face_id, entities, &index)
        .with_context(|| format!("face #{face_id} has unsupported loops"))?;
    if !loops
        .iter()
        .flat_map(|boundary| &boundary.edges)
        .any(|edge| edge.edge_id == edge_id)
    {
        bail!("edge #{edge_id} is not incident on face #{face_id}");
    }
    let surface_id = brep::face_surface(face_id, entities, &index)
        .with_context(|| format!("face #{face_id} has no surface"))?;
    let brep::SurfaceSupport::Plane(plane) = brep::surface_support(surface_id, entities, &index)
    else {
        bail!("face #{face_id} does not have an analytic planar carrier");
    };
    Ok(ReferencePlane {
        origin_mm: plane.origin_mm,
        normal: plane.normal,
    })
}

/// Establish a lower bound on the directed distance between two STEP edges.
///
/// The source and target may use different degrees, knots and parameter
/// speeds. A positive result is a geometric non-equivalence witness for
/// their full 3-D curves. No claim is made about the surfaces/PCURVEs or
/// whether the two edges have equivalent oriented topology.
///
/// # Errors
/// Returns an error when either edge lacks a usable untrimmed B-spline,
/// has ambiguous endpoint coverage, or the plane bracket is unsupported.
pub fn bound_complete_step_edges_separation_by_plane(
    source: StepEdge<'_>,
    target: StepEdge<'_>,
    plane: ReferencePlane,
    endpoint_tolerance_mm: f64,
    target_gap_mm: f64,
) -> Result<PlaneDirectedSeparationBounds> {
    let source_curve = complete_bspline_edge(source, endpoint_tolerance_mm)?;
    let target_curve = complete_bspline_edge(target, endpoint_tolerance_mm)?;
    let origin = Point3::new(plane.origin_mm[0], plane.origin_mm[1], plane.origin_mm[2]);
    let normal = Vector3::new(plane.normal[0], plane.normal[1], plane.normal[2]);
    bound_nurbs_curves_separation_by_plane(
        &source_curve,
        &target_curve,
        origin,
        normal,
        target_gap_mm,
    )
    .ok_or_else(|| anyhow!("whole-edge plane separation is numerically unsupported"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_endpoints_alone_do_not_certify_a_revisiting_curve() {
        let spline = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        assert!(!has_strict_chord_progress(
            &spline,
            &[1.0; 4],
            3,
            &[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ));
        let ordered = [
            [0.0, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.7, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ];
        assert!(has_strict_chord_progress(
            &ordered,
            &[1.0, 2.0, 0.5, 1.0],
            3,
            &[0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ));
        // An internal degree+1 knot repetition introduces a discontinuity:
        // an ordered list of controls cannot authorize two disconnected loci.
        assert!(!has_strict_chord_progress(
            &ordered,
            &[1.0; 4],
            1,
            &[0.0, 0.0, 0.5, 0.5, 1.0, 1.0],
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ));
    }

    #[test]
    fn manufacturer_cubic_has_full_domain_but_is_not_planar() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/production_radial_slot.step");
        let text = std::str::from_utf8(bytes)?;
        let source = step_p21::parser::parse(text)?;
        let [section] = source.data.as_slice() else {
            bail!("manufacturer fixture must have one DATA section");
        };
        let entities = &section.entities;
        let edge = StepEdge {
            id: 93481,
            entities,
        };
        let plane = incident_edge_plane(125603, edge.id, entities)?;
        let bounds =
            bound_complete_step_edges_separation_by_plane(edge, edge, plane, 1.0e-8, 1.0e-10)?;
        assert!(bounds.source.lower > 1.5706e-4, "{bounds:?}");
        assert_eq!(bounds.minimum_source_to_target_mm, 0.0);
        assert!(incident_edge_plane(125603, 11, entities).is_err());

        // The underlying source curve is still the same, but this broken
        // topological EDGE_CURVE now claims a different (zero-length) trim.
        // The full spline excursion cannot be substituted for that edge.
        let tampered = text.replacen(
            "#93481=EDGE_CURVE('NONE',#11689,#34532,#60612,.T.);",
            "#93481=EDGE_CURVE('NONE',#11689,#11689,#60612,.T.);",
            1,
        );
        assert_ne!(tampered, text);
        let altered = step_p21::parser::parse(&tampered)?;
        let [section] = altered.data.as_slice() else {
            bail!("changed fixture must have one DATA section");
        };
        assert!(
            complete_bspline_edge(
                StepEdge {
                    id: 93481,
                    entities: &section.entities,
                },
                1.0e-8,
            )
            .is_err()
        );
        Ok(())
    }
}

/// Resolve the source solid's declared STEP coordinate length unit.
///
/// The input geometry API consumes the source coordinates as written in STEP.
/// A bounding tolerance in millimetres is only meaningful after resolving the
/// unit from the solid's own B-rep shape representation context. A missing or
/// contradictory representation is an error, not an implicit mm assumption.
///
/// # Errors
/// Returns an error for absent, malformed or inconsistent source length units.
pub fn source_solid_length_scale_mm(solid_id: u64, entities: &[EntityInstance]) -> Result<f64> {
    use crate::step_document::StepSectionView;
    use crate::step_entities::representation_items_and_context;
    use crate::step_graph::simple_record;

    let index = build_index(entities);
    let section = StepSectionView::new(entities, &index);
    let mut scale: Option<f64> = None;
    for entity in entities {
        let Some(record) = simple_record(entity) else {
            continue;
        };
        if record.name != "ADVANCED_BREP_SHAPE_REPRESENTATION" {
            continue;
        }
        let Some((items, context_id)) = representation_items_and_context(entity) else {
            bail!("malformed B-rep shape representation");
        };
        if !items.contains(&solid_id) {
            continue;
        }
        let units = section
            .entity(context_id)
            .and_then(|entry| entry.record_named("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
            .and_then(|record| record.value().list())
            .ok_or_else(|| anyhow!("solid #{solid_id} has no valid unit context"))?;
        let unit_refs = units
            .get(0)
            .and_then(crate::step_document::StepValue::list)
            .filter(|_| units.len() == 1)
            .ok_or_else(|| anyhow!("solid #{solid_id} has malformed unit context entries"))?;
        let mut context_scale = None;
        for unit in unit_refs.iter() {
            let unit_id = unit
                .entity_ref()
                .ok_or_else(|| anyhow!("solid #{solid_id} has a malformed unit reference"))?;
            if section
                .entity(unit_id)
                .and_then(|entry| entry.record_named("LENGTH_UNIT"))
                .is_none()
            {
                continue;
            }
            let value = crate::units::length_unit_scale_mm(unit_id, &section)
                .ok_or_else(|| anyhow!("solid #{solid_id} has an unresolvable length unit"))?;
            if context_scale.replace(value).is_some() {
                bail!("solid #{solid_id} has multiple length units in one context");
            }
        }
        let value = context_scale
            .ok_or_else(|| anyhow!("solid #{solid_id} has no declared length unit"))?;
        if let Some(prior) = scale.replace(value)
            && prior.total_cmp(&value).is_ne()
        {
            bail!("solid #{solid_id} has conflicting length unit contexts");
        }
    }
    scale.ok_or_else(|| anyhow!("solid #{solid_id} has no B-rep shape representation"))
}

/// A whole-domain maximum-distance bracket for a B-spline on an incident plane.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize)]
pub struct SourcePlaneDeviation {
    pub lower_mm: f64,
    pub upper_mm: f64,
    pub arithmetic_margin_mm: f64,
    pub subdivisions: usize,
}

impl From<monstertruck_geometry::analytic_plane_deviation::PlaneDeviationBounds>
    for SourcePlaneDeviation
{
    fn from(value: monstertruck_geometry::analytic_plane_deviation::PlaneDeviationBounds) -> Self {
        Self {
            lower_mm: value.lower,
            upper_mm: value.upper,
            arithmetic_margin_mm: value.arithmetic_margin,
            subdivisions: value.subdivisions,
        }
    }
}

/// Undetermined is neither a clean source plane nor a certified defect.
/// Unsupported indicates a failure of the full-edge curve/domain proof.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SourcePlaneIncidenceStatus {
    WithinTolerance,
    ExceedsTolerance,
    Undetermined,
    Unsupported,
}

/// A real STEP edge/face incidence, not a geometric-proximity guess.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct SourcePlaneIncidence {
    pub face_id: u64,
    pub edge_id: u64,
    pub status: SourcePlaneIncidenceStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bound: Option<SourcePlaneDeviation>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

/// Audit source B-spline edges on authored planar faces.
///
/// Inspects complete B-spline incidences on literal STEP `PLANE` faces.
/// Other edge carriers are counted, not declared conformant. This does not
/// establish or admit an editable CAD reconstruction.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct SourcePlanarEdgeAudit {
    pub solid_id: u64,
    pub tolerance_mm: f64,
    pub coordinate_scale_mm: f64,
    pub planar_faces: usize,
    pub non_bspline_edge_incidences: usize,
    pub unresolved_curve_edge_incidences: usize,
    pub incidences: Vec<SourcePlaneIncidence>,
}

impl SourcePlanarEdgeAudit {
    /// Applies only to examined B-spline incidence loci, never whole solids.
    #[must_use]
    pub fn checked_bspline_incidences_are_within_tolerance(&self) -> bool {
        !self.incidences.is_empty()
            && self.unresolved_curve_edge_incidences == 0
            && self
                .incidences
                .iter()
                .all(|item| item.status == SourcePlaneIncidenceStatus::WithinTolerance)
    }
}

#[derive(Clone, Copy)]
struct PlaneIncidenceSettings {
    endpoint_tolerance: f64,
    source_target_gap: f64,
    coordinate_scale_mm: f64,
    tolerance_mm: f64,
}

fn check_bspline_planar_incidence(
    edge: StepEdge<'_>,
    index: &std::collections::HashMap<u64, usize>,
    face_id: u64,
    plane: brep::PlaneSupport,
    settings: PlaneIncidenceSettings,
) -> SourcePlaneIncidence {
    let (status, bound, reason) =
        match complete_bspline_edge_with_index(edge, index, settings.endpoint_tolerance) {
            Err(error) => (
                SourcePlaneIncidenceStatus::Unsupported,
                None,
                Some(error.to_string()),
            ),
            Ok(curve) => {
                let origin =
                    Point3::new(plane.origin_mm[0], plane.origin_mm[1], plane.origin_mm[2]);
                let normal = Vector3::new(plane.normal[0], plane.normal[1], plane.normal[2]);
                monstertruck_geometry::analytic_plane_deviation::bound_nurbs_plane_deviation(
                    &curve,
                    origin,
                    normal,
                    settings.source_target_gap,
                )
                .and_then(|bounds| {
                    let mut deviation: SourcePlaneDeviation = bounds.into();
                    deviation.lower_mm *= settings.coordinate_scale_mm;
                    deviation.upper_mm *= settings.coordinate_scale_mm;
                    deviation.arithmetic_margin_mm *= settings.coordinate_scale_mm;
                    (deviation.lower_mm.is_finite()
                        && deviation.upper_mm.is_finite()
                        && deviation.arithmetic_margin_mm.is_finite())
                    .then_some(deviation)
                })
                .map_or_else(
                    || {
                        (
                            SourcePlaneIncidenceStatus::Unsupported,
                            None,
                            Some("numerically unsupported whole-curve plane deviation".to_owned()),
                        )
                    },
                    |deviation| {
                        let status = if deviation.lower_mm > settings.tolerance_mm {
                            SourcePlaneIncidenceStatus::ExceedsTolerance
                        } else if deviation.upper_mm <= settings.tolerance_mm {
                            SourcePlaneIncidenceStatus::WithinTolerance
                        } else {
                            SourcePlaneIncidenceStatus::Undetermined
                        };
                        (status, Some(deviation), None)
                    },
                )
            }
        };
    SourcePlaneIncidence {
        face_id,
        edge_id: edge.id,
        status,
        bound,
        reason,
    }
}

/// Analyze every literal analytic planar face and its incident source B-spline
/// edge loci, without an external CAD kernel or a reconstructed candidate.
///
/// All plane faces must have resolvable topology. Ambiguous B-spline geometry
/// is reported as `Unsupported`, not silently dropped. `ExceedsTolerance` proves
/// an internal source inconsistency independently of later projection/export.
///
/// # Errors
/// Returns an error for missing or malformed source solid/plane topology,
/// non-finite/negative tolerances, or unsupported source plane coordinates.
pub fn audit_source_solid_planar_bspline_incidences(
    solid_id: u64,
    entities: &[EntityInstance],
    tolerance_mm: f64,
    target_gap_mm: f64,
    coordinate_scale_mm: f64,
) -> Result<SourcePlanarEdgeAudit> {
    if !tolerance_mm.is_finite() || tolerance_mm < 0.0 {
        bail!("source plane incidence tolerance must be nonnegative and finite");
    }
    if !target_gap_mm.is_finite() || target_gap_mm <= 0.0 {
        bail!("source plane incidence target gap must be positive and finite");
    }
    if !coordinate_scale_mm.is_finite() || coordinate_scale_mm <= 0.0 {
        bail!("source coordinate millimetre scale must be positive and finite");
    }
    let source_target_gap = target_gap_mm / coordinate_scale_mm;
    let endpoint_tolerance = 1.0e-7 / coordinate_scale_mm;
    if !source_target_gap.is_finite()
        || source_target_gap <= 0.0
        || !endpoint_tolerance.is_finite()
        || endpoint_tolerance <= 0.0
    {
        bail!("source unit scale cannot represent required precision");
    }
    let settings = PlaneIncidenceSettings {
        endpoint_tolerance,
        source_target_gap,
        coordinate_scale_mm,
        tolerance_mm,
    };
    let index = build_index(entities);
    let face_ids = brep::solid_face_ids(solid_id, entities, &index)
        .with_context(|| format!("solid #{solid_id} has unsupported face topology"))?;
    let mut result = SourcePlanarEdgeAudit {
        solid_id,
        tolerance_mm,
        coordinate_scale_mm,
        planar_faces: 0,
        non_bspline_edge_incidences: 0,
        unresolved_curve_edge_incidences: 0,
        incidences: Vec::new(),
    };
    for face_id in face_ids {
        let surface_id = brep::face_surface(face_id, entities, &index)
            .with_context(|| format!("face #{face_id} has no surface"))?;
        // A spline approximated as planar is not an authored STEP PLANE.
        // A malformed literal PLANE, however, is a hard failure, not a skip.
        let literal_plane = index
            .get(&surface_id)
            .and_then(|&ordinal| entities.get(ordinal))
            .and_then(crate::step_graph::simple_record)
            .is_some_and(|record| record.name == "PLANE");
        if !literal_plane {
            continue;
        }
        let brep::SurfaceSupport::Plane(plane) =
            brep::surface_support(surface_id, entities, &index)
        else {
            bail!("face #{face_id} contains a malformed analytic plane");
        };
        let loops = brep::face_loops(face_id, entities, &index)
            .with_context(|| format!("face #{face_id} has unsupported edge loops"))?;
        let mut edges = std::collections::BTreeMap::new();
        for loop_ in loops {
            for edge in loop_.edges {
                edges.entry(edge.edge_id).or_insert(edge);
            }
        }
        result.planar_faces += 1;
        for (edge_id, edge) in edges {
            match edge.support {
                CurveSupport::BSpline(_) => {}
                CurveSupport::Other { .. } => {
                    result.unresolved_curve_edge_incidences += 1;
                    continue;
                }
                _ => {
                    result.non_bspline_edge_incidences += 1;
                    continue;
                }
            }
            result.incidences.push(check_bspline_planar_incidence(
                StepEdge {
                    id: edge_id,
                    entities,
                },
                &index,
                face_id,
                plane,
                settings,
            ));
        }
    }
    result
        .incidences
        .sort_by_key(|item| (item.face_id, item.edge_id));
    Ok(result)
}

#[cfg(test)]
mod source_planar_audit_tests {
    use super::*;

    fn simple_source(interior_z: f64) -> String {
        format!(
            "ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('source'),'2;1');
FILE_NAME('','',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=MANIFOLD_SOLID_BREP('',#2);
#2=CLOSED_SHELL('',(#3));
#3=ADVANCED_FACE('',(#4),#5,.T.);
#4=FACE_OUTER_BOUND('',#6,.T.);
#5=PLANE('',#16);
#6=EDGE_LOOP('',(#7));
#7=ORIENTED_EDGE('',*,*,#8,.T.);
#8=EDGE_CURVE('',#9,#10,#11,.T.);
#9=VERTEX_POINT('',#12);
#10=VERTEX_POINT('',#13);
#11=B_SPLINE_CURVE_WITH_KNOTS('',3,(#12,#14,#15,#13),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.);
#12=CARTESIAN_POINT('',(0.,0.,0.));
#13=CARTESIAN_POINT('',(1.,0.,0.));
#14=CARTESIAN_POINT('',(0.25,0.2,{interior_z}));
#15=CARTESIAN_POINT('',(0.75,0.2,{interior_z}));
#16=AXIS2_PLACEMENT_3D('',#17,#18,#19);
#17=CARTESIAN_POINT('',(0.,0.,0.));
#18=DIRECTION('',(0.,0.,1.));
#19=DIRECTION('',(1.,0.,0.));
ENDSEC;
END-ISO-10303-21;"
        )
    }

    fn inspect(data: &str) -> Result<SourcePlanarEdgeAudit> {
        let exchange = step_p21::parser::parse(data)?;
        let [section] = exchange.data.as_slice() else {
            bail!("expected a single DATA section");
        };
        audit_source_solid_planar_bspline_incidences(1, &section.entities, 1.0e-5, 1.0e-10, 1.0)
    }

    #[test]
    fn authored_planes_are_checked_over_full_bspline_spans() -> Result<()> {
        let exact = inspect(&simple_source(0.0))?;
        assert_eq!(exact.planar_faces, 1);
        assert_eq!(exact.non_bspline_edge_incidences, 0);
        assert_eq!(exact.incidences.len(), 1);
        assert_eq!(exact.incidences[0].edge_id, 8);
        assert_eq!(
            exact.incidences[0].status,
            SourcePlaneIncidenceStatus::WithinTolerance
        );
        assert!(exact.checked_bspline_incidences_are_within_tolerance());

        let off_plane = inspect(&simple_source(0.002))?;
        assert_eq!(off_plane.incidences.len(), 1);
        assert_eq!(
            off_plane.incidences[0].status,
            SourcePlaneIncidenceStatus::ExceedsTolerance
        );
        assert!(
            off_plane.incidences[0]
                .bound
                .is_some_and(|b| b.lower_mm > 0.001)
        );
        assert!(!off_plane.checked_bspline_incidences_are_within_tolerance());

        // The same source coordinates are within tolerance when their unit
        // is a micrometre, rather than a millimetre.
        let exchange = step_p21::parser::parse(&simple_source(0.002))?;
        let [section] = exchange.data.as_slice() else {
            bail!("expected one DATA section");
        };
        let scaled = audit_source_solid_planar_bspline_incidences(
            1,
            &section.entities,
            1.0e-5,
            1.0e-10,
            0.001,
        )?;
        assert_eq!(scaled.coordinate_scale_mm, 0.001);
        assert_eq!(
            scaled.incidences[0].status,
            SourcePlaneIncidenceStatus::WithinTolerance
        );
        assert!(scaled.checked_bspline_incidences_are_within_tolerance());
        Ok(())
    }

    #[test]
    fn ambiguous_edges_and_unsupported_carriers_are_not_assumed_conformant() -> Result<()> {
        let invalid = simple_source(0.0).replace(
            "#15=CARTESIAN_POINT('',(0.75,0.2,0));",
            "#15=CARTESIAN_POINT('',(-0.75,0.2,0));",
        );
        // The exact numeric representation in the synthetic fixture is
        // generated by Display, so use a geometric control edit.
        let invalid = if invalid == simple_source(0.0) {
            simple_source(0.0).replace(
                "#15=CARTESIAN_POINT('',(0.75,0.2,0.0));",
                "#15=CARTESIAN_POINT('',(-0.75,0.2,0.0));",
            )
        } else {
            invalid
        };
        let refused = inspect(&invalid)?;
        assert_eq!(
            refused.incidences[0].status,
            SourcePlaneIncidenceStatus::Unsupported
        );
        assert!(refused.incidences[0].reason.is_some());

        let line = simple_source(0.0).replace(
            "#11=B_SPLINE_CURVE_WITH_KNOTS('',3,(#12,#14,#15,#13),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.);",
            "#11=LINE('',#12,#20);\n#20=VECTOR('',#19,1.);",
        );
        let unexamined = inspect(&line)?;
        assert_eq!(unexamined.planar_faces, 1);
        assert_eq!(unexamined.non_bspline_edge_incidences, 1);
        assert!(unexamined.incidences.is_empty());
        assert!(!unexamined.checked_bspline_incidences_are_within_tolerance());
        assert_eq!(unexamined.unresolved_curve_edge_incidences, 0);

        let malformed = simple_source(0.0).replace(
            "#11=B_SPLINE_CURVE_WITH_KNOTS('',3,(#12,#14,#15,#13),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.);",
            "#11=UNKNOWN_GEOMETRIC_CURVE('');",
        );
        let unknown = inspect(&malformed)?;
        assert_eq!(unknown.non_bspline_edge_incidences, 0);
        assert_eq!(unknown.unresolved_curve_edge_incidences, 1);
        assert!(unknown.incidences.is_empty());
        assert!(!unknown.checked_bspline_incidences_are_within_tolerance());
        Ok(())
    }

    #[test]
    fn real_radial_slot_audit_finds_original_off_plane_source_edge() -> Result<()> {
        let source = std::str::from_utf8(include_bytes!(
            "../validation/fixtures/production_radial_slot.step"
        ))?;
        let exchange = step_p21::parser::parse(source)?;
        let [section] = exchange.data.as_slice() else {
            bail!("expected one DATA section");
        };
        assert_eq!(source_solid_length_scale_mm(15385, &section.entities)?, 1.0);
        let centimeter_source =
            source.replacen("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)", 1);
        assert_ne!(centimeter_source, source);
        let centimeter_exchange = step_p21::parser::parse(&centimeter_source)?;
        let [centimeter_section] = centimeter_exchange.data.as_slice() else {
            bail!("expected one DATA section");
        };
        assert_eq!(
            source_solid_length_scale_mm(15385, &centimeter_section.entities)?,
            10.0,
        );

        let audit = audit_source_solid_planar_bspline_incidences(
            15385,
            &section.entities,
            1.0e-5,
            1.0e-10,
            1.0,
        )?;
        let edge = audit
            .incidences
            .iter()
            .find(|edge| edge.edge_id == 93481 && edge.face_id == 125603)
            .ok_or_else(|| anyhow!("radial source edge was not audited"))?;
        assert_eq!(edge.status, SourcePlaneIncidenceStatus::ExceedsTolerance);
        assert!(
            edge.bound
                .is_some_and(|bounds| bounds.lower_mm > 0.00015706)
        );
        assert!(!audit.checked_bspline_incidences_are_within_tolerance());
        Ok(())
    }
}
