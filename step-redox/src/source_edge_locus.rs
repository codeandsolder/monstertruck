//! Whole-domain STEP edge-locus comparison via a separating analytic plane.
//!
//! The input `EDGE_CURVE` must cover its entire referenced B-spline domain.
//! STEP permits topological trims of longer basis curves: those are rejected
//! here rather than treating an off-trim excursion as a real edge witness.
//! The geometric comparison is independent of knot/parameter correspondence.

use crate::brep::{self, BSplineSupport, CurveSupport, LineSupport};
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

/// Domain structure shared by both endpoint certificates. A full internal
/// multiplicity could introduce a discontinuity and is never accepted.
fn clamped_continuous_knots(degree: usize, knots: &[f64], weights: &[f64]) -> bool {
    if degree == 0 || degree > 64 {
        return false;
    }
    let count = weights.len();
    let order = degree + 1;
    if count < order || knots.len() != count + order {
        return false;
    }
    if !knots.iter().all(|knot| knot.is_finite()) || knots.windows(2).any(|w| w[0] > w[1]) {
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
    // A full internal multiplicity of degree + 1 can conceal a geometric
    // discontinuity. Endpoint matching must not authorize that union of loci.
    let domain_start = knots[degree];
    let domain_end = knots[count];
    if domain_start >= domain_end || knots[order] <= domain_start || knots[count - 1] >= domain_end
    {
        return false;
    }
    if knots.windows(degree + 1).any(|window| {
        window[0].total_cmp(&window[degree]).is_eq()
            && window[0] > domain_start
            && window[0] < domain_end
    }) {
        return false;
    }
    true
}

/// Show that the endpoint is strictly outside the convex hull of all
/// supplied non-endpoint controls. Any separating direction suffices;
/// candidate directions are not themselves the certificate.
fn endpoint_separated(endpoint: [f64; 3], controls: &[[f64; 3]], floor: f64) -> bool {
    use crate::math3::{add, dot, mul, norm, sub};
    let offsets = controls
        .iter()
        .map(|&point| sub(point, endpoint))
        .collect::<Vec<_>>();
    let Some(mut closest) = offsets
        .iter()
        .copied()
        .min_by(|&a, &b| norm(a).total_cmp(&norm(b)))
    else {
        return false;
    };
    // Frank-Wolfe minimum-norm point iterations only propose a direction.
    // A direct all-controls separating-plane inequality is the actual proof;
    // nonconvergence or near-zero separation refuses the trim.
    for _ in 0..128 {
        let length = norm(closest);
        if !length.is_finite() || length <= floor {
            return false;
        }
        let direction = mul(closest, length.recip());
        let Some(witness) = offsets
            .iter()
            .copied()
            .min_by(|&a, &b| dot(a, direction).total_cmp(&dot(b, direction)))
        else {
            return false;
        };
        if dot(witness, direction) > floor {
            return true;
        }
        let step = sub(witness, closest);
        let denominator = dot(step, step);
        if !denominator.is_finite() || denominator <= 0.0 {
            return false;
        }
        let fraction = (-dot(closest, step) / denominator).clamp(0.0, 1.0);
        if !fraction.is_finite() || fraction <= 0.0 {
            return false;
        }
        closest = add(closest, mul(step, fraction));
    }
    false
}

/// A nonmonotone curve may still have an unambiguous complete open trim.
/// Each genuinely positive knot span lies inside the convex hull of its
/// active controls. Exclude the opposite endpoint from every such hull;
/// exclude the matching endpoint from every hull except the endpoint span,
/// where its one clamped control has strictly less than unit basis weight
/// in the interior. This permits turning tangents without permitting an
/// additional visit to either endpoint. No sampled curve points are used.
fn unique_endpoint_visits(spline: &BSplineSupport, weights: &[f64]) -> bool {
    let points = &spline.control_points_mm;
    let degree = spline.degree;
    let knots = &spline.knots;
    if points.len() != weights.len()
        || !clamped_continuous_knots(degree, knots, weights)
        || points.iter().flatten().any(|v| !v.is_finite())
        || weights.iter().any(|&w| !w.is_finite() || w <= 0.)
    {
        return false;
    }
    let scale = points
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(1.0_f64, f64::max);
    let floor = 256. * f64::EPSILON * scale;
    let Some(&first) = points.first() else {
        return false;
    };
    let Some(&last) = points.last() else {
        return false;
    };
    (degree..points.len())
        .filter(|&i| knots[i] < knots[i + 1])
        .all(|i| {
            let start = i - degree;
            let first_control = usize::from(start == 0);
            let last_control = usize::from(i == points.len() - 1);
            endpoint_separated(first, &points[start + first_control..=i], floor)
                && endpoint_separated(last, &points[start..=i - last_control], floor)
        })
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
    if control_points.len() != weights.len() || !clamped_continuous_knots(degree, knots, weights) {
        return false;
    }
    let chord = [last.x - first.x, last.y - first.y, last.z - first.z];
    let len2 = chord[0].mul_add(chord[0], chord[1].mul_add(chord[1], chord[2] * chord[2]));
    if !len2.is_finite() || len2 <= 0.0 {
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

/// Canonical straight-edge support for recognizers, not a source-file edit.
/// Besides the finite locus certificate, require one unambiguous traversal
/// between the authored endpoints. A collinear but retracing spline must
/// not silently turn into an ordinary oriented edge of a feature.
pub(crate) fn certified_straight_spline(
    spline: &BSplineSupport,
    start: [f64; 3],
    end: [f64; 3],
    tolerance: f64,
) -> Option<LineSupport> {
    use monstertruck_geometry::analytic_line_segment::{SegmentLocus, certify_nurbs_segment_locus};
    let first = *spline.control_points_mm.first()?;
    let last = *spline.control_points_mm.last()?;
    let point = |p: [f64; 3]| Point3::new(p[0], p[1], p[2]);
    let unit_weights;
    let weights = if let Some(weights) = &spline.weights {
        weights.as_slice()
    } else {
        unit_weights = vec![1.0; spline.control_points_mm.len()];
        &unit_weights
    };
    if !has_strict_chord_progress(
        &spline.control_points_mm,
        weights,
        spline.degree,
        &spline.knots,
        point(first),
        point(last),
    ) {
        return None;
    }
    let rational = spline.rational_curve()?;
    if certify_nurbs_segment_locus(&rational, point(start), point(end), tolerance)
        != SegmentLocus::WithinTolerance
    {
        return None;
    }
    Some(LineSupport {
        origin_mm: start,
        direction: crate::math3::normalize(crate::math3::sub(end, start), tolerance)?,
    })
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
    complete_bspline_edge_with_index(
        step,
        &index,
        endpoint_tolerance_mm,
        EndpointProof::StrictChordProgress,
    )
}

#[derive(Clone, Copy)]
enum EndpointProof {
    StrictChordProgress,
    UniqueVisits,
}

fn complete_bspline_edge_with_index(
    step: StepEdge<'_>,
    index: &std::collections::HashMap<u64, usize>,
    endpoint_tolerance_mm: f64,
    proof: EndpointProof,
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
    complete_bspline_support_with_proof(spline, [start, end], endpoint_tolerance_mm, proof)
        .with_context(|| format!("edge #{}", step.id))
}

/// Construct the same full-domain, monotone, positive-weight spline from
/// an already resolved source edge (e.g. the native topology planner).
#[cfg(feature = "cad-kernel-monstertruck")]
pub(crate) fn complete_bspline_support(
    spline: brep::BSplineSupport,
    [start, end]: [[f64; 3]; 2],
    endpoint_tolerance_mm: f64,
) -> Result<NurbsCurve<Vector4>> {
    complete_bspline_support_with_proof(
        spline,
        [start, end],
        endpoint_tolerance_mm,
        EndpointProof::StrictChordProgress,
    )
}

fn complete_bspline_support_with_proof(
    spline: brep::BSplineSupport,
    [start, end]: [[f64; 3]; 2],
    endpoint_tolerance_mm: f64,
    proof: EndpointProof,
) -> Result<NurbsCurve<Vector4>> {
    if !endpoint_tolerance_mm.is_finite() || endpoint_tolerance_mm <= 0.0 {
        bail!("edge endpoint tolerance must be positive and finite");
    }
    let count = spline.control_points_mm.len();
    let unit_weights;
    let weights = if let Some(weights) = &spline.weights {
        weights.as_slice()
    } else {
        unit_weights = vec![1.0; count];
        &unit_weights
    };
    let first_control = spline
        .control_points_mm
        .first()
        .copied()
        .ok_or_else(|| anyhow!("edge has no B-spline controls"))?;
    let last_control = spline
        .control_points_mm
        .last()
        .copied()
        .ok_or_else(|| anyhow!("edge has no B-spline controls"))?;
    let monotone = has_strict_chord_progress(
        &spline.control_points_mm,
        weights,
        spline.degree,
        &spline.knots,
        Point3::new(first_control[0], first_control[1], first_control[2]),
        Point3::new(last_control[0], last_control[1], last_control[2]),
    );
    if !(monotone
        || matches!(proof, EndpointProof::UniqueVisits) && unique_endpoint_visits(&spline, weights))
    {
        bail!("edge lacks an unambiguous full-domain trim certificate");
    }
    let controls = spline
        .control_points_mm
        .into_iter()
        .zip(weights.iter().copied())
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
        bail!("edge has invalid active B-spline range");
    }
    let first_point = curve.evaluate(first);
    let last_point = curve.evaluate(last);
    if point_distance(first_point, [last_point.x, last_point.y, last_point.z])
        <= endpoint_tolerance_mm
    {
        bail!("edge has ambiguous coincident geometric endpoints");
    }
    let forward = point_distance(first_point, start) <= endpoint_tolerance_mm
        && point_distance(last_point, end) <= endpoint_tolerance_mm;
    let backward = point_distance(first_point, end) <= endpoint_tolerance_mm
        && point_distance(last_point, start) <= endpoint_tolerance_mm;
    if !forward && !backward {
        bail!("edge trims a subrange of its B-spline geometry");
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
    let brep::SurfaceSupport::Plane(plane) = brep::surface_support(
        surface_id,
        entities,
        &index,
        crate::source_accuracy::source_face_recovery_tolerance_mm(face_id, entities, &index),
    ) else {
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
    fn straight_spline_certificate_requires_complete_directed_segment() {
        let mut spline = BSplineSupport {
            control_points_mm: vec![[0., 0., 0.], [0.2, 0., 0.], [0.7, 0., 0.], [1., 0., 0.]],
            degree: 3,
            knots: vec![0., 0., 0., 0., 1., 1., 1., 1.],
            weights: Some(vec![1., 2., 0.5, 1.]),
        };
        let a = [0., 0., 0.];
        let b = [1., 0., 0.];
        let tolerance = 1.0e-7;
        assert!(certified_straight_spline(&spline, a, b, tolerance).is_some());
        assert!(certified_straight_spline(&spline, b, a, tolerance).is_some());
        assert!(certified_straight_spline(&spline, a, [0.9, 0., 0.], tolerance).is_none());
        spline.control_points_mm[1][1] = 0.001;
        assert!(certified_straight_spline(&spline, a, b, tolerance).is_none());
        spline.control_points_mm[1] = [0.8, 0., 0.];
        assert!(certified_straight_spline(&spline, a, b, tolerance).is_none());
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

/// A whole-domain maximum shortest-distance bracket for a B-spline on
/// its literal incident carrier, not an estimate of whole-body equivalence.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize)]
pub struct SourceCarrierDeviation {
    pub lower_mm: f64,
    pub upper_mm: f64,
    pub arithmetic_margin_mm: f64,
    pub subdivisions: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SourceCarrierKind {
    Plane,
    Cylinder,
    Cone,
}

#[derive(Clone, Copy)]
enum AuditedCarrier {
    Plane(brep::PlaneSupport),
    Cylinder(brep::CylinderSupport),
    Cone(brep::ConeSupport),
}

impl AuditedCarrier {
    const fn kind(self) -> SourceCarrierKind {
        match self {
            Self::Plane(_) => SourceCarrierKind::Plane,
            Self::Cylinder(_) => SourceCarrierKind::Cylinder,
            Self::Cone(_) => SourceCarrierKind::Cone,
        }
    }

    fn bound(self, curve: &NurbsCurve<Vector4>, gap: f64) -> Option<SourceCarrierDeviation> {
        use monstertruck_geometry::{
            analytic_plane_deviation::bound_nurbs_plane_deviation,
            analytic_radial_deviation::bound_nurbs_radial_deviation,
        };
        let point = |p: [f64; 3]| Point3::new(p[0], p[1], p[2]);
        let vector = |v: [f64; 3]| Vector3::new(v[0], v[1], v[2]);
        let radial = |origin, axis, radius, slope| {
            let bounds = bound_nurbs_radial_deviation(
                curve,
                point(origin),
                vector(axis),
                radius,
                slope,
                gap,
            )?;
            Some(SourceCarrierDeviation {
                lower_mm: bounds.lower,
                upper_mm: bounds.upper,
                arithmetic_margin_mm: bounds.arithmetic_margin,
                subdivisions: bounds.subdivisions,
            })
        };
        match self {
            Self::Plane(plane) => {
                let bounds = bound_nurbs_plane_deviation(
                    curve,
                    point(plane.origin_mm),
                    vector(plane.normal),
                    gap,
                )?;
                Some(SourceCarrierDeviation {
                    lower_mm: bounds.lower,
                    upper_mm: bounds.upper,
                    arithmetic_margin_mm: bounds.arithmetic_margin,
                    subdivisions: bounds.subdivisions,
                })
            }
            Self::Cylinder(cylinder) => radial(
                cylinder.axis_origin_mm,
                cylinder.axis,
                cylinder.radius_mm,
                0.0,
            ),
            Self::Cone(cone) => radial(
                cone.reference_origin_mm,
                cone.axis,
                cone.reference_radius_mm,
                cone.semi_angle_rad.tan(),
            ),
        }
    }
}

/// Undetermined is neither a conforming source incidence nor a certified defect.
/// Unsupported indicates a failure of the full-edge curve/domain proof.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SourceCarrierIncidenceStatus {
    WithinTolerance,
    ExceedsTolerance,
    Undetermined,
    Unsupported,
}

/// A real STEP edge/face incidence, not a geometric-proximity guess.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct SourceCarrierIncidence {
    pub face_id: u64,
    pub edge_id: u64,
    pub carrier: SourceCarrierKind,
    pub status: SourceCarrierIncidenceStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bound: Option<SourceCarrierDeviation>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

/// Audit source B-spline edges on authored planes, cylinders and cones.
///
/// Inspects complete B-spline incidences on supported literal STEP carriers.
/// Other edge carriers are counted, not declared conformant. This does not
/// establish or admit an editable CAD reconstruction.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct SourceSplineEdgeAudit {
    pub solid_id: u64,
    pub tolerance_mm: f64,
    pub coordinate_scale_mm: f64,
    pub planar_faces: usize,
    pub cylindrical_faces: usize,
    pub conical_faces: usize,
    pub unexamined_faces: usize,
    pub declared_uncertainty_mm: Option<f64>,
    pub non_bspline_edge_incidences: usize,
    pub unresolved_curve_edge_incidences: usize,
    pub incidences: Vec<SourceCarrierIncidence>,
}

impl SourceSplineEdgeAudit {
    /// Applies only to examined B-spline incidence loci, never whole solids.
    #[must_use = "The boolean result of `checked_bspline_incidences_are_within_tolerance` must be observed to use the predicate."]
    pub fn checked_bspline_incidences_are_within_tolerance(&self) -> bool {
        !self.incidences.is_empty()
            && self.unresolved_curve_edge_incidences == 0
            && self
                .incidences
                .iter()
                .all(|item| item.status == SourceCarrierIncidenceStatus::WithinTolerance)
    }
}

#[derive(Clone, Copy)]
struct CarrierIncidenceSettings {
    endpoint_tolerance: f64,
    source_target_gap: f64,
    coordinate_scale_mm: f64,
    tolerance_mm: f64,
}

fn check_bspline_carrier_incidence(
    edge: StepEdge<'_>,
    index: &std::collections::HashMap<u64, usize>,
    face_id: u64,
    carrier: AuditedCarrier,
    settings: CarrierIncidenceSettings,
) -> SourceCarrierIncidence {
    let (status, bound, reason) = match complete_bspline_edge_with_index(
        edge,
        index,
        settings.endpoint_tolerance,
        EndpointProof::UniqueVisits,
    ) {
        Err(error) => (
            SourceCarrierIncidenceStatus::Unsupported,
            None,
            Some(error.to_string()),
        ),
        Ok(curve) => carrier
            .bound(&curve, settings.source_target_gap)
            .and_then(|bounds| {
                let mut deviation = bounds;
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
                        SourceCarrierIncidenceStatus::Unsupported,
                        None,
                        Some("numerically unsupported whole-curve carrier deviation".to_owned()),
                    )
                },
                |deviation| {
                    let status = if deviation.lower_mm > settings.tolerance_mm {
                        SourceCarrierIncidenceStatus::ExceedsTolerance
                    } else if deviation.upper_mm <= settings.tolerance_mm {
                        SourceCarrierIncidenceStatus::WithinTolerance
                    } else {
                        SourceCarrierIncidenceStatus::Undetermined
                    };
                    (status, Some(deviation), None)
                },
            ),
    };
    SourceCarrierIncidence {
        face_id,
        edge_id: edge.id,
        carrier: carrier.kind(),
        status,
        bound,
        reason,
    }
}

/// Resolve literal source carriers, not analytic fits of authored splines.
fn audited_carrier(
    surface_id: u64,
    entities: &[EntityInstance],
    index: &std::collections::HashMap<u64, usize>,
) -> Result<Option<AuditedCarrier>> {
    let name = index
        .get(&surface_id)
        .and_then(|&i| entities.get(i))
        .and_then(crate::step_graph::simple_record)
        .map(|record| record.name.as_str());
    let carrier = match name {
        Some("PLANE" | "CYLINDRICAL_SURFACE" | "CONICAL_SURFACE") => {
            match (
                name,
                brep::surface_support(
                    surface_id,
                    entities,
                    index,
                    monstertruck_core::tolerance::NATIVE_STEP_RECOVERY_MM,
                ),
            ) {
                (Some("PLANE"), brep::SurfaceSupport::Plane(plane)) => AuditedCarrier::Plane(plane),
                (Some("CYLINDRICAL_SURFACE"), brep::SurfaceSupport::Cylinder(cylinder)) => {
                    AuditedCarrier::Cylinder(cylinder)
                }
                (Some("CONICAL_SURFACE"), brep::SurfaceSupport::Cone(cone)) => {
                    AuditedCarrier::Cone(cone)
                }
                _ => bail!("surface #{surface_id} contains a malformed analytic carrier"),
            }
        }
        _ => return Ok(None),
    };
    Ok(Some(carrier))
}

/// Analyze every supported literal planar, cylindrical or conical face and its incident source B-spline
/// edge loci, without an external CAD kernel or a reconstructed candidate.
///
/// All examined faces must have resolvable topology. Ambiguous B-spline geometry
/// is reported as `Unsupported`, not silently dropped. `ExceedsTolerance` proves
/// an internal source inconsistency independently of later projection/export.
///
/// # Errors
/// Returns an error for missing or malformed source solid/carrier topology,
/// non-finite/negative tolerances, or unsupported source carrier coordinates.
pub fn audit_source_solid_bspline_incidences(
    solid_id: u64,
    entities: &[EntityInstance],
    tolerance_mm: f64,
    target_gap_mm: f64,
    coordinate_scale_mm: f64,
) -> Result<SourceSplineEdgeAudit> {
    if !tolerance_mm.is_finite() || tolerance_mm < 0.0 {
        bail!("source carrier incidence tolerance must be nonnegative and finite");
    }
    if !target_gap_mm.is_finite() || target_gap_mm <= 0.0 {
        bail!("source carrier incidence target gap must be positive and finite");
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
    let settings = CarrierIncidenceSettings {
        endpoint_tolerance,
        source_target_gap,
        coordinate_scale_mm,
        tolerance_mm,
    };
    let index = build_index(entities);
    let face_ids = brep::solid_face_ids(solid_id, entities, &index)
        .with_context(|| format!("solid #{solid_id} has unsupported face topology"))?;
    let mut result = SourceSplineEdgeAudit {
        solid_id,
        tolerance_mm,
        coordinate_scale_mm,
        planar_faces: 0,
        cylindrical_faces: 0,
        conical_faces: 0,
        unexamined_faces: 0,
        declared_uncertainty_mm:
            crate::source_accuracy::declared_source_uncertainty_by_representation_item(
                entities, &index,
            )
            .get(&solid_id)
            .copied(),
        non_bspline_edge_incidences: 0,
        unresolved_curve_edge_incidences: 0,
        incidences: Vec::new(),
    };
    for face_id in face_ids {
        let surface_id = brep::face_surface(face_id, entities, &index)
            .with_context(|| format!("face #{face_id} has no surface"))?;
        let Some(carrier) = audited_carrier(surface_id, entities, &index)? else {
            result.unexamined_faces += 1;
            continue;
        };
        match carrier.kind() {
            SourceCarrierKind::Plane => result.planar_faces += 1,
            SourceCarrierKind::Cylinder => result.cylindrical_faces += 1,
            SourceCarrierKind::Cone => result.conical_faces += 1,
        }
        let loops = brep::face_loops(face_id, entities, &index)
            .with_context(|| format!("face #{face_id} has unsupported edge loops"))?;
        let mut edges = std::collections::BTreeMap::new();
        for loop_ in loops {
            for edge in loop_.edges {
                edges.entry(edge.edge_id).or_insert(edge);
            }
        }
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
            result.incidences.push(check_bspline_carrier_incidence(
                StepEdge {
                    id: edge_id,
                    entities,
                },
                &index,
                face_id,
                carrier,
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
mod source_carrier_audit_tests {
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

    fn inspect(data: &str) -> Result<SourceSplineEdgeAudit> {
        let exchange = step_p21::parser::parse(data)?;
        let [section] = exchange.data.as_slice() else {
            bail!("expected a single DATA section");
        };
        audit_source_solid_bspline_incidences(1, &section.entities, 1.0e-5, 1.0e-10, 1.0)
    }

    #[test]
    fn crossing_cylinder_source_splines_exceed_declared_uncertainty() -> Result<()> {
        let bytes = include_bytes!("../validation/fixtures/native_mixed_treatment_4a6356.step");
        let exchange = step_p21::parser::parse(std::str::from_utf8(bytes)?)?;
        let entities = &exchange.data[0].entities;
        let scale = source_solid_length_scale_mm(171, entities)?;
        let audit = audit_source_solid_bspline_incidences(171, entities, 1e-5, 1e-9, scale)?;
        assert_eq!(audit.declared_uncertainty_mm, Some(1e-5));
        for edge_id in [1049, 1063] {
            let incidences = audit
                .incidences
                .iter()
                .filter(|i| i.edge_id == edge_id)
                .collect::<Vec<_>>();
            assert_eq!(incidences.len(), 2);
            for incidence in incidences {
                assert_eq!(incidence.carrier, SourceCarrierKind::Cylinder);
                assert_eq!(
                    incidence.status,
                    SourceCarrierIncidenceStatus::ExceedsTolerance,
                    "{incidence:?}"
                );
                assert!(
                    incidence
                        .bound
                        .is_some_and(|bound| bound.lower_mm > 1e-4 && bound.upper_mm < 2e-4),
                    "{incidence:?}"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn nonmonotone_trim_is_audited_when_endpoint_hulls_are_separated() -> Result<()> {
        let bowed = simple_source(0.0).replace(
            "#14=CARTESIAN_POINT('',(0.25,0.2,0));",
            "#14=CARTESIAN_POINT('',(-0.25,0.2,0.));",
        );
        assert_ne!(bowed, simple_source(0.0));
        let audit = inspect(&bowed)?;
        assert_eq!(
            audit.incidences[0].status,
            SourceCarrierIncidenceStatus::WithinTolerance
        );
        Ok(())
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
            SourceCarrierIncidenceStatus::WithinTolerance
        );
        assert!(exact.checked_bspline_incidences_are_within_tolerance());

        let off_plane = inspect(&simple_source(0.002))?;
        assert_eq!(off_plane.incidences.len(), 1);
        assert_eq!(
            off_plane.incidences[0].status,
            SourceCarrierIncidenceStatus::ExceedsTolerance
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
        let scaled =
            audit_source_solid_bspline_incidences(1, &section.entities, 1.0e-5, 1.0e-10, 0.001)?;
        assert_eq!(scaled.coordinate_scale_mm, 0.001);
        assert_eq!(
            scaled.incidences[0].status,
            SourceCarrierIncidenceStatus::WithinTolerance
        );
        assert!(scaled.checked_bspline_incidences_are_within_tolerance());
        Ok(())
    }

    #[test]
    fn ambiguous_edges_and_unsupported_carriers_are_not_assumed_conformant() -> Result<()> {
        let invalid = simple_source(0.0)
            .replace(
                "#14=CARTESIAN_POINT('',(0.25,0.2,0));",
                "#14=CARTESIAN_POINT('',(-0.25,0.,0.));",
            )
            .replace(
                "#15=CARTESIAN_POINT('',(0.75,0.2,0));",
                "#15=CARTESIAN_POINT('',(-0.75,0.,0.));",
            );
        assert_ne!(invalid, simple_source(0.0));
        let refused = inspect(&invalid)?;
        assert_eq!(
            refused.incidences[0].status,
            SourceCarrierIncidenceStatus::Unsupported
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

        let audit =
            audit_source_solid_bspline_incidences(15385, &section.entities, 1.0e-5, 1.0e-10, 1.0)?;
        let edge = audit
            .incidences
            .iter()
            .find(|edge| edge.edge_id == 93481 && edge.face_id == 125603)
            .ok_or_else(|| anyhow!("radial source edge was not audited"))?;
        assert_eq!(edge.status, SourceCarrierIncidenceStatus::ExceedsTolerance);
        assert!(
            edge.bound
                .is_some_and(|bounds| bounds.lower_mm > 0.00015706)
        );
        assert!(!audit.checked_bspline_incidences_are_within_tolerance());
        Ok(())
    }
}
