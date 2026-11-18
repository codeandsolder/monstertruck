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
    if !endpoint_tolerance_mm.is_finite() || endpoint_tolerance_mm <= 0.0 {
        bail!("edge endpoint tolerance must be positive and finite");
    }
    let index = build_index(step.entities);
    let params = named_entity_parameters(step.id, "EDGE_CURVE", step.entities, &index)
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
    let start = vertex_point(start_id, step.entities, &index)
        .ok_or_else(|| anyhow!("edge #{} has invalid starting vertex", step.id))?;
    let end = vertex_point(end_id, step.entities, &index)
        .ok_or_else(|| anyhow!("edge #{} has invalid ending vertex", step.id))?;
    let CurveSupport::BSpline(spline) = brep::curve_support(curve_id, step.entities, &index) else {
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
