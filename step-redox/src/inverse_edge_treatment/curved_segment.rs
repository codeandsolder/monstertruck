//! Replace supplier B-spline approximations of an analytic plane/cylinder
//! intersection by the exact conic when extending an inverse planar chamfer.
//!
//! The fitted source curve and its adjoining cylinder are checked over every
//! knot span. An analytic ellipse is rebuilt with a new intersection endpoint,
//! preserving the original edge's `StableId` on both neighboring source faces.
//! This changes the *geometric carrier*, not the numerical patch tolerance.

use super::{
    CurveSupport, DIR_TOL, EdgeTreatmentEvidence, EndpointEdit, GEOM_TOL_MM, KernelFaceSplice,
    KernelNewEdge, KernelPatchCurve, KernelPatchEdgeRef, KernelPatchEdgeUse, KernelPatchVertexRef,
    KernelTopologyPatch, KernelTopologyProvenance, OrientedEdgeUse, ReducedFace, ReducedSolid,
    SurfaceSupport, add, carrier_run_neighbors, cross, distance, dot, edge_face_users, face,
    global_edge_endpoints, mul, norm, normalize, normalize_retained_edge_rewrites,
    plan_terminal_splices, sub, topology_edges, topology_face, unique_face_edges,
};
use crate::brep::{BSplineSupport, CylinderSupport, PlaneSupport};
use anyhow::{Context, Result, bail, ensure};
use monstertruck_core::tolerance::DEFAULT_STEP_RECOVERY_TOLERANCES;
use std::collections::BTreeSet;

const FIT_MM: f64 = crate::numeric::DEFAULT_FEATURE_HEAL_TOL_MM;
const SAMPLES_PER_KNOT_SPAN: u32 = 32;

#[derive(Clone, Copy)]
struct CylinderPlaneConic {
    center: [f64; 3],
    cosine: [f64; 3],
    sine: [f64; 3],
    cylinder: CylinderSupport,
    axis: [f64; 3],
    radial_u: [f64; 3],
    radial_v: [f64; 3],
}

impl CylinderPlaneConic {
    fn new(plane: PlaneSupport, cylinder: CylinderSupport) -> Result<Self> {
        let normal = normalize(plane.normal, DIR_TOL).context("degenerate plane normal")?;
        let axis = normalize(cylinder.axis, DIR_TOL).context("degenerate cylinder axis")?;
        let denominator = dot(normal, axis);
        ensure!(
            denominator.abs() >= DEFAULT_STEP_RECOVERY_TOLERANCES.source_rounding_arc_residual_rad,
            "plane parallel to cylindrical rulings: intersection is not a conic"
        );
        let radial_u = normalize(
            sub(
                cylinder.x_direction,
                mul(axis, dot(cylinder.x_direction, axis)),
            ),
            DIR_TOL,
        )
        .context("cylinder has no independent radial reference")?;
        let radial_v = cross(axis, radial_u);
        let shift = dot(normal, sub(plane.origin_mm, cylinder.axis_origin_mm)) / denominator;
        let center = add(cylinder.axis_origin_mm, mul(axis, shift));
        let cosine = mul(
            sub(radial_u, mul(axis, dot(normal, radial_u) / denominator)),
            cylinder.radius_mm,
        );
        let sine = mul(
            sub(radial_v, mul(axis, dot(normal, radial_v) / denominator)),
            cylinder.radius_mm,
        );
        ensure!(
            [center, cosine, sine]
                .into_iter()
                .flatten()
                .all(f64::is_finite)
                && norm(cross(cosine, sine)) > GEOM_TOL_MM,
            "plane/cylinder conic is degenerate"
        );
        Ok(Self {
            center,
            cosine,
            sine,
            cylinder,
            axis,
            radial_u,
            radial_v,
        })
    }

    fn angle(self, point: [f64; 3]) -> f64 {
        let radial = sub(point, self.cylinder.axis_origin_mm);
        dot(radial, self.radial_v).atan2(dot(radial, self.radial_u))
    }

    fn point(self, angle: f64) -> [f64; 3] {
        add(
            self.center,
            add(mul(self.cosine, angle.cos()), mul(self.sine, angle.sin())),
        )
    }

    fn radial_residual(self, point: [f64; 3]) -> f64 {
        let relative = sub(point, self.cylinder.axis_origin_mm);
        let radial = sub(relative, mul(self.axis, dot(relative, self.axis)));
        (norm(radial) - self.cylinder.radius_mm).abs()
    }

    fn fitted_to_source(self, spline: &BSplineSupport, plane: PlaneSupport) -> Result<()> {
        let (first, last) = crate::boundary_proof::bspline_domain(spline)
            .context("retained B-spline has no finite supported knot interval")?;
        let normal = normalize(plane.normal, DIR_TOL).context("plane is degenerate")?;
        let mut intervals = vec![first];
        for &knot in &spline.knots {
            if knot > *intervals.last().context("no previous knot")? && knot < last {
                intervals.push(knot);
            }
        }
        intervals.push(last);
        for span in intervals.windows(2) {
            for step in 0..=SAMPLES_PER_KNOT_SPAN {
                let t = span[0].mul_add(
                    1.0 - f64::from(step) / f64::from(SAMPLES_PER_KNOT_SPAN),
                    span[1] * f64::from(step) / f64::from(SAMPLES_PER_KNOT_SPAN),
                );
                let point = crate::boundary_proof::evaluate_bspline(spline, t)
                    .context("retained spline cannot be evaluated")?;
                let plane_residual = dot(normal, sub(point, plane.origin_mm)).abs();
                ensure!(
                    plane_residual <= FIT_MM && self.radial_residual(point) <= FIT_MM,
                    "retained spline differs from plane/cylinder intersection by more than 10 um"
                );
            }
        }
        Ok(())
    }

    fn extend_to(self, point_on_virtual_line: [f64; 3], line: [f64; 3]) -> Result<[f64; 3]> {
        let relative = sub(point_on_virtual_line, self.cylinder.axis_origin_mm);
        let radial = sub(relative, mul(self.axis, dot(relative, self.axis)));
        let radial_dir = sub(line, mul(self.axis, dot(line, self.axis)));
        let a = dot(radial_dir, radial_dir);
        ensure!(
            a > DIR_TOL,
            "virtual line has no transverse cylinder intersection"
        );
        let b = 2.0 * dot(radial, radial_dir);
        let c = self
            .cylinder
            .radius_mm
            .mul_add(-self.cylinder.radius_mm, dot(radial, radial));
        let discriminant = b.mul_add(b, -4.0 * a * c);
        ensure!(discriminant >= 0.0, "virtual line misses cylinder entirely");
        let root = discriminant.sqrt();
        let t = [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]
            .into_iter()
            .min_by(|a, b| a.abs().total_cmp(&b.abs()))
            .context("no intersection with cylinder")?;
        let snapped = add(point_on_virtual_line, mul(line, t));
        ensure!(
            distance(snapped, point_on_virtual_line) <= FIT_MM
                && self.radial_residual(snapped) <= GEOM_TOL_MM,
            "virtual segment requires more than 10 um of supplier fit"
        );
        ensure!(
            distance(self.point(self.angle(snapped)), snapped) <= GEOM_TOL_MM,
            "virtual line/cylinder intersection falls outside the source plane"
        );
        Ok(snapped)
    }

    fn trimmed(self, first: [f64; 3], second: [f64; 3]) -> Result<KernelPatchCurve> {
        let a = self.angle(first);
        let b = self.angle(second);
        ensure!(
            distance(self.point(a), first) <= GEOM_TOL_MM
                && distance(self.point(b), second) <= GEOM_TOL_MM,
            "new conic endpoints are not on exact plane/cylinder intersection"
        );
        let span =
            (b - a + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI;
        ensure!(
            span.abs() > DEFAULT_STEP_RECOVERY_TOLERANCES.source_rounding_quarter_turn_slack_rad
                && span.abs() < std::f64::consts::PI,
            "retained source conic arc has ambiguous or degenerate winding"
        );
        Ok(KernelPatchCurve::EllipseArc {
            center_mm: self.center,
            cosine_mm: self.cosine,
            sine_mm: self.sine,
            start_angle_rad: a,
            end_angle_rad: a + span,
        })
    }
}

struct CurvedNeighbor {
    edge: OrientedEdgeUse,
    vertex: u64,
    original: [f64; 3],
    conic: CylinderPlaneConic,
}

fn curved_neighbor(
    solid: &ReducedSolid,
    planar_carrier: &ReducedFace,
    edge: OrientedEdgeUse,
    touched_at_end: bool,
) -> Result<CurvedNeighbor> {
    let SurfaceSupport::Plane(plane) = planar_carrier.surface else {
        bail!(
            "curved segment carrier #{} is not planar",
            planar_carrier.id
        );
    };
    let CurveSupport::BSpline(ref spline) = edge.support else {
        bail!(
            "retained edge #{} is not a spline approximation",
            edge.edge_id
        );
    };
    let neighbor_id = edge_face_users(solid, edge.edge_id)
        .into_iter()
        .find(|&id| id != planar_carrier.id)
        .context("retained spline has no other source face")?;
    let SurfaceSupport::Cylinder(cylinder) = face(solid, neighbor_id)?.surface else {
        bail!("retained spline is not shared by a cylindrical carrier");
    };
    let conic = CylinderPlaneConic::new(plane, cylinder)?;
    conic.fitted_to_source(spline, plane)?;
    let (global_vertices, global_points) = global_edge_endpoints(&edge);
    // The original trimmed B-spline must describe this entire source edge.
    let first = crate::boundary_proof::evaluate_bspline(
        spline,
        crate::boundary_proof::bspline_domain(spline)
            .context("retained spline has no domain")?
            .0,
    )
    .context("invalid spline start")?;
    let last = crate::boundary_proof::evaluate_bspline(
        spline,
        crate::boundary_proof::bspline_domain(spline)
            .context("retained spline has no domain")?
            .1,
    )
    .context("invalid spline end")?;
    ensure!(
        (distance(first, global_points[0]) <= FIT_MM && distance(last, global_points[1]) <= FIT_MM)
            || (distance(first, global_points[1]) <= FIT_MM
                && distance(last, global_points[0]) <= FIT_MM),
        "retained curve is subtrimmed: cannot replace full source spline"
    );
    let (vertex, original) = if touched_at_end {
        (edge.end_vertex, edge.end_mm)
    } else {
        (edge.start_vertex, edge.start_mm)
    };
    ensure!(
        global_vertices.contains(&vertex),
        "source spline lost endpoint identity"
    );
    Ok(CurvedNeighbor {
        edge,
        vertex,
        original,
        conic,
    })
}

struct CarrierConicPlan {
    face_id: u64,
    shared: BTreeSet<u64>,
    previous: CurvedNeighbor,
    next: CurvedNeighbor,
    previous_target: usize,
    next_target: usize,
    snapped: [[f64; 3]; 2],
}

fn effective_carrier_neighbors(
    face: &ReducedFace,
    shared: &BTreeSet<u64>,
) -> Result<(OrientedEdgeUse, OrientedEdgeUse)> {
    // The native STEP importer traverses each FACE_BOUND relative to the
    // enclosing ADVANCED_FACE. Raw EDGE_LOOP order is not the kernel wire
    // order when the bound and surface senses differ.
    let mut oriented = face.clone();
    for wire in &mut oriented.loops {
        if face.surface_forward != wire.orientation {
            wire.edges.reverse();
            for edge in &mut wire.edges {
                std::mem::swap(&mut edge.start_vertex, &mut edge.end_vertex);
                std::mem::swap(&mut edge.start_mm, &mut edge.end_mm);
                edge.parameter_forward = !edge.parameter_forward;
            }
            wire.orientation = face.surface_forward;
        }
    }
    carrier_run_neighbors(&oriented, shared)
}

fn plan_carrier(
    solid: &ReducedSolid,
    carrier: &crate::edge_treatments::EdgeTreatmentCarrier,
    virtual_ends: [[f64; 3]; 2],
    direction: [f64; 3],
) -> Result<CarrierConicPlan> {
    let face = face(solid, carrier.face_id)?;
    let shared = carrier
        .shared_edge_ids
        .iter()
        .copied()
        .collect::<BTreeSet<_>>();
    let (prev, next) = effective_carrier_neighbors(face, &shared)?;
    let previous = curved_neighbor(solid, face, prev, true)?;
    let next = curved_neighbor(solid, face, next, false)?;
    let first_order =
        distance(previous.original, virtual_ends[0]) + distance(next.original, virtual_ends[1]);
    let reverse_order =
        distance(previous.original, virtual_ends[1]) + distance(next.original, virtual_ends[0]);
    let (previous_target, next_target) = if first_order <= reverse_order {
        (0, 1)
    } else {
        (1, 0)
    };
    let mut snapped = [[0.0; 3]; 2];
    snapped[previous_target] = previous
        .conic
        .extend_to(virtual_ends[previous_target], direction)?;
    snapped[next_target] = next.conic.extend_to(virtual_ends[next_target], direction)?;
    Ok(CarrierConicPlan {
        face_id: face.id,
        shared,
        previous,
        next,
        previous_target,
        next_target,
        snapped,
    })
}

fn curve_edit(
    neighbor: CurvedNeighbor,
    carrier_face_id: u64,
    new_vertex: usize,
    snapped: [f64; 3],
) -> Result<EndpointEdit> {
    let (verts, points) = global_edge_endpoints(&neighbor.edge);
    let (front, back) = if verts[0] == neighbor.vertex {
        (snapped, points[1])
    } else if verts[1] == neighbor.vertex {
        (points[0], snapped)
    } else {
        bail!("curved retained edge has no touched endpoint");
    };
    let replacement_curve = Some(neighbor.conic.trimmed(front, back)?);
    Ok(EndpointEdit {
        carrier_face_id,
        edge: neighbor.edge,
        old_vertex: neighbor.vertex,
        new_vertex,
        new_point: snapped,
        replacement_curve,
    })
}

/// Plan a valid exact conic reconstruction when the retained neighbors are
/// B-spline approximations to existing plane/cylinder intersection curves.
/// The original per-face source `StableId`s are preserved.
pub(super) fn plan_curved_segment_patch(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    virtual_start: [f64; 3],
    virtual_end: [f64; 3],
) -> Result<Option<KernelTopologyPatch>> {
    let mut has_spline = false;
    for carrier in &evidence.carriers {
        let planar_carrier = face(solid, carrier.face_id)?;
        let shared = carrier
            .shared_edge_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let (prev, next) = effective_carrier_neighbors(planar_carrier, &shared)?;
        has_spline |= matches!(prev.support, CurveSupport::BSpline(_))
            || matches!(next.support, CurveSupport::BSpline(_));
    }
    if !has_spline {
        return Ok(None);
    }

    let virtual_ends = [virtual_start, virtual_end];
    let direction = normalize(sub(virtual_end, virtual_start), DIR_TOL)
        .context("degenerate candidate virtual segment")?;
    let plans = evidence
        .carriers
        .iter()
        .map(|carrier| plan_carrier(solid, carrier, virtual_ends, direction))
        .collect::<Result<Vec<_>>>()?;
    let [first, second] = plans.as_slice() else {
        bail!("partial segment inverse does not have exactly two planar carriers");
    };
    let snapped = first.snapped;
    ensure!(
        (0..2).all(|i| distance(snapped[i], second.snapped[i]) <= GEOM_TOL_MM),
        "opposing plane/cylinder reconstructions do not agree on one sharp edge"
    );
    let treatment = face(solid, evidence.treatment_face_id)?;
    let treatment_edges = unique_face_edges(treatment)?;
    let mut interface_edges = BTreeSet::new();
    let mut face_splices = Vec::with_capacity(4);
    let mut edits = Vec::with_capacity(4);
    for plan in plans {
        ensure!(
            interface_edges.is_disjoint(&plan.shared),
            "duplicate carrier source interface"
        );
        interface_edges.extend(plan.shared.iter().copied());
        face_splices.push(KernelFaceSplice {
            face: topology_face(provenance, data_section, plan.face_id)?,
            remove_edges: topology_edges(provenance, data_section, plan.shared)?,
            replacement: vec![KernelPatchEdgeUse {
                edge: KernelPatchEdgeRef::New(0),
                orientation: plan.previous_target == 0,
            }],
        });
        edits.push(curve_edit(
            plan.previous,
            plan.face_id,
            plan.previous_target,
            snapped[plan.previous_target],
        )?);
        edits.push(curve_edit(
            plan.next,
            plan.face_id,
            plan.next_target,
            snapped[plan.next_target],
        )?);
    }
    ensure!(
        interface_edges.is_subset(&treatment_edges),
        "curved carrier claims an edge missing from the treatment"
    );
    let terminal_edges = treatment_edges
        .difference(&interface_edges)
        .copied()
        .collect::<BTreeSet<_>>();
    face_splices.extend(plan_terminal_splices(
        solid,
        evidence,
        data_section,
        provenance,
        &terminal_edges,
    )?);
    let rewrites = normalize_retained_edge_rewrites(&edits, data_section, provenance)?;
    Ok(Some(KernelTopologyPatch {
        remove_faces: vec![topology_face(
            provenance,
            data_section,
            evidence.treatment_face_id,
        )?],
        new_vertices: snapped.to_vec(),
        new_edges: vec![KernelNewEdge {
            vertices: (KernelPatchVertexRef::New(0), KernelPatchVertexRef::New(1)),
            curve: KernelPatchCurve::Line {
                start_mm: snapped[0],
                end_mm: snapped[1],
            },
        }],
        rewrite_edges: rewrites,
        face_splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(provenance, data_section, treatment_edges)?,
    }))
}
