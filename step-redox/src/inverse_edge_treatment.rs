//! Pure normalization of proven STEP edge-treatment evidence into kernel topology patches.
//!
//! This module is deliberately read-only with respect to [`ReducedSolid`].  It uses the
//! parser-side B-rep only to understand adjacency and bounded source geometry, then emits a
//! backend-neutral patch addressed entirely by the kernel `StableId` journal.  Mutation belongs to
//! [`crate::kernel_topology::SourceTopologyEditor`].

use crate::brep::{CurveSupport, FaceLoop, OrientedEdgeUse, SurfaceSupport};
use crate::edge_treatments::{EdgeTreatmentEvidence, EdgeTreatmentVirtualEdge};
use crate::kernel_source_load::{KernelTopologyId, KernelTopologyProvenance};
use crate::kernel_topology::{
    KernelFaceSplice, KernelNewEdge, KernelPatchCurve, KernelPatchEdgeRef, KernelPatchEdgeUse,
    KernelPatchVertexRef, KernelRetainedEdge, KernelTopologyPatch,
};
use crate::math3::{
    add, cross, distance, dot, mul, norm, normalize, point_to_unit_line_distance, sub,
};
use crate::reduced_brep::{ReducedFace, ReducedSolid};
use anyhow::{Context, Result, bail, ensure};
use std::collections::{BTreeMap, BTreeSet};

const GEOM_TOL_MM: f64 = crate::numeric::NUMERICAL_GEOMETRY_TOL_MM;
const DIR_TOL: f64 = 1.0e-10;

mod curved_segment;

#[derive(Debug, Clone)]
struct EndpointEdit {
    carrier_face_id: u64,
    edge: OrientedEdgeUse,
    old_vertex: u64,
    new_vertex: usize,
    new_point: [f64; 3],
    replacement_curve: Option<KernelPatchCurve>,
}

#[derive(Debug)]
struct SegmentCarrierPlan {
    splice: KernelFaceSplice,
    endpoint_edits: [EndpointEdit; 2],
}

/// Normalize one proven inverse edge treatment without mutating parser-side topology.
pub fn plan_edge_treatment_patch(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
) -> Result<KernelTopologyPatch> {
    ensure!(
        evidence.solid_id == solid.solid_id(),
        "edge-treatment evidence belongs to solid #{}, not #{}",
        evidence.solid_id,
        solid.solid_id()
    );
    ensure!(
        evidence.carrier_pair_unique,
        "edge-treatment face #{} has no unique carrier pair",
        evidence.treatment_face_id
    );

    match evidence.virtual_edge {
        EdgeTreatmentVirtualEdge::Segment { start_mm, end_mm } => {
            plan_segment_patch(solid, evidence, data_section, provenance, start_mm, end_mm)
        }
        EdgeTreatmentVirtualEdge::Circle {
            center_mm,
            normal,
            radius_mm,
        } => plan_circle_patch(
            solid,
            evidence,
            data_section,
            provenance,
            center_mm,
            normal,
            radius_mm,
        ),
    }
}

/// Partition physical treatments into connected families with the same exact
/// inverse circle. STEP exporters commonly split one fillet around its azimuth:
/// each `ADVANCED_FACE` is then only half of the actual inverse feature.
pub fn edge_treatment_families<'a>(
    solid: &ReducedSolid,
    evidence: &[&'a EdgeTreatmentEvidence],
) -> Vec<Vec<&'a EdgeTreatmentEvidence>> {
    let boundaries = evidence
        .iter()
        .map(|item| {
            face(solid, item.treatment_face_id)
                .ok()
                .map_or_else(BTreeSet::new, |face| {
                    face.loop_edges()
                        .map(|edge| edge.edge_id)
                        .collect::<BTreeSet<_>>()
                })
        })
        .collect::<Vec<_>>();
    let mut seen = vec![false; evidence.len()];
    let mut families = Vec::new();
    for index in 0..evidence.len() {
        if seen[index] {
            continue;
        }
        seen[index] = true;
        let mut stack = vec![index];
        let mut group = Vec::new();
        while let Some(current) = stack.pop() {
            group.push(evidence[current]);
            for next in 0..evidence.len() {
                if seen[next]
                    || evidence[current].kind != evidence[next].kind
                    || !same_virtual_circle(
                        &evidence[current].virtual_edge,
                        &evidence[next].virtual_edge,
                    )
                    || boundaries[current].is_disjoint(&boundaries[next])
                {
                    continue;
                }
                seen[next] = true;
                stack.push(next);
            }
        }
        group.sort_by_key(|item| item.treatment_face_id);
        families.push(group);
    }
    families
}

fn same_virtual_circle(
    first: &EdgeTreatmentVirtualEdge,
    second: &EdgeTreatmentVirtualEdge,
) -> bool {
    let (
        EdgeTreatmentVirtualEdge::Circle {
            center_mm: first_center,
            normal: first_normal,
            radius_mm: first_radius,
        },
        EdgeTreatmentVirtualEdge::Circle {
            center_mm: second_center,
            normal: second_normal,
            radius_mm: second_radius,
        },
    ) = (first, second)
    else {
        return false;
    };
    distance(*first_center, *second_center) <= GEOM_TOL_MM
        && (first_radius - second_radius).abs() <= GEOM_TOL_MM
        && normalize(*first_normal, DIR_TOL)
            .zip(normalize(*second_normal, DIR_TOL))
            .is_some_and(|(a, b)| dot(a, b).abs() >= 1.0 - DIR_TOL)
}

#[derive(Debug)]
struct CircleCarrierRun {
    face_id: u64,
    shell_forward: bool,
    shared: BTreeSet<u64>,
    first: OrientedEdgeUse,
    last: OrientedEdgeUse,
    previous: Option<OrientedEdgeUse>,
    next: Option<OrientedEdgeUse>,
    start_angle: f64,
    sweep: f64,
}

#[derive(Clone, Copy)]
struct CircleFrame {
    center: [f64; 3],
    normal: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
    radius: f64,
}

impl CircleFrame {
    fn from_edge(edge: &EdgeTreatmentVirtualEdge) -> Result<Self> {
        let EdgeTreatmentVirtualEdge::Circle {
            center_mm,
            normal,
            radius_mm,
        } = edge
        else {
            bail!("mixed linear and circular treatment family");
        };
        ensure!(
            radius_mm.is_finite() && *radius_mm > GEOM_TOL_MM,
            "invalid inverse circle radius"
        );
        let normal = normalize(*normal, DIR_TOL).context("invalid inverse circle normal")?;
        let x = deterministic_perpendicular(normal).context("inverse circle has no radial axis")?;
        Ok(Self {
            center: *center_mm,
            normal,
            x,
            y: cross(normal, x),
            radius: *radius_mm,
        })
    }

    fn angle(self, point: [f64; 3]) -> Result<f64> {
        let relative = sub(point, self.center);
        let radial = sub(relative, mul(self.normal, dot(relative, self.normal)));
        let length = norm(radial);
        ensure!(
            length.is_finite() && length > GEOM_TOL_MM,
            "circular edge endpoint is on the inverse axis"
        );
        Ok(dot(radial, self.y)
            .atan2(dot(radial, self.x))
            .rem_euclid(std::f64::consts::TAU))
    }

    fn point(self, angle: f64) -> [f64; 3] {
        add(
            self.center,
            mul(
                add(mul(self.x, angle.cos()), mul(self.y, angle.sin())),
                self.radius,
            ),
        )
    }
}

/// STEP's raw EDGE_LOOP traversal is not necessarily the imported kernel's
/// face-local traversal. In particular ADVANCED_FACE(...,.F.) reverses a
/// FACE_OUTER_BOUND(...,.T.) wire. Both line- and circle-based inverse
/// treatment planners must use the same oriented neighbor order.
fn native_oriented_face_edges(face: &ReducedFace, wire: &FaceLoop) -> Vec<OrientedEdgeUse> {
    if face.surface_forward == wire.orientation {
        return wire.edges.clone();
    }
    wire.edges
        .iter()
        .rev()
        .map(|edge| {
            let mut edge = edge.clone();
            std::mem::swap(&mut edge.start_vertex, &mut edge.end_vertex);
            std::mem::swap(&mut edge.start_mm, &mut edge.end_mm);
            edge.parameter_forward = !edge.parameter_forward;
            edge
        })
        .collect()
}

fn circle_carrier_run(
    face: &ReducedFace,
    shared: BTreeSet<u64>,
    frame: CircleFrame,
) -> Result<CircleCarrierRun> {
    let matching = face
        .loops
        .iter()
        .filter(|wire| wire.edges.iter().any(|edge| shared.contains(&edge.edge_id)))
        .collect::<Vec<_>>();
    let [wire] = matching.as_slice() else {
        bail!(
            "circular carrier face #{} has more than one treatment boundary wire",
            face.id
        );
    };
    ensure!(
        !wire.edges.is_empty() && shared.len() <= wire.edges.len(),
        "circular carrier has invalid run"
    );
    let oriented_edges = native_oriented_face_edges(face, wire);
    let count = oriented_edges.len();
    let starts = (0..count)
        .filter(|&i| {
            shared.contains(&oriented_edges[i].edge_id)
                && (shared.len() == count
                    || !shared.contains(&oriented_edges[(i + count - 1) % count].edge_id))
        })
        .collect::<Vec<_>>();
    let start = if shared.len() == count {
        0
    } else {
        let [start] = starts.as_slice() else {
            bail!(
                "circular carrier face #{} has disconnected treatment run",
                face.id
            );
        };
        *start
    };
    let run = (0..shared.len())
        .map(|i| &oriented_edges[(start + i) % count])
        .collect::<Vec<_>>();
    ensure!(
        run.iter().all(|edge| shared.contains(&edge.edge_id)),
        "circular carrier face #{} has noncontiguous interfaces",
        face.id
    );
    let sweep = circle_interface_sweep(face.id, &run, frame)?;
    let first = (*run[0]).clone();
    let last = (*run[run.len() - 1]).clone();
    let full = shared.len() == count;
    let start_angle = frame.angle(first.start_mm)?;
    let end_angle = frame.angle(last.end_mm)?;
    let expected = if full {
        std::f64::consts::TAU
    } else if sweep > 0.0 {
        (end_angle - start_angle).rem_euclid(std::f64::consts::TAU)
    } else {
        (start_angle - end_angle).rem_euclid(std::f64::consts::TAU)
    };
    ensure!(
        (expected - sweep.abs()).abs() <= 1.0e-6,
        "circular carrier interface angles do not match its exact curve traversal"
    );
    let previous = (!full).then(|| oriented_edges[(start + count - 1) % count].clone());
    let next = (!full).then(|| oriented_edges[(start + run.len()) % count].clone());
    Ok(CircleCarrierRun {
        face_id: face.id,
        shell_forward: face.surface_forward,
        shared,
        first,
        last,
        previous,
        next,
        start_angle,
        sweep,
    })
}

fn circle_interface_sweep(
    face_id: u64,
    run: &[&OrientedEdgeUse],
    frame: CircleFrame,
) -> Result<f64> {
    let mut sweep = 0.0_f64;
    for edge in run {
        let CurveSupport::Circle(circle) = edge.support else {
            bail!(
                "circular carrier face #{face_id} has a non-circular interface edge #{}",
                edge.edge_id
            );
        };
        ensure!(
            normalize(circle.normal, DIR_TOL)
                .is_some_and(|normal| dot(normal, frame.normal).abs() >= 1.0 - DIR_TOL)
                && point_to_unit_line_distance(circle.center_mm, frame.center, frame.normal)
                    <= GEOM_TOL_MM,
            "circular interface #{} does not share the inverse axis",
            edge.edge_id
        );
        let from = frame.angle(edge.start_mm)?;
        let to = frame.angle(edge.end_mm)?;
        let increasing = edge.parameter_forward == (dot(circle.normal, frame.normal) > 0.0);
        let span = if distance(edge.start_mm, edge.end_mm) <= GEOM_TOL_MM {
            std::f64::consts::TAU
        } else if increasing {
            (to - from).rem_euclid(std::f64::consts::TAU)
        } else {
            -(from - to).rem_euclid(std::f64::consts::TAU)
        };
        ensure!(
            sweep == 0.0 || sweep.signum() == span.signum(),
            "circular carrier has opposing interface windings"
        );
        sweep += span;
    }
    ensure!(
        sweep.abs() <= std::f64::consts::TAU + 1.0e-7,
        "circular carrier interface covers more than one turn"
    );
    Ok(sweep)
}

fn circle_carrier_supports_inverse(face: &ReducedFace, frame: CircleFrame) -> bool {
    match &face.surface {
        SurfaceSupport::Plane(plane) => {
            normalize(plane.normal, DIR_TOL)
                .is_some_and(|n| dot(n, frame.normal).abs() >= 1.0 - DIR_TOL)
                && dot(sub(frame.center, plane.origin_mm), plane.normal).abs() <= GEOM_TOL_MM
        }
        SurfaceSupport::Cylinder(cylinder) => {
            normalize(cylinder.axis, DIR_TOL).is_some_and(|axis| {
                dot(axis, frame.normal).abs() >= 1.0 - DIR_TOL
                    && point_to_unit_line_distance(frame.center, cylinder.axis_origin_mm, axis)
                        <= GEOM_TOL_MM
            }) && (cylinder.radius_mm - frame.radius).abs() <= GEOM_TOL_MM
        }
        SurfaceSupport::Cone(cone) => normalize(cone.axis, DIR_TOL).is_some_and(|axis| {
            dot(axis, frame.normal).abs() >= 1.0 - DIR_TOL
                && point_to_unit_line_distance(frame.center, cone.reference_origin_mm, axis)
                    <= GEOM_TOL_MM
                && (cone.semi_angle_rad.tan().mul_add(
                    dot(sub(frame.center, cone.reference_origin_mm), axis),
                    cone.reference_radius_mm,
                ) - frame.radius)
                    .abs()
                    <= GEOM_TOL_MM
        }),
        _ => false,
    }
}

fn inverse_circle_angles(frame: CircleFrame, runs: &[CircleCarrierRun]) -> Result<Vec<f64>> {
    let mut breaks = vec![0.0];
    for run in runs {
        for point in [run.first.start_mm, run.last.end_mm] {
            let angle = frame.angle(point)?;
            if !breaks.iter().any(|a| {
                frame.radius * (a - angle).sin().abs() <= GEOM_TOL_MM
                    && frame.radius * (a - angle).cos() >= frame.radius - GEOM_TOL_MM
            }) {
                breaks.push(angle);
            }
        }
    }
    breaks.sort_by(f64::total_cmp);
    // Any original segment is split further into bounded arcs without
    // interpolating its geometric locus: all pieces use the same exact circle.
    let mut bounded = Vec::new();
    for i in 0..breaks.len() {
        let start = breaks[i];
        let end = breaks
            .get(i + 1)
            .copied()
            .unwrap_or(breaks[0] + std::f64::consts::TAU);
        let span = end - start;
        let pieces = (1_u32..=4)
            .find(|&count| span <= f64::from(count).mul_add(std::f64::consts::FRAC_PI_2, 1.0e-10))
            .context("inverse circle has an invalid angular partition")?;
        for j in 0..pieces {
            bounded.push(
                (f64::from(j) / f64::from(pieces))
                    .mul_add(span, start)
                    .rem_euclid(std::f64::consts::TAU),
            );
        }
    }
    ensure!(
        bounded.len() >= 2 && bounded.len() <= 128,
        "inverse ring exceeds supported exact arc partition size"
    );
    Ok(bounded)
}

fn circle_arc_cycle(frame: CircleFrame, angles: &[f64]) -> (Vec<[f64; 3]>, Vec<KernelNewEdge>) {
    let vertices = angles.iter().map(|a| frame.point(*a)).collect::<Vec<_>>();
    let edges = angles
        .iter()
        .enumerate()
        .map(|(i, a)| KernelNewEdge {
            vertices: (
                KernelPatchVertexRef::New(i),
                KernelPatchVertexRef::New((i + 1) % angles.len()),
            ),
            curve: KernelPatchCurve::CircleArc {
                center_mm: frame.center,
                normal: frame.normal,
                x_direction: frame.x,
                radius_mm: frame.radius,
                start_angle_rad: *a,
                end_angle_rad: angles
                    .get(i + 1)
                    .copied()
                    .unwrap_or(angles[0] + std::f64::consts::TAU),
            },
        })
        .collect();
    (vertices, edges)
}

fn circle_run_replacement(
    run: &CircleCarrierRun,
    frame: CircleFrame,
    angles: &[f64],
    flip_full_loop: bool,
) -> Result<Vec<KernelPatchEdgeUse>> {
    ensure!(
        !flip_full_loop || run.previous.is_none(),
        "cannot reverse a retained partial circular wire"
    );
    let tolerance = GEOM_TOL_MM / frame.radius;
    let start = angles
        .iter()
        .position(|a| {
            ((a - run.start_angle).rem_euclid(std::f64::consts::TAU))
                .min((run.start_angle - a).rem_euclid(std::f64::consts::TAU))
                <= tolerance
        })
        .context("circular interface start has no shared vertex")?;
    let full = run.previous.is_none();
    let end_angle = frame.angle(run.last.end_mm)?;
    let end = if full {
        start
    } else {
        angles
            .iter()
            .position(|a| {
                ((a - end_angle).rem_euclid(std::f64::consts::TAU))
                    .min((end_angle - a).rem_euclid(std::f64::consts::TAU))
                    <= tolerance
            })
            .context("circular interface end has no shared vertex")?
    };
    let mut result = Vec::new();
    let mut at = start;
    loop {
        if result.len() >= angles.len() {
            bail!("circular run would reuse a replacement arc");
        }
        let increasing = (run.sweep > 0.0) != flip_full_loop;
        let next = if increasing {
            (at + 1) % angles.len()
        } else {
            (at + angles.len() - 1) % angles.len()
        };
        let edge = if increasing { at } else { next };
        result.push(KernelPatchEdgeUse {
            edge: KernelPatchEdgeRef::New(edge),
            orientation: increasing,
        });
        at = next;
        if at == end && (!full || result.len() == angles.len()) {
            break;
        }
    }
    let span = result
        .iter()
        .map(|use_| {
            let KernelPatchEdgeRef::New(i) = use_.edge else {
                unreachable!()
            };
            let from = angles[i];
            let to = angles
                .get(i + 1)
                .copied()
                .unwrap_or(angles[0] + std::f64::consts::TAU);
            to - from
        })
        .sum::<f64>();
    ensure!(
        (span - run.sweep.abs()).abs() <= 1.0e-6,
        "reconstructed circular run has wrong total angular span"
    );
    Ok(result)
}

fn circle_line_endpoint_edit(
    carrier_id: u64,
    edge: OrientedEdgeUse,
    old_vertex: u64,
    new_vertex: usize,
    new_point: [f64; 3],
) -> Result<EndpointEdit> {
    let CurveSupport::Line(_) = edge.support else {
        bail!(
            "circular carrier #{carrier_id} has a non-linear retained neighbor #{}",
            edge.edge_id
        );
    };
    let delta = sub(edge.end_mm, edge.start_mm);
    let length = norm(delta);
    ensure!(
        length > GEOM_TOL_MM
            && norm(cross(delta, sub(new_point, edge.start_mm))) / length <= GEOM_TOL_MM,
        "retained edge #{} does not extend onto inverse circle",
        edge.edge_id
    );
    Ok(EndpointEdit {
        carrier_face_id: carrier_id,
        edge,
        old_vertex,
        new_vertex,
        new_point,
        replacement_curve: None,
    })
}

/// A complete circular carrier has no retained endpoints, so its winding
/// can be selected from the native two-sided edge incidence. Partial carrier
/// runs have retained neighbors and therefore keep their original winding.
fn choose_two_sided_circle_runs(
    runs: &[CircleCarrierRun],
    frame: CircleFrame,
    angles: &[f64],
) -> Result<Vec<Vec<KernelPatchEdgeUse>>> {
    let full_runs = runs
        .iter()
        .enumerate()
        .filter_map(|(index, run)| run.previous.is_none().then_some(index))
        .collect::<Vec<_>>();
    ensure!(
        full_runs.len() <= 2,
        "inverse circle has more than two complete carriers"
    );
    for mask in 0..(1_usize << full_runs.len()) {
        let replacements = runs
            .iter()
            .enumerate()
            .map(|(index, run)| {
                let flip = full_runs
                    .iter()
                    .enumerate()
                    .any(|(bit, &selected)| selected == index && mask & (1 << bit) != 0);
                circle_run_replacement(run, frame, angles, flip)
            })
            .collect::<Result<Vec<_>>>()?;
        let mut uses = vec![[0_usize; 2]; angles.len()];
        for (run, replacement) in runs.iter().zip(&replacements) {
            for edge_use in replacement {
                let KernelPatchEdgeRef::New(index) = edge_use.edge else {
                    unreachable!()
                };
                // The runtime shell inverts edge use when its face
                // orientation is reversed. Count *effective* shell uses,
                // not the face-local wire directions.
                uses[index][usize::from(edge_use.orientation == run.shell_forward)] += 1;
            }
        }
        if uses.iter().all(|uses| *uses == [1, 1]) {
            return Ok(replacements);
        }
    }
    bail!("inverse ring is not a complete two-sided shared circular interface")
}

fn circle_family_carrier_patch(
    runs: &[CircleCarrierRun],
    frame: CircleFrame,
    angles: &[f64],
    data_section: usize,
    provenance: &KernelTopologyProvenance,
) -> Result<(Vec<KernelFaceSplice>, Vec<KernelRetainedEdge>)> {
    let mut splices = Vec::new();
    let mut edits = Vec::new();
    let replacements = choose_two_sided_circle_runs(runs, frame, angles)?;
    for (run, replacement) in runs.iter().zip(replacements) {
        splices.push(KernelFaceSplice {
            face: topology_face(provenance, data_section, run.face_id)?,
            remove_edges: topology_edges(provenance, data_section, run.shared.iter().copied())?,
            replacement,
        });
        if let (Some(previous), Some(next)) = (&run.previous, &run.next) {
            let start_idx = angles
                .iter()
                .position(|a| {
                    distance(frame.point(*a), frame.point(run.start_angle)) <= GEOM_TOL_MM
                })
                .context("missing inverse run start")?;
            let end_angle = frame.angle(run.last.end_mm)?;
            let end_idx = angles
                .iter()
                .position(|a| distance(frame.point(*a), frame.point(end_angle)) <= GEOM_TOL_MM)
                .context("missing inverse run end")?;
            edits.push(circle_line_endpoint_edit(
                run.face_id,
                previous.clone(),
                previous.end_vertex,
                start_idx,
                frame.point(angles[start_idx]),
            )?);
            edits.push(circle_line_endpoint_edit(
                run.face_id,
                next.clone(),
                next.start_vertex,
                end_idx,
                frame.point(angles[end_idx]),
            )?);
        }
    }
    Ok((
        splices,
        normalize_retained_edge_rewrites(&edits, data_section, provenance)?,
    ))
}

/// Replace a connected, full-azimuth family of partial toroidal/conical
/// treatment faces as one exact native topology transaction.
pub fn plan_circle_family_patch(
    solid: &ReducedSolid,
    family: &[&EdgeTreatmentEvidence],
    data_section: usize,
    provenance: &KernelTopologyProvenance,
) -> Result<KernelTopologyPatch> {
    ensure!(
        family.len() > 1,
        "partial circular family requires multiple source treatment faces"
    );
    let frame = CircleFrame::from_edge(&family[0].virtual_edge)?;
    let mut removed = BTreeSet::new();
    let mut treatment_edges = BTreeSet::new();
    let mut by_carrier = BTreeMap::<u64, BTreeSet<u64>>::new();
    for member in family {
        ensure!(
            member.carrier_pair_unique
                && same_virtual_circle(&member.virtual_edge, &family[0].virtual_edge)
                && member.kind == family[0].kind,
            "inconsistent circular treatment family"
        );
        ensure!(
            removed.insert(member.treatment_face_id),
            "duplicated circular treatment face"
        );
        treatment_edges.extend(unique_face_edges(face(solid, member.treatment_face_id)?)?);
        for carrier in &member.carriers {
            let entry = by_carrier.entry(carrier.face_id).or_default();
            for &edge in &carrier.shared_edge_ids {
                ensure!(
                    entry.insert(edge),
                    "circular treatment family assigns shared interface edge more than once"
                );
            }
        }
    }
    ensure!(
        by_carrier.keys().all(|id| !removed.contains(id)),
        "circular treatment overlaps its own carrier"
    );
    let interfaces = by_carrier
        .values()
        .flat_map(|edges| edges.iter().copied())
        .collect::<BTreeSet<_>>();
    ensure!(
        interfaces.is_subset(&treatment_edges),
        "circular family interface is not a treatment boundary"
    );
    for &edge in &treatment_edges {
        let users = edge_face_users(solid, edge);
        ensure!(
            users.len() == 2 && users.iter().any(|id| removed.contains(id)),
            "circular treatment edge #{edge} is not manifold in its source solid"
        );
        if interfaces.contains(&edge) {
            ensure!(
                users.iter().filter(|id| removed.contains(id)).count() == 1
                    && users
                        .iter()
                        .any(|id| by_carrier.get(id).is_some_and(|set| set.contains(&edge))),
                "circular treatment edge #{edge} has an invalid surviving carrier"
            );
        } else {
            ensure!(
                users.iter().all(|id| removed.contains(id)),
                "circular treatment family has a real terminal boundary at edge #{edge}"
            );
        }
    }
    let mut runs = Vec::new();
    for (face_id, shared) in by_carrier {
        let carrier = face(solid, face_id)?;
        ensure!(
            circle_carrier_supports_inverse(carrier, frame),
            "circular inverse circle is not on carrier face #{face_id}"
        );
        runs.push(circle_carrier_run(carrier, shared, frame)?);
    }
    let angles = inverse_circle_angles(frame, &runs)?;
    let (new_vertices, new_edges) = circle_arc_cycle(frame, &angles);
    let (face_splices, rewrite_edges) =
        circle_family_carrier_patch(&runs, frame, &angles, data_section, provenance)?;
    Ok(KernelTopologyPatch {
        remove_faces: removed
            .into_iter()
            .map(|id| topology_face(provenance, data_section, id))
            .collect::<Result<Vec<_>>>()?,
        new_vertices,
        new_edges,
        rewrite_edges,
        face_splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(provenance, data_section, treatment_edges)?,
    })
}

fn plan_segment_patch(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    virtual_start: [f64; 3],
    virtual_end: [f64; 3],
) -> Result<KernelTopologyPatch> {
    ensure!(
        distance(virtual_start, virtual_end) > GEOM_TOL_MM,
        "inverse segment for treatment face #{} is degenerate",
        evidence.treatment_face_id
    );
    if let Some(patch) = curved_segment::plan_curved_segment_patch(
        solid,
        evidence,
        data_section,
        provenance,
        virtual_start,
        virtual_end,
    )? {
        return Ok(patch);
    }
    let treatment = face(solid, evidence.treatment_face_id)?;
    let treatment_edges = unique_face_edges(treatment)?;
    let treatment_face = topology_face(provenance, data_section, evidence.treatment_face_id)?;
    let (carrier_plans, interface_edges) = collect_segment_carrier_plans(
        solid,
        evidence,
        data_section,
        provenance,
        virtual_start,
        virtual_end,
    )?;
    ensure!(
        interface_edges.is_subset(&treatment_edges),
        "carrier interface contains an edge absent from treatment face #{}",
        evidence.treatment_face_id
    );
    let endpoint_edits = carrier_plans
        .iter()
        .flat_map(|plan| plan.endpoint_edits.iter().cloned())
        .collect::<Vec<_>>();
    let rewrite_edges =
        normalize_retained_edge_rewrites(&endpoint_edits, data_section, provenance)?;
    let terminal_edges = treatment_edges
        .difference(&interface_edges)
        .copied()
        .collect::<BTreeSet<_>>();
    ensure!(
        !terminal_edges.is_empty(),
        "segment treatment face #{} has no terminal boundary",
        evidence.treatment_face_id
    );
    let mut face_splices = carrier_plans
        .into_iter()
        .map(|plan| plan.splice)
        .collect::<Vec<_>>();
    face_splices.extend(plan_terminal_splices(
        solid,
        evidence,
        data_section,
        provenance,
        &terminal_edges,
    )?);

    Ok(KernelTopologyPatch {
        remove_faces: vec![treatment_face],
        new_vertices: vec![virtual_start, virtual_end],
        new_edges: vec![KernelNewEdge {
            vertices: (KernelPatchVertexRef::New(0), KernelPatchVertexRef::New(1)),
            curve: KernelPatchCurve::Line {
                start_mm: virtual_start,
                end_mm: virtual_end,
            },
        }],
        rewrite_edges,
        face_splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(
            provenance,
            data_section,
            treatment_edges.iter().copied(),
        )?,
    })
}

fn collect_segment_carrier_plans(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    virtual_start: [f64; 3],
    virtual_end: [f64; 3],
) -> Result<(Vec<SegmentCarrierPlan>, BTreeSet<u64>)> {
    let mut plans = Vec::with_capacity(2);
    let mut interface_edges = BTreeSet::new();
    for carrier in &evidence.carriers {
        for &edge_id in &carrier.shared_edge_ids {
            ensure!(
                interface_edges.insert(edge_id),
                "treatment face #{} assigns source edge #{edge_id} to both carriers",
                evidence.treatment_face_id
            );
        }
        plans.push(plan_segment_carrier(
            face(solid, carrier.face_id)?,
            &carrier.shared_edge_ids,
            data_section,
            provenance,
            virtual_start,
            virtual_end,
        )?);
    }
    Ok((plans, interface_edges))
}

fn plan_terminal_splices(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    terminal_edges: &BTreeSet<u64>,
) -> Result<Vec<KernelFaceSplice>> {
    let mut terminal_by_face = BTreeMap::<u64, BTreeSet<u64>>::new();
    for &terminal_edge in terminal_edges {
        let users = edge_face_users(solid, terminal_edge);
        ensure!(
            users.len() == 2 && users.contains(&evidence.treatment_face_id),
            "terminal treatment edge #{terminal_edge} is not manifold around face #{}",
            evidence.treatment_face_id
        );
        let terminal_face = users
            .iter()
            .copied()
            .find(|&face_id| face_id != evidence.treatment_face_id)
            .context("terminal treatment edge has no surviving face")?;
        ensure!(
            !evidence
                .carriers
                .iter()
                .any(|carrier| carrier.face_id == terminal_face),
            "terminal treatment edge #{terminal_edge} folds back into a carrier face"
        );
        terminal_by_face
            .entry(terminal_face)
            .or_default()
            .insert(terminal_edge);
    }
    terminal_by_face
        .into_iter()
        .map(|(face_id, edges)| {
            let terminal_face = face(solid, face_id)?;
            ensure!(
                edges_form_one_cyclic_run(terminal_face, &edges),
                "terminal edges {edges:?} are not one contiguous run on face #{face_id}"
            );
            Ok(KernelFaceSplice {
                face: topology_face(provenance, data_section, face_id)?,
                remove_edges: topology_edges(provenance, data_section, edges.iter().copied())?,
                replacement: Vec::new(),
            })
        })
        .collect()
}

fn plan_segment_carrier(
    face: &ReducedFace,
    shared_edge_ids: &[u64],
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    virtual_start: [f64; 3],
    virtual_end: [f64; 3],
) -> Result<SegmentCarrierPlan> {
    ensure!(
        !shared_edge_ids.is_empty(),
        "carrier has no shared treatment edges"
    );
    let shared = shared_edge_ids.iter().copied().collect::<BTreeSet<_>>();
    ensure!(
        shared.len() == shared_edge_ids.len(),
        "carrier face #{} repeats a shared treatment edge",
        face.id
    );
    let (previous, next) = carrier_run_neighbors(face, &shared)?;
    let CurveSupport::Line(previous_line) = previous.support else {
        bail!(
            "carrier face #{} previous retained edge is not linear",
            face.id
        );
    };
    let CurveSupport::Line(next_line) = next.support else {
        bail!("carrier face #{} next retained edge is not linear", face.id);
    };
    let previous_direction = normalize(previous_line.direction, DIR_TOL)
        .context("previous retained edge direction is degenerate")?;
    let next_direction = normalize(next_line.direction, DIR_TOL)
        .context("next retained edge direction is degenerate")?;
    let (healed_start, healed_end) = [(virtual_start, virtual_end), (virtual_end, virtual_start)]
        .into_iter()
        .find(|(previous_point, next_point)| {
            point_to_unit_line_distance(
                *previous_point,
                previous_line.origin_mm,
                previous_direction,
            ) <= GEOM_TOL_MM
                && point_to_unit_line_distance(*next_point, next_line.origin_mm, next_direction)
                    <= GEOM_TOL_MM
        })
        .with_context(|| {
            format!(
                "proven inverse segment does not meet retained neighbors on carrier face #{}",
                face.id
            )
        })?;
    let forward = distance(healed_start, virtual_start) <= GEOM_TOL_MM;
    let (start_vertex, end_vertex) = if forward { (0, 1) } else { (1, 0) };

    Ok(SegmentCarrierPlan {
        splice: KernelFaceSplice {
            face: topology_face(provenance, data_section, face.id)?,
            remove_edges: topology_edges(provenance, data_section, shared.iter().copied())?,
            replacement: vec![KernelPatchEdgeUse {
                edge: KernelPatchEdgeRef::New(0),
                orientation: forward,
            }],
        },
        endpoint_edits: [
            EndpointEdit {
                carrier_face_id: face.id,
                old_vertex: previous.end_vertex,
                edge: previous,
                new_vertex: start_vertex,
                new_point: healed_start,
                replacement_curve: None,
            },
            EndpointEdit {
                carrier_face_id: face.id,
                old_vertex: next.start_vertex,
                edge: next,
                new_vertex: end_vertex,
                new_point: healed_end,
                replacement_curve: None,
            },
        ],
    })
}

fn carrier_run_neighbors(
    face: &ReducedFace,
    shared: &BTreeSet<u64>,
) -> Result<(OrientedEdgeUse, OrientedEdgeUse)> {
    let matches = face
        .loops
        .iter()
        .enumerate()
        .filter_map(|(loop_index, loop_)| {
            let count = loop_
                .edges
                .iter()
                .filter(|edge| shared.contains(&edge.edge_id))
                .count();
            (count != 0).then_some((loop_index, count))
        })
        .collect::<Vec<_>>();
    let [(loop_index, count)] = matches.as_slice() else {
        bail!(
            "carrier face #{} does not have one treatment boundary run",
            face.id
        );
    };
    ensure!(
        *count == shared.len(),
        "carrier face #{} does not contain every shared treatment edge",
        face.id
    );
    let loop_ = &face.loops[*loop_index];
    let oriented_edges = native_oriented_face_edges(face, loop_);
    let edge_count = oriented_edges.len();
    ensure!(
        shared.len() < edge_count.saturating_sub(1),
        "carrier face #{} leaves no retained neighbors around treatment run",
        face.id
    );
    let starts = (0..edge_count)
        .filter(|&index| {
            shared.contains(&oriented_edges[index].edge_id)
                && !shared.contains(&oriented_edges[(index + edge_count - 1) % edge_count].edge_id)
        })
        .collect::<Vec<_>>();
    let [run_start] = starts.as_slice() else {
        bail!(
            "carrier face #{} treatment edges are not one cyclic run",
            face.id
        );
    };
    let mut run_len = 0usize;
    while run_len < edge_count
        && shared.contains(&oriented_edges[(run_start + run_len) % edge_count].edge_id)
    {
        run_len += 1;
    }
    ensure!(
        run_len == shared.len(),
        "carrier face #{} treatment run is disconnected",
        face.id
    );
    Ok((
        oriented_edges[(run_start + edge_count - 1) % edge_count].clone(),
        oriented_edges[(run_start + run_len) % edge_count].clone(),
    ))
}

fn normalize_retained_edge_rewrites(
    edits: &[EndpointEdit],
    data_section: usize,
    provenance: &KernelTopologyProvenance,
) -> Result<Vec<KernelRetainedEdge>> {
    let mut by_edge = BTreeMap::<u64, Vec<&EndpointEdit>>::new();
    for edit in edits {
        by_edge.entry(edit.edge.edge_id).or_default().push(edit);
    }
    let mut out = Vec::with_capacity(by_edge.len());
    for (edge_id, edits) in by_edge {
        let representative = edits
            .first()
            .context("retained edge rewrite lost its representative use")?;
        let (vertices, mut points) = global_edge_endpoints(&representative.edge);
        let curve_replacement = edits
            .iter()
            .find_map(|edit| edit.replacement_curve.as_ref());
        match &representative.edge.support {
            CurveSupport::Line(_) => ensure!(
                curve_replacement.is_none(),
                "retained line #{edge_id} cannot use a curved replacement"
            ),
            CurveSupport::BSpline(_) => ensure!(
                curve_replacement.is_some()
                    && edits
                        .iter()
                        .all(|edit| edit.replacement_curve.as_ref() == curve_replacement),
                "retained spline #{edge_id} lacks a consistent analytic replacement"
            ),
            _ => bail!("retained edge #{edge_id} is not reconstructible"),
        }
        let mut changed = [false; 2];
        let mut new_refs = [None, None];
        for edit in edits {
            let endpoint = if vertices[0] == edit.old_vertex {
                0
            } else if vertices[1] == edit.old_vertex {
                1
            } else {
                bail!(
                    "retained edge #{edge_id} on carrier face #{} does not contain edited vertex #{}",
                    edit.carrier_face_id,
                    edit.old_vertex
                );
            };
            if changed[endpoint] {
                ensure!(
                    new_refs[endpoint] == Some(edit.new_vertex)
                        && distance(points[endpoint], edit.new_point) <= GEOM_TOL_MM,
                    "retained edge #{edge_id} receives conflicting endpoint edits"
                );
                continue;
            }
            changed[endpoint] = true;
            new_refs[endpoint] = Some(edit.new_vertex);
            points[endpoint] = edit.new_point;
        }
        let refs = [0usize, 1usize].map(|endpoint| {
            new_refs[endpoint].map_or_else(
                || {
                    topology_vertex(provenance, data_section, vertices[endpoint])
                        .map(KernelPatchVertexRef::Existing)
                },
                |new_index| Ok(KernelPatchVertexRef::New(new_index)),
            )
        });
        let [front, back] = refs;
        ensure!(
            distance(points[0], points[1]) > GEOM_TOL_MM,
            "retained edge #{edge_id} collapses after inverse treatment"
        );
        out.push(KernelRetainedEdge {
            edge: topology_edge(provenance, data_section, edge_id)?,
            vertices: (front?, back?),
            curve: curve_replacement
                .cloned()
                .unwrap_or(KernelPatchCurve::Line {
                    start_mm: points[0],
                    end_mm: points[1],
                }),
        });
    }
    Ok(out)
}

const fn global_edge_endpoints(edge: &OrientedEdgeUse) -> ([u64; 2], [[f64; 3]; 2]) {
    if edge.parameter_forward == edge.curve_same_sense {
        (
            [edge.start_vertex, edge.end_vertex],
            [edge.start_mm, edge.end_mm],
        )
    } else {
        (
            [edge.end_vertex, edge.start_vertex],
            [edge.end_mm, edge.start_mm],
        )
    }
}

/// Continue one retained straight edge to its intersection with the exact
/// inverse circle. The extension must leave the source edge at the endpoint
/// touching the removed treatment, never change its interior locus.
fn extended_line_circle_endpoint(
    edge: &OrientedEdgeUse,
    frame: CircleFrame,
    previous: bool,
) -> Result<[f64; 3]> {
    ensure!(
        matches!(edge.support, CurveSupport::Line(_)),
        "retained circular-treatment neighbor #{} is not straight",
        edge.edge_id
    );
    let direction = sub(edge.end_mm, edge.start_mm);
    let length = norm(direction);
    ensure!(length > GEOM_TOL_MM, "retained neighbor edge is degenerate");
    let from_center = sub(edge.start_mm, frame.center);
    let axial = dot(from_center, frame.normal);
    let axial_direction = dot(direction, frame.normal);
    let mut candidates = Vec::<f64>::new();
    if axial_direction.abs() > DIR_TOL * length {
        candidates.push(-axial / axial_direction);
    } else if axial.abs() <= GEOM_TOL_MM {
        let radial = sub(from_center, mul(frame.normal, axial));
        let radial_direction = sub(direction, mul(frame.normal, axial_direction));
        let a = dot(radial_direction, radial_direction);
        let b = 2.0 * dot(radial, radial_direction);
        let c = frame.radius.mul_add(-frame.radius, dot(radial, radial));
        let discriminant = b.mul_add(b, -4.0 * a * c);
        if a > DIR_TOL * DIR_TOL && discriminant >= 0.0 {
            let root = discriminant.sqrt();
            candidates.extend([(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]);
        }
    }
    let at = if previous { 1.0 } else { 0.0 };
    let parameter = candidates
        .into_iter()
        .filter(|t| {
            t.is_finite()
                && if previous {
                    *t >= at - 1.0e-9
                } else {
                    *t <= at + 1.0e-9
                }
        })
        .min_by(|a, b| (a - at).abs().total_cmp(&(b - at).abs()))
        .context("retained straight neighbor does not extend to inverse circle")?;
    let point = add(edge.start_mm, mul(direction, parameter));
    ensure!(
        dot(sub(point, frame.center), frame.normal).abs() <= GEOM_TOL_MM
            && (norm(sub(point, frame.center)) - frame.radius).abs() <= GEOM_TOL_MM,
        "line/circle extension does not belong to inverse carrier support"
    );
    Ok(point)
}

fn open_inverse_circle_arc(
    frame: CircleFrame,
    first: [f64; 3],
    last: [f64; 3],
    source_increasing: bool,
) -> Result<(Vec<[f64; 3]>, Vec<KernelNewEdge>)> {
    let a = frame.angle(first)?;
    let b = frame.angle(last)?;
    let (start, span) = if source_increasing {
        (a, (b - a).rem_euclid(std::f64::consts::TAU))
    } else {
        (b, (a - b).rem_euclid(std::f64::consts::TAU))
    };
    ensure!(
        span > GEOM_TOL_MM / frame.radius
            && span < std::f64::consts::TAU - GEOM_TOL_MM / frame.radius,
        "partial inverse circle must span a proper open arc"
    );
    let pieces = (1_u32..=4)
        .find(|&n| span <= f64::from(n).mul_add(std::f64::consts::FRAC_PI_2, 1.0e-10))
        .context("partial inverse circle has an invalid arc")?;
    let angles = (0..=pieces)
        .map(|i| f64::from(i).mul_add(span / f64::from(pieces), start))
        .collect::<Vec<_>>();
    let vertices = angles
        .iter()
        .map(|&angle| frame.point(angle))
        .collect::<Vec<_>>();
    let edges = angles
        .windows(2)
        .enumerate()
        .map(|(index, pair)| KernelNewEdge {
            vertices: (
                KernelPatchVertexRef::New(index),
                KernelPatchVertexRef::New(index + 1),
            ),
            curve: KernelPatchCurve::CircleArc {
                center_mm: frame.center,
                normal: frame.normal,
                x_direction: frame.x,
                radius_mm: frame.radius,
                start_angle_rad: pair[0],
                end_angle_rad: pair[1],
            },
        })
        .collect();
    Ok((vertices, edges))
}

fn partial_circle_carrier_splices(
    carrier_runs: &[(CircleCarrierRun, [f64; 3], [f64; 3])],
    new_vertices: &[[f64; 3]],
    edge_count: usize,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
) -> Result<(Vec<KernelFaceSplice>, Vec<EndpointEdit>)> {
    let first = new_vertices
        .first()
        .context("partial circle has no first vertex")?;
    let last = new_vertices
        .last()
        .context("partial circle has no last vertex")?;
    let mut face_splices = Vec::new();
    let mut edits = Vec::new();
    let mut uses = vec![[0_usize; 2]; edge_count];
    for (run, start, end) in carrier_runs {
        let increasing = if distance(*start, *first) <= GEOM_TOL_MM
            && distance(*end, *last) <= GEOM_TOL_MM
        {
            true
        } else if distance(*start, *last) <= GEOM_TOL_MM && distance(*end, *first) <= GEOM_TOL_MM {
            false
        } else {
            bail!("partial circle carrier does not span the exact shared arc");
        };
        let replacements = if increasing {
            (0..edge_count)
                .map(|index| KernelPatchEdgeUse {
                    edge: KernelPatchEdgeRef::New(index),
                    orientation: true,
                })
                .collect::<Vec<_>>()
        } else {
            (0..edge_count)
                .rev()
                .map(|index| KernelPatchEdgeUse {
                    edge: KernelPatchEdgeRef::New(index),
                    orientation: false,
                })
                .collect::<Vec<_>>()
        };
        for edge_use in &replacements {
            let KernelPatchEdgeRef::New(index) = edge_use.edge else {
                unreachable!()
            };
            uses[index][usize::from(edge_use.orientation == run.shell_forward)] += 1;
        }
        face_splices.push(KernelFaceSplice {
            face: topology_face(provenance, data_section, run.face_id)?,
            remove_edges: topology_edges(provenance, data_section, run.shared.iter().copied())?,
            replacement: replacements,
        });
        let previous = run.previous.as_ref().context("partial run lost previous")?;
        let next = run.next.as_ref().context("partial run lost next")?;
        let start_index = if increasing {
            0
        } else {
            new_vertices.len() - 1
        };
        let end_index = if increasing {
            new_vertices.len() - 1
        } else {
            0
        };
        edits.push(circle_line_endpoint_edit(
            run.face_id,
            previous.clone(),
            previous.end_vertex,
            start_index,
            *start,
        )?);
        edits.push(circle_line_endpoint_edit(
            run.face_id,
            next.clone(),
            next.start_vertex,
            end_index,
            *end,
        )?);
    }
    ensure!(
        uses.iter().all(|pair| *pair == [1, 1]),
        "partial inverse arc lacks two oppositely oriented shell uses"
    );
    Ok((face_splices, edits))
}

fn plan_partial_circle_patch(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    frame: CircleFrame,
    treatment_edges: &BTreeSet<u64>,
) -> Result<KernelTopologyPatch> {
    let mut carrier_runs = Vec::new();
    let mut interface_edges = BTreeSet::new();
    for carrier in &evidence.carriers {
        let face = face(solid, carrier.face_id)?;
        ensure!(
            circle_carrier_supports_inverse(face, frame),
            "partial circular carrier #{} does not support inverse circle",
            carrier.face_id
        );
        let shared = carrier
            .shared_edge_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        ensure!(
            !shared.is_empty() && shared.len() == carrier.shared_edge_ids.len(),
            "partial circular treatment carrier has invalid shared edges"
        );
        ensure!(
            shared.iter().all(|id| interface_edges.insert(*id)),
            "partial circular treatment has duplicate carrier edges"
        );
        let run = circle_carrier_run(face, shared, frame)?;
        let (Some(previous), Some(next)) = (&run.previous, &run.next) else {
            bail!("partial inverse circle carrier has no retained neighbors");
        };
        let start = extended_line_circle_endpoint(previous, frame, true)?;
        let end = extended_line_circle_endpoint(next, frame, false)?;
        carrier_runs.push((run, start, end));
    }
    ensure!(
        interface_edges.is_subset(treatment_edges),
        "partial circle carrier assigns an edge outside the treatment boundary"
    );
    let terminal_edges = treatment_edges
        .difference(&interface_edges)
        .copied()
        .collect::<BTreeSet<_>>();
    ensure!(
        terminal_edges.len() == 2
            && carrier_runs.len() == 2
            && carrier_runs
                .iter()
                .all(|(_, start, end)| distance(*start, *end) > GEOM_TOL_MM),
        "partial circle needs two terminal edges and two nondegenerate carrier arcs"
    );
    let (first_run, first_start, first_end) = &carrier_runs[0];
    let (_, second_start, second_end) = &carrier_runs[1];
    ensure!(
        (distance(*first_start, *second_start) <= GEOM_TOL_MM
            && distance(*first_end, *second_end) <= GEOM_TOL_MM)
            || (distance(*first_start, *second_end) <= GEOM_TOL_MM
                && distance(*first_end, *second_start) <= GEOM_TOL_MM),
        "partial circle carriers do not meet at the same inverse endpoints"
    );
    let (new_vertices, new_edges) =
        open_inverse_circle_arc(frame, *first_start, *first_end, first_run.sweep > 0.0)?;
    let (mut face_splices, edits) = partial_circle_carrier_splices(
        &carrier_runs,
        &new_vertices,
        new_edges.len(),
        data_section,
        provenance,
    )?;
    face_splices.extend(plan_terminal_splices(
        solid,
        evidence,
        data_section,
        provenance,
        &terminal_edges,
    )?);
    Ok(KernelTopologyPatch {
        remove_faces: vec![topology_face(
            provenance,
            data_section,
            evidence.treatment_face_id,
        )?],
        new_vertices,
        new_edges,
        rewrite_edges: normalize_retained_edge_rewrites(&edits, data_section, provenance)?,
        face_splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(
            provenance,
            data_section,
            treatment_edges.iter().copied(),
        )?,
    })
}

fn plan_circle_patch(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    center_mm: [f64; 3],
    normal: [f64; 3],
    radius_mm: f64,
) -> Result<KernelTopologyPatch> {
    ensure!(
        radius_mm.is_finite() && radius_mm > GEOM_TOL_MM,
        "inverse circle for treatment face #{} has invalid radius",
        evidence.treatment_face_id
    );
    let normal = normalize(normal, DIR_TOL).context("inverse circle normal is degenerate")?;
    let x_direction = deterministic_perpendicular(normal)
        .context("inverse circle cannot construct a deterministic radial axis")?;
    let treatment_edges = unique_face_edges(face(solid, evidence.treatment_face_id)?)?;
    let carrier_interfaces = evidence
        .carriers
        .iter()
        .flat_map(|carrier| carrier.shared_edge_ids.iter().copied())
        .collect::<BTreeSet<_>>();
    if carrier_interfaces.len() < treatment_edges.len() {
        return plan_partial_circle_patch(
            solid,
            evidence,
            data_section,
            provenance,
            CircleFrame::from_edge(&evidence.virtual_edge)?,
            &treatment_edges,
        );
    }
    let (face_splices, interface_edges) =
        plan_circle_carriers(solid, evidence, data_section, provenance, center_mm, normal)?;
    ensure!(
        interface_edges == treatment_edges,
        "circular treatment face #{} has terminal edges; partial-circle surgery is not normalized yet",
        evidence.treatment_face_id
    );
    let (new_vertices, new_edges) =
        circle_patch_geometry(center_mm, normal, x_direction, radius_mm);

    Ok(KernelTopologyPatch {
        remove_faces: vec![topology_face(
            provenance,
            data_section,
            evidence.treatment_face_id,
        )?],
        new_vertices,
        new_edges,
        rewrite_edges: Vec::new(),
        face_splices,
        face_boundary_replacements: Vec::new(),
        require_unused_edges: topology_edges(
            provenance,
            data_section,
            treatment_edges.iter().copied(),
        )?,
    })
}

fn plan_circle_carriers(
    solid: &ReducedSolid,
    evidence: &EdgeTreatmentEvidence,
    data_section: usize,
    provenance: &KernelTopologyProvenance,
    center_mm: [f64; 3],
    normal: [f64; 3],
) -> Result<(Vec<KernelFaceSplice>, BTreeSet<u64>)> {
    let mut interface_edges = BTreeSet::new();
    let mut face_splices = Vec::with_capacity(2);
    for carrier in &evidence.carriers {
        let carrier_face = face(solid, carrier.face_id)?;
        let shared = carrier
            .shared_edge_ids
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        ensure!(
            !shared.is_empty() && shared.len() == carrier.shared_edge_ids.len(),
            "circular carrier face #{} has invalid shared-edge set",
            carrier.face_id
        );
        for &edge_id in &shared {
            ensure!(
                interface_edges.insert(edge_id),
                "circular treatment face #{} assigns edge #{edge_id} to both carriers",
                evidence.treatment_face_id
            );
        }
        let loop_ = complete_carrier_loop(carrier_face, &shared)?;
        let forward = circular_loop_forward(loop_, center_mm, normal)?;
        face_splices.push(KernelFaceSplice {
            face: topology_face(provenance, data_section, carrier.face_id)?,
            remove_edges: topology_edges(provenance, data_section, shared.iter().copied())?,
            replacement: circle_replacement_uses(forward),
        });
    }
    Ok((face_splices, interface_edges))
}

fn complete_carrier_loop<'a>(
    face: &'a ReducedFace,
    shared: &BTreeSet<u64>,
) -> Result<&'a crate::brep::FaceLoop> {
    let matching = face
        .loops
        .iter()
        .filter(|loop_| {
            loop_.edges.len() == shared.len()
                && loop_
                    .edges
                    .iter()
                    .map(|edge| edge.edge_id)
                    .collect::<BTreeSet<_>>()
                    == *shared
        })
        .collect::<Vec<_>>();
    let [loop_] = matching.as_slice() else {
        bail!(
            "circular carrier face #{} does not expose the treatment as one complete loop",
            face.id
        );
    };
    Ok(loop_)
}

fn circle_replacement_uses(forward: bool) -> Vec<KernelPatchEdgeUse> {
    if forward {
        (0..4)
            .map(|edge| KernelPatchEdgeUse {
                edge: KernelPatchEdgeRef::New(edge),
                orientation: true,
            })
            .collect()
    } else {
        (0..4)
            .rev()
            .map(|edge| KernelPatchEdgeUse {
                edge: KernelPatchEdgeRef::New(edge),
                orientation: false,
            })
            .collect()
    }
}

fn circle_patch_geometry(
    center_mm: [f64; 3],
    normal: [f64; 3],
    x_direction: [f64; 3],
    radius_mm: f64,
) -> (Vec<[f64; 3]>, Vec<KernelNewEdge>) {
    let y_direction = cross(normal, x_direction);
    let quarter = std::f64::consts::FRAC_PI_2;
    let angles = [0.0, quarter, 2.0 * quarter, 3.0 * quarter];
    let vertices = angles
        .map(|angle| {
            let radial = [
                x_direction[0].mul_add(angle.cos(), y_direction[0] * angle.sin()),
                x_direction[1].mul_add(angle.cos(), y_direction[1] * angle.sin()),
                x_direction[2].mul_add(angle.cos(), y_direction[2] * angle.sin()),
            ];
            [
                radial[0].mul_add(radius_mm, center_mm[0]),
                radial[1].mul_add(radius_mm, center_mm[1]),
                radial[2].mul_add(radius_mm, center_mm[2]),
            ]
        })
        .to_vec();
    let ranges = [
        (0.0, quarter),
        (quarter, 2.0 * quarter),
        (2.0 * quarter, 3.0 * quarter),
        (3.0 * quarter, 4.0 * quarter),
    ];
    let edges = ranges
        .into_iter()
        .enumerate()
        .map(|(edge, (start_angle_rad, end_angle_rad))| KernelNewEdge {
            vertices: (
                KernelPatchVertexRef::New(edge),
                KernelPatchVertexRef::New((edge + 1) % 4),
            ),
            curve: KernelPatchCurve::CircleArc {
                center_mm,
                normal,
                x_direction,
                radius_mm,
                start_angle_rad,
                end_angle_rad,
            },
        })
        .collect();
    (vertices, edges)
}

fn circular_loop_forward(
    loop_: &crate::brep::FaceLoop,
    center_mm: [f64; 3],
    normal: [f64; 3],
) -> Result<bool> {
    for edge in &loop_.edges {
        let start = sub(edge.start_mm, center_mm);
        let end = sub(edge.end_mm, center_mm);
        let turn = dot(normal, cross(start, end));
        if turn.abs() > GEOM_TOL_MM * GEOM_TOL_MM {
            return Ok(turn > 0.0);
        }
        if let CurveSupport::Circle(circle) = edge.support {
            let alignment = dot(circle.normal, normal);
            if alignment.abs() > 1.0 - DIR_TOL {
                let parameter_forward = if edge.parameter_forward { 1.0 } else { -1.0 };
                return Ok(parameter_forward * alignment > 0.0);
            }
        }
    }
    bail!("circular carrier loop has no resolvable traversal direction")
}

fn deterministic_perpendicular(normal: [f64; 3]) -> Option<[f64; 3]> {
    let seed = if normal[0].abs() <= normal[1].abs() && normal[0].abs() <= normal[2].abs() {
        [1.0, 0.0, 0.0]
    } else if normal[1].abs() <= normal[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalize(cross(normal, seed), DIR_TOL)
}

fn face(solid: &ReducedSolid, face_id: u64) -> Result<&ReducedFace> {
    solid
        .faces()
        .find(|face| face.id == face_id)
        .with_context(|| format!("reduced solid #{} has no face #{face_id}", solid.solid_id()))
}

fn unique_face_edges(face: &ReducedFace) -> Result<BTreeSet<u64>> {
    let mut edges = BTreeSet::new();
    for edge in face.loop_edges() {
        ensure!(
            edges.insert(edge.edge_id),
            "face #{} uses source edge #{} more than once",
            face.id,
            edge.edge_id
        );
    }
    Ok(edges)
}

fn edge_face_users(solid: &ReducedSolid, edge_id: u64) -> BTreeSet<u64> {
    solid
        .faces()
        .filter(|face| face.loop_edges().any(|edge| edge.edge_id == edge_id))
        .map(|face| face.id)
        .collect()
}

fn edges_form_one_cyclic_run(face: &ReducedFace, edges: &BTreeSet<u64>) -> bool {
    let matches = face
        .loops
        .iter()
        .filter_map(|loop_| {
            let count = loop_
                .edges
                .iter()
                .filter(|edge| edges.contains(&edge.edge_id))
                .count();
            (count != 0).then_some((loop_, count))
        })
        .collect::<Vec<_>>();
    let [(loop_, count)] = matches.as_slice() else {
        return false;
    };
    if *count != edges.len() || edges.len() >= loop_.edges.len() {
        return false;
    }
    let edge_count = loop_.edges.len();
    let starts = (0..edge_count)
        .filter(|&index| {
            edges.contains(&loop_.edges[index].edge_id)
                && !edges.contains(&loop_.edges[(index + edge_count - 1) % edge_count].edge_id)
        })
        .count();
    starts == 1
}

fn topology_vertex(
    provenance: &KernelTopologyProvenance,
    data_section: usize,
    source_id: u64,
) -> Result<KernelTopologyId> {
    provenance
        .vertex_for_source_id(data_section, source_id)
        .with_context(|| format!("source vertex #{source_id} has no unique kernel StableId"))
}

fn topology_edge(
    provenance: &KernelTopologyProvenance,
    data_section: usize,
    source_id: u64,
) -> Result<KernelTopologyId> {
    provenance
        .edge_for_source_id(data_section, source_id)
        .with_context(|| format!("source edge #{source_id} has no unique kernel StableId"))
}

fn topology_edges(
    provenance: &KernelTopologyProvenance,
    data_section: usize,
    source_ids: impl IntoIterator<Item = u64>,
) -> Result<Vec<KernelTopologyId>> {
    source_ids
        .into_iter()
        .map(|source_id| topology_edge(provenance, data_section, source_id))
        .collect()
}

fn topology_face(
    provenance: &KernelTopologyProvenance,
    data_section: usize,
    source_id: u64,
) -> Result<KernelTopologyId> {
    provenance
        .face_for_source_id(data_section, source_id)
        .with_context(|| format!("source face #{source_id} has no unique kernel StableId"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SourceEntity;
    use crate::brep::{FaceLoop, LineSupport, PlaneSupport, SurfaceSupport};
    use crate::edge_treatments::{
        EdgeTreatmentCarrier, EdgeTreatmentEvidenceStrength, EdgeTreatmentKind,
    };
    use crate::kernel_source_load::{KernelFaceTopologyProvenance, KernelTopologyEntityProvenance};

    /// Native STEP face reversal must be applied to the *segment* carrier
    /// run as well as circular ones. This source has both forward and reversed
    /// carrier faces; every individual candidate and the simultaneous patch
    /// transaction must remain a one-shell, closed native solid.
    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_reversed_step_carriers_close_after_inverse_chamfer_replay() -> Result<()> {
        use crate::cad_kernel::monstertruck::MonstertruckKernel;
        use crate::kernel_source_load::SourceStepLoadInspector;
        use crate::kernel_topology::{SourceTopologyEditor, plan_source_edge_treatment_patches};

        let bytes =
            include_bytes!("../validation/fixtures/native_mixed_carrier_chamfer_06d018ae.step");
        let kernel = MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(bytes)?.exchange,
        )?;
        let inspected = kernel.inspect_source_document(&document)?;
        let loaded = inspected
            .solids
            .iter()
            .find(|solid| solid.source.id == 19)
            .context("20-face real source absent")?;
        let provenance = loaded
            .topology_provenance
            .as_ref()
            .context("source has no stable topology provenance")?;
        let section = &document.as_ref().data[loaded.source.data_section];
        let index = crate::step_graph::build_index(&section.entities);
        let evidence = crate::edge_treatments::detect_edge_treatment_evidence_for_solids_with_index(
            &section.entities,
            &index,
            &std::collections::HashSet::from([19]),
        );
        for face_ids in [
            vec![85],
            vec![172],
            vec![182],
            vec![382],
            vec![402],
            vec![452],
            vec![172, 182],
            vec![85, 172, 182, 382, 402, 452],
        ] {
            let selected = evidence
                .iter()
                .filter(|item| face_ids.contains(&item.treatment_face_id))
                .cloned()
                .collect::<Vec<_>>();
            ensure!(
                selected.len() == face_ids.len(),
                "lost real source treatments {face_ids:?}"
            );
            let patch = plan_source_edge_treatment_patches(
                document.as_ref(),
                loaded.source,
                &selected,
                provenance,
            )?;
            let edited = kernel
                .edit_source_topology(&document, loaded.source, &[patch])
                .with_context(|| format!("real source treatments {face_ids:?} are not closed"))?;
            let summary = kernel.inspect_edited_source(&edited)?;
            assert_eq!(summary.boundary_shells, 1, "{face_ids:?}");
            assert_eq!(summary.faces, 20 - face_ids.len(), "{face_ids:?}");
        }
        Ok(())
    }

    /// A real 81-face antenna body with three split-ring families, and a
    /// real 178-face perforated cylinder with a split fillet at each end.
    /// Source STEP representations are imported once; surgery operates
    /// entirely on Monstertruck's exact-trim native topology.
    #[cfg(feature = "cad-kernel-monstertruck")]
    fn verify_real_split_ring_native_transaction(
        input: &[u8],
        source_id: u64,
        treatment_ids: &[u64],
        original_faces: usize,
    ) -> Result<()> {
        use crate::cad_kernel::monstertruck::MonstertruckKernel;
        use crate::kernel_source_load::SourceStepLoadInspector;
        use crate::kernel_topology::{SourceTopologyEditor, plan_source_edge_treatment_patches};

        let kernel = MonstertruckKernel::project_standard();
        let document = kernel.resolve_source_exchange(
            monstertruck_io::step::load::ParsedStepSource::parse(input)?.exchange,
        )?;
        let source_load = kernel.inspect_source_document(&document)?;
        let loaded = source_load
            .solids
            .iter()
            .find(|item| item.source.id == source_id)
            .context("real split-ring fixture source is missing")?;
        let provenance = loaded
            .topology_provenance
            .as_ref()
            .context("real split-ring fixture has no stable topology provenance")?;
        let section = &document.as_ref().data[loaded.source.data_section];
        let index = crate::step_graph::build_index(&section.entities);
        let evidence =
            crate::edge_treatments::detect_edge_treatment_evidence_for_solids_with_index(
                &section.entities,
                &index,
                &std::collections::HashSet::from([source_id]),
            )
            .into_iter()
            .filter(|item| treatment_ids.contains(&item.treatment_face_id))
            .collect::<Vec<_>>();
        ensure!(
            evidence.len() == treatment_ids.len(),
            "fixture lost one or more split-ring treatments"
        );
        let patch = plan_source_edge_treatment_patches(
            document.as_ref(),
            loaded.source,
            &evidence,
            provenance,
        )?;
        assert_eq!(patch.remove_faces.len(), treatment_ids.len());
        let removed = patch.remove_faces.iter().copied().collect::<BTreeSet<_>>();
        let retained = provenance
            .faces
            .iter()
            .map(|face| face.topology_id)
            .filter(|id| !removed.contains(id))
            .collect::<BTreeSet<_>>();
        let edited = kernel.edit_source_topology(&document, loaded.source, &[patch])?;
        let actual = kernel.inspect_edited_source(&edited)?;
        assert_eq!(actual.boundary_shells, 1);
        assert_eq!(actual.faces, original_faces - treatment_ids.len());
        assert_eq!(
            actual
                .retained_source_topology
                .faces
                .iter()
                .map(|face| face.topology_id)
                .collect::<BTreeSet<_>>(),
            retained,
            "native surgery changed identities outside the selected treatment faces",
        );
        Ok(())
    }

    #[cfg(feature = "cad-kernel-monstertruck")]
    #[test]
    fn real_split_ring_families_apply_atomically_with_original_native_identity() -> Result<()> {
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_antenna.step"),
            4570,
            &[131, 1542, 463, 2621, 4057, 4159],
            81,
        )?;
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_antenna.step"),
            4570,
            &[2031, 3387],
            81,
        )?;
        // The same two open arcs must coexist atomically with the three
        // complete split-ring families. All neighboring source face IDs
        // outside these eight removed treatments remain unchanged.
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_antenna.step"),
            4570,
            &[131, 1542, 463, 2621, 2031, 3387, 4057, 4159],
            81,
        )?;
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_hollow_cylinder.step"),
            7849,
            &[4123],
            178,
        )?;
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_hollow_cylinder.step"),
            7849,
            &[7044],
            178,
        )?;
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_hollow_cylinder.step"),
            7849,
            &[3393, 3413, 7654, 7840],
            178,
        )?;
        verify_real_split_ring_native_transaction(
            include_bytes!("../validation/fixtures/native_split_ring_hollow_cylinder.step"),
            7849,
            &[3393, 3413, 4123, 7044, 7654, 7840],
            178,
        )
    }

    const DATA_SECTION: usize = 0;
    const SOLID_ID: u64 = 500;

    fn source(id: u64) -> SourceEntity {
        SourceEntity {
            data_section: DATA_SECTION,
            ordinal: usize::try_from(id).unwrap_or(0),
            id,
        }
    }

    fn kid(raw: u64) -> KernelTopologyId { KernelTopologyId::from_raw(raw) }

    fn edge(
        id: u64,
        vertices: (u64, u64),
        points: ([f64; 3], [f64; 3]),
        forward: bool,
    ) -> OrientedEdgeUse {
        let (start_vertex, end_vertex, start_mm, end_mm) = if forward {
            (vertices.0, vertices.1, points.0, points.1)
        } else {
            (vertices.1, vertices.0, points.1, points.0)
        };
        OrientedEdgeUse {
            oriented_edge_id: id + 10_000,
            edge_id: id,
            curve_id: id + 20_000,
            curve_same_sense: true,
            parameter_forward: forward,
            start_vertex,
            end_vertex,
            start_mm,
            end_mm,
            support: CurveSupport::Line(LineSupport {
                origin_mm: points.0,
                direction: sub(points.1, points.0),
            }),
        }
    }

    fn plane_face(id: u64, edges: Vec<OrientedEdgeUse>) -> ReducedFace {
        ReducedFace {
            id,
            surface: SurfaceSupport::Plane(PlaneSupport {
                origin_mm: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                max_residual_mm: 0.0,
            }),
            surface_forward: true,
            loops: vec![FaceLoop {
                bound_id: id + 30_000,
                loop_id: id + 40_000,
                outer: true,
                orientation: true,
                edges,
            }],
        }
    }

    const CHAMFER_POINTS: [[f64; 3]; 10] = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.5, 0.0],
        [1.5, 2.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 1.0],
        [2.0, 0.0, 1.0],
        [2.0, 1.5, 1.0],
        [1.5, 2.0, 1.0],
        [0.0, 2.0, 1.0],
    ];
    const CHAMFER_VERTEX_IDS: [u64; 10] = [101, 102, 103, 104, 105, 106, 107, 108, 109, 110];
    const CHAMFER_EDGE_IDS: [u64; 15] = [
        201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215,
    ];
    const CHAMFER_FACE_IDS: [u64; 7] = [301, 302, 303, 304, 305, 306, 307];
    const CHAMFER_EDGE_VERTICES: [(usize, usize); 15] = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 0),
        (5, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 5),
        (0, 5),
        (1, 6),
        (2, 7),
        (3, 8),
        (4, 9),
    ];

    fn chamfer_faces() -> Vec<ReducedFace> {
        let edge_use = |index: usize, forward: bool| {
            let (front, back) = CHAMFER_EDGE_VERTICES[index];
            edge(
                CHAMFER_EDGE_IDS[index],
                (CHAMFER_VERTEX_IDS[front], CHAMFER_VERTEX_IDS[back]),
                (CHAMFER_POINTS[front], CHAMFER_POINTS[back]),
                forward,
            )
        };
        let uses = |entries: &[(usize, bool)]| {
            entries
                .iter()
                .map(|&(index, forward)| edge_use(index, forward))
                .collect()
        };
        vec![
            plane_face(
                CHAMFER_FACE_IDS[0],
                uses(&[(0, true), (1, true), (2, true), (3, true), (4, true)]),
            ),
            plane_face(
                CHAMFER_FACE_IDS[1],
                uses(&[(9, false), (8, false), (7, false), (6, false), (5, false)]),
            ),
            plane_face(
                CHAMFER_FACE_IDS[2],
                uses(&[(0, false), (10, true), (5, true), (11, false)]),
            ),
            plane_face(
                CHAMFER_FACE_IDS[3],
                uses(&[(1, false), (11, true), (6, true), (12, false)]),
            ),
            plane_face(
                CHAMFER_FACE_IDS[4],
                uses(&[(2, false), (12, true), (7, true), (13, false)]),
            ),
            plane_face(
                CHAMFER_FACE_IDS[5],
                uses(&[(3, false), (13, true), (8, true), (14, false)]),
            ),
            plane_face(
                CHAMFER_FACE_IDS[6],
                uses(&[(4, false), (14, true), (9, true), (10, false)]),
            ),
        ]
    }

    fn chamfer_evidence() -> EdgeTreatmentEvidence {
        EdgeTreatmentEvidence {
            solid_id: SOLID_ID,
            treatment_face_id: CHAMFER_FACE_IDS[4],
            kind: EdgeTreatmentKind::PlanarChamfer,
            strength: EdgeTreatmentEvidenceStrength::AnalyticTangency,
            carrier_pair_unique: true,
            treatment_support: "plane".to_owned(),
            carriers: [
                EdgeTreatmentCarrier {
                    face_id: CHAMFER_FACE_IDS[3],
                    support: "plane".to_owned(),
                    shared_edge_ids: vec![CHAMFER_EDGE_IDS[12]],
                },
                EdgeTreatmentCarrier {
                    face_id: CHAMFER_FACE_IDS[5],
                    support: "plane".to_owned(),
                    shared_edge_ids: vec![CHAMFER_EDGE_IDS[13]],
                },
            ],
            inferred_radius_mm: None,
            inferred_setbacks_mm: Some([0.5, 0.5]),
            inferred_angle_rad: Some(std::f64::consts::FRAC_PI_4),
            carrier_angle_rad: std::f64::consts::FRAC_PI_2,
            virtual_edge: EdgeTreatmentVirtualEdge::Segment {
                start_mm: [2.0, 2.0, 0.0],
                end_mm: [2.0, 2.0, 1.0],
            },
            max_boundary_residual_mm: 0.0,
        }
    }

    fn topology_provenance(
        vertex_ids: &[u64],
        edge_ids: &[u64],
        face_ids: &[u64],
    ) -> KernelTopologyProvenance {
        KernelTopologyProvenance {
            vertices: vertex_ids
                .iter()
                .enumerate()
                .map(|(index, &id)| KernelTopologyEntityProvenance {
                    source: source(id),
                    topology_id: kid(1_000 + u64::try_from(index).unwrap_or(0)),
                })
                .collect(),
            edges: edge_ids
                .iter()
                .enumerate()
                .map(|(index, &id)| KernelTopologyEntityProvenance {
                    source: source(id),
                    topology_id: kid(2_000 + u64::try_from(index).unwrap_or(0)),
                })
                .collect(),
            faces: face_ids
                .iter()
                .enumerate()
                .map(|(index, &id)| KernelFaceTopologyProvenance {
                    source: None,
                    shell_member: Some(source(id)),
                    topology_id: kid(3_000 + u64::try_from(index).unwrap_or(0)),
                })
                .collect(),
        }
    }

    fn chamfer_fixture() -> (
        ReducedSolid,
        EdgeTreatmentEvidence,
        KernelTopologyProvenance,
        [u64; 10],
        [u64; 15],
        [u64; 7],
    ) {
        (
            ReducedSolid::from_test_faces(SOLID_ID, chamfer_faces()),
            chamfer_evidence(),
            topology_provenance(&CHAMFER_VERTEX_IDS, &CHAMFER_EDGE_IDS, &CHAMFER_FACE_IDS),
            CHAMFER_VERTEX_IDS,
            CHAMFER_EDGE_IDS,
            CHAMFER_FACE_IDS,
        )
    }

    #[test]
    fn chamfer_planner_emits_complete_identity_preserving_patch() -> Result<()> {
        let (solid, evidence, provenance, vertex_ids, edge_ids, face_ids) = chamfer_fixture();
        let patch = plan_edge_treatment_patch(&solid, &evidence, DATA_SECTION, &provenance)?;

        assert_eq!(patch.remove_faces, vec![kid(3_004)]);
        assert_eq!(patch.new_vertices, vec![[2.0, 2.0, 0.0], [2.0, 2.0, 1.0]]);
        assert_eq!(patch.new_edges.len(), 1);
        assert_eq!(
            patch.new_edges[0].vertices,
            (KernelPatchVertexRef::New(0), KernelPatchVertexRef::New(1))
        );
        assert!(matches!(
            patch.new_edges[0].curve,
            KernelPatchCurve::Line {
                start_mm: [2.0, 2.0, 0.0],
                end_mm: [2.0, 2.0, 1.0]
            }
        ));

        let rewritten = patch
            .rewrite_edges
            .iter()
            .map(|edge| edge.edge)
            .collect::<Vec<_>>();
        assert_eq!(
            rewritten,
            vec![kid(2_001), kid(2_003), kid(2_006), kid(2_008)]
        );
        let expected_unused = [edge_ids[2], edge_ids[7], edge_ids[12], edge_ids[13]]
            .into_iter()
            .map(|id| {
                provenance
                    .edge_for_source_id(DATA_SECTION, id)
                    .context("fixture edge provenance")
            })
            .collect::<Result<Vec<_>>>()?;
        assert_eq!(patch.require_unused_edges, expected_unused);

        let by_face = patch
            .face_splices
            .iter()
            .map(|splice| (splice.face, splice))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(by_face.len(), 4);
        for terminal in [face_ids[0], face_ids[1]] {
            let face_id = provenance
                .face_for_source_id(DATA_SECTION, terminal)
                .context("fixture terminal face provenance")?;
            assert_eq!(by_face[&face_id].replacement, Vec::new());
        }
        for carrier in [face_ids[3], face_ids[5]] {
            let face_id = provenance
                .face_for_source_id(DATA_SECTION, carrier)
                .context("fixture carrier face provenance")?;
            assert_eq!(by_face[&face_id].replacement.len(), 1);
            assert_eq!(
                by_face[&face_id].replacement[0].edge,
                KernelPatchEdgeRef::New(0)
            );
        }

        // The untouched endpoint of every rewritten source edge still refers to its original
        // imported vertex identity; only the tangency endpoint becomes a new derived vertex.
        let expected_untouched = [vertex_ids[1], vertex_ids[4], vertex_ids[6], vertex_ids[9]]
            .into_iter()
            .map(|id| {
                provenance
                    .vertex_for_source_id(DATA_SECTION, id)
                    .context("fixture vertex provenance")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        let actual_untouched = patch
            .rewrite_edges
            .iter()
            .flat_map(|edge| <[KernelPatchVertexRef; 2]>::from(edge.vertices))
            .filter_map(|reference| match reference {
                KernelPatchVertexRef::Existing(id) => Some(id),
                KernelPatchVertexRef::New(_) => None,
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(actual_untouched, expected_untouched);
        Ok(())
    }

    fn circular_edge(id: u64, vertex_id: u64, forward: bool) -> OrientedEdgeUse {
        let seam = [1.0, 0.0, 0.0];
        OrientedEdgeUse {
            oriented_edge_id: id + 10_000,
            edge_id: id,
            curve_id: id + 20_000,
            curve_same_sense: true,
            parameter_forward: forward,
            start_vertex: vertex_id,
            end_vertex: vertex_id,
            start_mm: seam,
            end_mm: seam,
            support: CurveSupport::Circle(crate::brep::CircleSupport {
                center_mm: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                x_direction: [1.0, 0.0, 0.0],
                radius_mm: 1.0,
            }),
        }
    }

    fn circular_face(id: u64, edge: OrientedEdgeUse) -> ReducedFace { plane_face(id, vec![edge]) }

    fn circular_fixture() -> (
        ReducedSolid,
        EdgeTreatmentEvidence,
        KernelTopologyProvenance,
    ) {
        let first_face = 401;
        let second_face = 402;
        let treatment_face = 403;
        let first_edge = 501;
        let second_edge = 502;
        let treatment = ReducedFace {
            id: treatment_face,
            surface: SurfaceSupport::Other { entity_id: 9999 },
            surface_forward: true,
            loops: vec![
                FaceLoop {
                    bound_id: 701,
                    loop_id: 702,
                    outer: true,
                    orientation: true,
                    edges: vec![circular_edge(first_edge, 601, false)],
                },
                FaceLoop {
                    bound_id: 703,
                    loop_id: 704,
                    outer: false,
                    orientation: true,
                    edges: vec![circular_edge(second_edge, 602, true)],
                },
            ],
        };
        let solid = ReducedSolid::from_test_faces(
            SOLID_ID,
            vec![
                circular_face(first_face, circular_edge(first_edge, 601, true)),
                circular_face(second_face, circular_edge(second_edge, 602, false)),
                treatment,
            ],
        );
        let evidence = EdgeTreatmentEvidence {
            solid_id: SOLID_ID,
            treatment_face_id: treatment_face,
            kind: EdgeTreatmentKind::ToroidalFillet,
            strength: EdgeTreatmentEvidenceStrength::AnalyticTangency,
            carrier_pair_unique: true,
            treatment_support: "torus".to_owned(),
            carriers: [
                EdgeTreatmentCarrier {
                    face_id: first_face,
                    support: "plane".to_owned(),
                    shared_edge_ids: vec![first_edge],
                },
                EdgeTreatmentCarrier {
                    face_id: second_face,
                    support: "cylinder".to_owned(),
                    shared_edge_ids: vec![second_edge],
                },
            ],
            inferred_radius_mm: Some(0.1),
            inferred_setbacks_mm: None,
            inferred_angle_rad: None,
            carrier_angle_rad: std::f64::consts::FRAC_PI_2,
            virtual_edge: EdgeTreatmentVirtualEdge::Circle {
                center_mm: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                radius_mm: 1.0,
            },
            max_boundary_residual_mm: 0.0,
        };
        let provenance = KernelTopologyProvenance {
            vertices: Vec::new(),
            edges: [first_edge, second_edge]
                .into_iter()
                .enumerate()
                .map(|(index, id)| KernelTopologyEntityProvenance {
                    source: source(id),
                    topology_id: kid(5_000 + u64::try_from(index).unwrap_or(0)),
                })
                .collect(),
            faces: [first_face, second_face, treatment_face]
                .into_iter()
                .enumerate()
                .map(|(index, id)| KernelFaceTopologyProvenance {
                    source: None,
                    shell_member: Some(source(id)),
                    topology_id: kid(6_000 + u64::try_from(index).unwrap_or(0)),
                })
                .collect(),
        };
        (solid, evidence, provenance)
    }

    #[test]
    fn circular_planner_replaces_illegal_single_edge_loops_with_shared_four_arc_cycle() -> Result<()>
    {
        let (solid, evidence, provenance) = circular_fixture();
        let patch = plan_edge_treatment_patch(&solid, &evidence, DATA_SECTION, &provenance)?;
        assert_eq!(patch.new_vertices.len(), 4);
        assert_eq!(patch.new_edges.len(), 4);
        assert_eq!(patch.rewrite_edges, Vec::new());
        assert_eq!(patch.face_splices.len(), 2);
        assert_eq!(
            patch.face_splices[0]
                .replacement
                .iter()
                .map(|edge| (edge.edge, edge.orientation))
                .collect::<Vec<_>>(),
            vec![
                (KernelPatchEdgeRef::New(0), true),
                (KernelPatchEdgeRef::New(1), true),
                (KernelPatchEdgeRef::New(2), true),
                (KernelPatchEdgeRef::New(3), true),
            ]
        );
        assert_eq!(
            patch.face_splices[1]
                .replacement
                .iter()
                .map(|edge| (edge.edge, edge.orientation))
                .collect::<Vec<_>>(),
            vec![
                (KernelPatchEdgeRef::New(3), false),
                (KernelPatchEdgeRef::New(2), false),
                (KernelPatchEdgeRef::New(1), false),
                (KernelPatchEdgeRef::New(0), false),
            ]
        );
        assert!(
            patch
                .new_edges
                .iter()
                .all(|edge| edge.vertices.0 != edge.vertices.1)
        );
        Ok(())
    }

    #[test]
    fn planner_rejects_incomplete_terminal_topology() -> Result<()> {
        let (mut solid, evidence, provenance, _, edge_ids, face_ids) = chamfer_fixture();
        let faces = solid
            .faces()
            .filter(|face| face.id != face_ids[0])
            .cloned()
            .collect();
        solid = ReducedSolid::from_test_faces(SOLID_ID, faces);
        let error = plan_edge_treatment_patch(&solid, &evidence, DATA_SECTION, &provenance)
            .err()
            .context("missing terminal face unexpectedly planned successfully")?;
        assert!(error.to_string().contains(&format!("#{}", edge_ids[2])));
        Ok(())
    }
}
