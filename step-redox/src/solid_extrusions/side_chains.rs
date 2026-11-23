//! Complete proof of a ruled side whose authored cap or generator curves may
//! be subdivided differently at opposite ends. Match *geometric chains*, not
//! the number of STEP edges representing them.

use super::{
    ANGLE_TOL_RAD, DIR_TOL, FaceInfo, GEOM_TOL_MM, OrientedEdgeUse, TAU, add, circle_edge_sweep,
    cross, distance, dot, norm, normalize, raw_loop_is_continuous, side_support_residual, sub,
    translated_profile_edge_residual,
};
use crate::brep::{CircleSupport, CurveSupport};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Region {
    Bottom,
    Top,
    Generator,
}

struct SideRuns<'a> {
    bottom: Vec<&'a OrientedEdgeUse>,
    top: Vec<&'a OrientedEdgeUse>,
    generators: [Vec<&'a OrientedEdgeUse>; 2],
}

pub(super) struct SideProofContext<'a> {
    bottom_edges: &'a HashSet<u64>,
    top_edges: &'a HashSet<u64>,
    bottom_loops: HashMap<u64, usize>,
    top_loops: HashMap<u64, usize>,
    bottom_roles: &'a [bool],
    top_roles: &'a [bool],
    loop_match: HashMap<usize, usize>,
    reverse_match: HashMap<usize, usize>,
    extrusion: [f64; 3],
}

impl<'a> SideProofContext<'a> {
    pub(super) fn new(
        bottom: &FaceInfo,
        top: &FaceInfo,
        bottom_edges: &'a HashSet<u64>,
        top_edges: &'a HashSet<u64>,
        extrusion: [f64; 3],
        bottom_roles: &'a [bool],
        top_roles: &'a [bool],
    ) -> Option<Self> {
        fn loop_map(face: &FaceInfo) -> Option<HashMap<u64, usize>> {
            let mut result = HashMap::new();
            for (index, wire) in face.loops.iter().enumerate() {
                for edge in &wire.edges {
                    if result.insert(edge.edge_id, index).is_some() {
                        return None;
                    }
                }
            }
            Some(result)
        }
        let bottom_loops = loop_map(bottom)?;
        let top_loops = loop_map(top)?;
        (bottom_loops.len() == bottom_edges.len()
            && top_loops.len() == top_edges.len()
            && bottom_roles.len() == bottom.loops.len()
            && top_roles.len() == top.loops.len()
            && !bottom_edges.is_empty()
            && !top_edges.is_empty())
        .then_some(Self {
            bottom_edges,
            top_edges,
            bottom_loops,
            top_loops,
            bottom_roles,
            top_roles,
            loop_match: HashMap::new(),
            reverse_match: HashMap::new(),
            extrusion,
        })
    }

    pub(super) fn complete(&self) -> bool {
        self.loop_match.len() == self.bottom_roles.len()
            && self.reverse_match.len() == self.top_roles.len()
            && self.loop_match.iter().all(|(&bottom, &top)| {
                self.bottom_roles[bottom] == self.top_roles[top]
                    && self.reverse_match.get(&top) == Some(&bottom)
            })
    }

    fn match_loops(&mut self, bottom: usize, top: usize) -> Option<()> {
        if self.bottom_roles.get(bottom)? != self.top_roles.get(top)?
            || self
                .loop_match
                .get(&bottom)
                .is_some_and(|&prior| prior != top)
            || self
                .reverse_match
                .get(&top)
                .is_some_and(|&prior| prior != bottom)
        {
            return None;
        }
        self.loop_match.insert(bottom, top);
        self.reverse_match.insert(top, bottom);
        Some(())
    }
}

fn side_cap_runs<'a>(side: &'a FaceInfo, context: &SideProofContext<'_>) -> Option<SideRuns<'a>> {
    if side.loops.len() != 1 || !raw_loop_is_continuous(&side.loop_edges) {
        return None;
    }
    let edges = &side.loop_edges;
    let regions = edges
        .iter()
        .map(|edge| {
            if context.bottom_edges.contains(&edge.edge_id) {
                Region::Bottom
            } else if context.top_edges.contains(&edge.edge_id) {
                Region::Top
            } else {
                Region::Generator
            }
        })
        .collect::<Vec<_>>();
    // Start at a run boundary so that the cyclic partition never splits the
    // last generator into the first and last vector entries.
    let start = (0..regions.len())
        .find(|&i| regions[i] != regions[(i + regions.len() - 1) % regions.len()])?;
    let mut groups = Vec::<(Region, Vec<&OrientedEdgeUse>)>::new();
    for step in 0..regions.len() {
        let index = (start + step) % regions.len();
        let kind = regions[index];
        if groups.last().is_none_or(|(last, _)| *last != kind) {
            groups.push((kind, Vec::new()));
        }
        groups.last_mut()?.1.push(&edges[index]);
    }
    if groups.len() != 4 {
        return None;
    }
    let mut bottom = None;
    let mut top = None;
    let mut generators = Vec::new();
    for (kind, group) in groups {
        match kind {
            Region::Bottom => {
                if bottom.replace(group).is_some() {
                    return None;
                }
            }
            Region::Top => {
                if top.replace(group).is_some() {
                    return None;
                }
            }
            Region::Generator => generators.push(group),
        }
    }
    let [first, second] = generators.try_into().ok()?;
    Some(SideRuns {
        bottom: bottom?,
        top: top?,
        generators: [first, second],
    })
}

fn run_ends(edges: &[&OrientedEdgeUse]) -> Option<[[f64; 3]; 2]> {
    Some([edges.first()?.start_mm, edges.last()?.end_mm])
}

/// Each straight piece advances along one geometric segment. This is enough
/// to establish complete coverage of that segment, irrespective of splits.
fn line_chain(edges: &[&OrientedEdgeUse]) -> Option<f64> {
    let [start, end] = run_ends(edges)?;
    let direction = normalize(sub(end, start))?;
    let mut max_residual = 0.0_f64;
    for edge in edges {
        let CurveSupport::Line(line) = &edge.support else {
            return None;
        };
        let support_axis = normalize(line.direction)?;
        let move_axis = sub(edge.end_mm, edge.start_mm);
        if dot(move_axis, direction) <= GEOM_TOL_MM
            || dot(support_axis, direction).abs() < 1.0 - DIR_TOL
        {
            return None;
        }
        for point in [edge.start_mm, edge.end_mm] {
            let deviation = norm(cross(sub(point, start), direction))
                .max(norm(cross(sub(point, line.origin_mm), support_axis)));
            if deviation > GEOM_TOL_MM {
                return None;
            }
            max_residual = max_residual.max(deviation);
        }
    }
    Some(max_residual)
}

struct CircleChain {
    support: CircleSupport,
    sweep: f64,
    max_residual: f64,
}

/// The signed sum of consecutive supported arcs is the exact complete arc.
/// Forbid reversals and repeat traversal; merely matching endpoints would
/// confuse the complementary semicircle with the intended trim.
fn circle_chain(edges: &[&OrientedEdgeUse]) -> Option<CircleChain> {
    let CurveSupport::Circle(reference) = &edges.first()?.support else {
        return None;
    };
    let mut sweep = 0.0_f64;
    let mut sign = None;
    let mut max_residual = 0.0_f64;
    for edge in edges {
        let CurveSupport::Circle(circle) = &edge.support else {
            return None;
        };
        let axis_cos = dot(reference.normal, circle.normal);
        let residual = distance(reference.center_mm, circle.center_mm)
            .max((reference.radius_mm - circle.radius_mm).abs())
            .max(norm(cross(reference.normal, circle.normal)) * reference.radius_mm);
        if residual > GEOM_TOL_MM || axis_cos.abs() < 1.0 - DIR_TOL {
            return None;
        }
        let angle = circle_edge_sweep(edge, *circle)? * axis_cos.signum();
        if angle.abs() <= ANGLE_TOL_RAD
            || sign.is_some_and(|positive| positive != angle.is_sign_positive())
        {
            return None;
        }
        sign = Some(angle.is_sign_positive());
        sweep += angle;
        if sweep.abs() > TAU + ANGLE_TOL_RAD {
            return None;
        }
        max_residual = max_residual.max(residual);
    }
    (sweep.abs() < TAU - ANGLE_TOL_RAD).then_some(CircleChain {
        support: *reference,
        sweep,
        max_residual,
    })
}

/// Support-level correspondence. Prefer individually certified arcs when both
/// authors used identical segmentation; otherwise prove the *whole* line or
/// circle chain, so different segment counts do not invalidate the same part.
fn translated_chains(
    bottom: &[&OrientedEdgeUse],
    top: &[&OrientedEdgeUse],
    extrusion: [f64; 3],
) -> Option<f64> {
    let [b0, b1] = run_ends(bottom)?;
    let [t0, t1] = run_ends(top)?;
    let endpoint_residual = distance(add(b0, extrusion), t1).max(distance(add(b1, extrusion), t0));
    if endpoint_residual > GEOM_TOL_MM {
        return None;
    }
    if bottom.len() == top.len() {
        let pairwise = bottom
            .iter()
            .zip(top.iter().rev())
            .try_fold(0.0_f64, |max, (a, b)| {
                let ends = distance(add(a.start_mm, extrusion), b.end_mm)
                    .max(distance(add(a.end_mm, extrusion), b.start_mm));
                let support = translated_profile_edge_residual(a, b, extrusion)?;
                (ends.max(support) <= GEOM_TOL_MM).then_some(max.max(ends.max(support)))
            });
        if let Some(residual) = pairwise {
            return Some(residual.max(endpoint_residual));
        }
    }
    if bottom
        .iter()
        .all(|e| matches!(e.support, CurveSupport::Line(_)))
        && top
            .iter()
            .all(|e| matches!(e.support, CurveSupport::Line(_)))
    {
        let b = line_chain(bottom)?;
        let t = line_chain(top)?;
        return Some(endpoint_residual.max(b).max(t));
    }
    let b = circle_chain(bottom)?;
    let t = circle_chain(top)?;
    let residual = distance(add(b.support.center_mm, extrusion), t.support.center_mm)
        .max((b.support.radius_mm - t.support.radius_mm).abs())
        .max(norm(cross(b.support.normal, t.support.normal)) * b.support.radius_mm);
    if residual > GEOM_TOL_MM || dot(b.support.normal, t.support.normal).abs() < 1.0 - DIR_TOL {
        return None;
    }
    let top_sign = dot(b.support.normal, t.support.normal).signum();
    let angle_delta = top_sign.mul_add(t.sweep, b.sweep).abs();
    let angular_residual = angle_delta * b.support.radius_mm;
    (angular_residual <= GEOM_TOL_MM && angle_delta <= ANGLE_TOL_RAD).then_some(
        endpoint_residual
            .max(residual)
            .max(angular_residual)
            .max(b.max_residual)
            .max(t.max_residual),
    )
}

/// A generator can itself have multiple collinear pieces. Their ordered union
/// must traverse the complete bottom-to-top axial interval exactly once.
fn generator_chain(edges: &[&OrientedEdgeUse], extrusion: [f64; 3]) -> Option<f64> {
    let [start, end] = run_ends(edges)?;
    let forward = distance(sub(end, start), extrusion);
    let reverse = distance(sub(start, end), extrusion);
    let endpoint_residual = forward.min(reverse);
    if endpoint_residual > GEOM_TOL_MM {
        return None;
    }
    Some(endpoint_residual.max(line_chain(edges)?))
}

/// Every cap edge is consumed exactly once; every side creates one consistent
/// bottom/top loop correspondence, independently of cap-edge subdivision.
pub(super) fn prove_side_face(
    side: &FaceInfo,
    context: &mut SideProofContext<'_>,
    used_bottom_edges: &mut HashSet<u64>,
    used_top_edges: &mut HashSet<u64>,
) -> Option<f64> {
    let runs = side_cap_runs(side, context)?;
    let bottom_loop = *context.bottom_loops.get(&runs.bottom.first()?.edge_id)?;
    let top_loop = *context.top_loops.get(&runs.top.first()?.edge_id)?;
    if !runs
        .bottom
        .iter()
        .all(|e| context.bottom_loops.get(&e.edge_id) == Some(&bottom_loop))
        || !runs
            .top
            .iter()
            .all(|e| context.top_loops.get(&e.edge_id) == Some(&top_loop))
    {
        return None;
    }
    let mut residual = translated_chains(&runs.bottom, &runs.top, context.extrusion)?;
    for edge in &runs.bottom {
        let on_carrier = side_support_residual(side, edge, context.extrusion)?;
        if on_carrier > GEOM_TOL_MM {
            return None;
        }
        residual = residual.max(on_carrier);
    }
    for generator in &runs.generators {
        residual = residual.max(generator_chain(generator, context.extrusion)?);
    }
    if runs
        .bottom
        .iter()
        .any(|e| used_bottom_edges.contains(&e.edge_id))
        || runs.top.iter().any(|e| used_top_edges.contains(&e.edge_id))
    {
        return None;
    }
    context.match_loops(bottom_loop, top_loop)?;
    used_bottom_edges.extend(runs.bottom.iter().map(|e| e.edge_id));
    used_top_edges.extend(runs.top.iter().map(|e| e.edge_id));
    Some(residual)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{CylinderSupport, FaceLoop, LineSupport, PlaneSupport, SurfaceSupport};

    fn point_id(p: [f64; 3]) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hash = std::collections::hash_map::DefaultHasher::new();
        p.map(f64::to_bits).hash(&mut hash);
        hash.finish()
    }

    fn line(id: u64, start: [f64; 3], end: [f64; 3]) -> OrientedEdgeUse {
        let Some(direction) = normalize(sub(end, start)) else {
            panic!("test segment must have positive length");
        };
        OrientedEdgeUse {
            oriented_edge_id: id + 1000,
            edge_id: id,
            curve_id: id + 2000,
            curve_same_sense: true,
            parameter_forward: true,
            start_vertex: point_id(start),
            end_vertex: point_id(end),
            start_mm: start,
            end_mm: end,
            support: CurveSupport::Line(LineSupport {
                origin_mm: start,
                direction,
            }),
        }
    }

    fn circle(id: u64, start: [f64; 3], end: [f64; 3], forward: bool) -> OrientedEdgeUse {
        let center_mm = [0., 0., start[2]];
        OrientedEdgeUse {
            oriented_edge_id: id + 1000,
            edge_id: id,
            curve_id: id + 2000,
            curve_same_sense: true,
            parameter_forward: forward,
            start_vertex: point_id(start),
            end_vertex: point_id(end),
            start_mm: start,
            end_mm: end,
            support: CurveSupport::Circle(CircleSupport {
                center_mm,
                normal: [0., 0., 1.],
                x_direction: [1., 0., 0.],
                radius_mm: 1.,
            }),
        }
    }

    fn face(id: u64, surface: SurfaceSupport, edges: Vec<OrientedEdgeUse>) -> FaceInfo {
        FaceInfo {
            id,
            surface,
            surface_forward: true,
            loops: vec![FaceLoop {
                bound_id: id + 100,
                loop_id: id + 200,
                outer: true,
                orientation: true,
                edges: edges.clone(),
            }],
            loop_edges: edges,
        }
    }

    fn plane(z: f64) -> SurfaceSupport {
        SurfaceSupport::Plane(PlaneSupport {
            origin_mm: [0., 0., z],
            normal: [0., 0., 1.],
            max_residual_mm: 0.,
        })
    }

    #[test]
    fn independently_segmented_caps_and_generators_are_proven_without_edge_pairs() {
        let bottom = face(
            10,
            plane(0.),
            vec![
                line(1, [0., 0., 0.], [0.5, 0., 0.]),
                line(2, [0.5, 0., 0.], [1., 0., 0.]),
                line(11, [1., 0., 0.], [1., 1., 0.]),
                line(12, [1., 1., 0.], [0., 1., 0.]),
                line(13, [0., 1., 0.], [0., 0., 0.]),
            ],
        );
        let top = face(
            20,
            plane(1.),
            vec![
                line(5, [0., 0., 1.], [1., 0., 1.]),
                line(14, [1., 0., 1.], [1., 1., 1.]),
                line(15, [1., 1., 1.], [0., 1., 1.]),
                line(16, [0., 1., 1.], [0., 0., 1.]),
            ],
        );
        let wall = vec![
            bottom.loop_edges[0].clone(),
            bottom.loop_edges[1].clone(),
            line(7, [1., 0., 0.], [1., 0., 0.4]),
            line(8, [1., 0., 0.4], [1., 0., 1.]),
            line(5, [1., 0., 1.], [0., 0., 1.]),
            line(9, [0., 0., 1.], [0., 0., 0.6]),
            line(10, [0., 0., 0.6], [0., 0., 0.]),
        ];
        let mut side = face(
            30,
            SurfaceSupport::Plane(PlaneSupport {
                origin_mm: [0., 0., 0.],
                normal: [0., 1., 0.],
                max_residual_mm: 0.,
            }),
            wall,
        );
        let bottom_ids = bottom.loop_edges.iter().map(|e| e.edge_id).collect();
        let top_ids = top.loop_edges.iter().map(|e| e.edge_id).collect();
        for _ in 0..side.loop_edges.len() {
            let mut context = SideProofContext::new(
                &bottom,
                &top,
                &bottom_ids,
                &top_ids,
                [0., 0., 1.],
                &[true],
                &[true],
            )
            .unwrap_or_else(|| panic!("cap mapping"));
            let mut used_bottom = HashSet::new();
            let mut used_top = HashSet::new();
            assert!(
                prove_side_face(&side, &mut context, &mut used_bottom, &mut used_top)
                    .is_some_and(|residual| residual <= GEOM_TOL_MM)
            );
            assert_eq!(used_bottom, HashSet::from([1, 2]));
            assert_eq!(used_top, HashSet::from([5]));
            assert!(context.complete());
            side.loop_edges.rotate_left(1);
        }
        side.loop_edges[0].end_mm[1] += 0.0001;
        let mut context = SideProofContext::new(
            &bottom,
            &top,
            &bottom_ids,
            &top_ids,
            [0., 0., 1.],
            &[true],
            &[true],
        )
        .unwrap_or_else(|| panic!("cap mapping"));
        assert!(
            prove_side_face(
                &side,
                &mut context,
                &mut HashSet::new(),
                &mut HashSet::new()
            )
            .is_none()
        );
    }

    #[test]
    fn differently_split_circular_caps_prove_only_the_same_signed_arc() {
        let a = [1., 0., 0.];
        let b = [0., 1., 0.];
        let c = [-1., 0., 0.];
        let top_of = |p: [f64; 3]| [p[0], p[1], 1.];
        let bottom = face(
            10,
            plane(0.),
            vec![
                circle(1, a, b, true),
                circle(2, b, c, true),
                circle(3, c, a, true),
            ],
        );
        let top = face(
            20,
            plane(1.),
            vec![
                circle(5, top_of(a), top_of(c), true),
                circle(6, top_of(c), top_of(a), true),
            ],
        );
        let mut side = face(
            30,
            SurfaceSupport::Cylinder(CylinderSupport {
                axis_origin_mm: [0., 0., 0.],
                axis: [0., 0., 1.],
                x_direction: [1., 0., 0.],
                radius_mm: 1.,
            }),
            vec![
                bottom.loop_edges[0].clone(),
                bottom.loop_edges[1].clone(),
                line(7, c, top_of(c)),
                circle(5, top_of(c), top_of(a), false),
                line(8, top_of(a), a),
            ],
        );
        let bottom_ids = bottom.loop_edges.iter().map(|e| e.edge_id).collect();
        let top_ids = top.loop_edges.iter().map(|e| e.edge_id).collect();
        let validate = |side: &FaceInfo| {
            let mut context = SideProofContext::new(
                &bottom,
                &top,
                &bottom_ids,
                &top_ids,
                [0., 0., 1.],
                &[true],
                &[true],
            )
            .unwrap_or_else(|| panic!("cap mapping"));
            prove_side_face(side, &mut context, &mut HashSet::new(), &mut HashSet::new())
        };
        assert!(validate(&side).is_some_and(|residual| residual <= GEOM_TOL_MM));
        side.loop_edges[3].parameter_forward = true;
        assert!(
            validate(&side).is_none(),
            "complementary half-circle must be rejected"
        );
    }
}
