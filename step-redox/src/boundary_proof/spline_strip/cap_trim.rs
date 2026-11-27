//! Independent full-face trim witness for source spline-cap precision audits.
//!
//! Proximity of a source curve to its surface row alone cannot certify a STEP
//! face: `EDGE_CURVE` may trim a longer basis, or the face might contain holes,
//! dangling links or off-surface generators. Admit two whole cap edge uses
//! and two completely covered straight generator chains, using the original
//! vertex identities and *every* edge of the authored face. This is a
//! geometric trim witness, not the strict whole-solid boundary equivalence.

use super::BSplineSupport;
use crate::{
    brep::{self, CurveSupport, OrientedEdgeUse},
    math3::{cross, distance, dot, norm, normalize, sub},
};

fn endpoints_match(edge: &OrientedEdgeUse, profile: &BSplineSupport, tolerance: f64) -> bool {
    let Some((&start, &end)) = profile
        .control_points_mm
        .first()
        .zip(profile.control_points_mm.last())
    else {
        return false;
    };
    distance(edge.start_mm, start)
        .max(distance(edge.end_mm, end))
        .min(distance(edge.start_mm, end).max(distance(edge.end_mm, start)))
        <= tolerance
}

fn run_between(
    edges: &[OrientedEdgeUse],
    first: usize,
    second: usize,
) -> Option<Vec<&OrientedEdgeUse>> {
    let mut result = Vec::new();
    let mut index = (first + 1) % edges.len();
    while index != second {
        result.push(&edges[index]);
        index = (index + 1) % edges.len();
        if result.len() >= edges.len() {
            return None;
        }
    }
    (!result.is_empty()).then_some(result)
}

fn verify_generator_run(
    edges: &[&OrientedEdgeUse],
    reference: [[f64; 3]; 2],
    tolerance: f64,
) -> bool {
    let Some(first) = edges.first() else {
        return false;
    };
    let Some(last) = edges.last() else {
        return false;
    };
    let (start, end) = if distance(first.start_mm, reference[0])
        .max(distance(last.end_mm, reference[1]))
        <= tolerance
    {
        reference.into()
    } else if distance(first.start_mm, reference[1]).max(distance(last.end_mm, reference[0]))
        <= tolerance
    {
        (reference[1], reference[0])
    } else {
        return false;
    };
    let delta = sub(end, start);
    let Some(direction) = normalize(delta, tolerance) else {
        return false;
    };
    let length = norm(delta);
    let mut prior = start;
    for edge in edges {
        if distance(prior, edge.start_mm) > tolerance {
            return false;
        }
        let step = sub(edge.end_mm, edge.start_mm);
        let projection = dot(step, direction);
        if projection <= 0.0 || norm(cross(step, direction)) > tolerance {
            return false;
        }
        for point in [edge.start_mm, edge.end_mm] {
            let along = dot(sub(point, start), direction);
            if along < -tolerance
                || along > length + tolerance
                || norm(cross(sub(point, start), direction)) > tolerance
            {
                return false;
            }
        }
        let line = match &edge.support {
            CurveSupport::Line(line) => *line,
            CurveSupport::BSpline(spline) => {
                let Some(line) = crate::source_edge_locus::certified_straight_spline(
                    spline,
                    edge.start_mm,
                    edge.end_mm,
                    tolerance,
                ) else {
                    return false;
                };
                line
            }
            _ => return false,
        };
        let Some(support_axis) = normalize(line.direction, tolerance) else {
            return false;
        };
        if norm(cross(support_axis, direction)) > tolerance {
            return false;
        }
        if [edge.start_mm, edge.end_mm]
            .into_iter()
            .any(|point| norm(cross(sub(point, line.origin_mm), support_axis)) > tolerance)
        {
            return false;
        }
        prior = edge.end_mm;
    }
    distance(prior, end) <= tolerance
}

/// The source surface is linear in its second parameter, so corresponding
/// row endpoints trace complete straight generator lines. Each generator
/// chain must trace exactly one of them, with no omissions or repeated edges.
pub(super) fn complete_ruled_face_trim(
    loops: &[brep::FaceLoop],
    rows: &[BSplineSupport; 2],
    cap_edge_ids: [u64; 2],
    tolerance: f64,
) -> bool {
    let [boundary] = loops else {
        return false;
    };
    let edges = &boundary.edges;
    if edges.len() < 4 || !tolerance.is_finite() || tolerance <= 0.0 {
        return false;
    }
    if edges
        .iter()
        .enumerate()
        .any(|(i, edge)| edge.end_vertex != edges[(i + 1) % edges.len()].start_vertex)
    {
        return false;
    }
    let indexes = cap_edge_ids.map(|id| edges.iter().position(|edge| edge.edge_id == id));
    let [Some(a), Some(b)] = indexes else {
        return false;
    };
    if a == b
        || !endpoints_match(&edges[a], &rows[0], tolerance)
        || !endpoints_match(&edges[b], &rows[1], tolerance)
    {
        return false;
    }
    let Some(first) = run_between(edges, a, b) else {
        return false;
    };
    let Some(second) = run_between(edges, b, a) else {
        return false;
    };
    let Some(generators) = rows[0]
        .control_points_mm
        .first()
        .copied()
        .zip(rows[1].control_points_mm.first().copied())
        .zip(
            rows[0]
                .control_points_mm
                .last()
                .copied()
                .zip(rows[1].control_points_mm.last().copied()),
        )
        .map(|((a, b), (c, d))| [(a, b).into(), (c, d).into()])
    else {
        return false;
    };
    (verify_generator_run(&first, generators[0], tolerance)
        && verify_generator_run(&second, generators[1], tolerance))
        || (verify_generator_run(&first, generators[1], tolerance)
            && verify_generator_run(&second, generators[0], tolerance))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::brep::{FaceLoop, LineSupport};

    fn vertex(z: f64, x: f64) -> [f64; 3] { [x, 0.0, z] }

    fn curve(z: f64) -> BSplineSupport {
        BSplineSupport {
            degree: 1,
            control_points_mm: vec![vertex(z, 0.0), vertex(z, 1.0)],
            knots: vec![0.0, 0.0, 1.0, 1.0],
            weights: None,
        }
    }

    fn edge(
        id: u64,
        from: u64,
        to: u64,
        start_mm: [f64; 3],
        end_mm: [f64; 3],
        support: CurveSupport,
    ) -> OrientedEdgeUse {
        OrientedEdgeUse {
            oriented_edge_id: id + 100,
            edge_id: id,
            curve_id: id + 200,
            curve_same_sense: true,
            parameter_forward: true,
            start_vertex: from,
            end_vertex: to,
            start_mm,
            end_mm,
            support,
        }
    }

    fn generator(id: u64, from: u64, to: u64, a: [f64; 3], b: [f64; 3]) -> OrientedEdgeUse {
        let direction = if b[2] > a[2] {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 0.0, -1.0]
        };
        edge(
            id,
            from,
            to,
            a,
            b,
            CurveSupport::Line(LineSupport {
                origin_mm: a,
                direction,
            }),
        )
    }

    fn complete_face() -> FaceLoop {
        let [a, b, c, d] = [
            vertex(0.0, 0.0),
            vertex(0.0, 1.0),
            vertex(1.0, 1.0),
            vertex(1.0, 0.0),
        ];
        FaceLoop {
            bound_id: 1,
            loop_id: 2,
            outer: true,
            orientation: true,
            edges: vec![
                edge(1, 1, 2, a, b, CurveSupport::BSpline(curve(0.0))),
                generator(2, 2, 5, b, vertex(0.4, 1.0)),
                generator(3, 5, 3, vertex(0.4, 1.0), c),
                edge(4, 3, 4, c, d, CurveSupport::BSpline(curve(1.0))),
                generator(5, 4, 6, d, vertex(0.6, 0.0)),
                generator(6, 6, 1, vertex(0.6, 0.0), a),
            ],
        }
    }

    #[test]
    fn complete_multisegment_generators_and_whole_caps_are_required() {
        let rows = [curve(0.0), curve(1.0)];
        let face = complete_face();
        assert!(complete_ruled_face_trim(
            std::slice::from_ref(&face),
            &rows,
            [1, 4],
            1.0e-7
        ));

        let mut shortened_cap = face.clone();
        shortened_cap.edges[0].end_mm[0] = 0.9;
        assert!(!complete_ruled_face_trim(
            &[shortened_cap],
            &rows,
            [1, 4],
            1.0e-7
        ));

        let mut unsupported_generator = face.clone();
        let CurveSupport::Line(line) = &mut unsupported_generator.edges[1].support else {
            unreachable!("test line support is fixed");
        };
        line.origin_mm[0] += 0.01;
        assert!(!complete_ruled_face_trim(
            &[unsupported_generator],
            &rows,
            [1, 4],
            1.0e-7
        ));

        let mut disconnected = face.clone();
        disconnected.edges[2].start_vertex += 100;
        assert!(!complete_ruled_face_trim(
            &[disconnected],
            &rows,
            [1, 4],
            1.0e-7
        ));

        let mut retraced = face.clone();
        retraced.edges[2].start_mm[2] = 1.2;
        assert!(!complete_ruled_face_trim(
            &[retraced],
            &rows,
            [1, 4],
            1.0e-7
        ));

        assert!(!complete_ruled_face_trim(
            &[face.clone(), face.clone()],
            &rows,
            [1, 4],
            1.0e-7
        ));
        assert!(!complete_ruled_face_trim(&[face], &rows, [1, 1], 1.0e-7));
    }
}
