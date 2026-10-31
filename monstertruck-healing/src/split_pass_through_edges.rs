use super::{
    BoundedCurve, CompressedEdgeUse, CompressedTrimmedShell, CurveParameter, Cut, Edge, EdgeIndex,
    HashMap, KERNEL_EPSILON, MetricSpace, Origin, ParameterDivision1D, ParametricCurve3D, Point3,
    SearchNearestParameter, Shell, SnapCurveEndpoints, Wire, append_unassigned_stable_ids,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use std::env;

#[derive(Clone)]
struct EdgeSplitInfo {
    splits: Vec<(f64, usize)>,
    indices: Vec<usize>,
}

fn debug_translated_focus_edge<C>(edge: &Edge<C>) -> bool
where C: ParametricCurve3D<Point = Point3> + BoundedCurve<Point = Point3> {
    let Ok(debug) = env::var("MT_BOOL_DEBUG_PASS_THROUGH") else {
        return false;
    };
    if debug != "translated" && debug != "origin" && debug != "all" {
        return false;
    }
    let curve = &edge.curve;
    let (t0, t1) = curve.range_tuple();
    let midpoint = curve.evaluate(f64::midpoint(t0, t1));
    if debug == "translated" {
        midpoint.x > 0.75
            && midpoint.x < 1.01
            && midpoint.y > 0.65
            && midpoint.y < 0.86
            && midpoint.z > -0.01
            && midpoint.z < 0.30
    } else if debug == "origin" {
        midpoint.x > -0.01
            && midpoint.x < 0.35
            && midpoint.y > 0.15
            && midpoint.y < 0.35
            && midpoint.z > -0.01
            && midpoint.z < 0.05
    } else {
        midpoint.x > 0.75
            && midpoint.x < 1.01
            && midpoint.y > 0.65
            && midpoint.y < 0.86
            && midpoint.z > -0.01
            && midpoint.z < 0.30
            || midpoint.x > -0.01
                && midpoint.x < 0.35
                && midpoint.y > 0.15
                && midpoint.y < 0.35
                && midpoint.z > -0.01
                && midpoint.z < 0.05
    }
}

fn sampled_curve_aabb<C>(edge: &Edge<C>, tol: f64) -> (Point3, Point3)
where C: ParametricCurve3D<Point = Point3>
        + BoundedCurve<Point = Point3>
        + ParameterDivision1D<Point = Point3> {
    let (_, samples) = edge.curve.parameter_division(
        edge.curve.range_tuple(),
        f64::max(
            tol,
            DEFAULT_NUMERICAL_TOLERANCES.healing_sampling_floor_factor * KERNEL_EPSILON,
        ),
    );
    samples.into_iter().fold(
        (
            Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
            Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
        ),
        |(min, max), point| {
            (
                Point3::new(min.x.min(point.x), min.y.min(point.y), min.z.min(point.z)),
                Point3::new(max.x.max(point.x), max.y.max(point.y), max.z.max(point.z)),
            )
        },
    )
}

fn split_parameters<C>(edge: &Edge<C>, vertices: &[Point3], tol: f64) -> Vec<(f64, usize)>
where C: ParametricCurve3D<Point = Point3>
        + BoundedCurve<Point = Point3>
        + ParameterDivision1D<Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3> {
    let debug_focus = debug_translated_focus_edge(edge);
    let tol2 = tol * tol;
    let (t0, t1) = edge.curve.range_tuple();
    let param_eps = f64::mul_add(
        (t1 - t0).abs(),
        DEFAULT_NUMERICAL_TOLERANCES.split_parameter_relative,
        DEFAULT_NUMERICAL_TOLERANCES.split_parameter_absolute,
    );
    let (min, max) = sampled_curve_aabb(edge, tol);
    let candidates: Vec<_> = vertices
        .iter()
        .copied()
        .enumerate()
        .filter(|(index, _)| *index != edge.vertices.0 && *index != edge.vertices.1)
        .filter(|(_, point)| {
            point.x >= min.x - tol
                && point.x <= max.x + tol
                && point.y >= min.y - tol
                && point.y <= max.y + tol
                && point.z >= min.z - tol
                && point.z <= max.z + tol
        })
        .filter_map(|(index, point)| {
            let t = edge.curve.search_nearest_parameter(point, None, 100)?;
            let distance2 = edge.curve.evaluate(t).distance2(point);
            Some((index, point, t, distance2))
        })
        .collect();
    if debug_focus {
        let curve = &edge.curve;
        let midpoint = curve.evaluate(f64::midpoint(t0, t1));
        tracing::debug!(
            "debug pass-through edge front={:?} back={:?} midpoint={midpoint:?} candidates={:?}",
            vertices[edge.vertices.0],
            vertices[edge.vertices.1],
            candidates,
        );
    }
    let mut splits: Vec<_> = candidates
        .into_iter()
        .filter_map(|(index, point, t, distance2)| {
            let on_curve = distance2 <= tol2;
            if debug_focus {
                tracing::debug!(
                    "debug pass-through candidate index={index} point={point:?} t={t} distance2={distance2} on_curve={on_curve}",
                );
            }
            (on_curve && t > t0 + param_eps && t < t1 - param_eps).then_some((t, index))
        })
        .collect();
    splits.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));
    splits.dedup_by(|lhs, rhs| {
        let same_param = (lhs.0 - rhs.0).abs() <= param_eps;
        let same_vertex = lhs.1 == rhs.1;
        same_param || same_vertex
    });
    if debug_focus {
        tracing::debug!("debug pass-through splits={splits:?}");
    }
    splits
}

fn split_edge<C>(
    edge_index: usize,
    edges: &mut Vec<Edge<C>>,
    vertices: &[Point3],
    splits: &[(f64, usize)],
) -> Vec<usize>
where
    C: BoundedCurve<Point = Point3> + Cut + SnapCurveEndpoints,
{
    let edge = &mut edges[edge_index];
    let back = edge.vertices.1;
    let first = splits[0];
    let mut tail = edge.curve.cut(first.0);
    edge.vertices.1 = first.1;
    edge.curve
        .snap_endpoints(vertices[edge.vertices.0], vertices[first.1]);
    let mut previous = first.1;
    let mut indices = vec![edge_index];
    splits.iter().skip(1).for_each(|(t, vertex)| {
        let next_tail = tail.cut(*t);
        tail.snap_endpoints(vertices[previous], vertices[*vertex]);
        edges.push(Edge {
            vertices: (previous, *vertex),
            curve: tail.clone(),
        });
        indices.push(edges.len() - 1);
        tail = next_tail;
        previous = *vertex;
    });
    tail.snap_endpoints(vertices[previous], vertices[back]);
    edges.push(Edge {
        vertices: (previous, back),
        curve: tail,
    });
    indices.push(edges.len() - 1);
    indices
}

fn replace_edge_uses(wire: &mut Wire, replacements: &HashMap<usize, EdgeSplitInfo>) {
    *wire = wire
        .iter()
        .copied()
        .flat_map(|edge| {
            replacements.get(&edge.index).map_or_else(
                || vec![edge],
                |split_info| {
                    if edge.orientation {
                        split_info
                            .indices
                            .iter()
                            .copied()
                            .map(|index| EdgeIndex {
                                index,
                                orientation: true,
                            })
                            .collect::<Vec<_>>()
                    } else {
                        split_info
                            .indices
                            .iter()
                            .rev()
                            .copied()
                            .map(|index| EdgeIndex {
                                index,
                                orientation: false,
                            })
                            .collect::<Vec<_>>()
                    }
                },
            )
        })
        .collect();
}

fn split_trim_curve<T>(
    trim_curve: T,
    edge_range: (f64, f64),
    orientation: bool,
    splits: &[(f64, usize)],
) -> Vec<T>
where
    T: BoundedCurve + Cut + Clone,
{
    let (edge_t0, edge_t1) = edge_range;
    let edge_len = edge_t1 - edge_t0;
    let (trim_t0, trim_t1) = trim_curve.range_tuple();
    let trim_scale = trim_t1 - trim_t0;
    if trim_scale.so_small() {
        return vec![trim_curve];
    }
    let param_eps = f64::mul_add(
        edge_len.abs(),
        DEFAULT_NUMERICAL_TOLERANCES.split_parameter_relative,
        DEFAULT_NUMERICAL_TOLERANCES.split_parameter_absolute,
    );
    let trim_eps = f64::mul_add(
        trim_scale.abs(),
        DEFAULT_NUMERICAL_TOLERANCES.split_parameter_relative,
        DEFAULT_NUMERICAL_TOLERANCES.split_parameter_absolute,
    );
    let mut local_splits: Vec<_> = splits
        .iter()
        .map(|(t, vertex)| {
            let alpha = if edge_len.so_small() {
                0.5
            } else {
                (*t - edge_t0) / edge_len
            };
            let local_alpha = if orientation { alpha } else { 1.0 - alpha };
            (f64::mul_add(trim_scale, local_alpha, trim_t0), *vertex)
        })
        .filter(|(t, _)| *t > trim_t0 + trim_eps && *t < trim_t1 - trim_eps)
        .collect();
    local_splits.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));
    local_splits.dedup_by(|lhs, rhs| {
        let same_param = (lhs.0 - rhs.0).abs() <= param_eps;
        let same_vertex = lhs.1 == rhs.1;
        same_param || same_vertex
    });
    let mut head = trim_curve;
    let mut segments = Vec::with_capacity(local_splits.len() + 1);
    for (t, _) in local_splits {
        let (head_t0, head_t1) = head.range_tuple();
        if t > head_t0 + trim_eps && t < head_t1 - trim_eps {
            let tail = head.cut(t);
            segments.push(head.clone());
            head = tail;
        }
    }
    segments.push(head);
    segments
}

fn replace_trimmed_edge_uses<T, C>(
    wire: &mut Vec<CompressedEdgeUse<T>>,
    edges: &[Edge<C>],
    replacements: &HashMap<usize, EdgeSplitInfo>,
) where
    T: BoundedCurve + Cut + Clone,
    C: BoundedCurve<Point = Point3>,
{
    *wire = wire
        .iter()
        .cloned()
        .flat_map(|edge_use| {
            replacements.get(&edge_use.index).map_or_else(
                || vec![edge_use.clone()],
                |split_info| {
                    let trim_segments = edge_use.trim_curve.clone().map(|trim_curve| {
                        split_trim_curve(
                            trim_curve,
                            edges[edge_use.index].curve.range_tuple(),
                            edge_use.orientation,
                            &split_info.splits,
                        )
                    });
                    let indices: Vec<_> = if edge_use.orientation {
                        split_info.indices.clone()
                    } else {
                        split_info.indices.iter().rev().copied().collect()
                    };
                    let trim_values: Vec<_> = trim_segments.map_or_else(
                        || vec![None; indices.len()],
                        |segments| segments.into_iter().map(Some).collect(),
                    );
                    indices
                        .into_iter()
                        .zip(trim_values)
                        .map(|(index, trim_curve)| CompressedEdgeUse {
                            index,
                            orientation: edge_use.orientation,
                            trim_curve,
                        })
                        .collect::<Vec<_>>()
                },
            )
        })
        .collect();
}

pub fn split_pass_through_edges<C, S>(shell: &mut Shell<Point3, C, S>, tol: f64)
where C: ParametricCurve3D<Point = Point3>
        + BoundedCurve<Point = Point3>
        + Cut
        + ParameterDivision1D<Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3>
        + SnapCurveEndpoints
        + Clone {
    let len = shell.edges.len();
    let replacements: HashMap<_, _> = (0..len)
        .filter_map(|edge_index| {
            let splits = split_parameters(&shell.edges[edge_index], &shell.vertices, tol);
            (!splits.is_empty()).then(|| {
                (
                    edge_index,
                    EdgeSplitInfo {
                        indices: split_edge(edge_index, &mut shell.edges, &shell.vertices, &splits),
                        splits,
                    },
                )
            })
        })
        .collect();
    append_unassigned_stable_ids(&mut shell.edge_stable_ids, shell.edges.len() - len);
    if replacements.is_empty() {
        return;
    }
    shell.faces.iter_mut().for_each(|face| {
        face.boundaries
            .iter_mut()
            .for_each(|wire| replace_edge_uses(wire, &replacements));
    });
}

pub fn split_pass_through_edges_trimmed<C, S, T>(
    shell: &mut CompressedTrimmedShell<Point3, C, S, T>,
    tol: f64,
) where
    C: ParametricCurve3D<Point = Point3>
        + BoundedCurve<Point = Point3>
        + Cut
        + ParameterDivision1D<Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3>
        + SnapCurveEndpoints
        + Clone,
    T: BoundedCurve + Cut + Clone,
{
    let len = shell.edges.len();
    let replacements: HashMap<_, _> = (0..len)
        .filter_map(|edge_index| {
            let splits = split_parameters(&shell.edges[edge_index], &shell.vertices, tol);
            (!splits.is_empty()).then(|| {
                (
                    edge_index,
                    EdgeSplitInfo {
                        indices: split_edge(edge_index, &mut shell.edges, &shell.vertices, &splits),
                        splits,
                    },
                )
            })
        })
        .collect();
    if replacements.is_empty() {
        return;
    }
    shell.faces.iter_mut().for_each(|face| {
        face.boundaries
            .iter_mut()
            .for_each(|wire| replace_trimmed_edge_uses(wire, &shell.edges, &replacements));
    });
}

/// Unify coincident duplicate edges after the pass-through split.
///
/// When both sides of a shared span get subdivided (each face's copy split
/// at the same passing-through vertices), the shell holds several edges with
/// the SAME vertex pair and geometry under different indices, and closedness
/// counting still sees every one as half-used (sphere-minus-cube, cell 8:
/// seven duplicate groups after the split). Only a SOUND pairing is unified:
/// each side used exactly once, with the two uses traversing the span in
/// OPPOSITE directions after vertex-order adjustment (adjacent faces across
/// the span). Same-direction coincident edges are legitimate separate
/// topology (e.g. doubled seams at tangency) and must survive -- unifying
/// them deterministically breaks the `or()` union of offset cubes
/// (`or_outputs_boolean_named_selections`).
pub fn dedup_coincident_pass_through_edges<C>(
    edges: &[Edge<C>],
    use_orientations: &[Vec<(usize, bool)>],
    tol: f64,
) -> HashMap<usize, (usize, bool)>
where
    C: ParametricCurve3D<Point = Point3> + BoundedCurve<Point = Point3>,
{
    let tol2 = tol * tol;
    let midpoint = |edge: &Edge<C>| {
        let (t0, t1) = edge.curve.range_tuple();
        edge.curve.evaluate(f64::midpoint(t0, t1))
    };
    let mut canonical_by_pair: HashMap<(usize, usize), Vec<usize>> = HashMap::default();
    let mut remap: HashMap<usize, (usize, bool)> = HashMap::default();
    for (index, edge) in edges.iter().enumerate() {
        let (a, b) = edge.vertices;
        if a == b {
            continue;
        }
        let key = (a.min(b), a.max(b));
        let mid = midpoint(edge);
        let bucket = canonical_by_pair.entry(key).or_default();
        let matched = bucket.iter().copied().find(|&canonical_index| {
            if midpoint(&edges[canonical_index]).distance2(mid) > tol2 {
                return false;
            }
            // Only unify a SOUND pairing: each side used exactly once, and
            // the two uses traverse the span in OPPOSITE directions after
            // vertex-order adjustment (adjacent faces across the span).
            // Same-direction coincident edges are legitimate separate
            // topology (e.g. doubled seams) and must survive.
            let canonical_uses = use_orientations
                .get(canonical_index)
                .map_or(&[][..], Vec::as_slice);
            let duplicate_uses = use_orientations.get(index).map_or(&[][..], Vec::as_slice);
            let (
                &[(canonical_face, canonical_orientation)],
                &[(duplicate_face, duplicate_orientation)],
            ) = (canonical_uses, duplicate_uses)
            else {
                return false;
            };
            // Unifying two uses within ONE face wire makes that wire
            // non-simple (TrimmedShell rejects it); only pair across faces.
            if canonical_face == duplicate_face {
                return false;
            }
            let reversed = edges[canonical_index].vertices != edges[index].vertices;
            let duplicate_effective = duplicate_orientation ^ reversed;
            duplicate_effective != canonical_orientation
        });
        match matched {
            Some(canonical_index) => {
                let reversed = edges[canonical_index].vertices != edge.vertices;
                if std::env::var("MT_BOOL_DEBUG_PASS_THROUGH_DEDUP").is_ok() {
                    let m = midpoint(edge);
                    tracing::debug!(
                        "debug pass-through-dedup unify dup={index} -> canon={canonical_index} reversed={reversed} pair={key:?} mid=({:.4},{:.4},{:.4})",
                        m.x,
                        m.y,
                        m.z,
                    );
                }
                remap.insert(index, (canonical_index, reversed));
            }
            None => bucket.push(index),
        }
    }
    remap
}

/// Apply a [`dedup_coincident_pass_through_edges`] remap to trimmed wires.
pub fn remap_trimmed_edge_uses_for_dedup<T>(
    wire: &mut [CompressedEdgeUse<T>],
    remap: &HashMap<usize, (usize, bool)>,
) {
    for edge_use in wire {
        if let Some((canonical, reversed)) = remap.get(&edge_use.index) {
            edge_use.index = *canonical;
            if *reversed {
                edge_use.orientation = !edge_use.orientation;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_geometry::prelude::Line;
    use monstertruck_topology::{StableId, compress::CompressedFace as Face};

    #[test]
    fn pass_through_split_preserves_existing_stable_ids_and_marks_new_edge_unassigned() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let p2 = Point3::new(2.0, 0.0, 0.0);
        let mut shell = Shell {
            vertices: vec![p0, p1, p2],
            edges: vec![Edge {
                vertices: (0, 2),
                curve: Line(p0, p2),
            }],
            faces: vec![Face {
                boundaries: vec![vec![EdgeIndex {
                    index: 0,
                    orientation: true,
                }]],
                orientation: true,
                surface: (),
            }],
            vertex_stable_ids: Some(vec![
                StableId::new(11),
                StableId::new(12),
                StableId::new(13),
            ]),
            edge_stable_ids: Some(vec![StableId::new(21)]),
            face_stable_ids: Some(vec![StableId::new(31)]),
        };

        split_pass_through_edges(&mut shell, 1.0e-6);

        assert_eq!(shell.edges.len(), 2);
        assert_eq!(shell.edges[0].vertices, (0, 1));
        assert_eq!(shell.edges[1].vertices, (1, 2));
        assert_eq!(shell.faces[0].boundaries[0].len(), 2);
        assert_eq!(
            shell.edge_stable_ids,
            Some(vec![StableId::new(21), StableId::UNASSIGNED])
        );
        assert_eq!(
            shell.vertex_stable_ids,
            Some(vec![
                StableId::new(11),
                StableId::new(12),
                StableId::new(13)
            ])
        );
        assert_eq!(shell.face_stable_ids, Some(vec![StableId::new(31)]));
    }
}
