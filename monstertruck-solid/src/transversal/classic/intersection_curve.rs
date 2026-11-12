//! Marching intersection-curve backend, descended from the published 0.3.2 kernel.
//!
//! It still extracts mesh-vs-mesh interference segments and re-samples them onto
//! both carrier surfaces, but endpoint stitching is now fail-closed around branch
//! ambiguity and uses explicit sampled validation for short linear gaps. The old
//! parameter-space wrapper remains intentionally absent; the loops store consumes
//! only the cleaned 3-D [`IntersectionCurve`] leader.

use monstertruck_core::{
    cgmath64::{EuclideanSpace, InnerSpace, MetricSpace, Point3},
    tolerance::{DEFAULT_NUMERICAL_TOLERANCES, KERNEL_EPSILON},
};
use monstertruck_geometry::prelude::*;
use monstertruck_meshing::prelude::*;
use num_traits::ToPrimitive;

use crate::transversal::polyline_construction::{
    construct_polylines, enumerate_endpoint_joins, merge_endpoint_join,
    merge_endpoint_join_with_bridge, remove_redundant_interference_strands,
    stitch_nearby_polylines,
};

type Polyline = PolylineCurve<Point3>;

struct FragmentValidation<'a, S> {
    surface0: &'a S,
    polygon0: &'a PolygonMesh,
    surface1: &'a S,
    polygon1: &'a PolygonMesh,
    tolerance: f64,
}

fn endpoint_and_inner(
    polyline: &Polyline,
    front: bool,
    minimum_span: f64,
) -> Option<(Point3, Point3)> {
    if polyline.len() < 2 {
        return None;
    }
    let endpoint = if front {
        polyline[0]
    } else {
        *polyline.last()?
    };
    let inner = if front {
        polyline
            .iter()
            .copied()
            .skip(1)
            .find(|point| endpoint.distance(*point) > minimum_span)?
    } else {
        polyline
            .iter()
            .copied()
            .rev()
            .skip(1)
            .find(|point| endpoint.distance(*point) > minimum_span)?
    };
    Some((endpoint, inner))
}

/// Track the exact UV branch along adjacent samples of one short SSI join.
/// Each candidate is still independently checked on the carrier; the hint
/// skips repeated global NURBS presearch when consecutive samples belong to
/// the same smooth trim. Fall back to the original unhinted projection whenever
/// local continuation cannot establish point-on-surface incidence.
fn point_is_on_surface<S>(
    surface: &S,
    point: Point3,
    tolerance: f64,
    hint: &mut Option<(f64, f64)>,
) -> bool
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let within_tolerance = |parameter: (f64, f64)| {
        surface.evaluate(parameter.0, parameter.1).distance(point) <= tolerance
    };
    let local = hint.and_then(|previous| {
        surface
            .search_parameter(point, Some(previous), 100)
            .filter(|parameter| within_tolerance(*parameter))
            .or_else(|| {
                surface
                    .search_nearest_parameter(point, Some(previous), 100)
                    .filter(|parameter| within_tolerance(*parameter))
            })
    });
    let parameter = local.or_else(|| {
        surface
            .search_parameter(point, None, 100)
            .filter(|parameter| within_tolerance(*parameter))
            .or_else(|| {
                surface
                    .search_nearest_parameter(point, None, 100)
                    .filter(|parameter| within_tolerance(*parameter))
            })
    });
    if let Some(parameter) = parameter {
        *hint = Some(parameter);
        true
    } else {
        false
    }
}

fn validated_linear_fragment_join<S>(
    first: &Polyline,
    first_front: bool,
    second: &Polyline,
    second_front: bool,
    validation: &FragmentValidation<'_, S>,
) -> bool
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let Some((a, a_inner)) = endpoint_and_inner(first, first_front, validation.tolerance) else {
        return false;
    };
    let Some((b, b_inner)) = endpoint_and_inner(second, second_front, validation.tolerance) else {
        return false;
    };

    let bridge = b - a;
    let bridge_length = bridge.magnitude();

    let outward_a = a - a_inner;
    let outward_b = b - b_inner;
    let span_a = outward_a.magnitude();
    let span_b = outward_b.magnitude();
    if span_a <= validation.tolerance || span_b <= validation.tolerance {
        return false;
    }

    // The missing interval must continue outward from both fragments and remain
    // on their endpoint tangents to within the kernel validation tolerance.
    let direction_ok = outward_a.dot(bridge) > 0.0 && outward_b.dot(-bridge) > 0.0;
    let line_error_a = outward_a.cross(bridge).magnitude() / span_a;
    let line_error_b = outward_b.cross(bridge).magnitude() / span_b;
    let line_ok = line_error_a <= validation.tolerance && line_error_b <= validation.tolerance;
    if !(direction_ok && line_ok) {
        return false;
    }

    // Validate the entire short bridge, not just its endpoints. Dense samples must
    // stay on both exact carrier surfaces and inside both trimmed polygon faces.
    // Refuse pathological gaps rather than spending unbounded work validating them.
    let intervals = (bridge_length / (2.0 * validation.tolerance))
        .ceil()
        .max(2.0);
    if !intervals.is_finite() || intervals > 4096.0 {
        return false;
    }
    let samples: Vec<_> = (0_u32..=4096)
        .map(f64::from)
        .take_while(|&index| index <= intervals)
        .map(|index| a + bridge * (index / intervals))
        .collect();
    // The carrier polygons are tessellations, not proof of exact incidence,
    // but their trim bounds can cheaply reject impossible candidates *before*
    // the costly global NURBS inverse search. This only reorders two required
    // predicates, preserving the exact admission criteria.
    let polygon_tolerance = 2.0 * validation.tolerance;
    if !validation
        .polygon0
        .neighborhood_include(&samples, polygon_tolerance)
        || !validation
            .polygon1
            .neighborhood_include(&samples, polygon_tolerance)
    {
        return false;
    }
    let mut hint0 = None;
    let mut hint1 = None;
    samples.iter().copied().all(|point| {
        point_is_on_surface(validation.surface0, point, validation.tolerance, &mut hint0)
            && point_is_on_surface(validation.surface1, point, validation.tolerance, &mut hint1)
    })
}

fn stitch_validated_linear_intersection_fragments<S>(
    mut polylines: Vec<Polyline>,
    surface0: &S,
    polygon0: &PolygonMesh,
    surface1: &S,
    polygon1: &PolygonMesh,
    tolerance: f64,
) -> Vec<Polyline>
where
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let validation = FragmentValidation {
        surface0,
        polygon0,
        surface1,
        polygon1,
        tolerance,
    };
    loop {
        let Some(mut joins) =
            enumerate_endpoint_joins(&polylines, |first, first_front, second, second_front| {
                validated_linear_fragment_join(
                    first,
                    first_front,
                    second,
                    second_front,
                    &validation,
                )
            })
        else {
            return polylines;
        };

        joins.sort_by(|left, right| {
            left.distance2
                .total_cmp(&right.distance2)
                .then_with(|| left.first.chain.cmp(&right.first.chain))
                .then_with(|| left.second.chain.cmp(&right.second.chain))
        });
        let Some(join) = joins.first().copied() else {
            return polylines;
        };

        // Several fragments may lie on one validated straight intersection, so an
        // endpoint can legitimately see more than one farther mate. Continue to
        // the closest one. A near-equal alternative sharing either endpoint is a
        // real branch ambiguity and remains fail-closed.
        let join_distance = join.distance2.sqrt();
        let ambiguous = joins.iter().skip(1).any(|other| {
            let join_endpoints = [join.first, join.second];
            let shares_endpoint =
                join_endpoints.contains(&other.first) || join_endpoints.contains(&other.second);
            let distance_delta = (other.distance2.sqrt() - join_distance).abs();
            shares_endpoint && distance_delta <= tolerance
        });
        if ambiguous {
            return polylines;
        }

        if merge_endpoint_join(&mut polylines, join).is_none() {
            return polylines;
        }
    }
}

/// Correct a short curved SSI gap using transverse slicing planes.
/// Both endpoints come from certified fragments, and intermediate predictor
/// points are corrected by the exact two-surface/one-plane Newton problem.
/// This never treats a straight chord as the true intersection geometry.
fn corrected_curved_fragment_bridge<S>(
    first: &Polyline,
    first_front: bool,
    second: &Polyline,
    second_front: bool,
    surface0: &S,
    surface1: &S,
    tolerance: f64,
) -> Option<Vec<Point3>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let (a, a_inner) = endpoint_and_inner(first, first_front, tolerance)?;
    let (b, b_inner) = endpoint_and_inner(second, second_front, tolerance)?;
    let span = b - a;
    let length = span.magnitude();
    let mesh_scale = DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor;
    if !length.is_finite() || length <= tolerance || length > mesh_scale * 16.0 {
        return None;
    }
    let tangent_a = a - a_inner;
    let tangent_b = b - b_inner;
    let ta = tangent_a.magnitude();
    let tb = tangent_b.magnitude();
    if ta <= tolerance || tb <= tolerance {
        return None;
    }
    let direction = span / length;
    // Both samples must approach this missing interval from opposite ends;
    // require positive transverse progress, not identical straight tangents.
    if tangent_a.dot(direction) / ta <= 0.5 || -tangent_b.dot(direction) / tb <= 0.5 {
        return None;
    }

    let leader = PolylineCurve(vec![a, b]);
    let intersection = IntersectionCurve::new(surface0.clone(), surface1.clone(), leader);
    let debug = std::env::var_os("MT_CURVED_SSI_TRACE").is_some();
    let seed = |surface: &S, point: Point3| {
        surface
            .search_parameter(point, None, 100)
            .or_else(|| surface.search_nearest_parameter(point, None, 100))
    };
    let mut hint0 = seed(surface0, a).or_else(|| {
        if debug {
            eprintln!("CURVED_BRIDGE_FAIL seed_surface0 at={a:?}");
        }
        None
    })?;
    let mut hint1 = seed(surface1, a).or_else(|| {
        if debug {
            eprintln!("CURVED_BRIDGE_FAIL seed_surface1 at={a:?}");
        }
        None
    })?;
    let reference_tol =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_surface_projection_factor * KERNEL_EPSILON;
    let intervals = (length / (mesh_scale * 0.5)).ceil().clamp(2.0, 128.0) as usize;
    let mut samples = Vec::with_capacity(intervals + 1);
    samples.push(a);
    for index in 1..intervals {
        let fraction = index as f64 / intervals as f64;
        let predictor = a + span * fraction;
        let (q, uv0, uv1) = intersection.search_intersection_in_plane(
            predictor,
            direction,
            Some(hint0),
            Some(hint1),
            100,
        ).or_else(||{
            if debug {eprintln!("CURVED_BRIDGE_FAIL plane_corrector step={index}/{intervals} a={a:?} b={b:?} predictor={predictor:?} hint0={hint0:?} hint1={hint1:?}");}
            None
        })?;
        let s0 = surface0.evaluate(uv0.x, uv0.y);
        let s1 = surface1.evaluate(uv1.x, uv1.y);
        if debug {
            eprintln!(
                "CURVED_BRIDGE_SAMPLE step={index}/{intervals} q={q:?} shift={} residual={}",
                q.distance(predictor),
                s0.distance(s1)
            );
        }
        if s0.distance(s1) > reference_tol
            || q.distance(predictor) > length * 0.5
            || !q.x.is_finite()
            || !q.y.is_finite()
            || !q.z.is_finite()
            || (q - a).dot(direction) < 0.0
            || (b - q).dot(direction) < 0.0
        {
            return None;
        }
        if let Some(previous) = samples.last() {
            if (*previous - q).dot(direction) >= 0.0 {
                return None;
            }
        }
        hint0 = (uv0.x, uv0.y);
        hint1 = (uv1.x, uv1.y);
        samples.push(q);
    }
    samples.push(b);
    // Independent endpoint tangent checks detect a solver jumping to another
    // intersection branch near a tangency or self-crossing.
    let start = samples[1] - a;
    let end = b - samples[samples.len() - 2];
    if debug {
        eprintln!(
            "CURVED_BRIDGE_END length={length} start_alignment={} end_alignment={}",
            start.dot(tangent_a) / (start.magnitude() * ta),
            -end.dot(tangent_b) / (end.magnitude() * tb)
        );
    }
    if start.dot(tangent_a) / (start.magnitude() * ta) <= 0.5
        || -end.dot(tangent_b) / (end.magnitude() * tb) <= 0.5
    {
        return None;
    }
    Some(samples)
}

/// Reconstruct short gaps on curved SSI carriers. Each accepted bridge
/// includes its exact corrected interior samples, rather than a midpoint
/// snap or an unsupported straight chord.
fn stitch_validated_curved_intersection_fragments<S>(
    mut polylines: Vec<Polyline>,
    surface0: &S,
    surface1: &S,
    tolerance: f64,
) -> Vec<Polyline>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    loop {
        let Some(mut joins) =
            enumerate_endpoint_joins(&polylines, |first, first_front, second, second_front| {
                corrected_curved_fragment_bridge(
                    first,
                    first_front,
                    second,
                    second_front,
                    surface0,
                    surface1,
                    tolerance,
                )
                .is_some()
            })
        else {
            return polylines;
        };
        joins.sort_by(|a, b| {
            a.distance2
                .total_cmp(&b.distance2)
                .then_with(|| a.first.chain.cmp(&b.first.chain))
                .then_with(|| a.second.chain.cmp(&b.second.chain))
        });
        let Some(join) = joins.first().copied() else {
            return polylines;
        };
        let closest = join.distance2.sqrt();
        if joins.iter().skip(1).any(|other| {
            let same = [join.first, join.second];
            (same.contains(&other.first) || same.contains(&other.second))
                && (other.distance2.sqrt() - closest).abs() <= tolerance
        }) {
            return polylines;
        }
        let Some(bridge) = corrected_curved_fragment_bridge(
            &polylines[join.first.chain],
            join.first.front,
            &polylines[join.second.chain],
            join.second.front,
            surface0,
            surface1,
            tolerance,
        ) else {
            return polylines;
        };
        if std::env::var_os("MT_CURVED_SSI_TRACE").is_some() {
            eprintln!(
                "CURVED_SSI_JOIN chains_before={} a={} b={} gap={} bridge_points={}",
                polylines.len(),
                join.first.chain,
                join.second.chain,
                join.distance2.sqrt(),
                bridge.len()
            );
        }
        if merge_endpoint_join_with_bridge(&mut polylines, join, &bridge, tolerance).is_none() {
            return polylines;
        }
    }
}

/// Locate a certified SSI point along the *interior* of an existing
/// already-corrected intersection fragment. A near 3D chord is only a
/// broad-phase candidate: the exact surface/surface/slicing-plane Newton
/// solution must coincide with the point and remain on this local branch.
fn locate_on_exact_fragment<S>(
    curve: &Polyline,
    point: Point3,
    surface0: &S,
    surface1: &S,
    tolerance: f64,
) -> Option<(usize, f64)>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let mut best = None::<(f64, usize, f64)>;
    let mut matched_progress = None::<f64>;
    for (index, segment) in curve.windows(2).enumerate() {
        let delta = segment[1] - segment[0];
        let length2 = delta.magnitude2();
        if length2 <= tolerance * tolerance {
            continue;
        }
        let along = (point - segment[0]).dot(delta) / length2;
        if !(-1.0e-8..=1.0 + 1.0e-8).contains(&along) {
            continue;
        }
        let foot = segment[0] + delta * along.clamp(0.0, 1.0);
        let distance = foot.distance(point);
        if distance > DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor {
            continue;
        }
        let projection = IntersectionCurve::new(
            surface0.clone(),
            surface1.clone(),
            PolylineCurve(vec![segment[0], segment[1]]),
        );
        let hint = |surface: &S| {
            surface
                .search_parameter(segment[0], None, 100)
                .or_else(|| surface.search_nearest_parameter(segment[0], None, 100))
        };
        let Some((q, uv0, uv1)) = projection.search_intersection_in_plane(
            point,
            delta,
            hint(surface0),
            hint(surface1),
            100,
        ) else {
            continue;
        };
        let residual = surface0
            .evaluate(uv0.x, uv0.y)
            .distance(surface1.evaluate(uv1.x, uv1.y));
        if residual > tolerance || q.distance(point) > tolerance {
            continue;
        }
        let progress = index as f64 + along.clamp(0.0, 1.0);
        // Nonadjacent repeated solutions mean self-overlap, not a unique splice.
        if matched_progress.is_some_and(|prior| (progress - prior).abs() > 1.5) {
            return None;
        }
        matched_progress = Some(progress);
        if best.is_none_or(|(old, _, _)| distance < old) {
            best = Some((distance, index, along.clamp(0.0, 1.0)));
        }
    }
    best.map(|(_, index, along)| (index, along))
}

/// Find two fragments that trace a common local subarc in the same direction.
/// The first fragment overruns the second's start; discard *only* that
/// certified duplicate suffix. No new geometric surface tolerance is needed:
/// candidates are corrected on both exact surfaces and sampled through the
/// overlap, not guessed from a coincident endpoint.
fn overlap_join_candidate<S>(
    first: &Polyline,
    second: &Polyline,
    surface0: &S,
    surface1: &S,
    tolerance: f64,
) -> Option<Vec<Point3>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let first_start = *first.first()?;
    let second_start = *second.first()?;
    let first_end = *first.last()?;
    let second_end = *second.last()?;
    let (cut_index, cut_fraction) =
        locate_on_exact_fragment(first, second_start, surface0, surface1, tolerance)?;
    let (overlap_index, overlap_fraction) =
        locate_on_exact_fragment(second, first_end, surface0, surface1, tolerance)?;
    if cut_index == 0 && cut_fraction <= tolerance
        || cut_index + 1 >= first.len()
        || overlap_index + 1 >= second.len()
        || first_start.distance(second_start) <= tolerance
        || first_end.distance(second_end) <= tolerance
    {
        return None;
    }
    let (Some(forward0), Some(forward1)) = (
        first
            .get(cut_index + 1)
            .map(|point| point - first[cut_index]),
        second.get(1).map(|point| point - second[0]),
    ) else {
        return None;
    };
    let n0 = forward0.magnitude();
    let n1 = forward1.magnitude();
    if n0 <= tolerance || n1 <= tolerance || forward0.dot(forward1) / (n0 * n1) < 0.8 {
        return None;
    }
    // The discarded first suffix must progress along the retained second
    // fragment, not cross it or return along a distinct branch.
    let mut previous_progress = 0.0_f64;
    let mut trace = Vec::with_capacity(first.len() - cut_index);
    trace.push(second_start);
    trace.extend_from_slice(&first[cut_index + 1..]);
    let mesh_step = DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor * 0.5;
    let mut samples = 0usize;
    let ic = IntersectionCurve::new(
        surface0.clone(),
        surface1.clone(),
        PolylineCurve(vec![second_start, first_end]),
    );
    for pair in trace.windows(2) {
        let delta = pair[1] - pair[0];
        let count = (delta.magnitude() / mesh_step).ceil().max(1.0) as usize;
        if samples + count > 4096 {
            return None;
        }
        for step in 1..=count {
            let fraction = step as f64 / count as f64;
            let point = if step == count {
                pair[1]
            } else {
                let predictor = pair[0] + delta * fraction;
                let seed = |surface: &S| {
                    surface
                        .search_parameter(pair[0], None, 100)
                        .or_else(|| surface.search_nearest_parameter(pair[0], None, 100))
                };
                let (exact, u0, u1) = ic.search_intersection_in_plane(
                    predictor,
                    delta,
                    seed(surface0),
                    seed(surface1),
                    100,
                )?;
                if surface0
                    .evaluate(u0.x, u0.y)
                    .distance(surface1.evaluate(u1.x, u1.y))
                    > tolerance
                    || exact.distance(predictor) > delta.magnitude() * 0.5
                {
                    return None;
                }
                exact
            };
            let (index, local) =
                locate_on_exact_fragment(second, point, surface0, surface1, tolerance)?;
            let progress = index as f64 + local;
            if progress + 1.0e-8 < previous_progress {
                return None;
            }
            previous_progress = progress;
            samples += 1;
        }
    }
    if previous_progress < overlap_index as f64 + overlap_fraction - 1.0e-5 {
        return None;
    }
    let prefix = PolylineCurve(first[..=cut_index].to_vec());
    let bridge = corrected_curved_fragment_bridge(
        &prefix, false, second, true, surface0, surface1, tolerance,
    )?;
    let mut result = prefix.0;
    result.extend_from_slice(&bridge[1..bridge.len() - 1]);
    result.extend(second.iter().copied());
    Some(result)
}

/// Splice only unique, same-direction exact overlap matches; an ambiguous
/// repeated locus is left unmodified for upstream fail-closed handling.
fn stitch_overlapping_intersection_fragments<S>(
    mut chains: Vec<Polyline>,
    surface0: &S,
    surface1: &S,
    tolerance: f64,
) -> Vec<Polyline>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    loop {
        let mut matches = Vec::new();
        for first in 0..chains.len() {
            for second in 0..chains.len() {
                if first == second {
                    continue;
                }
                for reverse_first in [false, true] {
                    for reverse_second in [false, true] {
                        let mut a = chains[first].clone();
                        let mut b = chains[second].clone();
                        if reverse_first {
                            a.0.reverse();
                        }
                        if reverse_second {
                            b.0.reverse();
                        }
                        if let Some(merged) =
                            overlap_join_candidate(&a, &b, surface0, surface1, tolerance)
                        {
                            // A and B can admit the same union in reversed order.
                            // Use canonical 3D endpoint order to admit one
                            // orientation while retaining genuine ambiguity.
                            let front = *merged.first().unwrap();
                            let back = *merged.last().unwrap();
                            if front
                                .x
                                .total_cmp(&back.x)
                                .then_with(|| front.y.total_cmp(&back.y))
                                .then_with(|| front.z.total_cmp(&back.z))
                                .is_lt()
                            {
                                matches.push((
                                    first,
                                    second,
                                    reverse_first,
                                    reverse_second,
                                    merged,
                                ));
                            }
                        }
                    }
                }
            }
        }
        if matches.len() != 1 {
            return chains;
        }
        let (first, second, _, _, merged) = matches.remove(0);
        if std::env::var_os("MT_CURVED_SSI_TRACE").is_some() {
            eprintln!(
                "SSI_OVERLAP_SPLICE before={} first={} second={} result_points={}",
                chains.len(),
                first,
                second,
                merged.len()
            );
        }
        for index in [first.max(second), first.min(second)] {
            chains.remove(index);
        }
        chains.push(PolylineCurve(merged));
    }
}

/// Re-sample a raw interference polyline onto both surfaces, returning a
/// cleaned polyline intersection curve. Mirrors 0.3.2
/// `IntersectionCurveWithParameters::try_new`, keeping only the 3D leader.
fn build_intersection_curve<S>(
    surface0: S,
    surface1: S,
    poly: Polyline,
) -> Option<IntersectionCurve<Polyline, S, S>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let ic = IntersectionCurve::new(surface0.clone(), surface1.clone(), poly);
    let raw = ic.leader().clone();
    let len = raw.len();
    if len < 2 {
        return None;
    }
    let recover_triple = |point: Point3| {
        ic.search_nearest_point(point, None, None, 100).or_else(|| {
            let p0 = ic
                .surface0()
                .search_parameter(point, None, 100)
                .or_else(|| ic.surface0().search_nearest_parameter(point, None, 100))?;
            let p1 = ic
                .surface1()
                .search_parameter(point, None, 100)
                .or_else(|| ic.surface1().search_nearest_parameter(point, None, 100))?;
            let q0 = ic.surface0().evaluate(p0.0, p0.1);
            let q1 = ic.surface1().evaluate(p1.0, p1.1);
            Some((q0.midpoint(q1), p0.into(), p1.into()))
        })
    };
    // The original rescanning algorithm discarded the exact UV solution at
    // every sampled SSI vertex. On NURBS supports this makes the *next* point
    // redo two global 50x50 parameter preselections, before a 4-variable solve.
    // Continue the same intersection branch using the previously certified UVs
    // as Newton hints. Fall back to the existing global solver if local
    // continuation cannot prove exact incidence near the sampled leader.
    let mut prior = None::<((f64, f64), (f64, f64))>;
    let mut polyline = PolylineCurve(Vec::with_capacity(len));
    let support_tol =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_surface_projection_factor * KERNEL_EPSILON;
    let leader_tol = DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor;
    for index in 0..len {
        if index == len - 1 && raw[0].near(&raw[index]) {
            let first = polyline[0];
            polyline.push(first);
            break;
        }
        let parameter = index.to_f64()?;
        let point = raw[index];
        let tangent = ic.leader().derivative(parameter);
        let continuation = prior
            .filter(|_| tangent.magnitude2() > KERNEL_EPSILON * KERNEL_EPSILON)
            .and_then(|(uv0, uv1)| {
                ic.search_intersection_in_plane(point, tangent, Some(uv0), Some(uv1), 100)
            })
            .filter(|(candidate, uv0, uv1)| {
                // Newton convergence alone is not a correctness certificate.
                // Both carriers must agree on the SAME 3-D point and the
                // candidate must remain on the mesh-identified local branch.
                ic.surface0()
                    .evaluate(uv0.x, uv0.y)
                    .distance(ic.surface1().evaluate(uv1.x, uv1.y))
                    <= support_tol
                    && candidate.distance(point) <= leader_tol
            });
        let (q, uv0, uv1) = continuation
            .or_else(|| ic.search_triple(parameter, 100))
            .or_else(|| recover_triple(point))?;
        prior = Some(((uv0.x, uv0.y), (uv1.x, uv1.y)));
        polyline.push(q);
    }
    Some(IntersectionCurve::new(surface0, surface1, polyline))
}

type IntersectionTuple<S> = (Polyline, IntersectionCurve<Polyline, S, S>);

/// Marching intersection curves between two trimmed surfaces, keyed by their
/// tessellated polygon meshes. Returns, per intersection branch, the raw
/// interference polyline and the cleaned polyline intersection curve.
pub(super) fn intersection_curves<S>(
    surface0: S,
    polygon0: &PolygonMesh,
    surface1: S,
    polygon1: &PolygonMesh,
    tol: f64,
) -> Option<Vec<IntersectionTuple<S>>>
where
    S: ParametricSurface3D
        + Clone
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>,
{
    let interferences = polygon0.extract_interference(polygon1);
    let stitch_tolerance = 4.0 * tol.max(KERNEL_EPSILON);
    let polylines =
        remove_redundant_interference_strands(construct_polylines(&interferences), KERNEL_EPSILON);
    let polylines = stitch_nearby_polylines(polylines, stitch_tolerance)?;
    // Mesh intersections are seeds, not certified incidences. Project all
    // fragments onto both exact surfaces before attempting cross-gap joins.
    // Keep kernel-precision verification strict; only the mesh seed may move
    // up to the underlying discretization error.
    let projection_tol = DEFAULT_NUMERICAL_TOLERANCES.mesh_seed_floor;
    let support_tol =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_surface_projection_factor * KERNEL_EPSILON;
    let exact_fragments = polylines
        .into_iter()
        .map(|raw| {
            let recovered =
                build_intersection_curve(surface0.clone(), surface1.clone(), raw.clone())?;
            let corrected = recovered.leader().clone();
            if corrected.len() != raw.len() {
                return None;
            }
            let valid = raw.iter().zip(corrected.iter()).all(|(before, after)| {
                if before.distance(*after) > projection_tol {
                    return false;
                }
                let on = |surface: &S| {
                    let uv = surface
                        .search_parameter(*after, None, 100)
                        .or_else(|| surface.search_nearest_parameter(*after, None, 100));
                    uv.is_some_and(|(u, v)| surface.evaluate(u, v).distance(*after) <= support_tol)
                };
                on(&surface0) && on(&surface1)
            });
            valid.then_some(corrected)
        })
        .collect::<Option<Vec<_>>>()?;
    // This pass repairs numerical fragmentation introduced by tessellation.
    // Its admission envelope is an explicit kernel-precision bound rather than
    // a threshold derived from the user-facing modeling-error budget.
    let validation_tolerance =
        DEFAULT_NUMERICAL_TOLERANCES.boolean_fragment_validation_factor * KERNEL_EPSILON;
    let polylines = stitch_validated_linear_intersection_fragments(
        exact_fragments,
        &surface0,
        polygon0,
        &surface1,
        polygon1,
        validation_tolerance,
    );
    let polylines = stitch_validated_curved_intersection_fragments(
        polylines,
        &surface0,
        &surface1,
        validation_tolerance,
    );
    let polylines = stitch_overlapping_intersection_fragments(
        polylines,
        &surface0,
        &surface1,
        validation_tolerance,
    );
    polylines
        .into_iter()
        // Mesh interference can contain isolated point contacts. They do not
        // define a one-dimensional trim and therefore cannot divide either
        // face. Keep real intersection branches strict, but discard these
        // zero-dimensional artifacts instead of letting one make the whole
        // Boolean look like an SSI failure.
        .filter(|polyline| polyline.len() >= 2)
        .map(|polyline| {
            let curve =
                build_intersection_curve(surface0.clone(), surface1.clone(), polyline.clone())?;
            Some((polyline, curve))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use monstertruck_modeling::Surface;

    fn xy_mesh(x0: f64, x1: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(x0, -1.0, 0.0),
                    Point3::new(x1, -1.0, 0.0),
                    Point3::new(x1, 1.0, 0.0),
                    Point3::new(x0, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
        )
        .expect("fixture mesh is valid")
    }

    fn disconnected_xy_mesh(gap: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(-1.0, -1.0, 0.0),
                    Point3::new(0.0, -1.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(-1.0, 1.0, 0.0),
                    Point3::new(gap, -1.0, 0.0),
                    Point3::new(2.0, -1.0, 0.0),
                    Point3::new(2.0, 1.0, 0.0),
                    Point3::new(gap, 1.0, 0.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]]),
        )
        .expect("fixture mesh is valid")
    }

    fn zx_mesh(x0: f64, x1: f64) -> PolygonMesh {
        PolygonMesh::new(
            StandardAttributes {
                positions: vec![
                    Point3::new(x0, 0.0, -1.0),
                    Point3::new(x1, 0.0, -1.0),
                    Point3::new(x1, 0.0, 1.0),
                    Point3::new(x0, 0.0, 1.0),
                ],
                ..Default::default()
            },
            Faces::from_iter(vec![[0, 1, 2], [0, 2, 3]]),
        )
        .expect("fixture mesh is valid")
    }

    #[test]
    fn meshed_ssi_fragments_project_to_exact_surfaces_before_bridge_validation() {
        let surface0 = Plane::xy();
        let surface1 = Plane::zx();
        let raw = [
            PolylineCurve(vec![
                Point3::new(-0.5, 0.0, 0.0),
                Point3::new(-0.25, 2.0e-5, -1.5e-5),
                Point3::new(-0.0002, 2.0e-5, -1.5e-5),
            ]),
            PolylineCurve(vec![
                Point3::new(0.0002, 2.0e-5, -1.5e-5),
                Point3::new(0.25, 1.8e-5, -1.7e-5),
                Point3::new(0.5, 0.0, 0.0),
            ]),
        ];
        let projected = raw
            .iter()
            .map(|poly| {
                build_intersection_curve(surface0, surface1, poly.clone())
                    .expect("planes have an exact common line")
                    .leader()
                    .clone()
            })
            .collect::<Vec<_>>();
        for poly in &projected {
            assert!(
                poly.iter()
                    .all(|point| point.y.abs() < 4.0e-6 && point.z.abs() < 4.0e-6)
            );
        }
        let stitched = stitch_validated_linear_intersection_fragments(
            projected,
            &surface0,
            &xy_mesh(-1.0, 1.0),
            &surface1,
            &zx_mesh(-1.0, 1.0),
            4.0e-6,
        );
        assert_eq!(
            stitched.len(),
            1,
            "corrected samples permit exact line joining"
        );
    }

    #[test]
    fn curved_intersection_bridge_retains_exact_circular_arc_samples() {
        // A plane intersects this rational quadratic cylinder in a true
        // circular arc: its gap MUST be filled by curved, not chordal points.
        let half = std::f64::consts::FRAC_PI_4;
        let weight = half.cos();
        let controls = [
            Vector4::new(1.0, 0.0, -1.0, 1.0),
            Vector4::new(weight, weight, -weight, weight),
            Vector4::new(0.0, 1.0, -1.0, 1.0),
        ];
        let net = controls
            .into_iter()
            .map(|p| vec![p, Vector4::new(p.x, p.y, -p.z, p.w)])
            .collect::<Vec<_>>();
        let cylinder = Surface::NurbsSurface(NurbsSurface::new(BsplineSurface::new_unchecked(
            (KnotVector::bezier_knot(2), KnotVector::bezier_knot(1)),
            net,
        )));
        let plane = Surface::Plane(Plane::xy());
        for theta in [0.2_f64, 0.28, 0.288, 0.37] {
            let p = Point3::new(theta.cos(), theta.sin(), 0.0);
            let uv = cylinder.search_nearest_parameter(p, None, 100).unwrap();
            assert!(
                cylinder.evaluate(uv.0, uv.1).distance(p) < 1.0e-7,
                "test arc samples must lie on the exact cylinder"
            );
        }
        let on_circle = |theta: f64| Point3::new(theta.cos(), theta.sin(), 0.0);
        let a = PolylineCurve(vec![on_circle(0.20), on_circle(0.24), on_circle(0.28)]);
        let b = PolylineCurve(vec![on_circle(0.288), on_circle(0.33), on_circle(0.37)]);
        let corrected =
            corrected_curved_fragment_bridge(&a, false, &b, true, &plane, &cylinder, 4.0e-6)
                .expect("two exact intersecting surfaces admit continuous curved gap");
        assert!(corrected.len() >= 3);
        assert_eq!(corrected[0], a.back());
        assert_eq!(*corrected.last().unwrap(), b.front());
        let middle = corrected[corrected.len() / 2];
        assert!((middle.x * middle.x + middle.y * middle.y - 1.0).abs() < 1e-6);
        // Chord midpoint would be strictly inside the circle.
        let chord_mid = a.back().midpoint(b.front());
        assert!(middle.distance(chord_mid) > 1e-6);
        let stitched =
            stitch_validated_curved_intersection_fragments(vec![a, b], &plane, &cylinder, 4.0e-6);
        assert_eq!(stitched.len(), 1);
        assert!(stitched[0].len() > 6);
    }

    #[test]
    fn overlapping_exact_curved_fragments_splice_without_retracing() {
        let plane = Plane::xy();
        let cylinder = Plane::zx();
        let left = PolylineCurve(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.005, 0.0, 0.0),
            Point3::new(0.009, 0.0, 0.0),
            Point3::new(0.012, 0.0, 0.0),
        ]);
        let right = PolylineCurve(vec![
            Point3::new(0.008, 0.0, 0.0),
            Point3::new(0.011, 0.0, 0.0),
            Point3::new(0.016, 0.0, 0.0),
            Point3::new(0.025, 0.0, 0.0),
        ]);
        let result =
            stitch_overlapping_intersection_fragments(vec![left, right], &plane, &cylinder, 4.0e-6);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].front(), Point3::new(0.0, 0.0, 0.0));
        assert_eq!(result[0].back(), Point3::new(0.025, 0.0, 0.0));
        assert!(result[0].windows(2).all(|s| s[1].x > s[0].x));
    }

    #[test]
    fn overlapping_splice_refuses_transverse_crossings() {
        let plane = Plane::xy();
        let second = Plane::zx();
        let left = PolylineCurve(vec![
            Point3::new(-0.01, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.01, 0.0, 0.0),
        ]);
        let crossing = PolylineCurve(vec![
            Point3::new(0.0, -0.01, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.01, 0.0),
        ]);
        let result = stitch_overlapping_intersection_fragments(
            vec![left, crossing],
            &plane,
            &second,
            4.0e-6,
        );
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn curved_intersection_bridge_rejects_diverging_fragments() {
        let plane = Plane::xy();
        let transverse = Plane::zx();
        let a = PolylineCurve(vec![
            Point3::new(-0.05, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
        ]);
        // Endpoint touches the analytic intersection but the fragment
        // diverges along a different direction: no continuation allowed.
        let b = PolylineCurve(vec![
            Point3::new(0.003, 0.0, 0.0),
            Point3::new(0.004, 0.01, 0.0),
        ]);
        assert!(
            corrected_curved_fragment_bridge(&a, false, &b, true, &plane, &transverse, 4.0e-6)
                .is_none()
        );
    }

    #[test]
    fn intersection_leader_continuation_certifies_dense_transverse_planes() {
        let leader = PolylineCurve(
            (0..=100)
                .map(|i| {
                    let x = -0.5 + f64::from(i) / 100.0;
                    Point3::new(x, 0.0, 0.0)
                })
                .collect(),
        );
        let recovered = build_intersection_curve(Plane::xy(), Plane::zx(), leader)
            .expect("transverse planes share an exact line");
        assert_eq!(recovered.leader().len(), 101);
        for (i, point) in recovered.leader().iter().enumerate() {
            assert!((point.x - (-0.5 + i as f64 / 100.0)).abs() < 1.0e-6);
            assert!(point.y.abs() < 1.0e-6);
            assert!(point.z.abs() < 1.0e-6);
        }
    }

    #[test]
    fn validated_linear_stitch_skips_tiny_endpoint_segments() {
        let tolerance = 4.0e-6;
        let gap = 5.0e-5;
        let fragments = vec![
            PolylineCurve(vec![
                Point3::new(-0.5, 0.0, 0.0),
                Point3::new(-2.0e-6, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
            ]),
            PolylineCurve(vec![
                Point3::new(gap, 0.0, 0.0),
                Point3::new(gap + 2.0e-6, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ]),
        ];
        let stitched = stitch_validated_linear_intersection_fragments(
            fragments,
            &Plane::xy(),
            &xy_mesh(-1.0, 2.0),
            &Plane::zx(),
            &zx_mesh(-1.0, 2.0),
            tolerance,
        );
        assert_eq!(stitched.len(), 1);
        assert_eq!(stitched[0].front(), Point3::new(-0.5, 0.0, 0.0));
        assert_eq!(stitched[0].back(), Point3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn validated_linear_stitch_does_not_cross_real_trim_gap() {
        let tolerance = 4.0e-6;
        let gap = 5.0e-5;
        let fragments = vec![
            PolylineCurve(vec![
                Point3::new(-0.5, 0.0, 0.0),
                Point3::new(0.0, 0.0, 0.0),
            ]),
            PolylineCurve(vec![Point3::new(gap, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]),
        ];
        let stitched = stitch_validated_linear_intersection_fragments(
            fragments,
            &Plane::xy(),
            &disconnected_xy_mesh(gap),
            &Plane::zx(),
            &zx_mesh(-1.0, 2.0),
            tolerance,
        );
        assert_eq!(stitched.len(), 2);
    }
}
