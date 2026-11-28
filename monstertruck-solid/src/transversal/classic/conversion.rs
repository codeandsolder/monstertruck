//! Convert temporary Boolean intersection curves back to target curve types.

use super::tolerance::BooleanTolerance;
use super::{
    AltCurveShell, Alternative, BoundaryCurve2D, BoundedCurve, BoundedSurface, BsplineCurve, Cut,
    InnerProduct, IntersectionCurve, Invertible, KnotVector, MetricSpace, NurbsCurve,
    ParameterCurve, ParametricCurve, Point3, PolylineCurve, SearchNearestParameter, ShapeOpsCurve,
    ShapeOpsSurface, Shell, SurfaceCurve, Vector4,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use num_traits::ToPrimitive;
mod cylinder;

/// Convert a polyline leader to an exactly equivalent degree-1 B-spline.
///
/// `PolylineCurve` uses integer segment parameters, so scale the uniform knot
/// vector to the same `[0, segment_count]` domain. This preserves both the
/// guide geometry and parameterization exactly.
fn polyline_leader_bspline(polyline: &PolylineCurve<Point3>) -> Option<BsplineCurve<Point3>> {
    let segment_count = polyline.len().checked_sub(1)?;
    if segment_count == 0 {
        return None;
    }
    let mut knots = KnotVector::uniform_knot(1, segment_count);
    let segment_count = segment_count.to_f64()?;
    knots.transform(segment_count, 0.0);
    Some(BsplineCurve::new_unchecked(
        knots,
        polyline.iter().copied().collect(),
    ))
}

#[derive(Clone, Copy)]
enum IsoParameter {
    U,
    V,
}

/// A trimmed exact iso-section uses its *own* spline parameters, not world-space
/// distances or the analytic carrier's possibly different parameter frame.
fn trim_homogeneous_curve(
    mut curve: NurbsCurve<Vector4>,
    start: f64,
    end: f64,
) -> Option<NurbsCurve<Vector4>> {
    let (range0, range1) = curve.range_tuple();
    let slack = BooleanTolerance::iso_parameter_slack((range0, range1))?;
    if !start.is_finite() || !end.is_finite() || (start - end).abs() <= slack {
        return None;
    }
    let reverse = start > end;
    let (lo, hi) = if reverse { (end, start) } else { (start, end) };
    if lo < range0 - slack || hi > range1 + slack {
        return None;
    }
    // A tiny out-of-range drift does not authorize extrapolation of the curve.
    let lo = lo.max(range0);
    let hi = hi.min(range1);
    if hi - lo <= slack {
        return None;
    }
    if lo > range0 + slack {
        curve = curve.cut(lo);
    }
    let (_, current_end) = curve.range_tuple();
    if hi < current_end - slack {
        let _tail = curve.cut(hi);
    }
    if reverse {
        curve.invert();
    }
    Some(curve)
}

fn homogeneous_iso_support_curve<S: ShapeOpsSurface>(
    surface: &S,
    iso: IsoParameter,
    constant: f64,
) -> Option<NurbsCurve<Vector4>> {
    let mut homogeneous = surface.try_into_homogeneous_bspline_surface()?;
    let (u_range, v_range) = homogeneous.range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range, v_range);
    let curve = match iso {
        IsoParameter::U => {
            let slack = BooleanTolerance::iso_parameter_slack(u_range)?;
            if !constant.is_finite() || constant < u0 - slack || constant > u1 + slack {
                return None;
            }
            let coordinate = constant.clamp(u0, u1);
            let row = if (coordinate - u0).abs() <= slack {
                0
            } else if (coordinate - u1).abs() <= slack {
                homogeneous.control_points().len().checked_sub(1)?
            } else {
                let _right = homogeneous.cut_u(coordinate);
                homogeneous.control_points().len().checked_sub(1)?
            };
            NurbsCurve::new(homogeneous.curve_v(row))
        }
        IsoParameter::V => {
            let slack = BooleanTolerance::iso_parameter_slack(v_range)?;
            if !constant.is_finite() || constant < v0 - slack || constant > v1 + slack {
                return None;
            }
            let coordinate = constant.clamp(v0, v1);
            let column = if (coordinate - v0).abs() <= slack {
                0
            } else if (coordinate - v1).abs() <= slack {
                homogeneous.control_points().first()?.len().checked_sub(1)?
            } else {
                let _right = homogeneous.cut_v(coordinate);
                homogeneous.control_points().first()?.len().checked_sub(1)?
            };
            NurbsCurve::new(homogeneous.curve_u(column))
        }
    };
    Some(curve)
}

/// An analytic revolution may use a radian angular coordinate while its exact
/// rational B-spline conversion uses a normalized full-turn knot span. The
/// geometric inverse is valid only after the full circular carrier has been
/// algebraically certified, and the other carrier is a transverse exact plane.
/// This avoids treating approximate or arbitrarily reparameterized splines as
/// a proven exact intersection.
fn certified_plane_cylinder_iso<S: ShapeOpsSurface>(
    cylinder: &S,
    plane: &S,
    iso: IsoParameter,
) -> bool {
    let Some(mut net) = cylinder.try_into_homogeneous_bspline_surface() else {
        return false;
    };
    // The cylinder certificate expects rational circle in U and linear sweep
    // in V; a RevolutionSurface uses the opposite axis order.
    if matches!(iso, IsoParameter::U) {
        net.swap_axes();
    }
    let Some((_, axis, _, _)) = super::loops_store::certified_circular_cylinder(&net) else {
        return false;
    };
    let plane_net = plane.try_into_homogeneous_bspline_surface();
    let Some(exact_plane) = super::intersection_curve::certified_plane(plane, plane_net.as_ref())
    else {
        return false;
    };
    exact_plane.normal().cross(&axis).norm() <= DEFAULT_NUMERICAL_TOLERANCES.angular
}

/// Guide samples come from tessellation, not geometric identity. Use only to
/// choose the correct connected branch of an already certified exact circle.
/// Test the *whole* candidate against the ordered guide, so coincident circle
/// endpoints cannot promote the complement (or an entire extra revolution).
fn iso_candidate_follows_guide(
    candidate: &NurbsCurve<Vector4>,
    guide: &PolylineCurve<Point3>,
    tolerances: BooleanTolerance,
) -> bool {
    let (start, end) = candidate.range_tuple();
    let mut previous_progress = (0_usize, 0.0_f64);
    let sample_count = DEFAULT_NUMERICAL_TOLERANCES
        .boolean_iso_validation_segments
        .max(32);
    (0..=sample_count).all(|i| {
        let fraction = f64::from(i) / f64::from(sample_count);
        let point = candidate.evaluate(f64::mul_add(end - start, fraction, start));
        let Some((progress, residual)) = guide
            .windows(2)
            .enumerate()
            .filter_map(|(index, segment)| {
                let delta = segment[1] - segment[0];
                let squared = delta.norm_squared();
                if !squared.is_finite() || squared <= 0.0 {
                    return None;
                }
                let t = ((point - segment[0]).dot(&delta) / squared).clamp(0.0, 1.0);
                let distance = point.distance(segment[0] + delta * t);
                Some(((index, t), distance))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
        else {
            return false;
        };
        // Avoid lossy conversion of large polyline indices to f64. Exact
        // segment indices and local interpolation parameters form an ordered
        // progress key without collapsing distinct distant guide segments.
        let forward = progress.0 > previous_progress.0
            || (progress.0 == previous_progress.0 && progress.1 + 1.0e-8 >= previous_progress.1);
        previous_progress = progress;
        forward && residual.is_finite() && residual <= tolerances.mesh_seed
    })
}

fn surface_parameter_with_hint<S: ShapeOpsSurface>(
    surface: &S,
    point: Point3,
    hint: Option<(f64, f64)>,
    tolerances: BooleanTolerance,
) -> Option<(f64, f64)> {
    let search = |search_hint| {
        let trials = DEFAULT_NUMERICAL_TOLERANCES.boolean_inverse_search_trials;
        surface
            .search_parameter(point, search_hint, trials)
            .or_else(|| surface.search_nearest_parameter(point, search_hint, trials))
    };
    let hinted = search(hint);
    if hint.is_some()
        && hinted.is_some_and(|(u, v)| {
            surface.evaluate(u, v).distance(point) > tolerances.surface_support_residual()
        })
    {
        return search(None);
    }
    hinted.or_else(|| hint.and_then(|_| search(None)))
}

fn parameters_on_surface_along<S: ShapeOpsSurface>(
    surface: &S,
    points: impl IntoIterator<Item = Point3>,
    tolerances: BooleanTolerance,
) -> Option<Vec<(f64, f64)>> {
    let mut hint = None;
    points
        .into_iter()
        .map(|point| {
            let parameter = surface_parameter_with_hint(surface, point, hint, tolerances)?;
            hint = Some(parameter);
            Some(parameter)
        })
        .collect()
}

fn point_on_surface_with_hint<S: ShapeOpsSurface>(
    surface: &S,
    point: Point3,
    hint: &mut Option<(f64, f64)>,
    tolerances: BooleanTolerance,
) -> bool {
    let Some(parameter) = surface_parameter_with_hint(surface, point, *hint, tolerances) else {
        return false;
    };
    *hint = Some(parameter);
    surface.evaluate(parameter.0, parameter.1).distance(point)
        <= tolerances.surface_support_residual()
}

fn prove_iso_candidate<S: ShapeOpsSurface>(
    candidate: &NurbsCurve<Vector4>,
    surface0: &S,
    surface1: &S,
    leader: &PolylineCurve<Point3>,
    tolerances: BooleanTolerance,
) -> bool {
    let (t0, t1) = candidate.range_tuple();
    if candidate.evaluate(t0).distance(leader.front()) > tolerances.boundary_residual()
        || candidate.evaluate(t1).distance(leader.back()) > tolerances.boundary_residual()
    {
        return false;
    }
    let segments = DEFAULT_NUMERICAL_TOLERANCES.boolean_iso_validation_segments;
    let points = (0..=segments)
        .map(|index| {
            let ratio = f64::from(index) / f64::from(segments);
            candidate.evaluate(f64::mul_add(t1, ratio, t0 * (1.0 - ratio)))
        })
        .collect::<Vec<_>>();
    let mut hint0 = None;
    if !points
        .iter()
        .copied()
        .all(|point| point_on_surface_with_hint(surface0, point, &mut hint0, tolerances))
    {
        return false;
    }
    let mut hint1 = None;
    points
        .into_iter()
        .all(|point| point_on_surface_with_hint(surface1, point, &mut hint1, tolerances))
}

fn leader_matches_iso_locus<S: ShapeOpsSurface>(
    surface: &S,
    params: &[(f64, f64)],
    leader: &PolylineCurve<Point3>,
    iso: IsoParameter,
    constant: f64,
    tolerances: BooleanTolerance,
) -> bool {
    params.len() == leader.len()
        && constant.is_finite()
        && params
            .iter()
            .copied()
            .zip(leader.iter().copied())
            .all(|((u, v), point)| {
                let (iso_u, iso_v) = match iso {
                    IsoParameter::U => (constant, v),
                    IsoParameter::V => (u, constant),
                };
                surface.evaluate(iso_u, iso_v).distance(point)
                    <= tolerances.surface_support_residual()
            })
}

pub(super) fn exact_iso_leader<S: ShapeOpsSurface>(
    intersection: &IntersectionCurve<PolylineCurve<Point3>, S, S>,
    tolerances: BooleanTolerance,
) -> Option<NurbsCurve<Vector4>> {
    let leader = intersection.leader();
    if leader.len() < 2 {
        return None;
    }
    for surface in [intersection.surface0(), intersection.surface1()] {
        let Some(params) = parameters_on_surface_along(surface, leader.iter().copied(), tolerances)
        else {
            continue;
        };
        for iso in [IsoParameter::U, IsoParameter::V] {
            let values = params
                .iter()
                .map(|&(u, v)| match iso {
                    IsoParameter::U => u,
                    IsoParameter::V => v,
                })
                .collect::<Vec<_>>();
            let value_count = values.len().to_f64()?;
            let constant = values.iter().sum::<f64>() / value_count;
            // Inverse parameter search on tightly trimmed B-spline surfaces can
            // wander by hundreds of ulps even for an exact iso-parameter
            // intersection. Prove constancy in 3D instead: every marching
            // sample must lie on the exact iso locus within the Boolean validation
            // tolerance. The candidate is independently checked against both
            // carrier surfaces below before it can replace the marching leader.
            if !leader_matches_iso_locus(surface, &params, leader, iso, constant, tolerances) {
                continue;
            }
            let start = match iso {
                IsoParameter::U => params.first()?.1,
                IsoParameter::V => params.first()?.0,
            };
            let end = match iso {
                IsoParameter::U => params.last()?.1,
                IsoParameter::V => params.last()?.0,
            };
            let Some(full) = homogeneous_iso_support_curve(surface, iso, constant) else {
                continue;
            };
            // First try the inexpensive same-parameter-frame path. It remains
            // valid for ordinary native NURBS and planar surfaces.
            if let Some(candidate) = trim_homogeneous_curve(full.clone(), start, end)
                && prove_iso_candidate(
                    &candidate,
                    intersection.surface0(),
                    intersection.surface1(),
                    leader,
                    tolerances,
                )
                && iso_candidate_follows_guide(&candidate, leader, tolerances)
            {
                return Some(candidate);
            }
            // Radian and rational-circle knot coordinates are not the same
            // parameter. Find the exact circle's endpoint parameters by
            // geometric inverse only for a structurally certified cylinder
            // cut by an exact axis-normal plane.
            let other = if std::ptr::eq(surface, intersection.surface0()) {
                intersection.surface1()
            } else {
                intersection.surface0()
            };
            if !certified_plane_cylinder_iso(surface, other, iso) {
                continue;
            }
            let Some((start, end)) = full
                .search_nearest_parameter(leader.front(), None, 100)
                .zip(full.search_nearest_parameter(leader.back(), None, 100))
            else {
                continue;
            };
            let Some(candidate) = trim_homogeneous_curve(full, start, end) else {
                continue;
            };
            let carrier_ok = prove_iso_candidate(
                &candidate,
                intersection.surface0(),
                intersection.surface1(),
                leader,
                tolerances,
            );
            let guide_ok = iso_candidate_follows_guide(&candidate, leader, tolerances);
            if carrier_ok && guide_ok {
                return Some(candidate);
            }
        }
    }
    None
}

/// Convert an Alternative-curve shell back to the target curve type.
///
/// Surface intersection curves already carry a cleaned marching leader. When
/// that leader is provably an iso-parameter curve of either exact carrier,
/// replace the mesh-derived polyline by the exact homogeneous section before
/// serialization. Regular crossing cylinders instead receive an independently
/// bounded cubic leader. Otherwise keep the leader as a degree-1 B-spline.
/// `SurfaceCurve` still carries both supporting surfaces in either case.
pub(super) fn altshell_to_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    altshell: &AltCurveShell<C, S>,
    tolerances: BooleanTolerance,
) -> Option<Shell<Point3, C, S>> {
    let converted = altshell.try_mapped(
        |p| Some(*p),
        |c| match c {
            Alternative::FirstType(c) => Some(c.clone()),
            Alternative::SecondType(ic) => {
                let surface0 = ic.surface0().clone();
                let surface1 = ic.surface1().clone();
                let boundary0: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                let boundary1: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                if let Some(nurbs) = cylinder::crossing_cylinder_leader(ic, tolerances)
                    .or_else(|| exact_iso_leader(ic, tolerances))
                {
                    Some(
                        SurfaceCurve::with_boundaries(
                            surface0, surface1, nurbs, boundary0, boundary1,
                        )
                        .into(),
                    )
                } else {
                    let bspline = polyline_leader_bspline(ic.leader())?;
                    Some(
                        SurfaceCurve::with_boundaries(
                            surface0, surface1, bspline, boundary0, boundary1,
                        )
                        .into(),
                    )
                }
            }
        },
        |s| Some(s.clone()),
    )?;
    for (source_face, converted_face) in altshell.iter().zip(converted.iter()) {
        for (source_wire, converted_wire) in source_face
            .absolute_boundaries()
            .iter()
            .zip(converted_face.absolute_boundaries())
        {
            for (source, target) in source_wire.iter().zip(converted_wire.iter()) {
                let Alternative::SecondType(intersection) = source.curve() else {
                    continue;
                };
                let start = target.absolute_front().point();
                let end = target.absolute_back().point();
                let on_carrier = |surface: &S, point: Point3| {
                    surface
                        .search_parameter(point, None, 100)
                        .or_else(|| surface.search_nearest_parameter(point, None, 100))
                        .is_some_and(|(u, v)| {
                            surface.evaluate(u, v).distance(point)
                                <= tolerances.surface_support_residual()
                        })
                };
                if ![start, end].into_iter().all(|point| {
                    on_carrier(intersection.surface0(), point)
                        && on_carrier(intersection.surface1(), point)
                }) {
                    continue;
                }
                let mut curve = target.curve();
                if curve.front().distance(start) > tolerances.fragment_residual()
                    || curve.back().distance(end) > tolerances.fragment_residual()
                {
                    continue;
                }
                curve.snap_endpoints(start, end);
                target.set_curve(curve);
            }
        }
    }
    Some(converted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiny_native_knot_domain_remains_trimmable() {
        let mut knots = KnotVector::uniform_knot(1, 1);
        knots.transform(1.0e-9, 0.0);
        let curve = NurbsCurve::new(BsplineCurve::new_unchecked(
            knots,
            vec![
                Vector4::new(0.0, 0.0, 0.0, 1.0),
                Vector4::new(1.0, 0.0, 0.0, 1.0),
            ],
        ));
        let trimmed = trim_homogeneous_curve(curve.clone(), 2.0e-10, 8.0e-10)
            .expect("a tiny but well-resolved parameter interval is valid");
        let (lo, hi) = trimmed.range_tuple();
        assert!((lo - 2.0e-10).abs() <= 1.0e-15);
        assert!((hi - 8.0e-10).abs() <= 1.0e-15);
        assert!(trim_homogeneous_curve(curve.clone(), -1.0e-8, 8.0e-10).is_none());
        assert!(trim_homogeneous_curve(curve, 2.0e-10, 2.0e-10 + 1.0e-16).is_none());
    }
}
