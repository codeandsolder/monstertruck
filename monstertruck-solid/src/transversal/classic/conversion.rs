//! Convert temporary Boolean intersection curves back to target curve types.

use super::*;

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
    knots.transform(segment_count as f64, 0.0);
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

fn trim_homogeneous_curve(
    mut curve: NurbsCurve<Vector4>,
    start: f64,
    end: f64,
) -> Option<NurbsCurve<Vector4>> {
    if !start.is_finite() || !end.is_finite() || (start - end).abs() <= TOLERANCE {
        return None;
    }
    let reverse = start > end;
    let (lo, hi) = if reverse { (end, start) } else { (start, end) };
    let (range0, range1) = curve.range_tuple();
    if lo < range0 - TOLERANCE || hi > range1 + TOLERANCE {
        return None;
    }
    if lo > range0 + TOLERANCE {
        curve = curve.cut(lo);
    }
    let (_, current_end) = curve.range_tuple();
    if hi < current_end - TOLERANCE {
        let _tail = curve.cut(hi);
    }
    if reverse {
        curve.invert();
    }
    Some(curve)
}

fn homogeneous_iso_curve<S: ShapeOpsSurface>(
    surface: &S,
    iso: IsoParameter,
    constant: f64,
    start: f64,
    end: f64,
) -> Option<NurbsCurve<Vector4>> {
    let mut homogeneous = surface.try_into_homogeneous_bspline_surface()?;
    let (u_range, v_range) = homogeneous.range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range, v_range);
    let curve = match iso {
        IsoParameter::U => {
            if constant < u0 - TOLERANCE || constant > u1 + TOLERANCE {
                return None;
            }
            let row = if (constant - u0).abs() <= TOLERANCE {
                0
            } else if (constant - u1).abs() <= TOLERANCE {
                homogeneous.control_points().len().checked_sub(1)?
            } else {
                let _right = homogeneous.cut_u(constant);
                homogeneous.control_points().len().checked_sub(1)?
            };
            NurbsCurve::new(homogeneous.curve_v(row))
        }
        IsoParameter::V => {
            if constant < v0 - TOLERANCE || constant > v1 + TOLERANCE {
                return None;
            }
            let column = if (constant - v0).abs() <= TOLERANCE {
                0
            } else if (constant - v1).abs() <= TOLERANCE {
                homogeneous.control_points().first()?.len().checked_sub(1)?
            } else {
                let _right = homogeneous.cut_v(constant);
                homogeneous.control_points().first()?.len().checked_sub(1)?
            };
            NurbsCurve::new(homogeneous.curve_u(column))
        }
    };
    trim_homogeneous_curve(curve, start, end)
}

fn surface_parameter_with_hint<S: ShapeOpsSurface>(
    surface: &S,
    point: Point3,
    hint: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    let search = |search_hint| {
        surface
            .search_parameter(point, search_hint, 100)
            .or_else(|| surface.search_nearest_parameter(point, search_hint, 100))
    };
    let hinted = search(hint);
    if hint.is_some()
        && hinted.is_some_and(|(u, v)| surface.evaluate(u, v).distance(point) > 4.0 * TOLERANCE)
    {
        return search(None);
    }
    hinted.or_else(|| hint.and_then(|_| search(None)))
}

fn parameters_on_surface_along<S: ShapeOpsSurface>(
    surface: &S,
    points: impl IntoIterator<Item = Point3>,
) -> Option<Vec<(f64, f64)>> {
    let mut hint = None;
    points
        .into_iter()
        .map(|point| {
            let parameter = surface_parameter_with_hint(surface, point, hint)?;
            hint = Some(parameter);
            Some(parameter)
        })
        .collect()
}

fn point_on_surface_with_hint<S: ShapeOpsSurface>(
    surface: &S,
    point: Point3,
    hint: &mut Option<(f64, f64)>,
) -> bool {
    let Some(parameter) = surface_parameter_with_hint(surface, point, *hint) else {
        return false;
    };
    *hint = Some(parameter);
    surface.evaluate(parameter.0, parameter.1).distance(point) <= 4.0 * TOLERANCE
}

fn prove_iso_candidate<S: ShapeOpsSurface>(
    candidate: &NurbsCurve<Vector4>,
    surface0: &S,
    surface1: &S,
    leader: &PolylineCurve<Point3>,
) -> bool {
    let (t0, t1) = candidate.range_tuple();
    if candidate.subs(t0).distance(leader.front()) > 8.0 * TOLERANCE
        || candidate.subs(t1).distance(leader.back()) > 8.0 * TOLERANCE
    {
        return false;
    }
    let points = (0..=32)
        .map(|index| {
            let ratio = index as f64 / 32.0;
            candidate.subs(t0 * (1.0 - ratio) + t1 * ratio)
        })
        .collect::<Vec<_>>();
    let mut hint0 = None;
    if !points
        .iter()
        .copied()
        .all(|point| point_on_surface_with_hint(surface0, point, &mut hint0))
    {
        return false;
    }
    let mut hint1 = None;
    points
        .into_iter()
        .all(|point| point_on_surface_with_hint(surface1, point, &mut hint1))
}

fn leader_matches_iso_locus<S: ShapeOpsSurface>(
    surface: &S,
    params: &[(f64, f64)],
    leader: &PolylineCurve<Point3>,
    iso: IsoParameter,
    constant: f64,
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
                surface.evaluate(iso_u, iso_v).distance(point) <= 4.0 * TOLERANCE
            })
}

pub(super) fn exact_iso_leader<S: ShapeOpsSurface>(
    intersection: &IntersectionCurve<PolylineCurve<Point3>, S, S>,
) -> Option<NurbsCurve<Vector4>> {
    let leader = intersection.leader();
    if leader.len() < 2 {
        return None;
    }
    for surface in [intersection.surface0(), intersection.surface1()] {
        let Some(params) = parameters_on_surface_along(surface, leader.iter().copied()) else {
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
            let constant = values.iter().sum::<f64>() / values.len() as f64;
            // Inverse parameter search on tightly trimmed B-spline surfaces can
            // wander by hundreds of ulps even for an exact iso-parameter
            // intersection. Prove constancy in 3D instead: every marching
            // sample must lie on the exact iso locus within the Boolean validation
            // tolerance. The candidate is independently checked against both
            // carrier surfaces below before it can replace the marching leader.
            if !leader_matches_iso_locus(surface, &params, leader, iso, constant) {
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
            let Some(candidate) = homogeneous_iso_curve(surface, iso, constant, start, end) else {
                continue;
            };
            if prove_iso_candidate(
                &candidate,
                intersection.surface0(),
                intersection.surface1(),
                leader,
            ) {
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
/// serialization. Otherwise keep the leader exactly as a degree-1 B-spline.
/// `SurfaceCurve` still carries both supporting surfaces in either case.
pub(super) fn altshell_to_shell<C: ShapeOpsCurve<S>, S: ShapeOpsSurface>(
    altshell: &AltCurveShell<C, S>,
) -> Option<Shell<Point3, C, S>> {
    altshell.try_mapped(
        |p| Some(*p),
        |c| match c {
            Alternative::FirstType(c) => Some(c.clone()),
            Alternative::SecondType(ic) => {
                let surface0 = ic.surface0().clone();
                let surface1 = ic.surface1().clone();
                let boundary0: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                let boundary1: Option<ParameterCurve<BoundaryCurve2D, S>> = None;
                if let Some(nurbs) = exact_iso_leader(ic) {
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
    )
}
