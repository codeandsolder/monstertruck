//! Convert temporary Boolean intersection curves back to target curve types.

use super::tolerance::BooleanTolerance;
use super::{
    AltCurveShell, Alternative, BoundaryCurve2D, BoundedCurve, BoundedSurface, BsplineCurve, Cut,
    IntersectionCurve, Invertible, KnotVector, MetricSpace, NurbsCurve, ParameterCurve,
    ParametricCurve, Point3, PolylineCurve, ShapeOpsCurve, ShapeOpsSurface, Shell, SurfaceCurve,
    Vector4,
};
use monstertruck_core::tolerance::DEFAULT_NUMERICAL_TOLERANCES;
use num_traits::ToPrimitive;

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

/// A trimmed exact iso-section uses native spline parameters, not world-space distances.
fn trim_homogeneous_curve(
    mut curve: NurbsCurve<Vector4>,
    start: f64,
    end: f64,
    tolerances: BooleanTolerance,
) -> Option<NurbsCurve<Vector4>> {
    let (range0, range1) = curve.range_tuple();
    let slack = tolerances.iso_parameter_slack((range0, range1))?;
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

fn homogeneous_iso_curve<S: ShapeOpsSurface>(
    surface: &S,
    iso: IsoParameter,
    constant: f64,
    start: f64,
    end: f64,
    tolerances: BooleanTolerance,
) -> Option<NurbsCurve<Vector4>> {
    let mut homogeneous = surface.try_into_homogeneous_bspline_surface()?;
    let (u_range, v_range) = homogeneous.range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range, v_range);
    let curve = match iso {
        IsoParameter::U => {
            let slack = tolerances.iso_parameter_slack(u_range)?;
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
            let slack = tolerances.iso_parameter_slack(v_range)?;
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
    trim_homogeneous_curve(curve, start, end, tolerances)
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
            let Some(candidate) =
                homogeneous_iso_curve(surface, iso, constant, start, end, tolerances)
            else {
                continue;
            };
            if prove_iso_candidate(
                &candidate,
                intersection.surface0(),
                intersection.surface1(),
                leader,
                tolerances,
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
    tolerances: BooleanTolerance,
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
                if let Some(nurbs) = exact_iso_leader(ic, tolerances) {
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
        let tolerances = BooleanTolerance::new(1.0e-8);
        let trimmed = trim_homogeneous_curve(curve.clone(), 2.0e-10, 8.0e-10, tolerances)
            .expect("a tiny but well-resolved parameter interval is valid");
        let (lo, hi) = trimmed.range_tuple();
        assert!((lo - 2.0e-10).abs() <= 1.0e-15);
        assert!((hi - 8.0e-10).abs() <= 1.0e-15);
        assert!(trim_homogeneous_curve(curve.clone(), -1.0e-8, 8.0e-10, tolerances).is_none());
        assert!(trim_homogeneous_curve(curve, 2.0e-10, 2.0e-10 + 1.0e-16, tolerances).is_none());
    }
}
