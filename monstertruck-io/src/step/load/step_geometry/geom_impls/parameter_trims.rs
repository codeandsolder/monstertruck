//! Parameter-space trims: the 2D boundary a 3D STEP curve has on a STEP
//! surface, found exactly where the pair admits a closed form and by
//! projected sampling otherwise.

use super::{
    BoundaryCurveFromSamples, BoundedCurve, BsplineCurve, Conic2D, Conic3D, Curve2D, Curve3D,
    CurveTrimRef, ElementarySurface, Ellipse, EuclideanSpace, ExactParameterBoundary2D, Hyperbola,
    InnerSpace, Invertible, KnotVector, Line, Matrix3, Matrix4, MetricSpace, NurbsCurve, Ordering,
    Origin, Parabola, ParameterBoundary2D, ParameterCurve, ParameterDivision1D, ParametricCurve,
    ParametricCurve3D, ParametricSurface, Plane, Point2, Point3, Processor, RevolutionSurface,
    SearchNearestParameter, SearchParameter, SearchParameterHint2D, SphericalSurface, SquareMatrix,
    StepConvertingError, StepParameterCurve, StepRevolutionSurface, Surface, SurfaceCurve3D,
    SurfaceCurveTrimRef, SweepSurface, TOLERANCE, Tolerance, ToroidalSurface, Transform,
    TrimmedCurve, UnitCircle, UnitHyperbola, UnitParabola, Vector3, Vector4, uv_clamp,
};

#[cfg(test)]
use crate::step::load::step_geometry::{Sphere, SurfaceCurveKind, SurfaceCurveRepresentation};
#[cfg(test)]
use monstertruck_geometry::prelude::{BsplineSurface, Rad, Torus};
#[cfg(test)]
use monstertruck_traits::Transformed;
#[cfg(test)]
use std::f64::consts::TAU;

fn sampled_parameter_boundary<C>(
    curve: &C,
    surface: &Surface,
    tolerance: f64,
) -> Option<Vec<Point2>>
where
    C: ParametricCurve3D + BoundedCurve + ParameterDivision1D<Point = Point3>,
{
    fn abs_diff(previous: f64) -> impl Fn(&f64, &f64) -> Ordering {
        let distance = move |value: &f64| f64::abs(*value - previous);
        move |lhs: &f64, rhs: &f64| distance(lhs).total_cmp(&distance(rhs))
    }

    /// `surface` and `axis` are carried for the `MT_STEP_DEBUG_UV_CLAMP` lens
    /// only (`uv_clamp.rs`); with the variable unset every `record_axis` call is
    /// a `OnceLock` read and a return, and no value the lens computes is read
    /// back here.
    fn normalize_axis(
        value: f64,
        previous: Option<f64>,
        period: Option<f64>,
        range: Option<(f64, f64)>,
        surface: &Surface,
        axis: uv_clamp::Axis,
        real_range: bool,
    ) -> Option<f64> {
        if !value.is_finite() {
            None
        } else if let Some(previous) = previous {
            if let Some(period) = period {
                (-2..=2)
                    .map(|index| f64::mul_add(f64::from(index), period, value))
                    .min_by(abs_diff(previous))
            } else if let Some(range) = range {
                let clamped = clamp_near_range(value, range);
                uv_clamp::record_axis(
                    surface,
                    axis,
                    real_range,
                    uv_clamp::Event::Clamp {
                        moved: clamped - value,
                        to_max: clamped == range.1,
                    },
                );
                Some(clamped)
            } else {
                uv_clamp::record_axis(surface, axis, real_range, uv_clamp::Event::Unranged);
                Some(value)
            }
        } else if let Some((min, max)) = range {
            if let Some(period) = period {
                let span = max - min;
                if span.so_small() {
                    uv_clamp::record_axis(
                        surface,
                        axis,
                        real_range,
                        uv_clamp::Event::Periodic {
                            wrapped: min != value,
                            clamped: 0.0,
                            degenerate: true,
                        },
                    );
                    Some(min)
                } else {
                    let mut normalized = f64::floor((value - min) / period).mul_add(-period, value);
                    if normalized < min {
                        normalized += period;
                    }
                    if normalized > max {
                        normalized -= period;
                    }
                    let clamped = normalized.clamp(min, max);
                    uv_clamp::record_axis(
                        surface,
                        axis,
                        real_range,
                        uv_clamp::Event::Periodic {
                            wrapped: normalized != value,
                            clamped: clamped - normalized,
                            degenerate: false,
                        },
                    );
                    Some(clamped)
                }
            } else {
                let clamped = clamp_near_range(value, (min, max));
                uv_clamp::record_axis(
                    surface,
                    axis,
                    real_range,
                    uv_clamp::Event::Clamp {
                        moved: clamped - value,
                        to_max: clamped == max,
                    },
                );
                Some(clamped)
            }
        } else {
            uv_clamp::record_axis(surface, axis, real_range, uv_clamp::Event::Unranged);
            Some(value)
        }
    }

    fn clamp_near_range(value: f64, (min, max): (f64, f64)) -> f64 {
        if value < min && min - value < TOLERANCE {
            min
        } else if value > max && value - max < TOLERANCE {
            max
        } else {
            value
        }
    }

    // THE CLAMP IS NOT ASKED ABOUT AN AXIS THAT REPORTS A PLACEHOLDER.
    //
    // `clamp_near_range` snaps a value that overshoots the reported range by
    // less than `TOLERANCE` onto the boundary. On a knot-bounded axis that is
    // exactly right and it is doing real work: the net carries no data past the
    // knot vector, so an unbounded Newton that walked one ULP-cloud past the end
    // must come back. On a PLACEHOLDER axis there is no boundary to come back
    // to -- a plane's `[0, 1]^2` and a loaded cylinder's or cone's axial `[0, 1]`
    // are the values `parameter_range` returns so `range_tuple()` has something
    // to return, over an unbounded direction whose parameter is a world-scale
    // distance. Clamping there writes a number the solver did not produce, onto
    // a line in parameter space that means nothing (see
    // `uv_clamp::reported_range_bounds_the_surface`, and the 40-line doc on
    // `monstertruck-geometry/src/specifieds/plane.rs::parameter_range` that
    // forbids treating the square as a domain). Same repair as spec 011's
    // `a4604cef` on the modeling twin, and the same reason: do not ask a
    // placeholder a question it cannot answer.
    //
    // MEASURED before the change with `MT_STEP_DEBUG_UV_CLAMP=1` (the lens in
    // `uv_clamp.rs`, driven by `uv_clamp_probe.rs`), over two populations: the
    // 15 in-repo fixtures at full depth (58,244 chains, 939,416 points) and 7 of
    // the 8 big-assembly corpus files sampled by lowest solid id (81 solids,
    // 35,447 chains, 862,318 points; Scania-Engine not sampled -- its parse
    // alone exceeds nextest's 20-minute kill). **The clamp is LIVE**, and the
    // split is what decides the shape of the fix:
    //
    // | axis                     | calls   | MOVED  | -> min | -> max | max |move| |
    // |--------------------------|--------:|-------:|-------:|-------:|-----------:|
    // | IN-REPO                  |         |        |        |        |            |
    // | plane u (placeholder)    | 490,007 |  1,324 |  1,302 |     22 |    7.43e-7 |
    // | plane v (placeholder)    | 490,007 |  1,582 |  1,524 |     58 |   2.78e-10 |
    // | cylinder v (placeholder) |  67,472 |  3,655 |  3,569 |     86 |   5.63e-11 |
    // | cone v (placeholder)     |   3,926 |    704 |    704 |      0 |   9.61e-11 |
    // | bspline u+v (knots)      | 128,664 | 24,144 | 12,105 | 12,039 |    8.60e-7 |
    // | nurbs u+v (knots)        | 624,366 | 64,708 | 36,896 | 27,812 |    9.30e-7 |
    // | CORPUS                   |         |        |        |        |            |
    // | plane u (placeholder)    | 291,687 |    842 |    828 |     14 |    7.06e-7 |
    // | plane v (placeholder)    | 291,687 |    768 |    748 |     20 |    9.58e-7 |
    // | cylinder v (placeholder) | 311,100 | 17,492 | 14,717 |  2,775 |    2.50e-7 |
    // | cone v (placeholder)     | 103,061 | 26,338 | 23,057 |  3,281 |   1.21e-13 |
    // | bspline u+v (knots)      |  86,630 |  5,670 |  3,335 |  2,335 |    9.99e-7 |
    // | nurbs u+v (knots)        |  84,588 |  7,050 |  4,170 |  2,880 |    9.99e-7 |
    //
    // In-repo, 92.4% of all movement (88,852 of 96,117) is against a REAL knot
    // vector and is untouched by this change. **On the corpus the ratio
    // INVERTS**: 45,440 of 58,160 moves -- 78% -- are against a placeholder.
    // The in-repo fixtures are B-spline-heavy scan data; real assemblies are
    // analytic-heavy, so anyone measuring only in-repo would have concluded this
    // was a 7% edge case.
    //
    // The `-> min` column is why the placeholder half looked harmless: most
    // moves snap onto `min = 0`, which is the surface's own parameter ORIGIN --
    // a plane's placement point, a revolution profile's start -- so the write is
    // ~1e-10 and reads as de-noising. It is NOT principled de-noising: the same
    // face's sample at u = 5.0 gets none, and the noise floor is only small
    // because the placement happens to sit on the trim. The other **6,090
    // corpus moves snap onto `max = 1`**, one arbitrary unit along an unbounded
    // direction, which is the same write with none of the excuse.
    //
    // AND THE SQUARE IS NOT A DOMAIN, MEASURED AS SUCH: of the raw solver
    // answers on a placeholder axis, the fraction sitting OUTSIDE the reported
    // range by more than `TOLERANCE` -- where `clamp_near_range` does nothing at
    // all -- is 437,258 of 490,007 planar u in-repo (89.2%) and 274,196 of
    // 291,687 on the corpus (94.0%), reaching `u = 5.0e5` outside `[0, 1]`. A
    // range that four out of five samples violate by five orders of magnitude
    // is not a bound, and clamping the one-in-a-thousand that lands within
    // 1e-6 of it is an accident of where the placement was written.
    //
    // Only NON-PERIODIC axes are dropped. No placeholder axis carries a period
    // today (a revolution's turn and a sphere's/torus's angles are all real
    // bounds, and a `Line` profile has no period), so the guard is inert --
    // it is written this way so that the periodic wrap arm below, which is
    // selected by `range.is_some()`, provably cannot lose its range. Pinned by
    // `no_placeholder_axis_carries_a_period`.
    //
    // WHAT THIS CHANGE COSTS, measured by per-chain BIT digests of the produced
    // parameter loops, before and after, over both populations: 3,147 of 58,268
    // in-repo chains change (5.4%). NOT ONE changes its point count, and not one
    // flips between `Some` and `None` -- the chains that refused still refuse
    // and the chains that answered still answer, with the same number of
    // samples. The largest movement of any chain's u/v extent is 7.43e-7, and
    // only 19 chains move by more than 1e-9; the mode is ~1e-16, the ULP cloud
    // around the placement origin that used to be snapped flat to `0.0`.
    // Placeholder moves go 7,265 -> 0 and knot moves hold at exactly 88,852.
    let (u_real, v_real) = uv_clamp::reported_range_bounds_the_surface(surface);
    let normalize_uv = |uv: Point2, previous: Option<(f64, f64)>| {
        let (period_u, period_v) = (surface.period_u(), surface.period_v());
        let (urange, vrange) = surface.try_range_tuple();
        uv_clamp::record_point();
        uv_clamp::record_reported_excess(surface, uv_clamp::Axis::U, u_real, uv.x, urange);
        uv_clamp::record_reported_excess(surface, uv_clamp::Axis::V, v_real, uv.y, vrange);
        let urange = urange.filter(|_| u_real || period_u.is_some());
        let vrange = vrange.filter(|_| v_real || period_v.is_some());
        Some(Point2::new(
            normalize_axis(
                uv.x,
                previous.map(|(u, _)| u),
                period_u,
                urange,
                surface,
                uv_clamp::Axis::U,
                u_real,
            )?,
            normalize_axis(
                uv.y,
                previous.map(|(_, v)| v),
                period_v,
                vrange,
                surface,
                uv_clamp::Axis::V,
                v_real,
            )?,
        ))
    };
    let points = curve
        .try_parameter_division(curve.range_tuple(), tolerance)?
        .1;
    let project = |point: Point3, hint: Option<(f64, f64)>| {
        surface
            .search_parameter(point, hint, 100)
            .or_else(|| surface.search_parameter(point, None, 100))
            .or_else(|| surface.search_nearest_parameter(point, hint, 100))
            .or_else(|| surface.search_nearest_parameter(point, None, 100))
            .map(|(u, v)| Point2::new(u, v))
    };
    // ONE PASS, NO RETRY -- and the absence is measured, not an oversight.
    //
    // This chain used to end in `.or_else(|| <the identical hinted scan>)`.
    // `project` and `normalize_uv` are deterministic functions of
    // `(point, hint)`, neither solver keeps state across calls, and the retry
    // re-seeded the hint from the same `None` over the same `points` -- so it
    // recomputed the same sequence and returned the same `None`. Every time.
    //
    // The twin in `monstertruck-modeling/src/geometry/`
    // (`sampled_parameter_boundary`) retries differently: it re-projects every
    // point UNHINTED, abandoning the chain. That IS a different computation, and
    // copying it here was the other candidate. Both were measured before either
    // was chosen (spec 011 open item 6), over the 15 in-repo fixtures at full
    // depth plus 7 of the 8 corpus files sampled by solid:
    //
    // | population                       | chains | reach the retry | same-hinted rescues | UNHINTED rescues |
    // |----------------------------------|-------:|----------------:|--------------------:|-----------------:|
    // | 15 in-repo fixtures, all solids  | 58,244 |              14 |                   0 |                0 |
    // | ROTOR-201NAL-Z7, all 33 solids   | 11,810 |               2 |                   0 |                0 |
    // | Rocky_House, 12 of 156           |  8,938 |              11 |                   0 |                0 |
    // | Cruise_Assembly, 12              |    995 |               0 |                   0 |                0 |
    // | UMC-500, 12 of 217               |  3,094 |              52 |                   0 |                0 |
    // | Ai-14R, 4                        |  1,978 |               0 |                   0 |                0 |
    // | NissanGT-R, 4                    |  6,108 |              16 |                   0 |                0 |
    // | Scania-8x4, 4 of 832             |  3,256 |              32 |                   0 |                0 |
    // | **total**                        | 94,423 |         **127** |               **0** |            **0** |
    //
    // So the site was LIVE (127 chains, 0.13%, reached it) and the retry rescued
    // none of them -- as the determinism argument requires. The finding that
    // decided the shape of the fix is the last column: the modeling twin's
    // unhinted retry would have rescued **none of the same 127 either**. The
    // asymmetry between the twins was real, but "the other one is different, so
    // copy it" is not the conclusion -- it would have replaced a retry that
    // cannot rescue with one that measurably does not. Both are dead here; only
    // one of them looks alive.
    //
    // The 127 are not a loss the retry was hiding: the chain refuses typed, and
    // its caller (`monstertruck-solid`'s `reattach_preserved_face_trims`) falls
    // through to its own `sampled_trim_segment`. If a rescue is ever wanted at
    // THIS level, the measurement says it must be a genuinely different solve --
    // not a different hint into the same one.
    uv_clamp::begin_chain();
    let boundary = points
        .iter()
        .copied()
        .scan(None, |hint, point| {
            let uv = project(point, *hint).and_then(|uv| normalize_uv(uv, *hint));
            *hint = uv.map(|uv| (uv.x, uv.y));
            Some(uv)
        })
        .collect();
    uv_clamp::end_chain();
    boundary
}

fn exact_parameter_curve_on(curve: &Curve3D, surface: &Surface) -> Option<StepParameterCurve> {
    match (curve, surface) {
        (Curve3D::ParameterCurve(curve), surface)
            if SurfaceCurve3D::same_surface(curve.surface().as_ref(), surface) =>
        {
            Some(curve.clone())
        }
        (Curve3D::SurfaceCurve(curve), surface) => curve
            .parameter_curve_reparameterized_on(surface)
            .or_else(|| exact_parameter_curve_on(curve.leader(), surface)),
        (Curve3D::Line(curve), surface) => exact_line_parameter_curve_on(curve, surface),
        (
            Curve3D::BsplineCurve(curve),
            Surface::ElementarySurface(ElementarySurface::Plane(plane)),
        ) => exact_bspline_parameter_curve_on_plane(curve, plane, surface),
        (
            Curve3D::NurbsCurve(curve),
            Surface::ElementarySurface(ElementarySurface::Plane(plane)),
        ) => exact_nurbs_parameter_curve_on_plane(curve, plane, surface),
        (Curve3D::Conic(curve), surface) => exact_conic_parameter_curve_on(curve, surface),
        (Curve3D::IntersectionCurve(curve), surface) => {
            exact_parameter_curve_on(curve.leader().as_ref(), surface)
        }
        _ => None,
    }
}

#[inline]
fn matrix4_is_affine(transform: &Matrix4) -> bool {
    // This is a structural classification, not a geometric-nearness test. A
    // perspective term of 5e-7 is "near zero" at the kernel tolerance but is
    // still projective, and treating it as affine would turn an approximation
    // into proof-bearing trim geometry.
    const STRUCTURAL_TOLERANCE: f64 = 1.0e-12;
    transform.x.w.abs() <= STRUCTURAL_TOLERANCE
        && transform.y.w.abs() <= STRUCTURAL_TOLERANCE
        && transform.z.w.abs() <= STRUCTURAL_TOLERANCE
        && (transform.w.w - 1.0).abs() <= STRUCTURAL_TOLERANCE
}

fn transformed_xy_support_matches_plane(transform: &Matrix4, plane: &Plane) -> bool {
    if !matrix4_is_affine(transform) {
        return false;
    }

    let origin = transform.transform_point(Point3::origin());
    let x = transform.transform_point(Point3::new(1.0, 0.0, 0.0));
    let y = transform.transform_point(Point3::new(0.0, 1.0, 0.0));
    let u = x - origin;
    let v = y - origin;
    let basis_scale = u.magnitude().max(v.magnitude()).max(1.0);
    if u.cross(v).magnitude() <= TOLERANCE * basis_scale * basis_scale {
        return false;
    }

    let normal = u.cross(v);
    let scale = origin
        .to_vec()
        .magnitude()
        .max(plane.origin().to_vec().magnitude())
        .max(basis_scale)
        .max(1.0);
    let Some(parameter) = plane.parameter(origin) else {
        return false;
    };
    parameter.z.abs() <= TOLERANCE * scale && parallel_scaled(normal, plane.normal())
}

fn projected_conic_transform_on_plane(transform: &Matrix4, plane: &Plane) -> Option<Matrix3> {
    let project = |point| {
        plane
            .parameter(transform.transform_point(point))
            .map(|parameter| Point2::new(parameter.x, parameter.y))
    };
    let origin = project(Point3::origin())?;
    let u_axis = project(Point3::new(1.0, 0.0, 0.0))? - origin;
    let v_axis = project(Point3::new(0.0, 1.0, 0.0))? - origin;
    Some(Matrix3::from_cols(
        Vector3::new(u_axis.x, u_axis.y, 0.0),
        Vector3::new(v_axis.x, v_axis.y, 0.0),
        Vector3::new(origin.x, origin.y, 1.0),
    ))
}

// Numeric defense in depth after an analytic/exact construction. This must never
// be used to promote a sampled candidate to the exact API on its own.
fn pcurve_sanity_matches_surface_curve<C>(curve: &C, trim: &Curve2D, surface: &Surface) -> bool
where C: ParametricCurve3D<Point = Point3> + BoundedCurve {
    let (t0, t1) = curve.range_tuple();
    [
        t0,
        3.0f64.mul_add(t0, t1) * 0.25,
        f64::midpoint(t0, t1),
        3.0f64.mul_add(t1, t0) * 0.25,
        t1,
    ]
    .into_iter()
    .all(|parameter| {
        let uv = trim.evaluate(parameter);
        surface
            .evaluate(uv.x, uv.y)
            .near(&curve.evaluate(parameter))
    })
}

fn exact_ellipse_parameter_curve_on_plane(
    curve: &Ellipse<Point3, Matrix4>,
    plane: &Plane,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    transformed_xy_support_matches_plane(curve.transform(), plane).then_some(())?;
    let trimmed = TrimmedCurve::new(UnitCircle::new(), curve.entity().range());
    let mut projected = Processor::with_transform(
        trimmed,
        projected_conic_transform_on_plane(curve.transform(), plane)?,
    );
    if !curve.orientation() {
        projected.invert();
    }
    let trim = Curve2D::Conic(Conic2D::Ellipse(projected));
    pcurve_sanity_matches_surface_curve(curve, &trim, surface)
        .then(|| StepParameterCurve::new(Box::new(trim), Box::new(surface.clone())))
}

fn exact_hyperbola_parameter_curve_on_plane(
    curve: &Hyperbola<Point3, Matrix4>,
    plane: &Plane,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    transformed_xy_support_matches_plane(curve.transform(), plane).then_some(())?;
    let trimmed = TrimmedCurve::new(UnitHyperbola::new(), curve.entity().range());
    let mut projected = Processor::with_transform(
        trimmed,
        projected_conic_transform_on_plane(curve.transform(), plane)?,
    );
    if !curve.orientation() {
        projected.invert();
    }
    let trim = Curve2D::Conic(Conic2D::Hyperbola(projected));
    pcurve_sanity_matches_surface_curve(curve, &trim, surface)
        .then(|| StepParameterCurve::new(Box::new(trim), Box::new(surface.clone())))
}

fn exact_parabola_parameter_curve_on_plane(
    curve: &Parabola<Point3, Matrix4>,
    plane: &Plane,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    transformed_xy_support_matches_plane(curve.transform(), plane).then_some(())?;
    let trimmed = TrimmedCurve::new(UnitParabola::new(), curve.entity().range());
    let mut projected = Processor::with_transform(
        trimmed,
        projected_conic_transform_on_plane(curve.transform(), plane)?,
    );
    if !curve.orientation() {
        projected.invert();
    }
    let trim = Curve2D::Conic(Conic2D::Parabola(projected));
    pcurve_sanity_matches_surface_curve(curve, &trim, surface)
        .then(|| StepParameterCurve::new(Box::new(trim), Box::new(surface.clone())))
}

fn line_parameter_curve(line: Line<Point2>, surface: &Surface) -> StepParameterCurve {
    StepParameterCurve::new(Box::new(Curve2D::Line(line)), Box::new(surface.clone()))
}

fn affine_uv_parameter_curve(
    range: (f64, f64),
    front: Point2,
    back: Point2,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    let (t0, t1) = range;
    if !t0.is_finite() || !t1.is_finite() || t0.near(&t1) {
        return None;
    }
    let (start, end, points, reversed) = if t0 < t1 {
        (t0, t1, vec![front, back], false)
    } else {
        (t1, t0, vec![back, front], true)
    };
    let mut curve =
        BsplineCurve::new(KnotVector::from(vec![start, start, end, end]), points).ok()?;
    if reversed {
        curve.invert();
    }
    Some(StepParameterCurve::new(
        Box::new(Curve2D::BsplineCurve(curve)),
        Box::new(surface.clone()),
    ))
}

#[derive(Clone, Copy, Debug)]
struct AnalyticCircleFrame {
    center: Point3,
    radius: f64,
    normal: Vector3,
    front: Point3,
    tangent: Vector3,
    range: (f64, f64),
}

#[inline]
fn scaled_near(left: f64, right: f64, scale: f64) -> bool {
    (left - right).abs() <= TOLERANCE * scale.abs().max(left.abs()).max(right.abs()).max(1.0)
}

#[inline]
fn parallel_scaled(left: Vector3, right: Vector3) -> bool {
    let scale = left.magnitude() * right.magnitude();
    scale > TOLERANCE && left.cross(right).magnitude() <= TOLERANCE * scale.max(1.0)
}

fn ellipse_circle_frame_in_processor<S>(
    curve: &Ellipse<Point3, Matrix4>,
    surface: &Processor<S, Matrix4>,
) -> Option<AnalyticCircleFrame> {
    if !matrix4_is_affine(curve.transform()) || !matrix4_is_affine(surface.transform()) {
        return None;
    }
    let inverse = surface.transform().invert()?;
    let center = inverse.transform_point(curve.transform().transform_point(Point3::origin()));
    let axis0 = inverse.transform_vector(curve.transform().transform_vector(Vector3::unit_x()));
    let axis1 = inverse.transform_vector(curve.transform().transform_vector(Vector3::unit_y()));
    let radius0 = axis0.magnitude();
    let radius1 = axis1.magnitude();
    let scale = radius0.max(radius1).max(1.0);
    if !radius0.is_finite()
        || !radius1.is_finite()
        || radius0 <= TOLERANCE
        || radius1 <= TOLERANCE
        || !scaled_near(radius0, radius1, scale)
        || axis0.dot(axis1).abs() > TOLERANCE * radius0 * radius1
    {
        return None;
    }
    let cross = axis0.cross(axis1);
    let cross_magnitude = cross.magnitude();
    if cross_magnitude <= TOLERANCE * radius0 * radius1 {
        return None;
    }
    let (t0, t1) = curve.range_tuple();
    Some(AnalyticCircleFrame {
        center,
        radius: f64::midpoint(radius0, radius1),
        normal: cross / cross_magnitude,
        front: inverse.transform_point(curve.front()),
        tangent: inverse.transform_vector(curve.derivative(t0)),
        range: (t0, t1),
    })
}

fn circle_support_matches(
    frame: AnalyticCircleFrame,
    center: Point3,
    radius: f64,
    normal: Vector3,
) -> bool {
    let scale = frame
        .radius
        .abs()
        .max(radius.abs())
        .max(frame.center.to_vec().magnitude())
        .max(center.to_vec().magnitude())
        .max(1.0);
    frame.center.distance(center) <= TOLERANCE * scale
        && scaled_near(frame.radius, radius.abs(), scale)
        && parallel_scaled(frame.normal, normal)
}

fn unit_angular_rate(curve_tangent: Vector3, surface_tangent: Vector3) -> Option<f64> {
    let curve_speed = curve_tangent.magnitude();
    let surface_speed = surface_tangent.magnitude();
    let scale = curve_speed.max(surface_speed).max(1.0);
    if curve_speed <= TOLERANCE
        || surface_speed <= TOLERANCE
        || !scaled_near(curve_speed, surface_speed, scale)
        || curve_tangent.cross(surface_tangent).magnitude()
            > TOLERANCE * curve_speed * surface_speed
    {
        return None;
    }
    Some(if curve_tangent.dot(surface_tangent) >= 0.0 {
        1.0
    } else {
        -1.0
    })
}

#[inline]
const fn processor_parameter(orientation: bool, u: f64, v: f64) -> Point2 {
    if orientation {
        Point2::new(u, v)
    } else {
        Point2::new(v, u)
    }
}

fn exact_ellipse_parameter_curve_on_revolution<C>(
    curve: &Ellipse<Point3, Matrix4>,
    revolution: &Processor<RevolutionSurface<C>, Matrix4>,
    surface: &Surface,
) -> Option<StepParameterCurve>
where
    C: ParametricCurve3D<Point = Point3> + BoundedCurve,
{
    let frame = ellipse_circle_frame_in_processor(curve, revolution)?;
    let entity = revolution.entity();
    let (profile_parameter, angle) =
        entity.search_parameter(frame.front, SearchParameterHint2D::None, 30)?;
    let profile_point = entity.entity_curve().evaluate(profile_parameter);
    let axis = entity.axis();
    let offset = profile_point - entity.origin();
    let center = entity.origin() + axis * offset.dot(axis);
    let radius = (profile_point - center).magnitude();
    if !circle_support_matches(frame, center, radius, axis)
        || !entity.evaluate(profile_parameter, angle).near(&frame.front)
    {
        return None;
    }
    let rate = unit_angular_rate(frame.tangent, entity.derivative_v(profile_parameter, angle))?;
    let angle1 = f64::mul_add(rate, frame.range.1 - frame.range.0, angle);
    let front = processor_parameter(revolution.orientation(), profile_parameter, angle);
    let back = processor_parameter(revolution.orientation(), profile_parameter, angle1);
    let parameter_curve = affine_uv_parameter_curve(frame.range, front, back, surface)?;
    pcurve_sanity_matches_surface_curve(curve, parameter_curve.curve().as_ref(), surface)
        .then_some(parameter_curve)
}

fn exact_ellipse_parameter_curve_on_sphere(
    curve: &Ellipse<Point3, Matrix4>,
    processor: &SphericalSurface,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    let frame = ellipse_circle_frame_in_processor(curve, processor)?;
    let sphere = &processor.entity().0;
    let (u, v) = sphere.search_parameter(frame.front, SearchParameterHint2D::None, 30)?;
    let z = Vector3::unit_z();
    let parallel_center = sphere.center() + z * (sphere.radius() * u.cos());
    let parallel_radius = (sphere.radius() * u.sin()).abs();
    if parallel_radius > TOLERANCE
        && circle_support_matches(frame, parallel_center, parallel_radius, z)
    {
        let rate = unit_angular_rate(frame.tangent, sphere.derivative_v(u, v))?;
        let v1 = f64::mul_add(rate, frame.range.1 - frame.range.0, v);
        let front = processor_parameter(processor.orientation(), u, v);
        let back = processor_parameter(processor.orientation(), u, v1);
        let parameter_curve = affine_uv_parameter_curve(frame.range, front, back, surface)?;
        return pcurve_sanity_matches_surface_curve(
            curve,
            parameter_curve.curve().as_ref(),
            surface,
        )
        .then_some(parameter_curve);
    }

    let radial = Vector3::new(v.cos(), v.sin(), 0.0);
    let meridian_normal = radial.cross(z);
    if !circle_support_matches(frame, sphere.center(), sphere.radius(), meridian_normal) {
        return None;
    }
    let rate = unit_angular_rate(frame.tangent, sphere.derivative_u(u, v))?;
    let u1 = f64::mul_add(rate, frame.range.1 - frame.range.0, u);
    if u1 < -TOLERANCE || u1 > std::f64::consts::PI + TOLERANCE {
        return None;
    }
    let u1 = u1.clamp(0.0, std::f64::consts::PI);
    let front = processor_parameter(processor.orientation(), u, v);
    let back = processor_parameter(processor.orientation(), u1, v);
    let parameter_curve = affine_uv_parameter_curve(frame.range, front, back, surface)?;
    pcurve_sanity_matches_surface_curve(curve, parameter_curve.curve().as_ref(), surface)
        .then_some(parameter_curve)
}

fn exact_ellipse_parameter_curve_on_torus(
    curve: &Ellipse<Point3, Matrix4>,
    processor: &ToroidalSurface,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    let frame = ellipse_circle_frame_in_processor(curve, processor)?;
    let torus = processor.entity();
    let (u, v) = torus.search_parameter(frame.front, SearchParameterHint2D::None, 30)?;
    let z = Vector3::unit_z();

    let major_radius = torus.small_radius().mul_add(v.cos(), torus.large_radius());
    let major_center = torus.center() + z * (torus.small_radius() * v.sin());
    if major_radius.abs() > TOLERANCE
        && circle_support_matches(frame, major_center, major_radius.abs(), z)
    {
        let rate = unit_angular_rate(frame.tangent, torus.derivative_u(u, v))?;
        let u1 = f64::mul_add(rate, frame.range.1 - frame.range.0, u);
        let front = processor_parameter(processor.orientation(), u, v);
        let back = processor_parameter(processor.orientation(), u1, v);
        let parameter_curve = affine_uv_parameter_curve(frame.range, front, back, surface)?;
        return pcurve_sanity_matches_surface_curve(
            curve,
            parameter_curve.curve().as_ref(),
            surface,
        )
        .then_some(parameter_curve);
    }

    let radial = Vector3::new(u.cos(), u.sin(), 0.0);
    let minor_center = torus.center() + radial * torus.large_radius();
    let minor_normal = radial.cross(z);
    if !circle_support_matches(frame, minor_center, torus.small_radius(), minor_normal) {
        return None;
    }
    let rate = unit_angular_rate(frame.tangent, torus.derivative_v(u, v))?;
    let v1 = f64::mul_add(rate, frame.range.1 - frame.range.0, v);
    let front = processor_parameter(processor.orientation(), u, v);
    let back = processor_parameter(processor.orientation(), u, v1);
    let parameter_curve = affine_uv_parameter_curve(frame.range, front, back, surface)?;
    pcurve_sanity_matches_surface_curve(curve, parameter_curve.curve().as_ref(), surface)
        .then_some(parameter_curve)
}

fn exact_line_parameter_curve_on_linear_revolution(
    curve: &Line<Point3>,
    revolution: &Processor<RevolutionSurface<Line<Point3>>, Matrix4>,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    curve
        .exact_parameter_boundary_2d(revolution)
        .map(|boundary| wrap_exact_line_boundary(boundary, surface))
}

fn exact_line_parameter_curve_on_step_revolution(
    curve: &Line<Point3>,
    revolution: &StepRevolutionSurface,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    let Curve3D::Line(profile) = revolution.entity().entity_curve() else {
        return None;
    };
    let mut typed = Processor::with_transform(
        RevolutionSurface::by_revolution(
            *profile,
            revolution.entity().origin(),
            revolution.entity().axis(),
        ),
        *revolution.transform(),
    );
    if !revolution.orientation() {
        typed.invert();
    }
    exact_line_parameter_curve_on_linear_revolution(curve, &typed, surface)
}

fn exact_line_parameter_curve_on_plane(
    curve: &Line<Point3>,
    plane: &Plane,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    let (front_u, front_v) = plane.search_parameter(curve.front(), None, 1)?;
    let (back_u, back_v) = plane.search_parameter(curve.back(), None, 1)?;
    let trim = Curve2D::Line(Line(
        Point2::new(front_u, front_v),
        Point2::new(back_u, back_v),
    ));
    pcurve_sanity_matches_surface_curve(curve, &trim, surface)
        .then(|| StepParameterCurve::new(Box::new(trim), Box::new(surface.clone())))
}

fn wrap_exact_line_boundary<S>(
    boundary: ParameterCurve<Line<Point2>, S>,
    surface: &Surface,
) -> StepParameterCurve {
    let (line, _) = boundary.decompose();
    line_parameter_curve(line, surface)
}

fn exact_line_parameter_curve_on(
    curve: &Line<Point3>,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    match surface {
        Surface::ElementarySurface(ElementarySurface::Plane(plane)) => {
            exact_line_parameter_curve_on_plane(curve, plane, surface)
        }
        Surface::ElementarySurface(
            ElementarySurface::CylindricalSurface(revolution)
            | ElementarySurface::ConicalSurface(revolution),
        ) => exact_line_parameter_curve_on_linear_revolution(curve, revolution, surface),
        Surface::BsplineSurface(bspline) => curve
            .exact_parameter_boundary_2d(bspline)
            .map(|boundary| wrap_exact_line_boundary(boundary, surface)),
        Surface::NurbsSurface(nurbs) => curve
            .exact_parameter_boundary_2d(nurbs)
            .map(|boundary| wrap_exact_line_boundary(boundary, surface)),
        Surface::SweepSurface(SweepSurface::RevolutionSurface(revolution)) => {
            exact_line_parameter_curve_on_step_revolution(curve, revolution, surface)
        }
        Surface::ElementarySurface(
            ElementarySurface::Sphere(_) | ElementarySurface::ToroidalSurface(_),
        )
        | Surface::SweepSurface(SweepSurface::ExtrusionSurface(_)) => None,
    }
}

fn exact_bspline_parameter_curve_on_plane(
    curve: &BsplineCurve<Point3>,
    plane: &Plane,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    // A B-spline is an affine combination of its control points. If every
    // control point lies on the plane, projecting those controls into the
    // plane's (u, v) basis preserves the curve exactly, including its knot
    // vector and parameterization.
    let control_points = curve
        .control_points()
        .iter()
        .copied()
        .map(|point| {
            let uv = plane.parameter(point)?;
            let uv = Point2::new(uv.x, uv.y);
            plane.evaluate(uv.x, uv.y).near(&point).then_some(uv)
        })
        .collect::<Option<Vec<_>>>()?;
    let trim =
        Curve2D::BsplineCurve(BsplineCurve::new(curve.knot_vector().clone(), control_points).ok()?);
    pcurve_sanity_matches_surface_curve(curve, &trim, surface)
        .then(|| StepParameterCurve::new(Box::new(trim), Box::new(surface.clone())))
}

fn exact_nurbs_parameter_curve_on_plane(
    curve: &NurbsCurve<Vector4>,
    plane: &Plane,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    // Plane coordinates are affine, so they can be applied directly in
    // homogeneous space: (x*w, y*w, z*w, w) -> (u*w, v*w, w). Keeping the
    // original knots and weights therefore preserves the rational curve exactly.
    let control_points = curve
        .control_points()
        .iter()
        .copied()
        .map(|point| {
            if !point.w.is_finite() || point.w.abs() <= TOLERANCE {
                return None;
            }
            let world = Point3::new(point.x / point.w, point.y / point.w, point.z / point.w);
            let uv = plane.parameter(world)?;
            let uv = Point2::new(uv.x, uv.y);
            plane
                .evaluate(uv.x, uv.y)
                .near(&world)
                .then_some(Vector3::new(uv.x * point.w, uv.y * point.w, point.w))
        })
        .collect::<Option<Vec<_>>>()?;
    let trim = Curve2D::NurbsCurve(NurbsCurve::new(
        BsplineCurve::new(curve.knot_vector().clone(), control_points).ok()?,
    ));
    pcurve_sanity_matches_surface_curve(curve, &trim, surface)
        .then(|| StepParameterCurve::new(Box::new(trim), Box::new(surface.clone())))
}

fn exact_conic_parameter_curve_on(
    curve: &Conic3D,
    surface: &Surface,
) -> Option<StepParameterCurve> {
    match (curve, surface) {
        (Conic3D::Ellipse(curve), Surface::ElementarySurface(ElementarySurface::Plane(plane))) => {
            exact_ellipse_parameter_curve_on_plane(curve, plane, surface)
        }
        (
            Conic3D::Hyperbola(curve),
            Surface::ElementarySurface(ElementarySurface::Plane(plane)),
        ) => exact_hyperbola_parameter_curve_on_plane(curve, plane, surface),
        (Conic3D::Parabola(curve), Surface::ElementarySurface(ElementarySurface::Plane(plane))) => {
            exact_parabola_parameter_curve_on_plane(curve, plane, surface)
        }
        (
            Conic3D::Ellipse(curve),
            Surface::ElementarySurface(
                ElementarySurface::CylindricalSurface(revolution)
                | ElementarySurface::ConicalSurface(revolution),
            ),
        ) => exact_ellipse_parameter_curve_on_revolution(curve, revolution, surface),
        (
            Conic3D::Ellipse(curve),
            Surface::ElementarySurface(ElementarySurface::Sphere(sphere)),
        ) => exact_ellipse_parameter_curve_on_sphere(curve, sphere, surface),
        (
            Conic3D::Ellipse(curve),
            Surface::ElementarySurface(ElementarySurface::ToroidalSurface(torus)),
        ) => exact_ellipse_parameter_curve_on_torus(curve, torus, surface),
        (
            Conic3D::Ellipse(curve),
            Surface::SweepSurface(SweepSurface::RevolutionSurface(revolution)),
        ) => exact_ellipse_parameter_curve_on_revolution(curve, revolution, surface),
        _ => None,
    }
}

impl ParameterBoundary2D<Surface> for Curve3D {
    fn parameter_boundary_2d(&self, surface: &Surface, tolerance: f64) -> Option<Vec<Point2>> {
        match self {
            Self::ParameterCurve(curve) => {
                if curve.surface().as_ref() == surface {
                    curve
                        .curve()
                        .try_parameter_division(curve.curve().range_tuple(), tolerance)
                        .map(|(_, points)| points)
                } else {
                    sampled_parameter_boundary(curve, surface, tolerance)
                }
            }
            Self::SurfaceCurve(curve) => curve
                .parameter_curve_reparameterized_on(surface)
                .and_then(|parameter_curve| {
                    parameter_curve
                        .curve()
                        .try_parameter_division(parameter_curve.curve().range_tuple(), tolerance)
                        .map(|(_, points)| points)
                })
                .or_else(|| sampled_parameter_boundary(curve.leader(), surface, tolerance)),
            Self::IntersectionCurve(curve) => {
                exact_parameter_curve_on(curve.leader().as_ref(), surface)
                    .and_then(|parameter_curve| {
                        parameter_curve
                            .curve()
                            .try_parameter_division(
                                parameter_curve.curve().range_tuple(),
                                tolerance,
                            )
                            .map(|(_, points)| points)
                    })
                    .or_else(|| {
                        sampled_parameter_boundary(curve.leader().as_ref(), surface, tolerance)
                    })
                    .or_else(|| {
                        curve
                            .leader()
                            .try_parameter_division(curve.range_tuple(), tolerance)?
                            .0
                            .into_iter()
                            .map(|t| {
                                let (_, uv0, uv1) = curve.search_triple(t, 100)?;
                                if curve.surface0().as_ref() == surface {
                                    Some(uv0)
                                } else if curve.surface1().as_ref() == surface {
                                    Some(uv1)
                                } else {
                                    None
                                }
                            })
                            .collect::<Option<Vec<_>>>()
                    })
            }
            Self::Line(_)
            | Self::Polyline(_)
            | Self::Conic(_)
            | Self::BsplineCurve(_)
            | Self::NurbsCurve(_) => sampled_parameter_boundary(self, surface, tolerance),
        }
    }
}

impl ExactParameterBoundary2D<Surface> for Curve3D {
    type BoundaryCurve = StepParameterCurve;

    fn exact_parameter_boundary_2d(&self, surface: &Surface) -> Option<Self::BoundaryCurve> {
        CurveTrimRef::new(self, surface).try_into().ok()
    }
}

fn curve2d_from_sampled_boundary(points: Vec<Point2>) -> Option<Curve2D> {
    if points.len() < 2 {
        None
    } else {
        let front = points.first().copied()?;
        let back = points.last().copied()?;
        let line = Line(front, back);
        let is_linear = points.iter().copied().all(|point| {
            line.search_nearest_parameter(point, None, 1)
                .is_some_and(|t| line.evaluate(t).near(&point))
        });
        if is_linear {
            Some(Curve2D::Line(line))
        } else {
            let denom = (points.len() - 1) as f64;
            let knot_vec = KnotVector::from(
                std::iter::once(0.0)
                    .chain((0..points.len()).map(|index| index as f64 / denom))
                    .chain(std::iter::once(1.0))
                    .collect::<Vec<_>>(),
            );
            Some(Curve2D::BsplineCurve(BsplineCurve::new_unchecked(
                knot_vec, points,
            )))
        }
    }
}

impl BoundaryCurveFromSamples<Surface> for StepParameterCurve {
    fn boundary_curve_from_samples(surface: &Surface, points: Vec<Point2>) -> Option<Self> {
        curve2d_from_sampled_boundary(points)
            .map(|curve| Self::new(Box::new(curve), Box::new(surface.clone())))
    }
}

impl<'a> TryFrom<SurfaceCurveTrimRef<'a>> for StepParameterCurve {
    type Error = StepConvertingError;

    fn try_from(value: SurfaceCurveTrimRef<'a>) -> std::result::Result<Self, Self::Error> {
        let curve = value.curve();
        let surface = value.surface();
        exact_parameter_curve_on(&Curve3D::SurfaceCurve(curve.clone()), surface)
            .ok_or_else(|| "STEP surface curve has no exact trim on the requested surface.".into())
    }
}

impl<'a> TryFrom<CurveTrimRef<'a>> for StepParameterCurve {
    type Error = StepConvertingError;

    fn try_from(value: CurveTrimRef<'a>) -> std::result::Result<Self, Self::Error> {
        exact_parameter_curve_on(value.curve(), value.surface())
            .ok_or_else(|| "STEP curve has no exact trim on the requested surface.".into())
    }
}

#[test]
fn near_affine_projective_conic_is_not_promoted_to_exact() {
    let perspective = 0.5 * TOLERANCE;
    let transform = Matrix4::from_cols(
        Vector4::new(1.0, 0.0, 0.0, perspective),
        Vector4::new(0.0, 1.0, 0.0, 0.0),
        Vector4::new(0.0, 0.0, 1.0, 0.0),
        Vector4::new(0.0, 0.0, 0.0, 1.0),
    );
    let ellipse = Processor::with_transform(
        TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, TAU)),
        transform,
    );
    let curve = Curve3D::Conic(Conic3D::Ellipse(ellipse));
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn tilted_planar_conic_is_not_exact_even_when_old_five_samples_are_near_plane() {
    let slope = 0.9 * TOLERANCE;
    let transform = Matrix4::from_cols(
        Vector4::new(1.0, 0.0, slope, 0.0),
        Vector4::new(0.0, 1.0, slope, 0.0),
        Vector4::new(0.0, 0.0, 1.0, 0.0),
        Vector4::new(0.0, 0.0, 0.0, 1.0),
    );
    let ellipse = Processor::with_transform(
        TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, TAU)),
        transform,
    );
    let curve = Curve3D::Conic(Conic3D::Ellipse(ellipse));
    let plane = Plane::xy();
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(plane));

    for parameter in [0.0, 0.25 * TAU, 0.5 * TAU, 0.75 * TAU, TAU] {
        assert!(
            plane
                .parameter(curve.evaluate(parameter))
                .expect("test curve point must have plane coordinates")
                .z
                .abs()
                <= TOLERANCE
        );
    }
    assert!(
        plane
            .parameter(curve.evaluate(0.125 * TAU))
            .expect("test curve point must have plane coordinates")
            .z
            .abs()
            > TOLERANCE
    );
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn isoparametric_circle_on_cylinder_has_analytic_exact_trim() {
    let radius = 2.0;
    let axis = Vector3::unit_z();
    let point = Point3::new(radius, 0.0, 0.0);
    let surface =
        Surface::ElementarySurface(ElementarySurface::CylindricalSurface(Processor::new(
            RevolutionSurface::by_revolution(Line(point, point + axis), Point3::origin(), axis),
        )));
    let curve = Curve3D::Conic(Conic3D::Ellipse(
        Processor::new(TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, TAU)))
            .transformed(Matrix4::from_scale(radius)),
    ));

    let trim = curve
        .exact_parameter_boundary_2d(&surface)
        .expect("a revolution isoparametric circle has an analytic UV line proof");
    let Curve2D::BsplineCurve(line) = trim.curve().as_ref() else {
        panic!("isoparametric revolution circle must become a parameter-preserving UV line");
    };
    assert_eq!(line.degree(), 1);
    assert!(line.front().x.near(&line.back().x));
    assert!(line.range_tuple().0.near(&0.0));
    assert!(line.range_tuple().1.near(&TAU));
}

#[test]
fn tilted_great_circle_on_sphere_is_not_claimed_as_exact_uv_line() {
    let sphere = SphericalSurface::new(super::super::Sphere(
        monstertruck_geometry::prelude::Sphere::new(Point3::origin(), 1.0),
    ));
    let surface = Surface::ElementarySurface(ElementarySurface::Sphere(sphere));
    let curve = Curve3D::Conic(Conic3D::Ellipse(
        Processor::new(TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, TAU)))
            .transformed(Matrix4::from_angle_y(Rad(std::f64::consts::FRAC_PI_4))),
    ));

    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn exact_nurbs_trim_on_plane_preserves_weights_and_curve() {
    let plane = Plane::xy();
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(plane));
    let knots = KnotVector::bezier_knot(2);
    let source_control_points = vec![
        Vector4::new(0.0, 0.0, 0.0, 1.0),
        Vector4::new(1.0, 2.0, 0.0, 2.0),
        Vector4::new(3.0, 0.0, 0.0, 3.0),
    ];
    let curve = Curve3D::NurbsCurve(NurbsCurve::new(
        BsplineCurve::new(knots.clone(), source_control_points.clone())
            .expect("valid rational test curve"),
    ));

    let trim = curve
        .exact_parameter_boundary_2d(&surface)
        .expect("planar NURBS must project exactly into the plane chart");
    let Curve2D::NurbsCurve(trim_curve) = trim.curve().as_ref() else {
        panic!("planar NURBS must remain a rational trim");
    };
    assert_eq!(trim_curve.knot_vector(), &knots);
    assert_eq!(
        trim_curve.control_points().len(),
        source_control_points.len()
    );
    for (trim_point, source_point) in trim_curve
        .control_points()
        .iter()
        .zip(source_control_points.iter())
    {
        assert_eq!(trim_point.z, source_point.w);
    }

    let (t0, t1) = curve.range_tuple();
    for fraction in [0.0, 0.17, 0.5, 0.83, 1.0] {
        let t = t0 + (t1 - t0) * fraction;
        let uv = trim_curve.evaluate(t);
        assert!(surface.evaluate(uv.x, uv.y).near(&curve.evaluate(t)));
    }
}

#[test]
fn off_plane_nurbs_is_not_claimed_as_exact_planar_trim() {
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let curve = Curve3D::NurbsCurve(NurbsCurve::new(
        BsplineCurve::new(
            KnotVector::bezier_knot(2),
            vec![
                Vector4::new(0.0, 0.0, 0.0, 1.0),
                Vector4::new(1.0, 1.0, 2.0e-3, 2.0),
                Vector4::new(3.0, 0.0, 0.0, 3.0),
            ],
        )
        .expect("valid off-plane rational test curve"),
    ));
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn exact_bspline_trim_on_plane_preserves_curve() {
    let plane = Plane::new(
        Point3::new(3.0, -2.0, 5.0),
        Point3::new(5.0, -2.0, 5.0),
        Point3::new(3.0, 1.0, 5.0),
    );
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(plane));
    let knots = KnotVector::bezier_knot(2);
    let curve = Curve3D::BsplineCurve(
        BsplineCurve::new(
            knots.clone(),
            vec![
                Point3::new(3.5, -1.5, 5.0),
                Point3::new(4.25, 0.5, 5.0),
                Point3::new(4.75, -0.25, 5.0),
            ],
        )
        .expect("valid planar B-spline test curve"),
    );

    let trim = curve
        .exact_parameter_boundary_2d(&surface)
        .expect("planar B-spline must project exactly into the plane chart");
    let Curve2D::BsplineCurve(trim_curve) = trim.curve().as_ref() else {
        panic!("planar B-spline must remain a B-spline trim");
    };
    assert_eq!(trim_curve.knot_vector(), &knots);
    assert_eq!(trim_curve.control_points().len(), 3);

    let (t0, t1) = curve.range_tuple();
    for fraction in [0.0, 0.17, 0.5, 0.83, 1.0] {
        let t = t0 + (t1 - t0) * fraction;
        let uv = trim_curve.evaluate(t);
        assert!(surface.evaluate(uv.x, uv.y).near(&curve.evaluate(t)));
    }
}

#[test]
fn zero_range_bspline_refuses_exact_planar_retrim_without_panicking() {
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let curve = Curve3D::BsplineCurve(BsplineCurve::new_unchecked(
        KnotVector::from(vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
    ));
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn zero_range_nurbs_refuses_exact_planar_retrim_without_panicking() {
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let curve = Curve3D::NurbsCurve(NurbsCurve::new(BsplineCurve::new_unchecked(
        KnotVector::from(vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        vec![
            Vector4::new(0.0, 0.0, 0.0, 1.0),
            Vector4::new(0.5, 0.0, 0.0, 1.0),
            Vector4::new(1.0, 0.0, 0.0, 1.0),
        ],
    )));
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn off_plane_bspline_is_not_claimed_as_exact_planar_trim() {
    let plane = Plane::xy();
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(plane));
    let curve = Curve3D::BsplineCurve(
        BsplineCurve::new(
            KnotVector::bezier_knot(2),
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.5, 1.0e-3),
                Point3::new(1.0, 0.0, 0.0),
            ],
        )
        .expect("valid off-plane B-spline test curve"),
    );
    assert!(curve.exact_parameter_boundary_2d(&surface).is_none());
}

#[test]
fn surface_curve_leader_rebuilds_trim_on_equivalent_cylinder() {
    let center = Point3::origin();
    let axis = Vector3::unit_z();
    let radius = 2.0;
    let p = center + radius * Vector3::unit_x();
    let mut cylinder = Processor::new(RevolutionSurface::by_revolution(
        Line(p, p + axis),
        center,
        axis,
    ));
    cylinder.invert();
    let surface = Surface::ElementarySurface(ElementarySurface::CylindricalSurface(cylinder));
    let leader = Curve3D::Line(Line(
        Point3::new(radius, 0.0, -3.0),
        Point3::new(radius, 0.0, 4.0),
    ));
    let curve = Curve3D::SurfaceCurve(SurfaceCurve3D::new(
        SurfaceCurveKind::SurfaceCurve,
        Box::new(leader),
        Vec::new(),
        SurfaceCurveRepresentation::Curve3D,
    ));

    let trim = curve
        .exact_parameter_boundary_2d(&surface)
        .expect("3D leader must be reusable on a differently parameterized carrier");
    let (u0, u1) = trim.range_tuple();
    for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let u = u0 + (u1 - u0) * fraction;
        let uv = trim.curve().evaluate(u);
        let point = surface.evaluate(uv.x, uv.y);
        assert!(curve.search_parameter(point, None, 100).is_some());
    }
}

#[test]
fn exact_line_trim_is_rebuilt_on_cylindrical_surface() {
    let center = Point3::origin();
    let axis = Vector3::unit_z();
    let radius = 2.0;
    let p = center + radius * Vector3::unit_x();
    let mut cylinder = Processor::new(RevolutionSurface::by_revolution(
        Line(p, p + axis),
        center,
        axis,
    ));
    cylinder.invert();
    let surface = Surface::ElementarySurface(ElementarySurface::CylindricalSurface(cylinder));
    let curve = Curve3D::Line(Line(
        Point3::new(radius, 0.0, -3.0),
        Point3::new(radius, 0.0, 4.0),
    ));

    let trim = curve
        .exact_parameter_boundary_2d(&surface)
        .expect("a cylinder generatrix has an exact linear UV trim");
    let (u0, u1) = trim.range_tuple();
    for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let u = u0 + (u1 - u0) * fraction;
        let uv = trim.curve().evaluate(u);
        let point = surface.evaluate(uv.x, uv.y);
        assert!(
            curve.search_parameter(point, None, 100).is_some(),
            "trim point {point:?} is not on the source edge",
        );
    }
}

#[test]
fn sampled_parameter_boundary_preserves_unbounded_cylinder_axis_parameter() {
    let center = Point3::new(0.0, 0.0, 68.0);
    let axis = Vector3::unit_z();
    let radius = 0.3;
    let p = center + radius * Vector3::unit_x();
    let mut cylinder = Processor::new(RevolutionSurface::by_revolution(
        Line(p, p + axis),
        center,
        axis,
    ));
    cylinder.invert();
    let surface = Surface::ElementarySurface(ElementarySurface::CylindricalSurface(cylinder));
    let curve = Line(
        Point3::new(radius, 0.0, -6.25),
        Point3::new(radius, 0.0, 6.25),
    );

    let boundary = sampled_parameter_boundary(&curve, &surface, 0.001).unwrap();
    let max_abs = boundary
        .iter()
        .flat_map(|uv| [uv.x.abs(), uv.y.abs()])
        .fold(0.0, f64::max);

    assert!(max_abs > 10.0);
}

/// Spec 012 U2. A planar trim sample that overshoots the plane's fictional
/// `[0, 1]` square keeps the value the solver produced.
///
/// `Plane::xy()`'s `u` and `v` ARE world `x` and `y`, so this line's far
/// endpoint sits at `u = 1 + 5e-7` -- inside `TOLERANCE` of the reported `max`,
/// which is what `clamp_near_range` used to snap to exactly `1.0`. There is no
/// domain edge at `u = 1` on a plane; the number is a placeholder.
#[test]
fn a_planar_trim_sample_past_the_fictional_unit_square_is_not_snapped_to_it() {
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let overshoot = 5.0e-7;
    let curve = Line(
        Point3::new(1.0 - overshoot, 0.5, 0.0),
        Point3::new(1.0 + overshoot, 0.5, 0.0),
    );

    let boundary = sampled_parameter_boundary(&curve, &surface, 1.0e-3).unwrap();
    let last = *boundary.last().unwrap();

    assert!(
        last.x > 1.0,
        "the sample sits past u = 1 and must stay there; got {}",
        last.x,
    );
    assert!(
        (last.x - (1.0 + overshoot)).abs() < 1.0e-12,
        "the projected u must be the solver's answer, not the placeholder \
         boundary; got {}",
        last.x,
    );
}

/// The other half of U2: on a KNOT-bounded surface the clamp is load-bearing and
/// is kept. The net carries no data past the knot vector, so a Newton answer one
/// ULP-cloud past the end has a real boundary to come back to.
///
/// This bilinear patch spans world `[0, 1]^2` in the `z = 0` plane, so a sample
/// at world `x = 1 + 5e-7` extrapolates to `u = 1 + 5e-7` -- and must come back
/// to exactly `1.0`.
#[test]
fn a_knot_bounded_surface_still_snaps_a_sample_that_overshoots_its_knots() {
    let knots = KnotVector::bezier_knot(1);
    let surface = Surface::BsplineSurface(BsplineSurface::new_unchecked(
        (knots.clone(), knots),
        vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)],
            vec![Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ],
    ));
    let overshoot = 5.0e-7;
    let curve = Line(
        Point3::new(1.0 - overshoot, 0.5, 0.0),
        Point3::new(1.0 + overshoot, 0.5, 0.0),
    );

    let boundary = sampled_parameter_boundary(&curve, &surface, 1.0e-3).unwrap();

    assert_eq!(
        boundary.last().unwrap().x,
        1.0,
        "a knot vector IS a bound and the clamp must still snap onto it.",
    );
}

/// The guard on the drop is inert, and this is why: no axis whose reported range
/// is a placeholder carries a period. If that ever stops being true the drop
/// would take the range away from the periodic wrap arm, which is selected by
/// `range.is_some()` -- so the invariant is pinned rather than assumed.
#[test]
fn no_placeholder_axis_carries_a_period() {
    let center = Point3::new(0.0, 0.0, 0.0);
    let axis = Vector3::unit_z();
    let profile = Point3::new(0.3, 0.0, 0.0);
    let mut cylinder = Processor::new(RevolutionSurface::by_revolution(
        Line(profile, profile + axis),
        center,
        axis,
    ));
    cylinder.invert();
    let cone = Processor::new(RevolutionSurface::by_revolution(
        Line(profile, profile + axis + Vector3::unit_x()),
        center,
        axis,
    ));

    let surfaces = [
        Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy())),
        Surface::ElementarySurface(ElementarySurface::CylindricalSurface(cylinder)),
        Surface::ElementarySurface(ElementarySurface::ConicalSurface(cone)),
        Surface::ElementarySurface(ElementarySurface::Sphere(Processor::new(Sphere(
            monstertruck_geometry::prelude::Sphere::new(center, 2.0),
        )))),
        Surface::ElementarySurface(ElementarySurface::ToroidalSurface(Processor::new(
            Torus::new_unchecked(center, 3.0, 1.0),
        ))),
    ];

    for surface in surfaces {
        let (u_real, v_real) = uv_clamp::reported_range_bounds_the_surface(&surface);
        assert!(
            u_real || surface.period_u().is_none(),
            "placeholder u axis with a period: {surface:?}",
        );
        assert!(
            v_real || surface.period_v().is_none(),
            "placeholder v axis with a period: {surface:?}",
        );
    }
}

#[test]
fn raw_conic_boundary_without_pcurve_uses_sampled_projection_at_safe_tolerance() {
    let curve = Curve3D::Conic(Conic3D::Ellipse(
        Processor::new(TrimmedCurve::new(UnitCircle::<Point3>::new(), (0.0, TAU)))
            .transformed(Matrix4::from_nonuniform_scale(100.0, 100.0, 100.0)),
    ));
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));

    let boundary = curve
        .parameter_boundary_2d(&surface, 1.0e-3)
        .expect("safe raw conic projection should produce a parameter boundary.");

    assert!(boundary.len() > 4);
    assert!(
        boundary
            .iter()
            .any(|point| point.distance2(Point2::new(100.0, 0.0)) < 1.0e-6)
    );
}

#[test]
fn raw_line_boundary_without_pcurve_uses_sampled_projection() {
    let curve = Curve3D::Line(Line(
        Point3::new(0.25, 0.5, 0.0),
        Point3::new(0.75, 0.5, 0.0),
    ));
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));
    let boundary = curve
        .parameter_boundary_2d(&surface, 1.0e-5)
        .expect("safe raw line projection should produce a parameter boundary.");

    assert!(boundary.len() >= 2);
    assert!(
        boundary
            .first()
            .is_some_and(|point| point.near(&Point2::new(0.25, 0.5)))
    );
    assert!(
        boundary
            .last()
            .is_some_and(|point| point.near(&Point2::new(0.75, 0.5)))
    );
}

#[test]
fn surface_curve_line_without_pcurve_converts_to_exact_pcurve() {
    let leader = Curve3D::Line(Line(
        Point3::new(0.25, 0.5, 0.0),
        Point3::new(0.75, 0.5, 0.0),
    ));
    let curve = Curve3D::SurfaceCurve(SurfaceCurve3D::new(
        SurfaceCurveKind::SurfaceCurve,
        Box::new(leader),
        Vec::new(),
        SurfaceCurveRepresentation::Curve3D,
    ));
    let surface = Surface::ElementarySurface(ElementarySurface::Plane(Plane::xy()));

    let boundary = StepParameterCurve::try_from(CurveTrimRef::new(&curve, &surface))
        .expect("surface curve leader should produce an exact parameter boundary.");

    match boundary.curve().as_ref() {
        Curve2D::Line(line) => {
            assert!(line.front().near(&Point2::new(0.25, 0.5)));
            assert!(line.back().near(&Point2::new(0.75, 0.5)));
        }
        curve => panic!("expected line boundary, got {curve:?}"),
    }
}

#[test]
fn surface_curve_line_without_pcurve_on_cylinder_is_exact_generator() {
    let axis = Vector3::unit_z();
    let center = Point3::origin();
    let point = Point3::new(1.0, 0.0, 0.0);
    let profile = Line(point, point + axis);
    let surface = Surface::ElementarySurface(ElementarySurface::CylindricalSurface(
        Processor::new(RevolutionSurface::by_revolution(profile, center, axis)),
    ));
    let angle = 0.7;
    let leader = Curve3D::Line(Line(
        surface.evaluate(0.0, angle),
        surface.evaluate(1.0, angle),
    ));
    let curve = Curve3D::SurfaceCurve(SurfaceCurve3D::new(
        SurfaceCurveKind::SurfaceCurve,
        Box::new(leader),
        Vec::new(),
        SurfaceCurveRepresentation::Curve3D,
    ));

    let boundary = StepParameterCurve::try_from(CurveTrimRef::new(&curve, &surface))
        .expect("cylinder generator has an analytic UV-line proof");
    assert!(matches!(boundary.curve().as_ref(), Curve2D::Line(_)));
    assert!(curve.exact_parameter_boundary_2d(&surface).is_some());
}

#[test]
fn raw_line_without_pcurve_on_cylinder_is_exact_generator() {
    let axis = Vector3::unit_z();
    let center = Point3::origin();
    let point = Point3::new(1.0, 0.0, 0.0);
    let profile = Line(point, point + axis);
    let surface = Surface::ElementarySurface(ElementarySurface::CylindricalSurface(
        Processor::new(RevolutionSurface::by_revolution(profile, center, axis)),
    ));
    let angle = 0.7;
    let curve = Curve3D::Line(Line(
        surface.evaluate(0.0, angle),
        surface.evaluate(1.0, angle),
    ));

    let boundary = StepParameterCurve::try_from(CurveTrimRef::new(&curve, &surface))
        .expect("cylinder generator has an analytic UV-line proof");
    assert!(matches!(boundary.curve().as_ref(), Curve2D::Line(_)));
    assert!(curve.exact_parameter_boundary_2d(&surface).is_some());
}

#[test]
fn cone_apex_generator_is_exact_in_both_edge_directions() {
    let axis = Vector3::unit_z();
    let center = Point3::origin();
    let profile = Line(center, Point3::new(1.0, 0.0, 1.0));
    let surface = Surface::ElementarySurface(ElementarySurface::ConicalSurface(Processor::new(
        RevolutionSurface::by_revolution(profile, center, axis),
    )));
    let angle = std::f64::consts::FRAC_PI_4;
    let apex = surface.evaluate(0.0, angle);
    let rim = surface.evaluate(1.0, angle);

    for line in [Line(apex, rim), Line(rim, apex)] {
        let curve = Curve3D::Line(line);
        let boundary = StepParameterCurve::try_from(CurveTrimRef::new(&curve, &surface))
            .expect("cone generator through the apex has an analytic UV-line proof");
        let Curve2D::Line(uv) = boundary.curve().as_ref() else {
            panic!("cone generator must become a UV line");
        };
        assert!((uv.front().y - angle).abs() <= TOLERANCE);
        assert!((uv.back().y - angle).abs() <= TOLERANCE);
        assert!(curve.exact_parameter_boundary_2d(&surface).is_some());
    }
}
