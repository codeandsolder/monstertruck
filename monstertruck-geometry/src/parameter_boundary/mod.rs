use crate::prelude::*;

/// Face-local 2D trim curve used to preserve better boundary structure than a
/// raw UV polyline.
#[derive(
    Clone,
    Debug,
    ParametricCurve,
    BoundedCurve,
    ParameterDivision1D,
    Cut,
    Invertible,
    SearchNearestParameterD1,
    SearchParameterD1,
)]
pub enum BoundaryCurve2D {
    /// Linear trim segment in UV.
    Line(Line<Point2>),
    /// Degree-1 B-spline trim through sampled UV points.
    BsplineCurve(BsplineCurve<Point2>),
    /// Rational trim curve in UV.
    NurbsCurve(NurbsCurve<Vector3>),
}

impl<S> ParameterBoundary2D<S> for Line<Point3> {}

impl<S> ParameterBoundary2D<S> for BsplineCurve<Point3> {}

impl<S> ParameterBoundary2D<S> for NurbsCurve<Vector4> {}

impl<C, S> ExactParameterBoundary2D<S> for ParameterCurve<C, S>
where
    C: ParametricCurve2D<Point = Point2>
        + BoundedCurve<Point = Point2>
        + SearchParameter<CurveParameter, Point = Point2>
        + SearchNearestParameter<CurveParameter, Point = Point2>
        + Cut
        + Invertible
        + Clone,
    S: ParametricSurface3D
        + SearchParameter<SurfaceParameter, Point = Point3>
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + Clone
        + PartialEq,
{
    type BoundaryCurve = Self;

    fn exact_parameter_boundary_2d(&self, surface: &S) -> Option<Self::BoundaryCurve> {
        (self.surface() == surface).then(|| self.clone())
    }
}

impl<S: Clone> BoundaryCurveFromSamples<S> for ParameterCurve<Line<Point2>, S> {
    fn boundary_curve_from_samples(surface: &S, points: Vec<Point2>) -> Option<Self> {
        let (front, rest) = points.split_first()?;
        let back = rest.last()?;
        Some(Self::new(Line(*front, *back), surface.clone()))
    }
}

impl<S: Clone> BoundaryCurveFromSamples<S> for ParameterCurve<BoundaryCurve2D, S> {
    fn boundary_curve_from_samples(surface: &S, points: Vec<Point2>) -> Option<Self> {
        let (front, rest) = points.split_first()?;
        let back = rest.last()?;
        {
            let line = Line(*front, *back);
            let is_linear = points.iter().copied().all(|point| {
                line.search_nearest_parameter(point, None, 1)
                    .is_some_and(|t| line.evaluate(t).near(&point))
            });
            Some(Self::new(
                if is_linear {
                    BoundaryCurve2D::Line(line)
                } else {
                    let denom = crate::usize_to_f64(points.len() - 1);
                    let knot_vec = KnotVector::from(
                        std::iter::once(0.0)
                            .chain(
                                (0..points.len()).map(|index| crate::usize_to_f64(index) / denom),
                            )
                            .chain(std::iter::once(1.0))
                            .collect::<Vec<_>>(),
                    );
                    BoundaryCurve2D::BsplineCurve(BsplineCurve::new_unchecked(knot_vec, points))
                },
                surface.clone(),
            ))
        }
    }
}

impl TryFrom<ParameterCurve<Line<Point2>, Plane>> for BsplineCurve<Point3> {
    type Error = ();

    fn try_from(
        value: ParameterCurve<Line<Point2>, Plane>,
    ) -> std::result::Result<Self, Self::Error> {
        let (line, plane) = value.decompose();
        Ok(Self::from(Line(
            plane.evaluate(line.0.x, line.0.y),
            plane.evaluate(line.1.x, line.1.y),
        )))
    }
}

impl<C, S> ParameterBoundary2D<S> for ParameterCurve<C, S>
where
    C: ParametricCurve2D + BoundedCurve + ParameterDivision1D<Point = Point2>,
    S: PartialEq,
{
    fn parameter_boundary_2d(&self, surface: &S, tolerance: f64) -> Option<Vec<Point2>> {
        if self.surface() == surface {
            self.curve()
                .try_parameter_division(self.curve().range_tuple(), tolerance)
                .map(|(_, points)| points)
        } else {
            None
        }
    }
}

impl<C, S> ParameterBoundary2D<S> for TrimmedCurve<C>
where C: ParameterBoundary2D<S>
{
    fn parameter_boundary_2d(&self, surface: &S, tolerance: f64) -> Option<Vec<Point2>> {
        self.curve().parameter_boundary_2d(surface, tolerance)
    }
}

impl<C, T, S> ParameterBoundary2D<S> for Processor<C, T> where Self: ParametricCurve3D {}

fn exact_line_boundary_on_plane(
    line: &Line<Point3>,
    plane: &Plane,
) -> Option<ParameterCurve<Line<Point2>, Plane>> {
    let (u0, v0) = plane.search_parameter(line.front(), None, 1)?;
    let (u1, v1) = plane.search_parameter(line.back(), None, 1)?;
    let boundary = ParameterCurve::new(Line(Point2::new(u0, v0), Point2::new(u1, v1)), *plane);
    boundary
        .evaluate(0.5)
        .near(&line.evaluate(0.5))
        .then_some(boundary)
}

fn exact_line_boundary_on_affine_surface<S>(
    line: &Line<Point3>,
    surface: &S,
) -> Option<ParameterCurve<Line<Point2>, S>>
where
    S: Clone + ParametricSurface3D,
{
    let (u_range, v_range) = surface.try_range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range?, v_range?);
    let p00 = surface.evaluate(u0, v0);
    let p10 = surface.evaluate(u1, v0);
    let p01 = surface.evaluate(u0, v1);
    let p11 = surface.evaluate(u1, v1);
    let midpoint = surface.evaluate(f64::midpoint(u0, u1), f64::midpoint(v0, v1));
    (p11.near(&(p10 + (p01 - p00))) && midpoint.near(&(p00 + ((p10 - p00) + (p01 - p00)) * 0.5)))
        .then_some(())?;
    let u_axis = p10 - p00;
    let v_axis = p01 - p00;
    let uu = u_axis.dot(u_axis);
    let uv = u_axis.dot(v_axis);
    let vv = v_axis.dot(v_axis);
    let uv_squared = uv * uv;
    let denominator = f64::mul_add(uu, vv, -uv_squared);
    (!denominator.so_small()).then_some(())?;
    let project = |point: Point3| {
        let delta = point - p00;
        let du = delta.dot(u_axis);
        let dv = delta.dot(v_axis);
        let su = f64::mul_add(dv, -uv, du * vv) / denominator;
        let sv = f64::mul_add(du, -uv, dv * uu) / denominator;
        let parameter = Point2::new(f64::mul_add(u1 - u0, su, u0), f64::mul_add(v1 - v0, sv, v0));
        surface
            .evaluate(parameter.x, parameter.y)
            .near(&point)
            .then_some(parameter)
    };
    let front = project(line.front())?;
    let middle = project(line.evaluate(0.5))?;
    let back = project(line.back())?;
    let boundary = Line(front, back);
    boundary
        .evaluate(0.5)
        .near(&middle)
        .then(|| ParameterCurve::new(boundary, surface.clone()))
}

fn exact_line_boundary_on_homogeneous_extrusion_surface(
    line: &Line<Point3>,
    surface: &NurbsSurface<Vector4>,
) -> Option<ParameterCurve<Line<Point2>, NurbsSurface<Vector4>>> {
    let extrusion = match surface.try_into_analytic_surface_kind()? {
        AnalyticSurfaceKind::HomogeneousExtrusion(extrusion) => extrusion,
        _ => None?,
    };
    let vector = extrusion.vector;
    let line_vector = line.back() - line.front();
    let vector2 = vector.magnitude2();
    let line_vector2 = line_vector.magnitude2();
    (!vector2.so_small() && !line_vector2.so_small()).then_some(())?;
    (line_vector.cross(vector).magnitude2()
        <= KERNEL_EPSILON * KERNEL_EPSILON * line_vector2 * vector2)
        .then_some(())?;
    let base_curve = NurbsCurve::new(extrusion.curve.clone());
    let curve_tolerance = KERNEL_EPSILON
        * (extrusion.curve_range.1 - extrusion.curve_range.0)
            .abs()
            .max(1.0);
    let curve_parameter = base_curve.search_nearest_parameter(line.front(), None, 30)?;
    let back_curve_parameter = base_curve
        .search_nearest_parameter(line.back(), Some(curve_parameter), 30)
        .or_else(|| base_curve.search_nearest_parameter(line.back(), None, 30))?;
    ((back_curve_parameter - curve_parameter).abs() <= curve_tolerance).then_some(())?;
    let base_point = base_curve.evaluate(curve_parameter);
    let extrusion_parameter = |point: Point3| {
        let factor = (point - base_point).dot(vector) / vector2;
        let parameter = f64::mul_add(
            extrusion.extrusion_range.1 - extrusion.extrusion_range.0,
            factor,
            extrusion.extrusion_range.0,
        );
        let uv = match (extrusion.curve_axis, extrusion.extrusion_axis) {
            (SurfaceParameterAxis::U, SurfaceParameterAxis::V) => {
                Point2::new(curve_parameter, parameter)
            }
            (SurfaceParameterAxis::V, SurfaceParameterAxis::U) => {
                Point2::new(parameter, curve_parameter)
            }
            _ => None?,
        };
        surface.evaluate(uv.x, uv.y).near(&point).then_some(uv)
    };
    let front = extrusion_parameter(line.front())?;
    let middle = extrusion_parameter(line.evaluate(0.5))?;
    let back = extrusion_parameter(line.back())?;
    let boundary = Line(front, back);
    boundary
        .evaluate(0.5)
        .near(&middle)
        .then(|| ParameterCurve::new(boundary, surface.clone()))
}

fn restricted_bspline_segment<P>(
    boundary: &BsplineCurve<P>,
    start: f64,
    end: f64,
) -> Option<BsplineCurve<P>>
where
    P: Clone + control_point::ControlPoint<f64> + Tolerance,
{
    if !start.is_finite() || !end.is_finite() || start.near(&end) {
        return None;
    }
    let (lo, hi, reversed) = if start < end {
        (start, end, false)
    } else {
        (end, start, true)
    };
    let (b0, b1) = boundary.range_tuple();
    if lo < b0 - KERNEL_EPSILON || hi > b1 + KERNEL_EPSILON {
        return None;
    }
    let lo = lo.clamp(b0, b1);
    let hi = hi.clamp(b0, b1);
    let mut segment = if lo.near(&b0) {
        boundary.clone()
    } else {
        boundary.clone().cut(lo)
    };
    if !hi.near(&b1) {
        drop(segment.cut(hi));
    }
    if reversed {
        segment.invert();
    }
    Some(segment)
}

fn compatible_polynomial_control_nets(
    curve: &BsplineCurve<Point3>,
    boundary: &BsplineCurve<Point3>,
) -> bool {
    let mut curve = curve.clone();
    let mut boundary = boundary.clone();
    curve.knot_normalize();
    boundary.knot_normalize();
    curve.syncro_degree(&mut boundary);
    curve.syncro_knots(&mut boundary);
    curve.control_points().len() == boundary.control_points().len()
        && curve
            .control_points()
            .iter()
            .zip(boundary.control_points())
            .all(|(left, right)| {
                let scale = left
                    .to_vec()
                    .magnitude()
                    .max(right.to_vec().magnitude())
                    .max(1.0);
                left.distance(*right) <= KERNEL_EPSILON * scale
            })
}

fn compatible_homogeneous_control_nets(
    curve: &BsplineCurve<Vector4>,
    boundary: &BsplineCurve<Vector4>,
) -> bool {
    let mut curve = curve.clone();
    let mut boundary = boundary.clone();
    curve.knot_normalize();
    boundary.knot_normalize();
    curve.syncro_degree(&mut boundary);
    curve.syncro_knots(&mut boundary);
    if curve.control_points().len() != boundary.control_points().len() {
        return false;
    }

    // Rational homogeneous control nets are unique for a fixed basis up to one
    // global non-zero scale.  Recover that scale from the largest-magnitude
    // component, then compare the complete synchronized nets.
    let mut pivot = None::<(f64, f64)>;
    for (left, right) in curve.control_points().iter().zip(boundary.control_points()) {
        for (l, r) in [
            (left.x, right.x),
            (left.y, right.y),
            (left.z, right.z),
            (left.w, right.w),
        ] {
            if r.abs() > pivot.map_or(0.0, |(_, best)| best.abs()) {
                pivot = Some((l, r));
            }
        }
    }
    let Some((left_pivot, right_pivot)) = pivot else {
        return false;
    };
    if right_pivot.abs() <= KERNEL_EPSILON {
        return false;
    }
    let ratio = left_pivot / right_pivot;
    if !ratio.is_finite() || ratio.abs() <= KERNEL_EPSILON {
        return false;
    }
    curve
        .control_points()
        .iter()
        .zip(boundary.control_points())
        .all(|(left, right)| {
            let scaled = *right * ratio;
            let scale = left.magnitude().max(scaled.magnitude()).max(1.0);
            (*left - scaled).magnitude() <= KERNEL_EPSILON * scale
        })
}

fn exact_bspline_boundary_segment(
    curve: &BsplineCurve<Point3>,
    boundary: &BsplineCurve<Point3>,
) -> Option<(f64, f64)> {
    let start = boundary.search_parameter(curve.front(), None, 100)?;
    let end = boundary
        .search_parameter(curve.back(), Some(start), 100)
        .or_else(|| boundary.search_parameter(curve.back(), None, 100))?;
    let segment = restricted_bspline_segment(boundary, start, end)?;
    compatible_polynomial_control_nets(curve, &segment).then_some((start, end))
}

fn exact_homogeneous_boundary_segment(
    curve: &BsplineCurve<Vector4>,
    boundary: &BsplineCurve<Vector4>,
) -> Option<(f64, f64)> {
    let curve_nurbs = NurbsCurve::new(curve.clone());
    let boundary_nurbs = NurbsCurve::new(boundary.clone());
    let start = boundary_nurbs.search_parameter(curve_nurbs.front(), None, 100)?;
    let end = boundary_nurbs
        .search_parameter(curve_nurbs.back(), Some(start), 100)
        .or_else(|| boundary_nurbs.search_parameter(curve_nurbs.back(), None, 100))?;
    let segment = restricted_bspline_segment(boundary, start, end)?;
    compatible_homogeneous_control_nets(curve, &segment).then_some((start, end))
}

fn exact_nurbs_boundary_segment(
    curve: &NurbsCurve<Vector4>,
    boundary: &NurbsCurve<Vector4>,
) -> Option<(f64, f64)> {
    exact_homogeneous_boundary_segment(curve.non_rationalized(), boundary.non_rationalized())
}

fn exact_boundary_line_on_surface<B>(
    last_u: usize,
    last_v: usize,
    column_curve: impl Fn(usize) -> B,
    row_curve: impl Fn(usize) -> B,
    ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64)),
    match_segment: impl Fn(&B) -> Option<(f64, f64)>,
) -> Option<Line<Point2>> {
    match_segment(&column_curve(0))
        .map(|(s0, s1)| Line(Point2::new(u0, s0), Point2::new(u0, s1)))
        .or_else(|| {
            match_segment(&column_curve(last_u))
                .map(|(s0, s1)| Line(Point2::new(u1, s0), Point2::new(u1, s1)))
        })
        .or_else(|| {
            match_segment(&row_curve(0))
                .map(|(s0, s1)| Line(Point2::new(s0, v0), Point2::new(s1, v0)))
        })
        .or_else(|| {
            match_segment(&row_curve(last_v))
                .map(|(s0, s1)| Line(Point2::new(s0, v1), Point2::new(s1, v1)))
        })
}

fn exact_bspline_boundary_on_surface(
    curve: &BsplineCurve<Point3>,
    surface: &BsplineSurface<Point3>,
) -> Option<ParameterCurve<Line<Point2>, BsplineSurface<Point3>>> {
    let (u_range, v_range) = surface.try_range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range?, v_range?);
    let last_u = surface.control_points().len().checked_sub(1)?;
    let last_v = surface.control_points().first()?.len().checked_sub(1)?;
    let boundary = exact_boundary_line_on_surface(
        last_u,
        last_v,
        |index| surface.curve_v(index),
        |index| surface.curve_u(index),
        ((u0, u1), (v0, v1)),
        |boundary| exact_bspline_boundary_segment(curve, boundary),
    )?;
    Some(ParameterCurve::new(boundary, surface.clone()))
}

fn exact_nurbs_boundary_on_surface(
    curve: &NurbsCurve<Vector4>,
    surface: &NurbsSurface<Vector4>,
) -> Option<ParameterCurve<Line<Point2>, NurbsSurface<Vector4>>> {
    let (u_range, v_range) = surface.try_range_tuple();
    let ((u0, u1), (v0, v1)) = (u_range?, v_range?);
    let last_u = surface.control_points().len().checked_sub(1)?;
    let last_v = surface.control_points().first()?.len().checked_sub(1)?;
    let boundary = exact_boundary_line_on_surface(
        last_u,
        last_v,
        |index| surface.curve_v(index),
        |index| surface.curve_u(index),
        ((u0, u1), (v0, v1)),
        |boundary| exact_nurbs_boundary_segment(curve, boundary),
    )?;
    Some(ParameterCurve::new(boundary, surface.clone()))
}

fn nearest_periodic_angle(value: f64, reference: f64) -> f64 {
    ((reference - value) / (2.0 * std::f64::consts::PI))
        .round()
        .mul_add(2.0 * std::f64::consts::PI, value)
}

fn exact_revolution_isoparametric_boundary_from_homogeneous<C>(
    curve: &BsplineCurve<Vector4>,
    front: Point3,
    back: Point3,
    surface: &RevolutionSurface<C>,
) -> Option<Line<Point2>>
where
    C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve,
{
    // Search only proposes the candidate branch and endpoint profile
    // parameters.  Acceptance below is structural: the complete input curve
    // must equal the corresponding rotated profile segment after both exact
    // homogeneous representations are put on the same spline basis.
    let (u0, raw_v0) = surface.search_parameter(front, SearchParameterHint2D::None, 30)?;
    let (u1, raw_v1) =
        surface.search_parameter(back, SearchParameterHint2D::Parameter(u0, raw_v0), 30)?;
    let (profile_start, profile_end) = surface.entity_curve().range_tuple();
    let front_singular = (surface.is_front_fixed() && u0.near(&profile_start))
        || (surface.is_back_fixed() && u0.near(&profile_end));
    let back_singular = (surface.is_front_fixed() && u1.near(&profile_start))
        || (surface.is_back_fixed() && u1.near(&profile_end));

    // The angular coordinate at a point on the revolution axis is undefined.
    // Do not compare an arbitrary spelling returned by search against a real
    // angle from the other endpoint. Singular endpoints contribute only u;
    // a non-singular endpoint chooses the angular branch. If both endpoints
    // are singular, an interior point proposes v and the complete control-net
    // equality below still supplies the proof.
    let angle = match (front_singular, back_singular) {
        (false, false) => {
            let v1 = nearest_periodic_angle(raw_v1, raw_v0);
            if (v1 - raw_v0).abs() > KERNEL_EPSILON {
                return None;
            }
            raw_v0
        }
        (true, false) => raw_v1,
        (false, true) => raw_v0,
        (true, true) => {
            let nurbs = NurbsCurve::new(curve.clone());
            let (t0, t1) = nurbs.range_tuple();
            let midpoint = nurbs.evaluate(f64::midpoint(t0, t1));
            surface
                .search_parameter(midpoint, SearchParameterHint2D::None, 30)?
                .1
        }
    };

    let mut segment = if let Some(mut continued) = surface
        .entity_curve()
        .try_into_homogeneous_bspline_curve_over((u0, u1))
    {
        if u0 > u1 {
            continued.invert();
        }
        continued
    } else {
        let profile = surface
            .entity_curve()
            .try_into_homogeneous_bspline_curve()?;
        restricted_bspline_segment(&profile, u0, u1)?
    };
    let origin = surface.origin().to_vec();
    let rotation = Matrix3::from_axis_angle(surface.axis(), Rad(angle));
    segment.control_points_mut().for_each(|point| {
        let xyz = Vector3::new(point.x, point.y, point.z);
        let rotated = rotation * (xyz - origin * point.w) + origin * point.w;
        point.x = rotated.x;
        point.y = rotated.y;
        point.z = rotated.z;
    });
    compatible_homogeneous_control_nets(curve, &segment)
        .then_some(Line(Point2::new(u0, angle), Point2::new(u1, angle)))
}

fn exact_revolution_isoparametric_boundary<C, D>(
    curve: &D,
    surface: &RevolutionSurface<C>,
) -> Option<ParameterCurve<Line<Point2>, RevolutionSurface<C>>>
where
    C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve,
    D: ParametricCurve3D<Point = Point3> + BoundedCurve + TryIntoHomogeneousBsplineCurve,
{
    let homogeneous = curve.try_into_homogeneous_bspline_curve()?;
    let boundary = exact_revolution_isoparametric_boundary_from_homogeneous(
        &homogeneous,
        curve.front(),
        curve.back(),
        surface,
    )?;
    Some(ParameterCurve::new(boundary, surface.clone()))
}

type ProcessorRevolutionBoundary<C> =
    ParameterCurve<Line<Point2>, Processor<RevolutionSurface<C>, Matrix4>>;

fn exact_processor_revolution_isoparametric_boundary<C, D>(
    curve: &D,
    surface: &Processor<RevolutionSurface<C>, Matrix4>,
) -> Option<ProcessorRevolutionBoundary<C>>
where
    C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve,
    D: ParametricCurve3D<Point = Point3> + BoundedCurve + TryIntoHomogeneousBsplineCurve,
{
    let inverse = surface.transform().invert()?;
    let mut homogeneous = curve.try_into_homogeneous_bspline_curve()?;
    homogeneous
        .control_points_mut()
        .for_each(|point| *point = inverse * *point);
    let mut boundary = exact_revolution_isoparametric_boundary_from_homogeneous(
        &homogeneous,
        inverse.transform_point(curve.front()),
        inverse.transform_point(curve.back()),
        surface.entity(),
    )?;
    if !surface.orientation() {
        boundary = Line(
            Point2::new(boundary.0.y, boundary.0.x),
            Point2::new(boundary.1.y, boundary.1.x),
        );
    }
    Some(ParameterCurve::new(boundary, surface.clone()))
}

impl ExactParameterBoundary2D<Plane> for Line<Point3> {
    type BoundaryCurve = ParameterCurve<Line<Point2>, Plane>;

    fn exact_parameter_boundary_2d(&self, surface: &Plane) -> Option<Self::BoundaryCurve> {
        exact_line_boundary_on_plane(self, surface)
    }
}

impl ExactParameterBoundary2D<BsplineSurface<Point3>> for Line<Point3> {
    type BoundaryCurve = ParameterCurve<Line<Point2>, BsplineSurface<Point3>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &BsplineSurface<Point3>,
    ) -> Option<Self::BoundaryCurve> {
        (surface.degrees() == (1, 1))
            .then(|| exact_line_boundary_on_affine_surface(self, surface))
            .flatten()
    }
}

impl ExactParameterBoundary2D<BsplineSurface<Point3>> for BsplineCurve<Point3> {
    type BoundaryCurve = ParameterCurve<Line<Point2>, BsplineSurface<Point3>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &BsplineSurface<Point3>,
    ) -> Option<Self::BoundaryCurve> {
        exact_bspline_boundary_on_surface(self, surface)
    }
}

impl ExactParameterBoundary2D<NurbsSurface<Vector4>> for Line<Point3> {
    type BoundaryCurve = ParameterCurve<Line<Point2>, NurbsSurface<Vector4>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &NurbsSurface<Vector4>,
    ) -> Option<Self::BoundaryCurve> {
        exact_line_boundary_on_homogeneous_extrusion_surface(self, surface).or_else(|| {
            (surface.degrees() == (1, 1))
                .then(|| exact_line_boundary_on_affine_surface(self, surface))
                .flatten()
        })
    }
}

impl ExactParameterBoundary2D<NurbsSurface<Vector4>> for NurbsCurve<Vector4> {
    type BoundaryCurve = ParameterCurve<Line<Point2>, NurbsSurface<Vector4>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &NurbsSurface<Vector4>,
    ) -> Option<Self::BoundaryCurve> {
        exact_nurbs_boundary_on_surface(self, surface)
    }
}

impl<C> ExactParameterBoundary2D<RevolutionSurface<C>> for Line<Point3>
where C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve
{
    type BoundaryCurve = ParameterCurve<Line<Point2>, RevolutionSurface<C>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &RevolutionSurface<C>,
    ) -> Option<Self::BoundaryCurve> {
        exact_revolution_isoparametric_boundary(self, surface)
    }
}

impl<C> ExactParameterBoundary2D<RevolutionSurface<C>> for BsplineCurve<Point3>
where C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve
{
    type BoundaryCurve = ParameterCurve<Line<Point2>, RevolutionSurface<C>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &RevolutionSurface<C>,
    ) -> Option<Self::BoundaryCurve> {
        exact_revolution_isoparametric_boundary(self, surface)
    }
}

impl<C> ExactParameterBoundary2D<RevolutionSurface<C>> for NurbsCurve<Vector4>
where C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve
{
    type BoundaryCurve = ParameterCurve<Line<Point2>, RevolutionSurface<C>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &RevolutionSurface<C>,
    ) -> Option<Self::BoundaryCurve> {
        exact_revolution_isoparametric_boundary(self, surface)
    }
}

impl<C> ExactParameterBoundary2D<Processor<RevolutionSurface<C>, Matrix4>> for Line<Point3>
where C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve
{
    type BoundaryCurve = ParameterCurve<Line<Point2>, Processor<RevolutionSurface<C>, Matrix4>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &Processor<RevolutionSurface<C>, Matrix4>,
    ) -> Option<Self::BoundaryCurve> {
        exact_processor_revolution_isoparametric_boundary(self, surface)
    }
}

impl<C> ExactParameterBoundary2D<Processor<RevolutionSurface<C>, Matrix4>> for BsplineCurve<Point3>
where C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve
{
    type BoundaryCurve = ParameterCurve<Line<Point2>, Processor<RevolutionSurface<C>, Matrix4>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &Processor<RevolutionSurface<C>, Matrix4>,
    ) -> Option<Self::BoundaryCurve> {
        exact_processor_revolution_isoparametric_boundary(self, surface)
    }
}

impl<C> ExactParameterBoundary2D<Processor<RevolutionSurface<C>, Matrix4>> for NurbsCurve<Vector4>
where C: Clone + ParametricCurve3D + BoundedCurve + TryIntoHomogeneousBsplineCurve
{
    type BoundaryCurve = ParameterCurve<Line<Point2>, Processor<RevolutionSurface<C>, Matrix4>>;

    fn exact_parameter_boundary_2d(
        &self,
        surface: &Processor<RevolutionSurface<C>, Matrix4>,
    ) -> Option<Self::BoundaryCurve> {
        exact_processor_revolution_isoparametric_boundary(self, surface)
    }
}

impl<C, S0, S1, T, S> ExactParameterBoundary2D<S> for SurfaceCurve<C, S0, S1, T, T>
where
    T: ParametricCurve3D<Point = Point3>
        + BoundedCurve<Point = Point3>
        + SearchNearestParameter<CurveParameter, Point = Point3>
        + Cut
        + Invertible
        + Clone,
    S0: PartialEq<S>,
    S1: PartialEq<S>,
{
    type BoundaryCurve = T;

    fn exact_parameter_boundary_2d(&self, surface: &S) -> Option<Self::BoundaryCurve> {
        if self.surface0() == surface {
            self.boundary0().cloned()
        } else if self.surface1() == surface {
            self.boundary1().cloned()
        } else {
            None
        }
    }
}

impl<C, S0, S1, T0, T1, S> ParameterBoundary2D<S> for SurfaceCurve<C, S0, S1, T0, T1>
where
    C: ParametricCurve3D + BoundedCurve,
    S0: ParametricSurface3D
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + PartialEq<S>,
    S1: ParametricSurface3D
        + SearchNearestParameter<SurfaceParameter, Point = Point3>
        + PartialEq<S>,
    T0: ParameterBoundary2D<S> + Clone,
    T1: ParameterBoundary2D<S> + Clone,
{
    fn parameter_boundary_2d(&self, surface: &S, tolerance: f64) -> Option<Vec<Point2>> {
        if self.surface0() == surface {
            self.boundary0()
                .and_then(|boundary| boundary.parameter_boundary_2d(surface, tolerance))
                .or_else(|| {
                    self.try_parameter_division(self.range_tuple(), tolerance)?
                        .0
                        .into_iter()
                        .map(|t| self.search_triple(t, 100).map(|(_, uv0, _)| uv0))
                        .collect()
                })
        } else if self.surface1() == surface {
            self.boundary1()
                .and_then(|boundary| boundary.parameter_boundary_2d(surface, tolerance))
                .or_else(|| {
                    self.try_parameter_division(self.range_tuple(), tolerance)?
                        .0
                        .into_iter()
                        .map(|t| self.search_triple(t, 100).map(|(_, _, uv1)| uv1))
                        .collect()
                })
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests;
