use monstertruck_geometry::prelude::*;
use std::f64::consts::PI;

use super::types::FilletCurve;

pub fn circle_arc(
    point: Point3,
    origin: Point3,
    axis: Vector3,
    angle: Rad<f64>,
    (w0, w1): (f64, f64),
) -> NurbsCurve<Vector4> {
    let origin = origin + (axis.dot(point - origin)) * axis;
    let diag = point - origin;
    let axis_trsf = Matrix4::from_cols(
        diag.extend(0.0),
        axis.cross(diag).extend(0.0),
        axis.extend(0.0),
        origin.to_homogeneous(),
    );
    let mut unit_curve = unit_circle_arc(angle, w0, w1);
    unit_curve.transform_by(axis_trsf);
    unit_curve
}

pub fn circle_arc_by_three_points(
    point0: Vector4,
    point1: Vector4,
    transit: Point3,
) -> NurbsCurve<Vector4> {
    let (w0, w1) = (point0.w, point1.w);
    let point0 = Point3::from_homogeneous(point0);
    let point1 = Point3::from_homogeneous(point1);
    let origin = circum_center(point0, point1, transit);
    let (vec0, vec1) = (point0 - transit, point1 - transit);
    let axis = vec1.cross(vec0).normalize();
    let angle = Rad(PI) - vec0.angle(vec1);
    circle_arc(point0, origin, axis, angle * 2.0, (w0, w1))
}

fn circum_center(pt0: Point3, pt1: Point3, pt2: Point3) -> Point3 {
    let (vec0, vec1) = (pt1 - pt0, pt2 - pt0);
    let (a2, ab, b2) = (vec0.dot(vec0), vec0.dot(vec1), vec1.dot(vec1));
    let det = ab.mul_add(-ab, a2 * b2);
    // Solving the 2×2 Gram system for the circumcenter coefficients gives
    // b²(a² - a·b) and a²(b² - a·b). Factoring the shared terms makes
    // the symmetry explicit and avoids the easy-to-misread product grouping.
    let u = (a2 - ab) * b2;
    let v = (b2 - ab) * a2;
    pt0 + u / (2.0 * det) * vec0 + v / (2.0 * det) * vec1
}

fn unit_circle_arc(angle: Rad<f64>, w0: f64, w1: f64) -> NurbsCurve<Vector4> {
    let p0 = Vector3::new(w0, 0.0, w0);
    let p1 = Vector3::new(Rad::cos(angle), Rad::sin(angle), 1.0) * w1;
    let (sin_half, cos_half) = (angle / 2.0).sin_cos();
    let middle_weight = f64::sqrt(w0 * w1);
    let c = Vector3::new(
        middle_weight * cos_half,
        middle_weight * sin_half,
        middle_weight * cos_half,
    );
    let mut curve = BsplineCurve::new_unchecked(
        KnotVector::bezier_knot(2),
        vec![
            Vector4::new(p0.x, p0.y, 0.0, p0.z),
            Vector4::new(c.x, c.y, 0.0, c.z),
            Vector4::new(p1.x, p1.y, 0.0, p1.z),
        ],
    );

    curve.add_knot(0.25);
    curve.add_knot(0.5);
    curve.add_knot(0.75);
    NurbsCurve::new(curve)
}

#[inline]
pub fn unit_circle_knot_vec() -> KnotVector { KnotVector::uniform_knot(2, 4) }

#[inline]
pub const fn number_of_cpts_of_unit_circle() -> usize { 6 }

#[test]
fn unit_circle_info() {
    let uc = unit_circle_arc(Rad(PI), 1.0, 1.0);
    assert_eq!(uc.knot_vector(), &unit_circle_knot_vec());
    assert_eq!(uc.control_points().len(), number_of_cpts_of_unit_circle());
}

#[cfg(test)]
proptest::proptest! {
    #[test]
    fn test_unit_circle(
        angle in 0.1f64..(1.5 * PI),
        w0 in 0.1f64..5.0,
        w1 in 0.1f64..5.0,
    ) {
        use proptest::prelude::*;
        let uc = unit_circle_arc(Rad(angle), w0, w1);
        const N: usize = 10;
        for i in 0..=N {
            let t = i as f64 / N as f64;
            let p = uc.evaluate(t).to_vec();
            let v = uc.derivative(t);
            prop_assert_near!(p.magnitude2(), 1.0, "{w0} {w1} {p:?} {angle}");
            prop_assert!(p.z.so_small2());
            prop_assert!(p.x * v.y - p.y * v.x > 0.0, "minus area {:?}", uc.control_point(1));
        }
        prop_assert_near!(uc.evaluate(0.0), Point3::new(1.0, 0.0, 0.0));
        prop_assert_near!(uc.evaluate(1.0), Point3::new(f64::cos(angle), f64::sin(angle), 0.0));
    }
}

fn interpolate_bezier(points: &[Vector4]) -> Option<BsplineCurve<Vector4>> {
    let n = points.len().checked_sub(1)?;
    if n == 0 {
        return None;
    }
    let parameter_points = points
        .iter()
        .enumerate()
        .map(|(i, p)| (crate::count_as_f64(i) / crate::count_as_f64(n), *p))
        .collect::<Vec<_>>();
    BsplineCurve::interpolate(KnotVector::bezier_knot(n), parameter_points).ok()
}

pub fn composite_line_bezier(
    line: Line<Point2>,
    surface: &NurbsSurface<Vector4>,
) -> Option<BsplineCurve<Vector4>> {
    let curve = ParameterCurve::new(line, surface.non_rationalized());
    let degree = surface.udegree() + surface.vdegree();
    if degree == 0 {
        return None;
    }
    let points = (0..=degree)
        .map(|i| curve.evaluate(crate::count_as_f64(i) / crate::count_as_f64(degree)))
        .collect::<Vec<_>>();
    interpolate_bezier(&points)
}

#[derive(Clone, Copy, Debug)]
pub struct RelaySphere {
    /// Center of fillet circle.
    #[expect(
        dead_code,
        reason = "retained with the relay-sphere state for fillet diagnostics and inspection"
    )]
    pub(super) center: Point3,
    /// Contact point of sphere and surface0, 3d coordinate and parameter.
    pub(super) contact0: (Point3, Point2),
    /// Contact point of sphere and surface1, 3d coordinate and parameter.
    pub(super) contact1: (Point3, Point2),
    /// Transit of fillet circle.
    pub(super) transit: Point3,
}

impl RelaySphere {
    /// Contact point on planes and sphere with radius `radius`.
    /// Returns `(center, contact_point0, contact_point1)`.
    fn contact_point(
        // point and derivation
        point_on_curve: (Point3, Vector3),
        // origin and normal
        plane0: (Point3, Vector3),
        // origin and normal
        plane1: (Point3, Vector3),
        radius: f64,
        flip_side: bool,
    ) -> Option<(Point3, Point3, Point3)> {
        let ((p, der), (p0, n0), (p1, n1)) = (point_on_curve, plane0, plane1);
        let n = n0.cross(n1);
        // Which side of the edge the sphere sits on is keyed off the seam
        // direction relative to the surface normals; `flip_side` selects the
        // opposite solution for shells wound the other way (swapping the
        // faces cannot flip it: that negates both the cross product and
        // `der`, leaving the sign invariant).
        let orientation = if flip_side { -1.0 } else { 1.0 };
        let sign = f64::signum(n.dot(der)) * orientation;
        let mat = Matrix3::from_cols(n, n0, n1).transpose();
        let vec = Vector3::new(
            n.dot(p.to_vec()),
            f64::mul_add(sign, -radius, n0.dot(p0.to_vec())),
            f64::mul_add(sign, -radius, n1.dot(p1.to_vec())),
        );
        let center = Point3::from_vec(mat.invert()? * vec);
        let q0 = center + sign * radius * n0;
        let q1 = center + sign * radius * n1;
        Some((center, q0, q1))
    }

    fn next_point(
        surface: &(impl ParametricSurface3D + SearchParameter<SurfaceParameter, Point = Point3>),
        parameters: (f64, f64),
        points: (Point3, Point3),
    ) -> Option<(Point3, (f64, f64))> {
        let first_tangent = surface.derivative_u(parameters.0, parameters.1);
        let second_tangent = surface.derivative_v(parameters.0, parameters.1);
        let displacement = points.1 - points.0;
        let first_norm = first_tangent.dot(first_tangent);
        let cross_term = first_tangent.dot(second_tangent);
        let second_norm = second_tangent.dot(second_tangent);
        let matrix = Matrix2::new(first_norm, cross_term, cross_term, second_norm);
        let projection = Vector2::new(
            first_tangent.dot(displacement),
            second_tangent.dot(displacement),
        );
        let correction = matrix.invert()? * projection;
        let updated = (parameters.0 + correction.x, parameters.1 + correction.y);
        Some((surface.evaluate(updated.0, updated.1), updated))
    }

    pub(super) fn generate(
        point_on_curve: (Point3, Vector3),
        surface0: &(
             impl ParametricSurface3D
             + SearchParameter<SurfaceParameter, Point = Point3>
             + SearchNearestParameter<SurfaceParameter, Point = Point3>
         ),
        surface1: &(
             impl ParametricSurface3D
             + SearchParameter<SurfaceParameter, Point = Point3>
             + SearchNearestParameter<SurfaceParameter, Point = Point3>
         ),
        radius: f64,
        flip_side: bool,
    ) -> Option<Self> {
        let (p, der) = point_on_curve;
        let (mut p0, mut p1) = (p, p);
        let (mut u0, mut v0) = surface0
            .search_nearest_parameter(p0, None, 50)
            .or_else(|| surface0.search_parameter(p0, None, 10))?;
        let (mut u1, mut v1) = surface1
            .search_nearest_parameter(p1, None, 50)
            .or_else(|| surface1.search_parameter(p1, None, 10))?;
        let mut converged = false;
        let mut center = Point3::origin();
        for _ in 0..100 {
            let (n0, n1) = (surface0.normal(u0, v0), surface1.normal(u1, v1));
            let (c, q0, q1) = Self::contact_point((p, der), (p0, n0), (p1, n1), radius, flip_side)?;
            if p0.near(&q0) && p1.near(&q1) {
                center = c;
                converged = true;
                break;
            }
            (p0, (u0, v0)) = Self::next_point(surface0, (u0, v0), (p0, q0))?;
            (p1, (u1, v1)) = Self::next_point(surface1, (u1, v1), (p1, q1))?;
        }
        if !converged {
            return None;
        }
        Some(Self {
            center,
            contact0: (p0, (u0, v0).into()),
            contact1: (p1, (u1, v1).into()),
            transit: center + radius * (p - center).normalize(),
        })
    }

    pub(super) fn fillet_wire(self) -> NurbsCurve<Vector4> {
        circle_arc_by_three_points(
            self.contact0.0.to_homogeneous(),
            self.contact1.0.to_homogeneous(),
            self.transit,
        )
    }
}

pub fn relay_spheres(
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    curve: &impl FilletCurve,
    division: usize,
    radius: impl Fn(f64) -> f64,
    extend: bool,
    flip_side: bool,
) -> Option<Vec<RelaySphere>> {
    let (t0, t1) = curve.range_tuple();
    let generator = move |i: isize| {
        let a = crate::signed_index_as_f64(i) / crate::count_as_f64(division);
        let delta = 0.25 / crate::count_as_f64(division);
        [0.0, -delta, delta, -2.0 * delta, 2.0 * delta]
            .into_iter()
            .find_map(|offset| {
                let b = a + offset;
                let t = b.mul_add(t1, (1.0 - b) * t0);
                RelaySphere::generate(
                    (curve.evaluate(t), curve.derivative(t)),
                    surface0,
                    surface1,
                    radius(b),
                    flip_side,
                )
            })
    };
    let division = division.cast_signed();
    let range = if extend {
        -1..=(division + 1)
    } else {
        0..=division
    };
    range.map(generator).collect()
}

#[expect(
    dead_code,
    reason = "reference construction retained for comparison with the active custom fillet implementation"
)]
fn almost_fillet_patch(rs0: RelaySphere, rs1: RelaySphere) -> BsplineSurface<Vector4> {
    let nurbs0 = rs0.fillet_wire();
    let nurbs1 = rs1.fillet_wire();
    let knot_vecs = (KnotVector::bezier_knot(1), nurbs0.knot_vector().clone());
    let control_points = vec![
        nurbs0.control_points().clone(),
        nurbs1.control_points().clone(),
    ];
    BsplineSurface::new_unchecked(knot_vecs, control_points)
}

pub fn expand_fillet(
    relay_spheres: &[RelaySphere],
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
) -> Option<NurbsSurface<Vector4>> {
    const LEN: usize = number_of_cpts_of_unit_circle();
    if relay_spheres.len() < 2 {
        return None;
    }
    let mut collectors = vec![CurveCollector::<BsplineCurve<_>>::Singleton; LEN];
    for (n, u) in relay_spheres.windows(2).enumerate() {
        let transit_line = Line(u[0].transit, u[1].transit);
        let mut bezier0 = {
            let line0 = Line(u[0].contact0.1, u[1].contact0.1);
            composite_line_bezier(line0, surface0)?
        };
        let mut bezier1 = {
            let line1 = Line(u[0].contact1.1, u[1].contact1.1);
            composite_line_bezier(line1, surface1)?
        };
        BsplineCurve::syncro_degree(&mut bezier0, &mut bezier1);

        let (cpts0, cpts1) = (bezier0.control_points(), bezier1.control_points());
        let iter = cpts0.iter().zip(cpts1).enumerate();
        let control_points = iter
            .map(|(i, (p0, p1))| {
                let t = crate::count_as_f64(i) / crate::count_as_f64(cpts0.len() - 1);
                let transit = transit_line.evaluate(t);
                circle_arc_by_three_points(*p0, *p1, transit)
            })
            .map(|wire| BsplineCurve::from(wire).destruct().1)
            .collect();
        let homo_surface = BsplineSurface::new_unchecked(
            (bezier0.knot_vector().clone(), unit_circle_knot_vec()),
            control_points,
        );

        for (i, collector) in collectors.iter_mut().enumerate() {
            let mut curve = homo_surface.curve_u(i);
            curve.knot_translate(crate::count_as_f64(n));
            collector.try_concat(&curve).ok()?;
        }
    }
    let curves = collectors
        .into_iter()
        .map(Option::<BsplineCurve<Vector4>>::from)
        .collect::<Option<Vec<_>>>()?;
    let knot_vector_v = curves[0].knot_vector().clone();
    let control_points = curves
        .into_iter()
        .map(|curve| curve.destruct().1)
        .collect::<Vec<_>>();

    let knot_vecs = (unit_circle_knot_vec(), knot_vector_v);
    let mut bsp_surface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    bsp_surface.knot_normalize();
    Some(NurbsSurface::new(bsp_surface))
}

pub fn asymmetric_chamfer_contacts(
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    curve: &impl FilletCurve,
    division: usize,
    setbacks_mm: [f64; 2],
    extend: bool,
    flip_side: bool,
) -> Option<Vec<RelaySphere>> {
    let [first_setback_mm, second_setback_mm] = setbacks_mm;
    let (t0, t1) = curve.range_tuple();
    let reference_radius = first_setback_mm.min(second_setback_mm);
    let generator = move |i: isize| {
        let a = crate::signed_index_as_f64(i) / crate::count_as_f64(division);
        let delta = 0.25 / crate::count_as_f64(division);
        [0.0, -delta, delta, -2.0 * delta, 2.0 * delta]
            .into_iter()
            .find_map(|offset| {
                let b = a + offset;
                let t = b.mul_add(t1, (1.0 - b) * t0);
                let point = curve.evaluate(t);
                let reference = RelaySphere::generate(
                    (point, curve.derivative(t)),
                    surface0,
                    surface1,
                    reference_radius,
                    flip_side,
                )?;
                let direction0 = (reference.contact0.0 - point).normalize();
                let direction1 = (reference.contact1.0 - point).normalize();
                let target0 = point + first_setback_mm * direction0;
                let target1 = point + second_setback_mm * direction1;
                let parameter0 = surface0
                    .search_parameter(target0, None, 20)
                    .or_else(|| surface0.search_nearest_parameter(target0, None, 50))?;
                let parameter1 = surface1
                    .search_parameter(target1, None, 20)
                    .or_else(|| surface1.search_nearest_parameter(target1, None, 50))?;
                let contact0 = surface0.evaluate(parameter0.0, parameter0.1);
                let contact1 = surface1.evaluate(parameter1.0, parameter1.1);
                if !contact0.near(&target0) || !contact1.near(&target1) {
                    return None;
                }
                Some(RelaySphere {
                    center: point,
                    contact0: (contact0, parameter0.into()),
                    contact1: (contact1, parameter1.into()),
                    transit: contact0.midpoint(contact1),
                })
            })
    };
    let division = division.cast_signed();
    let range = if extend {
        -1..=(division + 1)
    } else {
        0..=division
    };
    range.map(generator).collect()
}

pub fn expand_chamfer(
    relay_spheres: &[RelaySphere],
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
) -> Option<NurbsSurface<Vector4>> {
    const CHAMFER_CPTS: usize = 2; // degree-1 cross-section: 2 control points
    if relay_spheres.len() < 2 {
        return None;
    }
    let mut collectors = vec![CurveCollector::<BsplineCurve<_>>::Singleton; CHAMFER_CPTS];
    for (n, u) in relay_spheres.windows(2).enumerate() {
        let mut bezier0 = {
            let line0 = Line(u[0].contact0.1, u[1].contact0.1);
            composite_line_bezier(line0, surface0)?
        };
        let mut bezier1 = {
            let line1 = Line(u[0].contact1.1, u[1].contact1.1);
            composite_line_bezier(line1, surface1)?
        };
        BsplineCurve::syncro_degree(&mut bezier0, &mut bezier1);

        let (cpts0, cpts1) = (bezier0.control_points(), bezier1.control_points());
        // Degree-1 cross-section: straight line from contact0 to contact1
        let control_points = cpts0
            .iter()
            .zip(cpts1)
            .map(|(p0, p1)| BsplineCurve::new_unchecked(KnotVector::bezier_knot(1), vec![*p0, *p1]))
            .map(|wire| wire.destruct().1)
            .collect();
        let homo_surface = BsplineSurface::new_unchecked(
            (bezier0.knot_vector().clone(), KnotVector::bezier_knot(1)),
            control_points,
        );

        for (i, collector) in collectors.iter_mut().enumerate() {
            let mut curve = homo_surface.curve_u(i);
            curve.knot_translate(crate::count_as_f64(n));
            collector.try_concat(&curve).ok()?;
        }
    }
    let curves = collectors
        .into_iter()
        .map(Option::<BsplineCurve<Vector4>>::from)
        .collect::<Option<Vec<_>>>()?;
    let knot_vector_v = curves[0].knot_vector().clone();
    let control_points = curves
        .into_iter()
        .map(|curve| curve.destruct().1)
        .collect::<Vec<_>>();

    let knot_vecs = (KnotVector::bezier_knot(1), knot_vector_v);
    let mut bsp_surface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    bsp_surface.knot_normalize();
    Some(NurbsSurface::new(bsp_surface))
}

pub fn expand_ridge(
    relay_spheres: &[RelaySphere],
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
) -> Option<NurbsSurface<Vector4>> {
    const RIDGE_CPTS: usize = 3; // degree-1 piecewise: p0, transit, p1
    if relay_spheres.len() < 2 {
        return None;
    }
    let mut collectors = vec![CurveCollector::<BsplineCurve<_>>::Singleton; RIDGE_CPTS];
    for (n, u) in relay_spheres.windows(2).enumerate() {
        let transit_line = Line(u[0].transit, u[1].transit);
        let mut bezier0 = {
            let line0 = Line(u[0].contact0.1, u[1].contact0.1);
            composite_line_bezier(line0, surface0)?
        };
        let mut bezier1 = {
            let line1 = Line(u[0].contact1.1, u[1].contact1.1);
            composite_line_bezier(line1, surface1)?
        };
        BsplineCurve::syncro_degree(&mut bezier0, &mut bezier1);

        let (cpts0, cpts1) = (bezier0.control_points(), bezier1.control_points());
        let control_points = cpts0
            .iter()
            .zip(cpts1)
            .enumerate()
            .map(|(i, (p0, p1))| {
                let t = crate::count_as_f64(i) / crate::count_as_f64(cpts0.len() - 1);
                let transit = transit_line.evaluate(t);
                BsplineCurve::new_unchecked(
                    KnotVector::uniform_knot(1, 2),
                    vec![*p0, transit.to_homogeneous(), *p1],
                )
            })
            .map(|wire| wire.destruct().1)
            .collect();
        let homo_surface = BsplineSurface::new_unchecked(
            (
                bezier0.knot_vector().clone(),
                KnotVector::uniform_knot(1, 2),
            ),
            control_points,
        );

        for (i, collector) in collectors.iter_mut().enumerate() {
            let mut curve = homo_surface.curve_u(i);
            curve.knot_translate(crate::count_as_f64(n));
            collector.try_concat(&curve).ok()?;
        }
    }
    let curves = collectors
        .into_iter()
        .map(Option::<BsplineCurve<Vector4>>::from)
        .collect::<Option<Vec<_>>>()?;
    let knot_vector_v = curves[0].knot_vector().clone();
    let control_points = curves
        .into_iter()
        .map(|curve| curve.destruct().1)
        .collect::<Vec<_>>();

    let knot_vecs = (KnotVector::uniform_knot(1, 2), knot_vector_v);
    let mut bsp_surface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    bsp_surface.knot_normalize();
    Some(NurbsSurface::new(bsp_surface))
}

pub fn ridge_fillet_surface(
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    curve: &impl FilletCurve,
    division: usize,
    radius: impl Fn(f64) -> f64,
    extend: bool,
    flip_side: bool,
) -> Option<NurbsSurface<Vector4>> {
    let relay_spheres = relay_spheres(
        surface0, surface1, curve, division, radius, extend, flip_side,
    )?;
    expand_ridge(&relay_spheres, surface0, surface1)
}

fn map_profile_to_3d(
    p0: &Vector4,
    p1: &Vector4,
    transit: Point3,
    profile: &BsplineCurve<Point2>,
) -> BsplineCurve<Vector4> {
    let pt0 = Point3::from_homogeneous(*p0);
    let pt1 = Point3::from_homogeneous(*p1);
    let mid = pt0.midpoint(pt1);
    let toward_transit = transit - mid;
    let cpts: Vec<Vector4> = profile
        .control_points()
        .iter()
        .map(|cp| {
            let px = cp.x;
            let py = cp.y;
            let pos = pt0 + px * (pt1 - pt0) + py * toward_transit;
            let w = f64::mul_add(px, p1.w, (1.0 - px) * p0.w);
            Vector4::new(pos.x * w, pos.y * w, pos.z * w, w)
        })
        .collect();
    BsplineCurve::new_unchecked(profile.knot_vector().clone(), cpts)
}

pub fn expand_custom(
    relay_spheres: &[RelaySphere],
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    profile: &BsplineCurve<Point2>,
) -> Option<NurbsSurface<Vector4>> {
    let num_cpts = profile.control_points().len();
    if relay_spheres.len() < 2 {
        return None;
    }
    let mut collectors = vec![CurveCollector::<BsplineCurve<_>>::Singleton; num_cpts];
    for (n, u) in relay_spheres.windows(2).enumerate() {
        let transit_line = Line(u[0].transit, u[1].transit);
        let mut bezier0 = {
            let line0 = Line(u[0].contact0.1, u[1].contact0.1);
            composite_line_bezier(line0, surface0)?
        };
        let mut bezier1 = {
            let line1 = Line(u[0].contact1.1, u[1].contact1.1);
            composite_line_bezier(line1, surface1)?
        };
        BsplineCurve::syncro_degree(&mut bezier0, &mut bezier1);

        let (cpts0, cpts1) = (bezier0.control_points(), bezier1.control_points());
        let control_points = cpts0
            .iter()
            .zip(cpts1)
            .enumerate()
            .map(|(i, (p0, p1))| {
                let t = crate::count_as_f64(i) / crate::count_as_f64(cpts0.len() - 1);
                let transit = transit_line.evaluate(t);
                map_profile_to_3d(p0, p1, transit, profile)
            })
            .map(|wire| wire.destruct().1)
            .collect();
        let homo_surface = BsplineSurface::new_unchecked(
            (bezier0.knot_vector().clone(), profile.knot_vector().clone()),
            control_points,
        );

        for (i, collector) in collectors.iter_mut().enumerate() {
            let mut curve = homo_surface.curve_u(i);
            curve.knot_translate(crate::count_as_f64(n));
            collector.try_concat(&curve).ok()?;
        }
    }
    let curves = collectors
        .into_iter()
        .map(Option::<BsplineCurve<Vector4>>::from)
        .collect::<Option<Vec<_>>>()?;
    let knot_vector_v = curves[0].knot_vector().clone();
    let control_points = curves
        .into_iter()
        .map(|curve| curve.destruct().1)
        .collect::<Vec<_>>();

    let knot_vecs = (profile.knot_vector().clone(), knot_vector_v);
    let mut bsp_surface = BsplineSurface::new_unchecked(knot_vecs, control_points);
    bsp_surface.knot_normalize();
    Some(NurbsSurface::new(bsp_surface))
}

// One argument past the lint threshold: `radius: impl Fn` resists struct
// bundling and the sibling builders (`chamfer_fillet_surface`,
// `rolling_ball_fillet_surface`) share this argument shape, so a params struct
// here would either desync those signatures or ripple through `relay_spheres`.
#[expect(
    clippy::too_many_arguments,
    reason = "arguments mirror the two surfaces, relay spheres, and boundary derivatives without an artificial parameter bag"
)]
pub fn custom_fillet_surface(
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    curve: &impl FilletCurve,
    division: usize,
    radius: impl Fn(f64) -> f64,
    extend: bool,
    flip_side: bool,
    profile: &BsplineCurve<Point2>,
) -> Option<NurbsSurface<Vector4>> {
    let relay_spheres = relay_spheres(
        surface0, surface1, curve, division, radius, extend, flip_side,
    )?;
    expand_custom(&relay_spheres, surface0, surface1, profile)
}

pub fn chamfer_fillet_surface(
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    curve: &impl FilletCurve,
    division: usize,
    radius: impl Fn(f64) -> f64,
    extend: bool,
    flip_side: bool,
) -> Option<NurbsSurface<Vector4>> {
    let relay_spheres = relay_spheres(
        surface0, surface1, curve, division, radius, extend, flip_side,
    )?;
    expand_chamfer(&relay_spheres, surface0, surface1)
}

/// Orientation of `curve` should be compatible with orientation of `surface0`.
pub fn rolling_ball_fillet_surface(
    surface0: &NurbsSurface<Vector4>,
    surface1: &NurbsSurface<Vector4>,
    curve: &impl FilletCurve,
    division: usize,
    radius: impl Fn(f64) -> f64,
    extend: bool,
    flip_side: bool,
) -> Option<NurbsSurface<Vector4>> {
    let relay_spheres = relay_spheres(
        surface0, surface1, curve, division, radius, extend, flip_side,
    )?;
    expand_fillet(&relay_spheres, surface0, surface1)
}
