//! `PolyBoundary`: assembling projected pieces into closed, consistently
//! wound parameter-space loops, testing point inclusion, cutting isoparametric
//! curves, and feeding the constrained Delaunay triangulation.
//!
//! Split out of the module file so the source stays readable; the module path
//! is unchanged.

use super::{
    Cdt, FixedVertexHandle, HashGen, HashMap, InnerSpace, IsoparametricCurveOptions, Itertools,
    MetricSpace, Origin, Point2, Point3, PolyBoundaryPiece, PreMeshableSurface, SPoint2,
    SurfacePoint, TOLERANCE, Tolerance, Triangulation, UvKey, Vector2, align_full_period_open_pair,
    align_periodic_open_pair, close_open_periodic_curve_to_singular, compact_periodic_closed_loop,
    connect_edges, mesh, normalize_range, open_periodic_closed_loop, periodic_axis_full_span,
    periodic_open_axis_full_span, polyline_on_surface, surface_point_with_cache,
};
use mesh::spade_round;
use std::{f64::consts::TAU, iter, mem};

const ISOPARAM_BOUNDARY_SEARCH_ITERATIONS: usize = 32;

const fn end_points<T: Copy>(slice: &[T]) -> (T, T) { (slice[0], slice[slice.len() - 1]) }

#[derive(Debug, Clone)]
pub(in crate::tessellation::triangulation) struct PolyBoundary {
    pub(in crate::tessellation::triangulation) loops: Vec<Vec<SurfacePoint>>,
    /// UV-space axis-aligned bounding box for cheap rejection in `include()`.
    pub(in crate::tessellation::triangulation) uv_min: Point2,
    pub(in crate::tessellation::triangulation) uv_max: Point2,
}

impl Default for PolyBoundary {
    fn default() -> Self {
        Self {
            loops: Vec::new(),
            uv_min: Point2::new(f64::INFINITY, f64::INFINITY),
            uv_max: Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
        }
    }
}

fn loop_signed_area(curve: &[SurfacePoint]) -> f64 {
    curve
        .iter()
        .circular_tuple_windows()
        .fold(0.0, |sum, (p, q)| f64::mul_add(q.x + p.x, q.y - p.y, sum))
}

fn loop_orientation(curve: &[SurfacePoint]) -> bool { loop_signed_area(curve) > 0.0 }

fn normalize_closed_loop_winding(closed: &mut [Vec<SurfacePoint>]) {
    let outer_index = closed
        .iter()
        .enumerate()
        .max_by(|(_, lhs), (_, rhs)| {
            loop_signed_area(lhs)
                .abs()
                .total_cmp(&loop_signed_area(rhs).abs())
        })
        .map(|(index, _)| index);
    if let Some(outer_index) = outer_index {
        closed.iter_mut().enumerate().for_each(|(index, curve)| {
            let should_be_positive = index == outer_index;
            if loop_orientation(curve) != should_be_positive {
                curve.reverse();
            }
        });
    }
}

fn open_pair_connector_score(
    p0: SurfacePoint,
    p1: SurfacePoint,
    q0: SurfacePoint,
    q1: SurfacePoint,
) -> f64 {
    let p1q0 = p1.uv - q0.uv;
    let q1p0 = q1.uv - p0.uv;
    q1p0.y.mul_add(
        q1p0.y,
        q1p0.x
            .mul_add(q1p0.x, p1q0.y.mul_add(p1q0.y, p1q0.x * p1q0.x)),
    )
}

pub(in crate::tessellation::triangulation) fn boundary_segment_parameter(
    point: Point2,
    front: Point2,
    back: Point2,
) -> Option<f64> {
    let segment = back - front;
    let denom = segment.dot(segment);
    if denom.so_small() {
        None
    } else {
        let offset = point - front;
        let parameter = offset.dot(segment) / denom;
        let projected = front + segment * parameter;
        (parameter > 0.0 && parameter < 1.0 && projected.distance(point) <= 1.0e-9)
            .then_some(parameter)
    }
}

fn push_isoparam_uv(curve: &mut Vec<Point2>, uv: Point2) {
    if curve.last().is_none_or(|previous| !previous.near(&uv)) {
        curve.push(uv);
    }
}

fn push_finished_isoparam_curve(curves: &mut Vec<Vec<Point2>>, curve: &mut Vec<Point2>) {
    if curve.len() >= 2 {
        curves.push(mem::take(curve));
    } else {
        curve.clear();
    }
}

fn add_direct_cdt_constraint(
    triangulation: &mut Cdt,
    added_constraints: &mut usize,
    skipped_constraints: &mut usize,
    front: FixedVertexHandle,
    back: FixedVertexHandle,
) -> bool {
    if front == back || !triangulation.can_add_constraint(front, back) {
        *skipped_constraints += 1;
        false
    } else {
        let constraints = triangulation.add_constraint_and_split(front, back, |point| point);
        if constraints.is_empty() {
            *skipped_constraints += 1;
            false
        } else {
            *added_constraints += constraints.len();
            true
        }
    }
}

fn classify_boundary_pieces(
    pieces: Vec<PolyBoundaryPiece>,
    surface: &impl PreMeshableSurface,
) -> (Vec<Vec<SurfacePoint>>, Vec<Vec<SurfacePoint>>) {
    let (mut closed, mut open) = (Vec::new(), Vec::new());
    for PolyBoundaryPiece(mut curve) in pieces {
        compact_periodic_closed_loop(&mut curve, surface);
        if let Some((axis, period)) = periodic_axis_full_span(&curve, surface) {
            open.push(open_periodic_closed_loop(curve, surface, axis, period));
        } else if curve[0].uv.distance(curve[curve.len() - 1].uv) < 1.0e-3 {
            curve.pop();
            closed.push(curve);
        } else {
            open.push(curve);
        }
    }
    open.retain(|curve| curve.len() >= 2);
    closed.retain(|curve| curve.len() >= 2);
    (closed, open)
}

fn close_curve_via_corners(
    surface: &impl PreMeshableSurface,
    curve: Vec<SurfacePoint>,
    corner0: Point2,
    corner1: Point2,
    tolerance: f64,
    point_cache: &mut HashMap<UvKey, Point3>,
) -> Vec<SurfacePoint> {
    let (front, back) = end_points(&curve);
    let first_corner = surface_point_with_cache(surface, corner0, point_cache);
    let second_corner = surface_point_with_cache(surface, corner1, point_cache);
    let edge0 = polyline_on_surface(surface, back, first_corner, tolerance, point_cache);
    let edge1 = polyline_on_surface(surface, first_corner, second_corner, tolerance, point_cache);
    let edge2 = polyline_on_surface(surface, second_corner, front, tolerance, point_cache);
    connect_edges([edge0, edge1, edge2, curve])
}

fn close_single_open_boundary(
    mut curve: Vec<SurfacePoint>,
    surface: &impl PreMeshableSurface,
    tolerance: f64,
    point_cache: &mut HashMap<UvKey, Point3>,
) -> Option<Vec<SurfacePoint>> {
    if let Some((axis, _)) = periodic_open_axis_full_span(&curve, surface)
        && let Some(singular_curve) = close_open_periodic_curve_to_singular(
            curve.clone(),
            surface,
            axis,
            tolerance,
            point_cache,
        )
    {
        return Some(singular_curve);
    }
    let (Some((u0, u1)), Some((v0, v1))) = surface.try_range_tuple() else {
        return None;
    };
    let (front, back) = end_points(&curve);
    if front.x < back.x - TOLERANCE {
        normalize_range(&mut curve, 0, (u0, u1));
        Some(close_curve_via_corners(
            surface,
            curve,
            Point2::new(u1, v1),
            Point2::new(u0, v1),
            tolerance,
            point_cache,
        ))
    } else if back.x < front.x - TOLERANCE {
        normalize_range(&mut curve, 0, (u0, u1));
        Some(close_curve_via_corners(
            surface,
            curve,
            Point2::new(u0, v0),
            Point2::new(u1, v0),
            tolerance,
            point_cache,
        ))
    } else if front.y < back.y - TOLERANCE {
        normalize_range(&mut curve, 1, (v0, v1));
        Some(close_curve_via_corners(
            surface,
            curve,
            Point2::new(u0, v1),
            Point2::new(u0, v0),
            tolerance,
            point_cache,
        ))
    } else if back.y < front.y - TOLERANCE {
        normalize_range(&mut curve, 1, (v0, v1));
        Some(close_curve_via_corners(
            surface,
            curve,
            Point2::new(u1, v0),
            Point2::new(u1, v1),
            tolerance,
            point_cache,
        ))
    } else {
        None
    }
}

fn close_open_boundary_pair(
    mut curve0: Vec<SurfacePoint>,
    mut curve1: Vec<SurfacePoint>,
    surface: &impl PreMeshableSurface,
    tolerance: f64,
    point_cache: &mut HashMap<UvKey, Point3>,
) -> Vec<SurfacePoint> {
    let full_period_pair_aligned = align_full_period_open_pair(&mut curve0, &mut curve1, surface);
    let ((p0, p1), (q0, q1)) = (end_points(&curve0), end_points(&curve1));
    if !full_period_pair_aligned && !p0.x.near(&p1.x) && !q0.x.near(&q1.x) {
        if let Some(period) = surface.period_u() {
            align_periodic_open_pair(&curve0, &mut curve1, 0, period);
        } else if let (Some(urange), _) = surface.try_range_tuple() {
            normalize_range(&mut curve0, 0, urange);
            normalize_range(&mut curve1, 0, urange);
        }
    } else if !full_period_pair_aligned && !p0.y.near(&p1.y) && !q0.y.near(&q1.y) {
        if let Some(period) = surface.period_v() {
            align_periodic_open_pair(&curve0, &mut curve1, 1, period);
        } else if let (_, Some(vrange)) = surface.try_range_tuple() {
            normalize_range(&mut curve0, 1, vrange);
            normalize_range(&mut curve1, 1, vrange);
        }
    }
    let ((p0, p1), (q0, q1)) = (end_points(&curve0), end_points(&curve1));
    if open_pair_connector_score(p0, p1, q1, q0) + TOLERANCE
        < open_pair_connector_score(p0, p1, q0, q1)
    {
        curve1.reverse();
    }
    let ((p0, p1), (q0, q1)) = (end_points(&curve0), end_points(&curve1));
    let edge0 = polyline_on_surface(surface, p1, q0, tolerance, point_cache);
    let edge1 = polyline_on_surface(surface, q1, p0, tolerance, point_cache);
    connect_edges([curve0, edge0, curve1, edge1])
}

fn add_surface_domain_boundary_if_needed(
    closed: &mut Vec<Vec<SurfacePoint>>,
    surface: &impl PreMeshableSurface,
    tolerance: f64,
    point_cache: &mut HashMap<UvKey, Point3>,
) {
    if closed.iter().any(|curve| loop_orientation(curve)) {
        return;
    }
    let (Some((u0, u1)), Some((v0, v1))) = surface.try_range_tuple() else {
        return;
    };
    let corners = [
        surface_point_with_cache(surface, Point2::new(u0, v0), point_cache),
        surface_point_with_cache(surface, Point2::new(u1, v0), point_cache),
        surface_point_with_cache(surface, Point2::new(u1, v1), point_cache),
        surface_point_with_cache(surface, Point2::new(u0, v1), point_cache),
    ];
    let edge0 = polyline_on_surface(surface, corners[0], corners[1], tolerance, point_cache);
    let edge1 = polyline_on_surface(surface, corners[1], corners[2], tolerance, point_cache);
    let edge2 = polyline_on_surface(surface, corners[2], corners[3], tolerance, point_cache);
    let edge3 = polyline_on_surface(surface, corners[3], corners[0], tolerance, point_cache);
    closed.push(connect_edges([edge0, edge1, edge2, edge3]));
}

fn boundary_uv_bounds(closed: &[Vec<SurfacePoint>]) -> (Point2, Point2) {
    let (mut uv_min, mut uv_max) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for point in closed.iter().flatten() {
        uv_min.x = f64::min(uv_min.x, point.x);
        uv_min.y = f64::min(uv_min.y, point.y);
        uv_max.x = f64::max(uv_max.x, point.x);
        uv_max.y = f64::max(uv_max.y, point.y);
    }
    (uv_min, uv_max)
}

impl PolyBoundary {
    pub(in crate::tessellation::triangulation) fn new(
        pieces: Vec<PolyBoundaryPiece>,
        surface: &impl PreMeshableSurface,
        tolerance: f64,
    ) -> Self {
        let (mut closed, mut open) = classify_boundary_pieces(pieces, surface);
        let mut point_cache = HashMap::<UvKey, Point3>::default();
        match open.len() {
            1 => {
                if let Some(curve) =
                    close_single_open_boundary(open.remove(0), surface, tolerance, &mut point_cache)
                {
                    closed.push(curve);
                }
            }
            2 => closed.push(close_open_boundary_pair(
                open.remove(0),
                open.remove(0),
                surface,
                tolerance,
                &mut point_cache,
            )),
            _ => {}
        }
        if closed.len() == 1 && !loop_orientation(&closed[0]) {
            closed[0].reverse();
        }
        add_surface_domain_boundary_if_needed(&mut closed, surface, tolerance, &mut point_cache);
        normalize_closed_loop_winding(&mut closed);
        let (uv_min, uv_max) = boundary_uv_bounds(&closed);
        Self {
            loops: closed,
            uv_min,
            uv_max,
        }
    }

    /// Whether `c` is included in the domain with boundary = `self`.
    pub(in crate::tessellation::triangulation) fn include(&self, c: Point2) -> bool {
        // AABB early reject.
        if c.x < self.uv_min.x || c.x > self.uv_max.x || c.y < self.uv_min.y || c.y > self.uv_max.y
        {
            return false;
        }
        let t = TAU * HashGen::hash1(c);
        let r = Vector2::new(f64::cos(t), f64::sin(t));
        self.loops
            .iter()
            .flat_map(|vec| vec.iter().circular_tuple_windows())
            .try_fold(0_i32, move |counter, (p0, p1)| {
                let a = **p0 - c;
                let b = **p1 - c;
                let s0 = r.y.mul_add(-a.x, r.x * a.y); // v times a.
                let s1 = r.y.mul_add(-b.x, r.x * b.y); // v times b.
                let s2 = a.y.mul_add(-b.x, a.x * b.y); // a times b.
                let x = s2 / (s1 - s0);
                if x.so_small() && s0 * s1 < 0.0 {
                    None
                } else if x > 0.0 && s0 <= 0.0 && s1 > 0.0 {
                    Some(counter + 1)
                } else if x > 0.0 && s0 >= 0.0 && s1 < 0.0 {
                    Some(counter - 1)
                } else {
                    Some(counter)
                }
            })
            .is_some_and(|counter| counter > 0)
    }

    pub(in crate::tessellation::triangulation) fn isoparametric_curves<S: PreMeshableSurface>(
        &self,
        surface: &S,
        options: IsoparametricCurveOptions,
    ) -> Vec<Vec<Point3>> {
        if self.loops.is_empty()
            || options.samples_per_direction == 0
            || options.segments_per_curve == 0
            || self.uv_min.x >= self.uv_max.x
            || self.uv_min.y >= self.uv_max.y
        {
            Vec::new()
        } else {
            (0..=1)
                .flat_map(|axis| {
                    self.isoparametric_axis_values(axis, options.samples_per_direction)
                        .into_iter()
                        .flat_map(move |axis_value| {
                            self.isoparametric_uv_curves(
                                axis,
                                axis_value,
                                options.segments_per_curve,
                            )
                            .into_iter()
                            .map(|curve| {
                                curve
                                    .into_iter()
                                    .map(|uv| surface.evaluate(uv.x, uv.y))
                                    .collect::<Vec<_>>()
                            })
                        })
                })
                .collect()
        }
    }

    fn isoparametric_axis_values(&self, axis: usize, count: usize) -> Vec<f64> {
        let (min, max) = match axis {
            0 => (self.uv_min.x, self.uv_max.x),
            _ => (self.uv_min.y, self.uv_max.y),
        };
        (0..count)
            .map(|index| {
                let parameter =
                    (crate::count_as_f64(index) + 1.0) / (crate::count_as_f64(count) + 1.0);
                f64::mul_add(max - min, parameter, min)
            })
            .collect()
    }

    fn isoparametric_uv_curves(
        &self,
        constant_axis: usize,
        constant_value: f64,
        segments: usize,
    ) -> Vec<Vec<Point2>> {
        let variable_axis = 1 - constant_axis;
        let (min, max) = match variable_axis {
            0 => (self.uv_min.x, self.uv_max.x),
            _ => (self.uv_min.y, self.uv_max.y),
        };
        let samples = (0..=segments)
            .map(|index| {
                let parameter = crate::count_as_f64(index) / crate::count_as_f64(segments);
                let variable_value = f64::mul_add(max - min, parameter, min);
                match constant_axis {
                    0 => Point2::new(constant_value, variable_value),
                    _ => Point2::new(variable_value, constant_value),
                }
            })
            .collect::<Vec<_>>();
        self.trim_uv_samples(samples)
    }

    fn trim_uv_samples(&self, samples: Vec<Point2>) -> Vec<Vec<Point2>> {
        let (_, mut curves, mut current) = samples.into_iter().fold(
            (
                None::<(Point2, bool)>,
                Vec::<Vec<Point2>>::new(),
                Vec::<Point2>::new(),
            ),
            |(previous, mut curves, mut current), uv| {
                let included = self.include(uv);
                if let Some((previous_uv, previous_included)) = previous
                    && previous_included != included
                {
                    let crossing =
                        self.isoparametric_boundary_crossing(previous_uv, uv, previous_included);
                    push_isoparam_uv(&mut current, crossing);
                    if previous_included {
                        push_finished_isoparam_curve(&mut curves, &mut current);
                    }
                    if included {
                        push_isoparam_uv(&mut current, crossing);
                    }
                }
                if included {
                    push_isoparam_uv(&mut current, uv);
                }
                (Some((uv, included)), curves, current)
            },
        );
        push_finished_isoparam_curve(&mut curves, &mut current);
        curves
    }

    fn isoparametric_boundary_crossing(
        &self,
        front: Point2,
        back: Point2,
        front_included: bool,
    ) -> Point2 {
        let (mut included, mut excluded) = if front_included {
            (front, back)
        } else {
            (back, front)
        };
        (0..ISOPARAM_BOUNDARY_SEARCH_ITERATIONS).for_each(|_| {
            let middle = Point2::new(
                f64::midpoint(included.x, excluded.x),
                f64::midpoint(included.y, excluded.y),
            );
            if self.include(middle) {
                included = middle;
            } else {
                excluded = middle;
            }
        });
        Point2::new(
            f64::midpoint(included.x, excluded.x),
            f64::midpoint(included.y, excluded.y),
        )
    }

    /// Inserts points and adds constraint into triangulation.
    pub(in crate::tessellation::triangulation) fn insert_to(
        &self,
        triangulation: &mut Cdt,
        boundary_map: &mut HashMap<FixedVertexHandle, Point3>,
    ) -> (usize, usize, usize) {
        let poly2tri: Vec<_> = self
            .loops
            .iter()
            .flatten()
            .map(|pt| {
                let p = [spade_round(pt.x), spade_round(pt.y)];
                triangulation.insert(SPoint2::from(p)).ok().inspect(|&idx| {
                    boundary_map.insert(idx, pt.point);
                })
            })
            .collect();
        let split_vertices = triangulation
            .vertices()
            .map(|vertex| {
                let point = vertex.position();
                (vertex.fix(), Point2::new(point.x, point.y))
            })
            .collect::<Vec<_>>();
        let mut counter = 0;
        let mut added_constraints = 0usize;
        let mut skipped_constraints = 0usize;
        let mut add_constraint = |front: FixedVertexHandle, back: FixedVertexHandle| {
            if front == back {
                false
            } else {
                let front_point = triangulation.vertex(front).position();
                let back_point = triangulation.vertex(back).position();
                let front_point = Point2::new(front_point.x, front_point.y);
                let back_point = Point2::new(back_point.x, back_point.y);
                let mut chain = split_vertices
                    .iter()
                    .filter_map(|(handle, point)| {
                        if *handle == front || *handle == back {
                            None
                        } else {
                            boundary_segment_parameter(*point, front_point, back_point)
                                .map(|parameter| (parameter, *handle))
                        }
                    })
                    .collect::<Vec<_>>();
                chain.sort_by(|lhs, rhs| lhs.0.total_cmp(&rhs.0));
                let handles = iter::once(front)
                    .chain(chain.into_iter().map(|(_, handle)| handle))
                    .chain(iter::once(back))
                    .collect::<Vec<_>>();
                let mut did_handle = false;
                handles.windows(2).for_each(|window| {
                    if window[0] != window[1] {
                        did_handle = true;
                        add_direct_cdt_constraint(
                            triangulation,
                            &mut added_constraints,
                            &mut skipped_constraints,
                            window[0],
                            window[1],
                        );
                    }
                });
                did_handle
            }
        };
        self.loops.iter().map(Vec::len).for_each(|len| {
            let range = counter..counter + len;
            counter += len;
            let mut prev: Option<usize> = None;
            range.circular_tuple_windows().for_each(|(i, j)| {
                let Some(vj) = poly2tri[j] else { return };
                if let Some(p) = prev {
                    let Some(v) = poly2tri[p] else { return };
                    if add_constraint(v, vj) {
                        prev = None;
                    }
                } else {
                    let Some(vi) = poly2tri[i] else { return };
                    if !add_constraint(vi, vj) {
                        prev = Some(i);
                    }
                }
            });
        });
        (boundary_map.len(), added_constraints, skipped_constraints)
    }
}
